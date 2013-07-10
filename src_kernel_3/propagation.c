#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "struct.h"
#include "propagation.h"
#include "ray.h"
#include "RND_lyman.h"
#include "io.h"
#include "myrand.h"
//#include "propagation_kernel.h"


double LocalLyaH(double x, double a)
{
    double z, P, q, H;
    z = 0.0; P=0.0; q=0.0; H=0.0;

    z = (x*x - 0.855)/(x*x + 3.42);
    
    P = 5.674*(z*z*z*z) - 9.207*(z*z*z)  + 4.421*(z*z) + 0.1117*z;

    if(z<=0)
    {
	q = 0.0;
    }
    else
    {
	q = (1.0 + (21.0/(x*x)))*(a/(PI*(x*x + 1.0)))*P;
    }

    H = q*sqrt(PI);
    H = H + exp(-(x*x));
    return H;
}

double LocalLyaCrossSection(double x, double a){
    double sigma;
    double nu_doppler;


    nu_doppler = Lya_nu_line_width_CGS/(2.0*a);
    sigma = 0.4162*sqrt(PI)*(CHARGEELECTRON*CHARGEELECTRON)/(ELECTRONMASS*C_LIGHT*nu_doppler);
    sigma = sigma*LocalLyaH(x,a);
        
    return sigma;
}

void PropagateAllSetup(void)
/*makes the general problem setup from the input values*/
{
  double a, nu_doppler, ly_sigma_0, column_HI, v_thermal, ly_sigma_v, x_new;
  
  nu_doppler = CONSTANT_NU_DOPPLER*sqrt(All.Temperature/10000.0); /* in s^-1 */
  a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
  ly_sigma_0 =  LocalLyaCrossSection(0,a);
  
  v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
  x_new = 200.0*1.0E5/v_thermal;
  ly_sigma_v = LocalLyaCrossSection(x_new,a);

  column_HI = All.Tau / ly_sigma_0;
  All.SlabLength = All.Tau/(All.NumberDensityHI*ly_sigma_0);
  
#ifdef DEBUG
  fprintf(stdout, "cross section for delta_v=200km/s %e\n", ly_sigma_v);
  fprintf(stdout, "cross section %e\n", ly_sigma_0);
  fprintf(stdout, "ColumnDensinty %e\n", column_HI);
  fprintf(stdout, "Temp %g\n", All.Temperature);
  fprintf(stdout, "a %g\n", a);
  fprintf(stdout, "ly_sigma_0 %g\n", ly_sigma_0);
  fprintf(stdout, "Slab Length [cm] %g\n", All.SlabLength);
  fprintf(stdout, "maximum at (a*tau)^1/3 %g\n", pow(a*All.Tau,0.333));
  fprintf(stdout, "scatterings %g\n", 1.612*All.Tau);
#endif

}

    
void PropagateAll(void)
{
    int n_packages;
    int i;
    lyman_RT_photons *Ph;
    int n_iter;
    char FileName[MAX_FILENAME_SIZE];
    int status;

    /*update some geometrical values*/
    PropagateAllSetup();
    
    /*get the number of photons to propagate*/
    n_packages = All.BusyPhotons;
#ifdef DEBUG
    fprintf(stdout, "%d packages to propagate\n", n_packages);
#endif 

    /*create the packages*/
    Ph = PhotonListCreate(n_packages);
    
    /*initialize the packages*/
    PhotonListInitialize(Ph);
        
    /*propagate each package*/
    status = PropagatePackage(Ph->PosX, Ph->PosY, Ph->PosZ,
			      Ph->DirX, Ph->DirY, Ph->DirZ,
			      Ph->ScatterHI, Ph->x_out, Ph->Active, 
			      n_packages);
    
    /*free the memory*/
    PhotonFree(Ph);
}

void TestRND(void)
{
    double x_in;
    double r_travel;
    double k_in_photon[3];
    double n_HI;
    int status;
    int i;
    char FileName[MAX_FILENAME_SIZE];
    FILE *out;

    sprintf(FileName, "%s/%s_RND.proc.%d", All.OutputDir, All.OutputTestFile, ThisProc);
    if(!(out=fopen(FileName, "w"))){
	fprintf(stderr, "TestRND problem opening file %s\n", FileName);
	exit(0);
    }

    fprintf(out, "%d %e %e\n", N_POINTS_IN_TEST, 0.0, 0.0);
    for(i=0;i<N_POINTS_IN_TEST;i++){
	r_travel = RandFloatUnit();
	fprintf(out, "%e\n", r_travel);
    }

    fclose(out);
} 

void TestAll(void){
    if(All.TestParallelVel){
	RND_lyman_parallel_vel_test(All.Test_x,All.Test_a);
    }

    if(All.TestParallelVelFast){
	RND_lyman_parallel_vel_fast_test(All.Test_x,All.Test_a);
    }    
    /*
    if(All.TestFirstScatter){
	TestFirstScatter(All.Test_x, All.Test_a);
    }
    */

    if(All.TestRND){
	TestRND();
    }

    if(All.TestPerpVel){
       RND_lyman_perp_vel_test();
    }
}

