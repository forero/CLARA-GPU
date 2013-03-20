#ifndef STRUCT_H
#define STRUCT_H

/*Internal flag/constant difinitions*/
#define EPSILON_DOUBLE 1E-6
#define MAX_VEL_ITER 1000000
#define EXCEEDED_ITERATIONS 10
#define MAX_FILENAME_SIZE 1024

/*some units in (cgs) sistem*/
#define PI 3.14159265358979323846264338327
#define G_GRAVITY         6.672e-8
#define HUBBLE            3.2407789e-18 /* in h/sec */
#define HUBBLE_TIME       3.09e+17   /*in sec/h*/
#define C_LIGHT           2.9979e+10   /*cm/s*/
#define CM_PER_MPC        3.085678e24
#define HUBBLE            3.2407789e-18	/* in h/sec */
#define LOG10_L_SUN       22.583     /*log10(3.83E22) = L_bol_sun [in ergs/s] */
#define BOLTZMANN         1.3806e-16 /*cgs units*/
#define PROTONMASS        1.6726e-24 /*cgs units*/
#define GAMMA             1.66666 
#define CHARGEELECTRON	  4.8032E-10  /*E.S.E			H.Scheffler/Elsässer Bau und Physik der Galaxies*/
#define ELECTRONMASS	  9.109382616E-28			/*g				wikipedia*/
#define PLANCK            6.626075E-27 

/*Lyman alpha related constants*/
#define Lya_nu_center_CGS 2.466E15
#define Lya_nu_line_width_CGS 9.936E7 
#define Lya_lambda_line_center_CGS 121.6E-7 /* cm wikipedia */
#define CONSTANT_NU_DOPPLER 1.057E11

typedef struct global_setup
{
  /*input files and formats*/
  char InputDir[MAX_FILENAME_SIZE];  /* directory where the configuration file is*/
  
  /*output files*/
  char OutputDir[MAX_FILENAME_SIZE];            /* directory where all the outputs are written*/
  char OutputFile[MAX_FILENAME_SIZE];            /* file were the outputs are written*/
  char OutputTestFile[MAX_FILENAME_SIZE];            /* file were the tests are written*/

  /*output options*/
  int OutputInitList;
  int OutputFinalList;
  int OutputBinary;  
  int RandomSeed;
  
  /*define the problem*/
  int NeufeldSlab;
  int NeufeldCube;
  int ExpandingSphere;
  int RotatingSphere;
  int HomogeneousInit; /*specifies if the photons are to be homogeneously distributed over the volume*/
    
  /*define the tests*/
  int TestParallelVel;
  int TestParallelVelFast;
  int TestFirstScatter;
  int TestPerpVel;
  int TestRND;

  /*parameters controling the algorithm*/
  int NPackages;

  /*Define some physical characteristics of the problem*/
  float Temperature;
  float Tau;
  float InputFrequency;
  float Vmax;
  
  /*parameters for the dust model*/
  double GrainSize;         /*in cm*/
  double TauDust;
  double DustAbsorptionProb;   /*when interacting with dust the probability to be absorbed*/

  /*units*/
  double UnitMass_in_g;  
  double UnitVelocity_in_cm_per_s;
  double UnitLength_in_cm;
  double UnitLymanLuminosity;
  double UnitTime_in_s;
  double UnitDensity_in_cgs;
  double UnitEnergy_in_cgs;
    
  /*parameters to be used in case of test*/
  double Test_a;
  double Test_x;
} setup;

extern setup All;

#endif
