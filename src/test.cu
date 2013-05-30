__global__ void TestFirstScatter(FLOAT *x, setup *S, curandState *state){
    FLOAT x_in;
    FLOAT k_in_photon[3];
    int status;
    FLOAT a, nu_doppler;
    FLOAT tau_travel;
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    
    x_in = x[id];

    nu_doppler = CONSTANT_NU_DOPPLER*sqrt(S->Temperature/10000.0); /* in cm/s */
    a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
    RND_spherical(k_in_photon, state);    
    status = PropagateStep(&x_in, &(k_in_photon[0]), &tau_travel, &a, state, S);
    x[id] = x_in;
} 



extern "C" void DumpArray(char *filename, FLOAT *x, int n_points){
    FILE *out;
  int i;

  if(!(out=fopen(filename, "w"))){
    fprintf(stderr, "Problem opening file %s\n", filename);
    exit(1);
  }
  
  for(i=0;i<n_points;i++){
    fprintf(out, "%f\n", x[i]);

  }
  fclose(out);
}

/* This is the driver for all the tests to be run on the GPU*/
extern "C" void all_tests(void){
  FLOAT *x_aux;
  FLOAT *x_aux_d;
  setup *S;
  int n_aux;
  int m_seed;
  int nBlocks, blockSize;
  int i;
  char fileout[MAX_FILENAME_SIZE];
  cudaError_t cudaResult = cudaSuccess;
  //allocate and copy the global setup structure 
  cudaMalloc((void**) &S, sizeof(setup));
  cudaMemcpy(S, &All, sizeof(setup), cudaMemcpyHostToDevice);



  //allocate memory on host an device
  n_aux = SIZE_TEST_ARRAY;
  if(!(x_aux = (FLOAT *)malloc(n_aux * sizeof(FLOAT)))){
    fprintf(stderr, "Problem allocating the auxiliary array\n");
    exit(1);
  }
  cudaMalloc((void **) &x_aux_d, n_aux * sizeof(FLOAT));
  
  //initialization
  for(i=0;i<n_aux;i++){
    x_aux[i] = 0.0;
  }

  //from host to device
  cudaMemcpy(x_aux_d, x_aux, sizeof(FLOAT) * n_aux, cudaMemcpyHostToDevice);
 

  //block size and number of blocks
  blockSize = 512; // This is the number of threads inside a block
  nBlocks = (n_aux)/blockSize + (n_aux%blockSize == 0?0:1); // This is the number of blocks
  fprintf(stdout, "nBlocks for testing %d\n", nBlocks);

  //allocate memory for RNG states
  curandState *d_rngStates = 0;
  cudaResult = cudaMalloc((void **)&d_rngStates, blockSize * nBlocks * sizeof(curandState));

  // Initialise RNG
  m_seed = All.RandomSeed;
  initRNG<<<nBlocks, blockSize>>>(d_rngStates, m_seed);

  //Make all the tests
  if(All.TestFirstScatter){
    TestFirstScatter <<< nBlocks, blockSize >>> (x_aux_d, S, d_rngStates);    
    cudaMemcpy(x_aux, x_aux_d, sizeof(FLOAT) * n_aux, cudaMemcpyDeviceToHost);
    sprintf(fileout, "%s/%s", All.OutputDir, All.OutputTestFile);
    DumpArray(fileout, x_aux, n_aux);
  }

}
