 #include "nfft3conf.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
/*#ifdef HAVE_COMPLEX_H*/
#include <complex.h>
/*#endif*/


#include "nfft3util.h"
#include "nfft3.h"

static void our_test_nfft_2d(double* x, double* y, double* f, int nCoeffs, int nNodes, double* kx, double* ky, double* fhat, int iter)
{
  nfft_plan p;
  solver_plan_complex ip;
  int i, j;
  int nIterMax = iter;

  /** init a two dimensional plan */
  nfft_init_2d(&p,nCoeffs,nCoeffs,nNodes);


  for(i=0; i<nNodes; i++)
  {
	p.x[2*i]  = x[i];
	p.x[2*i+1] = y[i];
//	p.f[i]  = f[i];
  }


  /** precompute psi, the entries of the matrix B */
  if(p.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p);

  /* Initialize complex */
  //solver_init_advanced_double(&ip, (nfft_mv_plan_complex*)(&p), CGNR);
  solver_init_advanced_complex(&ip, (nfft_mv_plan_complex*)(&p), CGNR);
  printf("%d %d\n", p.M_total, nNodes);

  /* Copy from input */
  for(i=0; i<nNodes; i++)
  {
	ip.y[i] = f[i];
  }

  /* init some guess */
  for(i=0;i<p.N_total;i++)
        ip.f_hat_iter[i]=0.;

  /* inverse system solve */
  solver_before_loop_complex(&ip);
  for(i=0;i<nIterMax;i++)
  {
	solver_loop_one_step_complex(&ip);
	printf("Residual ||r||=%e,\n",sqrt(ip.dot_r_iter));
	if(sqrt(ip.dot_r_iter)<1e-30)	break;
  }
  printf("System solved.\n");
  /* Copy to output */
  for(i = 0; i < nCoeffs; i++)
  {
	for(j = 0; j < nCoeffs; j++)
	{
		ky[i+j*(nCoeffs)] = (double)(-nCoeffs)/2.0 + (double)i;
		kx[i+j*(nCoeffs)] = (double)(-nCoeffs)/2.0 + (double)j;
		fhat[i+j*(nCoeffs)] = cabs(ip.f_hat_iter[i+j*(nCoeffs)]);
	}
  }
  printf("Copied to output.\n");

  /* Finalize */
  solver_finalize_complex(&ip);
  nfft_finalize(&p);
  printf("Finalized!\n");


}

int main(int argc, char** argv)
{

  char filename[256];
  double* x, *y, *f;
  double* kx, *ky, *fhat;
  float testx, testy, testf;
  int readDataFromFile = 0;
  int nCoeffs;
  int nNodes;
  int i, j;
  int iter;
  FILE* fidFunction, *fidTransform, *fidInput;
  double xmin, xmax;


  if(argc==5)
  {
  	printf("Using an input file for the function to transform...\n");
  	nNodes  = atoi(argv[1]);
  	nCoeffs = atoi(argv[2]);
	strcpy(filename,argv[3]);
        iter = atoi(argv[4]);
	readDataFromFile = 1;
	printf("nNodes = %d\nnCoeffs = %d\ninput file = %s\n\n",nNodes,nCoeffs,filename);
  }
  else
  {
  	nNodes = 128;
  	nCoeffs = 128;
  	readDataFromFile = 0;
        iter = 100;
  }


  /* Allocate arrays */
  x    = malloc(nNodes*nNodes*sizeof(double));
  y    = malloc(nNodes*nNodes*sizeof(double));
  f    = malloc(nNodes*nNodes*sizeof(double));
  kx   = malloc(nCoeffs*nCoeffs*sizeof(double));
  ky   = malloc(nCoeffs*nCoeffs*sizeof(double));
  fhat = malloc(nCoeffs*nCoeffs*sizeof(double));


  if(readDataFromFile == 1)
  {
  	/* Open input file */
  	fidInput = fopen(filename,"r");

    /* Read data from file */
    for(i=0; i<nNodes*nNodes; i++)
    {
      //fscanf(fidInput,"%15.15E\t%15.15E",&x[i],&y[i]);
      fscanf(fidInput,"%f\t%f\t%f",&testx,&testy,&testf);
      x[i] = (double)testx;
      f[i] = (double)testf;
      y[i] = (double)testy;
    }

    /* Close input file */
    fclose(fidInput);

  }
  else
  {
    xmin = 0;
    xmax = 1;

  	/** Create equally spaced points in x **/
  	for(i=0; i<nNodes; i++)
  	{
	  for(j = 0; j < nNodes; j++)
	  {
	    x[i+j*(nNodes)] = xmin + (xmax-xmin)/(nNodes-1)*i;
	    y[i+j*(nNodes)] = xmin + (xmax-xmin)/(nNodes-1)*j;
	    f[i+j*(nNodes)] = 10*sin(2*PI*2*x[i+j*(nNodes-1)]);
	  }
  	}
  }


  /* Compute the FFT */
  our_test_nfft_2d(x,y,f,nCoeffs,nNodes*nNodes,kx, ky,fhat, iter);

  /* Write function to file */
  fidFunction  = fopen("function.txt",  "w");
  fidTransform = fopen("transform.txt", "w");
  for(i=0; i<nNodes*nNodes; i++)
  {
    fprintf(fidFunction,  "%E\t%E\t%E\n", x[i], y[i], f[i]);
  }
  for(i=0; i<nCoeffs*nCoeffs; i++)
  {
  	fprintf(fidTransform, "%E\t%E\t%E\n", kx[i], ky[i], fhat[i]);
  }
  fclose(fidFunction);
  fclose(fidTransform);

  /* Cleanup allocated arrays */
  free(x);
  free(y);
  free(f);
  free(kx);
  free(ky);
  free(fhat);



  return 1;
}
