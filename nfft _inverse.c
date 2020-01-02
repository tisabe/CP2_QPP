#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>


/*******************************************************************************
Compile with:
gcc -o nfft nfft.c -lgsl -lm -O3
*******************************************************************************/


/*******************************************************************************
void nfft(double complex *out, double complex *in, int N, int D)

It calculates the multidimensional DFT of "in" using the FFT radix2 GSL
implementation, and it stores it in the array pointed by "out".
IMPORTANT: This function assumes lexicographic ordering of the lattice
points!!!

If "ik" is the index associated to the lattice point with coordinates
(k[0],k[1],...,k[D-1]) and "ix" is the index associated to the lattice point
with coordinates (x[0],x[1],...,x[D-1]), then this function calculates

   out[ik] = sum_ix in[ix] * exp( - I 2 pi k.x / N )

In the above formula, k.x denotes the scalar product

   k.x = k[0]*x[0]+k[1]*x[1]+...+k[D-1]*x[D-1]

A description of the parameters follows.

double complex *out  Pointer to the array containing the DFT transform of "in"
                     (output). The function assumes that the array pointed by
                     "out" has been already allocated.

int N                Number of lattice points in each direction.

int D                Number of dimensions
*******************************************************************************/
void nfft(double complex *out, double complex *in, int N, int D)
{
   /****************************************************************************
   VOLUME = N^D
   ****************************************************************************/
   int VOLUME=1;
   for(int j=0;j<D;j++) VOLUME*=N;

   memcpy(out,in,sizeof(double complex)*VOLUME);

   int Nj=1;
   for(int j=0;j<D;j++)
   {
      /*************************************************************************
      At this point, Nj=N^j
      *************************************************************************/
      
      for(int n=0;n<VOLUME/N;n++)
      {
         /**********************************************************************
         In this loop, the variable "index" runs over all lattice points with
         fixed x[j]=0 (where x is the coordinate array). Notice that there are
         VOLUME/N of such points.
         **********************************************************************/
         int index=n%Nj+N*Nj*(n/Nj);
         gsl_complex_packed_array data=(double*)(out+index);
         gsl_fft_complex_radix2_invserse(data,Nj,N);
      }
      
      Nj*=N;
   }
}


/*******************************************************************************
It takes the index identifying a lattice point and returns its coordinates.

int *x          Pointer to the array of coordinates of the lattice point
                (output). The function assumes that x has been already
                allocated.
        
int index       Integer between 0 and N^D-1 which identifies a lattice point
                uniquely (index).
        
int N           Number of lattice points in each direction.

int D           Number of dimensions
*******************************************************************************/
static void index2coord(int *x, int index, int N, int D)
{
   for(int k=0;k<D;k++)
   {
      x[k]=index%N;
      index/=N;
   }
}


/*******************************************************************************
It compares the DFT calculated with the "nfft" function and some simple
self-made implementation.

int N           Number of lattice points in each direction.

int D           Number of dimensions
*******************************************************************************/
static void check_nfft(int N, int D)
{
   int VOLUME=1;
   for(int k=0;k<D;k++) VOLUME*=N;
   
   double complex *in,*gslout,*selfout;
   in=malloc(sizeof(double complex)*VOLUME*3);
   gslout=in+VOLUME;
   selfout=gslout+VOLUME;
   
   /****************************************************************************
   Generation of random input
   ****************************************************************************/
   for(int ix=0;ix<VOLUME;ix++)
   {
      in[ix]= ((1.0*rand())/RAND_MAX-0.5) + I*((1.0*rand())/RAND_MAX-0.5);
   }
   
   /****************************************************************************
   Inefficient and straightforward implementation of n-dimensional DFT (with the
   GSL definition).
   ****************************************************************************/
   double complex *z;
   z=malloc(sizeof(double complex)*N);
   for(int p=0;p<N;p++) z[p]=cexp(-I* (2.0*M_PI*p)/N);
   
   int *k,*x,kx;
   k=malloc(sizeof(int)*D*2);
   x=k+D;
   for(int ik=0;ik<VOLUME;ik++)
   {
      index2coord(k,ik,N,D);
      selfout[ik]=0.0;
      for(int ix=0;ix<VOLUME;ix++)
      {
         index2coord(x,ix,N,D);
         
         kx=0;
         for(int j=0;j<D;j++) kx+=k[j]*x[j];
         
         selfout[ik]+=in[ix]*z[kx%N];
      }
   }

   /****************************************************************************
   Calculation of the n-dimensional DFT with the GSL library
   ****************************************************************************/
   nfft(gslout,in,N,D);
   

   /****************************************************************************
   Comparison of the two implementation
   ****************************************************************************/
   double maxerr,d;
   maxerr=0.0;
   for(int ik=0;ik<VOLUME;ik++)
   {
      d=cabs(gslout[ik]-selfout[ik]);
      if(d>maxerr) maxerr=d;
   }
   
   printf("Comparison for N= %3d , D= %d : err/VOLUME= %e  (should be < 1e-15)\n",N,D,maxerr/VOLUME);
   
   
   free(in);
   free(z);
   free(k);
}


int main(int argc, char *argv[])
{
   printf(
   "\n"
   "This program contains the definition of the \"nfft\" function, which\n"
   "calculates the multidimensional Discrete Fourier Transform using the FFT\n"
   "implementation of the Gnu Scientific Library (GSL). You can extract this\n"
   "function and use it in your project.\n"
   "\n"
   "As a simple check of the implementation, this program compares the output\n"
   "of the \"nfft\" function with a simple but inefficient implementation of\n"
   "the DFT. The results of the comparison are printed below, for a selection\n"
   "of values of N and D.\n"
   "\n");
   
   check_nfft(4,1);
   check_nfft(8,1);
   check_nfft(16,1);
   check_nfft(8,2);
   check_nfft(16,2);
   check_nfft(8,3);
   check_nfft(16,3);
   check_nfft(8,4);
   
   printf("\n");
   
   return 0;
}