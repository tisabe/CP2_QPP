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
void nfft_inverse(double complex *out, double complex *in, int N, int D)

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
         gsl_fft_complex_radix2_inverse(data,Nj,N);
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


