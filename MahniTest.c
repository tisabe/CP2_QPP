#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <float.h>
#include <math.h>

#include "structs.h"
#include "nfft.h"
#include "vmath.h"

#define _USE_MATH_DEFINES

int main(){

  int N = 16;
  int D = 2;
  int L = ipow(N, D);

  double complex *x= malloc(L * sizeof(double complex));
  double complex *x2= malloc(L * sizeof(double complex));
  double complex *x3= malloc(L * sizeof(double complex));

  for(int i=0; i<L; i++){
    x[i] = 10.0;
  }
  nfft(x2, x, N, D);
  nfft_inverse(x3, x2, N, D);
  for( int i=0; i < L; i++){
		  printf("x[%d]=%f%+fi\n",i ,crealf(x3[i]), cimagf(x3[i]));
  }
  free(x);
  free(x2);
  free(x3);
}
