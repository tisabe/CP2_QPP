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

  for(int i=0; i<L; i++){
    x[i] = 10.0;
  }
  nfft(x, x, N, D);
  nfft_inverse(x, x, N, D);
  for( int i=0; i < L; i++){
		  printf("x[%d]=%f%+fi\n",i ,crealf(x[i]), cimagf(x[i]));
  }
  free(x);
}
