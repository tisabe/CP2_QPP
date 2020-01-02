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

<<<<<<< HEAD
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
=======
    
    double complex *x= malloc(L * sizeof(double complex));
double complex *x2= malloc(L * sizeof(double complex));
double complex *x3= malloc(L * sizeof(double complex));
    
    for(int i=0; i<L; i++){
        x[i] = 10.0;
    }
	nfft(x2, x, N, D);
	nfft_inverse(x3, x2, N, D);
for( int i=0; i < L; i++){
		printf("x3[i] \n");
	}
free(x);
free(x2);
free(x3);
>>>>>>> a6b4ff696a934fea1df68846c66244c18e6afa0d
}
