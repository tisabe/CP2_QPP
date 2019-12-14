#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include "vmath.h"


void matrix_generator(double complex *out, double complex *in, long int L){
  out[0] = 4.0*in[0] + 1.0*in[1];
  out[1] = 1.0*in[0] + 3.0*in[1];
}

void wikipedia_test(){
  //[[4,1],[1,3]]*[x1,x2]=[1,2]

  long int N = 2;
  unsigned int D = 2;
  int max_iter = 10;
  double tol = 0.0000001;
  long int L = ipow(N,D);

  double complex *out = malloc(L*sizeof(double complex));
  double complex *in = malloc(N*sizeof(double complex));
  //double complex *out_func = malloc(L*sizeof(double complex));
  //double complex *in_func = malloc(N*sizeof(double complex));

  in[0] = 1 + 0*I;
  in[1] = 2 + 0*I;

  cg(out, matrix_generator, in, max_iter, tol, L);

  printf("[%lf+%lf i, %lf+%lf i]\n",creal(out[0]),cimag(out[0]),creal(out[1]),cimag(out[1]));

  free(out);
  free(in);
  //free(out_func);
  //free(in_func);
}

int main(){
  wikipedia_test();
  return 0;
}
