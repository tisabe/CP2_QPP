#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include "vmath.h"

void cg(double complex *out, void (*f)(double complex *, double complex *, long int L), double complex *in, int max_iter, double tol, long int N, unsigned int D);

void matrix_generator(double complex *out, double complex *in, long int L){
  out[0] = 4;
  out[1] = 1;
  out[2] = 1;
  out[3] = 3;
}

void wikipedia_test(){
  //[[4,1],[1,3]]*[x1,x2]=[1,2]

  long int N = 2;
  unsigned int D = 2;
  int max_iter = 10;
  double tol = 0.0000001;
  long int L = ipow(N,D);

  double complex *out = malloc(L*sizeof(double complex));
  double complex *in = malloc(L*sizeof(double complex));
  double complex *out_func = malloc(L*sizeof(double complex));
  double complex *in_func = malloc(L*sizeof(double complex));

  in[0] = 1;
  in[1] = 2;

  cg(out,matrix_generator(out_func, in_func, L), in, max_iter, tol, N, D);

  printf("[%i, %i]\n",out[0],out[1]);

  free(out);
  free(in);
  free(out_func);
  free(in_func);
}

int main(){
  wikipedia_test();
  return 0;
}
