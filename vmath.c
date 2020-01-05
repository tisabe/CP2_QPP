/* this file will contain all math functions for vector operations
required in the QPP problem*/

#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "structs.h"

long int ipow(long int base, unsigned int exp) {
/*integer power function so math does not need to be included*/
/*peer reviewed by frohlofl*/
    int result = 1;
    for (unsigned int i=0; i<exp; i++){
      result *= base;
    }
    return result;
}

double complex dot_product(double complex *v, double complex *w, long int L) {
/*returns the complex dot product of two vectors v, w*/
  double complex sprod = 0.0;
  for(long int i=0; i<L; i++) {
    sprod+=(~(v[i]))*w[i];
  }
  return sprod;
}

void scalar_vec(double complex *out, double complex *vec, double complex s, long int L) {
/*returns the product of a scalar s and a vector vec*/
  for(long int i=0; i<L; i++) {
    out[i] = s*vec[i];
  }
}

void mul_compl_vec_real_scal(double complex *out, double complex *vec, double s, long int L) {
/*same as previous function but with real scalar*/
  for(long int i=0; i<L; i++){
    out[i] = s*vec[i];
  }
}

void mul_element(double complex *out, double complex *w, double complex *v, long int L) {
/*returns the elementwise multiplication of two complex vectors v, w*/
  for(long int i=0; i<L; i++) {
    out[i] = v[i]*w[i];
  }
}

void add_vec(double complex *out, double complex *w, double complex *v, long int L) {
/*returns the sum of two complex vectors v, w*/
  for(long int i=0; i<L; i++) {
    out[i] = v[i]+w[i];
  }
}

void sub_vec(double complex *out, double complex *v, double complex *w, long int L) {
/*returns the difference of two complex vectors v, w*/
  for(long int i=0; i<L; i++) {
    out[i] = v[i]-w[i];
  }
}

void assign_vec(double complex *out, double complex *in, long int L) {
/*copy the values of in to values of out*/
  for(long int i=0; i<L; i++) {
    out[i] = in[i];
  }
}

void set_zero(double complex *in, long int L) {
/*set all values of in to zero*/
  for(long int i=0; i<L; i++) {
    in[i] = 0.0;
  }
}

double abs_vec(double complex *in, long int L) {
  return sqrt(creal(dot_product(in, in, L)));
}

long int fact(int n){
    long int res = 1;
    for(int i=1; i<=n; i++){
        res *= i;
    }
    return res;
}

void mat_vec_mul(double complex *out, double complex *mat, double complex *vec, long int N){ //N is the length of the vector
    set_zero(out, N);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            out[i] += mat[N*i+j] * vec[j];
        }
    }
}

void hermtransp_quadr_matrix(double complex *out, double complex *in, long int N){ //N is the length of one size
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            out[i*N+j] = (~(in[j*N+i]));
        }
    }
}

void cg(double complex *out, void (*f)(double complex */*out*/, double complex */*in*/, parameters /*params*/), double complex *in, parameters params) {
/*this will perform the conjugate gradient algorithm to find x for y=f(x),
where f is a positive, "matrix-like" function. y is passed as in,
f as a function pointer *f and x is saved in out. The maximum number of
iterations is passed as max_iter and the maximum tolerance as tol

parameters:   input:
                out: location where to save the result x of the cg for y=f(x), double complex *
                (*f): function which to perform on input vectors, needs to be symmetric and positive-definite (in matrix representation)
                  out(f): location where to save output of functions, double complex *
                  in(f): input vector of function f, double complex *
                  L(f): length of input vector, long int
                in: input vector y for the cg, double complex *
                params: input parameters, struct parameters
*/
  long int L = params.L;
  double complex * x = malloc(L * sizeof(double complex));
  double complex * x_next = malloc(L * sizeof(double complex));
  double complex * r =  malloc(L * sizeof(double complex));
  double complex * r_next = malloc(L * sizeof(double complex));
  double complex * z = malloc(L * sizeof(double complex));
  double complex * d = malloc(L * sizeof(double complex));
  double complex * temp = malloc(L * sizeof(double complex));
  double complex alpha = 0.0;
  double complex beta = 0.0;
  int k = 0; // current iteration

  /*initialize values of arrays*/
  set_zero(z, L);
  set_zero(x, L);
  set_zero(x_next, L);
  set_zero(r_next, L);
  (*f)(z, x, params);
  sub_vec(r, in, z, L);
  assign_vec(d, r, L);

  while((k < params.max_iter) && (abs_vec(r, L) > params.tol)) {
    (*f)(z, d, params);
    alpha = dot_product(r_next, r_next, L)/dot_product(r, z, L);
    scalar_vec(temp, d, alpha, L);
    add_vec(x_next, x, temp, L);
    scalar_vec(temp, z, alpha, L);
    sub_vec(r_next, r, temp, L);
    beta = dot_product(r_next, r_next, L)/dot_product(r, r, L);
    scalar_vec(temp, d, beta, L);
    add_vec(d, r_next, temp, L);
    // update old values r and x
    assign_vec(x, x_next, L);
    assign_vec(r, r_next, L);
    k++;
  }
  if (k >= params.max_iter) {
    printf("Maximum number of iterations (%d) reached, aborting calculation at precision (relative error) %e\n", k, abs_vec(r,L));
  }/*else{
    printf("CG successfully terminated after %i iterations.\n", k);
  }*/

  assign_vec(out, x, L);

  free(x);
  free(x_next);
  free(r);
  free(r_next);
  free(z);
  free(d);
  free(temp);
}
