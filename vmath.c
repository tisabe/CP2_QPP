/* this file will contain all math functions for vector operations
required in the QPP problem*/

#include <complex.h>

long int ipow(long int base, unsigned int exp){
/*integer power function so math does not need to be included*/
/*peer reviewed by frohlofl*/
    int result = 1;
    for (unsigned int i=0; i<exp; i++){
      result *= base;
    }
    return result;
}

double complex dot_product(double complex *v, double complex *w, long int N, unsigned int D) {
/*returns the complex dot product of two vectors v, w*/
  double complex sprod;
  sprod=0.0;
  long int maxi = ipow(N,D);

  for(long int i=0; i<maxi; i++) {
    sprod+=(~(v[i]))*w[i];
  }
  return sprod;
}

void scalar_vec(double complex *out, double complex *vec, double complex *s, long int N, unsigned int D){
/*returns the product of a scalar s and a vector vec*/
  long int maxi = ipow(N,D);
  for(long int i=0; i<maxi; i++) {
    out[i] = s*vec[i];
  }
}

void add_vec(double complex *out, double complex *w double complex *v, long int N, unsigned int D) {
/*returns the sum of two complex vectors v, w*/
  long int maxi = ipow(N,D);

  for(long int i=0; i<maxi; i++) {
    out[i] = v[i]+w[i];
  }
}

void sub_vec(double complex *out, double complex *w double complex *v, long int N, unsigned int D) {
/*returns the difference of two complex vectors v, w*/
  long int maxi = ipow(N,D);

  for(long int i=0; i<maxi; i++) {
    out[i] = v[i]-w[i];
  }
}

void assign_vec(double complex *out, double complex *in, long int N, unsigned int D) {
/*copy the values of in to values of out*/
  long int maxi = ipow(N,D);

  for(long int i=0; i<maxi; i++) {
    out[i] = in[i];
  }
}

void cg(double complex *out, void (*f)(double complex *, double complex *), double complex *in, int max_iter, double tol, long int N, unsigned int D) {
/*this will perform the conjugate gradient algorithm to find x for y=f(x),
where f is a positive, "matrix-like" function. y is passed as in,
f as a function pointer *f and x is saved in out. The maximum number of
iterations is passed as max_iter and the maximum tolerance as tol*/
  long int maxi = ipow(N,D);
}
