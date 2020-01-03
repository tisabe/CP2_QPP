#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

#include "structs.h"
#include "geometry.h"
#include "hamiltonian.h"
#include "vmath.h"

void ran_vec(double complex *vec, long int L) {
  /*Generates a random vector of size L and saves it in *vec.*/
  for (long int i=0; i<L; i++) {
    vec[i] = (double)rand()/RAND_MAX; // sets to random value between 0.0 and 1.0
  }
}
void check_add_ran(long int N, unsigned int D) {
  /* This function checks the additivity of the hamiltonian with a random
  wavefunction and potential given axis length N and number of dimenions D.
  The relative error will be printed as a result.*/
  printf("Testing the additivity of the Hamiltonian with random vectors for N=%d, D=%d\n", N, D);
  long int L = ipow(N, D);
  double complex *v;
  double complex *w;
  double complex *res;
  double complex *ref;
  v = malloc(L*sizeof(double complex));
  w = malloc(L*sizeof(double complex));
  res = malloc(L*sizeof(double complex));
  ref = malloc(L*sizeof(double complex));
  double err = 0;

  parameters params;
  params.N = N;
  params.D = D;
  params.L = L;
  params.tol = 1.0;
  params.max_iter = 1;
  params.tauhat = 1.0;
  params.epsilon = hbar * omega;
  params.mhat = 1e-4;
  params.khat = params.mhat;
  params.a = 1.0;
  params.pot = malloc(params.L*sizeof(double complex));
  ran_vec(params.pot, params.L);

  add_vec(res, v, w, L);                  // res = v + w
  hamiltonian(v, v, params);              // v = H(v)
  hamiltonian(w, w, params);              // w = H(w)
  hamiltonian(res, res, params);          // res = H(res) = H(v + w)
  add_vec(ref, v, w, L);                  // ref = H(v) + H(w)
  sub_vec(ref, ref, res, L);              // ref = ref - res = H(v) + H(w) - H(v + w)
  err = abs_vec(ref, L)/abs_vec(res, L);  // err = abs(ref)
}

void check_hom_ran(long int N, unsigned int D) {
  /* This function checks the homogeneity of the hamiltonian with a random
  wavefunction and potential given axis length N and number of dimenions D.
  The relative error will be printed as a result.*/

}

void check_eigen(long int N, unsigned int D) {
  /*This function checks if standing waves are eigenfunctions for the
  hamiltonian with a potential equal to zero at all points. The periodic
  boundary conditions determine the frequency of the eigenfunctions.*/

}

int main() {



  return 0;
}
