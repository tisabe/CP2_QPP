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
  long int L = ipow(N, D);
  double complex *v;
  double complex *w;
  double complex *res;
  double complex *ref;
  v = malloc(L*sizeof(double complex));
  w = malloc(L*sizeof(double complex));
  res = malloc(L*sizeof(double complex));
  ref = malloc(L*sizeof(double complex));
  ran_vec(v, L);
  ran_vec(w, L);
  double err = 0;

  parameters params;
  params.N = N;
  params.D = D;
  params.L = L;
  params.mhat = 1e-4;
  params.pot = malloc(params.L*sizeof(double complex));
  ran_vec(params.pot, params.L);

  add_vec(res, v, w, L);                  // res = v + w
  hamiltonian(v, v, params);              // v = H(v)
  hamiltonian(w, w, params);              // w = H(w)
  hamiltonian(res, res, params);          // res = H(res) = H(v + w)
  add_vec(ref, v, w, L);                  // ref = H(v) + H(w)
  sub_vec(ref, ref, res, L);              // ref = ref - res = H(v) + H(w) - H(v + w)
  err = abs_vec(ref, L)/abs_vec(res, L);  // err = abs(ref)
  printf("N=%ld\t\tD=%d\tRelative error: %e\n", N, D, err);

  free(v);
  free(w);
  free(res);
  free(ref);
  free(params.pot);
}

void check_hom_ran(long int N, unsigned int D) {
  /* This function checks the homogeneity of the hamiltonian with a random
  wavefunction and potential given axis length N and number of dimenions D.
  The relative error will be printed as a result.*/
  long int L = ipow(N, D);
  double complex *v;
  double complex *w;
  double complex *res;
  double complex *ref;
  v = malloc(L*sizeof(double complex));
  res = malloc(L*sizeof(double complex));
  ref = malloc(L*sizeof(double complex));
  ran_vec(v, L);
  double err = 0;

  parameters params;
  params.N = N;
  params.D = D;
  params.L = L;
  params.mhat = 1e-4;
  params.pot = malloc(params.L*sizeof(double complex));
  ran_vec(params.pot, params.L);

  double complex s = (double)rand()/RAND_MAX + 1I*(double)rand()/RAND_MAX; // generates a random scalar multiplier

  scalar_vec(res, v, s, L);               // res = s*v
  hamiltonian(res, res, params);          // res = H(res) = H(s*v)
  hamiltonian(ref, v, params);            // ref = H(v)
  scalar_vec(ref, ref, s, L);           // ref = s*ref = s*H(v)
  sub_vec(ref, ref, res, L);              // ref = ref - res = s*H(v) - H(s*v)
  err = abs_vec(ref, L)/abs_vec(res, L);  // err = abs(ref)
  printf("N=%ld\t\tD=%d\tRelative error: %e\n", N, D, err);

  free(v);
  free(res);
  free(ref);
  free(params.pot);
}

void check_eigen(long int N, unsigned int D) {
  /*This function checks if standing waves are eigenfunctions for the
  hamiltonian with a potential equal to zero at all points. The periodic
  boundary conditions determine the frequency of the eigenfunctions.
  This test is based on the check_exp function in laplacian_test.c,
  the formulas are explained there. Here we include the additional factor
  -1/(2*mhat)*/
  long int L = ipow(N, D);
  double complex *v;
  double complex *w;
  double complex *res;
  double complex *ref;
  v = malloc(L*sizeof(double complex));
  res = malloc(L*sizeof(double complex));
  ref = malloc(L*sizeof(double complex));
  double err = 0;

  parameters params;
  params.N = N;
  params.D = D;
  params.L = L;
  params.mhat = 1e-4;
  params.pot = malloc(params.L*sizeof(double complex));
  for (long int i=0; i<params.L; i++) {
    params.pot[i] = 0.0;
  }

  long int *coord = malloc(D*sizeof(long int)); // coordinate array
	double complex *coord_c = malloc(D*sizeof(double complex)); // array to cast coord into double complex
  double complex *omega = malloc(D*sizeof(double complex)); // omega vector
  // initialze omega with random values from 0 to 2*pi
	for(unsigned int i=0; i<D; i++) {
		omega[i] = (rand()%N)*2*M_PI/N;
	}

	// calculate exponentials depending on coordinates and omega
	for(long int i=0; i<L; i++) {
		index2coord(coord, i, N, D);
		for(unsigned int d=0; d<D; d++) {
			coord_c[d] = (double complex) coord[d]; // cast coordinates to double complex and save in coor_c
		}
		v[i] = cexp(I*dot_product(coord_c, omega, D));
	}

  // calculate the analytical hamiltonian
  double complex temp_sum = 0.0 + 1I*0.0;
  for(long int i=0; i<L; i++) {
		temp_sum = 0;
		for(unsigned int d=0; d<D; d++) {
			temp_sum += cos(omega[d]) - 1;
		}
		ref[i] = (double complex)(-1 / (2*params.mhat))*2*v[i]*temp_sum; // factor of kinetic part included
	}

  hamiltonian(res, v, params);
  sub_vec(ref, ref, res, L);
  err = abs_vec(ref, L)/abs_vec(res, L);  // err = abs(ref)
  printf("N=%ld\t\tD=%d\tRelative error: %e\n", N, D, err);

  free(v);
  free(res);
  free(ref);
  free(params.pot);
  free(coord_c);
  free(coord);
}

int main() {
  unsigned int arr_D[5] = {1, 1, 1, 2, 3};
  long int arr_N[5] = {101, 10001, 1000001, 1001, 101};
  printf("Testing the additivity of the Hamiltonian with random vectors...\n");
  for (int i=0; i<5; i++) {
    check_add_ran(arr_N[i], arr_D[i]);
  }
  printf("Testing the homogeneity of the Hamiltonian with random vectors...\n");
  for (int i=0; i<5; i++) {
    check_hom_ran(arr_N[i], arr_D[i]);
  }
  printf("Testing the analytical solution of the Hamiltonian...\n");
  for (int i=0; i<5; i++) {
    check_eigen(arr_N[i], arr_D[i]);
  }

  return 0;
}
