#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "vmath.h"
#include "laplacian.h"
#include "geometry.h"


void check_constant() {
	/* This functions checks that the laplacian of a flattened array is 0
	as can be easily calculated analytically. */
	printf("Checking constant vector...\n");
	long int N = 1001;
	unsigned int D = 1;
	long int L = ipow(N, D);
	double complex *arr;
	double complex *res;
	double err;
	arr = malloc(L*sizeof(double complex));
	res = malloc(L*sizeof(double complex));
	printf("Size of double: %d\n", sizeof(double));

	for(long int i=0; i<L; i++) {
		arr[i] = 10.0;
	}
	laplacian(res, arr, N, D);
	err = abs_vec(res, L);
	printf("Absolute magnitude of Laplacian of flat array: %.2e\n", err);
	free(arr);
	free(res);
}

void check_sine() {
	/* In this function the analytical and discrete laplacians of a sum od sines
	are calculated and the relative error will be printed. */
	printf("Checking sum of sines...\n");
	long int N = 10001;
	unsigned int D = 1;
	long int L = ipow(N, D);
	printf("Length of array: %d\n", L);
	double complex *arr, *res, *ref;
	arr = malloc(L*sizeof(double complex)); // input array
	res = malloc(L*sizeof(double complex)); // output array of laplacian
	ref = malloc(L*sizeof(double complex)); // reference array with analytical result
	//printf("mem allocation complete");
	long int *coord;
	coord = malloc(D*sizeof(long int)); // coordinate array
	double omega = 2*M_PI/(N/2.0);
	double err;

	// calculate values of function f(n) = sin(w*n_1) + sin(w*n_2)+...+sin(w*n_(D-1))
	for(long int i=0; i<L; i++) {
		index2coord(coord, i, N, D);
		arr[i] = 0.0;
		for(unsigned int d=0; d<D; d++) {
			arr[i] += sin(omega*coord[d]);
		}
	}

	// calculate reference values (expected result of laplacian)
	for(long int i=0; i<L; i++) {
		index2coord(coord, i, N, D);
		ref[i] = 0.0;
		for(unsigned int d=0; d<D; d++) {
			ref[i] += pow(omega,2.0)*(-1)*sin(omega*coord[d]);
		}
	}

	// next the discrete laplacian is calculated
	laplacian(res, arr, N, D);
	sub_vec(res, ref, res, L);
	err = abs_vec(res, L)/abs_vec(ref, L);
	printf("Relative error magnitude of Laplacian of sine sum array: %.2e\n", err);

	free(arr);
	free(res);
	free(ref);
	free(coord);
}

void check_exp() {
	/* We want to calculate the discretized laplacian analytically for a function f(\vec{omega})=e^{vec{n}\dot\vec{omega}}.
	The expected result is Laplacian(f)(\vec{n})=2*e^{vec{n}\dot\vec{omega}*\sum^{D-1}_{i=0}cos(\omega_i)-1.
	The absolute difference between function result and analytical result is printed as a measure if the laplacian is computed as intended.*/
	printf("Checking analytical solution of exponential function e^(inw)...\n");
	long int N = 101;
	unsigned int D = 3;
	long int L = ipow(N, D);
	double complex *arr = malloc(L*sizeof(double complex)); // input array
	double complex *res = malloc(L*sizeof(double complex)); // output array of laplacian
	double complex *ref = malloc(L*sizeof(double complex)); // reference array with analytical result
	double complex *omega = malloc(D*sizeof(double complex)); // omega vector
	long int *coord = malloc(D*sizeof(long int)); // coordinate array
	double complex *coord_c = malloc(D*sizeof(double complex)); // array to cast coord into double complex
	double complex temp_sum = 0.0;
	double err;
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
		arr[i] = cexp(I*dot_product(coord_c, omega, D));
	}
	//calculate the laplacian with implemented function
	laplacian(res, arr, N, D);
	// calculate the analytical laplacian
	for(long int i=0; i<L; i++) {
		temp_sum = 0;
		for(unsigned int d=0; d<D; d++) {
			temp_sum += cos(omega[d]) - 1;
		}
		ref[i] = 2*arr[i]*temp_sum;
	}
	sub_vec(res, ref, res, L);
	err = abs_vec(res, L);
	printf("Absolute deviation from analytical solution: %.2e\n", err);

	free(arr);
	free(res);
	free(ref);
	free(omega);
	free(coord);
	free(coord_c);;
}

int main() {
	// here all the tests of the laplacian are called
  printf("Starting test of Laplacians\n");
	check_constant();
	check_sine();
	check_exp();
  return 0;
}
