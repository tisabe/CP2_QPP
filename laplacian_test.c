#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vmath.h"
#include "laplacian.h"
#include "indices.h"


void check_constant() {
	printf("Checking constant vector\n");
	long int N = 101;
	unsigned int D = 3;
	long int L = ipow(N, D);
	double complex *arr;
	double complex *res;
	double err;
	arr = malloc(L*sizeof(double complex));
	res = malloc(L*sizeof(double complex));
	for(long int i=0; i<L; i++) {
		arr[i] = 10.0;
	}
	laplacian(res, arr, N, D);
	err = abs_vec(res, L);
	printf("Absolute magnitude of Laplacian of flat array: %.2e\n", err);
}

void check_sine() {
	printf("Checking sum of sines\n");
	long int N = 101;
	unsigned int D = 3;
	long int L = ipow(N, D);
	printf("Length of array: %d\n", L);
	double complex *arr, *res, *ref;
	arr = malloc(L*sizeof(double complex)); // input array
	res = malloc(L*sizeof(double complex)); // output array of laplacian
	ref = malloc(L*sizeof(double complex)); // reference array with analytical result
	printf("mem allocation complete");
	long int *coord;
	coord = malloc(D*sizeof(long int)); // coordinate array
	double omega = 2*M_PI/(N/2.0);
	double err;
	
	// calculate values of function
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
	
	laplacian(res, arr, N, D);
	sub_vec(res, ref, res, L);
	err = abs_vec(res, L);
	printf("Absolute error magnitude of Laplacian of sine sum array: %.2e\n", err);
	
	free(arr);
	free(res);
	free(ref);
	free(coord);
}

int main(){
    printf("Starting test of Laplacians\n");
	check_constant();
	check_sine();
    return 0;
}














