/* **** This function calculates the discrete symmetric laplacian of the input wave function ****
Valid inputs:
        *out: output vector with laplacian of wave function		complex double array of size N**D
		*in: input vector of discretized wavefunction		    complex double array of size N**D
		N: length of each axis						            long int
		D: number of dimensions						            unsigned int

*/
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

#include "structs.h"
#include "geometry.h"
#include "vmath.h"

// Function assumes that *out and *in are already allocated as double complex arrays of length L=N**D

void laplacian(double complex *out, double complex *in, long int N, unsigned int D){
    /*Use the static array neighbours to save all next neighbours in this array.
    There are 2*D next neighbours per point in the lattice, all next neighbours are saved in this array sequentially
    starting with the nextneighbour in 1. dimension in -1 direction, then +1 direction, then 2. dimension -1 direction, ...*/
    static long int *neighbours = NULL;
    static long int prev_N, prev_D, L;

    //Check if geometry parameters have changed. If yes, free, reallocate and recalculate the next neighbours the array again
    if((neighbours == NULL) || (prev_N != N) || (prev_D != D)){
        free(neighbours);
        L = ipow(N, D);
        neighbours = malloc(2*D*L*sizeof(long int));
        nneighbour_init(neighbours, N, D);
        prev_N = N;
        prev_D = D;
    }

    //Loop over all points in the lattice
    for(long int i=0; i<L; i++){
        out[i] = 0; // initialize current element to zero
        // Loop over all dimensions
        for(int d=0; d<D; d++){
            // Calculate discrete second derivative according to d^2 psi/ dx^2 = (f(x+h)-2*f(x)+f(x-h))/h^2
            out[i] += in[neighbours[2*(i*D+d)]]-2*in[i]+in[neighbours[2*(i*D+d)+1]];
        }
    }
}

void laplacian_slow(double complex *out, double complex *in, long int N, unsigned int D){
    /*Use the static array neighbours to save all next neighbours in this array.
    There are 2*D next neighbours per point in the lattice, all next neighbours are saved in this array sequentially
    starting with the nextneighbour in 1. dimension in -1 direction, then +1 direction, then 2. dimension -1 direction, ...*/
    static long int *neighbours = NULL;
    static long int prev_N, prev_D, L;

    //Check if geometry parameters have changed. If yes, free, reallocate and recalculate the next neighbours the array again
    if((neighbours == NULL) || (prev_N != N) || (prev_D != D)){
        free(neighbours);
        L = ipow(N, D);
        neighbours = malloc(2*D*L*sizeof(long int));
        nneighbour_init(neighbours, N, D);
        prev_N = N;
        prev_D = D;
    }

    //Loop over all points in the lattice
    for(long int i=0; i<L; i++){
        out[i] = 0; // initialize current element to zero
        // Loop over all dimensions
        for(int d=0; d<D; d++){
            // Calculate discrete second derivative according to d^2 psi/ dx^2 = (f(x+h)-2*f(x)+f(x-h))/h^2
            out[i] += in[nneighbour(i,d,-1,N,D)]-2*in[i]+in[nneighbour(i,d,+1,N,D)];
        }
    }
}
