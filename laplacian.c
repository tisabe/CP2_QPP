#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h> // compile with -fopenmp

#include "nneighbour.h"
#include "vmath.h"

// Function assumes that *out and *in are already allocated as double complex arrays of length L=N**D

void laplacian(double complex *out, double complex *in, long int N, unsigned int D){
    static long int *neighbours = NULL;
    static long int prev_N, prev_D;
    static long int L;
    //Check if parameters for neighbours array have changed. If yes, calculate the array again
    if((neighbours == NULL) || (prev_N != N) || (prev_D != D)){
        free(neighbours);
        L = ipow(N, D);
        neighbours = malloc(2*D*L*sizeof(long int));
        nneighbour_init(neighbours, N, D);
        prev_N = N;
        prev_D = D;
    }

    //Loop over all indices
    for(int i=0; i<L; i++){
        // Loop over all dimensions
        for(int d=0; d<D; d++){
            // Calculate discrete second derivative according to d^2 psi/ dx^2 = (f(x+h)-2*f(x)+f(x-h))/h^2
            out[i] += (in[neighbours[2*(i*D+d)]]-2*in[i]+in[neighbours[2*(i*D+d)+1]]);
        }
    }
}

void laplacian_mp(double complex *out, double complex *in, long int N, unsigned int D){
    static long int *neighbours = NULL;
    static long int prev_N, prev_D;
    static long int L;
    //Check if parameters for neighbours array have changed. If yes, calculate the array again
    if((neighbours == NULL) || (prev_N != N) || (prev_D != D)){
        free(neighbours);
        L = ipow(N, D);
        neighbours = malloc(2*D*L*sizeof(long int));
        nneighbour_init(neighbours, N, D);
        //printf("neighbours calculated\n");
        prev_N = N;
        prev_D = D;
    }
    long int i;
    unsigned int d;
    #pragma omp parallel for reduction(+:out[:L])
    for(i=0; i<L; i++){
        // Loop over all dimensions
        for(d=0; d<D; d++){
            // Calculate discrete second derivative according to d^2 psi/ dx^2 = (f(x+h)-2*f(x)+f(x-h))/h^2
            out[i] += (in[neighbours[2*(i*D+d)]]-2*in[i]+in[neighbours[2*(i*D+d)+1]]);
        }
    }
}
