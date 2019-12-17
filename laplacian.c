#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#iclude "structs.h"
#include "nneighbour.h"
#include "vmath.h"

// Function assumes that *out and *in are already allocated as double complex arrays of length L=N**D

void laplacian(double complex *out, double complex *in, parameters params){
    static long int *neighbours = NULL;
    static long int prev_N, prev_D;

    //Check if parameters for neighbours array have changed. If yes, calculate the array again
    if((neighbours == NULL) || (prev_N != params.N) || (prev_D != params.D)){
        free(neighbours);
        params.L = ipow(params.N,params.D);
        neighbours = malloc(2*params.D*params.L*sizeof(long int));
        nneighbour_init(neighbours, params.N, params.D);
        prev_N = params.N;
        prev_D = params.D;
    }

    //Loop over all indices
    for(int i=0; i<params.L; i++){
        // Loop over all dimensions
        for(int d=0; d<params.D; d++){
            // Calculate discrete second derivative according to d^2 psi/ dx^2 = (f(x+h)-2*f(x)+f(x-h))/h^2
            out[i] += (in[neighbours[2*(i*params.D+d)]]-2*in[i]+in[neighbours[2*(i*params.D+d)+1]]);
        }
    }
}
