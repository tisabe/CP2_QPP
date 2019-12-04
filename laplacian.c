#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include "nneighbour.h"

// Function assumes that *out and *in are already allocated as double complex arrays of length L=N**D

void laplacian(double complex *out, long int *neighbours, double complex *in, long int N, unsigned int D, long int L){
    /*using the next neighbours saved in a static variable*/
    for(int i=0; i<L; i++){
        // Loop over all dimensions
        for(int d=0; d<D; d++){
            // Calculate discrete second derivative according to d^2 psi/ dx^2 = (f(x+h)-2*f(x)+f(x-h))/h^2
            out[i] += (in[neighbours[(i*2*D+2*d)]]-2*in[i]+in[neighbours[(i*2*D+2*d)+1]])/2.;
        }
    }
}
