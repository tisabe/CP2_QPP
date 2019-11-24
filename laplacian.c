#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include "nneighbour.h"

// Change if ipow() is moved
#include "indices.h"

// Function assumes that *out and *in are already allocated as double complex arrays of length N**D

void laplacian(double complex *out, double complex *in, long int N, unsigned int D){
    // Loop over all indices
    for(int i=0; i<ipow(N,D); i++){
        // Loop over all dimensions
        for(int d=0; d<D; d++){
            // Calculate discrete second derivative according to d^2 psi/ dx^2 = (f(x+h)-2*f(x)+f(x-h))/h^2
            out[i] += (in[nneighbour(i, d, 1, N, D)]-2*in[i]+in[nneighbour(i, d, -1, N, D)])/2.;
        }
    }
}
