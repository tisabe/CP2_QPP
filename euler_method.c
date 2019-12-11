#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "hamiltonian.h"

/*  Compute the time propagation according to the Euler method.
    Note: Normalisation of the wave function is not conserved!! */



/* ************      To be tested        ************ */

void euler_method(double complex *out, double complex *in, long int N, unsigned int D, long int L, double tauhat, double total_time, double m, double epsilon, int ext_potential_type, double parameter) {
    // Allocate memory for a temporary array of length L
    double complex *ham_out = malloc(L*sizeof(double complex));

    // Copy data from in array to out array (= 0 time steps calculated)
    for(long int i=0; i<L; i++){
        out[i] = in[i];
    }

    // For loop that is being executed until the simulation length specified in total_length (to be given in 1/E_H approx 2.4 * 10^{-17} s) is reached
    for(int t=0; (t * tauhat) < total_time; t++){
        // Calculate H(psi)
        hamiltonian(ham_out, out, N, D, L, m, epsilon, ext_potential_type, parameter);
        // psi^{q+1} = psi^q - i tau H(psi^q)
        for(long int i=0; i<L; i++){
            out[i] -= 1I * tauhat * ham_out[i];
        }
    }

    free(ham_out);
}
