#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "structs.h"
#include "hamiltonian.h"

/*  Compute the time propagation according to the Euler method.
    Note: Normalisation of the wave function is not conserved!! */



/* ************      To be tested        ************ */

void euler_method(double complex *out, double complex *in, parameters params) {
    // Allocate memory for a temporary array of length L
    double complex *ham_out = malloc(params.L*sizeof(double complex));

    // Calculate H(psi)
    hamiltonian(ham_out, out, params);

    // psi^{q+1} = psi^q - i tau H(psi^q)
    for(long int i=0; i<params.L; i++){
        out[i] = in[i] - 1I * params.tauhat * ham_out[i];
    }

    free(ham_out);
}
