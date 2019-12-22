#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <float.h>
#include <math.h>

#include "structs.h"
#include "hamiltonian.h"
#include "vmath.h"
#include "geometry.h"
#include "euler_method.h"

#define _USE_MATH_DEFINES

int main(){

    parameters params;

    double omega = 8.3e14;
    double mass_H = 1.7e-27;
    double hbar = 1.054571817e-34;

    params.N = 101;
    params.D = 1;
    params.L = ipow(params.N, params.D);
    params.tau = 0.1;
    params.tol = DBL_EPSILON;
    params.max_iter = 1000;
    params.total_time = 10;
    params.epsilon = hbar * omega;
    params.parameter = 1;
    params.ext_potential_type = 0;
    params.mhat = 1;

    double a = hbar * sqrt(2*params.mhat/(mass_H * params.epsilon));

    double complex *out = malloc(params.L * sizeof(double complex));
    double complex *start_wf = malloc(params.L * sizeof(double complex));
    double complex *x_space = malloc(params.L * sizeof(double complex));
    long int *coords = malloc(params.D * sizeof(long int));

    FILE *output_file;

    for(long int i=0; i<params.L; i++){
        index2coord(coords, i, params.N, params.D);
        x_space[i] = coords[0] * a;
    }

    for(long int i=0; i<params.L; i++){
        start_wf[i] = sqrt(sqrt(mass_H * omega / (M_PI * hbar))) * exp(-1/2. * mass_H * omega / hbar * x_space[i] * x_space[i]);
    }

    euler_method(out, start_wf, params);

    output_file = fopen("h_test_input.txt","w");

    for(long int i=0; i<params.L; i++){
        printf("%e%+e i\n", creal(out[i]), cimag(out[i]));
        fprintf(output_file, "%e%+e i\n", creal(start_wf[i]), cimag(start_wf[i]));
    }

    fclose(output_file);

    free(coords);
    free(out);
    free(start_wf);
    free(x_space);

    return 0;
}
