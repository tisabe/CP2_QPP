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
#include "crank_nicolson.h"
#include "observables.h"

#define _USE_MATH_DEFINES

int main(){

    parameters params;

    double omega = 8.3e14; //Hz
    double mass_H = 1.67e-27; //kg
    double hbar = 1.054571817e-34; //Js

    params.N = 1001;
    params.D = 1;
    params.L = ipow(params.N, params.D);
    params.tol = DBL_EPSILON;
    params.max_iter = 1000;

    params.tauhat = 1e-8;
    params.epsilon = hbar * omega;
    params.mhat = 5e-4;
    params.khat = params.mhat;
    params.a = hbar * sqrt(params.mhat / (mass_H*params.epsilon));

    double complex *out_wf = malloc(params.L * sizeof(double complex));
    double complex *start_wf = malloc(params.L * sizeof(double complex));
    double complex *x_space = malloc(params.L * sizeof(double complex));
    long int *coords = malloc(params.D * sizeof(long int));
    params.pot = malloc(params.L * sizeof(double complex));

    FILE *startwf_file;
    FILE *potential_file;
    FILE *norm_file;

    for(long int i=0; i<params.L; i++){
        index2coord(coords, i, params.N, params.D);
        x_space[i] = coords[0] * params.a;
    }

    for(long int i=0; i<params.L; i++){
        start_wf[i] = pow(sqrt(params.a),params.D) * sqrt(sqrt(mass_H * omega / (M_PI * hbar))) * exp(-1/2. * mass_H * omega / hbar * x_space[i] * x_space[i]);
    }

    parameters *p = &params;
    gen_pot_harmonic(p, omega);

    startwf_file = fopen("h_test_input.txt","w");
    for(long int i=0; i<params.L; i++){
        fprintf(startwf_file, "%e\n", creal(start_wf[i]));
    }
    fclose(startwf_file);

    norm_file = fopen("h_norm_output.txt","w");
    double norm_obs;

    for(long int t=0; (t * params.tauhat) < 1e-1; t++){
        if(t % 100000 == 0){
	    norm_obs = cabs(obs_norm(start_wf, params));
            printf("%li\t%e\n", t, norm_obs - 1);
	    fprintf(norm_file, "%li\t%e\n", t, norm_obs - 1);
        }
        euler_method(out_wf, start_wf, params);
        assign_vec(start_wf, out_wf, params.L);
    }

    fclose(norm_file);

    potential_file = fopen("h_potential.txt","w");

    for(long int i=0; i<params.L; i++){
        fprintf(potential_file, "%e\n", creal(params.pot[i]));
    }

    printf("a = %e", params.a);

    fclose(potential_file);

    free(coords);
    free(start_wf);
    free(x_space);
    free(params.pot);
    //free(out_wf);

    return 0;
}
