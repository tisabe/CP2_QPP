#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <float.h>
#include <math.h>

#include "structs.h"
#include "hamiltonian.h"
#include "vmath.h"
#include "geometry.h"
#include "integrators.h"
#include "observables.h"

#define _USE_MATH_DEFINES

int main(){

    parameters params;

    double omega = 8.3e14; //Hz
    double mass_H = 1.67e-27; //kg
    double hbar = 1.054571817e-34; //Js

    params.N = 5001;
    params.D = 1;
    params.L = ipow(params.N, params.D);
    params.tol = DBL_EPSILON;
    params.max_iter = params.L;

    params.tauhat = 1e-3;
    params.epsilon = hbar * omega;
    params.mhat = 1e-4;
    params.khat = params.mhat;
    params.a = hbar * sqrt(params.mhat / (mass_H*params.epsilon));

    double simulation_duration = 1e1;
    int number_time_steps = 100;

    double complex *out_wf = malloc(params.L * sizeof(double complex));
    double complex *start_wf = malloc(params.L * sizeof(double complex));
    double complex *x_space = malloc(params.L * sizeof(double complex));
    long int *coords = malloc(params.D * sizeof(long int));
    params.pot = malloc(params.L * sizeof(double complex));

    FILE *startwf_file;
    FILE *potential_file;
    FILE *norm_file;
    FILE *output_file;

    for(long int i=0; i<params.L; i++){
        index2coord(coords, i, params.N, params.D);
        x_space[i] = coords[0] * params.a;
    }

    for(long int i=0; i<(params.L); i++){
        start_wf[i] = pow(sqrt(params.a),params.D) * sqrt(sqrt(mass_H * omega / (M_PI * hbar))) * exp(-1/2. * mass_H * omega / hbar * pow(x_space[i] + 500*params.a, 2));
    }

    parameters *p = &params;
    gen_pot_harmonic(p, omega);

    startwf_file = fopen("hydr_input.txt","w");
    for(long int i=0; i<params.L; i++){
        fprintf(startwf_file, "%e\n", creal(start_wf[i]));
    }
    fclose(startwf_file);

    norm_file = fopen("hydr_norm_output.txt","w");
    output_file = fopen("hydr_output.txt","w");
    double norm_obs;

    for(long int t=0; (t * params.tauhat) < simulation_duration; t++){
        if(t % (long)(simulation_duration/(number_time_steps*params.tauhat)) == 0){
            norm_obs = cabs(obs_norm(start_wf, params));
            printf("%li\t%e\n", t, norm_obs - 1);
            fprintf(norm_file, "%li\t%e\n", t, norm_obs - 1);
            for(long int i=0; i < params.L; i++){
                fprintf(output_file,"%e\t",cabs(start_wf[i]));
            }
            fprintf(output_file,"\n");
        }
        step_cn(out_wf, start_wf, params);
        assign_vec(start_wf, out_wf, params.L);
    }

    fclose(norm_file);
    fclose(output_file);

    potential_file = fopen("hydr_potential.txt","w");

    fprintf(potential_file, "a=\t%e\teps=\t%e\n", params.a, params.epsilon);
    for(long int i=0; i<params.L; i++){
        fprintf(potential_file, "%e\n", creal(params.pot[i]));
    }

    printf("Simulation finished!\nPlease find the generated data in the output files\nand use the IPython script to generate an animation");

    fclose(potential_file);

    free(coords);
    free(start_wf);
    free(x_space);
    free(params.pot);
    //free(out_wf);

    return 0;
}
