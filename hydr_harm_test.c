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

    params.N = 1000;
    params.D = 1;
    params.L = ipow(params.N, params.D);
    params.tol = DBL_EPSILON;
    params.max_iter = params.L;

    params.tauhat = 1e-3;
    params.epsilon = hbar * omega;
    params.a = 2e-13;
    params.mhat = pow(params.a/hbar,2)*mass_H*params.epsilon;     //=mass_H*omega*a^2/(2*hbar)
    params.khat = params.mhat * pow(hbar*omega/params.epsilon,2); //=params.mhat

    double simulation_duration = 3e1;
    int number_time_steps = 100;
    double offset[3] = {50,50,50};

    double complex *out_wf = malloc(params.L * sizeof(double complex));
    double complex *start_wf = malloc(params.L * sizeof(double complex));
    long int *coords = malloc(params.D * sizeof(long int));
    params.pot = malloc(params.L * sizeof(double complex));

    FILE *potential_file;
    FILE *obs_file;
    FILE *output_file_real;
    FILE *output_file_imag;

    set_zero(start_wf, params.L);

    for(long int i=0; i<(params.L); i++){
        index2coord(coords, i, params.N, params.D);
        for(int j=0; j<params.D; j++){
            start_wf[i] += pow(sqrt(params.a),params.D) * sqrt(sqrt(params.mhat/M_PI)/params.a) * exp(-.5*params.mhat*pow((coords[j]-params.N/2+100),2));
        }
    }

    parameters *p = &params;
    gen_pot_harmonic(p, omega);

    obs_file = fopen("hydr_observables_output.txt","w");
    output_file_real = fopen("hydr_wavefunction_output_real.txt","w");
    output_file_imag = fopen("hydr_wavefunction_output_imag.txt","w");
    double norm_obs;

    printf("\nElapsed time\tNorm - 1\t<H>\t\t<x>\t\t<p>\n\n");

    for(long int t=0; (t * params.tauhat) < simulation_duration; t++){
        if(t % (long)(simulation_duration/(number_time_steps*params.tauhat)) == 0){
            norm_obs = cabs(obs_norm(start_wf, params));
            printf("%e\t%e\t%e\t%e\t%e\n", t*params.tauhat, norm_obs - 1, creal(obs_E(start_wf,params)), creal(obs_x(start_wf, 0, params)), creal(obs_p(start_wf, 0, params)));
            fprintf(obs_file, "%e\t%e\t%e\t%e\t%e\n", t*params.tauhat, norm_obs - 1, creal(obs_E(start_wf,params)), creal(obs_x(start_wf, 0, params)), creal(obs_p(start_wf, 0, params)));
            for(long int i=0; i < params.L; i++){
                fprintf(output_file_real,"%e\t",creal(start_wf[i]));
                fprintf(output_file_imag,"%e\t",cimag(start_wf[i]));
            }
            fprintf(output_file_real,"\n");
            fprintf(output_file_imag,"\n");
        }
        step_strang(out_wf, start_wf, params);
        assign_vec(start_wf, out_wf, params.L);
    }

    fclose(obs_file);
    fclose(output_file_real);
    fclose(output_file_imag);

    potential_file = fopen("hydr_potential.txt","w");

    fprintf(potential_file, "a=\t%e\teps=\t%e\n", params.a, params.epsilon);
    for(long int i=0; i<params.L; i++){
        fprintf(potential_file, "%e\n", creal(params.pot[i]));
    }

    printf("Simulation finished!\nPlease find the generated data in the output files\nand use the IPython script to generate an animation");

    fclose(potential_file);

    free(coords);
    free(start_wf);
    free(params.pot);
    //free(out_wf);

    return 0;
}
