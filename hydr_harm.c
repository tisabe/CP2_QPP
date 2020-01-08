/* ********** This program solves the time-dependent Schroedinger equation for a particle with a certain wave function in a harmonic potential ********** */
/* ********** The parameters of the program are set according to the harmonic approximation of two hydrogen atoms in a vibrating mode ********** */


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
#include "phi.h"

#define _USE_MATH_DEFINES

int main(){
    // Set simulation parameters
    parameters params;

    //Set external parameters which the simulation depends upon
    double omega = 8.3e14; //Hz
    double mass_H = 1.67e-27; //kg
    double hbar = 1.054571817e-34; //Js

    //Set parameters for discretising space
    printf("\n Please set the number of points (recommended 1000-5000): N = ");
    scanf("%li",&params.N);
    params.D = 1;
    params.L = ipow(params.N, params.D);

    //Set tolerance/ aim parameters for CG algorithm
    params.tol = DBL_EPSILON;
    params.max_iter = params.L;

    //Set energy scale
    params.epsilon = hbar * omega;

    //Set dimensionless simulation parameters
    printf("\n Please set tauhat (recommended 1e-7): tauhat = ");
    scanf("%lf",&params.tauhat); // = omega * tau (with epsilon = hbar*omega)
    printf("\n Please set the distance between two lattice points a (recommended 2e-13): a = ");
    scanf("%lf", &params.a);
    params.mhat = pow(params.a/hbar,2)*mass_H*params.epsilon;
    params.khat = params.mhat * pow(hbar*omega/params.epsilon,2);

    double simulation_duration;
    int number_timesteps_to_record;
    printf("\n Please set the simulation duration as multiple of omega (recommended 10): total_duration = ");
    scanf("%lf", &simulation_duration); //Set simulation duration as a multiple of omega
    printf("\n Please set the number of time steps to record (recommended 100): time_steps_to_record = ");
    scanf("%i", &number_timesteps_to_record); //Set the number of evenly spaced time steps that are to be recorded

    //Initialise arrays to store the start wave function, the wave function at the current iteration and the potential
    double complex *start_wf = malloc(params.L * sizeof(double complex));
    double complex *out_wf = malloc(params.L * sizeof(double complex));
    params.pot = malloc(params.L * sizeof(double complex));

    //Initialise txt files to store the potential, observables, real part of the wave functions, imaginary part of the wave functions
    FILE *potential_file;
    FILE *obs_file;
    FILE *output_file_real;
    FILE *output_file_imag;

    //Initialise an array to store coordinates
    long int *coords = malloc(params.D * sizeof(long int));

    //Set start wave function to be the ground state of the solution of the quantum harmonic oscillator
    int excitation;
    double offset, scaling;
    printf("\n Please set the parameters for the start wave function.\nPlease choose if you want the solution of the ground state or the first excited state (0 or 1): excitation = ");
    scanf("%i", &excitation);
    printf("\n Please set a scaling factor for the wave function (set 1 for eigenfunctions): scaling = ");
    scanf("%lf", &scaling);
    printf("\n Please set an offset of the wave function from the center (recommended N/10): offset = ");
    scanf("%lf", &offset);
    for(long int i=0; i<(params.L); i++){
        index2coord(coords, i, params.N, params.D);
        if(excitation == 0){
            start_wf[i] = sqrt(sqrt(scaling))*pow(sqrt(params.a),params.D) * sqrt(sqrt(params.mhat/M_PI)/params.a) * exp(-.5*scaling*params.mhat*pow((coords[0]-params.N/2+offset),2));
        }else if(excitation == 1){
            start_wf[i] = sqrt(sqrt(scaling))*2*sqrt(scaling)*sqrt(params.mhat/2)*(coords[0]-params.N/2+offset)*pow(sqrt(params.a),params.D) * sqrt(sqrt(params.mhat/M_PI)/params.a) * exp(-.5*scaling*params.mhat*pow((coords[0]-params.N/2+offset),2));
        }
    }

    //Set the potential to be the harmonic potential with eigenfrequency omega
    parameters *p = &params;
    gen_pot_harmonic(p, omega);

    //Open the txt files to store the observables, real part of the wave functions, imaginary part of the wave functions
    obs_file = fopen("hydr_observables_output.txt","w");
    output_file_real = fopen("hydr_wavefunction_output_real.txt","w");
    output_file_imag = fopen("hydr_wavefunction_output_imag.txt","w");

    printf("\nElapsed time\tNorm - 1\t<E>\t\t<x>\t\t<p>\n\n");

    //Iterate over t and have the for loop running until t*tauhat is greater or equal to simulation_duration
    for(long int t=0; (t * params.tauhat) < simulation_duration; t++){

        //If the condition below holds, the wave function as well as the observables normalisation - 1, energy, x expectation value and p expectation value are recorded in the output files
        if(t % (long)(simulation_duration/(number_timesteps_to_record*params.tauhat)) == 0){
            //Print all the observables to the console...
            printf("%e\t%e\t%e\t%e\t%e\n", t*params.tauhat, cabs(obs_norm(start_wf, params)) - 1, creal(obs_E(start_wf,params)), creal(obs_x(start_wf, 0, params)), creal(obs_p(start_wf, 0, params)));
            //... and save them to the obs_file as well
            fprintf(obs_file, "%e\t%e\t%e\t%e\t%e\n", t*params.tauhat, cabs(obs_norm(start_wf, params)) - 1, creal(obs_E(start_wf,params)), creal(obs_x(start_wf, 0, params)), creal(obs_p(start_wf, 0, params)));
            //Save the wave functions
            for(long int i=0; i < params.L; i++){
                fprintf(output_file_real,"%e\t",creal(start_wf[i]));
                fprintf(output_file_imag,"%e\t",cimag(start_wf[i]));
            }
            fprintf(output_file_real,"\n");
            fprintf(output_file_imag,"\n");
        }

        //Execute one step of the chosen integrator
        step_euler(out_wf, start_wf, params);
        //step_cn(out_wf, start_wf, params);
        //step_strang(out_wf, start_wf, params);

        //Assign the vector to the old out_vector to prepare for the next step
        assign_vec(start_wf, out_wf, params.L);
    }

    //Close output files
    fclose(obs_file);
    fclose(output_file_real);
    fclose(output_file_imag);

    //Save potential into a file
    potential_file = fopen("hydr_potential.txt","w");

    fprintf(potential_file, "a=\t%e\teps=\t%e\ttauhat=\t%e\n", params.a, params.epsilon, params.tauhat);
    for(long int i=0; i<params.L; i++){
        fprintf(potential_file, "%e\n", creal(params.pot[i]));
    }

    fclose(potential_file);

    printf("Simulation finished!\nPlease find the generated data in the output files\nand use the IPython script 'animation-1d' to generate an animation");

    free(coords);
    free(start_wf);
    free(params.pot);
    //free(out_wf);

    return 0;
}
