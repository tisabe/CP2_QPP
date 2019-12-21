#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <float.h>

#include "structs.h"
#include "euler_method.h"
#include "vmath.h"
#include "geometry.h"

#include "hermite_polynomial.h"

#ifndef USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#ifndef LEN
#define LEN 256
#endif // LEN

int main(){
    parameters params;
    params.max_iter = 100;
    params.tol = DBL_EPSILON;
    params.mhat = 1;
    params.epsilon = 27.4;

    int potential_type, poly_type, poly_order, integrator;
    double omega;
    FILE *file, *file_start_wf;

    printf("********  Solver for the Schroedinger equation with three different integrators  ********\n");
    /*printf("Please set the number of dimensions D = ");
    scanf("%ui",&params.D);
    printf("Please set the number of points in each dimension N = ");
    scanf("%li",&params.N);
    params.L = ipow(params.N,params.D);
    printf("That makes %li points in the lattice!\n",params.L);
    printf("Please select the integrator that you would like to use:\n\t1. Euler integrator (does not conserve probability)\n\t2. Crank Nicolson method\n\t3. Strang splitting method (work in progress)\n");
    printf("Integrator: ");
    scanf("%i",&integrator);
    printf("Please set a time step (in 2.4 * 10^-17 s): ");
    scanf("%lf",&params.tau);
    printf("Please set the duration that the simulation should run for: ");
    scanf("%lf",&params.total_time);
    printf("Please choose a potential:\n\t1. Harmonic potential in %u dimension(s)\n\t2. Box potential in %ui dimension(s) \n",params.D,params.D);
    printf("Potential: ");
    scanf("%i",&potential_type);
    printf("Please choose a parameter: ");
    scanf("%lf",&omega);
    params.parameter = omega;
    printf("Now please choose a wave function that you want to start with:\n\t1. Solutions to 1 dim Quantum HO\n");
    printf("Start wave function: ");
    scanf("%i",&poly_type);
    */

    params.D = 1;
    params.N = 1501;
    params.L = ipow(params.N, params.D);
    integrator = 1;
    params.tau = 0.05;
    params.total_time = 100;
    potential_type = 1;
    omega = 1;
    poly_type = 1;
    poly_order = 1;

    double complex *start_wf = malloc(params.L*sizeof(double complex));
    set_zero(start_wf, params.L);
    long int *coords = malloc(params.D*sizeof(long int));

    if(poly_type == 1){
        /*printf("Which order would you like? (1, 2 or 3) ");
        scanf("%i",poly_order);*/
        if(params.D == 1){
            double norm_const = 1/sqrt(ipow(2,(unsigned int) poly_order)*fact((int) poly_order))*sqrt(sqrt(omega/M_PI));
            double *hermite = malloc((poly_order+1)*sizeof(double));
            file_start_wf = fopen("main_start_wf.txt","w");

            for(long int i=0; i<params.L; i++){
                index2coord(coords, i, params.N, params.D);
                double *point = malloc(1*sizeof(double));
                point[0] = sqrt(omega)*coords[0];
                hermite = h_polynomial_value(1,poly_order,point);
                start_wf[i] = norm_const*exp(omega*pow(coords[0],2.)/2.)*hermite[poly_order];
                free(point);
                fprintf(file_start_wf, "%lf, %lf\n", creal(start_wf[i]), cimag(start_wf[i]));
            }
            fclose(file_start_wf);

            free(hermite);
        }else{
            printf("Invalid number of dimensions! Ongoing work at this point!\n");
        }
    }else{
        printf("That is not a valid input!\n");
    }
    double complex *output = malloc(params.L*sizeof(double complex));
    if(integrator == 1){
        euler_method(output, start_wf, params);

    }else if(integrator == 2){
        //Call Crank Nicolson
        printf("This integrator has not been implemented yet. Please choose another one.\n");
    }else if(integrator == 3){
        //Call Strang splitting
        printf("This integrator has not been implemented yet. Please choose another one.\n");
    }else{
        printf("That is not a valid input!\n");
    }

    file = fopen("main_output.txt","w");
    for(long int i=0; i<params.L; i++){
        fprintf(file, "%li, %lf, %lf\n", i, creal(output[i]), cimag(output[i]));
    }
    fclose(file);

    free(output);
    free(start_wf);
    free(coords);

    return 0;
}
