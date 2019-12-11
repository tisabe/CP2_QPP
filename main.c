#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "euler_method.h"
#include "vmath.h"
#include "indices.h"

#include "hermite_polynomial.h"

#define _USE_MATH_DEFINES

int main(){
    unsigned int D;
    long int N, L, integrator, poly_type, poly_order;
    int potential_type;
    double tau, total_time, omega;

    printf("********  Solver for the Schroedinger equation with three different integrators  ********\n");
    printf("Please set the number of dimensions D = ");
    scanf("%ui",&D);
    printf("Please set the number of points in each dimension N = ");
    scanf("%li",&N);
    L = ipow(N,D);
    printf("That makes %li points in the lattice!\n",L);
    printf("Please select the integrator that you would like to use:\n\t1. Euler integrator (does not conserve probability)\n\t2. Crank Nicolson method\n\t3. Strang splitting method (work in progress)\n");
    printf("Integrator: ");
    scanf("%i",&integrator);
    printf("Please set a time step (in 2.4 * 10^-17 s)");
    scanf("%lf",&tau);
    printf("Please set the duration that the simulation should run for:");
    scanf("%lf",&total_time);
    printf("Please choose a potential:\n\t1. Harmonic potential in %ui dimension(s)\n\t2. Box potential in %ui dimension(s)",D,D);
    printf("Potential: ");
    scanf("%i",&potential_type);
    printf("Please choose a parameter: ");
    scanf("%lf",&omega);
    printf("Now please choose a wave function that you want to start with:\n\t1. Solutions to 1 dim Quantum HO\n");
    scanf("%i",&poly_type);

    double complex *start_wf = malloc(L*sizeof(double complex));
    set_zero(start_wf, L);
    long int *coords = malloc(D*sizeof(long int));

    if(poly_type == 1){
        printf("Which order would you like? (1, 2 or 3) ");
        scanf("%i",poly_order);
        if(D == 1){
            double norm_const = 1/sqrt(ipow(2,(unsigned int) poly_order)*fact((int) poly_order))*sqrt(sqrt(omega/M_PI));
            double *hermite = malloc((poly_order+1)*sizeof(double));

            for(long int i=0; i<L; i++){
                index2coord(coords, i, N, D);
                double *point = malloc(1*sizeof(double));
                point[0] = sqrt(omega)*coords[0];
                hermite = h_polynomial_value(1,poly_order,point);
                start_wf[i] = norm_const*exp(omega*pow(coords[0],2.)/2.)*hermite[poly_order];
                free(point);
            }

            free(hermite);
        }else{
            printf("Invalid number of dimensions! Ongoing work at this point!\n");
        }
    }else{
        printf("That is not a valid input!\n");
    }
    if(integrator == 1){
        double complex *euler_out = malloc(L*sizeof(double complex));
        euler_method(euler_out, start_wf, N, D, L, tau, total_time, 1.0, 27.4, potential_type, omega);
    }else if(integrator == 2){
        //Call Crank Nicolson
        printf("This integrator has not been implemented yet. Please choose another one.\n");
    }else if(integrator == 3){
        //Call Strang splitting
        printf("This integrator has not been implemented yet. Please choose another one.\n");
    }else{
        printf("That is not a valid input!\n");
    }

    free(start_wf);
    free(coords);
}
