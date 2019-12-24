#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

#include "structs.h"
#include "geometry.h"
#include "laplacian.h"
#include "vmath.h"



/*calculates the kinetic part of the hamiltonian given the parameter m*/
void kinetic(double complex *out, double complex *in, parameters params){
  /*first the laplacian is calculated*/
  double complex *phi_laplacian= malloc(params.L*sizeof(double complex));
  laplacian(phi_laplacian, in, params.N, params.D);

  /*calculate the kinetic term for in */
  scalar_vec(out, phi_laplacian, (double complex)(-1/(2*params.mhat)), params.L); // multiply by the factor -1/(2*mhat)

  free(phi_laplacian);
}

void kinetic_exp(double complex *out, double complex *in, parameters params) {
  /* rewritten hamiltonian with structs by Tim, untested*/
  laplacian(out, in, params.N, params.D); // calculate the D-dimensional laplacian of in and store it in out
  scalar_vec(out, out, (double complex)(-1/(2*params.mhat)), params.L); // multiply by the factor -1/(2*mhat)
}

void gen_pot_harmonic(parameters *params, double omega){
  for (long int i=0; i<params->L; i++) { /*(*params.L and params->L are synonyms)*/
    long int *coordinates = malloc(params->D* sizeof(long int));
    index2coord(coordinates, i, params->N, params->D);
    double potential = 0;
    for (int j=0; j<params->D; j++) {
      potential += pow(coordinates[j]*params->a, 2);
    }
    potential=0.5*pow(omega, 2)*potential;

    params->pot[i] = (double complex) potential;

    free(coordinates);
  }
  scalar_vec(params->pot, params->pot, 1/params->epsilon, params->L);
}

void gen_pot_box(parameters *params, double height){
    for (long int i=0; i<params->L; i++) {
      long int *coordinates = malloc(params->D* sizeof(long int));
      index2coord(coordinates, i, params->N, params->D);
      double boxy=0;
      for (int j=0; j<params->D; j++) {
          if(fabs(coordinates[j])<params->N/4){
              boxy=0;
          }
          else{
              boxy=height;
              break;
          }
      }

      params->pot[i]=(double complex)boxy;

      free(coordinates);
    }
    scalar_vec(params->pot, params->pot, 1/params->epsilon, params->L);
}

void gen_pot_well(parameters *params, double height){
  for (long int i=0; i<params->L; i++) {
    long int *coordinates = malloc(params->D* sizeof(long int));
    index2coord(coordinates, i, params->N, params->D);
    double boxy=0;
    for (int j=0; j<params->D; j++) {
        if(coordinates[j]>params->N/4+params->N * (int) (coordinates[j]/params->N) && coordinates[j]< 3/4*params->N+params->N * (int) (coordinates[j]/params->N)){
            boxy=0;
        }
        else{
            boxy=height;
            break;
        }
    }

    params->pot[i];

    free(coordinates);
    }
    scalar_vec(params->pot, params->pot, 1/params->epsilon, params->L);
}

void hamiltonian(double complex *out, double complex *in, parameters params){
  /*calculate T*phi and store it in phi_kinetic*/
  double complex *phi_kinetic= malloc(params.L*sizeof(double complex));
  kinetic(phi_kinetic,in,params);

  /*Calculate V*phi and store it in out*/
  mul_element(out, params.pot, in, params.L);

  /*add both parts of the hamiltonian and save it in out */
  add_vec(out, out, phi_kinetic, params.L);

  free(phi_kinetic);
}

void hamiltonian_exp(double complex *out, double complex *in, parameters params){
  /* rewritten hamiltonian with structs by Tim, untested*/
  /*Calculate the kinetic part*/
  double complex *phi_kinetic= malloc(params.L*sizeof(double complex));

  kinetic_exp(phi_kinetic, in, params); // calculate T*phi and store it in phi_kinetic
  mul_element(out, params.pot, in, params.L); // calculate V*phi and store it in out
  add_vec(out, out, phi_kinetic, params.L);

  free(phi_kinetic);
}
