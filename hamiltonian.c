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
  for (int i=0; i<params.L; i++) {
    out[i] = -1/(2*params.mhat)*phi_laplacian[i];
  }

  free(phi_laplacian);
}

void kinetic_exp(double complex *out, double complex *in, parameters params) {
  /* rewritten hamiltonian with structs by Tim, untested*/
  laplacian(out, in, params.N, params.D); // calculate the D-dimensional laplacian of in and store it in out
  scalar_vec(out, out, (double complex)(-1/(2*params.mhat)), params.L); // multiply by the factor -1/(2*mhat)
}

double harmonic0(long int index, parameters params) {
    long int *coordinates = malloc(params.D* sizeof(long int));
    index2coord(coordinates,index, params.N, params.D);
    double potential = 0;
    for (int i=0; i<params.D; i++) {
        potential += pow(coordinates[i], 2);
  }
    potential=0.5*pow(params.parameter, 2)*potential;
    free(coordinates);
    return potential;
}

void harmonic(double *harmonicpotential, parameters params){
    for (int i=0; i<params.L; i++) {
        harmonicpotential[i]=harmonic0( i, params);
    }
}


double boxpotential0(long int index, parameters params) {
    long int *coordinates = malloc(params.D* sizeof(long int));
    index2coord(coordinates,index, params.N, params.D);
    long int boxy=0;
    for (int i=0; i<params.D; i++) {
        if(fabs(coordinates[i])<params.N/4){
            boxy=boxy;
        }else{
            boxy=params.parameter;
            break;
            }
    free(coordinates);
    return boxy;
  }
}

void box(double *boxpotential, parameters params){
    for (int i=0; i<params.L; i++) {
        boxpotential[i]=boxpotential0( i, params);
    }
}

long int potentialwell(parameters params) {
    long int *boxy= malloc(params.L * sizeof(long int));
    for (int i=0; i<params.D; i++) {
        if(i>params.N/4+params.N * (int) (i/params.N) && i< 3/4*params.N+params.N * (int) (i/params.N)){
            boxy[i]=0;
        }else{
            boxy[i]=params.parameter;
            }
    }
  return *boxy;
}

void well(double *wellpotential, parameters params){
  for (int i=0; i<params.L; i++) {
      wellpotential[i]=potentialwell(params);
  }
}

void hamiltonian(double complex *out, double complex *in, parameters params){
  /*Calculate the kinetic part*/
  double complex *phi_kinetic= malloc(params.L*sizeof(double complex));
  kinetic(phi_kinetic,in,params);

  /*Calculate the harmonic part*/
  double complex *phi_potential= malloc(params.L*sizeof(double complex));

  if (params.pot==NULL){
    free(params.pot);
    params.pot = malloc(params.L * sizeof(double complex));
    double *temp = malloc(params.L * sizeof(double));
    if (params.ext_potential_type==0) {
	     harmonic(temp, params);
	  }

    else if (params.ext_potential_type==1) {
	     box(temp, params);
	  }


    else if (params.ext_potential_type==2) {
	     well(temp, params);
	  }

    for(long int i=0; i<params.L; i++){
        params.pot[i] = (double complex) temp[i];
    }
    free(temp);
    scalar_vec(params.pot, params.pot, 1/params.epsilon, params.L);     //Error here as scalar_vec expects double complex values but params.pot is double
  }

  mul_element(phi_potential, params.pot, in, params.L);
  add_vec(out, phi_potential, phi_kinetic, params.L);

  free(phi_kinetic);
  free(phi_potential);
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
