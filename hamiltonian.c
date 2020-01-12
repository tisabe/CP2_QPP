#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <omp.h>

#include "structs.h"
#include "geometry.h"
#include "laplacian.h"
#include "vmath.h"

void kinetic(double complex *out, double complex *in, parameters params) {
  laplacian(out, in, params.N, params.D); // calculate the D-dimensional laplacian of in and store it in out
  mul_compl_vec_real_scal(out, out, -1 / (2*params.mhat), params.L); // multiply by the factor -1/(2*mhat)
}

static void kinetic_slow(double complex *out, double complex *in, parameters params) {
  laplacian_slow(out, in, params.N, params.D); // calculate the D-dimensional laplacian of in and store it in out
  mul_compl_vec_real_scal(out, out, -1 / (2*params.mhat), params.L); // multiply by the factor -1/(2*mhat)
}

/* here you chose one of the three potentials with params.pot, the remaining parameters are also declared in parameters*/
//calculate the harmonic potential
void gen_pot_harmonic(parameters *params){
  double potential;
  long int *coordinates = malloc(params->D* sizeof(long int));
//loop over every index and transform them into the coordinates
  for (long int i=0; i<params->L; i++) { /*(*params.L and params->L are synonyms)*/
    index2coord(coordinates, i, params->N, params->D);
    potential = 0.0;
    for (int j=0; j<params->D; j++) { //loop over the coordinate array to calculate its absolute value with offset N/2 to center the potential
      potential += pow(coordinates[j]-params->N/2, 2);
    }
    potential=0.5*(params->khat)*potential;

    params->pot[i] = (double complex) potential;
  }

  free(coordinates);
}
//calculate the box potential
void gen_pot_box(parameters *params, double height){
    double boxy = 0;
    for (long int i=0; i<params->L; i++) {
      long int *coordinates = malloc(params->D* sizeof(long int));
      index2coord(coordinates, i, params->N, params->D);

      //the following loop checks if one the coordinates for index i is greater than N/4, if this is true -> set potential at this point to height
      for (int j=0; j<params->D; j++) {
          if(fabs(coordinates[j]-params->N/2)<params->N/4){
              boxy=0;
          }
          else{
              boxy=height;
              break;
          }
      }

      params->pot[i]=(double complex) boxy;

      free(coordinates);
    }
    scalar_vec(params->pot, params->pot, 1/params->epsilon, params->L);
}

//calculate the potential well, similar to the boxpotenial with one dimension fewer
void gen_pot_well(parameters *params, double height){
  double boxy = 0;
  //this loop checks if the index is in the gap between N/4 and 3N/4 of its row, outside this gap the potenial is set to "height"
  for (long int i=0; i<params->L; i++) {
        if(i>params->N/4+params->N * (int) (i/params->N) && i< 3/4*params->N+params->N * (int) (i/params->N)){
            boxy=0;
        }
        else{
            boxy=height;
            break;
        }

		params->pot[i]=(double complex) boxy;

	}
  scalar_vec(params->pot, params->pot, 1/params->epsilon, params->L);
}

void hamiltonian(double complex *out, double complex *in, parameters params){
  /*calculate T*phi and store it in phi_kinetic*/
  static double complex *phi_kinetic = NULL;
  static long int lPrev;
  if((phi_kinetic == NULL) || (lPrev != params.L)){
    free(phi_kinetic);
    phi_kinetic =  malloc(params.L*sizeof(double complex));
    lPrev = params.L;
  }

  kinetic(phi_kinetic,in,params);

  /*Calculate V*phi and store it in out*/
  mul_element(out, params.pot, in, params.L);

  /*add both parts of the hamiltonian and save it in out */
  add_vec(out, out, phi_kinetic, params.L);
}

void hamiltonian_slow(double complex *out, double complex *in, parameters params){
  /*calculate T*phi and store it in phi_kinetic*/
  static double complex *phi_kinetic = NULL;
  static long int lPrev;
  if((phi_kinetic == NULL) || (lPrev != params.L)){
    free(phi_kinetic);
    phi_kinetic =  malloc(params.L*sizeof(double complex));
    lPrev = params.L;
  }

  kinetic_slow(phi_kinetic,in,params);

  /*Calculate V*phi and store it in out*/
  mul_element(out, params.pot, in, params.L);

  /*add both parts of the hamiltonian and save it in out */
  add_vec(out, out, phi_kinetic, params.L);
}

void hamiltonian_single_loop(double complex *out, double complex *in, parameters params) {
  // first part comes from laplacian
  static long int *neighbours = NULL;
  //static long int prev_N, prev_D, L;
  static long int prev_N, prev_D;

  //Check if geometry parameters have changed. If yes, free, reallocate and recalculate the next neighbours the array again
  if((neighbours == NULL) || (prev_N != params.N) || (prev_D != params.D)){
      free(neighbours);
      //L = ipow(params.N, params.D);
      neighbours = malloc(2*params.D*params.L*sizeof(long int));
      nneighbour_init(neighbours, params.N, params.D);
      prev_N = params.N;
      prev_D = params.D;
  }

  //Loop over all points in the lattice
  for(long int i=0; i<params.L; i++){
      out[i] = params.pot[i]*in[i]; //first initialize element with potential part
      // Loop over all dimensions
      for(int d=0; d<params.D; d++){
          // Calculate discrete second derivative according to d^2 psi/ dx^2 = (f(x+h)-2*f(x)+f(x-h))/h^2
          // and multiply with factor from kinetic energy
          out[i] += -1 / (2*params.mhat)*(in[neighbours[2*(i*params.D+d)]]-2*in[i]+in[neighbours[2*(i*params.D+d)+1]]);
      }
  }
}

void hamiltonian_parallel(double complex *out, double complex *in, parameters params) {
  /* This implementation of the hamiltonian should be the fastest.
  Parallel speedup depends on array length and number of available computing cores.*/
  // first part comes from laplacian
  static long int *neighbours = NULL;
  //static long int prev_N, prev_D, L;
  static long int prev_N, prev_D;

  //Check if geometry parameters have changed. If yes, free, reallocate and recalculate the next neighbours the array again
  if((neighbours == NULL) || (prev_N != params.N) || (prev_D != params.D)){
      free(neighbours);
      //L = ipow(params.N, params.D);
      neighbours = malloc(2*params.D*params.L*sizeof(long int));
      nneighbour_init(neighbours, params.N, params.D);
      prev_N = params.N;
      prev_D = params.D;
  }
  long int i;

  // the next line is a directive for the compiler to spread the loop over all available cores
  #pragma omp parallel for
  for(i=0; i<params.L; i++){
      out[i] = params.pot[i]*in[i]; //first initialize element with potential part
      // Loop over all dimensions
      for(int d=0; d<params.D; d++){
          // Calculate discrete second derivative according to d^2 psi/ dx^2 = (f(x+h)-2*f(x)+f(x-h))/h^2
          // and multiply with factor from kinetic energy
          out[i] += -1 / (2*params.mhat)*(in[neighbours[2*(i*params.D+d)]]-2*in[i]+in[neighbours[2*(i*params.D+d)+1]]);
      }
  }
}
