#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

#include "indices.h"
#include "laplacian.h"
#include "structs.h"
#include "vmath.h"



/*calculates the kinetic part of the hamiltonian given the parameter m*/
void kinetic(double complex *out, double complex *in, long int N, unsigned int D, long int L, double m){
  /*first the laplacian is calculated*/
  double complex *phi_laplacian= malloc(L*sizeof(double complex));
  laplacian(phi_laplacian, in, N, D);

  /*calculate the kinetic term for in */
  for (int i=0; i<L; i++) {
    out[i] = -1/(2*m)*phi_laplacian[i];
  }

  free(phi_laplacian);
}

void kinetic_exp(double complex *out, double complex *in, parameters params) {
  /* rewritten hamiltonian with structs by Tim, untested*/
  laplacian(out, in, params.N, params.D); // calculate the D-dimensional laplacian of in and store it in out
  scalar_vec(out, out, (double complex)(-1/(2*params.mhat)), params.L); // multiply by the factor -1/(2*mhat)
}

double harmonic0(long int index, double omega, long int N, unsigned int D) {
    long int * coordinates = malloc(D* sizeof(long int));
    index2coord(coordinates,index, N, D);
    double potential = 0;
    for (int i=0; i<D; i++) {
        potential += pow(coordinates[i], 2);
  }
    potential=0.5*pow(omega, 2)*potential;
    free(coordinates);
    return potential;
}
void harmonic(double *harmonicpotential, double omega, long int N, unsigned int D){
    for (int i=0; i<ipow(N, D); i++) {
        harmonicpotential[i]=harmonic0( i, omega,N, D);
    }
}



double boxpotential0(long int index, int height, long int N, unsigned int D) {
    long int * coordinates = malloc(D* sizeof(long int));
    index2coord(coordinates,index, N, D);
    long int boxy=0;
    for (int i=0; i<D; i++) {
        if(fabs(coordinates[i])<N/4){
            boxy=boxy;
        }else{
            boxy=height;
            break;
            }
    free(coordinates);
    return boxy;
  }
}
void box(double *boxpotential, int height, long int N, unsigned int D){
    for (int i=0; i<ipow(N, D); i++) {
        boxpotential[i]=boxpotential0( i, height, N, D);
    }
}



long int potentialwell(int height, long int N, unsigned int D) {
    long int * boxy= malloc(ipow(N, D) * sizeof(long int));
    for (int i=0; i<D; i++) {
        if(i>N/4+N * (int) (i/N) && i< 3/4*N+N * (int) (i/N)){
            boxy[i]=0;
        }else{
            boxy[i]=height;
            }
    }
  return *boxy;
}
void well(double *wellpotential, int height, long int N, unsigned int D){
  for (int i=0; i<ipow(N, D); i++) {
      wellpotential[i]=potentialwell(height, N, D);
  }
}

void hamiltonian(double complex *out, double complex *in, long int N, unsigned int D, long int L, double m, double epsilon, int ext_potential_type, double parameter){
  /*Calculate the kinetic part*/
  double complex *phi_kinetic= malloc(L*sizeof(double complex));
  kinetic(phi_kinetic,in,N,D,L,m);

  /*Calculate the harmonic part*/
  double *phi_potential= malloc(L*sizeof(double));

  if (ext_potential_type==0) {
	  harmonic(phi_potential,parameter,N,D);
	  }

  else if (ext_potential_type==1) {
	  box(phi_potential,parameter,N,D);
	  }


  else if (ext_potential_type==2) {
	  well(phi_potential,parameter,N,D);
	  }



  /*Calculate the Hamiltonian of in */
  for (int i=0; i<L; i++) {
        /*Calculates the hamiltonian of wwave vector in at index i  */
    out[i] = phi_kinetic[i]+1/epsilon*phi_potential[i]*in[i];
  }

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
