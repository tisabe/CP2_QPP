#include <stdlib.h>
#include <stdio.h>
#include "kinetic.h"
#include "potentials.h"

void hamiltonian(double complex *out, double complex *in, long int N, unsigned int D, long int L, double m, double epsilon, int ext_potential_type){
  /*Calculate the kinetic part*/
  double complex *phi_kinetic= malloc(L*sizeof(double complex));
  kinetic(phi_kinetic,in,N,D,L,m);

  /*Calculate the harmonic part*/
  double complex *phi_potential= malloc(L*sizeof(double complex));
  potential(phi_potential,in,N,D,L,epsilon,ext_potential_type);

  /*Calculte the Hamiltonian of in */
  for(int i=0; i<L; i++){
        /*Calculates the hamiltonian of wwave vector in at index i  */
    out[i] = phi_kinetic[i]+phi_harmonic[i];
  }

  free(phi_kinetic);
  free(phi_potential);
}
