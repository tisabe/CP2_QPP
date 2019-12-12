#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#include "kinetic.h"
#include "potentials.h"

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



  /*Calculte the Hamiltonian of in */
  for (int i=0; i<L; i++) {
        /*Calculates the hamiltonian of wwave vector in at index i  */
    out[i] = phi_kinetic[i]+1/epsilon*phi_potential[i]*in[i];
  }

  free(phi_kinetic);
  free(phi_potential);
}
