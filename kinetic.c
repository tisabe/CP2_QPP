#include <stdlib.h>
#include <stdio.h>
#include "laplacian.h"

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
