#include <complex.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "vmath.h"
#include "indices.h"
#include "hamiltonian.h"
#include "fft.h" // hypothetical

/*In this file functions to calculate different observables from a wavefunction
are defined. not finished, not tested*/

double complex obs_norm(double complex *in, parameters params) {
  double complex var_norm = 0.0;
  var_norm = dotproduct(in, in, params.L);

  free(var_norm);
  return var_norm;
}

double complex obs_E(double complex *in, parameters params) {
  double complex var_E = 0.0;
  double complex *in_H = malloc(params.L*sizeof(double complex));
  hamiltonian(in_H, in, params);
  var_norm = dot_product(in, in_H, params.L)/obs_norm(in, params);

  free(in_H);
  return var_norm/params.N;
}

double complex obs_x(double complex *in, unsigned int d, parameters params) {
  double complex var_x = 0.0;
  long int *coord = malloc(params.D*sizeof(double complex));
  double complex *in_x = malloc(params.L*sizeof(double complex));
  for(long int i=0; i<params.L; i++) {
    index2coord(coord, i, params.N, params.D);
    in_x[i] = coord[d]*in[i];
  }
  var_x = dot_product(in, in_x, params.L)/obs_norm(in, params);

  free(coord);
  free(in_x);
  return var_x/params.N;
}

double complex obs_p(double complex *in, unsigned int d, parameters params) {
  double complex var_p = 0.0;
  double complex *psi_k = malloc(params.L*sizeof(double complex));
  fft(psi_k, in, params.N, params.D);
  for(long int i=0; i<params.L; i++) {
    var_p += 0;
  }

  free(psi_k);
  return var_p/params.N;
}
