#include <complex.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "vmath.h"
#include "geometry.h"
#include "hamiltonian.h"
#include "nfft.h"
#include "laplacian.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

/*In this file functions to calculate different observables from a wavefunction
are defined. not finished, not tested*/

double complex obs_norm(double complex *in, parameters params) {
  /* Computes the normalization of a wavefunction. */
  double complex var_norm = 0.0;

  var_norm = dot_product(in, in, params.L);

  return var_norm;
}

double complex obs_E(double complex *in, parameters params) {
  /* Computes the energy expectation value for a given wavefunction and hamiltonian. */
  double complex var_E = 0.0;
  double complex *in_H = malloc(params.L*sizeof(double complex));
  hamiltonian(in_H, in, params);
  var_E = dot_product(in, in_H, params.L)/obs_norm(in, params);

  free(in_H);
  return var_E;
}

double complex obs_x(double complex *in, unsigned int d, parameters params) {
  /* Computes the average position of a wavefunction in dimension d given the
  number of dimensions D and length of axis N using the struct params.
  According to eq. 84.*/
  double complex var_x = 0.0;
  long int *coord = malloc(params.D*sizeof(double complex));
  for(long int i=0; i<params.L; i++) {
    index2coord(coord, i, params.N, params.D);
    var_x += (~(in[i]))*coord[d]*in[i]; // calculated relative to origin
  }

  free(coord);
  return var_x/obs_norm(in, params);
}

double complex obs_delta_x(double complex *in, parameters params) {
  /* Computes the position width of a wavefunction given the number of
  dimensions D and length of axis N using the struct params.
  According to eq. 85.*/
  double complex var_x = 0.0;
  long int *coord = malloc(params.D*sizeof(double complex));
  double abs_coord = 0.0; // squared absolute of coordinate vector
  for(long int i=0; i<params.L; i++) {
    index2coord(coord, i, params.N, params.D);
    for(unsigned int d=0; d<params.D; d++){
      abs_coord += pow((double)coord[d], 2);
    }
    var_x += (~(in[i]))*abs_coord*in[i]; // calculated relative to origin
  }
  var_x *= 1/obs_norm(in, params);
  for(unsigned int d=0; d<params.D; d++) {
    var_x -= pow(obs_x(in, d, params), 2);
  }
  var_x = sqrt(var_x);

  free(coord);
  return var_x;
}



double complex obs_p(double complex *in, unsigned int d, parameters params) {
  /* Computes the average momentum of a wavefunction in dimension d given the
  number of dimensions D and length of axis N using the struct params.
  According to eq. 89.*/
  double complex norm= obs_norm(in, params);
  double complex var_p = 0.0;
  double complex *psi_k = malloc(params.L*sizeof(double complex));
  long int *coordinate = malloc(params.D* sizeof(long int));
  nfft(psi_k, in, params.N, params.D);

  for(long int i=0; i<params.L; i++) {
    index2coord(coordinate, i, params.N, params.D);
    var_p += (~(psi_k[i])) * cexp(1I * M_PI/params.N * coordinate[d]) * 2 * sin(M_PI/params.N*coordinate[d]) * psi_k[i];
  }

  free(psi_k);
  return var_p/norm;
}

double complex obs_delta_p(double complex *in, parameters params) {
  /* Computes the momentum width of a wavefunction given the number of
  dimensions D and length of axis N using the struct params.
  According to eq. 90. */
  double complex norm= obs_norm(in, params);
  double complex var_p = 0.0;
  double complex *psi_lap = malloc(params.L*sizeof(double complex));
  long int *coordinate = malloc(params.D* sizeof(long int));

  laplacian(psi_lap, in, params.N, params.D);

  var_p = dot_product(in, psi_lap, params.L);
  var_p *= -1/obs_norm(in, params);
  for(unsigned int d=0; d<params.D; d++) {
    var_p -= pow(obs_p(in, d, params), 2);
  }
  var_p = sqrt(var_p);

  free(psi_lap);
  free(coordinate);
  return var_p/norm;
}
