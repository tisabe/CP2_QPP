#include <complex.h>
#include <stdlib.h>

#include "structs.h"
#include "vmath.h"
#include "hamiltonian.h"

static void f_eta(double complex *out_eta, double complex *in_eta, parameters params) {
  double complex *temp = malloc(params.L * sizeof(double complex));
  hamiltonian(temp, in_eta, params);
  hamiltonian(out_eta, temp, params);
  scalar_vec(out_eta, out_eta, (double complex) params.tauhat*params.tauhat/4, params.L);
  add_vec(out_eta, out_eta, in_eta, params.L);
  free(temp);
}

static void f_psi(double complex *out_psi, double complex *in_psi, parameters params) {
  double complex *temp = malloc(params.L * sizeof(double complex));
  hamiltonian(temp, in_psi, params);
  hamiltonian(out_psi, temp, params);
  scalar_vec(out_psi, out_psi, (-1)*1I*params.tauhat/2, params.L);
  add_vec(out_psi, out_psi, in_psi, params.L);
  free(temp);
}

void step_cn(double complex *out, double complex *in, parameters params) {
  double complex *temp = malloc(params.L * sizeof(double complex));
  cg(temp, f_eta, in, params);
  f_psi(out, temp, params);
  free(temp);
}
