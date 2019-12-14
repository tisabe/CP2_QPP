#include <complex.h>

#include "structs.h"
#include "vmath.h"
#include "hamiltonian.h"

static void f_eta(double complex *out_eta, double complex *in_eta, parameters params) {
  hamiltonian_exp(out_eta, in_eta, params); // hier muessen noch parameter eingefuegt werden, vielleicht als typedef struct
  hamiltonian_exp(out_eta, out_eta, params);
  scalar_vec(out_eta, out_eta, (double complex) params.tau*params.tau/4, params.L);
  add_vec(out_eta, out_eta, in_eta, params.L);
}

static void f_psi(double complex *out_psi, double complex *in_psi, parameters params) {
  hamiltonian_exp(out_psi, in_psi, params);
  hamiltonian_exp(out_psi, out_psi, params);
  scalar_vec(out_psi, out_psi, (-1)*1I*params.tau/2, params.L);
  add_vec(out_psi, out_psi, in_psi, params.L);
}

void step_cn(double complex *out, double complex *in, , parameters params) {
  cg(out, f_eta, in, params);
  f_psi(out, out, params); // maybe this does not like out, out. not tested yet
}
