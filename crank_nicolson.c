#include <complex.h>

#include "vmath.h"
#include "hamiltonian.h"

static void f_eta(double complex *out_eta, double complex *in_eta, long int L, double tau) {
  hamiltonian(out_eta, in_eta, ); // hier muessen noch parameter eingefuegt werden, vielleicht als typedef struct
  hamiltonian(out_eta, out_eta, );
  scalar_vec(out_eta, out_eta, (double complex) tau*tau/4, L);
  add_vec(out_eta, out_eta, in_eta, L);
}

static void f_psi(double complex *out_psi, double complex *in_psi, long int L, double tau) {
  hamiltonian(out_psi, in_psi, );
  hamiltonian(out_psi, out_psi, );
  scalar_vec(out_psi, out_psi, (-1)*1Im*tau/2, L);
  add_vec(out_psi, out_psi, in_psi, L);
}

void step_cn(double complex *out, double complex *in, long int N, unsigned int D, double tau, int max_iter, double tol) {
  long int L = ipow(N, D);
  cg(out, f_eta, in, max_iter, tol, L);
  f_psi(out, out, L, tau);
}
