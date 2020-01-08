#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "structs.h"
#include "vmath.h"
#include "nfft.h"
#include "hamiltonian.h"
#include "geometry.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define _USE_MATH_DEFINES

void step_euler(double complex *out, double complex *in, parameters params) {
    /*  Compute the time propagation according to the Euler method.
      Note: Normalisation of the wave function is not conserved!! */

    // Allocate memory for a temporary array of length L
    static double complex *ham_out = NULL;
    static long int lPrev;
    if((ham_out == NULL) || (lPrev != params.L)){
      free(ham_out);
      ham_out = malloc(params.L*sizeof(double complex));
      lPrev = params.L;
    }

    // Calculate H(psi)
    hamiltonian(ham_out, out, params);

    // psi^{q+1} = psi^q - i tau H(psi^q)
    for(long int i=0; i<params.L; i++){
        out[i] = in[i] - 1I * params.tauhat * ham_out[i];
    }
}

/*These are the functions required to run the crank-nicolson integration: */

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
  double complex *temp2 = malloc(params.L * sizeof(double complex));
  hamiltonian(temp, in_psi, params);
  hamiltonian(out_psi, temp, params);
  scalar_vec(out_psi, out_psi, (-1)*params.tauhat*params.tauhat/4, params.L);
  scalar_vec(temp2, temp, -1I*params.tauhat, params.L);
  add_vec(out_psi, out_psi, in_psi, params.L);
  add_vec(out_psi, out_psi, temp2, params.L);
  free(temp);
  free(temp2);
}

void step_cn(double complex *out, double complex *in, parameters params) {
  double complex *temp = malloc(params.L * sizeof(double complex));
  cg(temp, f_eta, in, params);
  f_psi(out, temp, params);
  free(temp);
}

void step_strang(double complex *out, double complex *in, parameters params) {
  /* Computes a time step of the wavefunction according to the
  Strang-Splitting-method, with the parameters stored in params. */
  long int N = params.N;
  unsigned int D = params.D;
  long int L = params.L;

  double sin_sum= 0;
  long int *coordinate = malloc(D* sizeof(long int));

  // https://www.gnu.org/software/gsl/doc/html/fft.html#c.gsl_fft_complex_forward

  /* calculate eta according to equation (75) */
  for (int i=0; i<L; i++) {
    in[i]=cexp(- 1I/2 * params.tauhat * params.pot[i]) * in[i];
  }

	/* calculate eta_dft according to equation (76) */
  nfft(in,in,N,D);

	/* calculate chi tilde (chi_dft) according to equation (77) */
	for (int i=0; i<L; i++) {
    for (int j=0; j<D; j++) {  /* calculating the sum in the exponential function */
      index2coord(coordinate, i, N, D);
      sin_sum += pow(sin(M_PI/N*coordinate[j]),2);
    }
    in[i]=cexp(- 1I * 2 * params.tauhat/params.mhat * sin_sum) * in[i];
  }

	/* calculate chi according to equation (78) */
	nfft_inverse(in,in,N,D);
  /* calculate psi(q+1) according to equation (79) */
  for (int i=0; i<L; i++) {
    out[i]=cexp(- 1I/2 * params.tauhat * params.pot[i]) * in[i];
  }

  free(coordinate);

}
