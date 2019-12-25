#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "structs.h"
#include "vmath.h"
#include "potentials.h"
#include "geometry.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define _USE_MATH_DEFINES

/*  Compute the time propagation according to the strang splitting method.
     */



/* ************      To be tested        ************ */

void strang_splitting(double complex *out, double complex *in, parameters params) {
  long int N = params.N;
  unsigned int D = params.D;
  long int L = params.L;

  double sin_sum= 0;
  long int *coordinate = malloc(D* sizeof(long int));

  // Das braucht die fft aus gnu aus https://www.gnu.org/software/gsl/doc/html/fft.html#c.gsl_fft_complex_forward
  gsl_fft_complex_wavetable * wavetable;
  gsl_fft_complex_workspace * workspace;

  wavetable = gsl_fft_complex_wavetable_alloc (L);
  workspace = gsl_fft_complex_workspace_alloc (L);

  /* calculate eta according to equation (75) */
  for (int i=0; i<L; i++) {
    in[i]=cexp(- 1I/2 * params.tau * params.pot[i]) * in[i];
  }

	/* calculate eta_dft according to equation (76) */
  gsl_fft_complex_forward (in, 1, L, wavetable, workspace);

	/* calculate chi tilde (chi_dft) according to equation (77) */
	for (int i=0; i<L; i++) {
    for (int j=0; j<D; j++) {  /* calculating the sum in the exponential function */
      index2coord(coordinate, i, N, D);
      sin_sum += pow(sin(M_PI/N*coordinate[j]),2);
    }
    in[i]=cexp(- 1I * 2 * params.tau/params.mhat * sin_sum) * in[i];
  }

	/* calculate chi according to equation (78) */
	gsl_fft_complex_inverse (in, 1, L, wavetable, workspace);

  /* calculate psi(q+1) according to equation (79) */
  for (int i=0; i<L; i++) {
    out[i]=cexp(- 1I/2 * tauhat * params.pot[i]) * chi_dft[i];
  }

  free(coordinate);
  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);
}
