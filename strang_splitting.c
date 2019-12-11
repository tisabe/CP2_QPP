#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "hamiltonian.h"        //Still to be implemented

/*  Compute the time propagation according to the strang splitting method.
     */



/* ************      To be tested        ************ */

void strang_method(double complex *out, double complex *in, long int N, unsigned int D, long int L, double tauhat, double mhat, int ext_potential_type, double total_time){
    
/*Calculate the different potentials*/
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

double complex *eta= malloc(ipow(N, D)*sizeof(double));
double complex *eta_dft= malloc(ipow(N, D)*sizeof(double));
double complex *chi_dft= malloc(ipow(N, D)*sizeof(double));
double *sin_sum= malloc(1*sizeof(double));
long int * coordinate = malloc(D* sizeof(long int));
double complex chi= inverse_fft(chi_dft);

double complex *psi= malloc(ipow(N, D)*sizeof(double));
psi=in;

/* Das braucht die fft aus gnu aus https://www.gnu.org/software/gsl/doc/html/fft.html#c.gsl_fft_complex_forward*/
gsl_fft_complex_wavetable * wavetable;
gsl_fft_complex_workspace * workspace;

wavetable = gsl_fft_complex_wavetable_alloc (L);
workspace = gsl_fft_complex_workspace_alloc (L);


for(int t=0; (t * tauhat) < total_time; t++){

	/* calculate eta according to equation (73) */
	for (int i=0; i<ipow(N, D); i++) {
        	eta[i]=cexp(- 1Im/2 * tauhat * phi_potential[i]) * psi[i]; 
    	return eta;
 	 }

	
	/* calculate eta_dft according to equation (74) */
	eta_dft=gsl_fft_complex_forward (eta, 1, L, wavetable, workspace);

	
	/* calculate chi tilde (chi_dft) according to equation (75) */
	for (int i=0; i<ipow(N, D); i++) {
		for (int j=0; j<D; j++) {  /* calculatin the sum in the exponential function */
			index2coord(coordinate,i, N, D);
			sin_sum += sin(pi/N*coordinate[j]) * sin(pi/N*coordinate[j]);
		return sin_sum;
		}

    	    chi_dft[i]=cexp(- 1Im * 2 * tauhat/mhat * sin_sum) * eta_dft[i]; 
  	  return chi_dft;
 	 }

	/* calculate chi according to equation (76) */
	chi=gsl_fft_complex_inverse (chi_dft, 1, L, wavetable, workspace);
	

	for (int i=0; i<ipow(N, D); i++) {
  	      psi[i]=cexp(- 1Im/2 * tauhat * phi_potential[i]) * chi[i]; 
  	  return psi;
 	 }

return psi;

free(eta);
free(eta_dft);
free(chi_dft);
free(chi);
free(sin_sum);
free(coordinate);

 gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);
}
}




