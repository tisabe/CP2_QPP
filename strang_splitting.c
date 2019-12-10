#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "hamiltonian.h"        //Still to be implemented

/*  Compute the time propagation according to the strang splitting method.
     */



/* ************      To be tested        ************ */

void strang_method(double complex *out, double complex *in, long int N, unsigned int D, long int L, double tauhat, double, mhat, double total_time){
    
/*Calculate the harmonic potential*/
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

complex double *eta= malloc(ipow(N, D)*sizeof(double));
double *eta_dft= malloc(ipow(N, D)*sizeof(double));
complex double *chi_dft= malloc(ipow(N, D)*sizeof(double));
double *sin_sum= malloc(1*sizeof(double));
long int * coordinate = malloc(D* sizeof(long int));
complex double chi= inverse_fft(chi_dft);

double complex psi=in

for(int t=0; (t * tauhat) < total_time; t++){

	/* calculate eta according to equation (73) */
	for (int i=0; i<ipow(N, D); i++) {
        	eta[i]=exp(- 1Im/2 * tauhat * phi_potential[i]) * psi[i]; /* exp() might not take complex double */
    	return eta;
 	 }

	
	/* calculate eta_dft according to equation (74) */
	eta_dft=fft(eta)...

	
	/* calculate chi tilde according to equation (75) */
	for (int i=0; i<ipow(N, D); i++) {
		for (int j=0; i<D; i++) {
			index2coord(coordinate,i, N, D);
			sin_sum += sin(pi/N*coordinate[j]) * sin(pi/N*coordinate[j]);
		return sin_sum;
		}

    	    chi_dft[i]=exp(- 1Im * 2 * tauhat/mhat * sin_sum) * eta_dft[i]; /* exp() might not take complex double */
  	  return chi_dft;
 	 }

	/* calculate chi according to equation (76) */
	chi=inverse_fft(chi_dft)
	

	for (int i=0; i<ipow(N, D); i++) {
  	      psi[i]=exp(- 1Im/2 * tauhat * phi_potential[i]) * chi[i]; /* exp() might not take complex double */
  	  return psi;
 	 }

return psi;

free(eta);
free(eta_dft);
free(chi_dft);
free(chi);
free(sin_sum);
free(coordinate);
}
}




