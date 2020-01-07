#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "structs.h"
#include "geometry.h"
#include "vmath.h"

void phi(double complex *out, parameters params, int type){
  long int *phi_coords;
  phi_coords = malloc(params.D*sizeof(long int));

	if (type==0){
		/* initialise wave function with ones */
		for (i=0;i<params.L:i++) {
			out[i]=1
		}
		/* wave function as product of cosines */
		for (i=0;i<params.L;i++) {
			index2coord(phi_coords, i, params.N, params.D);
			for (j=0;j<params.D;j++) {
				out[i]*=cosines(phi_coords[j]-params.N/2);
			}
		}
	}

	if (type==1) {
		/* wave function as square of coordinates with fixed imaginary parts*/
		for (i=0;i<params.L;i++) {
			index2coord(phi_coords, i, params.N, params.D);
			for (j=0;j<params.D;j++) {
				out[i]+=(phi_coords[j]-params.N/2)*(phi_coords[j]-params.N/2)*pow(I,j);
			}
		}
	}

	if(type==2){
		/* wave function as 1/abs(phi_coords) with fixed imaginary parts*/
		for (i=0;i<params.L;i++) {
			index2coord(phi_coords, i, params.N, params.D);
			for (j=0,j<params.D,j++){
				phi_coords[j]=(phi_coords[j]-params.N/2);
			}
			out[i]+=1/abs_vec(phi_coords,params.D)*pow(I,j);
		}
	}

	if(type==3){
		/* wave function as gaussian*/
		for (i=0;i<params.L;i++) {
			index2coord(phi_coords, i, params.N, params.D);
			long int *exponent = malloc(params.D*sizeof(long int));
			for (j=0,j<params.D,j++){
				phi_coords[j]=(phi_coords[j]-params.N/2);
			}
			sub_vec(exponent, phi_coords, params.phi0, params.L);
			gaussian_exp=dot_product(exponent, exponent, params.L)*(-1)/(params.sigma*params.sigma*2);

			out[i]=exp(gaussian_exp);
			free(exponent);
		}
	}

  double normalization_factor=abs_vec(out,params.L);
  for (i=0;i<params.L;i++) {
		/* normalization of the wave function */
		out[i]=out[i]/normalization_factor;
	}

	free(phi_coords);
}
