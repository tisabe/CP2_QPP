#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "indeces.h"
#include "vmath.h"

void phi(double complex *out, long int N, long int D, long int L, int type){
	long int *phi_coords;
  phi_coords = malloc(D*sizeof(long));

	if (type==0){
		/* initialise wave function with ones */
		for(i=0;i<L:i++){
			out[i]=1
		}
		/* wave function as product of cosines */
		for(i=0;i<L;i++){
			index2coord(phi_coords,i,N,D);
			for(j=0;j<D;j++){
				out[i]*=cosines(phi_coords[j]);
			}
		}
	}

	if (type==1){
		/* wave function as square of coordinates with fixed imaginary parts*/
		for(i=0;i<L;i++){
			index2coord(phi_coords,i,N,D);
			for(j=0;j<D;j++){
				out[i]+=phi_coords[j]*phi_coords[j]*pow(I,j);
			}
		}
	}

	if (type==2){
		/* wave function as 1/abs(phi_coords) with fixed imaginary parts*/
		for(i=0;i<L;i++){
			index2coord(phi_coords,i,N,D);
			out[i]+=1/abs_vec(phi_coords,D)*pow(I,j);
		}
	}

	double normalization_factor=abs_vec(out,L);
  for(i=0;i<L;i++){
		/* normalization of the wave function */
		out[i]=out[i]/normalization_factor;
	}

	free phi_coords
}
