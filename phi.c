#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
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
		/* wave function as product of sines */
		for(i=0;i<L;i++){
			index2coord(phi_coords,i,N,D);
			for(j=0;j<D;j++){
				out[i]*=sin(phi_coords[j]);
				}
			}
		double normalization_factor=abs_vec(out,L);
		for(i=0;i<L;i++){
			/* normalization of the wave function */
			out[i]=out[i]/normalization_factor;
			}
	if (type==1){
		/* wave function as square of coordinates */
		for(i=0;i<L;i++){
			index2coord(phi_coords,i,N,D);
			for(j=0;j<D;j++){
				out[i]+=j*j*ipow(I);
				}
		}
	if (type==2){
		
		}
	}


