#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "nfft.h"

void ran_vec(double complex *vec, long int L) {
  /*Generates a random vector of size L and saves it in *vec.*/
  for (long int i=0; i<L; i++) {
    vec[i] = (double)rand()/RAND_MAX; // sets to random value between 0.0 and 1.0
  }
}
  int main(){
      
      int N= 25600;
    int D = 1;
 long int L = pow(N, D);  
 
 double complex *vec = malloc(L* sizeof(double complex));
 double complex *vec_fft = malloc(L* sizeof(double complex));
 double complex *difference = malloc(L* sizeof(double complex));
 ran_vec(vec, L);
 
 
 
 nfft(vec_fft,vec,N,D);
 nfft_inverse(vec_fft,vec_fft,N,D);
 
 FILE *f;
 f=fopen("nfft_testfile.txt","w");
 for(int i=0; i<L; i++) {
		difference[i] = vec[i]-vec_fft[i];
        fprintf(f,"%e \n", difference[i] );
	}
fclose(f);
	
	free(vec);
	free(vec_fft);
    free(difference);
    
    double complex *vec_sinus = malloc(L* sizeof(double complex));

    
    for(int i=0; i<L; i++) {
		vec_sinus[i] = sin(2*M_PI/48000*500*i);
	}
	nfft(vec_sinus,vec_sinus,N,D);
    FILE *f2;
    f2=fopen("nfft_testfile2.txt","w");
 for(int i=0; i<L; i++) {
     
        fprintf(f2,"%e \n", vec_sinus[i] );
     
	}
fclose(f2);
  }
