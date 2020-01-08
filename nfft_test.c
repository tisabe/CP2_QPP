#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "nfft.h"
/* in this short file we make a plausibility test of the nfft to see if any fatal errors occur when we dont use the radix2  (in one dimension)*/
//function which generates a random vector "vec" of length L
void ran_vec(double complex *vec, long int L) {
  /*Generates a random vector of size L and saves it in *vec.*/
  for (long int i=0; i<L; i++) {
    vec[i] = (double)rand()/RAND_MAX; // sets to random value between 0.0 and 1.0
  }
}
  int main(){
      //declare paramters
      int N= 25600;
    int D = 1;
 long int L = pow(N, D);  
 
 double complex *vec = malloc(L* sizeof(double complex));
 double complex *vec_fft = malloc(L* sizeof(double complex));
 double complex *difference = malloc(L* sizeof(double complex));
 
 //generate random vector
 ran_vec(vec, L);
 //apply fft on the random vector
 nfft(vec_fft,vec,N,D);
 //apply the inverse fft 
 nfft_inverse(vec_fft,vec_fft,N,D);
 
 
 // calculate the difference between the original random vector and the on after applying the fft and inverse fft. we expect every element to be 0
 // write the result in .txt file "nfft_testfile.txt"
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
    
// generate a sin() function and apply the fft, the output in the textfile can be plotted to see if the receive desired frequency
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
