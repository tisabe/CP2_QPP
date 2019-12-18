#include "vmath.h"
#include "indices.h"
#include <complex.h>
#include <stdlib.h>
#include <math.h>

double harmonic0(long int index, double omega, long int N, unsigned int D) {
    long int * coordinates = malloc(D* sizeof(long int));
    index2coord(coordinates,index, N, D);
    double potential = 0;
    for (int i=0; i<D; i++) {
        potential += pow(coordinates[i], 2);
  }
    potential=0.5*pow(omega, 2)*potential;
    free(coordinates);
    return potential;
}
void harmonic(double *harmonicpotential, double omega, long int N, unsigned int D){
    for (int i=0; i<ipow(N, D); i++) {
        harmonicpotential[i]=harmonic0( i, omega,N, D);
    }
}



double boxpotential0(long int index, int height, long int N, unsigned int D) {
    long int * coordinates = malloc(D* sizeof(long int));
    index2coord(coordinates,index, N, D);
    long int boxy=0;
    for (int i=0; i<D; i++) {
        if(fabs(coordinates[i])<N/4){
            boxy=boxy;
        }else{
            boxy=height;
            break;
            }
    free(coordinates);
    return boxy;
  }
}
void box(double *boxpotential, int height, long int N, unsigned int D){
    for (int i=0; i<ipow(N, D); i++) {
        boxpotential[i]=boxpotential0( i, height, N, D);
    }
}



double potentialwell(int height, long int N, unsigned int D) {
    long int * boxy= malloc(ipow(N, D) * sizeof(long int));
    for (int i=0; i<D; i++) {
        if(i>N/4+N * (int) (i/N) && i< 3/4*N+N * (int) (i/N)){
            boxy[i]=0;
        }else{
            boxy[i]=height;
            }
    return boxy;
  }
}
void well(double *wellpotential, int height, long int N, unsigned int D){
  for (int i=0; i<ipow(N, D); i++) {
      wellpotential[i]=potentialwell(height, N, D);
  }
}
