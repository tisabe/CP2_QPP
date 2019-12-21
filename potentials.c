#include <complex.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "vmath.h"
#include "geometry.h"

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
    double boxy=0;
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
    for (long int i=0; i<ipow(N, D); i++) {
        boxpotential[i]=boxpotential0( i, height, N, D);
    }
}



double potentialwell(long int i, int height, long int N, unsigned int D) {
    double boxy;
      if(i>N/4+N * (int) (i/N) && i< 3/4*N+N * (int) (i/N)){
        boxy=0;
      }else{
         boxy=height;
      }
    return boxy;
}
void well(double *wellpotential, int height, long int N, unsigned int D){
  for (long int i=0; i<ipow(N, D); i++) {
      wellpotential[i]=potentialwell(i, height, N, D);
  }
}
