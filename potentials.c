#include  "vmath.h"
#include <complex.h>
#include <stdlib.h>
#include <math.h>

double harmonic0(long int index, int omega, int N, int D) {
    long int * coordinates = malloc(D* sizeof(long int));
    index2coord(coordinates,index, N, D);
    double potential = 0;
    for (int i=0; i<D; i++) {
        potential = potential+ipow(coordinates[i], 2);
  }
    potential=0.5*ipow(omega, 2)*potential;
    free(coordinates);
    return potential;
}

void harmonic(double * harmonicpotential, int omega, int N, int D){
    for (int i=0; i<ipow(N, D); i++) {
        harmonicpotential[i]=harmonic0( i, omega,N, D);
    }
}

double boxpotential0(long int index, int omega, int N, int D) {
    long int * coordinates = malloc(D* sizeof(long int));
    index2coord(coordinates,index, N, D);
    int count=1;
    long int boxy=0;
    for (int i=0; i<D; i++) {
        if(faps(coordinates[i])<N/4){
            boxy=boxy+boxy;
        else
            boxy=height;
            }
    return boxy;
