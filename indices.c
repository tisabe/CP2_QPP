/*totally unfinished and untested*/

#include <stdio.h>
#include "vmath.h"
/*
    int result = 1;

    for (int i=0; i<exp; i++){
      result *= base;
    }

    return result;
}*/

long int coord2index(long int *coord, long int N, unsigned int D){
  long int index = 0;
  long int shift = (N-1)/2; // coordinate shift amount to have coordinate origin at middle of axis

  if (((N+1)%2) != 0) {
    printf("Error: Length of axis N should be uneven");
    return -1;
  }

  for (int i=0; i<D; i++){
    index += (coord[i]+shift)*ipow(N,i);
  }

  return index;
}

void index2coord(long int *coord, long int index, long int N, unsigned int D){
  /*coordinate vector needs to be returned as pointer in C*/
  //static int coord[D];
  if (((N+1)%2) != 0) {
    printf("Error: Length of axis N should be uneven");
  }

  long int b = 0;
  long int shift = (N-1)/2; // coordinate shift amount to have coordinate origin at middle of axis

  for (int i=0; i<D; i++) {
    b = ipow(N, D-1-i);
    coord[D-1-i] = index/b - shift;
    index %= b;
  }
}
