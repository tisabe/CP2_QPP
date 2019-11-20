/*totally unfinished and untested*/

#include <stdio.h>

unsigned int coord2index(int coord[D], unsigned int N, unsigned int D){
  unsigned int index = 0;

  for (i=0; i<D; i++){
    index += coord[i]*ipow(N,i)
  }

  return index;
}

int * index2coord(int index, unsigned int N, unsigned int D){
  /*coordinate vector needs to be returned as pointer in C*/
  static int coord[D];
  int b = 0;

  for (i=0; i<D; i++) {
    b = ipow(N, D-1-i);
    coord[D-1-j] = index/b;
    index %= b;
  }

  return coord;
}

int ipow(int base, int exp){
/*integer power function so math does not need to be included*/
    int result = 1;

    for (i=0; i<exp; i++){
      result *= base;
    }

    return result;
}
