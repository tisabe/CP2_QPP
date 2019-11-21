/*totally unfinished and untested*/

#include <stdio.h>

long int coord2index(long int *coord, long int N, unsigned int D){
  long int index = 0;

  for (i=0; i<D; i++){
    index += coord[i]*ipow(N,i)
  }

  return index;
}

void index2coord(long int *coord, long int index, long int N, unsigned int D){
  /*coordinate vector needs to be returned as pointer in C*/
  //static int coord[D];
  long int b = 0;

  for (i=0; i<D; i++) {
    b = ipow(N, D-1-i);
    coord[D-1-i] = index/b;
    index %= b;
  }
}

long int ipow(long int base, unsigned int exp){
/*integer power function so math does not need to be included*/
/*peer reviewed by frohlofl*/
    int result = 1;

    for (i=0; i<exp; i++){
      result *= base;
    }

    return result;
}
