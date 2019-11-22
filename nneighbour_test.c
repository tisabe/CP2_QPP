#include <stdio.h>
#include "nneighbour.h"


int main() {
  printf("Testing the nneighbour function\n");
  printf("i:, \t nni:\n");

  long int N = 7;
  unsigned int D = 2;
  long int nni = 0;

  for (int i=0; i<20; i++){
    nni = nneighbour(i, 0, 1, N, D);
    printf("%d, \t %d\n", i, nni);
  }
  return 0;
}
