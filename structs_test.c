#include <stdio.h>
#include <complex.h>

#include "vmath.c"
#include "structs.h"

/*This is just a little test to see if the structs file works as intended*/

int main(int argc, char const *argv[]) {
  printf("Testing parameters struct...\n");
  unsigned int D = 3;
  long int N = 100;
  struct parameters params;
  params.N = N;
  printf("params.N = %d\n", N);
  params.D = D;
  printf("params.D = %d\n", D);
  params.pot = malloc(ipow(N,D)*sizeof(double complex));
  params.pot[0] = (double complex) 10.0;
  printf("params.pot[0] = %.2e\n", params.pot[0]);
  return 0;
}
