#include <stdio.h>
#include <stdlib.h>

#include "structs.h"
#include "vmath.h"

/*This is just a little test to see if the structs file works as intended
and as a demonstration how to use it.*/

void func(parameters f_params) {
  printf("Testing parameters passed to function...\n");
  f_params.pot[1]++;
  printf("params.pot[0] = %.2e\n", f_params.pot[1]);
}

int main() {
  printf("Testing parameters struct...\n");
  unsigned int D = 3;
  long int N = 100;
  parameters params;
  params.N = N;
  printf("params.N = %d\n", N);
  params.D = D;
  printf("params.D = %d\n", D);
  //params.pot = malloc(ipow(N,D)*sizeof(double complex));
  params.pot = malloc(ipow(N,D)*sizeof(double));
  params.pot[0] = 10.0;
  printf("params.pot[0] = %.2e\n", params.pot[0]);
  func(params);
  return 0;
}
