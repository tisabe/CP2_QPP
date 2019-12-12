#include <stdlib.h>
#include <complex.h>
#include "vmath.h" // only for ipow, may be improved
/*
Here all structs used in the QPP problem will be defined.
*/

struct parameters {
  long int N; // length of each axis, number of points along each axis
  unsigned int D; // number of dimensions
  double complex *pot = malloc(10000000*sizeof(double));
  //double complex *pot = malloc(ipow(N,D)*sizeof(double complex)); // array potential
  double tau; // length of time step
  double mhat; // natural parameter m^ = m*epsilon*a^2/hbar^2
  double tolerance; // numerical tolerance for solvers (atm. for to cg method)
  double max_iter; // maximum number of iterations to perform (atm. for to cg method)
};
