#ifndef STRUCTS_H
#define STRUCTS_H

#include <complex.h>

//typedef struct parameters_tag parameters;
struct parameters_tag;

typedef struct parameters_tag {
  long int N; // length of each axis, number of points along each axis
  unsigned int D; // number of dimensions
  long int L; // length of whole array
  double complex *pot; // space to store potential array
  double tauhat; // length of time step
  double mhat; // natural parameter m^ = m*epsilon*a^2/hbar^2
  double khat; // natural parameter k^ = m^ * hbar^2 *omega^2 / epsilon^2
  double tol; // numerical tolerance for solvers (atm. for to cg method)
  long int max_iter; // maximum number of iterations to perform (atm. for to cg method)
  double epsilon; // energy scale
  double a; // distance between two lattice points
} parameters /*name to reference struct*/;

#endif
