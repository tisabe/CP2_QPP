#ifndef STRUCTS_H
#define STRUCTS_H

#include <complex.h>

//typedef struct parameters_tag parameters;
struct parameters_tag;

typedef struct parameters_tag {
  long int N; // length of each axis, number of points along each axis
  unsigned int D; // number of dimensions
  long int L; // length of whole array
  double a; // length scale
  double complex *pot; // space to store potential array
  double tau; // length of time step
  double mhat; // natural parameter m^ = m*epsilon*a^2/hbar^2
  double tol; // numerical tolerance for solvers (atm. for to cg method)
  long int max_iter; // maximum number of iterations to perform (atm. for to cg method)
  double total_time; //Time that the simulation should run for
  double epsilon; //Energy scale
  int ext_potential_type; //External potential type
  double sigma; //Sigma in a gaussian distribution
  long int *phi0; //Center in a gaussian distribution
} parameters /*name to reference struct*/;

#endif
