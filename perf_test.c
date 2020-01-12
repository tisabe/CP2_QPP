#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

#include "structs.h"
#include "geometry.h"
#include "hamiltonian.h"
#include "vmath.h"


int main() {
  parameters params;
  params.N = 100000;
  params.D = 1;
  params.L = ipow(params.N, params.D);
  double offset = 0.0;
  double scaling = 1.0;
  double omega = 8.3e14; //Hz
  double mass_H = 1.67e-27; //kg
  double hbar = 1.054571817e-34; //Js
  params.epsilon = hbar * omega;
  params.mhat = pow(params.a/hbar,2)*mass_H*params.epsilon;
  params.khat = params.mhat * pow(hbar*omega/params.epsilon,2);

  double complex *start_wf = malloc(params.L * sizeof(double complex));
  double complex *out_wf = malloc(params.L * sizeof(double complex));
  long int *coords = malloc(params.D * sizeof(long int));
  params.pot = malloc(params.L * sizeof(double complex));

  parameters *p = &params;
  gen_pot_harmonic(p, omega);

  for(long int i=0; i<(params.L); i++){
    index2coord(coords, i, params.N, params.D);
    start_wf[i] = sqrt(sqrt(scaling))*pow(sqrt(params.a),params.D) * sqrt(sqrt(params.mhat/M_PI)/params.a)* exp(-.5*scaling*params.mhat*pow((coords[0]-params.N/2+offset),2));
  }

  double t_end, t_start;
  printf("N=%ld, D=%d\n", params.N, params.D);
  for(int i=0; i<10; i++){
    t_start = omp_get_wtime();
    hamiltonian(out_wf, start_wf, params);
    t_end = omp_get_wtime() - t_start;
    printf("%e\t", t_end);

    t_start = omp_get_wtime();
    hamiltonian_slow(out_wf, start_wf, params);
    t_end = omp_get_wtime() - t_start;
    printf("%e\t", t_end);

    t_start = omp_get_wtime();
    hamiltonian_single_loop(out_wf, start_wf, params);
    t_end = omp_get_wtime() - t_start;
    printf("%e\t", t_end);

    t_start = omp_get_wtime();
    hamiltonian_parallel(out_wf, start_wf, params);
    t_end = omp_get_wtime() - t_start;
    printf("%e\t", t_end);
    printf("\n");
  }


}
