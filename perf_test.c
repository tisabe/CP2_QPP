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

  double offset = 0.0;
  double scaling = 1.0;
  double omega = 8.3e14; //Hz
  double mass_H = 1.67e-27; //kg
  double hbar = 1.054571817e-34; //Js
  params.epsilon = hbar * omega;
  params.a = 2e-13;
  params.mhat = pow(params.a/hbar,2)*mass_H*params.epsilon;
  params.khat = params.mhat * pow(hbar*omega/params.epsilon,2);

  FILE *perf_output_file;
  perf_output_file = fopen("perf_test_output.txt","w");
  fprintf(perf_output_file, "N\tStatic nn array\tCalc nn every time\tSingle loop\tParallel\n");

  for(double N = 8; N <= pow(2,21); N *= sqrt(2)){
    params.N = (long int)N;
    params.D = 1;
    params.L = ipow(params.N, params.D);

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

    double t_duration, t_start;
    double reg_min = 1e10;
    double slow_min = 1e10;
    double single_loop_min = 1e10;
    double parallel_min = 1e10;

    printf("N=%ld, D=%d\n", params.N, params.D);
    //printf("Regular\t\tSlow\t\tSingle loop\tParallel\n\n");
    for(int i=0; i < 100; i++){
      t_start = omp_get_wtime();
      hamiltonian(out_wf, start_wf, params);
      t_duration = omp_get_wtime() - t_start;
      if(t_duration < reg_min){reg_min = t_duration;}
      //printf("%e\t", t_duration);

      t_start = omp_get_wtime();
      hamiltonian_slow(out_wf, start_wf, params);
      t_duration = omp_get_wtime() - t_start;
      if(t_duration < slow_min){slow_min = t_duration;}
      //printf("%e\t", t_duration);

      t_start = omp_get_wtime();
      hamiltonian_single_loop(out_wf, start_wf, params);
      t_duration = omp_get_wtime() - t_start;
      if(t_duration < single_loop_min){single_loop_min = t_duration;}
      //printf("%e\t", t_duration);

      t_start = omp_get_wtime();
      hamiltonian_parallel(out_wf, start_wf, params);
      t_duration = omp_get_wtime() - t_start;
      if(t_duration < parallel_min){parallel_min = t_duration;}
      //printf("%e\t", t_duration);
      //printf("\n");
    }

    printf("Minimum time:\n");
    printf("Static nn array\tCalc nn every time\tSingle loop\tParallel\n");
    printf("%e\t%e\t%e\t%e\n\n\n",reg_min,slow_min,single_loop_min,parallel_min);
    fprintf(perf_output_file, "%li\t%e\t%e\t%e\t%e\n",params.N,reg_min,slow_min,single_loop_min,parallel_min);

    free(start_wf);
    free(out_wf);
    free(coords);
    free(params.pot);
  }

  fclose(perf_output_file);
}
