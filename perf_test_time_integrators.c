#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <float.h>

#include "structs.h"
#include "geometry.h"
#include "hamiltonian.h"
#include "vmath.h"
#include "integrators.h"
#include "observables.h"

#define _USE_MATH_DEFINES

extern double time_ham_total;

int main() {
  double time_total;

  time_total = omp_get_wtime();

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
  params.tauhat = 1e-3;
  double simulation_duration = 1e1;

  FILE *perf_output_file;
  perf_output_file = fopen("perf_test_output.txt","w");
  fprintf(perf_output_file, "N\tStatic nn array\tCalc nn every time\tSingle loop\tParallel\n");

  params.N = 2000;
  params.D = 1;
  params.L = ipow(params.N, params.D);

  params.tol = DBL_EPSILON;
  params.max_iter = params.L;

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

  for(long int t=0; (t * params.tauhat) < simulation_duration; t++){
    if(t % (long)(1e4) == 0){
        printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n", t*params.tauhat, cabs(obs_norm(start_wf, params)) - 1, creal(obs_E(start_wf,params)), creal(obs_x(start_wf, 0, params)), creal(obs_p(start_wf, 0, params)), creal(obs_delta_x(start_wf, params)),creal(obs_delta_p(start_wf, params)));
    }
    step_cn_timed(out_wf, start_wf, params);
    assign_vec(start_wf, out_wf, params.L);
  }

  printf("\nTime in Hamiltonian = %lf s\n", time_ham_total);

  free(start_wf);
  free(out_wf);
  free(coords);
  free(params.pot);

  fclose(perf_output_file);

  time_total = omp_get_wtime() - time_total;

  printf("\nTotal time = %lf s\n",time_total);
}
