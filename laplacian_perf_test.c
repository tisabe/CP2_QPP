#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "vmath.h"
#include "laplacian.h"

int run_serial(long int N, unsigned int D) {
  double complex *in = malloc(ipow(N,D)*sizeof(double complex));
  double complex *out = malloc(ipow(N,D)*sizeof(double complex));

  clock_t start = clock(), diff;
  laplacian(out, in, N, D);
  diff = clock() - start;
  free(in);
  free(out);
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  return msec;
}

int run_parallel(long int N, unsigned int D) {
  double complex *in = malloc(ipow(N,D)*sizeof(double complex));
  double complex *out = malloc(ipow(N,D)*sizeof(double complex));

  clock_t start = clock(), diff;
  laplacian_mp(out, in, N, D);
  diff = clock() - start;
  free(in);
  free(out);
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  return msec;
}

int main() {
  unsigned int D = 1;
  printf("N\t serial time\t parallel time\n");
  int time_s;
  int time_p;
  //omp_set_dynamic(0);
  //omp_set_num_threads(omp_get_num_procs());
  for (long int n = 1001; n < 1000000; n += 10000) {
    time_s = run_serial(n, D);
    time_p = run_parallel(n, D);
    printf("%d\t %d\t %d\n", n, time_s, time_p);
  }
  return 0;
}
