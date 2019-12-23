#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "vmath.h"
#include "laplacian.h"

double run_serial(long int N, unsigned int D) {
  double complex *in = malloc(ipow(N,D)*sizeof(double complex));
  double complex *out = malloc(ipow(N,D)*sizeof(double complex));
  int n = 11;
  int total = 0;
  double msec;
  clock_t diff = 0;
  //int total = 0;
  for (int i=0; i<n; i++) {
    clock_t start = clock(), diff;
    laplacian(out, in, N, D);
    diff += clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    if(i>0){
      total += msec;
    }
  }
  free(in);
  free(out);
  return total/(double)(n-1);
}

double run_parallel(long int N, unsigned int D) {
  double complex *in = malloc(ipow(N,D)*sizeof(double complex));
  double complex *out = malloc(ipow(N,D)*sizeof(double complex));
  int n = 11;
  int total = 0;
  double msec;
  clock_t diff = 0;
  //int total = 0;
  for (int i=0; i<n; i++) {
    clock_t start = clock(), diff;
    laplacian_mp(out, in, N, D);
    diff += clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    if(i>0){
      total += msec;
    }
  }
  free(in);
  free(out);
  return total/(double)(n-1);
}

int main() {
  unsigned int D = 1;
  printf("N\t serial time\t parallel time\n");
  double time_s;
  double time_p;
  //omp_set_dynamic(0);
  //omp_set_num_threads(omp_get_num_procs());
  for (long int n = 1001; n < 1000000; n += 100000) {
    time_s = run_serial(n, D);
    time_p = run_parallel(n, D);
    printf("%d\t %f\t %f \t %f\n", n, time_s, time_p, time_s/time_p);
  }
  return 0;
}
