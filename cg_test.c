#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>

#include "structs.h"
#include "vmath.h"


void matrix_generator_exp(double complex *out, double complex *in, parameters params){
  // function rewritten by Tim for use of struct
  static double complex matrix[4] = {4.0, 1.0, 1.0, 3.0};
  mat_vec_mul(out, matrix, in, params.N);
}

void wikipedia_test_exp(){
  // function rewritten by Tim for use of struct
  // [[4,1],[1,3]]*[1/11,7/11]=[1,2]
  parameters params;
  params.N = 2;
  params.D = 2;
  params.max_iter = 100;
  params.tol = DBL_EPSILON;
  params.L = ipow(params.N,params.D);

  double complex *out = malloc(params.L*sizeof(double complex));
  double complex *in = malloc(params.N*sizeof(double complex));
  double complex *check_res = malloc(params.N*sizeof(double complex));
  double *diff = malloc(params.N*sizeof(double complex));

  in[0] = 1.0;
  in[1] = 2.0;

  cg(out, matrix_generator_exp, in, params);

  printf("[%lf+%lf i, %lf+%lf i]\n",creal(out[0]),cimag(out[0]),creal(out[1]),cimag(out[1]));

  matrix_generator_exp(check_res, out, params);
  for(int i=0; i<params.N; i++){
    diff[i] = cabs(check_res[i]-in[i]);
  }
  printf("Diff: [%e, %e]\n",diff[0],diff[1]);

  free(out);
  free(in);
  free(diff);
  free(check_res);
}


/*Include some good random number generator here*/

double random_double(){
    srand(time(NULL));
    return rand()/RAND_MAX;
}

/* This test creates a matrix filled with random complex numbers with cabs(z) <= 1/sqrt(2).
This matrix is then hermitian transposed and added to itself to make it a hermitian matrix with all entries cabs(z) <= 1.
Adding a diagonal matrix of N*eye(N) forces the matrix to be positive definite.*/

void random_matrix(double complex *out, double complex *in, parameters params){
    static double complex *matrix = NULL;
    static long int Nprev;

    //Intialise the random matrix once and keep it until the parameter N of the system changes
    if((matrix == NULL) || (Nprev != params.N)){
        free(matrix);
        matrix = malloc(params.L*sizeof(double complex));
        double complex *transp = malloc(params.L*sizeof(double complex));

        //Generate a random matrix filled with complex numbers of Re(z) in {-.5, .5} and Im(z) in {-.5,.5}.
        for(int i=0; i<params.L; i++){
            matrix[i] = (random_double() - 1) + (random_double() - 1) * I;
        }

        //Add the hermitian transposed of the random matrix to itself to make it hermitian
        hermtransp_quadr_matrix(transp, matrix, params.N);
        add_vec(out, transp, matrix, params.L);

        //Add the N*eye(N) to make the matrix positive definite
        for(int i=0; i<params.N; i++){
            out[i*params.N+i] += params.N;
        }

        free(transp);
        Nprev = params.N;
    }

    //Return the value of matrix*in matrix-vector-multiplication in the out vector
    mat_vec_mul(out, matrix, in, params.N);
}

void random_matrix_test(){

  parameters params;
  params.N = 5;
  params.D = 2;
  params.max_iter = 100;
  params.tol = DBL_EPSILON;
  params.L = ipow(params.N,params.D);

  double complex *out = malloc(params.L*sizeof(double complex));
  double complex *in = malloc(params.N*sizeof(double complex));
  double complex *check_res = malloc(params.N*sizeof(double complex));

  for(int i=0; i<params.N; i++){
    in[i] = random_double();
  }

  cg(out, random_matrix, in, params);

  printf("x = \n");
  for(int i=0; i<params.N; i++){
    printf("[%lf+%lf i]\n",creal(out[i]),cimag(out[i]));
  }

  random_matrix(check_res, out, params);

  printf("\n\nA*x-b = \n");
  for(int i=0; i<params.N; i++){
    printf("%e\n",cabs(check_res[i]-in[i]));
  }

  free(out);
  free(in);
  free(check_res);
}

int main(){
  wikipedia_test_exp();
  //random_matrix_test();
  return 0;
}
