#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>

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
    return (double)rand()/RAND_MAX;
}

/* This test creates a matrix filled with random complex numbers with cabs(z) <= 1/sqrt(2).
This matrix is then hermitian transposed and added to itself to make it a hermitian matrix with all entries cabs(z) <= 1.
Adding a diagonal matrix of N*eye(N) forces the matrix to be positive definite.*/

void random_matrix(double complex *out, double complex *in, parameters params){
    static double complex *matrix = NULL;
    static long int Nprev;

    //Initialise the random matrix once and keep it until the parameter N of the system changes
    if((matrix == NULL) || (Nprev != params.N)){
        free(matrix);
        matrix = malloc(ipow(params.N,2)*sizeof(double complex));
        double complex *transp = malloc(ipow(params.N,2)*sizeof(double complex));

        //Generate a random matrix filled with complex numbers of Re(z) in {-.5, .5} and Im(z) in {-.5,.5}.
        for(long int i=0; i<ipow(params.N,2); i++){
            matrix[i] = (random_double() - 0.5) + (random_double() - 0.5) * I;
        }

        for(int i=0; i<params.N; i++){
            matrix[params.N*i+i] = random_double() + 0.0 * I;
        }

        printf("Matrix generated\n");

        //Add the hermitian transposed of the random matrix to itself to make it hermitian
        // Proof: H = M + M^dagger
        //        H^dagger = (M + M^dagger)^dagger = M^dagger + (M^dagger)^dagger = M^dagger + M = H
        hermtransp_quadr_matrix(transp, matrix, params.N);

        printf("Matrix transposed\n");

        add_vec(matrix, transp, matrix, ipow(params.N,2));

        printf("Matrices added\n");

        //Add the N*eye(N) to make the matrix positive definite.
        //A diagonally dominant matrix with all diagonal entries > 0  are strictly positive definite
        //Diagonally dominant Matrix: |a_{ii}| > \sum_{j!=i} |a_{ij}| \forall i
        for(int i=0; i<params.N; i++){
            matrix[i*params.N+i] += params.N;
        }

        printf("Eye added\n\nCG starting\n");

        free(transp);
        Nprev = params.N;
    }

    //Return the value of matrix*in matrix-vector-multiplication in the out vector
    mat_vec_mul(out, matrix, in, params.N);
}

void random_matrix_test(){

  parameters params;
  params.N = 10000;
  params.D = 1;
  params.max_iter = 10000;
  params.tol = DBL_EPSILON;
  params.L = ipow(params.N,params.D);

  double complex *result_cg = malloc(params.N*sizeof(double complex));
  double complex *start_vec = malloc(params.N*sizeof(double complex));
  double complex *check_res = malloc(params.N*sizeof(double complex));

  double max_error = 0.0;

  for(int i=0; i<params.N; i++){
    start_vec[i] = (double complex)(rand()-0.5*RAND_MAX + (rand()-0.5*RAND_MAX) * I);
  }

  cg(result_cg, random_matrix, start_vec, params);

  random_matrix(check_res, result_cg, params);

  printf("\nMax_Error(A*x-b) = \n");
  for(int i=0; i<params.N; i++){
    if(cabs(check_res[i]-start_vec[i]) > max_error){
        max_error = cabs(check_res[i]-start_vec[i]);
    }
  }
  printf("%e\n",max_error);

  free(result_cg);
  free(start_vec);
  free(check_res);
}

int main(){
  srand(time(NULL));
  double time_spent = 0.0;
  clock_t begin = clock();
  //wikipedia_test_exp();
  random_matrix_test();
  clock_t end = clock();
  time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time elpased is %f seconds\n\n\n", time_spent);
  return 0;
}
