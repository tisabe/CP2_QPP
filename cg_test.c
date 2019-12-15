#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>

#include "structs.h"
#include "vmath.h"

/* This test creates a matrix filled with random complex numbers with cabs(z) <= 1/sqrt(2).
This matrix is then hermitian transposed and added to itself to make it a hermitian matrix with all entries cabs(z) <= 1.
Adding a diagonal matrix of N*eye(N) forces the matrix to be positive definite due to properties of strictly diagonally dominant matrices.*/

//Function to return random doubles between 0.0 and 1.0
double random_double(){
    return (double)rand()/RAND_MAX;
}

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

        /* Add the hermitian transposed of the random matrix to itself to make it hermitian
        Proof: H = M + M^dagger
               H^dagger = (M + M^dagger)^dagger = M^dagger + (M^dagger)^dagger = M^dagger + M = H */
        hermtransp_quadr_matrix(transp, matrix, params.N);

        printf("Matrix transposed\n");

        add_vec(matrix, transp, matrix, ipow(params.N,2));

        printf("Matrices added\n");

        /* Add the N*eye(N) to make the matrix positive definite.
        A diagonally dominant matrix with all diagonal entries > 0  are strictly positive definite.
        A diagonally dominant matrix is defined by: |a_{ii}| > \sum_{j!=i} |a_{ij}| \forall i */
        for(int i=0; i<params.N; i++){
            matrix[i*params.N+i] += params.N;
        }

        printf("N * Unit matrix added\n\nCG starting\n");

        free(transp);
        Nprev = params.N;
    }

    //Return the value of matrix*in matrix-vector-multiplication in the out vector
    mat_vec_mul(out, matrix, in, params.N);
}

void random_matrix_test(){
  //Set simulation parameters
  parameters params;

  //Set array size to run the test with
  printf("Please specify the length of the array to run the test with (memory scales with N**2). N = ");
  scanf("%li",&params.N);               //Choose array length
  params.D = 1;                         //Set number of dimensions to 1. Only meaningful for interpretation of the flat output array as a matrix.
  params.max_iter = 100;                //Set the maximum number of iterations of the cg algorithm
  params.tol = DBL_EPSILON;             //Set the tolerance level to quit the cg algorithm. Best possible: DBL_PRECISION (macro from float.h)
  params.L = ipow(params.N,params.D);

  //Allocate arrays for computation
  double complex *result_cg = malloc(params.N*sizeof(double complex));
  double complex *start_vec = malloc(params.N*sizeof(double complex));
  double complex *check_res = malloc(params.N*sizeof(double complex));
  //Variable to store maximum relative error of the result. See below
  double max_error = 0.0;

  //Generate a random vector b
  for(int i=0; i<params.N; i++){
    start_vec[i] = (double complex)(rand()-0.5*RAND_MAX + (rand()-0.5*RAND_MAX) * I);
  }

  //Call the conjugate gradient algorithm to calculate x in Ax=b. Store the result in result_cg
  cg(result_cg, random_matrix, start_vec, params);

  //Calculate the matrix vector product of the random matrix A with the solution vector x from cg
  random_matrix(check_res, result_cg, params);

  //Calculate the maximum relative error of the solution with respect to the analytical solution.
  //max_error = max( ((Ax)_i-b_i)/b_i )
  for(int i=0; i<params.N; i++){
    if(cabs(check_res[i]/start_vec[i]-1) > max_error){
        max_error = cabs(check_res[i]/start_vec[i]-1);
    }
  }
  printf("\nMaximum Relative Error: max(((Ax)_i-b_i)/b_i) = %e\n",max_error);

  free(result_cg);
  free(start_vec);
  free(check_res);
}

int main(){
  //Set seed for the random number generator
  srand(time(NULL));

  printf("\n*******Testing the CG algorithm with a randomly generated input*******\n\n");
  //Start testing CG with a random hermitian strictly positive definite matrix A and a random complex vector b to solve Ax=b for x
  random_matrix_test();

  return 0;
}
