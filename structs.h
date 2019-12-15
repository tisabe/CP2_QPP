//typedef struct parameters_tag parameters;
struct parameters_tag;

typedef struct parameters_tag {
  long int N; // length of each axis, number of points along each axis
  unsigned int D; // number of dimensions
  long int L; // length of whole array
  double *pot; // space to store potential array
  double tau; // length of time step
  double mhat; // natural parameter m^ = m*epsilon*a^2/hbar^2
  double tol; // numerical tolerance for solvers (atm. for to cg method)
  double max_iter; // maximum number of iterations to perform (atm. for to cg method)
} parameters /*name to reference struct*/;