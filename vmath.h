#include <complex.h>

long int ipow(long int base, unsigned int exp);
double complex dot_product(double complex *v, double complex *w, long int L);
void scalar_vec(double complex *out, double complex *vec, double complex *s, long int L);
void add_vec(double complex *out, double complex *w, double complex *v, long int L);
void sub_vec(double complex *out, double complex *w, double complex *v, long int L);
void assign_vec(double complex *out, double complex *in, long int L);
double abs_vec(double complex *in, long int L);
void cg(double complex *out, void (*f)(double complex *, double complex *, long int L), double complex *in, int max_iter, double tol, long int N, unsigned int D);
