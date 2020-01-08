#ifndef VMATH_H
#define VMATH_H

#include <complex.h>

long int ipow(long int base, unsigned int exp);
double complex dot_product(double complex *v, double complex *w, long int L);
void scalar_vec(double complex *out, double complex *vec, double complex s, long int L);
void mul_compl_vec_real_scal(double complex *out, double complex *vec, double s, long int L);
void mul_element(double complex *out, double complex *w, double complex *v, long int L);
void add_vec(double complex *out, double complex *w, double complex *v, long int L);
void sub_vec(double complex *out, double complex *w, double complex *v, long int L);
void assign_vec(double complex *out, double complex *in, long int L);
void set_zero(double complex *in, long int L);
double abs_vec(double complex *in, long int L);
long int fact(int n);
void mat_vec_mul(double complex *out, double complex *mat, double complex *vec, long int N);
void hermtransp_quadr_matrix(double complex *out, double complex *in, long int N);
void cg(double complex *out, void (*f)(double complex *, double complex *, parameters), double complex *in, parameters params);

#endif
