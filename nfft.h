#ifndef NFFT_H
#define NFFT_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

void nfft(double complex *out, double complex *in, int N, int D);
void nfft_inverse(double complex *out, double complex *in, int N, int D);

#endif
