#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include <complex.h>

void step_cn(double complex *out, double complex *in, parameters params);
void cn(double complex *out, double complex *in, parameters params);

#endif
