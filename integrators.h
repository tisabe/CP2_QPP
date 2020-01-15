#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include <complex.h>
#include "structs.h"

void step_euler(double complex *out, double complex *in, parameters params);
void step_euler_timed(double complex *out, double complex *in, parameters params);
void step_cn(double complex *out, double complex *in, parameters params);
void step_cn_timed(double complex *out, double complex *in, parameters params);
void step_strang(double complex *out, double complex *in, parameters params);;

#endif
