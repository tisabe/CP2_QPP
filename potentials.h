#ifndef POTENTIALS_H
#define POTENTIALS_H

void well(double *wellpotential, parameters params, double height);
void box(double *boxpotential, parameters params, double height);
void harmonic(double *harmonicpotential, parameters params, double omega);

#endif
