#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

void hamiltonian(double complex *out, double complex *in, parameters params);
void hamiltonian_timed(double complex *out, double complex *in, parameters params);
void hamiltonian_slow(double complex *out, double complex *in, parameters params);
void hamiltonian_single_loop(double complex *out, double complex *in, parameters params);
void hamiltonian_parallel(double complex *out, double complex *in, parameters params);
void gen_pot_harmonic(parameters *params, double omega);
void gen_pot_box(parameters *params, double height);
void gen_pot_well(parameters *params, double height);

#endif
