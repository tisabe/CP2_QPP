SOLVER FOR THE TIME-DEPENDENT SCHROEDINGER EQUATION


Please use the file 'hydr_harm.c' to generate executables that simulate the time evolution of a wave function in the harmonic potential between two hydrogen atoms.
Precompiled executables can be found in the folder 'executables'. For compiling new executables please use the 'Makefile'.
In the folder 'sample_outputs' you can find outputs which were generated using the hydr_harm_crank_nicolson with the parameters as shown in the readme in that folder.


The functions are structured the following way:

geometry.c:
	index2coord		calculates the coordinates that correspond to the index of a lattice point
	coord2index		calculates the index that corresponds to the coordinates of a lattice point
	nneighbour		calculates the indices of next neighbours of a point (index given) in a specified dimension and direction
	nneighbour_init		used to create an array that contains all the next neighbours

laplacian.c:
	laplacian		uses the array of next neighbours to calculate the symmetrical discrete laplacian of a given point (index)

vmath.c:
	miscancellous functions that perform mathematical operations (esp. vector operations)
	cg			implementation of the conjugate gradient algorithm that solves for the vector x in a matrix-vector equation Ax=b

hamiltonian.c:
	hamiltonian		calculates the Hamiltonian of an input wave function
	gen_pot_harmonic	used to generate a harmonic potential
	gen_pot_box		used to generate a quantum well potential

integrators.c:
	step_euler		calculates one time step with Euler method time integration
	step_cn			calculates one time step with Crank-Nicolson time integration
	step_strang		calculates one time step with Strang-splitting time integration

structs.h:
	defines a set of parameters that are given to function to make parameter passing less confusing

observables.c:
	obs_norm		calculates the normalisation of a wave function
	obs_E			calculates the energy expectation value of a wave function
	obs_x			calculates the spacial expectation value of a wave function in a specified dimension
	obs_p			calculates the momentum expectation value of a wave function in a specified dimension
	obs_delta_x		calculates the spacial width of a wave function
	obs_delta_p		calculates the momentum width of a wave function
