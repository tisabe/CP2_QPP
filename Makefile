indices_test: indices_test.c geometry.c vmath.c
	gcc indices_test.c -lm geometry.c vmath.c -o indices_test.exe -O2
	./indices_test.exe

laplacian_test: laplacian_test.c vmath.c laplacian.c geometry.c
	gcc laplacian_test.c vmath.c laplacian.c geometry.c -o laplacian_test.exe -O2
	./laplacian_test.exe
	
nneighbour_test: nneighbour_test.c geometry.c vmath.c
	gcc nneighbour_test.c geometry.c vmath.c -lm -o nneighbour_test.exe -O2
	./nneighbour_test.exe
	
cg_test: cg_test.c vmath.c
	gcc -lm cg_test.c vmath.c -o cg_test.exe -O2
	./cg_test.exe

mainfile: main.c integrators.c vmath.c hermite_polynomial.c geometry.c hamiltonian.c laplacian.c
	gcc -lm main.c integrators.c vmath.c hermite_polynomial.c geometry.c hamiltonian.c laplacian.c -o main.exe -O2
	./main.exe

hydr_harm: hydr_harm.c hamiltonian.c vmath.c geometry.c laplacian.c integrators.c observables.c nfft.c
	gcc hydr_harm.c -lm hamiltonian.c vmath.c geometry.c laplacian.c integrators.c observables.c nfft.c -L/usr/local/lib -lgsl -o hydr_harm.exe -O2
	./hydr_harm.exe

structs_test: structs_test.c vmath.c
	gcc structs_test.c vmath.c -lm -o structs_test.exe -O2
	./structs_test.exe