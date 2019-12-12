indices_test: indices_test.c indices.c
	gcc indices_test.c indices.c -o indices_test.exe -O2
	./indices_test.exe

laplacian_test: laplacian_test.c vmath.c laplacian.c indices.c nneighbour.c
	gcc laplacian_test.c vmath.c laplacian.c indices.c nneighbour.c -o laplacian_test.exe -O2
	./laplacian_test.exe
	
cg_test: cg_test.c vmath.c
	gcc cg_test.c vmath.c -o cg_test.exe -O2
	./cg_test.exe

mainfile: main.c euler_method.c vmath.c hermite_polynomial.c indices.c hamiltonian.c
	gcc -lm main.c euler_method.c vmath.c hermite_polynomial.c indices.c hamiltonian.c -o main.exe
	./main.exe
