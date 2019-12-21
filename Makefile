indices_test: indices_test.c geometry.c vmath.c
	gcc indices_test.c geometry.c vmath.c -o indices_test.exe -O2
	./indices_test.exe

laplacian_test: laplacian_test.c vmath.c laplacian.c geometry.c
	gcc laplacian_test.c vmath.c laplacian.c geometry.c -o laplacian_test.exe -O2
	./laplacian_test.exe
	
cg_test: cg_test.c vmath.c
	gcc -lm cg_test.c vmath.c -o cg_test.exe -O2
	./cg_test.exe

mainfile: main.c euler_method.c vmath.c hermite_polynomial.c geometry.c hamiltonian.c laplacian.c
	gcc -lm main.c euler_method.c vmath.c hermite_polynomial.c geometry.c hamiltonian.c laplacian.c -o main.exe -O2
	./main.exe