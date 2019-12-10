indices_test: indices_test.c indices.c
	gcc indices_test.c indices.c -o indices_test.exe -O2 -std=c99
	./indices_test.exe

laplacian_test: laplacian_test.c vmath.c laplacian.c indices.c nneighbour.c
	gcc laplacian_test.c vmath.c laplacian.c indices.c nneighbour.c -o laplacian_test.exe -O2 -std=c99
	./laplacian_test.exe