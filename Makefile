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

hydr_harm: hydr_harm.c hamiltonian.c vmath.c geometry.c laplacian.c integrators.c observables.c nfft.c
	gcc hydr_harm.c -lm hamiltonian.c vmath.c geometry.c laplacian.c integrators.c observables.c nfft.c -L/usr/local/lib -lgsl -o hydr_harm.exe -O2
	./hydr_harm.exe

nfft_test: nfft_test.c nfft.c
	gcc nfft_test.c -lm nfft.c -L/usr/local/lib -lgsl -o nfft_test.exe -O2
	./nfft_test.exe

make timed_perf_test: perf_test_time_integrators.c geometry.c hamiltonian.c vmath.c integrators.c observables.c nfft.c laplacian.c
	gcc perf_test_time_integrators.c -lm geometry.c hamiltonian.c vmath.c integrators.c observables.c nfft.c laplacian.c -fopenmp -L/usr/local/lib -lgsl -o perf_test_time_integrators.exe -O2
	./perf_test_time_integrators.exe

make strang_perf_test: perf_test_strang.c geometry.c hamiltonian.c vmath.c integrators.c observables.c nfft.c laplacian.c
	gcc perf_test_strang.c -lm geometry.c hamiltonian.c vmath.c integrators.c observables.c nfft.c laplacian.c -fopenmp -L/usr/local/lib -lgsl -o perf_test_strang.exe -O2
	./perf_test_strang.exe