#include <stdio.h>
#include <stdbool.h>

#include "structs.h"
#include "geometry.h"
#include "vmath.h"

void test_inverse_nn() {
	/* This function tests the invertibility of the next neighbour function,
	by taking the neighbour in one direction and the neighbour of that neighbour
	in the opposite direction, which should result in the index we started with.
	This is tested for a couple of dimensions and axis lengths and every time
	for all directions that are possible (D times).
	The result of the test will be printed out. */
	printf("Testing inversion of next neighbour function...\n");

	unsigned int arrD[4] = {1,2,3,4};
	long int arrN[2] = {5,11};
	long int maxIndex = 0;
	bool passed = true;

	for (int i=0; i<4; i++) {							// loop over dimensions
		printf("Testing D=%d ", arrD[i]);
		for (int j=0; j<2; j++) {						// loop over coordinates
			printf("Testing N=%d ", arrN[j]);
			maxIndex = ipow(arrN[j], arrD[i]);
			for (int ind=0; ind<maxIndex; ind++) {		// loop over indices
				for (int d=0; d<arrD[i]; d++){			// loop over axis
					long int nni = nneighbour(ind, d, 1, arrN[j], arrD[i]);
					long int inni = nneighbour(nni, d, -1, arrN[j], arrD[i]);
					if (ind != inni) {
						passed = false;
						printf("\n Error at index %d, nni: %d, inni: %d ", ind, nni, inni);
					}
				}

			}
		}
		printf("\n");
	}
	if (passed) {
		printf("Test passed!");
	} else {
		printf("Test not passed!");
	}
}

int main() {
	long int N = 11;
	unsigned int D = 2;
	long int nni = 0;
	long int nnii = 0;
	long int coord[D];

	/* For one specific example the indices and coordinates are printed out
	to check for inconsistencies. */
	printf("Testing the nneighbour function with N=%ld, D=%d\n", N, D);
	printf("i:, \t coord: \t nni: \t nnii: \n");

	for (int i=0; i<20; i++){
		nni = nneighbour(i, 0, 1, N, D);
		nnii = nneighbour(nni, 0, -1, N, D);
		index2coord(coord, i, N, D);
		printf("%d \t (%d,%d) \t \t %d \t %d \n", i, coord[0], coord[1], nni, nnii);
	}
	
	test_inverse_nn();
	return 0;
}
