#include <stdio.h>
#include <stdbool.h>
#include "nneighbour.h"
#include "indices.h" //for ipow, maybe seperate functions later

void test_inverse_nn() {
	printf("Testing inversion of next neighbour function\n");

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
	printf("Testing the nneighbour function\n");
	printf("i:, \t coord: \t nni: \t nnii: \n");

	long int N = 11;
	unsigned int D = 2;
	long int nni = 0;
	long int nnii = 0;
	long int coord[D];

	for (int i=0; i<20; i++){
		nni = nneighbour(i, 0, 1, N, D);
		nnii = nneighbour(nni, 0, -1, N, D);
		index2coord(coord, i, N, D);
		printf("%d \t (%d,%d) \t \t %d \t %d \n", i, coord[0], coord[1], nni, nnii);
	}
	test_inverse_nn();
	return 0;
}
