#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "indices.h"

// No real test routines implemented yet. Just looping over all indices and printing their coordinates.

void check_index2coord(){
    for(unsigned int D=1; D<=3; D++){
        srand(time(NULL));
        long int N = rand()%100;
        printf("%li\n\n",N);

        long int *coords = malloc(D*sizeof(long));

        for(int i=0; i<ipow(N,D); i++){
            index2coord(coords,i,N,D);
            for(int j=0; j<D; j++){
                printf("%li ",coords[j]);
            }
            printf("\n");
        }
        printf("\n\n");

        free(coords);
    }
    printf("Done!");
}

int main(int argc, char *argv[]){
    check_index2coord();
    return 0;
}
