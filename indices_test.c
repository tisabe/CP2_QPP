#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "indices.h"

// No real test routines implemented yet. Just looping over all indices and printing their coordinates.

void check_index2coord(){
    int not_passed = 0;
    long int N;
    printf("Set array length (recommended <= 1000, 0 for random value up to 1000): ");
    scanf("%li",&N);
    if(N==0){
        srand(time(NULL));
        N = (rand()/RAND_MAX)%1000;
        printf("N = %li\n",N);
    }

    for(unsigned int D=1; D<=3; D++){
        printf("Testing %u dimension(s)...\n",D);

        long int *coords_cur = malloc(D*sizeof(long));
        long int *coords_old = malloc(D*sizeof(long));
        long int *sum = malloc(D*sizeof(long));

        for(int i=0; i<3; i++){
            sum[i] = 0;
        }

        for(int i=0; i<ipow(N,D); i++){
            index2coord(coords_cur,i,N,D);
            for(int j=0; j<D; j++){
                if(!((coords_cur[j] >= coords_old[j]) || (coords_cur[j] == 0))){
                    not_passed = 1;
                }
                coords_old[j] = coords_cur[j];
                sum[j] += coords_cur[j];
            }
        }

        for(int j=0; j<D; j++){
            if(sum[j] != ipow(N,D)*(N-1)/2){
                not_passed = 1;
                printf("(%li, %li)\n",sum[j],ipow(N,D)*(N-1)/2);
            }
        }

        free(coords_cur);
        free(coords_old);
        free(sum);
    }
    if(not_passed == 0){
        printf("All incrementing and sum tests passed!");
    }else{
        printf("Incrementing test and/ or sum test not passed!");
    }
}

int main(int argc, char *argv[]){
    check_index2coord();
    return 0;
}
