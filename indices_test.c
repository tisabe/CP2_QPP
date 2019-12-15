#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "indices.h"
#include "vmath.h"

long int setN(){
  long int N;
  printf("Set array length. Must be uneven. (Recommended N < 1000, 0 for random value up to 999): ");
  scanf("%li",&N);
  // Set N to a random uneven number if 0 is selected. Then print the random value.
  if(N==0){
    while(N%2 == 0){
      srand(time(NULL));
      N = rand()%1000;
    }
    printf("N = %li\n",N);
  }
  return N;
}

void check_index2coord_inc_sum(long int N){
    int not_passed = 0;

    printf("\n\n*******Testing if coordinates increment and if sum over each axis is 0*******\n");

    for(unsigned int D=1; D<=3; D++){
        printf("\nTesting %u dimension(s)...",D);

        long int *coords_cur = malloc(D*sizeof(long));
        long int *coords_old = malloc(D*sizeof(long));
        long int *sum = malloc(D*sizeof(long));

        for(int i=0; i<3; i++){
            sum[i] = 0;
        }

        for(int i=0; i<ipow(N,D); i++){
            index2coord(coords_cur,i,N,D);
            for(int j=0; j<D; j++){
                if(!((coords_cur[j] == coords_old[j]+1) || (coords_cur[j] == coords_old[j]) || (coords_cur[j] == -(N-1)/2))){
                    not_passed = 1;
                    printf("Error coordinate: %li",coords_cur[j]);
                }
                coords_old[j] = coords_cur[j];
                sum[j] += coords_cur[j];
            }
        }

        for(int j=0; j<D; j++){
            if(sum[j] != 0){
                not_passed = 1;
                printf("Error sum: (%li, %li)\n",sum[j],0);
            }
        }

        if(not_passed == 0){
          printf("ok");
        }else{
          printf("error");
        }

        free(coords_cur);
        free(coords_old);
        free(sum);
    }
    if(not_passed == 0){
        printf("\nAll incrementing and sum tests passed!\n");
    }else{
        printf("\nIncrementing test and/ or sum test not passed!\n");
    }
}

void check_inverse(long int N){
  int not_passed = 0;

  printf("\n\n*******Testing if coord2index is inverse of index2coord*******\n");

  for(unsigned int D=1; D<=3; D++){
    printf("\nTesting %u dimension(s)...",D);

    long int index_new;
    long int *coordinates = malloc(D*sizeof(long));

    for(int i=0; i<ipow(N,D); i++){
      index2coord(coordinates, i, N, D);
      index_new = coord2index(coordinates, N, D);
      if(index_new != i){
        not_passed = 1;
        printf("Error index: %i\n",i);
      }
    }

    if(not_passed == 0){
      printf("ok");
    }else{
      printf("error");
    }

    free(coordinates);
  }
  if(not_passed == 0){
    printf("\nAll inversion tests passed!\n");
  }else{
    printf("\nInversion test not passed!\n");
  }
}

int main(int argc, char *argv[]){
    long int N = setN();
    check_index2coord_inc_sum(N);
    printf("\nNext test...");
    check_inverse(N);
    return 0;
}
