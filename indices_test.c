/* **** This function tests the index2coord and the coord2index function checking if the coordinates increment and if the sum over one axis equals the checksum (N-1)/2*N**D *** */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "geometry.h"
#include "vmath.h"


// Function that waits for a user input to set the number of points per dimension
long int setN(){
  long int N;
  printf("Set array length. (Recommended N < 1000, enter 0 for random value up to 999): ");
  scanf("%li",&N);

  // Set N to a random number if 0 is selected. Then print the random value.
  if(N==0){
    srand(time(NULL));
    N = rand()%1000;
    printf("\nN = %li\n",N);
  }
  return N;
}


// The program checks every index from 0 to N**D-1 whether its coordinates increment by 1, stay the same or are 0.
// Whenever a coordinate reaches N-1 the next index will have a 0 at this coordinate and the coordinate of the next dimension will be incremented by one.
// As this test is not strictly unambiguous, the checksum is calculated secondly to make sure that every coordinate will be reached exactly once.
// The checksum calculates as follows: If j_i (j=0,...,D-1; i=0,...,N-1) are the coordinates, checksum[j] = sum_{i=0}^{N**D-1} j_i with j = 0,...,D-1.
// For every j it should be checksum[j] = (N-1)/2*N**D
void check_index2coord_inc_sum(long int N){

    // Variable that is changed if any test does not pass
    int not_passed = 0;

    printf("\n\n******* Testing if coordinates increment and if checksum is (N-1)/2*N**D *******\n");

    // Loop over 1, 2, 3 dimensions
    for(unsigned int D=1; D<=3; D++){
        printf("\nTesting %u dimension(s)...",D);

        long int *coords_cur = malloc(D*sizeof(long));
        long int *coords_old = malloc(D*sizeof(long));
        long int *checksum = malloc(D*sizeof(long));

        // Array to store the checksum
        for(int i=0; i<3; i++){
            checksum[i] = 0;
        }

        // Loop over all indices
        for(int i=0; i<ipow(N,D); i++){
            index2coord(coords_cur,i,N,D); //Convert current index to coordinates
            for(int j=0; j<D; j++){ // Loop over all dimensions
                // Check condition if current coordinates are either incremented by 1 with respect to the previous coordinate, the same as the previous coordinate or 0.
                if(!((coords_cur[j] == coords_old[j]+1) || (coords_cur[j] == coords_old[j]) || (coords_cur[j] == 0))){
                    not_passed = 1;
                    printf("Error coordinate: %li",coords_cur[j]);
                }
                coords_old[j] = coords_cur[j];
                checksum[j] += abs(coords_cur[j]);
            }
        }

        // Check whether checksum is equal to the expected value
        for(int j=0; j<D; j++){
            if(checksum[j] != (ipow(N,D)*(N-1)/2)){
                not_passed = 1;
                printf("Error checksum: (%li, %li)\n",checksum[j],ipow(N,D)*(N-1)/2);
            }
        }

        if(not_passed == 0){
          printf("ok");
        }else{
          printf("error");
        }

        free(coords_cur);
        free(coords_old);
        free(checksum);
    }
    if(not_passed == 0){
        printf("\nAll incrementing and sum tests passed!\n");
    }else{
        printf("\nIncrementing test and/ or sum test not passed!\n");
    }
}

// Function that tests if coord2index is the inverse of index2coord.
// It loops over all indices from 0 to N**D-1 and checks if the inversion hold for every single index
void check_inverse(long int N){
  int not_passed = 0;

  printf("\n\n*******Testing if coord2index is inverse of index2coord*******\n");

  // Looping over all dimensions
  for(unsigned int D=1; D<=3; D++){
    printf("\nTesting %u dimension(s)...",D);

    long int index_new;
    long int *coordinates = malloc(D*sizeof(long));

    // Looping over all indices
    for(int i=0; i<ipow(N,D); i++){
      // Transform index to coordinates...
      index2coord(coordinates, i, N, D);
      // ...and back again
      index_new = coord2index(coordinates, N, D);
      // Check if the index after the conversion to coords and back to index matches the previous input index
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

// Main that calls the functions described above
int main(){
    long int N = setN();
    check_index2coord_inc_sum(N);
    printf("\nNext test...");
    check_inverse(N);
    return 0;
}
