#include <stdio.h>
#include <stdlib.h>

#include "structs.h"
#include "vmath.h"

long int coord2index(long int *coord, long int N, unsigned int D){
  long int index = 0;
  long int shift = (N-1)/2; // coordinate shift amount to have coordinate origin at middle of axis

  if (((N+1)%2) != 0) {
    printf("Error: Length of axis N should be uneven");
    return -1;
  }

  for (int i=0; i<D; i++){
    index += (coord[i]+shift)*ipow(N,i);
  }

  return index;
}

void index2coord(long int *coord, long int index, long int N, unsigned int D){
  /*coordinate vector needs to be returned as pointer in C*/
  //static int coord[D];
  if (((N+1)%2) != 0) {
    printf("Error: Length of axis N should be uneven");
  }

  long int b = 0;
  long int shift = (N-1)/2; // coordinate shift amount to have coordinate origin at middle of axis

  for (int i=0; i<D; i++) {
    b = ipow(N, D-1-i);
    coord[D-1-i] = index/b - shift;
    index %= b;
  }
}


/* Calculate a next neighbour of a reference point given by an index. Specify axis and direction.
Valid values:
    index           0, ..., N**D-1
    axis            0, ..., D
    dir             -1, +1
    N               as specified in the main program
    D               as specified in the main program, typically 1, 2 or 3                       */


/* ########## TESTED BY TIM - SUCCESS - SEE nneighbour_test.c ########### */

long int nneighbour(long int index, unsigned int axis, int dir, long int N, unsigned int D)
{
    /* Initialise pointer for next neighbour coords */
    long int nneighbour_index;
    long int *nneighbour_coords;
    nneighbour_coords = malloc(D*sizeof(long));

    //Get coords of reference point. In-/Decrement by 1 in axis-direction
    index2coord(nneighbour_coords,index,N,D);
    nneighbour_coords[axis] += dir;

    //Check for boundary conditions
    if(nneighbour_coords[axis] > (N-1)/2)               //Edits because of change of coordinate origin
        nneighbour_coords[axis] = -(N-1)/2;
    if(nneighbour_coords[axis] < -(N-1)/2)
        nneighbour_coords[axis] = (N-1)/2;

    nneighbour_index = coord2index(nneighbour_coords, N, D);

    free(nneighbour_coords);
    // Return index of next neighbour
    return nneighbour_index;
}

//Call nneihbour_init to create an array that contains all the next neighbours for each point in consecutive order

//long int *out array has to be of length 2*D*L (i.e. long int *out = malloc(2*D*L*sizeof(long int)) )

void nneighbour_init(long int *out, long int N, unsigned int D){
  for(int i=0; i<ipow(N,D); i++){
    for(int j=0; j<D; j++){
      out[i*2*D+2*j] = nneighbour(i,j,-1,N,D);
      out[i*2*D+2*j+1] = nneighbour(i,j,1,N,D);
    }
  }
}
