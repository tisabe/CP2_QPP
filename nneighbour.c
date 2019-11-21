/* Calculate a next neighbour of a reference point given by an index. Specify axis and direction.
Valid values:
    index           0, ..., N**D-1
    axis            0, ..., D
    dir             -1, +1
    N               as specified in the main program
    D               as specified in the main program, typically 1, 2 or 3                       */


/* ########## Still to be tested. ########### */



#include <stdio.h>
#include <stdlib.h>
/* #include "index2coord.h" */
/* #include "coord2index.h" */

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
    if(nneighbour_coords[axis] > N-1)
        nneighbour_coords[axis] = 0;
    if(nneighbour_coords[axis] < 0)
        nneighbour_coords[axis] = N-1;

    nneighbour_index = index2coord(nneighbour_coords);

    free(nneighbour_coords);
    // Return index of next neighbour
    return nneighbour_index;
}
