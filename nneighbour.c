/* Lots of work to be done. Not running yet. */






#include <stdio.h>
/* #include "index2coord.h" */
/* #include "coord2index.h" */

void addArrays(int *out, int *first, int *second){
    for(int i=0; i<3; i++){
        out[i] = first[i] + second[i];
    }
}

unsigned int nneighbour(unsigned int neighbour_index, unsigned int index, unsigned int dim, int dir, unsigned int N, unsigned int D)
{
    /* Initialise integer neighbour */
    int neighbour;

    int dir_array[D] = {0};                                     // Create an array of zeros of the number of dimensions
    dir_array[dim] = dir;                                       // Change the entry of the value specified by dim into dir (-1 or +1)


    //Check for boundary conditions here

    /* Convert to coord, add dir_array and then convert back to index. */
    neighbour = index2coord(result,addArrays(coord2index(index),dir_array));
    if(dim==0){
        neighbour = coord2index(index2coord(index) + [dir,0,0])
    }
    if(dim==1){
        neighbour = coord2index(index2coord(index) + [0,dir,0])
    }
    if(dim==2){
        neighbour = coord2index(index2coord(index) + [0,0,dir])
    }else{
        println("Limited to 3 dimensions!")
    }

    return neighbour_index;
}
