#ifndef GEOMETRY_H
#define GEOMETRY_H

long int coord2index(long int *coord, long int N, unsigned int D);
void index2coord(long int *coord, long int index, long int N, unsigned int D);
long int nneighbour(long int index, unsigned int axis, int dir, long int N, unsigned int D);
void nneighbour_init(long int *out, long int N, unsigned int D);

#endif
