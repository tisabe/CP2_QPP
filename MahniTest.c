#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <float.h>
#include <math.h>

#include "write_text.h"

#define _USE_MATH_DEFINES

int main(){

    int L = 5

    
    double complex *x= malloc(L * sizeof(double complex));
    
    for(long int i=0; i<L; i++){
        x[i] = i;
    }
	data2text(x, char *blub);
}
