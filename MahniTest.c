#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <float.h>
#include <math.h>

#include "write_text.h"

#define _USE_MATH_DEFINES

int main(){

    int L = 5;
	char *blub;

    
    int *x= malloc(L * sizeof(int));
    
    for(int i=0; i<L; i++){
        x[i] = i;
    }
	data2text(x, blub);
}
