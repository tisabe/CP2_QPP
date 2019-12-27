#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "structs.h"
#include "vmath.h"
#include "potentials.h"
#include "geometry.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define _USE_MATH_DEFINES

#include "strang_splitting.h"
#include "euler_method.h"
#include "crank_nicolson.h"

//create the data "fdata" with one of the three methods
.
.
.

// write the created data into a .txt file
// https://www.cprogramming.com/tutorial/cfileio.html

FILE *f= fopen("data.txt", "w");
written=fwrite(fdata, sizeof(fdata[0]), sizeof(fdata), f);
// Test if fopen fails
if (written == 0) {
    printf("Error during writing to file !");
}
fclose(f);
