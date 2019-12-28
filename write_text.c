#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


// function which writes an array "fdata" into a .txt file named "data"
	// https://www.cprogramming.com/tutorial/cfileio.html

void data2text(fdata) {
FILE *f= fopen("data.txt", "w");
written=fwrite(fdata, sizeof(fdata[0]), sizeof(fdata), f);
// Test if fopen fails
	if (written == 0) {
    printf("Error during writing to file !");
	}
	/* Alternative
	for(long int i=0; i < sizeof(fdata); i++){
		fprintf(f, "fdata[i] \n");
	} */
fclose(f);
}