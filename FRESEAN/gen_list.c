#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

int printTitle() {
		printf("THE PURPOSE OF THIS PROGRAM IS TO GENERATE \n");
		printf("A CUSTOM TOPOLOGY FROM CP2K POSITION DATA. \n");
        printf("Author:\n");
        printf(" Michael Sauer and Dr. Matthias Heyden\n");
        printf(" School of Molecular Sciences\n");
        printf(" Arizona State University\n");
        printf(" Tempe, AZ, USA\n");
        printf(" e-mail: masauer2@asu.edu/mheyden1@asu.edu\n");
        printf("\n");
        return 0;
}

int main(int argc, char* argv[]){
    if(argc != 4){
        printf("usage: ./gen_list nCorr timestep (ps) output_file\n");
		exit(0);
	}

    int nCorr;
    float dt, df;
    char fnOut[100];
    FILE *out;
	sscanf(argv[1],"%d",&nCorr);
	sscanf(argv[2],"%f",&dt);
    sscanf(argv[3],"%s",&fnOut);
    out = fopen(fnOut, "w");
    df = 33.3567/((2*nCorr-1)*dt);
    printf("Frequency Resolution: %f cm^-1\n", df);
    printf("Generating list of frequencies at %s\n", fnOut);
    for(int i = 0; i < nCorr; i++){
        fprintf(out, "%d %f\n", i, df*i);
    }
    fclose(out);
    printf("Closed output file.");
}