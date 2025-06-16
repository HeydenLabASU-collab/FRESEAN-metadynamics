#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/dataTypes.h"
#include "../include/alloc.h"
#include "../include/io.h"
#include "../include/fatal.h"
#include "../include/progress.h"

int printTitle() {
        printf("THE PURPOSE OF THIS PROGRAM IS TO AVERAGE A SET OF\n");
        printf("BINARY FILES, EACH CONTAINING A LIST OF SQUARE MATRICES.\n");
	printf("THE FORMAT IS EXPECTED TO BE THE *.mmat FORMAT USED\n");
	printf("AS OUTPUT OF:\n");
	printf("gen-modes.exe OR gen-modes_omp.exe\n");
        printf("Version 1.0: April 27, 2023\n");
        printf("Author:\n");
        printf(" Michael Sauer and Dr. Matthias Heyden\n");
        printf(" School of Molecular Sciences\n");
        printf(" Arizona State University\n");
        printf(" Tempe, AZ, USA\n");
        printf(" e-mail: masauer2@asu.edu/mheyden1@asu.edu\n");
        printf("\n");
        return 0;
}

int getLineFromCOM(FILE *in,char *buffer,int max)
{
        if(fgets(buffer,max,in)==NULL) fatal("input file incomplete\n");
        while(strncmp(buffer,"#",1)==0)
        {
                if(fgets(buffer,max,in)==NULL) fatal("input file incomplete\n");
        }
        return 0;
}

int printKeys() {
        printf("nFiles (number of files to average over)\nfnList (list of individual matrix files, newline delimited)\nfnOut (average matrix output file)\n");
        return 0;
}

int getInput(char *fnCOM, int *nFiles,char ***fnList,char *fnOut) {
        FILE *io;
        char buffer[300];
        int format;
        int i=1;
        int j;
        float tmp;

        saveOpenRead(&io,fnCOM);

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nFiles);
        printf("%2d -> read %20s : %d\n",i,"nFiles",nFiles[0]);i++;

        fnList[0]=(char**)save_malloc(nFiles[0]*sizeof(char*));
        for(j=0;j<nFiles[0];j++) {
                fnList[0][j]=(char*)save_malloc(300*sizeof(char));
                getLineFromCOM(io,buffer,300);
                sscanf(buffer,"%s",fnList[0][j]);
                printf("%2d -> read %17s[%d] : %s\n",i,"fnList",j,fnList[0][j]);
                i++;
        }

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnOut);
        printf("%2d -> read %20s : %s\n",i,"fnOut",fnOut);i++;

        fclose(io);
        return 0;
}

typedef struct {
        double *m;
        int nu,nv;
} t_mat;

int getMat(t_mat mat,int u,int v,double *muv) {
        int i;
        i=u*mat.nv+v;
        muv[0]=mat.m[i];
        return 0;
}

int setMat(t_mat mat,int u,int v,double muv) {
        int i;
        i=u*mat.nv+v;
        mat.m[i]=muv;
        return 0;
}

int main(int argc, char *argv[])
{
        FILE  *tin,*cin,*vin,*job,*out;
        XDR xdrin;
        char fnCOM[300],fnOut[300];
        char **fnList;
        int nFiles=1;
        int magic1,magic2;
        t_mat *mat;
        int i,j,k;
        int nDOF,nDOFtest;;
        int nMatPerFile;
        int u,v;
        double muu,muv,mvu,mvv;
        FILE *in;
        t_mat *matrixList;
	double tmp;

        printTitle();
        if(argc!=2) {
                printKeys();
                fatal("no input file specified\n");
        }

        strcpy(fnCOM,argv[1]);
        getInput(fnCOM,&nFiles,&fnList,fnOut);
        saveOpenReadBin(&in,fnList[0]);
        fread(&magic1,sizeof(int),1,in);
        fread(&nMatPerFile,sizeof(int),1,in);
        printf("%d %d\n", magic1, nMatPerFile);
	fclose(in);
        printf("%s\n", fnList[0]); 
	matrixList=(t_mat*)save_malloc(nMatPerFile*sizeof(t_mat));
        
	
	for(i=0;i<nFiles;i++) {
                saveOpenReadBin(&in,fnList[i]);
                printf(" - reading file: %s\n",fnList[i]);
		fflush(stdout);
		
		printf("Made it 1");
                // Reading in MMAT file (list of square matrices)
                fread(&magic1,sizeof(int),1,in); 
                fread(&nMatPerFile,sizeof(int),1,in);
                fread(&nDOF,sizeof(int),1,in);
                fread(&nDOFtest,sizeof(int),1,in);
		if(nDOF!=nDOFtest) fatal("encountered non square matrix\n");
                fread(&magic1,sizeof(int),1,in);
		printf("Made it 0");
                for(j=0;j<nMatPerFile;j++) {
                        if(i==0){
                                matrixList[j].m=(double*)save_malloc(nDOF*nDOF*sizeof(double));
				for(k=0;k<nDOF*nDOF;k++) {
					matrixList[j].m[k]=0.0;
				}
                        }
                        matrixList[j].nu=nDOF;
                        matrixList[j].nv=nDOF;
                }
                printf("Made it 1");
                for(j=0;j<nMatPerFile;j++) {
                        fread(&magic2,sizeof(int),1,in);
                        for(k=0;k<nDOF*nDOF;k++) {
                                fread(&tmp,sizeof(double),1,in);
				matrixList[j].m[k]+=tmp;
                        }
                        fread(&magic2,sizeof(int),1,in);
//progressBar(j,1,nMatPerFile);
                }
		printf("Made it 2");
                fclose(in);
        }
        
        for(j=0;j<nMatPerFile;j++) {
                for(k=0;k<nDOF*nDOF;k++) {
                        matrixList[j].m[k]/=nFiles;
                }
        }
                
        printf(" - writing output file: %s\n",fnOut);fflush(stdout);
        magic1=3*sizeof(int);
        magic2=nDOF*nDOF*sizeof(double);
        out=fopen(fnOut,"wb");
        fwrite(&magic1,sizeof(int),1,out);
        fwrite(&nMatPerFile,sizeof(int),1,out);
        fwrite(&nDOF,sizeof(int),1,out);
        fwrite(&nDOF,sizeof(int),1,out);
        fwrite(&magic1,sizeof(int),1,out);
        for(j=0;j<nMatPerFile;j++) {
                fwrite(&magic2,sizeof(int),1,out);
                for(k=0;k<nDOF*nDOF;k++) {
                        fwrite(&matrixList[j].m[k],sizeof(double),1,out);
                }
                fwrite(&magic2,sizeof(int),1,out);
		progressBar(j,1,nMatPerFile);
        }
        fclose(out);
        return 0;
}

