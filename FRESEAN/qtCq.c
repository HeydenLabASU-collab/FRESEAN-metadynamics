#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "../include/dataTypes.h"
#include "../include/alloc.h"
#include "../include/fatal.h"
#include "../include/io.h"

int printTitle() {
    printf("THE PURPOSE OF THIS PROGRAM IS TO COMPUTE THE\n");
    printf("PRODUCT Q^T.C(x).Q AS A FUNCTION OF x, WHERE Q\n");
    printf("IS A MATRIX OF ORTHONORMAL VECTORS (AS COLUMNS),\n");
    printf("Q^T IS ITS TRANSPOSED FORM AND C(x) IS A SET OF\n");
    printf("SQUARE MATRICES THAT VARY WITH x.\n");
    printf("Version 1.0: Aug. 20, 2025\n");
    printf("Author:\n");
    printf(" Dr. Matthias Heyden\n");
    printf(" School of Molecular Sciences\n");
    printf(" Arizona State University\n");
    printf(" Tempe, AZ, USA\n");
    printf(" e-mail: mheyden1@asu.edu\n");
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
    printf("fnC (*.mmat-format)\n");
    printf("fnQ^T (*.mmat-format)\n");
    printf("selQset\nfirst\nlast\nx0\ndx\nfnOut\n");
    return 0;
}

int getInput(char *fnCOM,char *fnCx,char *fnQt,int *selQset,
                int *first,int *last,double *x0,double *dx,
                char *fnOut) {
    FILE *io;
    char buffer[300];
    int format;
    int i=1;
    float tmp;

    saveOpenRead(&io,fnCOM);

    getLineFromCOM(io,buffer,300);
    sscanf(buffer,"%s",fnCx);
    printf("%2d -> read %20s : %s\n",i,"fnCx",fnCx);i++;

    getLineFromCOM(io,buffer,300);
    sscanf(buffer,"%s",fnQt);
    printf("%2d -> read %20s : %s\n",i,"fnQt",fnQt);i++;
    
    getLineFromCOM(io,buffer,300);
    sscanf(buffer,"%d",selQset);
    printf("%2d -> read %20s : %d\n",i,"selQset",selQset[0]);i++;

    getLineFromCOM(io,buffer,300);
    sscanf(buffer,"%d",first);
    printf("%2d -> read %20s : %d\n",i,"first",first[0]);i++;

    getLineFromCOM(io,buffer,300);
    sscanf(buffer,"%d",last);
    printf("%2d -> read %20s : %d\n",i,"last",last[0]);i++;

    getLineFromCOM(io,buffer,300);
    sscanf(buffer,"%f",&tmp);
    x0[0]=tmp;
    printf("%2d -> read %20s : %f\n",i,"x0",tmp);i++;

    getLineFromCOM(io,buffer,300);
    sscanf(buffer,"%f",&tmp);
    dx[0]=tmp;
    printf("%2d -> read %20s : %f\n",i,"dx",tmp);i++;

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
    FILE  *inCx,*inQt,*out;
    char fnCOM[300],fnCx[300],fnQt[300],fnOut[300],fnOut2[300];;
    int i,j,k,l;
    int selQset,first,last;
    double x0,dx;
    int nx;
    t_mat Qt;
    t_mat Cx;
    double p;
    int n1,n2,n3,n4,bc1,bc2;
    int n;
    double tmp1;
    double a,b,c;
    double **res;

    printTitle();

    /*PARSE COMMAND LINE INPUT HERE*/
    if(argc!=2) {
        printKeys();
        fatal("no input file specified\n");
    }
    strcpy(fnCOM,argv[1]);
    getInput(fnCOM,fnCx,fnQt,&selQset,&first,&last,&x0,&dx,fnOut);

    n=last-first+1;
    if(n<1) fatal("first must be smaller than last\n");

    inQt=fopen(fnQt,"rb");
    fread(&bc1,sizeof(int),1,inQt);
    fread(&nx,sizeof(int),1,inQt);
    fread(&n1,sizeof(int),1,inQt);
    fread(&n2,sizeof(int),1,inQt);
    fread(&bc2,sizeof(int),1,inQt);
    if(selQset>nx) {
        printf("%s only has %d sets, set selQset=%d not available\n",fnQt,nx,selQset);
        exit(1);
    }
    if(bc1!=bc2) {
        printf("format error in %s\n",fnQt);
        exit(1);
    }
    if(last>n1) {
        printf("%s does not have %d vectors\n",fnQt,last);
        exit(1);
    }
    Qt.m=(double*)save_malloc(n1*n2*sizeof(double));
    Qt.nu=n1;
    Qt.nv=n2;

    for(i=0;i<selQset-1;i++) {
        fread(&bc1,sizeof(int),1,inQt);
        for(j=0;j<n1*n2;j++) {
            fread(&tmp1,sizeof(double),1,inQt);
        }
        fread(&bc2,sizeof(int),1,inQt);
        if(bc1!=bc2) {
            printf("format error in %s\n",fnQt);
            exit(1);
        }
    }
    for(i=0;i<1;i++) {
        fread(&bc1,sizeof(int),1,inQt);
        for(j=0;j<n1*n2;j++) {
            fread(&tmp1,sizeof(double),1,inQt);
            Qt.m[j]=tmp1;
        }
        fread(&bc2,sizeof(int),1,inQt);
        if(bc1!=bc2) {
            printf("format error in %s\n",fnQt);
            exit(1);
        }
    }
    fclose(inQt);

    inCx=fopen(fnCx,"rb");
    fread(&bc1,sizeof(int),1,inCx);
    fread(&nx,sizeof(int),1,inCx);
    fread(&n3,sizeof(int),1,inCx);
    fread(&n4,sizeof(int),1,inCx);
    fread(&bc2,sizeof(int),1,inCx);
    if(bc1!=bc2) {
        printf("format error in %s\n",fnCx);
        exit(1);
    }
    if(n3!=n4) {
        printf("incompatible dimensions in %s: n1 %d n2 %d\n",fnCx,n3,n4);
        exit(1);
    }
    if(n2!=n3) {
        printf("incompatible dimensions in %s and %s: Q^T = %d x %d and C(x) = %d x %d\n",fnQt,fnCx,n1,n2,n3,n4);
        exit(1);
    }
    Cx.m=(double*)save_malloc(n3*n3*sizeof(double));
    Cx.nu=n3;
    Cx.nv=n3;

    res=(double**)save_malloc(nx*sizeof(double*));
    for(i=0;i<nx;i++) {
        res[i]=(double*)save_malloc(n*sizeof(double));
        for(j=0;j<n;j++) {
            res[i][j]=0.0;
        }
    }

    for(i=0;i<nx;i++) {
        fread(&bc1,sizeof(int),1,inCx);
        for(j=0;j<n3*n3;j++) {
            fread(&tmp1,sizeof(double),1,inCx);
            Cx.m[j]=tmp1;
        }
        fread(&bc2,sizeof(int),1,inCx);
        if(bc1!=bc2) {
            printf("format error in %s\n",fnCx);
            exit(1);
        }
        #pragma omp parallel for private(k,l,p,a,b,c)
        for(j=first-1;j<last;j++) {
            for(k=0;k<n3;k++) {
                p=0.0;
                for(l=0;l<n3;l++) {
                    getMat(Cx,k,l,&a);
                    getMat(Qt,j,l,&b);
                    p+=a*b;
                }
                getMat(Qt,j,k,&c);
                p*=c;
                res[i][j]+=p;
            }
        }
	printf("%c%4d/%4d",(char)13,i+1,nx);fflush(stdout);
    }
    printf("\n");fflush(stdout);
    fclose(inCx);

    sprintf(fnOut2,"QtCQ_Qset%d_%d-%d_%s.dat",selQset,first,last,fnOut);
    printf(" - writing output in file: %s\n",fnOut2);
    out=fopen(fnOut2,"w");
    /*printing times*/
    for(j=0;j<nx;j++) {
            fprintf(out," %e",x0+j*dx);
    }
    fprintf(out,"\n");
    for(i=first-1;i<last;i++) {
            for(j=0;j<nx;j++) {
                    fprintf(out," %e",res[j][i]);
            }
            fprintf(out,"\n");
    }
    fclose(out);
    printf(" - closed file: %s\n",fnOut2);

    return 0;
}

