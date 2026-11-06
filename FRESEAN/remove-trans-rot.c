#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "../include/dataTypes.h"
#include "../include/alloc.h"
#include "../include/top.h"
#include "../include/io.h"
#include "../include/mol.h"
#include "../include/fatal.h"
#include "../include/job.h"
#include "../include/select.h"
#include "../include/geo.h"
#include "../include/mol.h"
#include "../include/align.h"
#include "../include/progress.h"

int printTitle() {
        printf("THE PURPOSE OF THIS PROGRAM IS TO REMOVE\n");
        printf("TRANSLATIONAL AND ROTATIONAL COMPONENTS FROM\n");
        printf("A FREQUENCY-DEPENDENT CROSS-CORRELATION MATRIX\n");
        printf("OF SQRT(M)-WEIGHTED VELOCITIES.\n");
        printf("Version 1.0: November 6, 2026\n");
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
        printf("fnTop (topology file [.mtop])\n");
        printf("fnRef (reference coordinate file)\nfnMatrix (cross-correlation matrix file)\n");
        printf("nCorr (number of points in time correlation function)\nbinaryMatrix (0 for human readable, 1 for binary matrix)\n");
        printf("fnOut (output matrix file)\n");
        return 0;
}

int getInput(char *fnCOM,char *fnTop,
                char *fnRef,char *fnMatrix,
                int *nCorr,int *binaryMatrix,
                char *fnOut) {
        FILE *io;
        char buffer[300];
        int format;
        int i=1;
        float tmp;

        saveOpenRead(&io,fnCOM);

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnTop);
        printf("%2d -> read %20s : %s\n",i,"fnTop",fnTop);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnRef);
        printf("%2d -> read %20s : %s\n",i,"fnRef",fnRef);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnMatrix);
        printf("%2d -> read %20s : %s\n",i,"fnMatrix",fnMatrix);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nCorr);
        printf("%2d -> read %20s : %d\n",i,"nCorr",nCorr[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",binaryMatrix);
        printf("%2d -> read %20s : %d\n",i,"binaryMatrix",binaryMatrix[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnOut);
        printf("%2d -> read %20s : %s\n",i,"fnOut",fnOut);i++;

        fclose(io);
        return 0;
}

int getRefTransRot(char *fnRef,int nAtoms,t_atom *atoms,double **evTrRot)
{
        int i;
        FILE *ref;
        char buffer[200];
        int nRef;
        float xr,yr,zr;
        real *masses;
        real totMass=0.0;

        ref=fopen(fnRef,"r");
        fgets(buffer,200,ref);
        fgets(buffer,200,ref);
        sscanf(buffer,"%d",&nRef);
        if(nRef!=nAtoms) {
                sprintf(buffer,"nAtoms in %s (%d) and matrix (%d) don't match!\n",fnRef,nRef,nAtoms);
                fatal(buffer);
        }
        allocReals(masses,nAtoms);
        refCrd[0]=(rvec*)save_malloc(nAtoms*sizeof(rvec));
        refCOM[0].x=0.0; refCOM[0].y=0.0; refCOM[0].z=0.0;
        for(i=0;i<nRef;i++) {
                fgets(buffer,200,ref);
                sscanf(&buffer[20],"%f %f %f",&xr,&yr,&zr);
                w_rls[0][i]=atoms[grp.atoms[i]].mass;
                refCrd[0][i][0]=10.0*xr;
                refCrd[0][i][1]=10.0*yr;
                refCrd[0][i][2]=10.0*zr;
                refCOM[0].x+=w_rls[0][i]*10.0*xr;
                refCOM[0].y+=w_rls[0][i]*10.0*yr;
                refCOM[0].z+=w_rls[0][i]*10.0*zr;
                totMass+=w_rls[0][i];
        }
        fclose(ref);
        refCOM[0].x/=totMass;
        refCOM[0].y/=totMass;
        refCOM[0].z/=totMass;

        for(i=0;i<nRef;i++)
        {
                refCrd[0][i][0]-=refCOM[0].x;
                refCrd[0][i][1]-=refCOM[0].y;
                refCrd[0][i][2]-=refCOM[0].z;
        }
        x[0]=(rvec*)save_malloc(nRef*sizeof(rvec));

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
        FILE *vin,*out;
        char fnCOM[100],fnTop[100],fnMatrix[100],fnOut[100];
        char fnOut2[100];
        t_atom *atoms;
        t_vecR *siteCrds;
        t_mol *mols;
        int nAtoms,nMols,nSites,format;
        int nFrames=0;
        real dt=1.0;
        real time=0.0;
        t_box box;
        int wc=0;
        int i=0;
        int j=0;
        int k,l,m;
        rvec *x;
        real *w_rls;
        rvec *refCrd;
        t_vecR refCOM;
        int bufferpos;
/*MODULAR*/
        double *sqm;
        char fnRef[100];
        int nDOF;
        int nCorr;
        int binaryMatrix=0;
        t_mat mat;
        int magic1,magic2;
        double **rcu,**rcv;
/*MODULAR*/

        printTitle();

/*PARTIALLY MODULAR*/
        /*PARSE COMMAND LINE INPUT HERE*/
        if(argc!=2) {
                printKeys();
                fatal("no input file specified\n");
        }
        strcpy(fnCOM,argv[1]);
        getInput(fnCOM,fnTop,
                fnRef,fnMatrix,
                &nCorr,&binaryMatrix,
                fnOut);
/*PARTIALLY MODULAR*/

        /*READ TOPOLOGY HERE*/
        readTOP(fnTop,&atoms,&nAtoms,&siteCrds,&nSites,&mols,&nMols);
        /*mols[nMols] is used to store bonds, angles, dihedrals and impropers between residues*/

/*MODULAR*/
        nDOF=3*nAtoms;
        sqm=(double*)save_malloc(nDOF*sizeof(double));
        k=0;
        for(i=0;i<nAtoms;i++) {
                sqm[k  ]=sqrt(atoms[i].mass);
                sqm[k+1]=sqrt(atoms[i].mass);
                sqm[k+2]=sqrt(atoms[i].mass);
                k+=3;
        }
 
        mat.m=(double*)save_malloc(nDOF*nDOF*sizeof(double));
        mat.nu=nDOF;
        mat.nv=nDOF;
/*MODULAR*/

        l=0;
        printf(" - reading input trajectory\n");fflush(stdout);
        /*HAPPY COMPUTING!!!*/
        for(i=0;i<nCorr;i++)
        {
                
                printf("step %d of %d%c",i+1,nCorr,(char)13); fflush(stdout);
        }
        printf("\n");

        printf(" - writing frequency-dependent velocity cross correlation matrix to output file\n");fflush(stdout);
        if(binaryMatrix!=0) {
                printf(" - writing cross correlation matrix in binary format\n");
                magic1=3*sizeof(int);
                magic2=nDOF*nDOF*sizeof(double);
                sprintf(fnOut2,"matrix_%s.mmat",fnOut);
                out=fopen(fnOut2,"wb");
                printf("   -> opening file: %s\n",fnOut2);fflush(stdout);
                fwrite(&magic1,sizeof(int),1,out);
                fwrite(&nCorr,sizeof(int),1,out);
                fwrite(&nDOF,sizeof(int),1,out);
                fwrite(&nDOF,sizeof(int),1,out);
                fwrite(&magic1,sizeof(int),1,out);
                for(i=0;i<nCorr;i++) {
                        fwrite(&magic2,sizeof(int),1,out);
                        for(j=0;j<nDOF*nDOF;j++) {
                                fwrite(&mat[i].m[j],sizeof(double),1,out);
                        }
                        fwrite(&magic2,sizeof(int),1,out);
                }
                fclose(out);
        } else {
                printf(" - writing cross correlation matrix in ASCII format\n");
                magic1=3*sizeof(int);
                magic2=nDOF*nDOF*sizeof(double);
                sprintf(fnOut2,"matrix_%s.mmat",fnOut);
                out=fopen(fnOut2,"w");
                fprintf(out,"%d\n%d %d %d\n%d\n",magic1,nCorr,nDOF,nDOF,magic1);
                for(i=0;i<nCorr;i++) {
                        fprintf(out,"%d\n",magic2);
                        for(j=0;j<nDOF*nDOF;j++) {
                                fprintf(out," %e",mat[i].m[j]);
                        }
                        fprintf(out,"\n%d\n",magic2);
                }
                fclose(out);
        }
        printf("   -> closing file: %s\n",fnOut2);fflush(stdout);
        printf(" - completed writing frequency-dependent velocity cross correlation matrix to file\n");fflush(stdout);

        return 0;
}

