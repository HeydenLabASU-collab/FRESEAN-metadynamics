#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <omp.h>
#include "../include/dataTypes.h"
#include "../include/alloc.h"
#include "../include/top.h"
#include "../include/io.h"
#include "../include/fatal.h"
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
        printf("binaryMatrix output (0 for human readable, 1 for binary matrix)\n");
        printf("fnOut (output matrix file)\n");
        return 0;
}

int getInput(char *fnCOM,char *fnTop,
                char *fnRef,char *fnMatrix,
                int *binaryMatrix,
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
        sscanf(buffer,"%d",binaryMatrix);
        printf("%2d -> read %20s : %d\n",i,"binaryMatrix",binaryMatrix[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnOut);
        printf("%2d -> read %20s : %s\n",i,"fnOut",fnOut);i++;

        fclose(io);
        return 0;
}

int readGROref(char *fnGRO,int nAtoms,rvec *crd)
{
        FILE *gro;
        char buffer[200];
        int i;
        float xr,yr,zr;
        int nAtomsGRO;

        gro=fopen(fnGRO,"r");
        fgets(buffer,200,gro);
        fgets(buffer,200,gro);
        sscanf(buffer,"%d",&nAtomsGRO);
        if(nAtomsGRO!=nAtoms) {
                sprintf(buffer,"nAtoms in %s (%d) and matrix (%d) don't match!\n",fnGRO,nAtomsGRO,nAtoms);
                fatal(buffer);
        }
        for(i=0;i<nAtoms;i++) {
                fgets(buffer,200,gro);
                sscanf(&buffer[20],"%f %f %f",&xr,&yr,&zr);
                crd[i][0]=10.0*xr;
                crd[i][1]=10.0*yr;
                crd[i][2]=10.0*zr;
        }
        fclose(gro);
        return 0;
}

int calcCOM(int nAtoms,real *m,real totMass,rvec *crd,rvec com)
{
        int i;

        com[0]=0.0; com[1]=0.0; com[2]=0.0;
        for(i=0;i<nAtoms;i++) {
                com[0]+=m[i]*crd[i][0];
                com[1]+=m[i]*crd[i][1];
                com[2]+=m[i]*crd[i][2];
        }
        com[0]/=totMass;
        com[1]/=totMass;
        com[2]/=totMass;

        return 0;
}

int calcRotAxisAndInertia(int nAtoms,real *m,rvec *crd,rvec *axes)
{
        double tensor[9];
        gsl_matrix *matrix;
        gsl_matrix *eVec;
        gsl_vector *eVal;
        gsl_eigen_symmv_workspace *ws;

        int i,j;
        /* compute elements of inertia tensor */
        for(i=0;i<9;i++) {
                tensor[i]=0.0;
        }
        for(i=0;i<nAtoms;i++) {
                tensor[0]+=m[i]*(crd[i][1]*crd[i][1]+crd[i][2]*crd[i][2]);
                tensor[4]+=m[i]*(crd[i][0]*crd[i][0]+crd[i][2]*crd[i][2]);
                tensor[8]+=m[i]*(crd[i][0]*crd[i][0]+crd[i][1]*crd[i][1]);
                tensor[1]-=m[i]*crd[i][0]*crd[i][1];
                tensor[2]-=m[i]*crd[i][0]*crd[i][2];
                tensor[5]-=m[i]*crd[i][1]*crd[i][2];
        }
        /* fill symmetric elements */
        tensor[3]=tensor[1];
        tensor[6]=tensor[2];
        tensor[7]=tensor[5];

        /* allocate GSL matrices, vectors and workspace */
        matrix=(gsl_matrix*)gsl_matrix_alloc(3,3);
        eVal=(gsl_vector*)gsl_vector_alloc(3);
        eVec=(gsl_matrix*)gsl_matrix_alloc(3,3);
        ws=(gsl_eigen_symmv_workspace*)gsl_eigen_symmv_alloc(12);
        /* fill matrix with inertia tensor elements */
        for(i=0;i<3;i++) {
            for(j=0;j<3;j++) {
                gsl_matrix_set(matrix,i,j,tensor[i*3+j]);
            }
        }
        /* compute eigenvalues and eigenvectors */
        gsl_eigen_symmv(matrix,eVal,eVec,ws);
        /* sort eigenvalues and eigenvectors (large to small) */
        gsl_eigen_symmv_sort(eVal,eVec,GSL_EIGEN_SORT_VAL_DESC);

        /* store eigenvectors as rotation axes */
        for(i=0;i<3;i++) {
            for(j=0;j<3;j++) {
                axes[i][j]=gsl_matrix_get(eVec,i,j);
            }
        }

        gsl_matrix_free(matrix);
        gsl_vector_free(eVal);
        gsl_matrix_free(eVec);
        gsl_eigen_symmv_free(ws);
        return 0;
}

int cross(rvec a,rvec b,rvec c)
{
        c[0]=a[1]*b[2]-a[2]*b[1];
        c[1]=a[2]*b[0]-a[0]*b[2];
        c[2]=a[0]*b[1]-a[1]*b[0];
        return 0;
}

int getRefTransRot(char *fnRef,int nAtoms,t_atom *atoms,double **evTrRot)
{
        int i,j;
        FILE *ref;
        char buffer[200];
        float xr,yr,zr;
        real *m;
        real *sqm;
        rvec *crd;
        rvec com;
        rvec axes[3];
        double totMass=0.0;
        rvec tmp;
        double norm;

        crd=(rvec*)save_malloc(nAtoms*sizeof(rvec));
        readGROref(fnRef,nAtoms,crd);

        allocReals(&m,nAtoms);
        allocReals(&sqm,nAtoms);
        for(i=0;i<nAtoms;i++) {
                m[i]=atoms[i].mass;
                sqm[i]=sqrt(m[i]);
                totMass+=m[i];
        }
        calcCOM(nAtoms,m,totMass,crd,com);

        for(i=0;i<nAtoms;i++)
        {
                crd[i][0]-=com[0];
                crd[i][1]-=com[1];
                crd[i][2]-=com[2];
        }
        calcRotAxisAndInertia(nAtoms,m,crd,axes);

        /* generate idealized translation and rotation modes */
        for(i=0;i<nAtoms;i++) {
                /* rigid body translation in x */
                evTrRot[0][3*i+0]=sqm[i]*1.0;
                evTrRot[0][3*i+1]=sqm[i]*0.0;
                evTrRot[0][3*i+2]=sqm[i]*0.0;
                /* rigid body translation in y */
                evTrRot[1][3*i+0]=sqm[i]*0.0;
                evTrRot[1][3*i+1]=sqm[i]*1.0;
                evTrRot[1][3*i+2]=sqm[i]*0.0;
                /* rigid body translation in z */
                evTrRot[2][3*i+0]=sqm[i]*0.0;
                evTrRot[2][3*i+1]=sqm[i]*0.0;
                evTrRot[2][3*i+2]=sqm[i]*1.0;
                /* rigid body rotation around inertia axis 1 */
                cross(axes[0],crd[i],tmp);
                evTrRot[3][3*i+0]=sqm[i]*tmp[0];
                evTrRot[3][3*i+1]=sqm[i]*tmp[1];
                evTrRot[3][3*i+2]=sqm[i]*tmp[2];
                /* rigid body rotation around inertia axis 2 */
                cross(axes[1],crd[i],tmp);
                evTrRot[4][3*i+0]=sqm[i]*tmp[0];
                evTrRot[4][3*i+1]=sqm[i]*tmp[1];
                evTrRot[4][3*i+2]=sqm[i]*tmp[2];
                /* rigid body rotation around inertia axis 3 */
                cross(axes[2],crd[i],tmp);
                evTrRot[5][3*i+0]=sqm[i]*tmp[0];
                evTrRot[5][3*i+1]=sqm[i]*tmp[1];
                evTrRot[5][3*i+2]=sqm[i]*tmp[2];
        }
        /* normalize idealized modes */
        for(i=0;i<6;i++) {
                norm=0;
                for(j=0;j<3*nAtoms;j++) {
                        norm+=evTrRot[i][j]*evTrRot[i][j];
                }
                norm=sqrt(norm);
                for(j=0;j<3*nAtoms;j++) {
                        evTrRot[i][j]/=norm;
                }
        }
        return 0;
}

typedef struct {
        double *m;
        int nu,nv;
} t_mat;

double getMat(t_mat mat,int u,int v) {
        int i;
        i=u*mat.nv+v;
        return mat.m[i];
}

int setMat(t_mat mat,int u,int v,double muv) {
        int i;
        i=u*mat.nv+v;
        mat.m[i]=muv;
        return 0;
}

int main(int argc, char *argv[])
{
        FILE *in,*out;;
        char fnCOM[100],fnTop[100],fnMatrix[100],fnOut[100];
        t_atom *atoms;
        t_vecR *siteCrds;
        t_mol *mols;
        int nAtoms,nMols,nSites,format;
        int nFrames=0;
        real dt=1.0;
        real time=0.0;
        t_box box;
        int wc=0;
        int i,j,k,l;
        rvec *x;
        real *w_rls;
        rvec *refCrd;
        t_vecR refCOM;
        int bufferpos;
/*MODULAR*/
        double **evTrRot;
        char fnRef[100];
        int nDOF;
        int nCorr;
        int n1,n2;
        int bc1,bc2;
        int binaryMatrix=0;
        t_mat mat,matFilter,res;
        int magic1,magic2;
        double **rcu,**rcv;
        double tmp;
        double a,b,c,p;
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
                &binaryMatrix,
                fnOut);
/*PARTIALLY MODULAR*/

        /*READ TOPOLOGY HERE*/
        readTOP(fnTop,&atoms,&nAtoms,&siteCrds,&nSites,&mols,&nMols);
        /*mols[nMols] is used to store bonds, angles, dihedrals and impropers between residues*/

/*MODULAR*/
        /* generate idealized translation and rotation modes */
        nDOF=3*nAtoms;
        evTrRot=(double**)save_malloc(6*sizeof(double*));
        for(i=0;i<6;i++) {
                evTrRot[i]=(double*)save_malloc(nDOF*sizeof(double));
        }
        getRefTransRot(fnRef,nAtoms,atoms,evTrRot);
 
        /* allocate a single matrix */
        mat.m=(double*)save_malloc(nDOF*nDOF*sizeof(double));
        mat.nu=nDOF;
        mat.nv=nDOF;
        /* allocate a single matrix */
        res.m=(double*)save_malloc(nDOF*nDOF*sizeof(double));
        res.nu=nDOF;
        res.nv=nDOF;
        /* allocate filter matrix */
        matFilter.m=(double*)save_malloc(nDOF*nDOF*sizeof(double));
        matFilter.nu=nDOF;
        matFilter.nv=nDOF;
        /* fill filter matrix */
        for(i=0;i<nDOF*nDOF;i++) {
                matFilter.m[i]=0.0;
        }
        for(i=0;i<nDOF;i++) {
                setMat(matFilter,i,i,1.0);
        }
        for(i=0;i<6;i++) {
                for(j=0;j<nDOF;j++) {
                        for(k=0;k<nDOF;k++) {
                                setMat(matFilter,j,k,
                                        getMat(matFilter,j,k)-evTrRot[i][j]*evTrRot[i][k]
                                );
                        }
                }
        }

        /* prepare reading binary input matrix */
        in=fopen(fnMatrix,"rb");
        fread(&bc1,sizeof(int),1,in);
        fread(&nCorr,sizeof(int),1,in);
        fread(&n1,sizeof(int),1,in);
        fread(&n2,sizeof(int),1,in);
        fread(&bc2,sizeof(int),1,in);
        if(bc2!=bc1) {
                printf("format error in file: %s\n",fnMatrix);
                exit(1);
        }
        if(n1!=n2) {
                printf("matrices in file %s are not square: %d x %d\n",fnMatrix,n1,n2);
                exit(1);
        }

        /* prepare writing output matrix */
        magic1=3*sizeof(int);
        magic2=nDOF*nDOF*sizeof(double);
        if(binaryMatrix!=0) {
                out=fopen(fnOut,"wb");
                fwrite(&magic1,sizeof(int),1,out);
                fwrite(&nCorr,sizeof(int),1,out);
                fwrite(&nDOF,sizeof(int),1,out);
                fwrite(&nDOF,sizeof(int),1,out);
                fwrite(&magic1,sizeof(int),1,out);
        } else {
                out=fopen(fnOut,"w");
                fprintf(out,"%d\n%d %d %d\n%d\n",magic1,nCorr,nDOF,nDOF,magic1);
        }

        /* reading correlation matrices from file */
        printf(" Reading and filtering correlation matrices from file %s\n",fnMatrix);
        printf(" Writing filtered correlation matrices to file %s\n",fnOut);
        for(i=0;i<nCorr;i++) {
                fread(&bc1,sizeof(int),1,in);
                for(j=0;j<n1;j++) {
                        for(k=0;k<n2;k++) {
                                fread(&tmp,sizeof(double),1,in);
                                setMat(mat,j,k,tmp);
                        }
                }
                fread(&bc2,sizeof(int),1,in);
                if(bc2!=bc1) {
                        printf("format error in file: %s\n",fnMatrix);
                        exit(1);
                }
                /* filter out translational and rotational components */
                for(j=0;j<nDOF;j++) {
                        for(k=0;k<nDOF;k++) {
                                p=0.0;
                                for(l=0;l<nDOF;l++) {
                                        a=getMat(mat,j,l);
                                        b=getMat(matFilter,l,k);
                                        p+=a*b;
                                }
                                setMat(res,j,k,p);
                        }
                }
                for(j=0;j<nDOF*nDOF;j++) {
                        mat.m[j]=res.m[j];
                }
                for(j=0;j<nDOF;j++) {
                        for(k=0;k<nDOF;k++) {
                                p=0.0;
                                for(l=0;l<nDOF;l++) {
                                        a=getMat(matFilter,j,l);
                                        b=getMat(mat,l,k);
                                        p+=a*b;
                                }
                                setMat(res,j,k,p);
                        }
                }
                /* write filtered matrix to output file */
                if(binaryMatrix!=0) {
                        fwrite(&magic2,sizeof(int),1,out);
                        for(j=0;j<nDOF*nDOF;j++) {
                                fwrite(&res.m[j],sizeof(double),1,out);
                        }
                        fwrite(&magic2,sizeof(int),1,out);
                } else {
                        fprintf(out,"%d\n",magic2);
                        for(j=0;j<nDOF*nDOF;j++) {
                                fprintf(out," %e",res.m[j]);
                        }
                        fprintf(out,"\n%d\n",magic2);
                }
                progressBar(i,1,nCorr);
        }
        fclose(in);
        fclose(out);
/*MODULAR*/

        return 0;
}

