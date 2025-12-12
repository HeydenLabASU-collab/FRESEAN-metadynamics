#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int offset;
    int nAtoms;
    double *masses;
    double totMass;
    float COM[3];
} t_residue;

typedef struct {
    int align;
    int centerCOM;
    int subCOMvel;
    int placeCOMInBox;
    int rotate;
    int rotVel;
    int nAtomsRef;
    double *massesRef;
    double totMassRef;
    float *crdRefFix;
    float COMrefFix[3];
    float *crdRef;
    int nAtomsSys;
    float *crdSys;
    float *velSys;
    int nResSys;
    t_residue *resSys;
    float box[9];
} t_align;

void *save_calloc(char *name,char *file,int line,
                  unsigned nelem,unsigned elsize)
{
    void *p;
    char err[300];

    p=NULL;
    if ((nelem==0)||(elsize==0))
        p=NULL;
    else {
        if ((p=calloc((size_t)nelem,(size_t)elsize))==NULL) {
            fprintf(stderr,"calloc for %s (nelem=%d, elsize=%d, file %s line %d)\n",name,nelem,elsize,file,line);
        }
    }
    return p;
}

void save_free(char *name,char *file,int line, void *ptr)
{
    if (ptr != NULL)
        free(ptr);
}

#define snew(ptr,nelem) (ptr)=save_calloc(#ptr,__FILE__,__LINE__,\
                        (nelem),sizeof(*(ptr)))

#define sfree(ptr) save_free(#ptr,__FILE__,__LINE__,(ptr))

typedef float matrix[3][3];

static inline void clear_mat(matrix a)
{
    const float nul=0.0;

    a[0][0]=a[0][1]=a[0][2]=nul;
    a[1][0]=a[1][1]=a[1][2]=nul;
    a[2][0]=a[2][1]=a[2][2]=nul;
}

static inline void oprod(float *a,float *b,float *c)
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
    a[k][l]=h+s*(g-h*tau);

void jacobi(double **a,int n,double d[],double **v,int *nrot)
{
    int j,i;
    int iq,ip;
    double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

    snew(b,n);
    snew(z,n);
    for (ip=0; ip<n; ip++) {
        for (iq=0; iq<n; iq++) v[ip][iq]=0.0;
        v[ip][ip]=1.0;
 }
    for (ip=0; ip<n;ip++) {
        b[ip]=d[ip]=a[ip][ip];
        z[ip]=0.0;
    }
    *nrot=0;
    for (i=1; i<=50; i++) {
        sm=0.0;
        for (ip=0; ip<n-1; ip++) {
            for (iq=ip+1; iq<n; iq++)
                sm += fabs(a[ip][iq]);
        }
        if (sm == 0.0) {
            sfree(z);
            sfree(b);
            return;
        }
        if (i < 4)
            tresh=0.2*sm/(n*n);
        else
            tresh=0.0;
        for (ip=0; ip<n-1; ip++) {
            for (iq=ip+1; iq<n; iq++) {
                g=100.0*fabs(a[ip][iq]);
                if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
                        && fabs(d[iq])+g == fabs(d[iq]))
                    a[ip][iq]=0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h=d[iq]-d[ip];
                    if (fabs(h)+g == fabs(h))
                        t=(a[ip][iq])/h;
                    else {
                        theta=0.5*h/(a[ip][iq]);
                        t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0) t = -t;
                    }
                    c=1.0/sqrt(1+t*t);
                    s=t*c;
                    tau=s/(1.0+c);
                    h=t*a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq]=0.0;
                    for (j=0; j<ip; j++) {
                        ROTATE(a,j,ip,j,iq)
	                }
                    for (j=ip+1; j<iq; j++) {
                        ROTATE(a,ip,j,j,iq)
                    }
                    for (j=iq+1; j<n; j++) {
                        ROTATE(a,ip,j,iq,j)
                    }
                    for (j=0; j<n; j++) {
                        ROTATE(v,j,ip,j,iq)
                    }
                    ++(*nrot);
                }
            }
        }
        for (ip=0; ip<n; ip++) {
            b[ip] +=    z[ip];
            d[ip]    =    b[ip];
            z[ip]    =    0.0;
        }
    }
    fprintf(stderr,"Error: Too many iterations!!!\n -> jacobi\n -> align.c\n");
}

void calc_fit_R(int natoms,double *w_rls,float *xp,float *x,matrix R)
{
    int c,r,n,j,m,i,irot;
    double **omega,**om;
    double d[2*3],xnr,xpc;
    matrix vh,vk,u;
    float mn;
    int index;
    float max_d;

    snew(omega,2*3);
    snew(om,2*3);
    for(i=0; i<2*3; i++) {
        snew(omega[i],2*3);
        snew(om[i],2*3);
    }
    
    for(i=0; i<2*3; i++) {
        d[i]=0;
        for(j=0; j<2*3; j++) {
            omega[i][j]=0;
            om[i][j]=0;
        }
    }
    
    /*calculate the matrix U*/
    clear_mat(u);
    for(n=0;(n<natoms);n++) {
        if ((mn = w_rls[n]) != 0.0) {
            for(c=0; (c<3); c++) {
	            xpc=xp[n*3+c];
	            for(r=0; (r<3); r++) {
	                xnr=x[n*3+r];
	                u[c][r]+=mn*xnr*xpc;
            	}
            }
        }
    }        
    
    /*construct omega*/
    /*omega is symmetric -> omega==omega' */
    for(r=0; r<2*3; r++) {
        for(c=0; c<=r; c++) {
            if (r>=3 && c<3) {
                omega[r][c]=u[r-3][c];
                omega[c][r]=u[r-3][c];
            } else {
                omega[r][c]=0;
                omega[c][r]=0;
            }
        }
    }
    
    /*determine h and k*/
    jacobi(omega,2*3,d,om,&irot);
    /*real     **omega = input matrix a[0..n-1][0..n-1] must be symmetric
     *int         natoms = number of rows and columns
     *real            NULL = d[0]..d[n-1] are the eigenvalues of a[][]
     *real             **v = v[0..n-1][0..n-1] contains the vectors in columns
     *int            *irot = number of jacobi rotations
     */
    
    if (irot==0) printf("IROT=0\n -> jacobi\n -> align.c");
    
    index=0; /* For the compiler only */

    /* Copy only the first two eigenvectors */    
    for(j=0; j<2; j++) {
        max_d=-1000;
        for(i=0; i<2*3; i++) {
            if (d[i]>max_d) {
                max_d=d[i];
                index=i;
            }
        }
        d[index]=-10000;
        for(i=0; i<3; i++) {
            vh[j][i]=1.41421356237*om[i][index];
            vk[j][i]=1.41421356237*om[i+3][index];
        }
    }
    /* Calculate the last eigenvector as the outer-product of the first two.
     * This insures that the conformation is not mirrored and
     * prevents problems with completely flat reference structures.
     */    
    oprod(vh[0],vh[1],vh[2]);
    oprod(vk[0],vk[1],vk[2]);

    /*determine R*/
    for(r=0; r<3; r++) {
        for(c=0; c<3; c++) {
            R[r][c] = vk[0][r]*vh[0][c] +
	                vk[1][r]*vh[1][c] +
	                vk[2][r]*vh[2][c];
        }
    }

    for(i=0; i<2*3; i++) {
        sfree(omega[i]);
        sfree(om[i]);
    }
    sfree(omega);
    sfree(om);
}

int residueCOM(t_residue *res,float *crd) {
    int i;
    float *crdPtr;
    float *velPtr;

    crdPtr=crd+(3*res->offset);
    for(i=0;i<3;i++) {
        res->COM[i]=0.0;
    }
    for(i=0;i<res->nAtoms;i++) {
        res->COM[0]+=res->masses[i]*crdPtr[3*i];
        res->COM[1]+=res->masses[i]*crdPtr[3*i+1];
        res->COM[2]+=res->masses[i]*crdPtr[3*i+2];
    }
    for(i=0;i<3;i++) {
        res->COM[i]/=res->totMass;
    }
    return 0;
}

int generalCOM(int nAtoms,double *masses,double totMass,float *crdVel,float *COM) {
    int i;
    for(i=0;i<3;i++) {
        COM[i]=0.0;
    }
    for (i=0;i<nAtoms;i++) {
        COM[0]+=masses[i]*crdVel[3*i];
        COM[1]+=masses[i]*crdVel[3*i+1];
        COM[2]+=masses[i]*crdVel[3*i+2];
    }
    for(i=0;i<3;i++) {
        COM[i]/=totMass;
    }
    return 0;
}

int subVec(int nAtoms,float *crdVel,float *sub) {
    int i;
    for(i=0;i<nAtoms;i++) {
        crdVel[3*i]-=sub[0];
        crdVel[3*i+1]-=sub[1];
        crdVel[3*i+2]-=sub[2];
    }
    return 0;
}

int addVec(int nAtoms,float *crdVel,float *add) {
    int i;
    for(i=0;i<nAtoms;i++) {
        crdVel[3*i]+=add[0];
        crdVel[3*i+1]+=add[1];
        crdVel[3*i+2]+=add[2];
    }
    return 0;
}

int alignRot(t_align *align,int doVel) {
    int j,c;
    matrix R;
    float x_old[3],v_old[3];
    
    /* Calculate the rotation matrix R */
    calc_fit_R(align->nAtomsRef,align->massesRef,align->crdRefFix,align->crdRef,R);

    /* rotate all coords */
    for(j=0;j<align->nAtomsSys;j++) {
        x_old[0]=align->crdSys[3*j];
        x_old[1]=align->crdSys[3*j+1];
        x_old[2]=align->crdSys[3*j+2];
        align->crdSys[3*j]=0.0;
        for(c=0;c<3;c++) align->crdSys[3*j]+=R[0][c]*x_old[c];
        align->crdSys[3*j+1]=0.0;
        for(c=0;c<3;c++) align->crdSys[3*j+1]+=R[1][c]*x_old[c];
        align->crdSys[3*j+2]=0.0;
        for(c=0;c<3;c++) align->crdSys[3*j+2]+=R[2][c]*x_old[c];
    }
    if(doVel!=0) {
        /* rotate all velocities */
        for(j=0;j<align->nAtomsSys;j++) {
            v_old[0]=align->velSys[3*j];
            v_old[1]=align->velSys[3*j+1];
            v_old[2]=align->velSys[3*j+2];
            align->velSys[3*j]=0.0;
            for(c=0; c<3; c++) align->velSys[3*j]+=R[0][c]*v_old[c];
            align->velSys[3*j+1]=0.0;
            for(c=0; c<3; c++) align->velSys[3*j+1]+=R[1][c]*v_old[c];
            align->velSys[3*j+2]=0.0;
            for(c=0; c<3; c++) align->velSys[3*j+2]+=R[2][c]*v_old[c];
        }
    }
    /* rotate box */
    for(j=0;j<3;j++) {
        x_old[0]=align->box[j*3+0];
        x_old[1]=align->box[j*3+1];
        x_old[2]=align->box[j*3+2];
        align->box[j*3+0]=0.0;
        for(c=0; c<3; c++) align->box[j*3+0]+=R[0][c]*x_old[c];
        align->box[j*3+1]=0.0;
        for(c=0; c<3; c++) align->box[j*3+1]+=R[1][c]*x_old[c];
        align->box[j*3+2]=0.0;
        for(c=0; c<3; c++) align->box[j*3+2]+=R[2][c]*x_old[c];
    }

    return 0;
}

int copyFloatArray(float *dest,float *ori,int n) {
    int i;
    for(i=0;i<n;i++) {
        dest[i]=ori[i];
    }
    return 0;
}

int copyDoubleArray(double *dest,double *ori,int n) {
    int i;
    for(i=0;i<n;i++) {
        dest[i]=ori[i];
    }
    return 0;
}

int setup(t_align *align,int nAtomsRef,double *massesRef,float *refCrd,int nAtomsSys,double *massesSys,int nRes,int *nAtomsPerRes,float *box) {
    int i,j;
    int offset=0;
    double *massPtr;

    align->nAtomsRef=nAtomsRef;
    /* prep reference atom coordinates */
    align->crdRefFix=(float*)malloc(3*nAtomsRef*sizeof(float));
    copyFloatArray(align->crdRefFix,refCrd,3*nAtomsRef);
    /* prep reference atom masses */
    align->massesRef=(double*)malloc(nAtomsRef*sizeof(double));
    copyDoubleArray(align->massesRef,massesRef,nAtomsRef);
    align->totMassRef=0.0;
    for(i=0;i<nAtomsRef;i++) {
        align->totMassRef+=massesRef[i];
    }
    generalCOM(nAtomsRef,massesRef,align->totMassRef,refCrd,align->COMrefFix);
    subVec(nAtomsRef,refCrd,align->COMrefFix);

    align->nAtomsSys=nAtomsSys;
    align->resSys=(t_residue*)malloc(nRes*sizeof(t_residue));
    align->nResSys=nRes;
    for(i=0;i<nRes;i++) {
        align->resSys[i].offset=offset;
        align->resSys[i].nAtoms=nAtomsPerRes[i];
        align->resSys[i].masses=(double*)malloc(nAtomsPerRes[i]*sizeof(double));
        massPtr=massesSys+offset;
        copyDoubleArray(align->resSys[i].masses,massPtr,nAtomsPerRes[i]);
        align->resSys[i].totMass=0.0;
        for(j=0;j<nAtomsPerRes[i];j++) {
            align->resSys[i].totMass+=massPtr[j];
        }
        offset+=nAtomsPerRes[i];
    }
    copyFloatArray(align->box,box,9);
    return 0;
}

int alignCrdVel(t_align *align,float *crdRef,float *velRef,float *crdSys,float *velSys,float *box) {
    int i;
    float refCOM[3];
    float refCOMvel[3];
    float halfBox[3];
    float pbcJump[0];
    float *crdPtr;

    align->crdRef=crdRef;
    align->crdSys=crdSys;
    align->velSys=velSys;
    copyFloatArray(align->box,box,9);
    
    if(align->align!=0) {
        generalCOM(align->nAtomsRef,align->massesRef,align->totMassRef,crdRef,refCOM);
        subVec(align->nAtomsRef,crdRef,refCOM);
        subVec(align->nAtomsSys,crdSys,refCOM);
        if(align->subCOMvel!=0) {
            generalCOM(align->nAtomsRef,align->massesRef,align->totMassRef,velRef,refCOMvel);
            subVec(align->nAtomsSys,velSys,refCOMvel);
        }
        if(align->placeCOMInBox!=0) {
            if(align->resSys==NULL) {
                fprintf(stderr,"Error: align->sysRes is NULL\n -> align\n -> align.c\n");
                return 1;
            }
            halfBox[0]=box[0]/2;
            halfBox[1]=box[4]/2;
            halfBox[2]=box[8]/2;
            #pragma omp parallel for private(i,crdPtr,pbcJump)
            for(i=0;i<align->nResSys;i++) {
                pbcJump[0]=0.0;
                pbcJump[1]=0.0;
                pbcJump[2]=0.0;
                residueCOM(&align->resSys[i],crdSys);
                if(align->resSys[i].COM[0]<-halfBox[0]) pbcJump[0]+=box[0];
                if(align->resSys[i].COM[0]> halfBox[0]) pbcJump[0]-=box[0];
                if(align->resSys[i].COM[1]<-halfBox[1]) pbcJump[1]+=box[4];
                if(align->resSys[i].COM[1]> halfBox[1]) pbcJump[1]-=box[4];
                if(align->resSys[i].COM[2]<-halfBox[2]) pbcJump[2]+=box[8];
                if(align->resSys[i].COM[2]> halfBox[2]) pbcJump[2]-=box[8];
                crdPtr=crdSys+(3*align->resSys[i].offset);
                addVec(align->resSys[i].nAtoms,crdPtr,pbcJump);
                addVec(1,align->resSys[i].COM,pbcJump);
            }
        }
        if(align->rotate!=0) {
            alignRot(align,align->rotVel);
        }
        if(align->centerCOM==0) {
            addVec(align->nAtomsSys,crdSys,align->COMrefFix);
        }
    }
    return 0;
}

int alignCrd(t_align *align,float *crdRef,float *crdSys,float *box) {
    int i;
    float refCOM[3];
    float halfBox[3];
    float pbcJump[0];
    float *crdPtr;

    align->crdRef=crdRef;
    align->crdSys=crdSys;
    copyFloatArray(align->box,box,9);
    
    if(align->align!=0) {
        generalCOM(align->nAtomsRef,align->massesRef,align->totMassRef,crdRef,refCOM);
        subVec(align->nAtomsRef,crdRef,refCOM);
        subVec(align->nAtomsSys,crdSys,refCOM);
        if(align->placeCOMInBox!=0) {
            if(align->resSys==NULL) {
                fprintf(stderr,"Error: align->sysRes is NULL\n -> align\n -> align.c\n");
                return 1;
            }
            halfBox[0]=box[0]/2;
            halfBox[1]=box[4]/2;
            halfBox[2]=box[8]/2;
            #pragma omp parallel for private(i,crdPtr,pbcJump)
            for(i=0;i<align->nResSys;i++) {
                pbcJump[0]=0.0;
                pbcJump[1]=0.0;
                pbcJump[2]=0.0;
                residueCOM(&align->resSys[i],crdSys);
                if(align->resSys[i].COM[0]<-halfBox[0]) pbcJump[0]+=box[0];
                if(align->resSys[i].COM[0]> halfBox[0]) pbcJump[0]-=box[0];
                if(align->resSys[i].COM[1]<-halfBox[1]) pbcJump[1]+=box[4];
                if(align->resSys[i].COM[1]> halfBox[1]) pbcJump[1]-=box[4];
                if(align->resSys[i].COM[2]<-halfBox[2]) pbcJump[2]+=box[8];
                if(align->resSys[i].COM[2]> halfBox[2]) pbcJump[2]-=box[8];
                crdPtr=crdSys+(3*align->resSys[i].offset);
                addVec(align->resSys[i].nAtoms,crdPtr,pbcJump);
                addVec(1,align->resSys[i].COM,pbcJump);
            }
        }
        if(align->rotate!=0) {
            alignRot(align,0);
        }
        if(align->centerCOM==0) {
            addVec(align->nAtomsSys,crdSys,align->COMrefFix);
        }
    }
    return 0;
}