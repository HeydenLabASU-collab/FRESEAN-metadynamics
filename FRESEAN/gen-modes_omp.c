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
        printf("THE PURPOSE OF THIS PROGRAM IS TO COMPUTE\n");
        printf("A CROSS-CORRELATION MATRIX OF WEIGHTED VELOCITIES\n");
        printf("FOR ADDITIONAL PROCESSING (BINARY OR ASCII FORMAT).\n");
        printf("THE ANALYSIS IS PERFORMED FOR A SELECTED GROUP OF\n");
        printf("ATOMS. THIS GROUP OF ATOMS SHOULD ALSO BE USED\n");
        printf("FOR TRANS/ROT FITTING. THIS OPERATION SHOULD BE\n");
        printf("PERFORMED BY THIS PROGRAM TO ENSURE THAT ALL\n");
        printf("ROTATIONS ARE ALSO APPLIED TO VELOCITIES!\n");
        printf("OPTIONALLY, THIS PROGRAM ALSO COMPUTES\n");
        printf("GENERALIZED NORMAL MODES TO DESCRIBE VIBRATIONS\n");
        printf("OF THE ATOMS IN THE SELECTED GROUP.\n");
        printf("IF YOU USE THIS OPTION, READ & CITE:\n");
        printf("G. Mathias & M. Baer\n");
        printf("J. Chem. Theory Comput. 2011, 7, 2028-2039.\n");
        printf("\n");
        printf("Version 1.1: May 3, 2022\n");
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
        printf("fnTop (topology file [.mtop])\nfnCrd\nfnVel (if format xyz,crd,dcd)\nfnJob (atom group file [.job])\n");
        printf("nRead (Number of frames to read)\nanalysisInterval (dFrame)\n");
        printf("fnRef (reference coordinate file)\nalignGrp (atom group to align trajectory to, from job file)\nanalyzeGrp (atom group to perform analysis on, from job file)\nwrap (atom group to wrap trajectory on, from job file)\n");
        printf("nCorr (number of points in time correlation function)\nwinSigma (width of smoothing funciton)\nbinaryMatrix (0 for human readable, 1 for binary matrix)\ndoGenModes (0 for no, 1 for yes)\n");
        printf("convergence (convergence value for generalized normal modes)\nmaxIter (maximum iterations for generalized normal modes)\nfnOut (output matrix file)\n");
        return 0;
}

int getInput(char *fnCOM,char *fnTop,char *fnCrd,char *fnVel,char *fnJob,
                int *nRead,int *analysisInterval,
                char *fnRef,int *alignGrp,int *analyzeGrp,int *wrap,
                int *nCorr,real *winSigma,int *binaryMatrix,int *doGenModes,
                real *convergence,int *maxIter,char *fnOut,
                int needVelo) {
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
        sscanf(buffer,"%s",fnCrd);
        getFormat(fnCrd,&format);
        printf("%2d -> read %20s : %s  ",i,"fnCrd",fnCrd);i++;
        if(format==1) printf("(xyz)\n");
        else if(format==2) printf("(crd)\n");
        else if(format==3) printf("(gro)\n");
        else if(format==4) printf("(trr)\n");
        else if(format==5) printf("(dcd)\n");
        if(needVelo==1 && (format==1 || format==2 || format==5)) {
                fgets(buffer,300,io);
                sscanf(buffer,"%s",fnVel);
                printf("%2d -> read %20s : %s\n",i,"fnVel",fnVel);i++;
        }

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnJob);
        printf("%2d -> read %20s : %s\n",i,"fnJob",fnJob);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nRead);
        printf("%2d -> read %20s : %d\n",i,"nRead",nRead[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",analysisInterval);
        printf("%2d -> read %20s : %d\n",i,"analysisInterval",analysisInterval[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnRef);
        printf("%2d -> read %20s : %s\n",i,"fnRef",fnRef);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",alignGrp);
        printf("%2d -> read %20s : %d\n",i,"alignGrp",alignGrp[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",analyzeGrp);
        printf("%2d -> read %20s : %d\n",i,"analyzeGrp",analyzeGrp[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",wrap);
        printf("%2d -> read %20s : %d\n",i,"wrap",wrap[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nCorr);
        printf("%2d -> read %20s : %d\n",i,"nCorr",nCorr[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        winSigma[0]=tmp;
        printf("%2d -> read %20s : %f\n",i,"winSigma",winSigma[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",binaryMatrix);
        printf("%2d -> read %20s : %d\n",i,"binaryMatrix",binaryMatrix[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",doGenModes);
        printf("%2d -> read %20s : %d\n",i,"doGenModes",doGenModes[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        convergence[0]=tmp;
        printf("%2d -> read %20s : %f\n",i,"convergence",convergence[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",maxIter);
        printf("%2d -> read %20s : %d\n",i,"maxIter",maxIter[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnOut);
        printf("%2d -> read %20s : %s\n",i,"fnOut",fnOut);i++;

        fclose(io);
        return 0;
}

int readTRJbuffer(int *pos,real *time,real *timebuffer,t_box *box,t_box *boxbuffer,int nAtoms,t_atom *atoms,t_atom **atomsbuffer,int nSites,t_vecR *siteCrds,t_vecR **siteCrdsbuffer)
{
        int i,j;
        i=pos[0];
        time[0]=timebuffer[i];
        vecRCpy(boxbuffer[i].a,&box[0].a);
        vecRCpy(boxbuffer[i].b,&box[0].b);
        vecRCpy(boxbuffer[i].c,&box[0].c);
        box[0].ortho=boxbuffer[i].ortho;
        vecRCpy(boxbuffer[i].trans.c1,&box[0].trans.c1);
        vecRCpy(boxbuffer[i].trans.c2,&box[0].trans.c2);
        vecRCpy(boxbuffer[i].trans.c3,&box[0].trans.c3);
        for(j=0;j<nAtoms;j++) {
                vecRCpy(atomsbuffer[i][j].crd,&atoms[j].crd);
                vecRCpy(atomsbuffer[i][j].vel,&atoms[j].vel);
        }
        for(j=0;j<nSites;j++) {
                vecRCpy(siteCrdsbuffer[i][j],&siteCrds[j]);
        }
        pos[0]++;
        return 0;
}

int getAlignRefSpec(char *fnRef,t_grp grp,t_atom *atoms,t_box box,real **w_rls,rvec **refCrd,t_vecR *refCOM,rvec **x)
{
        int i;
        FILE *ref;
        char buffer[200];
        int nRef;
        float xr,yr,zr;
        real totMass=0.0;

        ref=fopen(fnRef,"r");
        fgets(buffer,200,ref);
        fgets(buffer,200,ref);
        sscanf(buffer,"%d",&nRef);
        if(nRef!=grp.nAtoms) {
                sprintf(buffer,"nAtoms in %s (%d) and group %s (%d) don't match!\n",fnRef,nRef,grp.name,grp.nAtoms);
                fatal(buffer);
        }
        allocReals(w_rls,grp.nAtoms);
        refCrd[0]=(rvec*)save_malloc(grp.nAtoms*sizeof(rvec));
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

int jacobiTransform(t_mat mat,int u,int v,double c,double s,double *rcu,double *rcv) {
        int i;
        double tmp;
        double muu,muv,mvu,mvv;
        getMat(mat,u,u,&muu);
        getMat(mat,u,v,&muv);
        getMat(mat,v,u,&mvu);
        getMat(mat,v,v,&mvv);

        for(i=0;i<mat.nv;i++) {
                getMat(mat,u,i,&rcu[i]);
                getMat(mat,v,i,&rcv[i]);
        }
        for(i=0;i<mat.nv;i++) {
                /*testing for i!=u && i!=v too expensive*/
                /*we'll overwrite the uu,uv,vu,vv elements later*/
                tmp=rcu[i]*c-rcv[i]*s;
                setMat(mat,u,i,tmp);
                tmp=rcv[i]*c+rcu[i]*s;
                setMat(mat,v,i,tmp);
        }
        for(i=0;i<mat.nv;i++) {
                getMat(mat,i,u,&rcu[i]);
                getMat(mat,i,v,&rcv[i]);
        }
        for(i=0;i<mat.nu;i++) {
                /*testing for i!=u && i!=v too expensive*/
                /*we'll overwrite the uu,uv,vu,vv elements later*/
                tmp=rcu[i]*c-rcv[i]*s;
                setMat(mat,i,u,tmp);
                tmp=rcv[i]*c+rcu[i]*s;
                setMat(mat,i,v,tmp);
        }
        tmp=c*(muu*c-mvu*s)-s*(muv*c-mvv*s);
        setMat(mat,u,u,tmp);
        tmp=s*(mvu*c+muu*s)+c*(mvv*c+muv*s);
        setMat(mat,v,v,tmp);
        tmp=s*(muu*c-mvu*s)+c*(muv*c-mvv*s);
        setMat(mat,u,v,tmp);
        tmp=c*(mvu*c+muu*s)-s*(mvv*c+muv*s);
        setMat(mat,v,u,tmp);
        return 0;
}

int jacobiMult(t_mat mat,int u,int v,double c,double s,double *rcu,double *rcv) {
        int i;
        double tmp;

        for(i=0;i<mat.nu;i++) {
                getMat(mat,u,i,&rcu[i]);
                getMat(mat,v,i,&rcv[i]);
        }
        for(i=0;i<mat.nu;i++) {
                tmp=rcu[i]*c-rcv[i]*s;
                setMat(mat,u,i,tmp);
                tmp=rcv[i]*c+rcu[i]*s;
                setMat(mat,v,i,tmp);
        }
        return 0;
}

int main(int argc, char *argv[])
{
        FILE  *tin,*cin,*vin,*job,*out;
        XDR xdrin;
        char fnCOM[100],fnTop[100],fnCrd[100],fnVel[100],fnJob[100],fnOut[100];
        char fnOut2[100];
        t_atom *atoms;
        t_vecR *siteCrds;
        t_mol *mols;
        t_grp *grps;
        int nAtoms,nMols,nSites,format;
        int nFrames=0;
        real dt=1.0;
        real time=0.0;
        t_box box;
        t_grpDef grpDef;
        int needVelo=1;
        int veloFreq=1;
        int wc=0;
        int i=0;
        int j=0;
        int k,l,m;
        int nRead=5000;
        int nTimesAnalysis=1000;
        int analysisInterval=10;
        int nSample=10;
        int alignGrp,analyzeGrp;
        int wrap;
        rvec *x;
        real *w_rls;
        rvec *refCrd;
        t_vecR refCOM;
        real *timebuffer;
        t_box *boxbuffer;
        t_atom **atomsbuffer;
        t_vecR **siteCrdsbuffer;
        int bufferpos;
/*MODULAR*/
        real convergence;
        int maxIter;
        double *sqm;
        char fnRef[100];
        int nDOF;
        int nTimes;
        int nCorr;
        real winSigma;
        int binaryMatrix=0;
        int doGenModes=0;
        fftw_complex **sqmv;
        fftw_plan *FTlong,*FTshort,*invFTlong;
        fftw_complex **array,**array2;
        double *window;
        t_mat *mat;
        double FTlongNorm, FTshortNorm, winNorm, df;
        t_mat eigVec,old;
        double offDiagOld,offDiagNew;
        double delta;
        int jacobiCount;
        int u,v;
        double muu,muv,mvu,mvv;
        double chi1,chi2,chi3,eta,etaSq;
        double tkl[4],secD[4];
        double c,s,t;
        double t2,t4;
        int magic1,magic2;
        int nThreads,thread;
        double **rcu,**rcv;
/*MODULAR*/

        printTitle();

        timebuffer=(real*)save_malloc(100*sizeof(real));
        boxbuffer=(t_box*)save_malloc(100*sizeof(t_box));
        atomsbuffer=(t_atom**)save_malloc(100*sizeof(t_atom*));
        siteCrdsbuffer=(t_vecR**)save_malloc(100*sizeof(t_vecR*));

/*PARTIALLY MODULAR*/
        /*PARSE COMMAND LINE INPUT HERE*/
        if(argc!=2) {
                printKeys();
                fatal("no input file specified\n");
        }
        strcpy(fnCOM,argv[1]);
        getInput(fnCOM,fnTop,fnCrd,fnVel,fnJob,
                &nRead,&analysisInterval,
                fnRef,&alignGrp,&analyzeGrp,&wrap,
                &nCorr,&winSigma,&binaryMatrix,&doGenModes,
                &convergence,&maxIter,fnOut,
                needVelo);
        if(2*nCorr-1>=nRead) fatal("nRead must be greater than 2*nCorr-1\n");
/*PARTIALLY MODULAR*/

        /*READ TOPOLOGY HERE*/
        readTOP(fnTop,&atoms,&nAtoms,&siteCrds,&nSites,&mols,&nMols);
        for(i=0;i<100;i++) {
                atomsbuffer[i]=(t_atom*)save_malloc(nAtoms*sizeof(t_atom));
                siteCrdsbuffer[i]=(t_vecR*)save_malloc(nSites*sizeof(t_vecR));
        }
        /*mols[nMols] is used to store bonds, angles, dihedrals and impropers between residues*/

        /*READ JOB DATA HERE*/
        readJob(fnJob,&grpDef);

        /*SELECT STATIC GROUPS FROM TOPOLOGY*/
        allocGrps(grpDef,nAtoms,nMols,&grps);
        makeStaticGrps(atoms,nAtoms,mols,nMols,grpDef.statGrpDef,grpDef.nStat,grps);

        /*PREP CRD AND/OR VEL READ HERE*/
        prepTRJinput(fnCrd,fnVel,&cin,&vin,&xdrin,&format,&veloFreq,&time,&dt,&nFrames,needVelo,nAtoms,argc);

/*MODULAR*/
        #pragma omp parallel
        {
                nThreads=omp_get_num_threads();
        }
        nDOF=3*grps[analyzeGrp].nAtoms;
        nTimes=nRead/analysisInterval;
        sqm=(double*)save_malloc(nDOF*sizeof(double));
        k=0;
        for(i=0;i<grps[analyzeGrp].nAtoms;i++) {
                j=grps[analyzeGrp].atoms[i];
                sqm[k  ]=sqrt(atoms[j].mass);
                sqm[k+1]=sqrt(atoms[j].mass);
                sqm[k+2]=sqrt(atoms[j].mass);
                k+=3;
        }
        sqmv=(fftw_complex**)save_malloc(nDOF*sizeof(fftw_complex*));
        for(i=0;i<nDOF;i++) {
                sqmv[i]=(fftw_complex*)save_malloc(nTimes*sizeof(fftw_complex));
                for(j=0;j<nTimes;j++) {
                        sqmv[i][j][0]=0.0;
                        sqmv[i][j][1]=0.0;
                }
        }
        array=(fftw_complex**)save_malloc(nThreads*sizeof(fftw_complex*));
        array2=(fftw_complex**)save_malloc(nThreads*sizeof(fftw_complex*));
        FTlong=(fftw_plan*)save_malloc(nThreads*sizeof(fftw_plan));
        FTshort=(fftw_plan*)save_malloc(nThreads*sizeof(fftw_plan));
        invFTlong=(fftw_plan*)save_malloc(nThreads*sizeof(fftw_plan));
        for(i=0;i<nThreads;i++) {
                array[i]=(fftw_complex*)fftw_malloc(nTimes*sizeof(fftw_complex));
                array2[i]=(fftw_complex*)fftw_malloc((2*nCorr-1)*sizeof(fftw_complex));
                FTlong[i]=fftw_plan_dft_1d(nTimes,array[i],array[i],FFTW_FORWARD,FFTW_ESTIMATE);
                FTshort[i]=fftw_plan_dft_1d(2*nCorr-1,array2[i],array2[i],FFTW_FORWARD,FFTW_ESTIMATE);
                invFTlong[i]=fftw_plan_dft_1d(nTimes,array[i],array[i],FFTW_BACKWARD,FFTW_ESTIMATE);
        }
        window=(double*)save_malloc(nCorr*sizeof(double));
        mat=(t_mat*)save_malloc(nCorr*sizeof(t_mat));
        for(i=0;i<nCorr;i++) {
                mat[i].m=(double*)save_malloc(nDOF*nDOF*sizeof(double));
                mat[i].nu=nDOF;
                mat[i].nv=nDOF;
        }
        eigVec.m=(double*)save_malloc(nDOF*nDOF*sizeof(double));
        eigVec.nu=nDOF;
        eigVec.nv=nDOF;
        old.m=(double*)save_malloc(nDOF*nDOF*sizeof(double));
        old.nu=nDOF;
        old.nv=nDOF;
        rcu=(double**)save_malloc(nThreads*sizeof(double*));
        rcv=(double**)save_malloc(nThreads*sizeof(double*));
        for(i=0;i<nThreads;i++) {
                rcu[i]=(double*)save_malloc(nDOF*sizeof(double));
                rcv[i]=(double*)save_malloc(nDOF*sizeof(double));
        }
/*MODULAR*/

        l=0;
        printf(" - reading input trajectory\n");fflush(stdout);
        /*HAPPY COMPUTING!!!*/
        for(i=0;i<nRead;i++)
        {
                if(i%100==0) {
                        for(j=0;j<100;j++) {
                                readTRJ(i+j,format,cin,vin,&xdrin,dt,&timebuffer[j],&boxbuffer[j],nAtoms,nSites,needVelo,veloFreq,wc,atomsbuffer[j],siteCrdsbuffer[j]);
                        }
                        bufferpos=0;
                }
                readTRJbuffer(&bufferpos,&time,timebuffer,&box,boxbuffer,nAtoms,atoms,atomsbuffer,nSites,siteCrds,siteCrdsbuffer);
                if(i==1) {
                        dt=(timebuffer[1]-timebuffer[0]);
                        df=33.3564/((2*nCorr-1)*dt*analysisInterval);
                }

                if(i==0 && alignGrp!=-1) { 
                        getAlignRefSpec(fnRef,grps[alignGrp],atoms,box,&w_rls,&refCrd,&refCOM,&x);
                }
                if(i%analysisInterval==0) {
                        if(alignGrp!=-1) {
                                alignGroupVel(&grps[alignGrp],atoms,nAtoms,&box,w_rls,refCrd,refCOM,x);
                        }
                        if(wrap==1 || wrap ==2) wrapTrajectoryMol0(atoms,mols,nMols,box);
                        getMolCOM(atoms,mols,nMols,needVelo);
                        if(wrap==2) wrapTrajectoryMolCOM(atoms,mols,nMols,box);
                        getMolDM(atoms,mols,nMols,needVelo);
                }
                
                if(i%analysisInterval==0)
                {
                        updateGrps(atoms,nAtoms,mols,nMols,box,grpDef,grps);
/*MODULAR*/
                        /*store the sqrt(mass) weighted velocities in array*/
                        /*we switch from counting atoms to degrees of freedom*/
                        m=0;
                        for(j=0;j<grps[analyzeGrp].nAtoms;j++) {
                                k=grps[analyzeGrp].atoms[j];
                                sqmv[m  ][l][0]=sqm[m  ]*atoms[k].vel.x;
                                sqmv[m  ][l][1]=0.0;
                                sqmv[m+1][l][0]=sqm[m+1]*atoms[k].vel.y;
                                sqmv[m+1][l][1]=0.0;
                                sqmv[m+2][l][0]=sqm[m+2]*atoms[k].vel.z;
                                sqmv[m+2][l][1]=0.0;
                                m+=3;
                        }
                        l++;
/*MODULAR*/
                }
                printf("step %d of %d%c",i+1,nRead,(char)13); fflush(stdout);
        }
        printf("\n");
        closeInput(format,needVelo,cin,vin,&xdrin);
        printf(" - closed input trajectory\n");fflush(stdout);

        printf(" - preparing Gaussian window function\n");fflush(stdout);
        /*preparing Gaussian window function for smoothin*/
        /*winSigma is user defined parameter and given in wavenumbers (cm^-1)*/
        /*we could prepare this analytically, but it is done only once*/
        /*hence we define the window function in the frequency domain and FT into time domain*/
        winNorm=1.0/sqrt(2*3.14159*winSigma*winSigma);
        printf("   -> time resolution is %6.3f ps\n",dt*analysisInterval);fflush(stdout);
        printf("   -> frequency resolution is %6.3f cm^-1\n",df);fflush(stdout);
        FTshortNorm=1.0/sqrt(2*nCorr-1);
        /*first half of window function*/
        for(i=0;i<nCorr;i++) {
                if(i*df>7*winSigma) {
                        array2[0][i][0]=0.0;
                } else {
                        array2[0][i][0]=winNorm*exp(-1.0*(i*df)*(i*df)/(2*winSigma*winSigma));
                }
                array2[0][i][1]=0.0;
        }
        /*now creating second half to be exactly symmetric*/
        j=1;
        for(i=nCorr;i<2*nCorr-1;i++) {
                array2[0][i][0]=array2[0][nCorr-j][0];
                array2[0][i][1]=array2[0][nCorr-j][1];
                j++;
        }
        /*FT of window function from frequency domain into time domain*/
        fftw_execute(FTshort[0]);
        /*storing window function in its own array, no imaginary parts due to symmetry*/
        for(i=0;i<nCorr;i++) {
                window[i]=FTshortNorm*array2[0][i][0];
        }
        printf(" - Gaussian window ready\n");fflush(stdout);

        printf(" - perform FT for all velocity time series\n");fflush(stdout);
        /*FT of all velocity time series*/
        FTlongNorm=1.0/sqrt(nTimes);
        #pragma omp parallel for private(thread,j)
        for(i=0;i<nDOF;i++) {
                thread=omp_get_thread_num();
                for(j=0;j<nTimes;j++) {
                        array[thread][j][0]=sqmv[i][j][0];
                        array[thread][j][1]=sqmv[i][j][1];
                }
                fftw_execute(FTlong[thread]);
                for(j=0;j<nTimes;j++) {
                        sqmv[i][j][0]=FTlongNorm*array[thread][j][0];
                        sqmv[i][j][1]=FTlongNorm*array[thread][j][1];
                }
		//progressBar(i,1,nDOF);
        }
        printf(" - all velocities are now in the frequency domain\n");fflush(stdout);

        printf(" - computing frequency dependent matrix of velocity cross correlations\n");fflush(stdout);
        /*computing the frequency dependent matrix of time cross correlations*/
        /*in frequency domain = product with complex conjugate*/
        for(i=0;i<nDOF;i++) {
                /*need to do this only for one half of matrix*/
                /*we'll enforce symmetry*/
                #pragma omp parallel for private(thread,k,l)
                for(j=i;j<nDOF;j++) {
                        thread=omp_get_thread_num();
                        for(k=0;k<nTimes;k++) {
                                array[thread][k][0]=sqmv[i][k][0]*sqmv[j][k][0]+sqmv[i][k][1]*sqmv[j][k][1];
                                array[thread][k][1]=sqmv[j][k][0]*sqmv[i][k][1]-sqmv[i][k][0]*sqmv[j][k][1];
                        }
                        fftw_execute(invFTlong[thread]);
                        /*now we have the time cross correlation in time domain*/
                        /*imaginary components are gone, but the functionis not symmetric*/
                        /*functions for ij and ji are mirror images of each other*/
                        /*by averaging the mirror images, we obtain a symmetric correlation function*/
                        /*and symmetric matrix, hence we need to store the result twice for ij and ji*/
                        /*we need only nCorr elements though*/
                        /*we directly multiply with the previously prepared window function*/
                        array2[thread][0][0]=window[0]*FTlongNorm*array[thread][0][0];
                        array2[thread][0][1]=0.0;
                        for(k=1;k<nCorr;k++) {
                                l=nTimes-k;
                                array2[thread][k][0]=window[k]*FTlongNorm*(array[thread][k][0]+array[thread][l][0])/2.0;
                                array2[thread][k][1]=0.0;
                        }
                        /*now we generate the symmetric second half*/
                        l=1;
                        for(k=nCorr;k<2*nCorr-1;k++) {
                                array2[thread][k][0]=array2[thread][nCorr-l][0];
                                array2[thread][k][1]=0.0;
                                l++;
                        }
                        /*the FT gives us the smoothed spectrum of the time cross correlations*/
                        fftw_execute(FTshort[thread]);
                        /*now we store the result in the cros correlation matrix*/
                        for(k=0;k<nCorr;k++) {
                                setMat(mat[k],i,j,array2[thread][k][0]);
                                setMat(mat[k],j,i,array2[thread][k][0]);
                        }
                }
                printf("%c   DOF x DOF : %4d x %4d",(char)13,i+1,j+1);
        }
        printf("\n");
        printf(" - computed the frequency-dependent velocity cross correlation matrix\n");fflush(stdout);

        //k=0;
        //for(i=0;i<nDOF;i++) {
        //        for(j=0;j<nDOF;j++) {
        //                printf(" %7.4f",mat[0].m[k]);
        //                k++;
        //        }
        //        printf("\n");
        //}
        //printf("\n");

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
        
        if(doGenModes!=0) {
                printf(" - initializing eigenvector matrix (unit matrix)\n");fflush(stdout);
                for(i=0;i<nDOF;i++) {
                        for(j=0;j<nDOF;j++) {
                                setMat(old,i,j,0.0);
                                setMat(eigVec,i,j,0.0);
                        }
                        setMat(old,i,i,1.0);
                        setMat(eigVec,i,i,1.0);
                }
                printf(" - eigenvector matrix initialized\n");fflush(stdout);
                delta=2*convergence;
                jacobiCount=0;
                offDiagOld=0.0;
                for(k=0;k<nCorr;k++) {
                        for(u=0;u<nDOF;u++) {
                                for(v=0;v<nDOF;v++) {
                                        if(u!=v) {
                                                l=u*nDOF+v;
                                                offDiagOld+=fabs(mat[k].m[l]);
                                        }
                                }
                        }
                }
                offDiagOld/=(nDOF*nDOF-nDOF);
                printf(" - initial average off-diagonal sum offDiag: %e\n",offDiagOld);
                printf(" - beginning Jacobi sweeps\n");fflush(stdout);
                while(delta>convergence && jacobiCount<maxIter) {
                        if(jacobiCount>0) {
                                sprintf(fnOut2,"eigVec_%s.mat",fnOut);
                                out=fopen(fnOut2,"w");
                                k=0;
                                for(i=0;i<nDOF;i++) {
                                        for(j=0;j<nDOF;j++) {
                                                fprintf(out," %12.5e",eigVec.m[k]);
                                                k++;
                                        }
                                        fprintf(out,"\n");
                                }
                                fclose(out);
                                sprintf(fnOut2,"eigVal_%s.dat",fnOut);
                                out=fopen(fnOut2,"w");
                                fprintf(out,"#df= %6.3f cm^-1  nf= %5d  nDOF= %5d\n",df,nCorr,nDOF);
                                for(i=0;i<nDOF;i++) {
                                        for(j=0;j<nCorr;j++) {
                                                getMat(mat[j],i,i,&muu);
                                                fprintf(out," %12.5e",muu);
                                        }
                                        fprintf(out,"\n");
                                }
                                fclose(out);
                                printf("   wrote output files: eigVec_%s.mat & eigVal_%s.dat\n",fnOut,fnOut);
                        }
                        for(u=0;u<nDOF-1;u++) {
                                for(v=u+1;v<nDOF;v++) {
                                        chi1=0.0;
                                        chi2=0.0;
                                        chi3=0.0;
                                        for(i=0;i<nCorr;i++) {
                                                getMat(mat[i],u,u,&muu);
                                                getMat(mat[i],u,v,&muv);
                                                getMat(mat[i],v,u,&mvu);
                                                getMat(mat[i],v,v,&mvv);
                                                chi1+=muv*muv;
                                                chi2+=muv*(muu-mvv);
                                                chi3+=(muu-mvv)*(muu-mvv);
                                        }
                                        eta=(4*chi1-chi3)/chi2;
                                        etaSq=eta*eta;
                                        i=0;
                                        for(k=1;k<=2;k++) {
                                                for(l=1;l<=2;l++) {
                                                        tkl[i]=0.25*(-1.0*eta+pow(-1,k)*sqrt(16+etaSq)+pow(-1,l)*sqrt(2.0)*sqrt(16+etaSq-pow(-1,k)*eta*sqrt(16+etaSq)));
                                                        t2=tkl[i]*tkl[i];
                                                        t4=t2*t2;
                                                        secD[i]=-2*(chi2/pow(t2+4,4))*(2*tkl[i]*(t4-14*t2+9)+eta*(3*t4-8*t2+1));
                                                        i++;
                                                }
                                        }
                                        //for(i=0;i<4;i++) {
                                        //               printf("t[%d] %f sD[%d] %f\n",i,tkl[i],i,secD[i]);fflush(stdout);
                                        //}
                                        i=0;
                                        while((fabs(tkl[i])>1 || secD[i]<0) && i<4) i++;
                                        if(i==4) fatal("Jacobi minimization failed\n");
                                        t=tkl[i];
                                        c=1.0/sqrt(t*t+1);
                                        s=t*c;
                                        //printf("t %f c %f s %f\n",t,c,s);fflush(stdout);
                                        jacobiMult(eigVec,u,v,c,s,rcu[0],rcv[0]);
                                        #pragma omp parallel for private(thread)
                                        for(k=0;k<nCorr;k++) {
                                                thread=omp_get_thread_num();
                                                jacobiTransform(mat[k],u,v,c,s,rcu[thread],rcv[thread]);
                                        }
                                        //printf("transformed matrix 0\n");
                                        //for(i=0;i<nDOF;i++) {
                                        //        for(j=0;j<nDOF;j++) {
                                        //                getMat(mat[0],i,j,&muu);
                                        //                printf(" %7.4f",muu);
                                        //        }
                                        //        printf("\n");
                                        //}
                                        //exit(1);
                                        printf("%c   DOF x DOF : %4d x %4d",(char)13,u+1,v+1);
                                }
                        }
                        /*to monitor changes in eigenvector matrix*/
                        //delta=0.0;
                        //for(i=0;i<nDOF*nDOF;i++) {
                        //        delta+=fabs(eigVec.m[i]-old.m[i]);
                        //}
                        //delta/=nDOF*nDOF;
                        offDiagNew=0.0;
                        for(k=0;k<nCorr;k++) {
                                for(u=0;u<nDOF;u++) {
                                        for(v=0;v<nDOF;v++) {
                                                if(u!=v) {
                                                        l=u*nDOF+v;
                                                        offDiagNew+=fabs(mat[k].m[l]);
                                                }
                                        }
                                       }
                        }
                        offDiagNew/=(nDOF*nDOF-nDOF);
                        delta=offDiagOld-offDiagNew;
                        offDiagOld=offDiagNew;
                        for(i=0;i<nDOF*nDOF;i++) {
                                old.m[i]=eigVec.m[i];
                        }
                        jacobiCount++;
                        printf("\n");
                        printf(" - completed Jacobi sweep %4d: --> offDiag: %e  --> delta = %f\n",jacobiCount,offDiagNew,delta);
                }
        }

        return 0;
}

