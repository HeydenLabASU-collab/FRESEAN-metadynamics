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

int printTitle() {
        printf("***********************************************\n");
        printf("NEW CODE: TESTING REQUIRED\n");
        printf("***********************************************\n\n");
        printf("THE PURPOSE OF THIS PROGRAM IS TO PROJECT\n");
        printf("A TRAJECTORY OF ATOMIC DISPLACEMENTS RELATIVE\n");
        printf("TO A REFERENCE STRUCTURE ONTO A SET OF\n");
        printf("NORMAL COORDINATES.\n");
        printf("THE REFERENCE STRUCTURE IS THE SAME ONE USED\n");
        printf("FOR THE ROTATIONAL ALIGNMENT.\n");
        printf("Version 0.9: April 26, 2023\n");
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
        printf("fnTop\nfnCrd\nfnVel (if format xyz,crd,dcd)\nfnJob\n");
        printf("nRead\nanalysisInterval\n");
        printf("fnRefAlign\nfnRefDisp\nalignGrp\nanalyzeGrp\nwrap\n");
        printf("fnEigVec\nfreqSel\nmodeStart\nmodeEnd\nfnOut\n");
        return 0;
}

int getInput(char *fnCOM,char *fnTop,char *fnCrd,char *fnVel,char *fnJob,
                int *nRead,int *analysisInterval,
                char *fnRefAlign,char *fnRefDisp,int *alignGrp,int *analyzeGrp,int *wrap,
                char *fnEigVec,char *fnOut,
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
        sscanf(buffer,"%s",fnRefAlign);
        printf("%2d -> read %20s : %s\n",i,"fnRefAlign",fnRefAlign);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnRefDisp);
        printf("%2d -> read %20s : %s\n",i,"fnRefDisp",fnRefDisp);i++;

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
        sscanf(buffer,"%s",fnEigVec);
        printf("%2d -> read %20s : %s\n",i,"fnEigVec",fnEigVec);i++;

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

int getAlignRefSpec(char *fnRef,t_grp grp,t_atom *atoms,real **w_rls,rvec **refCrd,t_vecR *refCOM,rvec **x)
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

int getDispRefSpec(char *fnRef,t_grp grp,t_atom *atoms,t_vecR **refCrd,t_vecR *refCOM) {
        int i;
        FILE *ref;
        char buffer[200];
        int nRef;
        float xr,yr,zr;
        double m,totMass;

        totMass=0.0;
        ref=fopen(fnRef,"r");
        fgets(buffer,200,ref);
        fgets(buffer,200,ref);
        sscanf(buffer,"%d",&nRef);
        if(nRef!=grp.nAtoms) {
                sprintf(buffer,"nAtoms in %s (%d) and group %s (%d) don't match!\n",fnRef,nRef,grp.name,grp.nAtoms);
                fatal(buffer);
        }
        refCrd[0]=(t_vecR*)save_malloc(grp.nAtoms*sizeof(t_vecR));
        refCOM[0].x=0.0; refCOM[0].y=0.0; refCOM[0].z=0.0;
        for(i=0;i<nRef;i++) {
                fgets(buffer,200,ref);
                sscanf(&buffer[20],"%f %f %f",&xr,&yr,&zr);
                m=atoms[grp.atoms[i]].mass;
                refCrd[0][i].x=10.0*xr;
                refCrd[0][i].y=10.0*yr;
                refCrd[0][i].z=10.0*zr;
                refCOM[0].x+=m*10.0*xr;
                refCOM[0].y+=m*10.0*yr;
                refCOM[0].z+=m*10.0*zr;
                totMass+=m;
        }
        fclose(ref);
        refCOM[0].x/=totMass;
        refCOM[0].y/=totMass;
        refCOM[0].z/=totMass;

        for(i=0;i<nRef;i++)
        {
                refCrd[0][i].x-=refCOM[0].x;
                refCrd[0][i].y-=refCOM[0].y;
                refCrd[0][i].z-=refCOM[0].z;
        }

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
        FILE  *tin,*cin,*vin,*job,*io;
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
        int needVelo=0;
        int veloFreq=1;
        int wc=0;
        int i=0;
        int j=0;
        int k,l,m,n;
        int nRead=5000;
        int nTimesAnalysis=1000;
        int analysisInterval=10;
        int nSample=10;
        int alignGrp,analyzeGrp;
        int wrap;
        rvec *x;
        real *w_rls;
        rvec *refCrd;;
        t_vecR refCOM;
        real *timebuffer;
        t_box *boxbuffer;
        t_atom **atomsbuffer;
        t_vecR **siteCrdsbuffer;
        int bufferpos;
        char fnRefAlign[100];
/*MODULAR*/
        t_vecR *refCrdDisp;
        t_vecR refCOMDisp;
        char fnRefDisp[100];
        char fnEigVec[100];
        int nDOF;
        int nTimes;
        int nCorr,n1,n2;
        double *sqm;
        double **eigVec;
        int nProj;
        double **proj;
        double tmp1,tmp2,tmp3;
        float floatRead;
        int freqSel;
        int modeStart,modeEnd;
        int bc1,bc2;
	char line[100];
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
                fnRefAlign,fnRefDisp,&alignGrp,&analyzeGrp,&wrap,
                fnEigVec,fnOut,
                needVelo);
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
        nDOF=3*grps[analyzeGrp].nAtoms;
        nProj=1;
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
        eigVec=(double**)save_malloc(nProj*sizeof(double*));
        for(i=0;i<nProj;i++) {
                eigVec[i]=(double*)save_malloc(nDOF*sizeof(double));
        }
        io=fopen(fnEigVec,"rb");
        
	fgets(line, sizeof(line), io);
    printf("%s", line);    
	char tmp[10], ele[10], resName[10];
	int atomNum, resNum;
	double throw1, throw2;
    char buffer[300];
for(i=0;i<nProj;i++) {             
	    for(j=0;j<nDOF;j+=3) {
			fgets(buffer,300,io);
            sscanf(buffer, "%s%d%s%s%d %lf %lf %lf %lf %lf", tmp, &atomNum, ele, resName, &resNum, &tmp1, &tmp2, &tmp3, &throw1, &throw2);
            eigVec[i][j]=tmp1;
			eigVec[i][j+1]=tmp2;
			eigVec[i][j+2]=tmp3;
}
}
        fclose(io);
	
        proj=(double**)save_malloc((nProj+1)*sizeof(double*));
        for(i=0;i<nProj+1;i++) {
                proj[i]=(double*)save_malloc(nTimes*sizeof(double));
                for(j=0;j<nTimes;j++) {
                        proj[i][j]=0.0;
                }
        }
        l=0;
/*MODULAR*/

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

                if(i==0 && alignGrp!=-1) {
                        getAlignRefSpec(fnRefAlign,grps[alignGrp],atoms,&w_rls,&refCrd,&refCOM,&x);
                }
                if(i==0) {
                        getDispRefSpec(fnRefDisp,grps[analyzeGrp],atoms,&refCrdDisp,&refCOMDisp);
                }
                if(i%analysisInterval==0) {
                        if(alignGrp!=-1) {
                                alignGroup(&grps[alignGrp],atoms,nAtoms,&box,w_rls,refCrd,refCOM,x);
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
                        for(j=0;j<grps[analyzeGrp].nAtoms;j++) {
                                m=j*3;
                                k=grps[analyzeGrp].atoms[j];
                                tmp1=sqm[m  ]*(atoms[k].crd.x-refCrdDisp[j].x);
                                tmp2=sqm[m+1]*(atoms[k].crd.y-refCrdDisp[j].y);
                                tmp3=sqm[m+2]*(atoms[k].crd.z-refCrdDisp[j].z);
                                #pragma omp parallel for
                                for(n=0;n<nProj;n++) {
                                        proj[n][l]+=tmp1*eigVec[n][m  ];
                                        proj[n][l]+=tmp2*eigVec[n][m+1];
                                        proj[n][l]+=tmp3*eigVec[n][m+2];
                                }
                        }
                        proj[nProj][l]=time;
                        l++;
/*MODULAR*/
                }
                printf("step %d of %d%c",i+1,nRead,(char)13); fflush(stdout);
        }
        printf("\n");
        closeInput(format,needVelo,cin,vin,&xdrin);
        printf(" - closed input trajectory\n");fflush(stdout);

        sprintf(fnOut2,"proj_freq1_mode7.dat",fnOut);
        printf(" - writing eigenvector projection of weighted displacements in file: %s\n",fnOut2);
        io=fopen(fnOut2,"w");
        /*printing times*/
        for(j=0;j<nTimes;j++) {
                fprintf(io," %e",proj[nProj][j]);
        }
        fprintf(io,"\n");
        for(i=0;i<nProj;i++) {
                for(j=0;j<nTimes;j++) {
                        fprintf(io," %e",proj[i][j]);
                }
                fprintf(io,"\n");
        }
        fclose(io);
        printf(" - closed file: %s\n",fnOut2);

        return 0;
}

