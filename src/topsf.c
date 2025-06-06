#include <stdio.h>
#include <string.h>
#include "../include/fatal.h"
#include "../include/dataTypes.h"
#include "../include/alloc.h"
#include "../include/io.h"

typedef struct {
        int     atomNr;
        int     resNr;
        char    resName[6];
        char    atomName[6];
        char    atomType[12];
        float   charge;
        float   mass;
} t_atomLine;

typedef struct {
	char name[100];
	int nMol;
	int nRes;
	int *resNAtoms;
	int nAtoms;
	t_atomLine *atoms;
	int nBonds;
	/*the last element is used as a tag*/
	int (*bonds)[3];
	int nAngles;
	int (*angles)[4];
	int nDih;
	int (*dih)[5];
	int nImp;
	int (*imp)[5];
	int cnt;
} t_itp;

int initITP(t_itp *itpInfo) {
	itpInfo->nMol=0;
	itpInfo->nRes=0;
	itpInfo->nAtoms=0;
	itpInfo->nBonds=0;
	itpInfo->nAngles=0;
	itpInfo->nDih=0;
	itpInfo->nImp=0;
	itpInfo->cnt=0;
	return 0;
}

int getLine(FILE *in,char *buffer,int *lineCnt)
{
	int eof=0;

        if(fgets(buffer,300,in)==NULL) eof=1;
	else {
        	lineCnt[0]++;
        	while(eof==0 && (strncmp(buffer,";",1)==0 || strncmp(buffer,"#",1)==0))
       		{
        	        if(fgets(buffer,300,in)==NULL) eof=1;
        	        lineCnt[0]++;
        	}
	}
        return eof;
}

int goToSection(char *fn,FILE *in,char *keyword, int *lineCnt,int safe)
{
        char buffer[300];
	int eof=0;

	eof=getLine(in,buffer,lineCnt);
        while(eof==0 && strstr(buffer,keyword)==NULL) eof=getLine(in,buffer,lineCnt);
	if(eof==1 && safe==1) {
		sprintf(buffer,"section '%s' not found in file %s\n",keyword,fn);
		fatal(buffer);
	}
        return eof;
}

int main(int argc,char *argv[]) {
	char buffer[300],tmp[300];
	FILE *TOP,*ITP,*FFnb,*OUT;
	char fnTOP[300],fnFFnb[300],fnOUT[300],(*fnITP)[300];
	char mol[100];
	int lineCntTOP,lineCntITP,lineCntFFnb;
	int ITPcnt=0;
	int nITP;
	t_itp *itpInfo;
	int i,j,k,l,m,n;
	char dum1[100],dum2[100],dum3[100],dum4[100],dum5[100],dum6[100];
	int resCur,resLast,cgnr;
	char chains[26][2] = {"A","B","C","D","E","F","G","H","I","J","K","L","M",
		"N","O","P","Q","R","S","T","U","V","W","X","Y","Z"};
	int chainCnt=0;
	int nAtoms,nMols,nMolTypes;
	int start,end;
	char **types;
	float *sigmaLJ,*epsLJ;
	int nTypes;
	int OPLS=-1;
	char ffAtomType[12];
	int foundType;

	/*first we count the number of molecule types specified in the GMX top file*/
	printf("GMX topology filename: ");
	fgets(buffer,300,stdin);
	sscanf(buffer,"%s",fnTOP);
	saveOpenRead(&TOP,fnTOP);lineCntTOP=0;
	goToSection(fnTOP,TOP,"[ molecules ]",&lineCntTOP,1);
	printf("found section '[ molecules ]' in line %d of file %s\n",lineCntTOP,fnTOP);
	printf("molecules in system:\n");
	while(getLine(TOP,buffer,&lineCntTOP)==0) {
		sscanf(buffer,"%s %d",mol,&i);
		printf("%-20s%10d\n",mol,i);
		ITPcnt++;
	}
	nITP=ITPcnt;ITPcnt=0;
	if(nITP>26) fatal("SORRY: maximum of 26 itp files (chains) currently supported\n");
	fnITP=save_malloc(nITP*sizeof(*fnITP));
	itpInfo=(t_itp*)save_malloc(nITP*sizeof(t_itp));
	
	rewind(TOP);lineCntTOP=0;

	/*now we re-read the moleculetypes and open the cooresponding itp file*/
	/*we will go thorugh each itp file first to count the number of
		residues
		atoms
		bonds
		angles
		dihedrals
		impropers
	a bit cumbersome, but that way we can check how much memory we will need*/
	goToSection(fnTOP,TOP,"[ molecules ]",&lineCntTOP,1);
	while(getLine(TOP,buffer,&lineCntTOP)==0) {
		initITP(&itpInfo[ITPcnt]);
		sscanf(buffer,"%s %d",itpInfo[ITPcnt].name,&itpInfo[ITPcnt].nMol);
		printf("top/itp filename for molecule %s: ",itpInfo[ITPcnt].name);
		fgets(buffer,300,stdin);
		sscanf(buffer,"%s",fnITP[ITPcnt]);
		saveOpenRead(&ITP,fnITP[ITPcnt]);lineCntITP=0;
		goToSection(fnITP[ITPcnt],ITP,"[ moleculetype ]",&lineCntITP,1);
		getLine(ITP,buffer,&lineCntITP);
		sscanf(buffer,"%s",tmp);
		while(strcmp(tmp,itpInfo[ITPcnt].name)!=0) {
			goToSection(fnITP[ITPcnt],ITP,"[ moleculetype ]",&lineCntITP,1);
			getLine(ITP,buffer,&lineCntITP);
			sscanf(buffer,"%s",tmp);
		}
		printf("-found moleculetype '%s' in line %d of file %s\n",
			itpInfo[ITPcnt].name,lineCntITP,fnITP[ITPcnt]);
		goToSection(fnITP[ITPcnt],ITP,"[ atoms ]",&lineCntITP,1);
		resLast=-1;
		while(getLine(ITP,buffer,&lineCntITP)==0 && strcmp(buffer,"\n")!=0 && strcmp(buffer,"\t\n")!=0 
			&& strcmp(buffer," \n")!=0 && strcmp(buffer,"  \n")!=0 && strncmp(buffer,"[",1)!=0) {
			itpInfo[ITPcnt].nAtoms++;
			sscanf(buffer,"%s %s %d",dum1,dum2,&resCur);
			if(resCur!=resLast) {
				itpInfo[ITPcnt].nRes++;
				resLast=resCur;
			}
		}
		printf("--found %d residues for moleculetype %s\n",itpInfo[ITPcnt].nRes,itpInfo[ITPcnt].name);
		printf("--found %d atoms for moleculetype %s\n",itpInfo[ITPcnt].nAtoms,itpInfo[ITPcnt].name);
		itpInfo[ITPcnt].atoms=(t_atomLine*)save_malloc(itpInfo[ITPcnt].nAtoms*sizeof(t_atomLine));
		itpInfo[ITPcnt].resNAtoms=(int*)save_malloc(itpInfo[ITPcnt].nRes*sizeof(int));
		rewind(ITP);
		/*go back to atom section of the target moleuletype*/
		goToSection(fnITP[ITPcnt],ITP,"[ moleculetype ]",&lineCntITP,1);
                getLine(ITP,buffer,&lineCntITP);
                sscanf(buffer,"%s",tmp);
                while(strcmp(tmp,itpInfo[ITPcnt].name)!=0) {
                        goToSection(fnITP[ITPcnt],ITP,"[ moleculetype ]",&lineCntITP,1);
                        getLine(ITP,buffer,&lineCntITP);
                        sscanf(buffer,"%s",tmp);
                }
                goToSection(fnITP[ITPcnt],ITP,"[ atoms ]",&lineCntITP,1);
		/*now read atom section again to count atoms per residue*/
                resLast=-1;
		i=-1;
		while(getLine(ITP,buffer,&lineCntITP)==0 && strcmp(buffer,"\n")!=0 && strcmp(buffer,"\t\n")!=0 
			&& strcmp(buffer," \n")!=0 && strcmp(buffer,"  \n")!=0 && strncmp(buffer,"[",1)!=0) {
                        sscanf(buffer,"%s %s %d",dum1,dum2,&resCur);
                        if(resCur!=resLast) {
                                resLast=resCur;
				i++;
				itpInfo[ITPcnt].resNAtoms[i]=1;
                        } else {
				itpInfo[ITPcnt].resNAtoms[i]++;	
			}
                }

		if(goToSection(fnITP[ITPcnt],ITP,"[ bonds ]",&lineCntITP,0)==0) {
			while(getLine(ITP,buffer,&lineCntITP)==0 && strcmp(buffer,"\n")!=0 && strcmp(buffer,"\t\n")!=0 
				&& strcmp(buffer," \n")!=0 && strcmp(buffer,"  \n")!=0 && strncmp(buffer,"[",1)!=0) {
				itpInfo[ITPcnt].nBonds++;
			}
		}
		printf("--found %d bonds for moleculetype %s\n",itpInfo[ITPcnt].nBonds,itpInfo[ITPcnt].name);
		itpInfo[ITPcnt].bonds=save_malloc(itpInfo[ITPcnt].nBonds*sizeof(*itpInfo[ITPcnt].bonds));

		if(goToSection(fnITP[ITPcnt],ITP,"[ angles ]",&lineCntITP,0)==0) {
                        while(getLine(ITP,buffer,&lineCntITP)==0 && strcmp(buffer,"\n")!=0 && strcmp(buffer,"\t\n")!=0 
				&& strcmp(buffer," \n")!=0 && strcmp(buffer,"  \n")!=0 && strncmp(buffer,"[",1)!=0) {
                                itpInfo[ITPcnt].nAngles++;
                        }
                }
                printf("--found %d angles for moleculetype %s\n",itpInfo[ITPcnt].nAngles,itpInfo[ITPcnt].name);
		itpInfo[ITPcnt].angles=save_malloc(itpInfo[ITPcnt].nAngles*sizeof(*itpInfo[ITPcnt].angles));

		if(goToSection(fnITP[ITPcnt],ITP,"[ dihedrals ]",&lineCntITP,0)==0) {
                        while(getLine(ITP,buffer,&lineCntITP)==0 && strcmp(buffer,"\n")!=0 && strcmp(buffer,"\t\n")!=0 
				&& strcmp(buffer," \n")!=0 && strcmp(buffer,"  \n")!=0 && strncmp(buffer,"[",1)!=0) {
                                itpInfo[ITPcnt].nDih++;
                        }
                }
                printf("--found %d dihedrals for moleculetype %s\n",itpInfo[ITPcnt].nDih,itpInfo[ITPcnt].name);
		itpInfo[ITPcnt].dih=save_malloc(itpInfo[ITPcnt].nDih*sizeof(*itpInfo[ITPcnt].dih));

		/*gromacs calls the impoper section dihedrals too, so this is not a typo*/
		if(goToSection(fnITP[ITPcnt],ITP,"[ dihedrals ]",&lineCntITP,0)==0) {
                        while(getLine(ITP,buffer,&lineCntITP)==0 && strcmp(buffer,"\n")!=0 && strcmp(buffer,"\t\n")!=0 
				&& strcmp(buffer," \n")!=0 && strcmp(buffer,"  \n")!=0 && strncmp(buffer,"[",1)!=0) {
                                itpInfo[ITPcnt].nImp++;
                        }
                }
                printf("--found %d impropers for moleculetype %s\n",itpInfo[ITPcnt].nImp,itpInfo[ITPcnt].name);
		itpInfo[ITPcnt].imp=save_malloc(itpInfo[ITPcnt].nImp*sizeof(*itpInfo[ITPcnt].imp));
	
		fclose(ITP);
		ITPcnt++;
	}
	fclose(TOP);

	/*now we will re-read the itp files and store the data in allocated arrays*/
	for(ITPcnt=0;ITPcnt<nITP;ITPcnt++) {
		saveOpenRead(&ITP,fnITP[ITPcnt]);lineCntITP=0;
                goToSection(fnITP[ITPcnt],ITP,"[ moleculetype ]",&lineCntITP,1);
                getLine(ITP,buffer,&lineCntITP);
                sscanf(buffer,"%s",tmp);
                while(strcmp(tmp,itpInfo[ITPcnt].name)!=0) {
                        goToSection(fnITP[ITPcnt],ITP,"[ moleculetype ]",&lineCntITP,1);
                        getLine(ITP,buffer,&lineCntITP);
                        sscanf(buffer,"%s",tmp);
                }
		goToSection(fnITP[ITPcnt],ITP,"[ atoms ]",&lineCntITP,1);
		for(i=0;i<itpInfo[ITPcnt].nAtoms;i++) {
			getLine(ITP,buffer,&lineCntITP);
			j=sscanf(buffer,"%d %s %d %s %s %d %f %f",
				&itpInfo[ITPcnt].atoms[i].atomNr,
				itpInfo[ITPcnt].atoms[i].atomType,
				&itpInfo[ITPcnt].atoms[i].resNr,
				itpInfo[ITPcnt].atoms[i].resName,
				itpInfo[ITPcnt].atoms[i].atomName,
				&cgnr,
				&itpInfo[ITPcnt].atoms[i].charge,
				&itpInfo[ITPcnt].atoms[i].mass
			);
			if(j<7) {
				sprintf(buffer,"error in reading atom info in line %d of file %s\n",
					lineCntITP,fnITP[ITPcnt]);
				fatal(buffer);
			}
			if(j!=8) {
				if(strncmp(itpInfo[ITPcnt].atoms[i].atomType,"H",1)==0) {
					itpInfo[ITPcnt].atoms[i].mass=1.008;
				} else if(strncmp(itpInfo[ITPcnt].atoms[i].atomType,"C",1)==0) {
                                        itpInfo[ITPcnt].atoms[i].mass=12.011;
                                } else if(strncmp(itpInfo[ITPcnt].atoms[i].atomType,"O",1)==0) {
                                        itpInfo[ITPcnt].atoms[i].mass=15.9994;
                                } else if(strncmp(itpInfo[ITPcnt].atoms[i].atomType,"N",1)==0) {
                                        itpInfo[ITPcnt].atoms[i].mass=14.007;
                                } else if(strncmp(itpInfo[ITPcnt].atoms[i].atomType,"S",1)==0) {
                                        itpInfo[ITPcnt].atoms[i].mass=32.06;
                                } else itpInfo[ITPcnt].atoms[i].mass=0.0;
			}
			/*printf("%d %s %d %s %s %f %f\n",
                                itpInfo[ITPcnt].atoms[i].atomNr,
                                itpInfo[ITPcnt].atoms[i].atomType,
                                itpInfo[ITPcnt].atoms[i].resNr,
                                itpInfo[ITPcnt].atoms[i].resName,
                                itpInfo[ITPcnt].atoms[i].atomName,
                                itpInfo[ITPcnt].atoms[i].charge,
                                itpInfo[ITPcnt].atoms[i].mass
                        );*/
		}
		if(itpInfo[ITPcnt].nBonds!=0) {
			goToSection(fnITP[ITPcnt],ITP,"[ bonds ]",&lineCntITP,1);
			for(i=0;i<itpInfo[ITPcnt].nBonds;i++) {
                        	getLine(ITP,buffer,&lineCntITP);
                        	j=sscanf(buffer,"%d %d",
					&itpInfo[ITPcnt].bonds[i][0],
					&itpInfo[ITPcnt].bonds[i][1]
				);
				if(j!=2) {
                                	sprintf(buffer,"error in reading bond info in line %d of file %s\n",
                                        	lineCntITP,fnITP[ITPcnt]);
                                	fatal(buffer);
                        	}
				/*printf("%d %d\n",
					itpInfo[ITPcnt].bonds[i][0],
					itpInfo[ITPcnt].bonds[i][1]
				);*/
			}
		}
		if(itpInfo[ITPcnt].nAngles!=0) {
                        goToSection(fnITP[ITPcnt],ITP,"[ angles ]",&lineCntITP,1);
                        for(i=0;i<itpInfo[ITPcnt].nAngles;i++) {
                                getLine(ITP,buffer,&lineCntITP);
                                j=sscanf(buffer,"%d %d %d",
					&itpInfo[ITPcnt].angles[i][0],
					&itpInfo[ITPcnt].angles[i][2],
					&itpInfo[ITPcnt].angles[i][1]
				);
				if(j!=3) {
                                        sprintf(buffer,"error in reading angle info in line %d of file %s\n",
                                                lineCntITP,fnITP[ITPcnt]);
                                        fatal(buffer);
                                }
                                /*printf("%d %d %d\n",
					itpInfo[ITPcnt].angles[i][0],
					itpInfo[ITPcnt].angles[i][1],
					itpInfo[ITPcnt].angles[i][2]
				);*/
                        }
                }
		if(itpInfo[ITPcnt].nDih!=0) {
                        goToSection(fnITP[ITPcnt],ITP,"[ dihedrals ]",&lineCntITP,1);
                        for(i=0;i<itpInfo[ITPcnt].nDih;i++) {
                                getLine(ITP,buffer,&lineCntITP);
                                j=sscanf(buffer,"%d %d %d %d",
                                        &itpInfo[ITPcnt].dih[i][0],
                                        &itpInfo[ITPcnt].dih[i][1],
                                        &itpInfo[ITPcnt].dih[i][2],
					&itpInfo[ITPcnt].dih[i][3]
                                );
				if(j!=4) {
                                        sprintf(buffer,"error in reading dihedral info in line %d of file %s\n",
                                                lineCntITP,fnITP[ITPcnt]);
                                        fatal(buffer);
                                }
                                /*printf("%d %d %d %d\n",
                                        itpInfo[ITPcnt].dih[i][0],
                                        itpInfo[ITPcnt].dih[i][1],
                                        itpInfo[ITPcnt].dih[i][2],
					itpInfo[ITPcnt].dih[i][3]
                                );*/
                        }
                }
		if(itpInfo[ITPcnt].nImp!=0) {
                        goToSection(fnITP[ITPcnt],ITP,"[ dihedrals ]",&lineCntITP,1);
                        for(i=0;i<itpInfo[ITPcnt].nImp;i++) {
                                getLine(ITP,buffer,&lineCntITP);
                                j=sscanf(buffer,"%d %d %d %d",
                                        &itpInfo[ITPcnt].imp[i][0],
                                        &itpInfo[ITPcnt].imp[i][1],
                                        &itpInfo[ITPcnt].imp[i][2],
                                        &itpInfo[ITPcnt].imp[i][3]
                                );
				if(j!=4) {
                                        sprintf(buffer,"error in reading improper info in line %d of file %s\n",
                                                lineCntITP,fnITP[ITPcnt]);
                                        fatal(buffer);
                                }
                                /*printf("%d %d %d %d\n",
                                        itpInfo[ITPcnt].imp[i][0],
                                        itpInfo[ITPcnt].imp[i][1],
                                        itpInfo[ITPcnt].imp[i][2],
                                        itpInfo[ITPcnt].imp[i][3]
                                );*/
                        }
                }
		fclose(ITP);
	}

	/*now we sort out the atom types present in the system*/
	nTypes=0;
	for(i=0;i<nITP;i++) nTypes+=itpInfo[i].nAtoms;
	types=(char**)save_malloc(nTypes*sizeof(char*));
	sigmaLJ=(real*)save_malloc(nTypes*sizeof(real));
	epsLJ=(real*)save_malloc(nTypes*sizeof(real));
	for(i=0;i<nTypes;i++) types[i]=(char*)save_malloc(12*sizeof(char));
	n=0;
	for(i=0;i<nITP;i++) {
		for(j=0;j<itpInfo[i].nAtoms;j++) {
			if(n==0) {
				strcpy(types[n],itpInfo[i].atoms[j].atomType);
				n++;
			} else {
				/*move backward through the list of atom types 
				to see if we encountered this atomType before*/
				m=n-1;
				while(m>=0 && strcmp(itpInfo[i].atoms[j].atomType,types[m])!=0) m--;
				/*current atom type not found list => append*/
				if(m<0) {
					strcpy(types[n],itpInfo[i].atoms[j].atomType);
					n++;
				}
			}
		}
	}
	nTypes=n;
	printf("found %d atom types\n",nTypes);
	for(i=0;i<nTypes;i++) {
		printf("%12s ",types[i]);
		if((i+1)%5==0) printf("\n");
	}
	printf("\n");
	/*now we have to find the parameters in a nonbonded force field parameter file*/
	printf("requiring epislon and sigma in kJ/(mol*nm^2) and nm\n");
	printf("GMX non-bonded ff file: ");
	fgets(buffer,300,stdin);
	sscanf(buffer,"%s",fnFFnb);
	saveOpenRead(&FFnb,fnFFnb);lineCntFFnb=0;
        goToSection(fnFFnb,FFnb,"[ atomtypes ]",&lineCntFFnb,1);
	for(i=0;i<nTypes;i++) {
		foundType=0;
		while(foundType==0 && getLine(FFnb,buffer,&lineCntFFnb)==0) {
			sscanf(buffer,"%s",ffAtomType);
			if(OPLS==-1) {  /*this means this has not been set yet*/
				if(strncmp("opls",ffAtomType,4)==0) OPLS=1;
				else OPLS=0;
			}
			if(strcmp(ffAtomType,types[i])==0) foundType=1;
		}
		if(foundType==0) {
			sprintf(buffer,"did not find atom type %s in file %s\n",types[i],fnFFnb);
			fatal(buffer);
		}
		if(OPLS==1) sscanf(buffer,"%s %s %s %s %s %s %f %f",dum1,dum2,dum3,dum4,dum5,dum6,&sigmaLJ[i],&epsLJ[i]);
		else sscanf(buffer,"%s %s %s %s %s %f %f",dum1,dum2,dum3,dum4,dum5,&sigmaLJ[i],&epsLJ[i]);
		/*convert sigma from 'nm' to 'A'*/
		sigmaLJ[i]=10.0*sigmaLJ[i];
		printf("%12s %15.5e %15.5e\n",types[i],sigmaLJ[i],epsLJ[i]);
		rewind(FFnb);lineCntFFnb=0;
		goToSection(fnFFnb,FFnb,"[ atomtypes ]",&lineCntFFnb,1);
	}
	fclose(FFnb);

	/*now we assemble the data to write out the converted file*/
        nAtoms=0;
        for(i=0;i<nITP;i++) {
                nAtoms+=itpInfo[i].nMol*itpInfo[i].nAtoms;
        }
        nMols=0;
        for(i=0;i<nITP;i++) {
                nMols+=itpInfo[i].nMol*itpInfo[i].nRes;
        }
        nMolTypes=0;
        for(i=0;i<nITP;i++) {
                if(itpInfo[i].nRes==1) {
                        nMolTypes++;
                } else {
                        nMolTypes+=itpInfo[i].nMol*itpInfo[i].nRes;
                }
        }
        printf("SUMMARY:\n");
        printf("nAtoms:     %d\n",nAtoms);
        printf("nMols:      %d\n",nMols);
        printf("nMolTypes:  %d\n",nMolTypes);
	printf("output file: ");
	fgets(buffer,300,stdin);
	sscanf(buffer,"%s",fnOUT);
	OUT=fopen(fnOUT,"w");
	fprintf(OUT,"PSF CMAP\n\n");
	fprintf(OUT," %7d !NTITLE\n",3);
	fprintf(OUT,"REMARKS generated from Gromacs topology\n");
	fprintf(OUT,"REMARKS don't use for simulations, analyis only\n");
	fprintf(OUT,"REMARKS several terms might be missing\n\n");
	fprintf(OUT," %7d !NATOM\n",nAtoms);
	n=1;
	for(i=0;i<nITP;i++) {
		for(j=0;j<itpInfo[i].nMol;j++) {
			for(k=0;k<itpInfo[i].nAtoms;k++) {
				fprintf(OUT," %7d %s %5d %6s  %-4s %-4s %10.6f  %12.4f  %10d\n",
					n,
					chains[i],
					itpInfo[i].atoms[k].resNr+(j*itpInfo[i].nRes),
					itpInfo[i].atoms[k].resName,
					itpInfo[i].atoms[k].atomName,
					itpInfo[i].atoms[k].atomType,
					itpInfo[i].atoms[k].charge,
					itpInfo[i].atoms[k].mass,
					0
				);
				n++;
			}
		}
	}
	fprintf(OUT,"\n");
	n=0;
	for(i=0;i<nITP;i++) {
		n+=itpInfo[i].nMol*itpInfo[i].nBonds;
	}
	fprintf(OUT," %7d !NBOND: bonds\n",n);
	n=0;
	m=0;
	for(i=0;i<nITP;i++) {
		for(j=0;j<itpInfo[i].nMol;j++) {
			for(k=0;k<itpInfo[i].nBonds;k++) {
				fprintf(OUT," %7d %7d",
					itpInfo[i].bonds[k][0]+n,
					itpInfo[i].bonds[k][1]+n
				);
				m++;
				if(m%4==0) fprintf(OUT,"\n");
			}
			n+=itpInfo[i].nAtoms;
		}
	}
	if(m%4!=0) fprintf(OUT,"\n");
	fprintf(OUT,"\n");

	n=0;
        for(i=0;i<nITP;i++) {
                n+=itpInfo[i].nMol*itpInfo[i].nAngles;
        }
        fprintf(OUT," %7d !NTHETA: angles\n",n);
	n=0;
        m=0;
        for(i=0;i<nITP;i++) {
                for(j=0;j<itpInfo[i].nMol;j++) {
                        for(k=0;k<itpInfo[i].nAngles;k++) {
                                fprintf(OUT," %7d %7d %7d",
					itpInfo[i].angles[k][0]+n,
					itpInfo[i].angles[k][2]+n,
					itpInfo[i].angles[k][1]+n
				);
                                m++;
                                if(m%3==0) fprintf(OUT,"\n");
                        }
                        n+=itpInfo[i].nAtoms;
                }
        }
        if(m%3!=0) fprintf(OUT,"\n");
        fprintf(OUT,"\n");

	n=0;
        for(i=0;i<nITP;i++) {
                n+=itpInfo[i].nMol*itpInfo[i].nDih;
        }
	fprintf(OUT," %7d !NPHI: dihedrals\n",n);
	n=0;
	m=0;
	for(i=0;i<nITP;i++) {
	        for(j=0;j<itpInfo[i].nMol;j++) {
	                for(k=0;k<itpInfo[i].nDih;k++) {
	                        fprintf(OUT," %7d %7d %7d %7d",
	                                itpInfo[i].dih[k][0]+n,
	                                itpInfo[i].dih[k][1]+n,
	                                itpInfo[i].dih[k][2]+n,
					itpInfo[i].dih[k][3]+n
	                        );
	                        m++;
	                        if(m%2==0) fprintf(OUT,"\n");
	                }
	                n+=itpInfo[i].nAtoms;
	        }
	}
	if(m%2!=0) fprintf(OUT,"\n");
	fprintf(OUT,"\n");

	n=0;
	for(i=0;i<nITP;i++) {
	        n+=itpInfo[i].nMol*itpInfo[i].nImp;
	}
	fprintf(OUT," %7d !NIMPHI: impropers\n",n);
	n=0;
	m=0;
	for(i=0;i<nITP;i++) {
	        for(j=0;j<itpInfo[i].nMol;j++) {
	                for(k=0;k<itpInfo[i].nImp;k++) {
	                        fprintf(OUT," %7d %7d %7d %7d",
	                                itpInfo[i].imp[k][0]+n,
	                                itpInfo[i].imp[k][1]+n,
	                                itpInfo[i].imp[k][2]+n,
	                                itpInfo[i].imp[k][3]+n
	                        );
	                        m++;
	                        if(m%2==0) fprintf(OUT,"\n");
	                }
	                n+=itpInfo[i].nAtoms;
	        }
	}
	if(m%2!=0) fprintf(OUT,"\n");
	fprintf(OUT,"\n");

	fprintf(OUT," %7d !NDON: donors\n\n",0);
	fprintf(OUT," %7d !NACC: acceptors\n\n",0);
	fprintf(OUT," %7d !NNB: donors\n\n",0);
	for(i=0;i<nAtoms;i++) {
		fprintf(OUT," %7d",0);
		if((i+1)%8==0) fprintf(OUT,"\n");
	}
	if(i%8!=0) fprintf(OUT,"\n");
	fprintf(OUT,"\n");
	fprintf(OUT," %7d %7d !NGRP\n %7d %7d %7d\n\n",0,0,0,0,0);

	fprintf(OUT," %7d !NCRTERM: cross-terms\n\n",0);
	
	fclose(OUT);
	return 0;
}
