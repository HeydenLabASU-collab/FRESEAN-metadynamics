#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
	printf("THE PURPOSE OF THIS PROGRAM IS TO GENERATE\n");
	printf("A COARSE-GRAINED REPRESENTATION OF AN\n");
	printf("ALL -ATOM TRAJECTORY.\n");
	printf("THE COARSE-GRAINED REPRESENTATION IS\n");
	printf("DEFINED AS TWO CENTER OF MASS BEADS FOR.\n");
	printf("EACH AMINO ACID (ONE BACKBONE/ONE SIDECHAIN).\n");
	printf("Version 1: November 1, 2023\n");
	printf("Author:\n");
	printf("Michael A. Sauer\n");
	printf("School of Molecular Sciences\n");
	printf("Arizona State University\n");
	printf("Tempe, AZ, USA\n");
	printf("e-mail: masauer2@asu.edu\n");
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
        printf("fnTop (topology file [.mtop]) \nfnCrd\nfnVel (if format xyz,crd,dcd)\n");
        printf("fnJob (provided [.job] file)\analyzeGrp (atom group in job file to perform analysis on)\n");
		printf("nRead (Number of frames to read)\nnSample (dFrame)\nfnOutTraj (Output coarse-grain trajectory [.gro])\n");
		printf("fnOutTopol (Output coarse-grain topology [.mtop])\n");
        return 0;
}

int getInput(char *fnCOM,char *fnTop,char *fnCrd,char *fnVel, char *fnJob, int *analyzeGrp, 
                int *nRead, int *nSample, 
                char *fnOutTraj, char *fnOutTopol, int needVelo) {
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
        sscanf(buffer,"%s", fnJob);
        printf("%2d -> read %20s : %s\n",i,"fnJob",fnJob);i++;

		getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d", analyzeGrp);
        printf("%2d -> read %20s : %d\n",i,"analyzeGrp",analyzeGrp[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nRead);
        printf("%2d -> read %20s : %d\n",i,"nRead",nRead[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nSample);
        printf("%2d -> read %20s : %d\n",i,"nSample",nSample[0]);i++;;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnOutTraj);
        printf("%2d -> read %20s : %s\n",i,"fnOutTraj",fnOutTraj);i++;

		getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnOutTopol);
        printf("%2d -> read %20s : %s\n",i,"fnOutTopol",fnOutTopol);i++;

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


// Print the backbone bead to the file
void printBackbone(FILE *filewriter, t_atom atom, real *COMPOS, real *COMVEL, int *beadNum){
	fprintf(filewriter, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", atom.resNr + 1, atom.resName, "BACK", *beadNum, COMPOS[0]/10, COMPOS[1]/10, COMPOS[2]/10, COMVEL[0]/10, COMVEL[1]/10, COMVEL[2]/10);
	(*beadNum)++;
}

// Print the sidechain bead to the file
void printSidechain(FILE *filewriter, t_atom atom, real *COMPOS, real *COMVEL, int *beadNum){
	fprintf(filewriter, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", atom.resNr + 1, atom.resName, "SIDE", *beadNum, COMPOS[0]/10, COMPOS[1]/10, COMPOS[2]/10, COMVEL[0]/10, COMVEL[1]/10, COMVEL[2]/10);
	(*beadNum)++;
}

// Take in the position and velocities of a single atom and
// add to the center of mass of the bead.
// Add atom mass to mass of bead.
int assignCOM(real *COMPOS, real *COMVEL, real *MASS, t_atom atom){
	COMPOS[0] += atom.mass * atom.crd.x;
	COMPOS[1] += atom.mass * atom.crd.y;
	COMPOS[2] += atom.mass * atom.crd.z;
	COMVEL[0] += atom.mass * atom.vel.x;
	COMVEL[1] += atom.mass * atom.vel.y;
	COMVEL[2] += atom.mass * atom.vel.z;
	*MASS += atom.mass;
}

// Divide total position and velocity by the mass
int divideCOM(real *COMPOS, real *COMVEL, real MASS){
	for(int i = 0; i < 3; i++){
		COMPOS[i] /= MASS;
		COMVEL[i] /= MASS;
	}
}

// We will write the bead to the file and then reset the values
// for the next bead.
int resetCOM(real *COMPOS, real *COMVEL, real *MASS){
	for(int i = 0; i < 3; i++){
		COMPOS[i] = 0.0;
		COMVEL[i] = 0.0;
	}
	*MASS = 0.0;
}

// If we have a capping residue or a glycine residue, we do
// not need to consider a sidechain bead.
// NOTE: CHECK TOPOLOGY FOR OTHER ATOMS THAT HAVE THIS PROPERTY
int isBackboneOnly(char *resName){
	if(strcmp(resName, "ACE") == 0 || strcmp(resName, "NME") == 0 || strcmp(resName, "GLY") == 0){
		return 1;
	}
	else{
		return 0;
	}
	
}

// Check if the residue is C or N terminus
// NOTE: CURRENTLY UNUSED BUT COULD BE USEFUL IN THE FUTURE
int isNorCterminus(int resNum, int lastRes){
	if(resNum == 0){
		return 1;
	}
	else if(resNum == lastRes){
		return 2;
	}
	else{
		return 0;
	}
	
}

// Check to see if the current atom is in the backbone of the amino acid
// We will use the names from the topology.
// This includes differences in the backbone atoms when we have a c or n 
// terminus. 
int inBackbone(t_atom atom){
	if(strcmp(atom.atomName, "CA") == 0 || 
	   strcmp(atom.atomName, "OC1") == 0 || 
	   strcmp(atom.atomName, "OC2") == 0 ||
	   strcmp(atom.atomName, "C") == 0 || 
	   strcmp(atom.atomName, "O") == 0 || 
	   strcmp(atom.atomName, "HA") == 0 || 
	   strcmp(atom.atomName, "N") == 0 || 
	   strcmp(atom.atomName, "H") == 0 || 
	   strcmp(atom.atomName, "H1") == 0|| 
	   strcmp(atom.atomName, "H2") == 0|| 
	   strcmp(atom.atomName, "H3") == 0){
	   return 1;
	}
	else{
		return 0;
	}
	
}
int printReferenceFile(t_atom *atom, int nAtoms, t_grp grp){
	int nRes = 0;
	printf("Num atoms in grp: %d\n", grp.nAtoms);
	for(int j=0; j<grp.nAtoms;j++) {	
		if(j == grp.nAtoms - 1){
			nRes = atom[grp.atoms[j]].resNr;
		}
	}
	nRes = nRes + 1;
	printf("Number of Residues: %d\n", nRes);
	return nRes;
}

// If we are at the first frame, we will read in all of the atoms and coarse grain 
// the first frame seperately.
// Then we will use this to output a custom .mtop file.
// NOTE: Is there a better way to do this? Currently unoptimized but this only happens once.
int printCustomTopology(FILE *topol, t_atom *atom, int nRes, t_grp grp){

		printf("Number of Residues: %d\n", nRes);

		int curr_index = 0;
		int nBeads = 0;

		double BBMASS[nRes];
		double SCMASS[nRes];
		char *NAMES[nRes];
		real COM_POS_BB[3] = {0.0, 0.0, 0.0};
    	real COM_POS_SC[3] = {0.0, 0.0, 0.0};
    	real COM_VEL_BB[3] = {0.0, 0.0, 0.0};
   		real COM_VEL_SC[3] = {0.0, 0.0, 0.0};
    	real BB_MASS = 0.0;
    	real SC_MASS = 0.0;
		printf("Data for custom topology: %d\n", grp.nAtoms);
		for(int j=0;j<grp.nAtoms;j++){	

			curr_index = atom[grp.atoms[j]].resNr;
			printf("--------\nCurrent Index: %d -%d - %s\n", j, curr_index, atom[grp.atoms[j]].resName);
			/*
			Check the residue name to see if this atom is part of a residue with no side chain
			This includes any capping residues or glycine

			If it is one of these, we only need to assign the backbone center of mass position and velocity
			*/
			if(isBackboneOnly(atom[grp.atoms[j]].resName)){
				assignCOM(COM_POS_BB, COM_VEL_BB, &BB_MASS, atom[grp.atoms[j]]);
			}

			/*
			If the current atom contains both a backbone and a sidechain, then we need to identify whether
			the current atom belongs to the sidechain OR to the backbone and update the center of mass
			accordingly
			*/
			else{
				if(inBackbone(atom[grp.atoms[j]])){
					assignCOM(COM_POS_BB, COM_VEL_BB, &BB_MASS, atom[grp.atoms[j]]);
				}
				else{
					assignCOM(COM_POS_SC, COM_VEL_SC, &SC_MASS, atom[grp.atoms[j]]);
				}
			}
			
			/*
			 * Condition to check if we are switching residues or if we are at the last residue
			*/
			printf("Target: %d\nCurrent: %d\n", grp.nAtoms - 1, j);
			if((j == grp.nAtoms - 1) || (atom[grp.atoms[j+1]].resNr != curr_index)){
				printf("Here for iteration %d\n", curr_index);
				divideCOM(COM_POS_BB, COM_VEL_BB, BB_MASS);
				divideCOM(COM_POS_SC, COM_VEL_SC, SC_MASS);
					
				BBMASS[curr_index] = BB_MASS;
				SCMASS[curr_index] = SC_MASS;
				NAMES[curr_index] = atom[grp.atoms[j]].resName;

				resetCOM(COM_POS_BB, COM_VEL_BB, &BB_MASS);
				resetCOM(COM_POS_SC, COM_VEL_SC, &SC_MASS);		
				printf("Current res name of atom at selection %d: %s\n", j, atom[grp.atoms[j]].resName);
			}
		}
		printf("Number of Residues: %d\n", nRes);
		for(int l=0;l<3;l++){	
				printf("Name of arr: %f\n", BBMASS[l]);
		}
		// nBeads = Number of beads in our trajectory.
		for(int k = 0; k < nRes; k++){

			// Is Backbone mass non-zero? If so, there is a backbone bead.
			if(BBMASS[k] != 0){
				nBeads++;
			}

			// Is sidechain mass non-zero? If sp, there is a sidechain bead.
			if(SCMASS[k] != 0){
				nBeads++;
			}
		}
		
		// Output the header 
		// nAtoms = number of beads
		// nSites = 0?
		// nMols = number of residues
		// nMolTypes = number of residues
		fprintf(topol, "nAtoms:\t\t%d\n", nBeads);
		fprintf(topol, "nSites:\t\t%d\n", 0);
		fprintf(topol, "nMols:\t\t%d\n", nRes);
		fprintf(topol, "nMolTypes:\t%d\n", nRes);
	
		// Information for each residue
		// [molecule] = what molecule does the residue belong to? For single proteins, this is one.
		// molName = name of the residue
		// chain = chain number (NOTE: HARD CODED! CHANGE FOR COMPLEXES)
		// [atoms] = containing the bead name and the masses of each bead
		for(int k = 0; k < nRes; k++){
			fprintf(topol, "[molecule]\t%d\n", 1);
			fprintf(topol, "molName\t\t%s\n", NAMES[k]);
			fprintf(topol, "chain\t\t%s\n", "A");
			if(SCMASS[k] != 0){
				fprintf(topol, "[atoms]\t\t%d\n", 2);
				fprintf(topol, "%s\t%.4f\t%.4f\t0.0 %10.5e %10.5e\n", "BACK", BBMASS[k], 0.0, 0.0, 0.0);
				fprintf(topol, "%s\t%.4f\t%.4f\t0.0 %10.5e %10.5e\n", "SIDE", SCMASS[k], 0.0, 0.0, 0.0);
			}
			else{
				fprintf(topol, "[atoms]\t\t%d\n", 1);
				fprintf(topol, "%s\t%.4f\t%.4f\t0.0 %10.5e %10.5e\n", "BACK", BBMASS[k], 0.0, 0.0, 0.0);
			}

			// Not used for FRESEAN.
			fprintf(topol, "[sites]\t\t0\n");
			fprintf(topol, "[bonds]\t\t0\n");
			fprintf(topol, "[angles]\t\t0\n");
			fprintf(topol, "[dihedrals]\t\t0\n");
			fprintf(topol, "[impropers]\t\t0\n");
		}
		printf("Number of Beads: %d\n", nBeads);
		// Not used for FRESEAN.
		fprintf(topol, "[inter_bonds]\t\t0\n");
		fprintf(topol, "[inter_angles]\t\t0\n");
		fprintf(topol, "[inter_dihedrals]\t\t0\n");
		fprintf(topol, "[inter_impropers]\t\t0\n");

		printf("Number of Residues: %d\n", nRes);
		return nBeads;
}

int main(int argc, char *argv[])
{
	FILE *tin, *cin, *vin, *job, *out, *topol;
	XDR xdrin;
	char fnCOM[100],fnTop[100],fnJob[100],fnCrd[100],fnVel[100],fnOutTraj[100],fnOutTopol[100];
	t_atom *atoms;
	t_vecR *siteCrds;
	t_mol *mols;
	t_grp *grps;
	int analyzeGrp,nAtoms,nMols,nSites,format;
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
	int k,l,m,n,s;
	int nRead=5000;
	int nTimesAnalysis=1000;
	int analysisInterval=10;
	int nSample=10;
	
	real *timebuffer;
    t_box *boxbuffer;
    t_atom **atomsbuffer;
    t_vecR **siteCrdsbuffer;
    int bufferpos;


	printTitle();

    timebuffer=(real*)save_malloc(100*sizeof(real));
    boxbuffer=(t_box*)save_malloc(100*sizeof(t_box));
    atomsbuffer=(t_atom**)save_malloc(100*sizeof(t_atom*));
    siteCrdsbuffer=(t_vecR**)save_malloc(100*sizeof(t_vecR*));

	if(argc!=2) {
        printKeys();
        fatal("no input file specified\n");
    }

    strcpy(fnCOM,argv[1]);
	getInput(fnCOM,fnTop,fnCrd,fnVel,fnJob,&analyzeGrp,
                &nRead, &nSample, 
               	fnOutTraj, fnOutTopol, needVelo);
	
	/*READ TOPOLOGY HERE*/
	readTOP(fnTop,&atoms,&nAtoms,&siteCrds,&nSites,&mols,&nMols);
	for(i=0;i<100;i++) {
        atomsbuffer[i]=(t_atom*)save_malloc(nAtoms*sizeof(t_atom));
        siteCrdsbuffer[i]=(t_vecR*)save_malloc(nSites*sizeof(t_vecR));
    }

	/*READ JOB DATA HERE*/
	readJob(fnJob,&grpDef);

	/*SELECT STATIC GROUPS FROM TOPOLOGY*/
	allocGrps(grpDef,nAtoms,nMols,&grps);
	makeStaticGrps(atoms,nAtoms,mols,nMols,grpDef.statGrpDef,grpDef.nStat,grps);

	/*PREP CRD AND/OR VEL READ HERE*/
	prepTRJinput(fnCrd,fnVel,&cin,&vin,&xdrin,&format,&veloFreq,&time,&dt,&nFrames,needVelo,nAtoms,argc);

    out = fopen(fnOutTraj, "w"); // File we will output the coarse grain trajectory to 
    topol = fopen(fnOutTopol, "w"); // File we will output the custom .mtop file to

	/*
	Define the arrays to store the center of mass and total mass for each residue
	*/
    real COM_POS_BACKBONE[3] = {0.0, 0.0, 0.0};
    real COM_POS_SIDECHAIN[3] = {0.0, 0.0, 0.0};
    real COM_VEL_BACKBONE[3] = {0.0, 0.0, 0.0};
    real COM_VEL_SIDECHAIN[3] = {0.0, 0.0, 0.0};
    real BACKBONE_MASS = 0.0;
    real SIDECHAIN_MASS = 0.0;

	// Store how many residues there are in the structure and how many beads we have per frame
	int nRes = 0;
	int nBeads = 0;
	

	// Iterate over the topology
    for(i=0;i<nRead;i++){
		int beadCounter = 1;
		if(i%100==0) {
       		for(j=0;j<100;j++) {
            	readTRJ(i+j,format,cin,vin,&xdrin,dt,&timebuffer[j],&boxbuffer[j],nAtoms,nSites,needVelo,veloFreq,wc,atomsbuffer[j],siteCrdsbuffer[j]);
        	}
        	bufferpos=0;
    	}
		
    	readTRJbuffer(&bufferpos,&time,timebuffer,&box,boxbuffer,nAtoms,atoms,atomsbuffer,nSites,siteCrds,siteCrdsbuffer);
	    updateGrps(atoms,nAtoms,mols,nMols,box,grpDef,grps);
		
		/* Here we output the first frame of our trajectory as a reference file  
		   Define the custom topology file (requires COM calculation) */
		if(i == 0){
			nRes = printReferenceFile(atoms, nAtoms, grps[analyzeGrp]);
			printf("NRes value: %d\n", nRes);
			printf("Custom Topology Input: %d\n%d\n%d\n", nRes, nBeads, grps[analyzeGrp].nAtoms);
			nBeads = printCustomTopology(topol, atoms, nRes, grps[analyzeGrp]);
			printf("Value of nBeads: %d\n", nBeads);
			fclose(topol);
		}

		// Keep track of the current residue number so we know when to write out a new bead
		int curr_index;

		// Header for each frame
		fprintf(out, "coarse grained file\n");
		fprintf(out, "%d\n", nBeads);
		
		for(j=0;j<grps[analyzeGrp].nAtoms;j++) {

		curr_index = atoms[grps[analyzeGrp].atoms[j]].resNr;
		
		/*
		Check the residue name to see if this atom is part of a residue with no side chain
		This includes any capping residues or glycine

		If it is one of these, we only need to assign the backbone center of mass position and velocity
		*/
		if(isBackboneOnly(atoms[grps[analyzeGrp].atoms[j]].resName)){
			assignCOM(COM_POS_BACKBONE, COM_VEL_BACKBONE, &BACKBONE_MASS, atoms[grps[analyzeGrp].atoms[j]]);
		}

		/*
		If the current atom contains both a backbone and a sidechain, then we need to identify whether
		the current atom belongs to the sidechain OR to the backbone and update the center of mass
		accordingly
		*/
		else{
			if(inBackbone(atoms[grps[analyzeGrp].atoms[j]])){
				assignCOM(COM_POS_BACKBONE, COM_VEL_BACKBONE, &BACKBONE_MASS, atoms[grps[analyzeGrp].atoms[j]]);
			}
			else{
				assignCOM(COM_POS_SIDECHAIN, COM_VEL_SIDECHAIN, &SIDECHAIN_MASS, atoms[grps[analyzeGrp].atoms[j]]);
			}
		}

		/*
		If the next atom is part of a different residue, we cut off our center of mass calculation
		and output the values to the file in .gro format

		Afterwards, we reset our center of mass position, velocity, and mass
		*/
		
		if(j == grps[analyzeGrp].nAtoms - 1 || atoms[grps[analyzeGrp].atoms[j+1]].resNr != curr_index){
			if(isBackboneOnly(atoms[grps[analyzeGrp].atoms[j]].resName) == 0){
				
				divideCOM(COM_POS_BACKBONE, COM_VEL_BACKBONE, BACKBONE_MASS);
				divideCOM(COM_POS_SIDECHAIN, COM_VEL_SIDECHAIN, SIDECHAIN_MASS);
				printBackbone(out, atoms[grps[analyzeGrp].atoms[j]], COM_POS_BACKBONE, COM_VEL_BACKBONE, &beadCounter);
				printSidechain(out, atoms[grps[analyzeGrp].atoms[j]], COM_POS_SIDECHAIN, COM_VEL_SIDECHAIN, &beadCounter);
				resetCOM(COM_POS_BACKBONE, COM_VEL_BACKBONE, &BACKBONE_MASS);
				resetCOM(COM_POS_SIDECHAIN, COM_VEL_SIDECHAIN, &SIDECHAIN_MASS);
			}
					
			else{
		
				divideCOM(COM_POS_BACKBONE, COM_VEL_BACKBONE, BACKBONE_MASS);
				printBackbone(out, atoms[grps[analyzeGrp].atoms[j]], COM_POS_BACKBONE, COM_VEL_BACKBONE, &beadCounter);
				resetCOM(COM_POS_BACKBONE, COM_VEL_BACKBONE, &BACKBONE_MASS);
			}
			
		}
	}
	
	fprintf(out,"%.5f %.5f %.5f\n", boxbuffer[0].a.x/10, boxbuffer[0].b.y/10, boxbuffer[0].c.z/10);
	printf("step %d of %d%c",i+1,nRead,(char)13); fflush(stdout);
   	}
	printf("\n");
	closeInput(format,needVelo,cin,vin,&xdrin);
	fclose(out);
	
	return 0;

}
