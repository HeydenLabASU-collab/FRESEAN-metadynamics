#include <stdio.h>
#include <stdlib.h>

typedef struct t_branch {
    int node;
    struct t_branch* links;
    int nLinks;
} t_branch;

typedef struct t_trees {
    t_branch *trees;
    int nTrees;
} t_trees;

int getNextRoot(int nAtoms,int *atomTags) {
    int i;
    int node=-1;
    for(i=0;i<nAtoms;i++) {
        if(atomTags[i]==0) {
            node=i;
            atomTags[i]=1;
            return node;
        }
    }
    return node;
}

int buildTree(    t_branch *tree,
        int nAtoms,int *atomTags,
        int nBonds,int *bondTags,
        int *bondList,int *bondStart) {
    int i,b0,b1;
    int node;
    int tmp[10]; /*temporary storage*/
    int count=0; /*number of bonds found for current node*/
    int continuous=1;  /*track if any bonds were skipped*/
    t_branch *next;
    int error;

    node=tree->node;
    i=bondStart[0];
    while(i<nBonds && count<10) {
        b0=bondList[i*2+0];
        b1=bondList[i*2+1];
    
        if((b0<0 || b0>=nAtoms) || (b1<0 || b1>=nAtoms)) {
            FILE *err;
            err=fopen("error.log","a");
            fprintf(err,"ERROR: function wrap()\n aroms %d %d in bond %d out of range (0-%d)\n",b0,b1,i,nAtoms);
            fclose(err);
            return 1;
        }
        if(b0==node && atomTags[b1]==0) {
            tmp[count]=b1;
            count++;
            atomTags[b1]=1;
            bondTags[i]=1;
        }
        if(b1==node && atomTags[b0]==0) {
            tmp[count]=b0;
            count++;
            atomTags[b0]=1;
            bondTags[i]=1;
        }
        if(atomTags[b0]==1 && atomTags[b1]==1) {
            bondTags[i]=1;
        }
        if(continuous==1 && bondTags[i]==0) {
            continuous=0;
            bondStart[0]=i;
        }
        i++;
        if(continuous==0 && (b0>node && b1>node)) {
            /*bondList is assumed to be sorted*/
            /*therefore we can break the loop now*/
            break;
        }
    }
    if(count==10) {
        FILE *err;
        err=fopen("error.log","a");
        fprintf(err,"ERROR: function wrap()\n found >=10 bonds for atom %d\n",node);
        fclose(err);
        return 1;
    }
    tree->nLinks=count;
    tree->links=(t_branch*)malloc(count*sizeof(t_branch));
    for(i=0;i<count;i++) {
        tree->links[i].node=tmp[i];
    }
    if(count==0) {
        return 0;
    } else {
        for(i=0;i<count;i++) {
            next=&(tree->links[i]);
            error=buildTree(next,nAtoms,atomTags,nBonds,bondTags,bondList,bondStart);
            if(error==1) return 1;
        }
        return 0;
    }
}

int buildTrees(t_trees *trees,int nAtoms,int *atomTags,int nBonds,int *bondTags,int *bondList,float *masses) {
    int cnt=0;
    int node;
    int bondStart=0;
    t_branch *tree;
    int error;
    t_branch *tmp;
    int i;

    trees->trees=(t_branch*)malloc(nAtoms*sizeof(t_branch));
    node=getNextRoot(nAtoms,atomTags);
    while(node!=-1) {
        tree=&(trees->trees[cnt]);
        tree->node=node;
        error=buildTree(tree,nAtoms,atomTags,nBonds,bondTags,bondList,&bondStart);
        if(error==1) return 1;
        if(tree->nLinks==0 && masses[node]==0.0) {
            /*found isolated atoms with zero mass == dummy atom*/
            /*ensure that this is not the first tree (not allowed = error)*/
            if(cnt==0) return 1;
            /*will pretend there is a bond with first atom of previous molecule/tree*/
            /*first we create a copy of all branches connected to root of last tree*/
            /*we reserve space for one extra branch (+1) to connect the isolated atom*/
            tmp=(t_branch*)malloc((trees->trees[cnt-1].nLinks+1)*sizeof(t_branch));
            for(i=0;i<trees->trees[cnt-1].nLinks;i++) {
                tmp[i].node=trees->trees[cnt-1].links[i].node;
                tmp[i].links=trees->trees[cnt-1].links[i].links; /*pointers connecting to downstream branches*/
                tmp[i].nLinks=trees->trees[cnt-1].links[i].nLinks;
            }
            tmp[i].node=node;
            tmp[i].nLinks=0;
            free(trees->trees[cnt-1].links);
            trees->trees[cnt-1].links=tmp;
            trees->trees[cnt-1].nLinks+=1;
        } else {
            cnt++;
        }
        node=getNextRoot(nAtoms,atomTags);
    }
    trees->nTrees=cnt;
    return 0;
}

int unwrapTree(t_branch branch,float *coords,float *box) {
    int i,m;
    int node;
    float nodeCrd[3];
    int atomIdx;
    float atomCrd[3];
    float link[3];
    int nLinks;
    int error;

    node=branch.node;
    for(m=0;m<3;m++) {
        nodeCrd[m]=coords[3*node+m];
    }
    nLinks=branch.nLinks;
    for(i=0;i<nLinks;i++) {
        atomIdx=branch.links[i].node;
        for(m=0;m<3;m++) {
            atomCrd[m]=coords[3*atomIdx+m];
            link[m]=atomCrd[m]-nodeCrd[m];
            while(link[m]>box[m]/2.0) {
                link[m]-=box[m];
            }
            while(link[m]<-1.0*box[m]/2.0) {
                link[m]+=box[m];
            }
            coords[3*atomIdx+m]=nodeCrd[m]+link[m];
        }
        if(branch.links[i].nLinks!=0) {
            error=unwrapTree(branch.links[i],coords,box);
            if(error==1) return 1;
        }
    }
    return 0;
}

int unwrap(t_trees trees,float *coords,float *box) {
    int i,j;
    int error;
    for(i=0;i<trees.nTrees;i++) {
        error=unwrapTree(trees.trees[i],coords,box);
        if(error==1) return 1;
    }
    return 0;
}
