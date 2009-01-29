#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "global.h"
#include "tree.h"
#include "evolve.h"
#include "models.h"
#include "mutate.h"


#define HAMMINGDISTANCE 1

int numBases, numTaxa;
#if HAMMINGDISTANCE
int writenote, writeHammingdistance, writeHammingdistanceTree;
#endif
#if HAMMINGDISTANCE
FILE *hammingtree_fv;
FILE *interntree_fv;
FILE *internlength_fv;
FILE *anc_fv;
#endif
int neighborFile;
int verboseancestor;
char ancFilename[256];


/* prototypes */
void EvolveSequencesZ(TTree *tree, double scale, char *ancestor,double **QT2ij);
void EvolveNodeZ(TNode *anc, TNode *des, double scale,double **QT2ij);

void EvolveSequencesInd(TTree *tree, double scale, char *ancestor,double *Qij);
void EvolveNodeInd(TNode *anc, TNode *des, double scale,double *Qij);

void WriteSequences(FILE *fv, TTree *tree);
void WriteAncestralSequences(FILE *fv, TTree *tree);
void WriteAncestralSequencesNode(FILE *fv, TTree *tree, int *nodeNo, TNode *des);
#if HAMMINGDISTANCE
void WriteHammingdistanceSeq(FILE *fv, TTree *tree);
void WriteHammingdistanceNode(FILE *fv, TTree *tree, int *nodeNo, TNode *des, TNode *ancestral, double length, int *count0branch, int parentNo, int direction);
#endif

/* functions */

void EvolveSequencesZ(TTree *tree, double scale, char *ancestor,double **QT2ij)
{
	if (ancestor==NULL){
		/*if(random) {RandomSequence(tree->root->sequence);}
		else {*/
		/*RandomSequence_independent(tree->root->sequence);*/
		RandomSequenceMarkov1(tree->root->sequence);
		/*}*/
	}	
	else SetSequence(tree->root->sequence, ancestor);
	
	EvolveNodeZ(tree->root, tree->root->branch1, scale, QT2ij);
	EvolveNodeZ(tree->root, tree->root->branch2, scale, QT2ij);
	if (!tree->rooted)
		EvolveNodeZ(tree->root, tree->root->branch0, scale,QT2ij);
}



void EvolveSequencesInd(TTree *tree, double scale, char *ancestor,double *Qij)
{
	if (ancestor==NULL){
	     /* RandomSequence(tree->root->sequence);*/
    RandomSequence_independent(tree->root->sequence);
	}else SetSequence(tree->root->sequence, ancestor);
	   
	EvolveNodeInd(tree->root, tree->root->branch1, scale, Qij);
	EvolveNodeInd(tree->root, tree->root->branch2, scale, Qij);
	if (!tree->rooted)
		EvolveNodeInd(tree->root, tree->root->branch0, scale,Qij);
}

void EvolveNodeZ(TNode *anc, TNode *des, double scale,double **QT2ij)
{
	double len;
	
	
	if (scale<0.0)
		len=des->length0;
	else
		len=des->length0*scale;
		
	memcpy(des->sequence, anc->sequence, seqlen);
	
	if(fastqx==1)
	  MutateSequenceZqx(des->sequence, len,QT2ij);
	else
	  MutateSequenceZ(des->sequence, len,QT2ij);
	
	if (des->tipNo==-1) {
		EvolveNodeZ(des, des->branch1, scale, QT2ij);
		EvolveNodeZ(des, des->branch2, scale, QT2ij);
	}
	
	
}

void EvolveNodeInd(TNode *anc, TNode *des, double scale,double *Qij)
{
	double len;
	if (scale<0.0)
		len=des->length0;
	else
		len=des->length0*scale;
		
	memcpy(des->sequence, anc->sequence, seqlen);
	
	if(fastqx==1)
     MutateSequenceIndqx(des->sequence, len,Qij);
	else 
	 MutateSequenceInd(des->sequence, len,Qij);
	
	if (des->tipNo==-1) {
		EvolveNodeInd(des, des->branch1, scale, Qij);
		EvolveNodeInd(des, des->branch2, scale, Qij);
	}
}




/************************************** WRITE ********************************/
void WriteSequences(FILE *fv, TTree *tree)
{
	int i, j;
	char *P;


	if(FastaAln==0 && Clustal==0) fprintf(fv, " %d %d\n", numTaxa, numBases);
	if(Clustal==1) fprintf(fv, "CLUSTAL \n");
    relaxedPhylip=1;
	for (i=0; i<tree->numTips; i++) {
	    if(FastaAln==1){
		fprintf(fv, ">%s\n", tree->names[i]);
		}
		else if (relaxedPhylip)
			fprintf(fv, "%s ", tree->names[i]);
			/* ........ */
		else {
			j=0;
			P=tree->names[i];
			while (j<10 && *P) {
				fputc(*P, fv);
				j++;
				P++;
			}
			while (j<10) {
				fputc(' ', fv);
				j++;
			}
		}
		if(printfile==1 && verboseancestor==1){
	             sprintf(ancFilename, "%s_anc.sissi", filename);
	             if((anc_fv=fopen(ancFilename, "w"))==NULL)
	             fprintf(stderr, "anc.txt could not open\n");
	        }	
        P=tree->tips[i]->sequence;
		for (j=0; j<seqlen; j++) {
			fputc(*P, fv);
			if(printfile==1 && verboseancestor==1) fputc(*P, anc_fv);
			P++;
		}
		fputc('\n', fv);
		
	}
	if(FastaAln==1) fputc('\n', fv);
}


#if HAMMINGDISTANCE
void WriteHammingdistanceSeq(FILE *fv, TTree *tree){
	int  n;
	int count0branch=0;
	int myNo=0;

	
	
	
	if(writeHammingdistanceTree==1){
	   if((hammingtree_fv=fopen("hammingdistance.htree","a"))==NULL)
	      fprintf(fv, "hammingdistance.tree could not open\n");
	   if((interntree_fv=fopen("intern.htree","a"))==NULL)
	      fprintf(fv, "intern.htree could not open\n");	
	   if((internlength_fv=fopen("internlength.htree","a"))==NULL)
	      fprintf(fv, "internlength.htree could not open\n");	 	    
	}	    
	  	
	n = numTaxa + 1;
	myNo=n;
	
	if(writeHammingdistance){
	  fprintf(fv, "\n\nHAMMINGDTSTANCE");
	  fputc('\n', fv);    
	  fprintf(fv, "           ->    node   | length | hammingdistance \n");
	}  
	
	if(writeHammingdistanceTree==1){
    	fprintf(hammingtree_fv,  "(");
	    fprintf(interntree_fv,   "(");
		fprintf(internlength_fv, "(");
	}
	
    if(!tree->rooted){
	   WriteHammingdistanceNode(fv,tree, &n, tree->root->branch0, tree->root,tree->root->length0,&count0branch, myNo, 0);
	   
       if(writeHammingdistanceTree==1){
     	   fprintf(hammingtree_fv,  ",");
	       fprintf(interntree_fv,   ",");
		   fprintf(internlength_fv, ",");
	   }	   
	}   
	
	WriteHammingdistanceNode(fv,tree, &n, tree->root->branch1,tree->root,tree->root->length1,&count0branch, myNo, 1);
	
	if(writeHammingdistanceTree==1){
     	fprintf(hammingtree_fv,  ",");
	    fprintf(interntree_fv,   ",");
		fprintf(internlength_fv, ",");
    }	
	
	WriteHammingdistanceNode(fv,tree, &n, tree->root->branch2,tree->root,tree->root->length2,&count0branch, myNo, 2);
	
	if(writeHammingdistanceTree==1){
        fprintf(hammingtree_fv,  ")");
	    fprintf(interntree_fv,   ")");
		fprintf(internlength_fv, ")");
	    if(tree->rooted){
	       fprintf(hammingtree_fv,   "%d;\n",myNo);
	       fprintf(interntree_fv,    "%d;\n",myNo);
		   fprintf(internlength_fv, "%d;\n",myNo);
		   
	    }  else  {
		   fprintf(hammingtree_fv,  ";\n");
           fprintf(interntree_fv,   ";\n");
		   fprintf(internlength_fv, ";\n");
	}	   
	}  
	
   	fputc('\n', fv);
	if(count0branch != 0 || writenote==1){
	 fprintf(fv, "NOTE: %d branches without observed differences (Hammingdistance=0)\n", count0branch);
	 fputc('\n', fv); 
	} 
	
	if(writeHammingdistanceTree==1){
    	fclose(hammingtree_fv);
	    fclose(interntree_fv);
		fclose(internlength_fv);
    }
}

void WriteHammingdistanceNode(FILE *fv, TTree *tree, int *nodeNo, TNode *des, TNode *ancestral, double length,int *count0branch, int parentNo, int direction){
     int  j;
	 int  countbranchhammingdis;
	 char *P= des->sequence;
	 char *R= ancestral->sequence;
	 int  myNo;
	 
	 
	 
	 
	
	 if(des->tipNo==-1){
		(*nodeNo)++;  	
	    myNo= *nodeNo;
        	   
	    countbranchhammingdis=0;
		
				 
		 
				
	    for (j=0; j<seqlen; j++) {
	    if(*P!= *R) countbranchhammingdis++; 
		  P++;
		  R++;
	    }
		if(countbranchhammingdis==0)(*count0branch)++; 
		
		if(writeHammingdistance)
	    	fprintf(fv, "%10d ->%10d | %3.4f | %10d\n",parentNo, *nodeNo,length, countbranchhammingdis);
		
		if(writeHammingdistanceTree==1){
     		fprintf(hammingtree_fv,  "(");
	    	fprintf(interntree_fv,   "(");
			fprintf(internlength_fv, "(");
		}
			
	
		WriteHammingdistanceNode(fv,tree, nodeNo, des->branch1,des,des->length1,count0branch, myNo,1);
		
		if(writeHammingdistanceTree==1){
     		fprintf(hammingtree_fv,  ",");
	    	fprintf(interntree_fv,   ",");
			fprintf(internlength_fv, ",");
		}
			
		WriteHammingdistanceNode(fv,tree, nodeNo, des->branch2,des,des->length2,count0branch, myNo,2);
		
	    if(writeHammingdistanceTree==1){
    		fprintf(hammingtree_fv,  ")");
		    fprintf(interntree_fv,   ")");
			fprintf(internlength_fv, ")");
		    fprintf(hammingtree_fv,    "%d:%d", myNo,countbranchhammingdis);
		    fprintf(interntree_fv,        "%d", myNo);
			fprintf(internlength_fv, "%d:%.2f", myNo,length);
        }
				
	} else{
		countbranchhammingdis=0;
	    for (j=0; j<seqlen; j++) {
		  if(*P!=*R) countbranchhammingdis++; 
		  P++;
		  R++;
	    }
		
		if(countbranchhammingdis==0) (*count0branch)++;
		
		if(writeHammingdistance) 
	    	fprintf(fv, "%10d ->%10s | %3.4f | %10d\n", parentNo, tree->names[des->tipNo],length, countbranchhammingdis);
		
		if(writeHammingdistanceTree==1){
    		if(direction==1) fprintf(hammingtree_fv,  "%s:%d",tree->names[des->tipNo],countbranchhammingdis);
	    	if(direction==2) fprintf(hammingtree_fv,  "%s:%d",tree->names[des->tipNo],countbranchhammingdis);
	    	if(direction==1) fprintf(interntree_fv,   "%s",tree->names[des->tipNo]);
		    if(direction==2) fprintf(interntree_fv,   "%s",tree->names[des->tipNo]);
			if(direction==1) fprintf(internlength_fv, "%s:%.2f",tree->names[des->tipNo],length);
		    if(direction==2) fprintf(internlength_fv, "%s:%.2f",tree->names[des->tipNo],length);
		}
    }
}
#endif

 /*seq-gen*/
void WriteAncestralSequences(FILE *fv, TTree *tree)
{
	int j, n;
	char *P;

	if (!tree->rooted)
		n = (2 * numTaxa) - 3;
	else
		n = (2 * numTaxa) - 2;
	fprintf(fv, " %d %d\n", n, numBases);

	n = numTaxa + 1;
	
	fprintf(fv, "%d\t", n);
	P=tree->root->sequence;
	for (j=0; j<numBases; j++) { 
		fputc(*P, fv);
		P++;
	}
	fputc('\n', fv);
	
	if (!tree->rooted)
		WriteAncestralSequencesNode(fv, tree, &n, tree->root->branch0);
	WriteAncestralSequencesNode(fv, tree, &n, tree->root->branch1);
	WriteAncestralSequencesNode(fv, tree, &n, tree->root->branch2);
}
 /*seq-gen*/
void WriteAncestralSequencesNode(FILE *fv, TTree *tree, int *nodeNo, TNode *des)
{
	int j;
	char *P = des->sequence;
	
	if (des->tipNo==-1) {
		(*nodeNo)++;
		
		fprintf(fv, "%d\t", *nodeNo);
		for (j=0; j<numBases; j++) {
			fputc(*P, fv);
			P++;
		}
		fputc('\n', fv);
		
		WriteAncestralSequencesNode(fv, tree, nodeNo, des->branch1);
		WriteAncestralSequencesNode(fv, tree, nodeNo, des->branch2);
	} else {
		fprintf(fv, "%s\t", tree->names[des->tipNo]);
		for (j=0; j<numBases; j++) {
			fputc(*P, fv);
			P++;
		}
		fputc('\n', fv);
	}

}





