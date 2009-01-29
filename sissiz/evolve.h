#ifndef _EVOLVE_H_
#define _EVOLVE_H_

#include "tree.h"
#include "mutate.h"


#define HAMMINGDISTANCE 1

extern int numBases, numTaxa;
#if HAMMINGDISTANCE
extern int writenote, writeHammingdistance, writeHammingdistanceTree;
extern int neighborFile, verboseancestor;
#endif

/* prototypes */
void EvolveSequencesZ(TTree *tree, double scale, char *ancestor,double **QT2ij);
void EvolveNodeZ(TNode *anc, TNode *des, double scale,double **QT2ij);

void EvolveSequencesInd(TTree *tree, double scale, char *ancestor,double *Qij);
void EvolveNodeInd(TNode *anc, TNode *des, double scale,double *Qij);


void WriteSequences(FILE *fv, TTree *tree);                                         /*seq-gen*/
void WriteAncestralSequences(FILE *fv, TTree *tree);                                /*seq-gen*/
void WriteAncestralSequencesNode(FILE *fv, TTree *tree, int *nodeNo, TNode *des);   /*seq-gen*/

#if HAMMINGDISTANCE
void WriteHammingdistanceSeq(FILE *fv, TTree *tree);
void WriteHammingdistanceNode(FILE *fv, TTree *tree, int *nodeNo,TNode *des, TNode *ancestral, double length,int *count0branch, int parentNo, int direction);
#endif

#endif /* _EVOLVE_H_ */
