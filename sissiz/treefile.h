#ifndef _TREEFILE_H_
#define _TREEFILE_H_

/* This contains the structures, TTree and TNode which */
/* can be edited to add more elements                  */

#include "tree.h"

TTree *NewTree();
void DisposeTree(TTree *tree);
void FreeTree(TTree *tree);
void WriteAvailInfo();

int CountTrees(FILE *fv);
void ReadTree(FILE *fv, TTree *tree, int treeNum, int numNames, char **names, 
				int *outNumSites, double *outRelRate);
int IsTreeAvail(FILE *fv);
void WriteTree(FILE *fv, TTree *tree);

void UnrootRTree(TTree *tree);
void RerootUTree(TTree *tree, int tip);

/*Moved prototypes from treefile.c here */

TNode *NewNode(TTree *tree);
void InitTree(TTree *tree);
void CheckCapacity(TTree *tree, int required);
void DisposeNode(TNode *aNode);
void DisposeTreeNodes(TTree *tree);
void FreeNodes(void);

char ReadToNextChar(FILE *fv);
void ReadUntil(FILE *fv, char stopChar, char *what);
TNode *ReadTip(FILE *fv, char ch, TTree *tree, int numNames, char **names);
TNode *ReadNode(FILE *fv, TTree *tree, int numNames, char **names, int detectPolytomies);
TNode *ReadBranch(FILE *fv, TTree *tree, int numNames, char **names);
void WriteNode(FILE *fv, TTree *tree, TNode *node);



#endif /* _TREEFILE_H_ */

