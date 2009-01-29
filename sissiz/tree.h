#ifndef _TREE_H_
#define _TREE_H_

#define MAX_NAME_LEN 256

typedef struct TNode TNode;
struct TNode {
	TNode *branch0, *branch1, *branch2, *next;
	double length0, length1, length2, param;
	int tipNo;
	
	char *sequence;
};


typedef struct TTree TTree;
struct TTree {
	int rooted, lengths;
	TNode *root, *nodeList;
	int numTips, numNodes;
	double totalLength;
	char **names;
	TNode **tips;
	int capacity;
};

#endif /* _TREE_H_ */

