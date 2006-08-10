#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dataDef.h"
#include "ToMatrix.h"



/********************************************************************************************************
* FUNCTION: "ensembleToMatrix"
*
* ARGUMENTS:  DandO
*             nodePositions     for each node gives: node number, parent node, child nodes, tree number
*             nodeParameters    for each node gives: prediction, changeRSS, upper and lower limits
*             treeStart         treeStart[i] gives the start node of i-th tree 
*             treeEnd           treeEnd[i]   gives the last node (in nodePositions and nodeParameters) of tree "i"
*
* RETURNS: 
*
* PURPOSE:  This writes an ensemble of trees to a couple of vectors   
*  
* NOTE:    In use with R, the vectors (not DandO) must be allocated from R 
*          (max length from "depth" and "number trees")
********************************************************************************************************/

void ensembleToMatrix(struct DatAndOpt *DandO,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes)
{
int i,num,startIndx,endIndx;
struct treeList *treeListPtr;

/*********************************************
* Compute total number of nodes (in ensemble)
*********************************************/
treeListPtr=DandO->firstTree;
num=0;
for(i=0;i<DandO->numberOfTrees;i++)
{

num+=treeListPtr->numberOfNodes;
treeListPtr=treeListPtr->nextTree;
}
 
*totNumNodes=num;

treeListPtr=DandO->firstTree;

startIndx=0; /* start position of i-th tree*/

for(i=0;i<DandO->numberOfTrees;i++)
{

endIndx=startIndx+treeListPtr->numberOfNodes-1;  /* end position of i-th tree*/

treeStart[i]=startIndx;
treeEnd[i]=endIndx;

treeToVectors(treeListPtr->topNode,startIndx,endIndx,nodePositions,nodeParameters,DandO,*totNumNodes,treeListPtr->numberOfNodes);

treeListPtr=treeListPtr->nextTree;
startIndx=endIndx+1;
} /* for(i=0;i<DandO->numberOfTrees;i++)*/



}
















/********************************************************************************************************
* FUNCTION: "treeToVectors"
*
* ARGUMENTS:   rootNode            pointer to root node 
*              treeMatrix          matrix to which various node info is written, each row corresponds a node
*              lowerLimit          matrix of lower limits of p-dimensional rectangles defining nodes
*              upperLimit          matrix of upper limits of p-dimensional rectangles defining nodes
*             
* RETURNS:  void- info written to matrices  
*                 
* PURPOSE: 
**********************************************************************************************************/
void treeToVectors(struct node* rootNode,int startIndx,int endIndx,int *nodePositions,double *nodeParameters,struct DatAndOpt *DandO,int totalNumberNodes,int NumNodesTree)
{
int AddToMatrix,rowIndx,lastNodeLead;
struct unSplitList *currNode;
struct unSplitList *prevNode;
struct unSplitList *lastNode;

  
currNode=(struct unSplitList*) malloc(sizeof(struct unSplitList));

currNode->pnode=rootNode;
lastNode=currNode;
currNode->nextUnsplit=NULL;

lastNodeLead=0;
AddToMatrix=1;rowIndx=0;
while(AddToMatrix)
{
 
/*****************
* Node number
******************/
nodePositions[rowIndx+startIndx]=currNode->pnode->nodeNum;
/********************
* Parent Node number
*********************/
nodePositions[rowIndx+startIndx+totalNumberNodes]=currNode->pnode->parentNodeNum;
/*************************
* Left child Node position
***************************/
nodePositions[rowIndx+startIndx+2*totalNumberNodes]=-99;



/*************************
* Right child Node position
***************************/
nodePositions[rowIndx+startIndx+3*totalNumberNodes]=-99;

/**********************************************
* Index of variable node was split on (if any) 
***********************************************/
nodePositions[rowIndx+startIndx+4*totalNumberNodes]=currNode->pnode->splitVar;


/*************************
* split point of node
***************************/
nodeParameters[rowIndx+startIndx]=currNode->pnode->splitPoint;

/*************************
* prediction at node
***************************/
nodeParameters[rowIndx+startIndx+1*totalNumberNodes]=currNode->pnode->theta;


/*****************************
* change in RSS due to split
******************************/
nodeParameters[rowIndx+startIndx+2*totalNumberNodes]=currNode->pnode->changeRSS;

/****************************************
* Estimated number of data points in node
****************************************/
nodeParameters[rowIndx+startIndx+3*totalNumberNodes]=currNode->pnode->numInNode;





if(currNode->pnode->splitVar>-10)  /* if not split, then splitVar==-99 */
{


/* add left child to list */
lastNode->nextUnsplit=(struct unSplitList*) malloc(sizeof(struct unSplitList));

lastNode=lastNode->nextUnsplit;
lastNodeLead++; 
nodePositions[rowIndx+startIndx+2*totalNumberNodes]=lastNodeLead+startIndx;

/***********************************************************************************
*  "lastNodeLead" counts how many nodes from current node being written to matrix to 
*   to last node on "unSplitList"
***********************************************************************************/
lastNode->pnode=currNode->pnode->leftChild;
lastNode->nextUnsplit=NULL;

/* add right child to list */
lastNode->nextUnsplit=(struct unSplitList*) malloc(sizeof(struct unSplitList));
lastNode=lastNode->nextUnsplit;
lastNodeLead++; 
nodePositions[rowIndx+startIndx+3*totalNumberNodes]=lastNodeLead+startIndx;

lastNode->pnode=currNode->pnode->rightChild;
lastNode->nextUnsplit=NULL;
}
rowIndx++;


prevNode=currNode;
currNode=currNode->nextUnsplit;




if(NumNodesTree<rowIndx)
AddToMatrix=0;


if(currNode==NULL)
AddToMatrix=0;


free(prevNode); 
} /* matches "while(AddToMatrix)"*/



}








