/********************************************************************************************************
*  
*  File contains functions for free memory allocated to trees, forests and nodes
*  
*
*
*  
*                  
**********************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dataDef.h"
#include "freeTree.h"
#include "gbev.h"
#include "independentPred.h"


/********************************************************************************************************
* FUNCTION: "freeDatAndOpt"
*
* ARGUMENTS: pointer to DatAndOpt structure.
*
*
* RETURNS:  void- function frees all allocated memory associated with a DatAndOpt structure
*                  
**********************************************************************************************************/
void freeDatAndOpt(struct DatAndOpt *DandO)
{
int i,go;

free(DandO->oobIndicator);

/* Constraint points*/
if(DandO->useConstraints==1)
{
free(DandO->fitConstraintPoints);
free(DandO->constraintPoints);
free(DandO->nodeConstPoints);
}


free(DandO->w);
free(DandO->yRes);
free(DandO->minW);
free(DandO->maxW);
free(DandO->wcWeights);

/* Frees various quantities in DandO*/
divFreeDandO(DandO);


/*****************************************
* Free space for measurement error model
******************************************/
for(i=0;i<(DandO->n);i++)
free(DandO->ratio_pw[i]);

free(DandO->ratio_pw);


for(i=0;i<(DandO->NumberComponents*DandO->n);i++)
free(DandO->condMeanComp[i]);

free(DandO->condMeanComp);


for(i=0;i<DandO->NumberComponents;i++)
free(DandO->cholCov[i]);

free(DandO->cholCov);

for(i=0;i<(DandO->n*DandO->mc);i++)
free(DandO->x_mc[i]);

free(DandO->x_mc);

free(DandO->node_mc);
free(DandO->yPred);
if(DandO->SplitFunction==3)
{
free(DandO->pred_mc); 
free(DandO->predLink_mc);
free(DandO->yPredLink);
}




go=1;
if(go==1)
{

freeForest(DandO->firstTree,DandO->numberOfTrees);
}

}


/********************************************************************************************************
* FUNCTION: "freeForest"
*
* ARGUMENTS: firstTreePtr  pointer to first tree in forest to be deleted
*
*
* RETURNS:  void- function frees the forest whose first tree is pointed by firstTreePtr
*                  
**********************************************************************************************************/
void freeForest(struct treeList *firstTreePtr,int numberTrees)
{
int i;
struct treeList *treeListPtr;
struct treeList *treeListPtrOld;

treeListPtr=firstTreePtr;

for(i=0;i<numberTrees;i++)
{

freeTree(treeListPtr->topNode);
treeListPtrOld=treeListPtr;

treeListPtr=treeListPtr->nextTree;

free(treeListPtrOld);

}

}




/********************************************************************************************************
* FUNCTION: "freeTree"
*
* ARGUMENTS: ptrTopNode  pointer to the top node of tree to be deleted
*
*
* RETURNS:  void- function frees the tree whose top node is pointed to by "ptrTopNode"
*                  
**********************************************************************************************************/
void freeTree(struct node *ptrTopNode)
{
int go;
struct unSplitList *currUnsplit;
struct unSplitList *prevUnsplit;
struct unSplitList *lastUnsplit;


currUnsplit=(struct unSplitList*) malloc(sizeof(struct unSplitList));

currUnsplit->pnode=ptrTopNode;
currUnsplit->nextUnsplit=NULL;
lastUnsplit=currUnsplit;

go=1;
while(go)
{

if(currUnsplit->pnode->splitVar >-10)  
{
/* If node hasnt been split (ie terminal node), then "splitVar=-99"*/
lastUnsplit->nextUnsplit=(struct unSplitList*) malloc(sizeof(struct unSplitList));
lastUnsplit->nextUnsplit->pnode=currUnsplit->pnode->leftChild;
lastUnsplit=lastUnsplit->nextUnsplit;


lastUnsplit->nextUnsplit=(struct unSplitList*) malloc(sizeof(struct unSplitList));
lastUnsplit->nextUnsplit->pnode=currUnsplit->pnode->rightChild;
lastUnsplit=lastUnsplit->nextUnsplit;
lastUnsplit->nextUnsplit=NULL;

}

freeNode(currUnsplit->pnode);
prevUnsplit=currUnsplit;
currUnsplit=currUnsplit->nextUnsplit;

free(prevUnsplit);

if(currUnsplit==NULL)
go=0;

} /* matches while(go) */


}





/********************************************************************************************************
* FUNCTION: "freeNode"
*
* ARGUMENTS: ptrNode  pointer to a node
*
*
* RETURNS:  void- function frees the node pointed to by "ptrNode"
*                  
**********************************************************************************************************/
void freeNode(struct node *ptrNode)
{

/* first free all elements pointed to within node.*/


free(ptrNode->probInNode);


free(ptrNode->upperLim);


free(ptrNode->lowerLim);

free(ptrNode);
ptrNode=NULL;

}




/********************************************************************************************************
* FUNCTION: "freeDatAndOpt"
*
* ARGUMENTS: pointer to DatAndOpt structure.
*
*
* RETURNS:  void- function frees all allocated memory associated with a DatAndOpt structure
*                  
**********************************************************************************************************/
void freeDatAndOpt2(struct DatAndOpt *DandO)
{
int i;

free(DandO->oobIndicator);

/* Constraint points*/
if(DandO->useConstraints==1)
{
free(DandO->fitConstraintPoints);
free(DandO->constraintPoints);
free(DandO->nodeConstPoints);
}


free(DandO->w);
free(DandO->yRes);
free(DandO->minW);
free(DandO->maxW);
free(DandO->wcWeights);


/* Frees various quantities in DandO */
divFreeDandO(DandO);
divFreePred(DandO);


/*****************************************
* Free space for measurement error model
******************************************/
for(i=0;i<(DandO->n);i++)
free(DandO->ratio_pw[i]);

free(DandO->ratio_pw);


for(i=0;i<(DandO->NumberComponents*DandO->n);i++)
free(DandO->condMeanComp[i]);

free(DandO->condMeanComp);


for(i=0;i<DandO->NumberComponents;i++)
free(DandO->cholCov[i]);

free(DandO->cholCov);

for(i=0;i<(DandO->n*DandO->mc);i++)
free(DandO->x_mc[i]);

free(DandO->x_mc);

free(DandO->node_mc);
free(DandO->yPred);

if(DandO->SplitFunction==3)
{
free(DandO->pred_mc); 
free(DandO->predLink_mc);
free(DandO->yPredLink);
}



}


