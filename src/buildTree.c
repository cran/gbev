#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dataDef.h"
#include "splitter.h"
#include "Initialize.h"
#include "buildTree.h"
#include "nodeprob.h"
#include  "splitter_binary.h"
#include "auxFunctions.h"

#define FUDGE 0.0001
/* #define DEBUG */

/********************************************************************************************************
* FUNCTION: "buildTree"
*
* ARGUMENTS: void
*            
*             
* RETURNS: pointer to root node of tree
* 
********************************************************************************************************/





struct node* buildTree(struct DatAndOpt *DandO)
{
int i; /* used for debugging*/
int NodesToSplit;
int nodeNum,nodeCounter;
double h;




struct unSplitList *currUnsplit;
struct unSplitList *prevUnsplit;
struct unSplitList *lastUnsplit;


struct node *rootNode;
rootNode=DandO->nodeVEC[0];
currUnsplit=(struct unSplitList*) malloc(sizeof(struct unSplitList));

for(i=0;i<((DandO->n)*(DandO->mc));i++)
DandO->node_mc[i]=1;





/*******************************************
* Create and initialize root node
*******************************************/
InitiRootNode(rootNode,DandO);
 


/*******************************************
* Compute new residual vector
*******************************************/
rootNode->priorRSS=0.0;

if(DandO->SplitFunction==2){
for(i=0;i<DandO->n;i++){
h=DandO->yPred[i]*(1.0-DandO->yPred[i])+FUDGE;

DandO->yRes[i]=DandO->yRes[i]-(rootNode->theta)*(rootNode->probInNode[i])/h;
}
}

if(DandO->SplitFunction!=2){
for(i=0;i<DandO->n;i++)
DandO->yRes[i]=DandO->yRes[i]-(rootNode->theta)*(rootNode->probInNode[i]);
}

for(i=0;i<DandO->n;i++){
if(DandO->oobIndicator[i]==0){
rootNode->priorRSS=rootNode->priorRSS+(DandO->yRes[i])*(DandO->yRes[i])*DandO->weights[i];
}
}
#ifdef DEBUG
Rprintf(" RSS after rootnode=%lf \n",rootNode->priorRSS);
#endif 


/*****************************************
* Initialize unSplit list
*****************************************/
currUnsplit->pnode=rootNode;
currUnsplit->nextUnsplit=NULL;
lastUnsplit=currUnsplit;
prevUnsplit=currUnsplit;
nodeNum=1;
currUnsplit->pnode->nodeNum=nodeNum;
nodeCounter=0;




/* Enter loop that splits nodes */
NodesToSplit=1;
while(NodesToSplit)
{

/*****************************************
* Split node
*****************************************/
splitterDriver(currUnsplit->pnode,DandO->splitFinal,DandO->nodeProbMethod,DandO);



#ifdef DEBUG 
Rprintf(" splitFinal: changeRSS=%lf \n",DandO->splitFinal->changeRSS);
Rprintf(" splitFinal: splitVar=%d \n",DandO->splitFinal->splitVar);
Rprintf(" splitFinal: splitPoint=%lf \n",DandO->splitFinal->splitPoint);
Rprintf(" splitFinal: thetaL=%lf \n",DandO->splitFinal->thetaL);
Rprintf(" splitFinal: thetaR=%lf \n",DandO->splitFinal->thetaR);
Rprintf(" splitFinal: numL=%lf \n",DandO->splitFinal->numL);
Rprintf(" splitFinal: numR=%lf \n",DandO->splitFinal->numR);
#endif /* DEBUG */

if(DandO->successful_split)
{
/*****************************************
* Create left child node
*****************************************/
nodeCounter++;
currUnsplit->pnode->leftChild=DandO->nodeVEC[nodeCounter];
writeLeftChild(currUnsplit->pnode->leftChild,currUnsplit->pnode,DandO->splitFinal,DandO);

lastUnsplit=addToUnsplit(lastUnsplit,currUnsplit->pnode->leftChild,DandO);
nodeNum++;
currUnsplit->pnode->leftChild->nodeNum=nodeNum;


/*****************************************
* Create right child node
*****************************************/
nodeCounter++;
currUnsplit->pnode->rightChild=DandO->nodeVEC[nodeCounter];
writeRightChild(currUnsplit->pnode->rightChild,currUnsplit->pnode,DandO->splitFinal,DandO);
lastUnsplit=addToUnsplit(lastUnsplit,currUnsplit->pnode->rightChild,DandO);
nodeNum++;
currUnsplit->pnode->rightChild->nodeNum=nodeNum;


if((currUnsplit!=NULL)){

updateNodeMC(currUnsplit->pnode,DandO->splitFinal,DandO);


}

} /* matches "if(DandO->successful_split)"*/
else{

}

prevUnsplit=currUnsplit;
currUnsplit=currUnsplit->nextUnsplit;

if(currUnsplit==NULL)
NodesToSplit=0;

free(prevUnsplit);
}  
 
DandO->numNodes=nodeNum;



/* return root node address */
return(rootNode);

}






/********************************************************************************************************
* FUNCTION: "mallocNode"
*
* ARGUMENTS: void
*            
*             
* RETURNS: pointer to a node with all members allocated memory
*
* PURPOSE: 
* 
********************************************************************************************************/



struct node* mallocNode(struct DatAndOpt *DandO)
{
struct node *pnode;

 
pnode=(struct node*) malloc(sizeof(struct node));
pnode->probInNode=(double*) malloc((DandO->n)*sizeof(double));
pnode->upperLim=(double*) malloc((DandO->p)*sizeof(double));
pnode->lowerLim=(double*) malloc((DandO->p)*sizeof(double));

return(pnode);

}




/********************************************************************************************************
* FUNCTION: "writeLeftChild"
*
* ARGUMENTS:  
*            
*             
* RETURNS: void- writes the results of a split corresponding the left child node to
*                a node structure
*
* PURPOSE: 
* 
**********************************************************************************************************/

void writeLeftChild(struct node *LChild,struct node *parent,struct split *splitRes,struct DatAndOpt *DandO)
{
int i;

/* Need fill in some info */
parent->splitVar=splitRes->splitVar;
parent->splitPoint=splitRes->splitPoint;
parent->changeRSS=splitRes->changeRSS;



LChild->parent=parent;
parent->leftChild=LChild;
LChild->depth=parent->depth+1;
LChild->numInNode=splitRes->numL;
LChild->theta=splitRes->thetaL;
LChild->parentNodeNum=parent->nodeNum;

LChild->changeRSS=-99;
LChild->splitVar=-99;
LChild->splitPoint=-99;

for(i=0;i<(DandO->n);i++)
LChild->probInNode[i]=splitRes->pL[i];


for(i=0;i<(DandO->p);i++)
{
LChild->lowerLim[i]=parent->lowerLim[i];
LChild->upperLim[i]=parent->upperLim[i];
}

LChild->upperLim[parent->splitVar]=parent->splitPoint;


}




/********************************************************************************************************
* FUNCTION: "writeRightChild"
*
* ARGUMENTS:  
*            
*             
* RETURNS: void- writes the results of a split corresponding the right child node to
*                a node structure
*
* PURPOSE: 
* 
**********************************************************************************************************/

void writeRightChild(struct node *RChild,struct node *parent,struct split *splitRes,struct DatAndOpt *DandO)
{
int i;

RChild->parent=parent;
parent->rightChild=RChild;
RChild->depth=parent->depth+1;
RChild->numInNode=splitRes->numR;
RChild->theta=splitRes->thetaR;
RChild->parentNodeNum=parent->nodeNum;


RChild->changeRSS=-99;
RChild->splitVar=-99;
RChild->splitPoint=-99;


for(i=0;i<(DandO->n);i++)
RChild->probInNode[i]=splitRes->pR[i];


for(i=0;i<(DandO->p);i++)
{
RChild->lowerLim[i]=parent->lowerLim[i];
RChild->upperLim[i]=parent->upperLim[i];
}

RChild->lowerLim[parent->splitVar]=parent->splitPoint;


}


 

/********************************************************************************************************
* FUNCTION: "addToUnsplit"
*
* ARGUMENTS:  LastListElement  pointer to last elemtent 
*             pnode            pointer to new node
*             
* RETURNS: address to last element of unsplit list   
*                 
*
* PURPOSE: This function returns a pointer to the last element on the unsplit list,
*          if a "node" is added to the list a pointer to this is returned,
*          otherwise the former last-element is returned
**********************************************************************************************************/
struct unSplitList* addToUnsplit(struct unSplitList *LastListElement,struct node *pnode,struct DatAndOpt *DandO)
{ 
struct unSplitList *uS=NULL;

uS=LastListElement;

if((pnode->depth<(DandO->maxDepth))&&(pnode->numInNode>=(DandO->minSplit)))
{
/* OK to split this node */
uS=(struct unSplitList*) malloc(sizeof(struct unSplitList));
LastListElement->nextUnsplit=uS;
uS->pnode=pnode;
uS->nextUnsplit=NULL;
}

return(uS);

}


/********************************************************************************************************
* FUNCTION: "treeToMatrix"
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
void treeToMatrix(struct node* rootNode,double **treeMatrix,double **lowerLimit,double **upperLimit,struct DatAndOpt *DandO)
{
int AddToMatrix,rowIndx;
struct unSplitList *currNode;
struct unSplitList *prevNode;
struct unSplitList *lastNode;




  
currNode=(struct unSplitList*) malloc(sizeof(struct unSplitList));
currNode->pnode=rootNode;
lastNode=currNode;


AddToMatrix=1;rowIndx=0;
while(AddToMatrix)
{
/* add currNode->pnode to matrices */
nodeToArray(currNode->pnode,*(treeMatrix+rowIndx),*(lowerLimit+rowIndx),*(upperLimit+rowIndx),DandO);


rowIndx++;
if(currNode->pnode->leftChild!=NULL)
{
/* add left child to list */ 
lastNode->nextUnsplit=(struct unSplitList*) malloc(sizeof(struct unSplitList));

lastNode=lastNode->nextUnsplit;
lastNode->pnode=currNode->pnode->leftChild;
lastNode->nextUnsplit=NULL;
}

if(currNode->pnode->rightChild!=NULL)
{
/* add right child to list*/
lastNode->nextUnsplit=(struct unSplitList*) malloc(sizeof(struct unSplitList));
lastNode=lastNode->nextUnsplit;
lastNode->pnode=currNode->pnode->rightChild;
lastNode->nextUnsplit=NULL;
}

prevNode=currNode;
currNode=currNode->nextUnsplit;


if(currNode==NULL)
AddToMatrix=0;

free(prevNode);
} /* matches "while(AddToMatrix)" */


}




/********************************************************************************************************
* FUNCTION: "nodeToArray"
*
* ARGUMENTS:   pnode               pointer to node to be written to row of various matrices 
*              nodeArray           vector to which various node info is written
*              lowerLimit          vector of lower limits of p-dimensional rectangles defining the node
*              upperLimit          vector of upper limits of p-dimensional rectangles defining the node
*             
* RETURNS:  void- info on node written to row of 3 different matrices
*                 
* PURPOSE: This is a help function for the function "treeToMatrix"
**********************************************************************************************************/
void nodeToArray(struct node* pnode,double *nodeArray,double *lowerLimit,double *upperLimit,struct DatAndOpt *DandO)
{
int i;

 
nodeArray[0]=(double) pnode->nodeNum;

nodeArray[1]=(double) pnode->parentNodeNum;

nodeArray[2]=(double) pnode->depth;

nodeArray[3]=(double) pnode->splitVar;

nodeArray[4]=(double) pnode->splitPoint;

nodeArray[5]=(double) pnode->changeRSS;

nodeArray[6]=(double) pnode->theta;

nodeArray[7]=(double) pnode->numInNode;


for(i=0;i<DandO->p;i++)
{
lowerLimit[i]=pnode->lowerLim[i];
upperLimit[i]=pnode->upperLim[i];
}



}









/********************************************************************************************************
* FUNCTION: "buildFirstTree"
*
* ARGUMENTS:   
**********************************************************************************************************/
struct node* buildFirstTree(struct DatAndOpt *DandO)
{
struct node *rootNode;
int i;
double h=0.0;

rootNode=DandO->nodeVEC[0];

DandO->numNodes=1;
rootNode->parent=NULL;
rootNode->leftChild=NULL;
rootNode->rightChild=NULL;
rootNode->depth=1;
rootNode->nodeNum=1;
rootNode->splitVar=-99;
rootNode->parentNodeNum=-99;
rootNode->theta=meanW(DandO->y,DandO->n,DandO->weights,h);

rootNode->numInNode=(double) DandO->n;
DandO->firstTreeBuild=0;

rootNode->splitPoint=-99;      

rootNode->priorRSS=0.0;
rootNode->postRSS=0.0;

for(i=0;i<DandO->n;i++){
rootNode->priorRSS+=(DandO->y[i]*DandO->y[i]*DandO->weights[i]);
rootNode->postRSS+=((DandO->y[i]-rootNode->theta)*(DandO->y[i]-rootNode->theta)*DandO->weights[i]);
}  
rootNode->changeRSS=rootNode->priorRSS-rootNode->postRSS;         


return(rootNode);

}

#undef FUDGE
#undef DEBUG
