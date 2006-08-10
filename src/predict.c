#include <math.h>
#include <stdio.h>
#include "dataDef.h"
#include "buildTree.h"
#include "Initialize.h"
#include "predict.h"
#include "sampleMCdata.h"





/********************************************************************************************************
* FUNCTION: "predictForest"
*
* ARGUMENTS: x        covariate vector, the reponse of which is to be predicted
*           topNode   top node of tree, from which x is to be predicted.
*
*
* RETURNS:  returns predicted response at x
*                  
**********************************************************************************************************/


double predictForest(double *x,struct treeList *firstTree,int numberTrees)
{
double pred,px;
struct treeList *treeListPtr;
int i;

double ht;

pred=0;


treeListPtr=firstTree;

for(i=0;i<numberTrees;i++)
{
ht=predictTree(x,treeListPtr->topNode);



pred=pred+ht;

treeListPtr=treeListPtr->nextTree;

}

px=pred;

return(px);
}


void predictForest2(double *x,struct treeList *firstTree,int numberTrees,double *result)
{
double pred;
struct treeList *treeListPtr;
int i;

double ht;

pred=0;


treeListPtr=firstTree;

for(i=0;i<numberTrees;i++)
{
ht=predictTree(x,treeListPtr->topNode);



pred=pred+ht;

treeListPtr=treeListPtr->nextTree;

}

*result=pred;
}


/********************************************************************************************************
* FUNCTION: "predictTree"
*
* ARGUMENTS: x        covariate vector, the reponse of which is to be predicted
*           topNode   top node of tree, from which x is to be predicted.
*
*
* RETURNS:  returns predicted response at x
*                  
**********************************************************************************************************/


double predictTree(double *x,struct node *topNode)
{
int findTerminalNode;
struct node *nodePtr;
double pred;
pred=0;

nodePtr=topNode;
findTerminalNode=1;

while(findTerminalNode)
{

if((nodePtr->splitVar)< -10) /* if node is terminal, then no splitvariable, and splitVar=-99 */
{
/*found terminal node */ 
pred=nodePtr->theta;
findTerminalNode=0;
}
else{

if(x[nodePtr->splitVar]<=(nodePtr->splitPoint))
{ /* go to left child */
nodePtr=nodePtr->leftChild;
}
else{
nodePtr=nodePtr->rightChild;
}


}

} /* matches "while(findTerminalNode)" */


return(pred);
}






/********************************************************************************************************
* FUNCTION: "predictForestObs"
*
* ARGUMENTS: x        covariate vector, the reponse of which is to be predicted
*           topNode   top node of tree, from which x is to be predicted.
*
*
* RETURNS:  returns predicted response at x
*                  
**********************************************************************************************************/


double predictForestObs(double **x_mc,struct treeList *firstTree,int numberTrees,int mc)
{



double pred;
struct treeList *treeListPtr;
int i,k;



pred=0;




for(k=0;k<mc;k++)
{

treeListPtr=firstTree;
for(i=0;i<numberTrees;i++)
{

pred=pred+predictTree(x_mc[k],treeListPtr->topNode);

treeListPtr=treeListPtr->nextTree;

}

}

pred=pred/((double) mc);
return(pred);
}

