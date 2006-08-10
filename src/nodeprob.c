#include <R.h>
#include <math.h>
#include <stdio.h>
#include "dataDef.h"
#include "nodeprob.h"
#include "auxFunctions.h"


/********************************************************************************************************
* FUNCTION: "computeNodeProb"
*
* ARGUMENTS: sPoint         point to compute split
*            varIndx        index of covariate to be split
*            pnode          address of node to be split
*            yRes           =y-fit.y , where fit.y is fitted model upto current split 
*            splitRes       address of a split structure to store results
*
* RETURNS: void- writes Prob(in Left child) and Prob(in Right child) into "splitRes" structure
* 
* NOTE:    THE FUNCTION ASSUMES THAT SPACE FOR THESE TWO VECTORS HAS ALREADY BEEN RESERVED 
********************************************************************************************************/

void computeNodeProb(double sP,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,int method,struct DatAndOpt *DandO)
{

switch (method){

case 1: /* compute probability using MC-sampling */
nodeProbMC(sP,varIndx,pnode,yRes,splitRes,DandO);
break;

case 2: /* compute probability using numerical integration, not implemented yet */
break;


} 

} 




/********************************************************************************************************
* FUNCTION: "nodeProbMC"
*
* ARGUMENTS: sPoint         point to compute split
*            varIndx        index of covariate to be split
*            pnode          address of node to be split
*            yRes           =y-fit.y , where fit.y is fitted model upto current split 
*            splitRes       address of a split structure to store results
*
* RETURNS: void- writes Prob(in Left child) and Prob(in Right child) into "splitRes" structure
*          
* 
* NOTE:    THE FUNCTION ASSUMES THAT SPACE FOR THESE TWO VECTORS HAS ALREADY BEEN RESERVED 
*  
* DETAILS: Computes these node probabilities using MC-sampling. The relevant MC-sample is in
*          the global DatAndOpt structure DandO. The function makes use of the variable
*          "DandO.node_mc",   DandO.node_mc[j] is the node-number corresponding the j-th
*          MC-sampled data point, for the current state of the tree.        
********************************************************************************************************/


void nodeProbMC(double sP,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,struct DatAndOpt *DandO)
{

int i,j,numLeft,numRight;



for(i=0;i<DandO->n;i++)
{

numLeft=0;
numRight=0;
for(j=0;j<DandO->mc;j++)
{

if(DandO->node_mc[i*DandO->mc+j]==pnode->nodeNum) /* pnode should point to node being split*/
{

if(varIndx<(DandO->NumberME)){  
/* asks if covariate measured with error, the error measured covariates */
/* are assumed to be placed first in "DandO.w" */

if(DandO->x_mc[i*DandO->mc+j][varIndx]<=sP){
numLeft++;
}
else{
numRight++;
}
} 

#ifndef DEBUG
else{ /* corresponds if covariate not measured with error*/

if(DandO->w[varIndx][i]<=sP){
numLeft++;
}
else{
numRight++;
}

} /* else... */
#endif

} /* if(DandO->node_mc[i*DandO->mc+j]==pnode->nodeNum) */

} /* matches "for(j=0;j.." */

splitRes->pL[i]=((double) numLeft)/((double) DandO->mc);
splitRes->pR[i]=((double) numRight)/((double) DandO->mc);


} /* matches "for(i=0;i<DandO.n;i++)" */



} /* end-of "nodeProbMC"-function */





/********************************************************************************************************
* FUNCTION: "updateNodeMC"
*
* ARGUMENTS: parent
*            splitRes
*
* RETURNS: void- writes updates DandO.node_mc
*          
* 
* NOTE:    Updates if "DandO.nodeProbMethod==1", else falls through
*          and "DandO.nodeProbMethod==1" signifies that Crude MC is used to compute node-belonging probs.
*  
* DETAILS:    
********************************************************************************************************/

void updateNodeMC(struct node *parent,struct split *splitRes,struct DatAndOpt *DandO)
{
int i,j,indx;
double splitValue;
int splitVariable;
int parentNodeNumber;


parentNodeNumber=parent->nodeNum;
splitValue=splitRes->splitPoint;
splitVariable=splitRes->splitVar;

if(DandO->nodeProbMethod==1)
{

for(i=0;i<(DandO->n);i++)
{
for(j=0;j<(DandO->mc);j++)
{

indx=i*(DandO->mc)+j;
if(DandO->node_mc[indx]==parentNodeNumber)
{

/*********************************
* If an error measured covariate
**********************************/
if(splitVariable<DandO->NumberME)
{

if(DandO->x_mc[indx][splitVariable]<=splitValue){
DandO->node_mc[indx]=parent->leftChild->nodeNum;
}
else{
DandO->node_mc[indx]=parent->rightChild->nodeNum;
}

} /* if(splitVariable<DandO->pME) */
/******************************
* If covariate error free
*******************************/
else{

if(DandO->w[splitVariable][i]<=splitValue){
DandO->node_mc[indx]=parent->leftChild->nodeNum;
}
else{
DandO->node_mc[indx]=parent->rightChild->nodeNum;
}

} /* else.. */


}
}
}


} 

}





