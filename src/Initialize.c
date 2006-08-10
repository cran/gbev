#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dataDef.h"
#include "auxFunctions.h"
#include "Initialize.h"
#include "splitter_binary.h"
#include "splitterLS.h"




#define FUDGE 0.01



/****************************************************************************************************
* FUNCTION:   InitRootNode
* 
* ARGUMENTS:  n             number of observations
*
* RETURNS :   void
*
*
* PURPOSE:   This function initializes the root node of a tree.  
*             
****************************************************************************************************/
void InitiRootNode(struct node *pnode,struct DatAndOpt *DO)
{
int i;


/***********************
* Initialize root node 
************************/
pnode->depth=1;
pnode->parent=NULL;

pnode->numInNode=DO->n;
pnode->parentNodeNum=0;
pnode->splitVar=-70;

for(i=0;i<DO->n;i++)
pnode->probInNode[i]=1.0;


for(i=0;i<DO->p;i++)
{
pnode->lowerLim[i]=DO->minW[i]-.001;
pnode->upperLim[i]=DO->maxW[i]+.001;
}


/******************************************************************
* Initialize root node according to different regression types
*******************************************************************/
switch (DO->SplitFunction){

case 1:   /* Least squares regression */
pnode->theta=meanW(DO->yRes,DO->n,DO->weights,DO->lambda);  
break;

case 2:  /* Binary regression using log-likelihood */
InitRootNode_binary(pnode,DO);
break;



}





}


#undef FUDGE

void initializeFit(struct DatAndOpt *DandO)
{
switch (DandO->SplitFunction){

case 1:  /* least squares boosting */
initializeLS(DandO);
break;

case 2: /* binary regression using log-likelihood */
initialize_binary(DandO);
break;



}

}



