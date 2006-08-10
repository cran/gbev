#include <R.h>
#include <math.h>
#include <stdio.h>
#include "dataDef.h"
#include "splitter.h"
#include "splitterLS.h"
#include "splitter_binary.h"
#include "auxFunctions.h"

#define SQRT2 1.41421356237









/********************************************************************************************************
* FUNCTION: "splitterDriver"
*
* ARGUMENTS: pnode          address of node to be split
*            splitFinal     address of "split"-structure to which the results of split are written
*            nodeProbMethod controls which method is used to compute the node-belonging probabilities
*
* RETURNS:   void- all results are written into "splitRes" structure
*
* PURPOSE:   The function splits the node "pnode" (using specifications in "DandO". 
*            The results are written to a "split structure" named splitFinal 
********************************************************************************************************/


void splitterDriver(struct node *pnode,struct split *splitFinal,int nodeProbMethod,struct DatAndOpt *DandO)
{
 
switch (DandO->SplitFunction){


case 1: /* least squares splitting */
splitterLS(pnode,splitFinal,nodeProbMethod,DandO);
break;

case 2: /* binary regression using log-likelihood */
splitter_binary(pnode,splitFinal,nodeProbMethod,DandO); 
break;

default: 
splitterLS(pnode,splitFinal,nodeProbMethod,DandO);
break;

}


}  










