#include <stdio.h>
#include <stdlib.h>
#include "dataDef.h"
#include "partialDepend.h"



#define SPLIT_POINT 0
#define THETA 1
#define CHANGE_RSS 2
#define NUMBER_IN_NODE 3

#define NODE_NUMBER 0
#define PARENT_NODE_NUMBER 1
#define LEFT_CHILD 2
#define RIGHT_CHILD 3
#define SPLIT_VAR 4



/********************************************************************************************************
* FUNCTION: "partialDependenceR"
*
* ARGUMENTS:  
*             x                 vector of values to compute partial dependence. Is (nPred*pd) large.
*             nPred             number of points (of length "pd") at which to eval partial dependence
*             pd                number of variables computing partial dependence for
*             p                 total number of covariates
*             partDepVar        vector of indexes (length "pd") for partial dependence variables
*             startEnsemble     Tree index from which to start model "averaging"
*             endEnsemble       Tree index from which to end model "averaging"
*             nodePositions     for each node gives: node number, parent node, child nodes, tree number
*             nodeParameters    for each node gives: prediction, changeRSS, upper and lower limits
*             treeStart         treeStart[i] gives the start node of i-th tree 
*             treeEnd           treeEnd[i]   gives the last node (in nodePositions and nodeParameters) of tree "i"
*             totNumNodes,
*             res
*
* RETURNS: 
*
* PURPOSE:  Writes partial dependence of points in "x" to "res" 
*  
* NOTE:     Called from R.
********************************************************************************************************/


void partialDependenceR(double *x,int *nPred,int *pd,int *p,int *partDepVar,int *startEnsemble,int *endEnsemble,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes,double *res)
{

int j;

for(j=0;j<*nPred;j++)
res[j]=partialDependence(&(x[j*(*pd)]),pd,p,partDepVar,startEnsemble,endEnsemble,nodePositions,nodeParameters,treeStart,treeEnd,totNumNodes);


}

 



/********************************************************************************************************
* FUNCTION: "partialDependence"
*
* ARGUMENTS:  
*             x                 point to compute partial dependence
*             pd                number of variables computing partial dependence for
*             p                 total number of covariates
*             partDepVar        vector of indexes (length "pd") for partial dependence variables
*             startEnsemble     Tree index from which to start model "averaging"
*             endEnsemble       Tree index from which to end model "averaging"
*             nodePositions     for each node gives: node number, parent node, child nodes, tree number
*             nodeParameters    for each node gives: prediction, changeRSS, upper and lower limits
*             treeStart         treeStart[i] gives the start node of i-th tree 
*             treeEnd           treeEnd[i]   gives the last node (in nodePositions and nodeParameters) of tree "i"
*
* RETURNS:   partial dependence of point "x"
*
* PURPOSE:  
*  
* NOTE:   
********************************************************************************************************/


double partialDependence(double *x,int *pd,int *p,int *partDepVar,int *startEnsemble,int *endEnsemble,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes)
{
int i;
double partDep;

partDep=0.0;

for(i=((*startEnsemble)-1);i<((*endEnsemble));i++)
{
partDep=partDep+partDepTree(x,pd,partDepVar,nodePositions,nodeParameters,treeStart[i],treeEnd[i],*totNumNodes);
}

return(partDep);
}



/********************************************************************************************************
* FUNCTION: "partDepTree"
*
* ARGUMENTS:  
*             x                 point to compute partial dependence
*             pd                number of variables computing partial dependence for
*             p                 total number of covariates
*             partDepVar        vector of indexes (length "pd") for partial dependence variables
*             nodePositions     for each node gives: node number, parent node, child nodes, tree number
*             nodeParameters    for each node gives: prediction, changeRSS, upper and lower limits
*             treeStart         gives the start node of tree 
*             treeEnd           gives the last node (in nodePositions and nodeParameters) of tree 
*             totNumNodes       (adress of) total number of nodes describing ensemble
*
* RETURNS:   partial dependence of point "x" for tree starting at "treeStart" and ending at "treeEnd"
*
* PURPOSE:   
*  
* NOTE:    
********************************************************************************************************/

double partDepTree(double *x,int *pd,int *partDepVar,int *nodePositions,double *nodeParameters,int treeStart,int treeEnd,int totNumNodes)
{
double res,p,valAtSplitPoint=0.0;
int *nodes;
int moreNodes,i,splitVar,childNodeIndx;
int partDepVariable;
int nodeIndx;
int nodeEndIndx;
double *nodeWeights;

res=0;

/* node positions necesary to visit (in argument vectors) */
nodes=(int*) malloc(sizeof(int)*(treeEnd-treeStart+1));
/* weights at node positions necesary to visit */
nodeWeights=(double*) malloc(sizeof(double)*(treeEnd-treeStart+1));

nodeIndx=0;
nodes[nodeIndx]=treeStart;
nodeWeights[nodeIndx]=1.0;
nodeEndIndx=0;

moreNodes=1;
while(moreNodes)
{
splitVar=nodePositions[nodes[nodeIndx]+SPLIT_VAR*totNumNodes];

if(splitVar<-10){
/*************************************************
* Terminal node, add to partial dependence value
**************************************************/
res+=nodeWeights[nodeIndx]*nodeParameters[nodes[nodeIndx]+THETA*totNumNodes];
}
else{
/*************************************************
* Internal node, add to "visit list"
**************************************************/


/* check if variable split on is member of partial dependence set */
partDepVariable=0;
for(i=0;i<(*pd);i++){

if(splitVar==partDepVar[i]){
partDepVariable=1;
valAtSplitPoint=x[i]; /* value of part-dep point of split-variable */
}

}


/* if split on one of part.depend variables add only next (appropriate) node to list, dont update weights */
if(partDepVariable==1)
{


if(valAtSplitPoint<=nodeParameters[nodes[nodeIndx]+SPLIT_POINT*totNumNodes]){ 
/* go left */
nodeEndIndx++;
nodes[nodeEndIndx]=nodePositions[nodes[nodeIndx]+LEFT_CHILD*totNumNodes];
nodeWeights[nodeEndIndx]=nodeWeights[nodeIndx];

}
else{
/* go right */
nodeEndIndx++;
nodes[nodeEndIndx]=nodePositions[nodes[nodeIndx]+RIGHT_CHILD*totNumNodes];
nodeWeights[nodeEndIndx]=nodeWeights[nodeIndx];

}

} /* if(partDepVariable==1) */



/* if not split on one of part.depend variables add next two nodes to list, update weights */
if(partDepVariable==0)
{ 
nodeEndIndx++;
childNodeIndx=nodePositions[nodes[nodeIndx]+LEFT_CHILD*totNumNodes];
nodes[nodeEndIndx]=childNodeIndx;
p=nodeParameters[childNodeIndx+NUMBER_IN_NODE*totNumNodes]/nodeParameters[nodes[nodeIndx]+NUMBER_IN_NODE*totNumNodes];
nodeWeights[nodeEndIndx]=nodeWeights[nodeIndx]*p;

nodeEndIndx++;
childNodeIndx=nodePositions[nodes[nodeIndx]+RIGHT_CHILD*totNumNodes];
nodes[nodeEndIndx]=childNodeIndx;
nodeWeights[nodeEndIndx]=nodeWeights[nodeIndx]*(1-p);
}


} /* else{...internal node... } */
 
if(nodeIndx==nodeEndIndx)
moreNodes=0; 

nodeIndx++;
} /* while(moreNodes) */


free(nodes);
free(nodeWeights);

return(res);
}





#undef SPLIT_POINT  
#undef THETA  
#undef CHANGE_RSS  
#undef NUMBER_IN_NODE  

#undef  NODE_NUMBER  
#undef PARENT_NODE_NUMBER  
#undef LEFT_CHILD  
#undef RIGHT_CHILD  
#undef SPLIT_VAR  

