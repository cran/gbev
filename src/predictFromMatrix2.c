#include <stdio.h>
#include "dataDef.h"
#include "predictFromMatrix2.h"


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
* FUNCTION: "predictEnsembleR2_mult"
*
* ARGUMENTS:  
*             x                 vector of values to compute partial dependence. Is (nPred*pd) large.
*             nPred             number of points to be predicted 
*             p                 total number of covariates
*             startEnsemble     Tree index from which to start model "averaging"
*             endEnsemble       Tree index from which to end model "averaging"
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

void predictEnsembleR2_mult(double *x,int *nPred,int *p,int *startEnsemble,
                         int *endEnsemble,int *nEndTrees,int *nodePositions,double *nodeParameters,
                         int *treeStart,int *treeEnd,int *totNumNodes,double *res)
{

 
int j,k,sEnsemble,eEnsemble;
double h;

for(k=0;k<(*nEndTrees);k++)
{

if(k==0){
sEnsemble=*startEnsemble;
eEnsemble=endEnsemble[k];
}
if(k>0){
sEnsemble=endEnsemble[k-1]+1;
eEnsemble=endEnsemble[k];
}

for(j=0;j<*nPred;j++)
{ 


if(k==0)
res[j+(*nPred)*k]=0;


if(k>0)
res[j+(*nPred)*k]=res[j+(*nPred)*(k-1)];



h=predictEnsemble2(&(x[j*(*p)]),&sEnsemble,&eEnsemble,nodePositions,nodeParameters,treeStart,treeEnd,totNumNodes);

res[j+(*nPred)*k]+=h;

}
}

}







/********************************************************************************************************
* FUNCTION: "predictEnsemble"
*
* ARGUMENTS:  
*             x                 vector of values to compute partial dependence. Is (nPred*pd) large.
*             nPred             number of points to be predicted 
*             p                 total number of covariates
*             startEnsemble     Tree index from which to start model "averaging"
*             endEnsemble       Tree index from which to end model "averaging"
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

void predictEnsembleR2(double *x,int *nPred,int *p,int *startEnsemble,int *endEnsemble,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes,double *res)
{

 
int j;



for(j=0;j<*nPred;j++)
{ 


res[j]=predictEnsemble2(&(x[j*(*p)]),startEnsemble,endEnsemble,nodePositions,nodeParameters,treeStart,treeEnd,totNumNodes);

}


}







/********************************************************************************************************
* FUNCTION: "predictEnsemble"
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

double predictEnsemble2(double *x,int *startEnsemble,int *endEnsemble,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes)
{


double pred;
 
int i;
pred=0;



/* NOTE: startEnsemble=1 corresponds to the first tree, and so on. */

for(i=((*startEnsemble)-1);i<((*endEnsemble));i++)
{

pred=pred+predictTreeFromVec2(x,nodePositions,nodeParameters,treeStart[i],treeEnd[i],*totNumNodes);
}


return(pred);
}



double predictTreeFromVec2(double *x,int *nodePositions,double *nodeParameters,int treeStart,int treeEnd,int totalNumberNodes)
{
int FindTerminalNode;
int indx;
double pred;

pred=0;

indx=treeStart; /* index corresponding root node of tree */

FindTerminalNode=1;
while(FindTerminalNode)
{
if(nodePositions[indx+SPLIT_VAR*totalNumberNodes]<-10) {
/* this gives the index of split variable, =-99 if no split, so */
/* <-10 node is terminal */
pred=nodeParameters[indx+THETA*totalNumberNodes]; 
FindTerminalNode=0;
}
else{
/* current node was split, go to off-spring*/

if(x[nodePositions[indx+SPLIT_VAR*totalNumberNodes]]<=nodeParameters[indx+SPLIT_POINT*totalNumberNodes]){
/* go to left child-node*/
indx=nodePositions[indx+LEFT_CHILD*totalNumberNodes];
}
else{
/* go to left right-node*/
indx=nodePositions[indx+RIGHT_CHILD*totalNumberNodes];
}

}
} /* matches while(FindTerminalNode) */
return(pred);
}













/********************************************************************************************************
* FUNCTION: "predictEnsembleAll"
*
* ARGUMENTS:  
*             x                 vector of values to compute partial dependence. Is (nPred*pd) large.
*             nPred             number of points to be predicted 
*             p                 total number of covariates
*             startEnsemble     Tree index from which to start model "averaging"
*             endEnsemble       Tree index from which to end model "averaging"
*             nodePositions     for each node gives: node number, parent node, child nodes, tree number
*             nodeParameters    for each node gives: prediction, changeRSS, upper and lower limits
*             treeStart         treeStart[i] gives the start node of i-th tree 
*             treeEnd           treeEnd[i]   gives the last node (in nodePositions and nodeParameters) of tree "i"
*
* RETURNS: 
*
* PURPOSE:  This function returns the predictions of "x" for all sizes of the ensemble.
*           ie if there are "nTrees" in ensemble, a prediction for only using the first
*           tree, only the first 2 trees, the first 3 trees, and so on are provided
*  
* NOTE:    In use with R, the vectors (not DandO) must be allocated from R 
*          (max length from "depth" and "number trees")
********************************************************************************************************/

void predictEnsembleAll(double *x,int *nPred,int *p,int *startEnsemble,int *endEnsemble,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes,double *res)
{

 
int j,i,nTrees,endTree,k;

/* number of trees*/
nTrees=(*endEnsemble)-(*startEnsemble)+1;

for(j=0;j<*nPred;j++)
{ 

for(i=(*startEnsemble-1);i<(*endEnsemble);i++)
{
endTree=i+1;
k=i-(*startEnsemble-1);

if(k==0){
res[j*nTrees+k]=predictTreeFromVec2(&(x[j*(*p)]),nodePositions,nodeParameters,treeStart[i],treeEnd[i],*totNumNodes);
}
else{
res[j*nTrees+k]=res[j*nTrees+k-1]+predictTreeFromVec2(&(x[j*(*p)]),nodePositions,nodeParameters,treeStart[i],treeEnd[i],*totNumNodes);
}

}
}


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
