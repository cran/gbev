#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dataDef.h"
#include "buildTree.h"
#include "Initialize.h"
#include "sampleMCdata.h"

#include "auxFunctions.h"
#include "freeTree.h"
#include "ToMatrix.h"
#include "computeGradient.h"
#include "independentPred.h"
#include "gbev.h"

/*#define DEBUG*/

void InitializeSampleData(struct DatAndOpt *DandO);
void AllocSampleData(struct DatAndOpt *DandO);





/********************************************************************************************************
* FUNCTION: "gbevLWC"
*
* ARGUMENTS:   y         response vector
*              w         covariate vector
*              weights   weight vector
*              n         number of observations
*              p         number of covariates
*              pME       number of covariates that are measured with error, these must be placed first in "w"
*        numberOfBoosts  number of boosting iterations
*              mc        number of monte carlo samples
*              numComp   number of mixture components in marginal of X
*              compProb  a vector, where compProb[(i-1)*numComp+j] "probability w_i from component j"
*              condMean  vector of length "n*numComp*pME". condMean[(i-1)*numComp*pME+pME*(j-1)+k] is the 
*                        conditional mean of x[i,k] conditional on it coming from j-th mixture component.
*              cholCov   vector of length "numComp*pME*pME". This stores the cholesky factorization of the 
*                        variance covariance matrix for the density of X conditional on it coming from given component.                         
*                        cholCov[(j-1)*pME*pME+1....pME] is the first row of the cholesky factorization matrix "C"
*                        such that C%*%t(C) gives the desired conditional variance-covariance matrixs 
* RETURNS:  void-  
*                 
* PURPOSE: Function provides interface between R and C, i.e. this is the function called from R to 
*          build a tree. To get the estimated tree in matrix form, another function is called
*          A tree is built under assumption that the true covariate vector is mixture normal, and 
*          the measurement error is multivariate gaussian, node probabilities are computed using MC approx 
**********************************************************************************************************/



void gbevLWC(double *y,double *w,double *weights,
                        int *n,int *p,int *pME,int *nboost,
                        double *lambda,int *maxDepth,int *m,int *minSplit,
                        int *minBucket,int *sPoints,int *splitFunc,int *predUpdate,int *mc,
                        int *numComp,double *compProb,double *condMean,
                        double *cholCov,int *numNodes,int *nodePositions,
                        double *nodeParameters,
                        int *treeStart,int *treeEnd,int *totNumNodes,double *oobError,
                        int *recordBoosts,double *intermediatePred,
                        int *maxSplitAttempts,int *ranSeed)
{
int i,rb,maxNumNodes;
struct treeList *treeListPtr;
struct node *topNode;
struct DatAndOpt DandO;
int treeStartIndx,treeEndIndx;

int j,treeCount;

/*******************************************************************
* Need call prior to using R-functions for random variable generation
********************************************************************/
GetRNGstate();

*numNodes=87;



/****************************************
* Data info     
****************************************/
DandO.n=*n;
DandO.p=*p;
DandO.NumberComponents=*numComp;
DandO.NumberME=*pME;
DandO.mc=*mc;
DandO.weights=weights;
DandO.y=y;
DandO.wcWeights=(double*) malloc(DandO.n*sizeof(double));
DandO.yRes=(double*) malloc(DandO.n*sizeof(double));
DandO.yPred=(double*) malloc((DandO.n)*sizeof(double));

DandO.firstTreeBuild=1;



/****************************************
* TREE-building options    
****************************************/
DandO.lambda=*lambda;
DandO.m=*m;              
DandO.minSplit=*minSplit;         
DandO.minBucket=*minBucket;       
DandO.maxDepth=*maxDepth;         
DandO.sPoints=*sPoints;         
DandO.numberOfTrees=(*nboost);
DandO.oobIndicator=(int*) malloc(DandO.n*sizeof(int));
DandO.sizeBag=(int) DandO.n; /* Need pass as parameter */
DandO.oobError=oobError;
DandO.SplitFunction=*splitFunc;
DandO.maxSplitAttempts=*maxSplitAttempts;
DandO.nodeProbMethod=1; /* so far, there is only one such method */
DandO.idum=-(*ranSeed);
DandO.Indep_y_Update=*predUpdate; 


divAllocDandO(&DandO);
divAllocPred(&DandO);  /* allocate space for computing independent update of yPred, independent of fitting */

maxNumNodes=(int) pow(2.0,(double) (DandO.maxDepth));
maxNumNodes=maxNumNodes-1;
*totNumNodes=(maxNumNodes*DandO.numberOfTrees);




/******************************************************************
* Determine if split method requires the use of constraint points
*******************************t**********************************/
DandO.useConstraints=0;

 

/***********************************************************
* Allocate storage for MC data
***********************************************************/
AllocSampleData(&DandO);
InitializeSampleDensities(DandO.n,DandO.NumberME,DandO.mc,DandO.NumberComponents,compProb,condMean,cholCov,&DandO);



/***********************************************************
* Allocate a matrix to store covariates  
***********************************************************/
DandO.w=(double**) malloc(DandO.p*sizeof(double*));
for(i=0;i<DandO.p;i++)
{
DandO.w[i]=&(w[i*DandO.n]);
}



/******************************************
* minimum and maximum of covariate vectors    
*******************************************/
DandO.minW=repD(DandO.p,1.0); /* allocates and initializes to 1.0 */
DandO.maxW=repD(DandO.p,1.0);
 

for(i=0;i<DandO.p;i++)
{
DandO.minW[i]=minD(DandO.w[i],DandO.n);
DandO.maxW[i]=maxD(DandO.w[i],DandO.n);
}


 

/****************************************
*  initialize treeList (seem to need these next four lines, dunno why....)
****************************************/
DandO.firstTree=(struct treeList*) malloc(sizeof(struct treeList));
treeListPtr=DandO.firstTree;
treeListPtr->prevTree=NULL;
treeListPtr->nextTree=NULL;


 

 

 
treeCount=0;
/****************************************
*  Enter  Boosting loop
****************************************/
rb=0; /* 'rb' is used to detemine if prediction at iteration is to be saved*/
treeStartIndx=0;  

sampleMCdata(&DandO);

for(i=0;i<DandO.n;i++)
DandO.oobIndicator[i]=0;


for(i=0;i<*nboost;i++)
{
if(i==0){
sampleMCdata(&DandO);
initializeFit(&DandO);
}
else{
compute_gradient(&DandO,treeCount);
}

treeCount++;


if(DandO.firstTreeBuild==1){
topNode=buildFirstTree(&DandO);
}
else{
topNode=buildTree(&DandO);
}



/***************************************************************************
* Update "yPred" using independent sample, i.e independent of fitting sample
****************************************************************************/
if(DandO.Indep_y_Update==1){
sampleMCdataPred(&DandO);
}
yPredUpdate(&DandO);




/************************
* Write tree to vectors
*************************/
treeEndIndx=treeStartIndx+DandO.numNodes-1;

treeToVectors(topNode,treeStartIndx,treeEndIndx,nodePositions,nodeParameters,
              &DandO,*totNumNodes,DandO.numNodes);
treeStart[treeCount-1]=treeStartIndx;
treeEnd[treeCount-1]=treeEndIndx;
treeStartIndx=treeEndIndx+1;


if(recordBoosts[rb]==(treeCount-1)){
for(j=0;j<DandO.n;j++)
intermediatePred[j+(rb*DandO.n)]=DandO.yPred[j];
rb++;
}


} /* End boosting loop */




/******************************
* Free allocated space
******************************/
freeDatAndOpt2(&DandO);


/*******************************************************************
* Need call after using R-functions for random variable generation
********************************************************************/
PutRNGstate();


}











/***********************************************************
*  
* "divAllocDandO": Allocates various quantities for DandO
*
************************************************************/


void divAllocDandO(struct DatAndOpt *DandO)
{
int nMaxNodes,i;
nMaxNodes=(int) pow(2.0,(double) (DandO->maxDepth+1));
nMaxNodes=nMaxNodes-1;



DandO->splitRes=(struct split*) malloc(sizeof(struct split)); 
DandO->splitRes->pL=(double*) malloc(DandO->n*sizeof(double));
DandO->splitRes->pR=(double*) malloc(DandO->n*sizeof(double));
 
DandO->splitFinal=(struct split*) malloc(sizeof(struct split));  
DandO->splitFinal->pL=(double*) malloc(DandO->n*sizeof(double));
DandO->splitFinal->pR=(double*) malloc(DandO->n*sizeof(double));

DandO->res=(int*) malloc(DandO->p*sizeof(int));                 
DandO->probLeft=(double*) malloc(DandO->n*sizeof(double));
DandO->probRight=(double*) malloc(DandO->n*sizeof(double));

/* Allocate nodes to be used in all of the tree fitting */
DandO->nodeVEC=(struct node**) malloc(nMaxNodes*sizeof(struct node*)); 

for(i=0;i<nMaxNodes;i++)
DandO->nodeVEC[i]=(struct node*) mallocNode(DandO);


}



void divFreeDandO(struct DatAndOpt *DandO)
{
int nMaxNodes,i;
nMaxNodes=(int) pow(2.0,(double) (DandO->maxDepth+1));
nMaxNodes=nMaxNodes-1;


free(DandO->splitRes->pL);
free(DandO->splitRes->pR); 
free(DandO->splitRes);
 
free(DandO->splitFinal->pL); 
free(DandO->splitFinal->pR); 
free(DandO->splitFinal); 

free(DandO->res);                 
free(DandO->probLeft); 
free(DandO->probRight); 


/* Allocate nodes to be used in all of the tree fitting*/

for(i=0;i<nMaxNodes;i++)
freeNode(DandO->nodeVEC[i]);


free(DandO->nodeVEC);


}



/*#undef DEBUG*/

