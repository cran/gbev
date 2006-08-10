#include <R.h>
#include <math.h>
#include <stdio.h>
#include "dataDef.h"
#include "auxFunctions.h"
#include "nodeprob.h"
#include "splitterLS.h"

/*#define DEBUG*/



/********************************************************************************************************
* FUNCTION: "splitterLS"
*
* ARGUMENTS: pnode          address of node to be split
*            splitFinal     address of "split"-structure to which the results of split are written
*
* RETURNS:   void- all results are written into "splitRes" structure
*
* PURPOSE:   The function splits the node "pnode" (using specifications in "DandO". 
*            The results are written to a "split structure" named splitFinal 
********************************************************************************************************/

void splitterLS(struct node *pnode,struct split *splitFinal,int nodeProbMethod,struct DatAndOpt *DandO)
{


int i,j,covar,first_split_complete,successful_split=0;
int FindingValidPoint;
double sPoint;
struct split splitRes;

int countAttempts;



splitRes.pL=DandO->probLeft;
splitRes.pR=DandO->probRight;

/*********************************************
* Compute rss prior to split
*********************************************/
pnode->priorRSS=0.0;
for(i=0;i<DandO->n;i++){

if(DandO->oobIndicator[i]==0)
pnode->priorRSS=pnode->priorRSS+(DandO->yRes[i]*DandO->yRes[i])*DandO->weights[i];

}


/*******************************************
* Remove fit correponding node to be split
*******************************************/
for(i=0;i<DandO->n;i++)
DandO->yRes[i]=DandO->yRes[i]+(pnode->theta)*(pnode->probInNode[i]);



/************************************************
* Sample covariates to be searched in splitting
**************************************************/

sampleWOR(DandO->p,DandO->m,DandO->res,&(DandO->idum));


DandO->successful_split=0;  /* is set to "1" if node is successfully split.*/
first_split_complete=0;
/*******************************************
* Find best of "m" covariates
*******************************************/
for(i=0;i<DandO->m;i++)  
{
/* sample covariate*/
  
covar=DandO->res[i];

for(j=0;j<DandO->sPoints;j++)
{
FindingValidPoint=1;

countAttempts=0;

while(FindingValidPoint)
{
countAttempts++;
/**********************************************
* Sample point uniformly from permitted region
***********************************************/
sPoint=pnode->lowerLim[covar]+(pnode->upperLim[covar]-pnode->lowerLim[covar])*unif_rand();
splitRes.splitVar=covar;
splitRes.splitPoint=sPoint;
/**********************************************
* Compute fit at sampled point and covariate
***********************************************/
splitPointLS(sPoint,covar,pnode,DandO->yRes,&splitRes,nodeProbMethod,DandO);



/**********************************************
* Determine whether split is valid
***********************************************/
if((splitRes.numL>=DandO->minBucket)&&(splitRes.numR>=DandO->minBucket)){
FindingValidPoint=0;
successful_split=1;  
}

if(countAttempts>DandO->maxSplitAttempts){
FindingValidPoint=0;  
successful_split=0;
}

}  /* matches "while(FindingValidPoint)" */


if(first_split_complete==0){
if(successful_split==1){
copySplit(splitFinal,&splitRes,DandO);
DandO->successful_split=1;  /* means: atleast one successful split of current node has been made.*/
first_split_complete=1;
}
}
else{
if(successful_split==1){
if(splitRes.changeRSS>splitFinal->changeRSS) {
copySplit(splitFinal,&splitRes,DandO);
}
}
}


}  /* matches for(j=0;j<DandO.sPoints;j++)*/

}  /* matches for(i=0;i<m;i++)*/



/*******************************************
* Update fit 
*******************************************/

for(i=0;i<DandO->n;i++)
DandO->yRes[i]=DandO->yRes[i]-splitFinal->thetaL*(splitFinal->pL[i])-splitFinal->thetaR*(splitFinal->pR[i]);





}






/********************************************************************************************************
* FUNCTION: "splitPointLS"
*
* ARGUMENTS: sPoint         point to compute split
*            varIndx        index of covariate to be split
*            pnode          address of node to be split
*            yRes           =y-fit.y , where fit.y is fitted model upto current split 
*            splitRes       address of a split structure to store results
*            nodeProbMethod control parameter that determines method for computing node-belonging probabilities
*
* RETURNS:   void- all results are written into "splitRes" structure
*
* PURPOSE: This function computes a split, given a point at which to split, and a covariate to split on.
*          The predictions of the left and right child nodes of the split and the probabilites of 
*          falling in nodes, as well as the change in RSS are all computed, these are written to the 
*          "split structure"  splitRes. 
********************************************************************************************************/

void splitPointLS(double sPoint,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,int nodeProbMethod,struct DatAndOpt *DandO)
{

int i;


computeNodeProb(sPoint,varIndx,pnode,yRes,splitRes,nodeProbMethod,DandO);

computeThetaLS(sPoint,varIndx,pnode,yRes,splitRes,DandO);
computeChangeRSS_LS(pnode,yRes,splitRes,DandO);

/* compute expected number of obs in each node*/
splitRes->numL=0.0;
splitRes->numR=0.0;
for(i=0;i<DandO->n;i++)
{
splitRes->numL=splitRes->numL+splitRes->pL[i]*DandO->weights[i];
splitRes->numR=splitRes->numR+splitRes->pR[i]*DandO->weights[i];
} 

}







/********************************************************************************************************
* FUNCTION: "computeThetaLS"
*
* ARGUMENTS: sPoint         point to compute split
*            varIndx        index of covariate to be split
*            pnode          address of node to be split
*            yRes           =y-fit.y , where fit.y is fitted model upto current split 
*            splitRes       address of a split structure to store results
*
* RETURNS: void- writes theta(LeftChild) and theta(RightChild) into "splitRes" structure
*
* NOTE:    Function must be run after "computePrNode"
*  
********************************************************************************************************/

void computeThetaLS(double sPoint,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,struct DatAndOpt *DandO)
{

/*
* This routine entails computing the entries of the following 2x2 matrix equation and solving it. 
*
*   | a |  |c  d||thetaL|
*   |   |= |    ||      |
*   | b |  |e  f||thetaR|
* 
*
*   a=sum(pL*yRes*w)+DandO.lambda*pnode->theta      c= sum(pL*w*pL)+DandO.lambda     d= sum(pL*w*pR) 
*   b=sum(pR*yRes*w)+DandO.lambda*pnode->theta      e= sum(pR*w*pL)                  f= sum(pR*w*pR)+DandO.lambda
*
*
*  The solution is:
*
*   1   |f  -d|| a |  |thetaL|
* ----- |     ||   |= |      |
* fc-ed |-e  c|| b |  |thetaR|  
*  
*/


double a=0,b=0,c=0,d=0,e=0,f=0;
double y,pL,pR,w,detInv;
int i;

for(i=0;i<DandO->n;i++)
{

if(DandO->oobIndicator[i]==0)
{
pL=splitRes->pL[i];
pR=splitRes->pR[i];
y=yRes[i];
w=DandO->weights[i];

/*
ch1+=w;
ch2+=y;
ch3+=pL;
ch4+=pR;
*/


a=a+pL*w*y;
b=b+pR*w*y;
c=c+pL*w*pL;
d=d+pL*w*pR;
e=e+pL*w*pR;
f=f+pR*w*pR;
} /* if(DandO->oobIndicator[i]==0)*/

}


/*
c=c-DandO->lambda;
d=d+DandO->lambda;
e=e+DandO->lambda;
f=f-DandO->lambda;
*/
a=a+(DandO->lambda)*pnode->theta;
b=b+(DandO->lambda)*pnode->theta;
c=c+DandO->lambda;

f=f+DandO->lambda;

detInv=1/(f*c-e*d);




splitRes->thetaL=(f*a-d*b)*detInv;
splitRes->thetaR=(c*b-e*a)*detInv;


} 










/********************************************************************************************************
* FUNCTION: "computeThetaLS_2"
*
* ARGUMENTS: sPoint         point to compute split
*            varIndx        index of covariate to be split
*            pnode          address of node to be split
*            yRes           =y-fit.y , where fit.y is fitted model upto current split 
*            splitRes       address of a split structure to store results
*
* RETURNS: void- writes theta(LeftChild) and theta(RightChild) into "splitRes" structure
*
* NOTE:    Function must be run after "computePrNode"
*  
********************************************************************************************************/

void computeThetaLS_2(double sPoint,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,struct DatAndOpt *DandO)
{

/*
* This routine entails computing the entries of the following 2x2 matrix equation and solving it. 
*
*   | a |  |c  d||thetaL|
*   |   |= |    ||      |
*   | b |  |e  f||thetaR|
* 
*
*   a=sum(pL*yRes*w)     c= sum(pL*w*pL)     d= sum(pL*w*pR) 
*   b=sum(pR*yRes*w)     e= sum(pR*w*pL)     f= sum(pR*w*pR)
*
*
*  The solution is:
*
*   1   |f  -d|| a |  |thetaL|
* ----- |     ||   |= |      |
* fc-ed |-e  c|| b |  |thetaR|  
*  
*/


double a=0,b=0,c=0,d=0,e=0,f=0;
double y,pL,pR,w,detInv;
int i;

for(i=0;i<DandO->n;i++)
{

if(DandO->oobIndicator[i]==0)
{
pL=splitRes->pL[i];
pR=splitRes->pR[i];
y=yRes[i];
w=DandO->weights[i];

a=a+pL*w*y;
b=b+pR*w*y;
c=c+pL*w*pL;
d=d+pL*w*pR;
e=e+pL*w*pR;
f=f+pR*w*pR;
} 
}




detInv=1/(f*c-e*d);

splitRes->thetaL=(f*a-d*b)*detInv/DandO->lambda;
splitRes->thetaR=(c*b-e*a)*detInv/DandO->lambda;


} 










/********************************************************************************************************
* FUNCTION: "computeChangeRSS_LS"
*
* ARGUMENTS: 
*            pnode          address of node to be split
*            yRes           =y-fit.y , where fit.y is fitted model upto current split 
*                            NOTE, however, fit.y MUST NOT contain the prediction of the node
*                            to be split.
*            splitRes       address of a split structure to store results
*
* RETURNS: void- writes change in RSS to splitRes, a split structure
* 
* NOTE:    THIS FUNCTION ONLY MAKES SENSE TO RUN AFTER THE "computePrNode" AND "computeTheta"
*          FUNCTIONS, WONT WORK IF NOT
* 
********************************************************************************************************/

void computeChangeRSS_LS(struct node *pnode,double *yRes,struct split *splitRes,struct DatAndOpt *DandO)
{
int i;double rss,thL,thR,resid;

thL=splitRes->thetaL;
thR=splitRes->thetaR;
rss=0.0;
for(i=0;i< DandO->n;i++)
{

if(DandO->oobIndicator[i]==0)
{
resid=DandO->yRes[i]-thL*(splitRes->pL[i])-thR*(splitRes->pR[i]);
rss=rss+(resid*resid)*DandO->weights[i];
}

}

splitRes->changeRSS=pnode->priorRSS-rss;

} 


/********************************************************************************************************
* FUNCTION: "initializeLS"
*
* ARGUMENTS: 
*           
* RETURNS: void- 
* 
* PURPOSE: Initializes the gradient when using squared error-loss 
* 
********************************************************************************************************/


void initializeLS(struct DatAndOpt *DandO)
{
int i;


for(i=0;i<(DandO->n);i++){

DandO->yPred[i]=0.0;
DandO->yRes[i]=DandO->y[i]-DandO->yPred[i];
}





}


