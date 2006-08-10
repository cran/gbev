/*********************************************************************************************
* This file contains functions for the independent update of DandO->yPred, independent, that is,
* of the sample used in the fitting.
*
*********************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

#include "dataDef.h"
#include "auxFunctions.h"
#include "predict.h"
#include "sampleMCdata.h"


#define FUDGE 0.0001

/*#define DEBUG*/


void divAllocPred(struct DatAndOpt *DandO)
{
int i;


DandO->splitPred=(struct split*) malloc(sizeof(struct split)); 
DandO->splitPred->pL=(double*) malloc(DandO->n*sizeof(double));
DandO->splitPred->pR=(double*) malloc(DandO->n*sizeof(double));
 

DandO->x_mc_pred=(double**) malloc(DandO->n*DandO->mc*sizeof(double*));

for(i=0;i<(DandO->n*DandO->mc);i++)
DandO->x_mc_pred[i]=(double*) malloc(DandO->NumberME*sizeof(double));


DandO->node_mc_pred=(int*) malloc(DandO->n*sizeof(int));       

}


void divFreePred(struct DatAndOpt *DandO)
{
int i;

free(DandO->splitPred->pL);
free(DandO->splitPred->pR);
free(DandO->splitPred);


for(i=0;i<(DandO->n*DandO->mc);i++)
free(DandO->x_mc_pred[i]); 

free(DandO->x_mc_pred);

free(DandO->node_mc_pred);

}






/********************************************************************************************************************
* 
*   "sampleMCdataPred" samples data and places in "x_mc_pred" which is used to get independent update of
*                      yPred, independent that is of the sample used in the fitting.
* 
********************************************************************************************************************/

void sampleMCdataPred(struct DatAndOpt *DandO)
{
int i,j,k,m,indx;

int *compVec;
double *normVec;

compVec=(int*) malloc(DandO->mc*sizeof(int)); 
normVec=(double*) malloc(DandO->NumberME*sizeof(double));  /* "DandO.NumberME"= # covariates with error */ 


for(i=0;i<DandO->n;i++)
{

/* sample DandO.mc components for individ "i", place in "compVec" */
sample(DandO->mc,DandO->NumberComponents,DandO->ratio_pw[i],compVec,&(DandO->idum));


for(j=0;j<DandO->mc;j++)
{

/* generate N(0,1)'s */
for(k=0;k<DandO->NumberME;k++){

normVec[k]=norm_rand();
#ifdef DEBUG
Rprintf(" sampleMCdataPred-function: normVec[k]=%lf \n",normVec[k]);
#endif /*DEBUG*/

}

indx=i*DandO->mc+j; /* row number in DandO.x_mc */

for(k=0;k<DandO->NumberME;k++)
{
DandO->x_mc_pred[indx][k]=DandO->condMeanComp[i*DandO->NumberComponents+compVec[j]][k]; /* conditional mean, given component="compVec[j]" */

/* multiply "normVec" by cholesky factorization of var-covar matrix of the conditional density */
/* of the data at the sampled component */
for(m=0;m<DandO->NumberME;m++)
{
DandO->x_mc_pred[indx][k]+=DandO->cholCov[compVec[j]][m+(DandO->NumberME*k)]*normVec[m];
}
}

}
 
#ifdef DEBUG
#endif
}

free(normVec);
free(compVec);

}








void yPredUpdate(struct DatAndOpt *DandO)
{
int i,j,indx,n,mc,p,pME,k;
double pred,dmc;
double *xp;



xp=(double*) malloc((DandO->p)*sizeof(double));

n=DandO->n;
mc=DandO->mc;
p=DandO->p;
pME=DandO->NumberME;

dmc=(double) mc;


if(DandO->Indep_y_Update==1)
{
for(i=0;i<n;i++){


for(k=pME;k<p;k++)
xp[k]=DandO->w[k][i];

pred=0;
for(j=0;j<mc;j++){
indx=i*mc+j;

for(k=0;k<pME;k++)
xp[k]=DandO->x_mc_pred[indx][k];


pred+=predictTree(xp,DandO->nodeVEC[0]);
}
DandO->yPred[i]+=pred/dmc;

if(DandO->SplitFunction==2)
{
if(DandO->yPred[i]>(1.0-FUDGE))
DandO->yPred[i]=(1.0-FUDGE);

if(DandO->yPred[i]<FUDGE)
DandO->yPred[i]=FUDGE;
}


} 
}

if(DandO->Indep_y_Update==0)
{
 
for(i=0;i<n;i++){


for(k=pME;k<p;k++)
xp[k]=DandO->w[k][i];

pred=0;
for(j=0;j<mc;j++){
indx=i*mc+j;

for(k=0;k<pME;k++)
xp[k]=DandO->x_mc[indx][k];


pred+=predictTree(xp,DandO->nodeVEC[0]);
}
DandO->yPred[i]+=pred/dmc;

if(DandO->SplitFunction==2)
{
if(DandO->yPred[i]>(1.0-FUDGE))
DandO->yPred[i]=(1.0-FUDGE);

if(DandO->yPred[i]<FUDGE)
DandO->yPred[i]=FUDGE;
}

} 

}

free(xp);

#ifdef DEBUG
for(i=0;i<n;i++){
Rprintf(" DandO->yPred[%d]=%lf \n",i,DandO->yPred[i]);
}
#endif /*DEBUG*/

}

#undef FUDGE
#undef DEBUG

