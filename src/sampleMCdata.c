/********************************************************************************************************************
* 
*   
* File contains function "sampleMCdata". 
*
* This function samples MC-data from the conditional density of X given W (the "true" covarite given the "observed").
* The sampling is to be done prior to a tree-build. The storage for the sampled data is pointed to by
* structure "DatAndOpt", which also contains parameter for number of MC-samples ("DatAndOpt.mc"). 
* "DatAndOpt" contains a number of other parameters needed to generate data.
*
* NOTE: tis assumed that storage for the MC-data has already been allocated. 
* 
********************************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

#include "dataDef.h"
#include "sampleMCdata.h"
#include "auxFunctions.h"

/*#define DEBUG*/


void sampleMCdata(struct DatAndOpt *DandO)
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
Rprintf(" sampleMCdata-function: normVec[k]=%lf \n",normVec[k]);
#endif /*DEBUG*/ 

}


indx=i*DandO->mc+j; /* row number in DandO.x_mc */

for(k=0;k<DandO->NumberME;k++)
{
DandO->x_mc[indx][k]=DandO->condMeanComp[i*DandO->NumberComponents+compVec[j]][k]; /* conditional mean, given component="compVec[j]" */

/* multiply "normVec" by cholesky factorization of var-covar matrix of the conditional density */
/* of the data at the sampled component */
for(m=0;m<DandO->NumberME;m++)
{
DandO->x_mc[indx][k]+=DandO->cholCov[compVec[j]][m+(DandO->NumberME*k)]*normVec[m];
}
}

}
 
#ifdef DEBUG
#endif
}

free(normVec);
free(compVec);

}


/********************************************************************************************************
* FUNCTION: "AllocSampleData"
*
* ARGUMENTS:   
*             
* RETURNS:   
*                 
* PURPOSE:  This function allocates space to store sampling data, used in MC computations of node probability 
*           belonging. 
*
**********************************************************************************************************/
void AllocSampleData(struct DatAndOpt *DandO)
{

int i,k;

/*********************
* "DandO.ratio_pw" 
**********************/ 
DandO->ratio_pw=(double**) malloc(DandO->n*sizeof(double*));

for(i=0;i<(DandO->n);i++)
DandO->ratio_pw[i]=(double*) malloc(DandO->NumberComponents*sizeof(double));



/*************************
* "DandO.condMeanComp" 
**************************/ 
DandO->condMeanComp=(double**) malloc(DandO->NumberComponents*DandO->n*sizeof(double*));

for(i=0;i<(DandO->NumberComponents*DandO->n);i++)
DandO->condMeanComp[i]=(double*) malloc(DandO->NumberME*sizeof(double));

/*************************
* "DandO.cholCov" 
**************************/ 
DandO->cholCov=(double**) malloc(DandO->NumberComponents*sizeof(double*));

for(k=0;k<DandO->NumberComponents;k++)
DandO->cholCov[k]=(double*) malloc(DandO->NumberME*DandO->NumberME*sizeof(double));


/*************************
* "DandO.cholCov" 
**************************/ 
DandO->x_mc=(double**) malloc(DandO->n*DandO->mc*sizeof(double*));

for(i=0;i<(DandO->n*DandO->mc);i++)
DandO->x_mc[i]=(double*) malloc(DandO->NumberME*sizeof(double));


/*************************
* "DandO.node_mc"  
**************************/ 
DandO->node_mc=(int*) malloc(DandO->n*DandO->mc*sizeof(int));

}


/********************************************************************************************************
* FUNCTION: "InitializeSampleDensities"
*
* ARGUMENTS:  n          number of observations  
*             pME        number of covariates measured with error
*             mc         number of mc-samples
*             numComp    number of components in mixture model for covariates measured with error
*             compProb   a vector, where compProb[(i-1)*numComp+j] "probability w_i from component j"
*             condMean   vector of length "n*numComp*pME". condMean[(i-1)*numComp*pME+pME*(j-1)+k] is the 
*                        conditional mean of x[i,k] conditional on it coming from j-th mixture component.
*             cholCov    vector of length "numComp*pME*pME". This stores the cholesky factorization of the 
*                        variance covariance matrix for the density of X conditional on it coming from given component.                         
*                        cholCov[(j-1)*pME*pME+1....pME] is the first row of the cholesky factorization matrix "C"
*                        such that C%*%t(C) gives the desired conditional variance-covariance matrixs 

*
* 
* RETURNS:   
*                 
* PURPOSE:  This function takes vectors from R and writes their contents into vectors used in sampling
*           MC data for computing node probabilities.
**********************************************************************************************************/
void InitializeSampleDensities(int n,int pME,int mc,int numComp,double *compProb,double *condMean,double *cholCov,struct DatAndOpt *DandO)
{
int i,j,k;

for(i=0;i<n;i++){
for(j=0;j<numComp;j++){
DandO->ratio_pw[i][j]=compProb[i*(numComp)+j];
}
}

for(j=0;j<numComp;j++){
for(i=0;i<pME;i++){ 
for(k=0;k<pME;k++){
DandO->cholCov[j][i*pME+k]=cholCov[j*(pME*pME)+i*pME+k];
}
}
}

for(i=0;i<n;i++){
for(j=0;j<numComp;j++){
for(k=0;k<pME;k++){
DandO->condMeanComp[i*numComp+j][k]=condMean[i*numComp*pME+j*pME+k];
}
}
}

}

/********************************************************************************************************
* FUNCTION: "InitializeSampleData"
*
* ARGUMENTS:   
*             
* RETURNS:   
*                 
* PURPOSE:  This function allocates space to store sampling data, used in MC computations of node probability 
*           belonging. It also...
*
**********************************************************************************************************/

void InitializeSampleData(struct DatAndOpt *DandO)
{
int i,j,k;

/* Allocate space for quantities*/

for(i=0;i<(DandO->mc*DandO->n);i++)
DandO->node_mc[i]=1; 


for(i=0;i<(DandO->NumberComponents*DandO->n);i++)
{
for(k=0;k<DandO->NumberME;k++)
DandO->condMeanComp[i][k]=0;
}


for(i=0;i<(DandO->n);i++)
{
for(k=0;k<DandO->NumberComponents;k++)
DandO->ratio_pw[i][k]=1.0;
}

for(k=0;k<DandO->NumberComponents;k++)
{
for(i=0;i<DandO->NumberME;i++){
for(j=0;j<DandO->NumberME;j++){
DandO->cholCov[k][i*DandO->NumberME+j]=0;
if(i==j)
{
DandO->cholCov[k][i*DandO->NumberME+j]=.7;
}

}
}
}


} /* end-of "InitializeSampleData()" */


/********************************************************************************************************
* FUNCTION: "sample"
*
* ARGUMENTS:   sampleSize  sample size
*              n           length of prob
*              prob        vector of probabilities
*              results     vector to store results
*              idum        random seed
*             
* RETURNS:  results[j]=j with probability prob[j]
*                 
* PURPOSE: 
*
**********************************************************************************************************/


void sample(int sampleSize,int n,double* prob,int* results,long *idum)
{
int i,go,count;
double unif;
double* cumulative;

cumulative=(double*) malloc((n+1)*sizeof(double));

cumulative[0]=0;
for(i=1;i<=n;i++)
{
cumulative[i]=cumulative[i-1]+prob[i-1];
}
cumulative[0]=-.1;
cumulative[n]=1.1;

for(i=0;i<sampleSize;i++)
{



unif=unif_rand();
#ifdef DEBUG 
Rprintf(" sample-function: unif=%lf \n",unif);
#endif /*DEBUG*/ 
go=1;count=0;
while(go)
{

if((unif>cumulative[count])&&(unif<=cumulative[count+1]))
{
results[i]=count;
go=0;
}
count++;
}

}


free(cumulative);
}





/*
* "sampleMCdataTEST" is used to test whether sampling in "sampleMCdata" is correct.
*                    This function is identical, only that it writes the results to 
*                    the vector "res" which is allocated in R.
*/


void sampleMCdataTEST(struct DatAndOpt *DandO,double *res)
{
int i,j,k,m,indx;

int *compVec;
double *normVec;

compVec=(int*) malloc(DandO->mc*sizeof(int)); 
normVec=(double*) malloc(DandO->NumberME*sizeof(double));  /* "DandO.NumberME"= # covariates with error */ 

Rprintf(" Allocated \n");

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
Rprintf(" sampleMCdataTEST-function: normVec[k]=%lf \n",normVec[k]);
#endif /*DEBUG*/ 
}


indx=i*DandO->mc+j; /* row number in DandO.x_mc */

for(k=0;k<DandO->NumberME;k++)
{
DandO->x_mc[indx][k]=DandO->condMeanComp[i*DandO->NumberComponents+compVec[j]][k]; /* conditional mean, given component="compVec[j]" */

/* multiply "normVec" by cholesky factorization of var-covar matrix of the conditional density */
/* of the data at the sampled component */
for(m=0;m<DandO->NumberME;m++)
{
DandO->x_mc[indx][k]+=DandO->cholCov[compVec[j]][m+(DandO->NumberME*k)]*normVec[m];
}
}

}
 
#ifdef DEBUG
#endif
}

free(normVec);
free(compVec);

Rprintf(" got this far \n");
 

if(1)
{

for(i=0;i<DandO->n;i++){
for(j=0;j<DandO->mc;j++){
indx=i*DandO->mc+j;

for(k=0;k<DandO->NumberME;k++){
res[indx*(DandO->NumberME)+k]=DandO->x_mc[indx][k];
}

}
}

} /* matches if(0)*/

}


#undef DEBUG


