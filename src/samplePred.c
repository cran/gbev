#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "predictFromMatrix2.h"
#include "dataDef.h"
#include "sampleMCdata.h"
#include "auxFunctions.h"




/********************************************************************************************************
* FUNCTION: "samplePredData"
*
* ARGUMENTS:  n          number of observations 
*             p                 total number of covariates
*             startEnsemble     Tree index from which to start model "averaging"
*             endEnsemble       Tree index from which to end model "averaging"
*             nodePositions     for each node gives: node number, parent node, child nodes, tree number
*             nodeParameters    for each node gives: prediction, changeRSS, upper and lower limits
*             treeStart         treeStart[i] gives the start node of i-th tree 
*             treeEnd           treeEnd[i]   gives the last node (in nodePositions and nodeParameters) of tree "i"
 
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


void samplePredData(int *method,double *w,int *nw,int *mc,int *p,int *pME,
                    int *numComp,double *compProb,double *condMean,double *cholCov,
                    int *startEnsemble,int *endEnsemble,int *num_TreeEnd,int *nodePositions,double *nodeParameters,
                    int *treeStart,int *treeEnd,int *totNumNodes,double *res)
{


int i,j,k,m,indx,indx1,indx2,indxComp,indxCompCov;
long idum;
int *nPred;
int nPredV;
int *compVec;
int sEnsemble,eEnsemble;
double *normVec;
double *x;
double **x_mc;
double *res_mc;

idum=-8471244;

nPredV=2;
nPred=&nPredV;

compVec=(int*) malloc((*mc)*sizeof(int)); 
normVec=(double*) malloc((*pME)*sizeof(double));  /* "DandO.NumberME"= # covariates with error */ 
x_mc=(double**) malloc(((*mc)*(*nw))*sizeof(double*)); 

for(i=0;i<((*mc)*(*nw));i++)
x_mc[i]=(double*) malloc((*pME)*sizeof(double));


for(i=0;i<(*nw);i++)
{

/* sample "*mc" components for individ "i", place in "compVec" */

sample(*mc,*numComp,&compProb[(i*(*pME))],compVec,&idum);


for(j=0;j<(*mc);j++)
{

/* generate N(0,1)'s */
for(k=0;k<(*pME);k++){

normVec[k]=norm_rand();
}


indx=i*(*mc)+j; /* row number in x_mc */

/* "indxComp" indexes where the "compVec[j]" component conditional mean starts*/
indxComp=i*(*pME)*(*numComp)+(*pME)*compVec[j]; 

/* "indxCompCov" indexes where the cholesky factorization of component"compVec[j]" starts*/
indxCompCov=compVec[j]*(*pME)*(*pME);

for(k=0;k<(*pME);k++)
{


x_mc[indx][k]=condMean[indxComp+k]; /* conditional mean, given component="compVec[j]" */

/* multiply "normVec" by cholesky factorization of var-covar matrix of the conditional density */
/* of the data at the sampled component */
for(m=0;m<(*pME);m++)
{
x_mc[indx][k]+=cholCov[indxCompCov+k*(*pME)+m]*normVec[m];
}
}

}
 

}






/***********************************************************************************************
* x is now filled up with sampled values for error measured covariates, and the error free ones
*************************************************************************************************/
x=(double*) malloc(((*p)*(*mc)*(*nw))*sizeof(double));

if(x==NULL)
Rprintf("Couldnt allocate x \n");



for(i=0;i<(*nw);i++){


indx1=i*(*p)*(*mc);
indx2=i*(*mc);
for(j=0;j<(*mc);j++)
{


for(k=0;k<(*pME);k++){
x[indx1+j*(*p)+k]=x_mc[indx2+j][k];


}
if((*p)>(*pME))
{
for(k=(*pME);k<(*p);k++)
x[indx1+j*(*p)+k]=w[i*(*p)+k];
} 

}
}

/* Now predict at sampled values */
nPredV=(*mc)*(*nw);


res_mc=(double*) malloc(((*mc)*(*nw))*sizeof(double));

if(res_mc==NULL)
Rprintf("Couldnt allocate res_mc \n");




/*If there are multiple "treeEnd" values, then these are looped through, starting here */


/* If "(*num_TreeEnd)>1" and logistic regression is used, this function won't work correctly*/

if(((*num_TreeEnd)>1)&&((*method)>2))
Rprintf(" C-Function 'samplePredData', does not work with multiple endpoints for methods other than L2-boost \n"); 

for(k=0;k<(*num_TreeEnd);k++)
{


if(k==0)
predictEnsembleR2(x,&nPredV,p,startEnsemble,&endEnsemble[k],nodePositions,nodeParameters,treeStart,treeEnd,totNumNodes,res_mc);

if(k>0)
{
sEnsemble=(endEnsemble[k-1]+1);
eEnsemble=endEnsemble[k];

predictEnsembleR2(x,&nPredV,p,&sEnsemble,&eEnsemble,nodePositions,nodeParameters,treeStart,treeEnd,totNumNodes,res_mc);
}


for(i=0;i<(*nw);i++){

if(k==0)
res[i+(*nw)*k]=0;

if(k>0)
res[i+(*nw)*k]=res[i+(*nw)*(k-1)];

for(j=0;j<(*mc);j++){

if((*method)<3) /* least-squares regression */
res[i+(*nw)*k]+=res_mc[i*(*mc)+j];


if((*method)==3) /* logistic regression */
res[i+(*nw)*k]+=1/(1+exp(-res_mc[i*(*mc)+j]));

if((*method)>3)
Rprintf(" C-function 'samplePred.c' has not been implemented for this method. \n");

}


}


} /* matches "for(k=0;k<(*num_TreeEnd);k++)"*/

for(i=0;i<((*nw));i++){
for(k=0;k<(*num_TreeEnd);k++){
res[i+(*nw)*k]=res[i+(*nw)*k]/((double) (*mc));
}
}

free(res_mc);
for(i=0;i<((*mc)*(*nw));i++)
free(x_mc[i]);

free(x_mc);
free(x);
#ifdef DEBUGGER
#endif

free(normVec);
free(compVec);



}






















/********************************************************************************************************
* FUNCTION: "sampleData"
*
* ARGUMENTS:  n          number of observations 
*             p                 total number of covariates
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
* PURPOSE:  This function samples data from its conditional density (given observed covariates),
*           and places results in "res", for use in R
**********************************************************************************************************/


void sampleData(double *w,int *nw,int *mc,int *p,int *pME,int *numComp,double *compProb,double *condMean,double *cholCov,double *res)
{


int i,j,k,m,indx,indx1,indx2,indxComp,indxCompCov;
long idum;
int *compVec;
double *normVec;
double **x_mc;

idum=-8471244;


compVec=(int*) malloc((*mc)*sizeof(int)); 
normVec=(double*) malloc((*pME)*sizeof(double));  /* "DandO.NumberME"= # covariates with error */ 
x_mc=(double**) malloc(((*mc)*(*nw))*sizeof(double*)); 

for(i=0;i<((*mc)*(*nw));i++)
x_mc[i]=(double*) malloc((*pME)*sizeof(double));


for(i=0;i<(*nw);i++)
{

/* sample "*mc" components for individ "i", place in "compVec" */

sample(*mc,*numComp,&compProb[(i*(*pME))],compVec,&idum);


for(j=0;j<(*mc);j++)
{

/* generate N(0,1)'s */
for(k=0;k<(*pME);k++){

normVec[k]=norm_rand();
}


indx=i*(*mc)+j; /* row number in x_mc */

/* "indxComp" indexes where the "compVec[j]" component conditional mean starts*/
indxComp=i*(*pME)*(*numComp)+(*pME)*compVec[j]; 

/* "indxCompCov" indexes where the cholesky factorization of component"compVec[j]" starts*/
indxCompCov=compVec[j]*(*pME)*(*pME);

for(k=0;k<(*pME);k++)
{


x_mc[indx][k]=condMean[indxComp+k]; /* conditional mean, given component="compVec[j]" */

/* multiply "normVec" by cholesky factorization of var-covar matrix of the conditional density */
/* of the data at the sampled component */
for(m=0;m<(*pME);m++)
{
x_mc[indx][k]+=cholCov[indxCompCov+k*(*pME)+m]*normVec[m];
}
}

}
 

}






/***********************************************************************************************
* "res" will now filled up with sampled values for error measured covariates, and the error free ones
*************************************************************************************************/



for(i=0;i<(*nw);i++){


indx1=i*(*p)*(*mc);
indx2=i*(*mc);
for(j=0;j<(*mc);j++)
{


for(k=0;k<(*pME);k++){
res[indx1+j*(*p)+k]=x_mc[indx2+j][k];


}
if((*p)>(*pME))
{
for(k=(*pME);k<(*p);k++)
res[indx1+j*(*p)+k]=w[i*(*p)+k];
} 

}
}




for(i=0;i<((*mc)*(*nw));i++)
free(x_mc[i]);

free(x_mc);

#ifdef DEBUGGER
#endif

free(normVec);
free(compVec);



}


/********************************************************************************************************
* FUNCTION: "samplePredData"
*
* ARGUMENTS:  n          number of observations 
*             p                 total number of covariates
*             startEnsemble     Tree index from which to start model "averaging"
*             endEnsemble       Tree index from which to end model "averaging"
*             nodePositions     for each node gives: node number, parent node, child nodes, tree number
*             nodeParameters    for each node gives: prediction, changeRSS, upper and lower limits
*             treeStart         treeStart[i] gives the start node of i-th tree 
*             treeEnd           treeEnd[i]   gives the last node (in nodePositions and nodeParameters) of tree "i"
 
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


void samplePredDataTest(int *method,double *w,int *nw,int *mc,int *p,int *pME,
                    int *numComp,double *compProb,double *condMean,double *cholCov,
                    int *startEnsemble,int *endEnsemble,int *num_TreeEnd,int *nodePositions,double *nodeParameters,
                    int *treeStart,int *treeEnd,int *totNumNodes,double *res,double *x)
{


int i,j,k,m,indx,indx1,indx2,indxComp,indxCompCov;
long idum;
int *nPred;
int nPredV;
int *compVec;
int sEnsemble,eEnsemble;
double *normVec;

double **x_mc;
double *res_mc;

idum=-8471244;

nPredV=2;
nPred=&nPredV;

compVec=(int*) malloc((*mc)*sizeof(int)); 
normVec=(double*) malloc((*pME)*sizeof(double));  /* "DandO.NumberME"= # covariates with error */ 
x_mc=(double**) malloc(((*mc)*(*nw))*sizeof(double*)); 

for(i=0;i<((*mc)*(*nw));i++)
x_mc[i]=(double*) malloc((*pME)*sizeof(double));


for(i=0;i<(*nw);i++)
{

/* sample "*mc" components for individ "i", place in "compVec" */

sample(*mc,*numComp,&compProb[(i*(*pME))],compVec,&idum);


for(j=0;j<(*mc);j++)
{

/* generate N(0,1)'s */
for(k=0;k<(*pME);k++){

normVec[k]=norm_rand();
}


indx=i*(*mc)+j; /* row number in x_mc */

/* "indxComp" indexes where the "compVec[j]" component conditional mean starts*/
indxComp=i*(*pME)*(*numComp)+(*pME)*compVec[j]; 

/* "indxCompCov" indexes where the cholesky factorization of component"compVec[j]" starts*/
indxCompCov=compVec[j]*(*pME)*(*pME);

for(k=0;k<(*pME);k++)
{


x_mc[indx][k]=condMean[indxComp+k]; /* conditional mean, given component="compVec[j]" */

/* multiply "normVec" by cholesky factorization of var-covar matrix of the conditional density */
/* of the data at the sampled component */
for(m=0;m<(*pME);m++)
{
x_mc[indx][k]+=cholCov[indxCompCov+k*(*pME)+m]*normVec[m];
}
}

}
 

}






/***********************************************************************************************
* x is now filled up with sampled values for error measured covariates, and the error free ones
*************************************************************************************************/


for(i=0;i<(*nw);i++){


indx1=i*(*p)*(*mc);
indx2=i*(*mc);
for(j=0;j<(*mc);j++)
{


for(k=0;k<(*pME);k++){
x[indx1+j*(*p)+k]=x_mc[indx2+j][k];


}
if((*p)>(*pME))
{
for(k=(*pME);k<(*p);k++)
x[indx1+j*(*p)+k]=w[i*(*p)+k];
} 

}
}

/* Now predict at sampled values */
nPredV=(*mc)*(*nw);


res_mc=(double*) malloc(((*mc)*(*nw))*sizeof(double));

if(res_mc==NULL)
Rprintf("Couldnt allocate res_mc \n");




/* If there are multiple "treeEnd" values, then these are looped through, starting here */


/* If "(*num_TreeEnd)>1" and logistic regression is used, this function won't work correctly*/

if(((*num_TreeEnd)>1)&&((*method)>2))
Rprintf(" C-Function 'samplePredData', does not work with multiple endpoints for methods other than L2-boost \n"); 

for(k=0;k<(*num_TreeEnd);k++)
{


if(k==0)
predictEnsembleR2(x,&nPredV,p,startEnsemble,&endEnsemble[k],nodePositions,nodeParameters,treeStart,treeEnd,totNumNodes,res_mc);

if(k>0)
{
sEnsemble=(endEnsemble[k-1]+1);
eEnsemble=endEnsemble[k];

predictEnsembleR2(x,&nPredV,p,&sEnsemble,&eEnsemble,nodePositions,nodeParameters,treeStart,treeEnd,totNumNodes,res_mc);
}


for(i=0;i<(*nw);i++){

if(k==0)
res[i+(*nw)*k]=0;

if(k>0)
res[i+(*nw)*k]=res[i+(*nw)*(k-1)];

for(j=0;j<(*mc);j++){

if((*method)<3) /* least-squares regression */
res[i+(*nw)*k]+=res_mc[i*(*mc)+j];


if((*method)==3) /* logistic regression */
res[i+(*nw)*k]+=1/(1+exp(-res_mc[i*(*mc)+j]));

if((*method)>3)
Rprintf(" C-function 'samplePred.c' has not been implemented for this method. \n");

}


}


}  

for(i=0;i<((*nw));i++){
for(k=0;k<(*num_TreeEnd);k++){
res[i+(*nw)*k]=res[i+(*nw)*k]/((double) (*mc));
}
}

free(res_mc);
for(i=0;i<((*mc)*(*nw));i++)
free(x_mc[i]);

free(x_mc);

#ifdef DEBUGGER
#endif

free(normVec);
free(compVec);



}














