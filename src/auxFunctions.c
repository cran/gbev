#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dataDef.h"
#include "auxFunctions.h"

/*#define DEBUG*/

#define SQRT2 1.41421356237




/********************************************************************************************************
* FUNCTION: "min2num"
*
* ARGUMENTS: a
*            b
*             
* RETURNS: minimum of two numbers a and b
* 
********************************************************************************************************/

double min2num(double a,double b)
{
 

if(a<=b)
{
return(a);
}
else{
return(b);
}

}


/********************************************************************************************************
* FUNCTION: "max2num"
*
* ARGUMENTS: a
*            b
*             
* RETURNS: maximum of two numbers a and b
* 
********************************************************************************************************/

double max2num(double a,double b)
{
 

if(a>=b)
{
return(a);
}
else{
return(b);
}


}



/********************************************************************************************************
* FUNCTION: "bootsample"
*
* ARGUMENTS: n         number of observations 
*            vec       where to write results, this is a vector of doubles, of length "n" (not size)
*            idum      seed of random number generator
*
* RETURNS:  void- results are written to "vec"
* 
* PURPOSE:  The function draws a bootstrap sample, that is a sample with replacement is drawn from the
*           the cases of the original data. The results are written to the (NOTE) double vector "vec",
*           and the results are stored as counts of the number of times given data points are 
*           included in bootstrap sample. This form is convenient when passing the results as weights
*            
********************************************************************************************************/
void bootsample(int n,double *vec,long *idum)
{
int i,aux;
double unif,nAux;

nAux=(double) n;

for(i=0;i<n;i++)
vec[i]=0.0;


/* Sample a vector of uniforms*/
for(i=0;i<n;i++)
{

unif=unif_rand();
#ifdef DEBUG
Rprintf("bootstrap-function: unif=%lf \n",unif);
#endif /*DEBUG*/

aux=(int) nAux*unif;
vec[aux]=vec[aux]+1.0;
}


}





/********************************************************************************************************
* FUNCTION: "sampleOne"
*
* ARGUMENTS: n             number of observations 
*            idum          seed of random number generator
*            removedIndx   index of length n (0's and 1's) of whether observation has been removed, i.e cannot 
*                          be sampled.
*            numRemoved    number removed.
*
* RETURNS:  returns on integer sampled from 1 to n, but not from set with removedIndx=1.
* 
* PURPOSE:  This is a help function for performing sampling without replacement
* 
********************************************************************************************************/
int sampleOne(int n,long *idum,int *removedIndx,int numRemoved)
{
int *sampleVec;
int count,i;double nAux;

sampleVec=(int*) malloc((n-numRemoved)*sizeof(int));

count=0;
for(i=0;i<n;i++)
{
if(removedIndx[i]==0)
{
sampleVec[count]=i;
count++;
}
}

nAux=(double) (n-numRemoved);

nAux=nAux*unif_rand();
#ifdef DEBUG
Rprintf("sampleOne-function: unif=%lf \n",nAux);
#endif /*DEBUG*/


count=(int) nAux;
count=sampleVec[count];
free(sampleVec);
return(count);
}




/********************************************************************************************************
* FUNCTION: "sampleWOR"
*
* ARGUMENTS: n      the length of the index vector being sampled wo replacement
*            size   number of elements sampled without replacement
*            res    a vector to store the results of the sampling
*            idum   a random seed
*
* RETURNS:   void- the sampled elements are placed in "res"
*
* PURPOSE: This function samples "size" elements without replacement from the 
*          vector (0,...,(n-1)). The results are placed in "res" 
********************************************************************************************************/
void sampleWOR(int n,int size,int *res,long *idum)
{
int i,k,j,nElements;
int *sampVec=(int*) malloc(n*sizeof(int));


double nA,h;

for(i=1;i<=n;i++)
sampVec[i-1]=i-1;

nElements=n;
nA=(double) n;
 
for(i=0;i<size;i++)
{

h=nA*unif_rand(); 
#ifdef DEBUG
Rprintf("sampleWOR-function: unif=%lf \n",h);
#endif /*DEBUG*/

k=(int) h;
res[i]=sampVec[k];

for(j=(k+1);j<nElements;j++)
{
sampVec[j-1]=sampVec[j];
}
nElements=nElements-1;
nA=nA-1;
}

free(sampVec);

}



/********************************************************************************************************
* FUNCTION: "copySplit"
*
* ARGUMENTS: splitF
*            splitS
*
* RETURNS:  void- function copies 
* 
* PURPOSE:  This function copies the elements in splitS into splitF,
*           both split-structures
*
********************************************************************************************************/
void copySplit(struct split *splitF,struct split *splitS,struct DatAndOpt *DandO)
{
int i;

splitF->changeRSS=splitS->changeRSS;
splitF->splitVar=splitS->splitVar;
splitF->splitPoint=splitS->splitPoint;
splitF->thetaL=splitS->thetaL;
splitF->thetaR=splitS->thetaR;
splitF->numL=splitS->numL;
splitF->numR=splitS->numR;

for(i=0;i<(DandO->n);i++)
{
splitF->pL[i]=splitS->pL[i];
splitF->pR[i]=splitS->pR[i];
}

}








/********************************************************************************************************
* FUNCTION: "sampleBAG"
*
* ARGUMENTS: n      the length of the index vector being sampled wo replacement
*            size   number of elements sampled without replacement
*            weights a vector of weights to be used in fitting, will be set to 1 if in BAG, 0 if outside
*       oobIndicator same as 1-"weights", an integer vector
*            res    a vector to store the results of the sampling
*            idum   a random seed
*
* RETURNS:   void- the sampled elements are placed in "res"
*
* PURPOSE: This function samples the "BAG" used for fitting of a tree
********************************************************************************************************/

void sampleBAG(int n,int size,double *weights,int *oobIndicator,long *idum)
{
int *res;
int i;
res=(int*) malloc(n*sizeof(int));

sampleWOR(n,size,res,idum);

for(i=0;i<n;i++)
{

oobIndicator[i]=1;
}

for(i=0;i<size;i++)
{

oobIndicator[res[i]]=0; /* data used in fitting, not outside bag */
}



free(res);
}




/********************************************************************************************************
* FUNCTION: "computeOOBerror"
*
* ARGUMENTS: n      the length of the index vector being sampled wo replacement
*            size   number of elements sampled without replacement
*            weights a vector of weights to be used in fitting, will be set to 1 if in BAG, 0 if outside
*       oobIndicator same as 1-"weights", an integer vector
*            res    a vector to store the results of the sampling
*            idum   a random seed
*
* RETURNS:   void- the sampled elements are placed in "res"
*
* PURPOSE: This function samples the "BAG" used for fitting of a tree
********************************************************************************************************/

double computeOOBerror(int n,int size,int *oobIndicator,double *yPred,double *yPrevPred,double *y)
{
double res;
int i;
double count;

count=0.0;
res=0;
for(i=0;i<n;i++)
{

if(oobIndicator[i]==1){
res+=y[i]*(log(yPred[i]+.0001)-log(yPrevPred[i]+.0001))+(1-y[i])*(log(1-yPred[i]+.0001)-log(1-yPrevPred[i]+.0001));

count=count+1.0;

}

}

res=res/count;

return(res);

}



/********************************************************************************************************
* FUNCTION: "computeChangeOOBerror"
*
* ARGUMENTS: n      the length of the index vector being sampled wo replacement
*            size   number of elements sampled without replacement
*            weights a vector of weights to be used in fitting, will be set to 1 if in BAG, 0 if outside
*       oobIndicator same as 1-"weights", an integer vector
*            res    a vector to store the results of the sampling
*            idum   a random seed
*
* RETURNS:   void- the sampled elements are placed in "res"
*
* PURPOSE: This function samples the "BAG" used for fitting of a tree
********************************************************************************************************/

double computeChangeOOBerror(int n,int size,int *oobIndicator,double *yResCurrent,double *yResPrevious)
{
double res,count;
int i;


res=0;count=0;
for(i=0;i<n;i++)
{

if(oobIndicator[i]==1)
{
count++;
res+=yResPrevious[i]*yResPrevious[i]-yResCurrent[i]*yResCurrent[i];
}

}

res=res/(count);

return(res);

}





/*
* FUNCTION: "gaussRanVec" 
* 
* ARGUMENTS: n        length of generated vector
*
* RETURNS:  Pointer to the start of a vector of N(0,1)  
*
*
*/

double* gaussRanVec(int n,long *idum)
{
double *res; 
int i;

res=(double*) malloc(n*sizeof(double));

for(i=0;i<n;i++)
{
res[i]=norm_rand();
#ifdef DEBUG
Rprintf(" gaussRanVec-function: res[i]=%lf \n",res[i]);
#endif /*DEBUG*/

}

return(res);
}

/****************************************************************************************
* FUNCTION: "repD" 
* 
* ARGUMENTS: n        length of generated vector
*            s        value placed in each entry of returned vector
*
* RETURNS:  Pointer to the start of a vector of "doubles" of length n, and all entries equal to "s"  
*
* NOTE:   This is the equivalent of the the R function "rep"
*
****************************************************************************************/

double* repD(int n,double s)
{
double *res; 
int i;

res=(double*) malloc(n*sizeof(double));

for(i=0;i<n;i++)
{
res[i]=s;
}

return(res);
}

/****************************************************************************************
* FUNCTION: "repI" 
* 
* ARGUMENTS: n        length of generated vector
*            s        value placed in each entry of returned vector
*
* RETURNS:  Pointer to the start of a vector of "integers" of length n, and all entries equal to "s"  
*
* NOTE:   This is the equivalent of the the R function "rep"
*
****************************************************************************************/

int* repI(int n,double s)
{
int *res; 
int i;

res=(int*) malloc(n*sizeof(int));

for(i=0;i<n;i++)
{
res[i]=s;
}

return(res);
}



/****************************************************************************************
* FUNCTION: "MultByConstD" 
* 
* ARGUMENTS: vec      pointer to double vector
*            n        length of generated vector
*            s        value to multiply by
*          
*
* RETURNS:  void  
*
* PURPOSE:  multiplies the elements in "vec" by the constant "s", "vec" must be doubles
*
* 
*
****************************************************************************************/

void MultByConstD(double *vec,int n,double s)
{
 
int i;

for(i=0;i<n;i++)
{
vec[i]=s*vec[i];
}

}



/****************************************************************************************
* FUNCTION: "AddConstD" 
* 
* ARGUMENTS: vec      pointer to double vector
*            n        length of generated vector
*            s        value to multiply by
*          
*
* RETURNS:  void  
*
* PURPOSE:  adds the constant "s" constant to the elements in "vec", "vec" must be doubles
*
* 
*
****************************************************************************************/

void AddConstD(double *vec,int n,double s)
{
 
int i;

for(i=0;i<n;i++)
{
vec[i]=vec[i]+s;
}

}


/****************************************************************************************
* FUNCTION: "mean" 
* 
* ARGUMENTS: vec      pointer to double vector
*            n        length of generated vector
*                     
*          
* RETURNS:  mean of "vec"
* 
*
****************************************************************************************/

double mean(double *vec,int n)
{
 
int i;
double res=0.0;

for(i=0;i<n;i++)
{
res=res+vec[i];
}

return(res/((double) n));
}


/****************************************************************************************
* FUNCTION: "meanW" 
* 
* ARGUMENTS: vec      pointer to double vector
*            n        length of generated vector
*            weights  pointer to vector of weights
*          
* RETURNS:  weighted mean of "vec"
* 
*
****************************************************************************************/

double meanW(double *vec,int n,double *weights,double lambda)
{
 
int i;
double res=0.0;
double sumW=0.0;

for(i=0;i<n;i++)
{
sumW=sumW+weights[i];
res=res+vec[i]*weights[i];
}

return(res/(sumW+lambda));
}


/****************************************************************************************
* FUNCTION:  "EquateVectorsD" 
* 
* ARGUMENTS: vecW     pointer to double vector to write over
*            vecR     pointer to double vector to write "vecW" over with
*            n        length of vectors         
*          
* RETURNS:   void
*
* PURPOSE:   The function sets the elements of "vecW" equal to the elements of "vecR" 
* 
*
****************************************************************************************/

void EquateVectorsD(double *vecW,double *vecR,int n)
{
 
int i;

for(i=0;i<n;i++)
{
vecW[i]=vecR[i];
}

}




/****************************************************************************************
* FUNCTION:  "minD" 
* 
* ARGUMENTS: vec     pointer to double vector of which to find minimum
*            n       length of vectors         
*          
* RETURNS:   minimum of "vec"
*
* PURPOSE:   The function finds and returns the minimum of a double vector, "vec" 
* 
*
****************************************************************************************/

double minD(double *vec,int n)
{
int i;double min;

min=vec[0];

for(i=1;i<n;i++)
{

if(vec[i]<min)
min=vec[i];

}

return(min);
}



/****************************************************************************************
* FUNCTION:  "maxD" 
* 
* ARGUMENTS: vec     pointer to double vector of which to find maximum
*            n       length of vectors         
*          
* RETURNS:   minimum of "vec"
*
* PURPOSE:   The function finds and returns the maximum of a double vector, "vec" 
* 
*
****************************************************************************************/

double maxD(double *vec,int n)
{
int i;
double max;

max=vec[0];

for(i=1;i<n;i++)
{

if(vec[i]>max)
max=vec[i];

}

return(max);
}


