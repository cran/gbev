#include <R.h>
#include <math.h>
#include <stdio.h>
#include "dataDef.h"
#include "computeGradient.h"
#include "sampleMCdata.h"


void compute_gradient(struct DatAndOpt *DandO,int numberOfTrees)
{

switch (DandO->SplitFunction){


case 1:  /* least squares boosting */

gradientLS_A(DandO);
break;



case 2: /* log-likelihood boosting */

gradientLS_A_WC_logit(DandO);
break;





}

}



void gradientLS_A(struct DatAndOpt *DandO)
{
int i;

sampleMCdata(DandO);

for(i=0;i<(DandO->n);i++)
DandO->yRes[i]=DandO->y[i]-DandO->yPred[i];


}

void gradientLS_A_WC_logit(struct DatAndOpt *DandO)
{
int i;
sampleMCdata(DandO);

for(i=0;i<(DandO->n);i++)
DandO->yRes[i]=DandO->y[i]-DandO->yPred[i];

}




