/********************************************************************************************************
* File contains declarations of some help functions 
*
* 
*
*
********************************************************************************************************/



double min2num(double a,double b);
double max2num(double a,double b);
void bootsample(int n,double *vec,long *idum);
int sampleOne(int n,long *idum,int *removedIndx,int numRemoved);
void sampleWOR(int n,int size,int *res,long *idum);
void copySplit(struct split *splitF,struct split *splitS,struct DatAndOpt *DandO);
void sampleBAG(int n,int size,double *weights,int *oobIndicator,long *idum);
double computeOOBerror(int n,int size,int *oobIndicator,double *yPred,double *yPrevPred,double *y);



double* gaussRanVec(int n,long *idum);
double* repD(int n,double s);
int* repI(int n,double s);
void MultByConstD(double *vec,int n,double s);
void AddConstD(double *vec,int n,double s);
double mean(double *vec,int n);
double meanW(double *vec,int n,double *weights,double lambda);
void EquateVectorsD(double *vecW,double *vecR,int n);
double maxD(double *vec,int n);
double minD(double *vec,int n);








