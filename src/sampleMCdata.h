float gasdev(long *idum);
void sample(int sampleSize,int n,double* prob,int* results,long *idum);
void sampleMCdata(struct DatAndOpt *DandO);
void InitializeSampleData(struct DatAndOpt *DandO);
void AllocSampleData(struct DatAndOpt *DandO);
void InitializeSampleDensities(int n,int pME,int mc,int numComp,double *compProb,double *condMean,double *cholCov,struct DatAndOpt *DandO);
void sampleMCdataTEST(struct DatAndOpt *DandO,double *res);



