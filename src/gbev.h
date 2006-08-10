void gbevLWC(double *y,double *w,double *weights,
                        int *n,int *p,int *pME,int *nboost,
                        double *lambda,int *maxDepth,int *m,int *minSplit,
                        int *minBucket,int *sPoints,int *splitFunc,int *predUpdate,int *mc,
                        int *numComp,double *compProb,double *condMean,
                        double *cholCov,int *numNodes,int *nodePositions,
                        double *nodeParameters,
                        int *treeStart,int *treeEnd,int *totNumNodes,double *oobError,
                        int *recordBoosts,double *intermediatePred,
                        int *maxSplitAttempts,int *ranSeed);

void divAllocDandO(struct DatAndOpt *DandO);
void divFreeDandO(struct DatAndOpt *DandO);

