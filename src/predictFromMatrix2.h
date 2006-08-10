double predictEnsemble2(double *x,int *startEnsemble,int *endEnsemble,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes);
double predictTreeFromVec2(double *x,int *nodePositions,double *nodeParameters,int treeStart,int treeEnd,int totalNumberNodes);
void predictEnsembleR2(double *x,int *nPred,int *p,int *startEnsemble,int *endEnsemble,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes,double *res);
void predictEnsembleR2_mult(double *x,int *nPred,int *p,int *startEnsemble,int *endEnsemble,int *nEndTrees,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes,double *res);


