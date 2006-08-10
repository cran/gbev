double predictTree(double *x,struct node *topNode);
double predictForest(double *x,struct treeList *firstTree,int numberTrees);
double predictForestObs(double **x_mc,struct treeList *firstTree,int numberTrees,int mc);
void predictForest2(double *x,struct treeList *firstTree,int numberTrees,double *res);





