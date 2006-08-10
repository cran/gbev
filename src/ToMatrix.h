void ensembleToMatrix(struct DatAndOpt *DandO,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes);
void treeToVectors(struct node* rootNode,int startIndx,int endIndx,int *nodePositions,double *nodeParameters,struct DatAndOpt *DandO,int totalNumberNodes,int NumNodesTree);
void printTreeVec(int *nodePositions,double *nodeParameters,struct DatAndOpt *DandO,int totalNumberNodes);




