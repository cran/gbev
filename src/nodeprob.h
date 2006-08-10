
void computeNodeProb(double sP,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,int method,struct DatAndOpt *DandO);
void nodeProbMC(double sP,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,struct DatAndOpt *DandO);
void updateNodeMC(struct node *parent,struct split *splitRes,struct DatAndOpt *DandO);

