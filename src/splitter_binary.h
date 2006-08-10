void splitter_binary(struct node *pnode,struct split *splitFinal,int nodeProbMethod,struct DatAndOpt *DandO);
void splitPoint_binary(double sPoint,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,int nodeProbMethod,struct DatAndOpt *DandO);
void computeChangeRSS_binary(struct node *pnode,double *yRes,struct split *splitRes,struct DatAndOpt *DandO);
void computeTheta_binary(double sPoint,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,struct DatAndOpt *DandO);
void initialize_binary(struct DatAndOpt *DandO);
void InitRootNode_binary(struct node *pnode,struct DatAndOpt *DandO);

