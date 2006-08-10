void splitterLS(struct node *pnode,struct split *splitFinal,int nodeProbMethod,struct DatAndOpt *DandO);
void splitPointLS(double sPoint,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,int nodeProbMethod,struct DatAndOpt *DandO);
void computeThetaLS(double sPoint,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,struct DatAndOpt *DandO);
void computeChangeRSS_LS(struct node *pnode,double *yRes,struct split *splitRes,struct DatAndOpt *DandO);
void initializeLS(struct DatAndOpt *DandO);
void computeThetaLS_2(double sPoint,int varIndx,struct node *pnode,double *yRes,struct split *splitRes,struct DatAndOpt *DandO);





