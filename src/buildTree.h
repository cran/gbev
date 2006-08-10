struct node* buildTree(struct DatAndOpt *DO);

struct node* mallocNode(struct DatAndOpt *DO);
void writeLeftChild(struct node *LChild,struct node *parent,struct split *splitRes,struct DatAndOpt *DO);
void writeRightChild(struct node *RChild,struct node *parent,struct split *splitRes,struct DatAndOpt *DO);
struct unSplitList* addToUnsplit(struct unSplitList *LastListElement,struct node *pnode,struct DatAndOpt *DO);
void treeToMatrix(struct node* rootNode,double **treeMatrix,double **lowerLimit,double **upperLimit,struct DatAndOpt *DandO);
void nodeToArray(struct node* pnode,double *nodeArray,double *lowerLimit,double *upperLimit,struct DatAndOpt *DO);
void printTreeMatrix(double **treeMatrix,double **lowerLimit,double **upperLimit,struct DatAndOpt *DO);
struct node* buildFirstTree(struct DatAndOpt *DandO);




