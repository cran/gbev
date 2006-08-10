

/*
*  "DatAndOpt" structure contains all data and options for growing trees, random-forest, and boosted-trees 
*
*
*/

struct DatAndOpt{

/* DATA */

double *y;             /* response   */
double *yRes;          /* used to store "y-fit.y" for boosted tree growing, and simply "y" for random forrest */  
double **w;            /* covariates */
double *yPred;         /* model prediction, is p(y[i]=1|w[i][]) for 0-1 response    */
double *yPredLink;     /* model prediction on fitted scale, is logit for 0-1 response    */
double *minW;          /* minimum of w */
double *maxW;          /* maximum of w */
int n;                 /* number obs */
int p;                 /* number covar */
int *cat;              /* cat[i] is number of categories if *x[i] is categorical, 0 otherwise */
int *oobIndicator;     /* used for stochastic boosting, oobIndicator[i]=1 if i-th obs not used in fitting, 0 otherwise */
int sizeBag;           /* size of "Bag" used in fitting */
double *oobError;      /* oobError[i] contains the out-of-bag error rate at i-th boosting iteration */
double *weights;       /* weights used in estimation */
double **condMean;     /* conditional means of error free covariates (X) given "w" */
double *condSD;        /* standard dev. of X|W  */
int numNodes;          /* number of nodes */

int NumberME;          /* number covariates measured with error, tis assumed that these are placed */
                       /* at beginning of w, i.e w[i][0:(NumberME-1)] */
double *wcWeights;     /* used to store weights for working-covariates   */


/* CONSTRAINT INFO */
int numberConstraintPoints;   /* exactly that */
double* fitConstraintPoints;  /* fitted function values at "constraintPoints" */
double** constraintPoints;    /* points where 0-1 constraints are imposed, */
                              /* if "splitFunction" (see below) corresponds such a method */
int useConstraints;           /* is 1 if splitFunction method uses constraints, 0 otherwise */
int *nodeConstPoints;         /* nodeConstPoints[i] gives node-number of "constraintPoints[][i]" */
                              /* for current tree-build */


/* RANDOM SEED */
long idum;

/* MEASUREMENT ERROR INFO */
int *meIndicator;       /* 0 if error-free, 1 if measured with error */
double *stdME;          /* vector of length "p", giving standard-deviation of measurement error */ 
double *stdX;           /* vector of length "p", giving standard-deviation of "latent" covariates */





/* MEASUREMENT ERROR INFO, MIXTURE MODEL */
int NumberComponents;     /* number of components in mixture density, i.e model for density of X, error free covariate */ 
double **ratio_pw;        /* ration_pw[i][j]=p_j(w_i)*p[j]/sum(p_k(w_i)*p[k]); where p_j(w_i) is j-th component density and p[j] corresponding mixture probability*/
double **condMeanComp;    /* condMeanComp[i*NumberComp+j][] contains conditional mean vector */
                          /*                                given component "j" for individual "i"*/
double **cholCov;         /* cholesky factorization of conditional density var-covar matrix, all components, i.e this needs to  */
                          /* stored s.a.: cholCov[j][0..NumberME-1] is first row of Cj */
                          /*              cholCov[j][NumberME..2*NumberME-1] is second row of Cj */
                          /* and Cj is such that the var-covar of j-th component density is Cj%*%t(Cj) */


/****************************************************************************************************************/
/* Note: future implementation should allow mixture component modelling of measurement error as well            */
/*       also, this density could vary across individuals, say variance depends on observed covariates          */
/*       This can all be implemented by only changing the computation of the node probabilities                 */
/*       that is, making an addition to "ComputeNodeProb()" function (i.e another option in "switch"-statement) */
/****************************************************************************************************************/


/* STORAGE FOR COMPUTING NODE-BELONGING PROBABILITY */
/* (ie data used by "ComputeNodeProb" to compute probability of true covariate vector falling in node */

/* For Monte-Carlo computation of "NODE-BELONGING PROBABILITY" */
int mc;                 /* number of mc-samples per observation */
double **x_mc;          /* x_mc[(1+(i-1)*mc):(i*mc),] monte-carlo data corresponding i-th observation */
                        /* this matrix must be filled prior to tree-building (if MC-calc is employed) */
int *node_mc;           /* node_mc[j] is node number which "x_mc[j][]" belongs to */ 
double *pred_mc;        /* pred_mc[j] is p(y=1| x_mc[j][],z[i]), used in computing "working covariate" */ 
                        /* for logitBoost */
double *predLink_mc;    /* same as pred_mc but on "link" scale */

/* TREE INFO*/
struct node *rootNode;   /* pointer to the first node in a tree */
    
/* First tree in treeList, used for procedures with multiple trees, boosting and RF*/
struct treeList *firstTree; 
int numberOfTrees;


/* OPTIONS */

double lambda;       /* regularization parameter */
int m;               /* number sampled covariates among which to find split, random forest notation */
double minSplit;        /* minimum number of obs in node to attempt a split */
double minBucket;       /* minimum number of obs in a terminal node */
int maxDepth;           /* maximum depth of tree */
int sPoints;         /* number of sampled points to attempt split */
int nodeProbMethod;  /* determines which method for computing node-belonging probabilities is used in */
                     /* function "computeNodeProb()"*/
int SplitFunction;   /* selects split function to use, like LS, LS_WC, and in future some others like binomLogLike */
int gradient_type;   /* If "gradient_type==1" then use sequentially updated gradient */
                     /* If "gradient_type==2" then use Monte-Carlo updated gradient */
                     /* If "gradient_type==3" then use Bias-adjusted gradient, if applicable */
int maxSplitAttempts; /* maximum number of attempts at finding a valid split point of node, for a given covariate*/
/* DEBUG INFO */
int track;

/* Additional Algo control parameters */
int firstTreeBuild;


/* Storage space for tree building */

struct split *splitRes;    /* Stores the results of an attempted split */
struct split *splitFinal;  /* Stores the results of a final split */
int *res;                 /* Stores a vector of sampled covariate indexes for determining split */
double *probLeft;          /* used for storing probability of belonging to left side of split node */
double *probRight;         /* used for storing probability of belonging to right side of split node */


struct node **nodeVEC; 

/* Split values */
int successful_split;    /* Stores the results of an attempted split */

/* Storage space for evaluation of independent "yPred" */
struct split *splitPred; /* Stores the results of a final split */
double **x_mc_pred;      /* x_mc_pred[(1+(i-1)*mc):(i*mc),] monte-carlo data corresponding i-th observation */
                         /* this matrix must to update yPred fit using sample independent of the one in fitting */
int *node_mc_pred;       /* node_mc_pred[j] is node number which "x_mc_pred[j][]" belongs to */ 
int Indep_y_Update;      /* if Indep_y_Update==1 then independent sample used to update yPred */
                         /* if Indep_y_Update==0 then same sample as used in fitting */



};



/*
*  "node" structure describes nodes of trees 
*
*
*/

struct node{

int depth;
struct node *parent;
struct node *leftChild;
struct node *rightChild;

int nodeNum;                 /* just that, a unique integer identifying the node */
int parentNodeNum;           /* just that, a unique integer identifying the parent node number */


int splitVar;             /* variable used to split node */
double splitPoint;        
double changeRSS;         /* change in RSS due to split*/
double priorRSS;          /* RSS prior to split */
double postRSS;           /* RSS after split */

double *probInNode;       /* vector giving probabilities of being in current node for all observations */
double theta;             /* prediction at node */

double numInNode;         /* expected number of observations falling in node, tis: sum(probInNode) */

double *upperLim;         /* vector of length "p", giving upper bound on rectangle defining node */
double *lowerLim;         /* vector of length "p", giving lower bound on rectangle defining node */

};



/*
*  "forest" used for forming a linked list of root nodes of all elements in a forest
*
*
*/

struct treeList{

struct treeList *nextTree;
struct treeList *prevTree;

/* Top node of current tree */
struct node *topNode;
int numberOfNodes;

};



/*
*  "split" stores info about a split 
*
*
*/

struct split{

double changeRSS;
int splitVar;
double splitPoint;

double *pL;     /* Probability of obs in Left node, alternatively used to store "working covariates" */
double *pR;     /* Probability of obs in Right node, alternatively used to store "working covariates" */

double thetaL;     /* prediction in Left node */
double thetaR;     /* prediction  in Right node */

double numL;      /* expected number of observations in Left node */
double numR;      /* expected number of observations in Right node */


};



/*
*  "unSplitList" elements of a list for storing nodes to be split in
*                  tree-building process.
*
*/

struct unSplitList{

struct node *pnode;

struct unSplitList *nextUnsplit;

};

