#include <math.h>
#include "mex.h"
#define EPSILON 2.220446e-16
/**
 * CBCD size 1 with general constraints
 **/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /** input args
     * [1] Matrix A
     * [2] Vettor b
     * [3] dimension d
     * [4] max iters
     * [5] precision of KKT conditions
     * here the lower and upper bounds are 
     * set by the input parameters
     * if the number of inputs is 5, then it has default constraint [0,1]
     * if the number of inputs is 7, then the last two are constraints
     * if lower is greater than upper, then there will be no constraints
     * [6] lower bound l
     * [7] upper bound u
     * [8] initial value of out_x
     */
    // first check the number of input values
    if((nrhs!=5)&&(nrhs!=8)){mexErrMsgTxt("Input Value Error");  return;}
    // [1]
    double *in_A;
    mwIndex *irs;// for sparse matrix
    mwIndex *jcs;// for sparse matrix
    // [2]
    double *in_b;
    // [3]
    int in_d;
    // [4]
    int in_max_iter;
    // [5]
    double in_precision;
    // [6]
    double lower;
    // [7]
    double upper;
    // [8]
    double in_init;
    /** ouput args
     * [1] the minimizer out_x
     * [2] the KKT condition residual
     *     the length will be adjusted in the end,
     *     and final output is *r
     */
    // [1]
    double *out_x;
    // [2]
    double residual;
    // parameters in the function
    int i,j,epoch;//loop
    double df;
    // get input args
    in_A = mxGetPr(prhs[0]);if(in_A==NULL){mexErrMsgTxt("pointer in_A is null");  return;}
    irs = mxGetIr(prhs[0]);if(irs==NULL){mexErrMsgTxt("pointer irs is null");  return;}
    jcs = mxGetJc(prhs[0]);if(jcs==NULL){mexErrMsgTxt("pointer jcs is null");  return;}
    in_b = mxGetPr(prhs[1]);if(in_b==NULL){mexErrMsgTxt("pointer in_b is null");  return;}
    in_d = mxGetScalar(prhs[2]);if(in_d==NULL){mexErrMsgTxt("pointer in_d is null");  return;}
    in_max_iter = mxGetScalar(prhs[3]);if(in_max_iter==NULL){mexErrMsgTxt("max_iter can not be 0");  return;}
    in_precision = mxGetScalar(prhs[4]);if(in_max_iter==NULL){mexErrMsgTxt("precisionr can not be 0");  return;}
    // set the bounds according to the input
    if (nrhs==8){
        lower = mxGetScalar(prhs[5]);//if(lower==NULL){mexErrMsgTxt("pointer lower is null");  return;};
        upper = mxGetScalar(prhs[6]);//if(upper==NULL){mexErrMsgTxt("pointer upper is null");  return;};
        in_init = mxGetScalar(prhs[7]);
        /* make sure that the upper bound is larger than the lower bound
         * if the lower bound is greater, 
         * then we solve for an unconstrained problem, with init 0
         */
        if (lower>=upper){mexPrintf("Bounds Error, Results without constraints\n");in_init=0;}
        /* if the in_init is out of the bound, set it as lower */
        else if ((in_init<lower)||(in_init>upper)){in_init=lower;}
    }
    else if (nrhs==5){
        // if the constraints are not set,
        // the default bound is [0,1]
        // the default initial is 0
        lower = 0;
        upper = 1;
        in_init = 0;
    }
    /* Non-Zero elements, is the value of last entry of jcs
     * lengths of in_A and irs are both NZmax
     * length of jcs is in_d + 1, and the last entry of jcs has value NZmax
    */
    int NZmax = jcs[in_d];
    mexPrintf("CBCD size 1.cpp...Sparsity = %.5f.\n",NZmax/double((in_d*in_d)));
    // allocate output, and init as all 0s
    plhs[0] = mxCreateDoubleMatrix(in_d,1,mxREAL);
    out_x = mxGetPr(plhs[0]);if(out_x==NULL){mexErrMsgTxt("pointer out_x is null");  return;} 
    for (i=0;i<in_d;i++){
        out_x[i] = 0;
    }
    // pre-allocate output of residual, length as max_iter
    double* out_r=new double[in_max_iter]; if(out_r==NULL){mexErrMsgTxt("pointer out_r is null");  return;} 
    // allocate gradient, will delete later
    double* grad=new double[in_d];  if(grad==NULL){mexErrMsgTxt("pointer grad is null");  return;} 
    // allocate diagonal, for solving small block problem
    double* diag_A=new double[in_d];  if(diag_A==NULL){mexErrMsgTxt("pointer diag_A is null");  return;} 
    /*grad and residual of init loop
     *out_x is initialized as in_init
     *in the loop the elements from the diagonal are also extracted
    */
    // init gradient as 0s
    // extract i th element in diagonal and below
    for (i=0;i<in_d;i++){
        grad[i]=0;
        // i th diag_A[]
        diag_A[i]=0;
        for (j=jcs[i];j<jcs[i+1];j++){
            if (irs[j]==i){
                diag_A[i]=in_A[j];
                //mexPrintf("diag_A[%d]=%.5f;\n",i,in_A[j]);
            }
        }
    }
    // update gradient
    for (j=0;j<in_d;j++){
        for (i=jcs[j];i<jcs[j+1];i++){
            grad[irs[i]] += in_A[i]*in_init;
        }
    }
    // get the residual of KKT condition
    residual = 0;
    // if with constraints
    if (lower<upper){
        for (i=0;i<in_d;i++){
            // i th residual
            df = grad[i]-in_b[i];
            if (in_init<=lower+2*EPSILON){
                if (df<0){
                    residual += df*df;
                }
            }
            else if (in_init>=upper-2*EPSILON){    
                if (df>0){
                    residual += df*df;
                }
            }
            else {
            residual += df*df;
            }
        }
    }
    // if without constraints
    else{
        for (i=0;i<in_d;i++){
            // i th residual
            df = grad[i]-in_b[i];
            residual += df*df;
        }
    }          
    residual = sqrt(residual);
    out_r[0] = residual;
    mexPrintf("init:     0, residual=%.15f\n",residual);
    epoch=1;
    /* if the bounds are defined, and lower<upper
     * we take them as [lower,upper]
     * then do the following
     */
    if (lower<upper){
        while ((residual>in_precision)&&(epoch<in_max_iter)){
            for (i=0;i<in_d;i++){
                /*calc temporal grad
                 *Since A is symmetric, A(i,:)*x(i) can be 
                 *replaced by A(:,i)*x(i)
                */
                for (j=jcs[i];j<jcs[i+1];j++){
                    grad[irs[j]] -= in_A[j]*out_x[i];
                }
                // descent
                out_x[i] = (in_b[i]-grad[i])/diag_A[i];
                // bounds
                if (out_x[i]>upper){
                    out_x[i] = upper;
                }
                else if (out_x[i]<lower){
                    out_x[i] = lower;
                }
                // update temporal grad
                for (j=jcs[i];j<jcs[i+1];j++){
                    grad[irs[j]] += in_A[j]*out_x[i];
                }
            }
            // init gradient as 0s
            for (i=0;i<in_d;i++){
                grad[i]=0;
            }
            // update true gradient
            for (j=0;j<in_d;j++){
                for (i=jcs[j];i<jcs[j+1];i++){
                    grad[irs[i]] += in_A[i]*out_x[j];
                }
            }
            // get the residual
            residual = 0;
            for (i=0;i<in_d;i++){
                // i th residual
                df = grad[i]-in_b[i];
                if (out_x[i]<=lower+2*EPSILON){
                    if (df<0){
                        residual += df*df;
                    }
                }
                else if (out_x[i]>=upper-2*EPSILON){    
                    if (df>0){
                        residual += df*df;
                    }
                }
                else {
                residual += df*df;
                }
            }
            residual = sqrt(residual);
            out_r[epoch] = residual;
            //mexPrintf("epoch:%5d, residual=%.15f\n",epoch,residual);
            epoch++;
        }
        
    }
    /* if the bounds are defined, but lower>=upper
     * we take them as unconstrained
     * then do the following
     */
    else if (lower>=upper){
        while ((residual>in_precision)&&(epoch<in_max_iter)){
            for (i=0;i<in_d;i++){
                /*calc temporal grad
                 *Since A is symmetric, A(i,:)*x(i) can be 
                 *replaced by A(:,i)*x(i)
                */
                for (j=jcs[i];j<jcs[i+1];j++){
                    grad[irs[j]] -= in_A[j]*out_x[i];
                }
                // descent
                // without bounds
                out_x[i] = (in_b[i]-grad[i])/diag_A[i];
                // update temporal grad
                for (j=jcs[i];j<jcs[i+1];j++){
                    grad[irs[j]] += in_A[j]*out_x[i];
                }
            }
            // init gradient as 0s
            for (i=0;i<in_d;i++){
                grad[i]=0;
            }
            // update true gradient
            for (j=0;j<in_d;j++){
                for (i=jcs[j];i<jcs[j+1];i++){
                    grad[irs[i]] += in_A[i]*out_x[j];
                }
            }
            // get the residual
            residual = 0;
            for (i=0;i<in_d;i++){
                // i th residual
                df = grad[i]-in_b[i];
                residual += df*df;
            }
            residual = sqrt(residual);
            out_r[epoch] = residual;
            //mexPrintf("epoch:%5d, residual=%.15f\n",epoch,residual);
            epoch++;
        }
    }
    plhs[1] = mxCreateDoubleMatrix(epoch,1,mxREAL);
    double* r = mxGetPr(plhs[1]);if(r==NULL){mexErrMsgTxt("pointer r is null");  return;}
    for (i=0;i<epoch;i++){
        r[i]=out_r[i];
    }
    delete grad; delete out_r; delete diag_A;
    mexPrintf("epoch:%5d, residual=%.15f\nEnd of CBCD size 1.cpp\n",epoch-1,residual);
}