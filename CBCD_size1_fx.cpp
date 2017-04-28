#include <math.h>
#include "mex.h"
/**
 * CBCD size 1 with function value output
 **/
#define EPSILON 2.220446e-16
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //input args
    //here the lower and upper bounds are set in the program
    //not set by the input parameters
    double *in_A;
    mwIndex *irs;// for sparse matrix
    mwIndex *jcs;// for sparse matrix
    double *in_b;
    int in_d;
    int in_max_iter;
    //ouput args
    double *out_x;// the minimizer
    //parameters in the function
    int i,j,epoch;//loop
    double residual, df;
    //get input args
    in_A = mxGetPr(prhs[0]);if(in_A==NULL){mexErrMsgTxt("pointer in_A is null");  return;}
    irs = mxGetIr(prhs[0]);if(irs==NULL){mexErrMsgTxt("pointer irs is null");  return;}
    jcs = mxGetJc(prhs[0]);if(jcs==NULL){mexErrMsgTxt("pointer jcs is null");  return;}
    in_b = mxGetPr(prhs[1]);if(in_b==NULL){mexErrMsgTxt("pointer in_b is null");  return;}
    in_d = mxGetScalar(prhs[2]);if(in_d==NULL){mexErrMsgTxt("pointer in_d is null");  return;}
    in_max_iter = mxGetScalar(prhs[3]);if(in_max_iter==NULL){mexErrMsgTxt("pointer in_max_iter is null");  return;}
    /* Non-Zero elements, is the value of last entry of jcs
     * lengths of in_A and irs are both NZmax
     * length of jcs is in_d + 1, and the last entry of jcs has value NZmax
    */
    int NZmax = jcs[in_d];
    mexPrintf("CBCD size 1.cpp...Sparsity = %.5f.\n",NZmax/double((in_d*in_d)));
    //allocate output, and init as all 0s
    plhs[0] = mxCreateDoubleMatrix(in_d,1,mxREAL);
    out_x = mxGetPr(plhs[0]);if(out_x==NULL){mexErrMsgTxt("pointer out_x is null");  return;} 
    for (i=0;i<in_d;i++){
        out_x[i] = 0;
    }
    //pre-allocate output of function value, length as max_iter
    double* out_fx=new double[in_max_iter]; if(out_fx==NULL){mexErrMsgTxt("pointer out_fx is null");  return;} 
    //allocate gradient, will delete later
    double* grad=new double[in_d];  if(grad==NULL){mexErrMsgTxt("pointer grad is null");  return;} 
    //allocate diagonal, for solving small block problem
    double* diag_A=new double[in_d];  if(diag_A==NULL){mexErrMsgTxt("pointer diag_A is null");  return;} 
    /*grad and residual of init in one loop
     *out_x is initialized as all 0s, so grad is 0 vector
     *residual calculation is simplified
     *in the loop the elements from the diagonal are also extracted
    */
    residual = 0;
    for (i=0;i<in_d;i++){
        //for initial x=0, g[i]=0
        grad[i]=0;
        //i th residual
        df = -in_b[i];
        if (df<0){
            residual += df*df;
        }
        //i th diag_A[]
        diag_A[i]=0;
        for (j=jcs[i];j<jcs[i+1];j++){
            if (irs[j]==i){
                diag_A[i]=in_A[j];
                //mexPrintf("diag_A[%d]=%.5f;\n",i,in_A[j]);
            }
        }
    }
    residual = sqrt(residual);
    out_fx[0] = 0;
    mexPrintf("init:     0, residual=%.15f\n",residual);
    epoch=1;
    while ((residual>1E-13)&&(epoch<in_max_iter)){
        for (i=0;i<in_d;i++){
            /*calc temporal grad
             *Since A is symmetric, A(i,:)*x(i) can be 
             *replaced by A(:,i)*x(i)
            */
            for (j=jcs[i];j<jcs[i+1];j++){
                grad[irs[j]] -= in_A[j]*out_x[i];
            }
            //descent
            out_x[i] = (in_b[i]-grad[i])/diag_A[i];
            //bounds
            if (out_x[i]>1){
                out_x[i] = 1;
            }
            else if (out_x[i]<0){
                out_x[i] = 0;
            }
            //update temporal grad
            for (j=jcs[i];j<jcs[i+1];j++){
                grad[irs[j]] += in_A[j]*out_x[i];
            }
        }
        //init gradient as 0s
        for (i=0;i<in_d;i++){
            grad[i]=0;
        }
        //update true gradient
        for (j=0;j<in_d;j++){
            for (i=jcs[j];i<jcs[j+1];i++){
                grad[irs[i]] += in_A[i]*out_x[j];
            }
        }
        //get the residual and function value
        residual = 0;
        out_fx[epoch] = 0;
        for (i=0;i<in_d;i++){
            // i th residual
            df = grad[i]-in_b[i];
            if (out_x[i]<=0+2*EPSILON){
                if (df<0){
                    residual += df*df;
                }
            }
            else if (out_x[i]>=1-2*EPSILON){    
                if (df>0){
                    residual += df*df;
                }
            }
            else {
            residual += df*df;
            }
            // function value
            // fx = 0.5*<x,grad>-<b,x>
            out_fx[epoch] += (df-in_b[i])/2.0*out_x[i];
        }
        residual = sqrt(residual);
        //mexPrintf("epoch:%5d, residual=%.15f\n",epoch,residual);
        epoch++;
    }
    //allocate output pointer for f(x)
    plhs[1] = mxCreateDoubleMatrix(epoch,1,mxREAL);
    double* fx = mxGetPr(plhs[1]);if(fx==NULL){mexErrMsgTxt("pointer fx is null");  return;}
    for (i=0;i<epoch;i++){
        fx[i]=out_fx[i];
    }
    delete grad; delete out_fx; delete diag_A;
    mexPrintf("epoch:%5d, residual=%.15f\nEnd of CBCD size 1.cpp\n",epoch-1,residual);
}