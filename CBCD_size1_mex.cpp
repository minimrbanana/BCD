#include <math.h>
#include "mex.h"
//#include <time.h>

#define EPSILON 2.220446e-16
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //input args
    //here the lower and upper bounds are set in the program
    //not set by the input parameters
    double *in_A;
    double *in_b;
    int in_d;
    int in_max_iter;
    //ouput args
    double *out_x;// the minimizer
    //parameters in the function
    int i,j,epoch,iid;//loop
    double residual, df;
    //get input args
    in_A = mxGetPr(prhs[0]);if(in_A==NULL){mexErrMsgTxt("pointer in_A is null");  return;}
    in_b = mxGetPr(prhs[1]);if(in_b==NULL){mexErrMsgTxt("pointer in_b is null");  return;}
    in_d = mxGetScalar(prhs[2]);if(in_d==NULL){mexErrMsgTxt("pointer in_d is null");  return;}
    in_max_iter = mxGetScalar(prhs[3]);if(in_max_iter==NULL){mexErrMsgTxt("pointer in_max_iter is null");  return;}
    mexPrintf("CBCD size 1.cpp...input args get\n");
    //allocate output, and init as all 0s
    plhs[0] = mxCreateDoubleMatrix(in_d,1,mxREAL);
    out_x = mxGetPr(plhs[0]);if(out_x==NULL){mexErrMsgTxt("pointer out_x is null");  return;} 
    for (i=0;i<in_d;i++){
        out_x[i] = 0;
    }
    //pre-allocate output of residual, length as max_iter
    double* out_r=new double[in_max_iter]; if(out_r==NULL){mexErrMsgTxt("pointer out_r is null");  return;} 
    //allocate gradient, will delete later
    double* grad=new double[in_d];  if(grad==NULL){mexErrMsgTxt("pointer grad is null");  return;} 
    //grad and residual of init in one loop
    //out_x is initialized as all 0s, so grad is 0 vector
    //residual calculation is simplified
    residual = 0;
    for (i=0;i<in_d;i++){
        //for initial x=0, g[i]=0
        grad[i]=0;
        ////iid=i*in_d;
        // i th residual
        df = -in_b[i];
        if (df<0){
            residual += df*df;
        }
    }
    residual = sqrt(residual);
    out_r[0] = residual;
    mexPrintf("init:     0, residual=%.15f\n",residual);
    epoch=1;
    while ((residual>1E-13)&&(epoch<in_max_iter)){
        for (i=0;i<in_d;i++){
            //calc temporal grad
            //use iid as temporal i*in_d
            iid = i*in_d;
            for (j=0;j<in_d;j++){
                grad[j] -= in_A[iid+j]*out_x[i];
            }
            //descent
            out_x[i] = (in_b[i]-grad[i])/in_A[iid + i];
            //bounds
            if (out_x[i]>1){
                out_x[i] = 1;
            }
            else if (out_x[i]<0){
                out_x[i] = 0;
            }
            //update temporal grad
            for (j=0;j<in_d;j++){
                grad[j] += in_A[iid+j]*out_x[i];
            }
        }
        //update true gradient and residual in one loop
        residual = 0;
        for (i=0;i<in_d;i++){
            grad[i]=0;
            iid=i*in_d;
            // i th gradient
            for (int j=0;j<in_d;j++){
                grad[i] += in_A[iid + j]*out_x[j];
            }
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
        }
        residual = sqrt(residual);
        out_r[epoch] = residual;
        //mexPrintf("epoch:%5d, residual=%.15f\n",epoch,residual);
        epoch++;
    }
    plhs[1] = mxCreateDoubleMatrix(epoch,1,mxREAL);
    double* r = mxGetPr(plhs[1]);if(r==NULL){mexErrMsgTxt("pointer r is null");  return;}
    for (i=0;i<epoch;i++){
        r[i]=out_r[i];
    }
    delete grad; delete out_r;
    mexPrintf("epoch:%5d, residual=%.15f\nEnd of CBCD size 1.cpp\n",epoch-1,residual);
}