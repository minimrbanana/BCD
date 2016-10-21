#include <math.h>
#include "mex.h"
//#include <time.h>

#define EPSILON 2.220446e-16
double fval_mex(double *A,double *b,int d,double *x);
double residual_mex(double *A,double *b,int d,double *x,double *grad);
double *grad_mex(double *A,double *x,int d, double *grad);
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
   int i,j,epoch;//loop
   double residual,xtemp;
   //get input args
   in_A = mxGetPr(prhs[0]);if(in_A==NULL){mexErrMsgTxt("pointer in_A is null");  return;}
   in_b = mxGetPr(prhs[1]);if(in_b==NULL){mexErrMsgTxt("pointer in_b is null");  return;}
   in_d = mxGetScalar(prhs[2]);if(in_d==NULL){mexErrMsgTxt("pointer in_d is null");  return;}
   in_max_iter = mxGetScalar(prhs[3]);if(in_max_iter==NULL){mexErrMsgTxt("pointer in_max_iter is null");  return;}
   mexPrintf("CBCD size 1...input args get\n");
   //allocate output, and init as all 0s
   plhs[0] = mxCreateDoubleMatrix(in_d,1,mxREAL);
   out_x = mxGetPr(plhs[0]);if(out_x==NULL){mexErrMsgTxt("pointer out_x is null");  return;} else{mexPrintf("pointer out_x created\n");} 
   for (i=0;i<in_d;i++){
       out_x[i] = 0;
   }
   //pre-allocate output of residual, length as max_iter
   double* out_r=new double[in_max_iter]; if(out_r==NULL){mexErrMsgTxt("pointer out_r is null");  return;} else{mexPrintf("pointer out_r created\n");} 
   //allocate gradient, will delete later
   double* grad=new double[in_d];  if(grad==NULL){mexErrMsgTxt("pointer grad is null");  return;} else{mexPrintf("pointer grad created\n");} 
   //grad and residual of init
   grad = grad_mex(in_A,out_x,in_d,grad);
   residual = residual_mex(in_A, in_b, in_d, out_x, grad);
   out_r[0] = residual;
   mexPrintf("init:     0, residual=%.15f\n",residual);
   epoch=1;
   while ((residual>1E-13)&&(epoch<=in_max_iter)){
       for (i=0;i<in_d;i++){
           //calc temporal grad
           for (j=0;j<in_d;j++){
               grad[j] -= in_A[i*in_d+j]*out_x[i];
           }
           //descent
           out_x[i] = (in_b[i]-grad[i])/in_A[i*in_d + i];
           //bounds
           if (out_x[i]>1){
               out_x[i] = 1;
           }
           else if (out_x[i]<0){
               out_x[i] = 0;
           }
           //update temporal grad
           for (j=0;j<in_d;j++){
               grad[j] += in_A[i*in_d+j]*out_x[i];
           }
       }
       //update true gradient
       grad = grad_mex(in_A,out_x,in_d,grad);
       //residual
       residual = residual_mex(in_A,in_b,in_d,out_x,grad);
       out_r[epoch] = residual;
       mexPrintf("epoch:%5d, residual=%.15f\n",epoch,residual);
       epoch++;
   }
   plhs[1] = mxCreateDoubleMatrix(epoch,1,mxREAL);
   double* r = mxGetPr(plhs[1]);if(r==NULL){mexErrMsgTxt("pointer r is null");  return;}
   for (i=0;i<epoch;i++){
       r[i]=out_r[i];
   }
   delete grad; delete out_r;
   mexPrintf("End of CBCD size 1\n");
}
double fval_mex(double *A,double *b,int d,double *x){
    double y;//output function value, y in R^(iter*1)
    int i,j;
    double Aix;
    y = 0;
    for (i=0;i<d;i++){
        Aix = 0;
        for (j=0;j<d;j++){
            Aix = Aix + A[i*d + j]*x[j];
        }
        y = y + (Aix/2-b[i]) * x[i];
    }
    return y;
}
double residual_mex(double *A,double *b,int d,double *x,double *grad){
    double r;//output residual, scalar
    int i,j;
    double df;
    r = 0;
    for (i=0;i<d;i++){
        df = grad[i]-b[i];
        if (x[i]<=0+2*EPSILON){
            if (df<0){
                r = r + df*df;
            }
        }
        else if (x[i]>=1-2*EPSILON){    
            if (df>0){
                r = r + df*df;
            }
        }
        else {
            r = r + df*df;
        }
    }
    r = sqrt(r);
    return r;
}
double *grad_mex(double *A,double *x,int d, double *grad){
    for (int i=0;i<d;i++){
        grad[i]=0;
        for (int j=0;j<d;j++){
            grad[i]=grad[i]+A[i*d + j]*x[j];
        }
    }
    return grad;
}