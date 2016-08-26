#include "mex.h"
#include "WELL512a.c"
double fval_mex(double *A,double *b,int d,double *x);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   //input args
   double *in_A;
   double *in_b;
   int in_d;
   double *in_lower;
   double *in_upper;
   int in_max_iter;
   //ouput args
   double *out_x;// the minimizer
   double *out_y;// function value of each iter
   //parameters in the function
   int i,j,k;
   double Aix,residual,temp;
   //get input args
   in_A = mxGetPr(prhs[0]);
   in_b = mxGetPr(prhs[1]);
   in_d = mxGetScalar(prhs[2]);
   in_lower = mxGetPr(prhs[3]);
   in_upper = mxGetPr(prhs[4]);
   in_max_iter = mxGetScalar(prhs[5]);
   mexPrintf("input args get\n");
   //allocate output, and init as all 0s
   plhs[0] = mxCreateDoubleMatrix(in_d,1,mxREAL);
   out_x = mxGetPr(plhs[0]);
   for (i=0;i<in_d;i++){
       out_x[i] = 0;
   }
   plhs[1] = mxCreateDoubleMatrix(in_max_iter,1,mxREAL);
   out_y = mxGetPr(plhs[1]);
   for (i=0;i<in_max_iter;i++){
       out_y[i] = 0;
   }
   for (k=0;k<in_max_iter;k++){
       i = k%in_d;
       temp = WELLRNG512a();
       //calc vector A(i,:)*x
       Aix = 0;
       for (j=0;j<in_d;j++){
           Aix = Aix + in_A[i*in_d + j]*out_x[j];
           //mexPrintf("iter:%5d, i=%5d, j=%5d, Aix=%.8f\n",k+1,i,j,Aix);
       }
       //descent
       out_x[i] = out_x[i] - (Aix-in_b[i])/in_A[i*in_d + i];
       //bounds
       if (out_x[i]>in_upper[i]){
           out_x[i] = in_upper[i];
       }
       if (out_x[i]<in_lower[i]){
           out_x[i] = in_lower[i];
       }
       //residual
       residual=0;
       for (i=0;i<in_d;i++){
           Aix = 0;
           for (j=0;j<in_d;j++){
               Aix = Aix + in_A[i*in_d + j]*out_x[j];
           }
           residual = residual + pow(Aix-in_b[i],2);
       }
       residual = sqrt(residual);
       out_y[k] = fval_mex(in_A, in_b, in_d, out_x);
       mexPrintf("iter:%5d, residual=%.8f, fval:%.8f\n",k+1,temp,out_y[k]);
   }
}
double fval_mex(double *A,double *b,int d,double *x){
    double y;//output function value y in R^(iter*1)
    int i,j,k;
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