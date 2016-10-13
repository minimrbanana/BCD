#include <math.h>
#include "mex.h"

#define EPSILON 2.220446e-16
double fval_mex(double *A,double *b,int d,double *x);
double residual_mex(double *A,double *b,int d,double *x,double *grad);
double *grad_mex(double *A,double *x,int d, double *grad);
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
   //parameters in the function
   int i,j,epoch;//loop
   double residual;
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
   double *out_y = NULL;// function value of each iter
   out_y = (double*)malloc(sizeof(double)*in_d);
   out_y[0] = fval_mex(in_A, in_b, in_d, out_x);
   for (i=1;i<=in_max_iter;i++){
       out_y[i] = 0;
   }
   // residual of init
   double *grad=NULL;
   grad = (double*)malloc(sizeof(double)*in_d);
   grad = grad_mex(in_A,out_x,in_d,grad);
   residual = 1;
   mexPrintf("epoch:    0, residual=%.15f, fval:%.8f\n",residual,out_y[0]);
   epoch=1;
   double A2[2][2];
   double b2[2];
   bool FLAG;
   double a21b1_a11,a12b2_a22,detA2;
   while ((residual>1E-13)&&(epoch<=in_max_iter)){
       for (i=0;i<in_d-1;i=i+2){
           //calc temporal grad
           for (j=0;j<in_d;j++){
               grad[j] = grad[j]-in_A[j*in_d+i]*out_x[i]-in_A[j*in_d+i+1]*out_x[i+1];
           }
           // update x(i)
           // define size 2 block
           A2[0][0]=in_A[i*in_d+i];     A2[0][1]=in_A[i*in_d+i+1];
           A2[1][0]=in_A[(i+1)*in_d+i]; A2[1][1]=in_A[(i+1)*in_d+i+1];
           b2[0]   =in_b[i]-grad[i];
           b2[1]   =in_b[i+1]-grad[i+1];
           // decission tree
           FLAG = false;
           a21b1_a11 = A2[1][0]*b2[0]/A2[0][0];
           a12b2_a22 = A2[0][1]*b2[1]/A2[1][1];
           detA2 = A2[0][0]*A2[1][1]-A2[0][1]*A2[1][0];
           // first assume x2=0
           if (b2[0]<=0 && b2[1]<=0){//case 1
                   out_x[i]=0;
                   out_x[i+1]=0;
                   FLAG = true;//mexPrintf("case1\n");
           }
           else if (b2[0]>0 && b2[0]<A2[0][0] && b2[1]<=a21b1_a11){//case 2
                   out_x[i]=b2[0]/A2[0][0];
                   out_x[i+1]=0;
                   FLAG = true;//mexPrintf("case2\n");
           }
           else if(b2[0]>=A2[0][0] && b2[1]<=A2[1][0]){//case 3
                   out_x[i]=1;
                   out_x[i+1]=0;
                   FLAG = true;//mexPrintf("case3\n");
           }
           // x2~=0, assume x2=1
           else if (b2[0]<=a12b2_a22 && b2[1]>0 && b2[1]<A2[1][1]){//case 4
                   out_x[i]=0;
                   out_x[i+1]=b2[1]/A2[1][1];
                   FLAG = true;//mexPrintf("case4\n");
           }
           else if(b2[0]>=a12b2_a22 + detA2/A2[1][1] && b2[1]>A2[1][0] && b2[1]<A2[1][0]+A2[1][1]){//case 6
                   out_x[i]=1;
                   out_x[i+1]=(b2[1]-A2[1][0])/A2[1][1];
                   FLAG = true;//mexPrintf("case6\n");
           }
           // x2~=0 & x2~=1 x2 in (0,1)
           else if (b2[0]<=A2[0][1] && b2[1]>=A2[1][1]){//case 7
                   out_x[i]=0;
                   out_x[i+1]=1;
                   FLAG = true;//mexPrintf("case7\n");
           }
           else if (b2[0]>A2[0][1] && b2[0]<A2[0][1]+A2[0][0] && b2[1]>=a21b1_a11+detA2/A2[0][0]){//case 8
                   out_x[i]=(b2[0]-A2[0][1])/A2[0][0];
                   out_x[i+1]=1;
                   FLAG = true;//mexPrintf("case8\n");
           }
           else if(b2[0]>=A2[0][0]+A2[0][1] && b2[1]>=A2[1][0]+A2[1][1]){//case 9
                   out_x[i]=1;
                   out_x[i+1]=1;
                   FLAG = true;//mexPrintf("case9\n");
           }
           else{//case 5
                   out_x[i]=(A2[1][1]*b2[0]-A2[0][1]*b2[1])/detA2;
                   out_x[i+1]=(A2[0][0]*b2[1]-A2[1][0]*b2[0])/detA2;
                   FLAG = true;
           }
           //update temporal grad
           for (j=0;j<in_d;j++){
               grad[j] = grad[j]+in_A[j*in_d+i]*out_x[i]+in_A[j*in_d+i+1]*out_x[i+1];
           }
       }
       // if mod(in_d,2)==1
       if (in_d%2==1){
           i=in_d-1;
           //calc temporal grad
           for (j=0;j<in_d;j++){
               grad[j] = grad[j] - in_A[j*in_d+i]*out_x[i];
           }
           //descent
           out_x[i] = (in_b[i]-grad[i])/in_A[i*in_d + i];
           //bounds
           if (out_x[i]>in_upper[i]){
               out_x[i] = in_upper[i];
           }
           if (out_x[i]<in_lower[i]){
               out_x[i] = in_lower[i];
           }
           //update temporal grad
           for (j=0;j<in_d;j++){
               grad[j] = grad[j] + in_A[j*in_d+i]*out_x[i];
           }
       }
       grad = grad_mex(in_A,out_x,in_d,grad);
       //residual
       residual = residual_mex(in_A,in_b,in_d,out_x,grad);
       out_y[epoch] = fval_mex(in_A, in_b, in_d, out_x);
       mexPrintf("epoch:%5d, residual=%.15f, fval:%.15f\n",epoch,residual,out_y[epoch]);
       epoch++;
   }
   plhs[1] = mxCreateDoubleMatrix(epoch,1,mxREAL);
   double* y = mxGetPr(plhs[1]);
   for (i=0;i<epoch;i++){
       y[i]=out_y[i];
   }
   delete grad; delete out_y;
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