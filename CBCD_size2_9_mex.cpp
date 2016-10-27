#include <math.h>
#include "mex.h"
#include <time.h>

#define EPSILON 2.220446e-16
double fval_mex(double *A,double *b,int d,double *x);
double residual_mex(double *A,double *b,int d,double *x,double *grad);
double* grad_mex(double *A,double *x,int d, double *grad);
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
   double residual;
   //get input args
   in_A = mxGetPr(prhs[0]);if(in_A==NULL){mexErrMsgTxt("pointer in_A is null");  return;}
   in_b = mxGetPr(prhs[1]);if(in_b==NULL){mexErrMsgTxt("pointer in_b is null");  return;}
   in_d = mxGetScalar(prhs[2]);if(in_d==NULL){mexErrMsgTxt("pointer in_d is null");  return;}
   in_max_iter = mxGetScalar(prhs[3]);if(in_max_iter==NULL){mexErrMsgTxt("pointer in_max_iter is null");  return;}
   mexPrintf("CBCD size 2...input args get\n");
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
   //grad and residual of init
   grad = grad_mex(in_A,out_x,in_d,grad);
   residual = residual_mex(in_A, in_b, in_d, out_x, grad);
   out_r[0] = residual;
   mexPrintf("init:     0, residual=%.15f\n",residual);
   epoch=1;
   // parameters for size 2 matrix
   double a11,a12,a21,a22;
   double b1,b2;
   bool FLAG;
   double a21b1_a11,a12b2_a22,detA2;
   clock_t time0,time1,time2,time3,time4;////////time//////////
   double dt1=0; double dt2=0; double dt3=0; double dt4=0;////////time//////////
   while ((residual>1E-13)&&(epoch<=in_max_iter)){
       //dt1=0;dt2=0;dt3=0;dt4=0;////////time//////////
       for (i=0;i<in_d-1;i=i+2){
           //calc temporal grad
           time0 = clock();////////time//////////
           for (j=0;j<in_d;j++){
               grad[j] = grad[j]-in_A[i*in_d+j]*out_x[i]-in_A[(i+1)*in_d+j]*out_x[i+1];
           }
           time1 = clock();////////time//////////
           // update x(i)
           // define size 2 block
           a11=in_A[i*in_d+i  ];     a12=in_A[(i+1)*in_d+i  ];
           a21=in_A[i*in_d+i+1];     a22=in_A[(i+1)*in_d+i+1];
           b1 =in_b[i]  -grad[i];
           b2 =in_b[i+1]-grad[i+1];
           // decission tree
           FLAG = false;
           a21b1_a11 = a21*b1/a11;
           a12b2_a22 = a12*b2/a22;
           detA2 = a11*a22-a12*a21;
           // first assume x2=0
           if (b1<=0 && b2<=0){//case 1
                   out_x[i]=0;
                   out_x[i+1]=0;
                   FLAG = true;//mexPrintf("case1\n");
           }
           else if (b1>0 && b1<a11 && b2<=a21b1_a11){//case 2
                   out_x[i]=b1/a11;
                   out_x[i+1]=0;
                   FLAG = true;//mexPrintf("case2\n");
           }
           else if(b1>=a11 && b2<=a21){//case 3
                   out_x[i]=1;
                   out_x[i+1]=0;
                   FLAG = true;//mexPrintf("case3\n");
           }
           // x2~=0, assume x2=1
           else if (b1<=a12b2_a22 && b2>0 && b2<a22){//case 4
                   out_x[i]=0;
                   out_x[i+1]=b2/a22;
                   FLAG = true;//mexPrintf("case4\n");
           }
           else if(b1>=a12b2_a22 + detA2/a22 && b2>a21 && b2<a21+a22){//case 6
                   out_x[i]=1;
                   out_x[i+1]=(b2-a21)/a22;
                   FLAG = true;//mexPrintf("case6\n");
           }
           // x2~=0 & x2~=1 x2 in (0,1)
           else if (b1<=a12 && b2>=a22){//case 7
                   out_x[i]=0;
                   out_x[i+1]=1;
                   FLAG = true;//mexPrintf("case7\n");
           }
           else if (b1>a12 && b1<a12+a11 && b2>=a21b1_a11+detA2/a11){//case 8
                   out_x[i]=(b1-a12)/a11;
                   out_x[i+1]=1;
                   FLAG = true;//mexPrintf("case8\n");
           }
           else if(b1>=a11+a12 && b2>=a21+a22){//case 9
                   out_x[i]=1;
                   out_x[i+1]=1;
                   FLAG = true;//mexPrintf("case9\n");
           }
           else{//case 5
                   out_x[i  ]=(a22*b1-a12*b2)/detA2;
                   out_x[i+1]=(a11*b2-a21*b1)/detA2;
                   FLAG = true;
           }
           time2 = clock();////////time//////////
           //update temporal grad
           for (j=0;j<in_d;j++){
               grad[j] = grad[j]+in_A[i*in_d+j]*out_x[i]+in_A[(i+1)*in_d+j]*out_x[i+1];
           }
           time3 = clock();////////time//////////
           dt1+=(double)(time1-time0);
           dt2+=(double)(time2-time1);
           dt3+=(double)(time3-time2);///(CLOCKS_PER_SEC);
           dt4+=(double)(time3-time0);
       }
       if (in_d%2==1){
           i=in_d-1;
           //calc temporal grad
           for (j=0;j<in_d;j++){
               grad[j] = grad[j] - in_A[i*in_d+j]*out_x[i];
           }
           //descent
           out_x[i] = (in_b[i]-grad[i])/in_A[i*in_d + i];
           //bounds
           if (out_x[i]>1){
               out_x[i] = 1;
           }
           if (out_x[i]<0){
               out_x[i] = 0;
           }
           //update temporal grad
           for (j=0;j<in_d;j++){
               grad[j] = grad[j] + in_A[i*in_d+j]*out_x[i];
           }
       }
       //update true gradient
       grad = grad_mex(in_A,out_x,in_d,grad);
       //residual
       residual = residual_mex(in_A,in_b,in_d,out_x,grad);
       out_r[epoch] = residual;
       if (epoch%5==1){
           mexPrintf("epoch:%5d, residual=%.15f\n",epoch,residual);
       }
       epoch++;
   }
   plhs[1] = mxCreateDoubleMatrix(epoch,1,mxREAL);
   double* r = mxGetPr(plhs[1]);if(r==NULL){mexErrMsgTxt("pointer r is null");  return;}
   for (i=0;i<epoch;i++){
       r[i]=out_r[i];
   }
   delete grad; delete out_r;
   mexPrintf("dt1 = %.5f, dt2 = %.5f, dt3 = %.5f, dt4 = %.5f\n",dt1,dt2,dt3,dt4);
   mexPrintf("End of CBCD size 2\n");
   
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