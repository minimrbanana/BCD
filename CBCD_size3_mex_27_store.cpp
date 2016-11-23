#include <math.h>
#include "mex.h"
#include <time.h>

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
   double residual;
   //get input args
   in_A = mxGetPr(prhs[0]);if(in_A==NULL){mexErrMsgTxt("pointer in_A is null");  return;}
   in_b = mxGetPr(prhs[1]);if(in_b==NULL){mexErrMsgTxt("pointer in_b is null");  return;}
   in_d = mxGetScalar(prhs[2]);if(in_d==NULL){mexErrMsgTxt("pointer in_d is null");  return;}
   in_max_iter = mxGetScalar(prhs[3]);if(in_max_iter==NULL){mexErrMsgTxt("pointer in_max_iter is null");  return;}
   mexPrintf("CBCD size 3...input args get\n");
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
   // parameters for size 3 matrix
   double a11,a12,a13,a21,a22,a23,a31,a32,a33;
   double b1,b2,b3,x1,x2,x3,detA;
   // FLAG stores the label of last choice of each block
   int* labels=new int[in_d];  if(labels==NULL){mexErrMsgTxt("pointer labels is null");  return;} ;
   for (i=0;i<in_d;i++){
       labels[i]=0;
   }
   int FLAG;
   clock_t time0,time1,time2,time3,time4;////////time//////////
   double dt1=0; double dt2=0; double dt3=0; double dt4=0;////////time//////////
   while ((residual>1E-13)&&(epoch<=in_max_iter)){
       //double dt1=0; double dt2=0; double dt3=0; double dt4=0;////////time//////////
       for (i=0;i<in_d-2;i=i+3){
           //calc temporal grad
           time0 = clock();////////time//////////
           for (j=0;j<in_d;j++){
               grad[j] = grad[j]-in_A[i*in_d+j]*out_x[i]-in_A[(i+1)*in_d+j]*out_x[i+1]-in_A[(i+2)*in_d+j]*out_x[i+2];
           }
           time1 = clock();////////time//////////
           // update x(i)
           // define size 3 block
           a11=in_A[i*in_d+i  ]; a12=in_A[(i+1)*in_d+i  ]; a13=in_A[(i+2)*in_d+i  ];
           a21=in_A[i*in_d+i+1]; a22=in_A[(i+1)*in_d+i+1]; a23=in_A[(i+2)*in_d+i+1];
           a31=in_A[i*in_d+i+2]; a32=in_A[(i+1)*in_d+i+2]; a33=in_A[(i+2)*in_d+i+2];
           b1 =in_b[i]  -grad[i];
           b2 =in_b[i+1]-grad[i+1];
           b3 =in_b[i+2]-grad[i+2];
           // decission tree
           FLAG = labels[i];
           switch (FLAG){
               case 1: {
                   if (b1<=0 && b2<=0 && b3<=0){
                        out_x[i]=0;
                        out_x[i+1]=0;
                        out_x[i+2]=0;
                        FLAG=1;
                   }
                   else {FLAG=0;}
               }
               case 2: {
                   if (a13>=b1 && a23>=b2 && a33<=b3){
                        out_x[i]=0;
                        out_x[i+1]=0;
                        out_x[i+2]=1;
                        FLAG=2;
                   }
                   else {FLAG=0;}
               }
               case 3: {
                   if (a12>=b1 && a22<=b2 && a32>=b3){
                        out_x[i]=0;
                        out_x[i+1]=1;
                        out_x[i+2]=0;
                        FLAG=3;
                   }
                   else {FLAG=0;}
               }
               case 4: {
                   if (a12+a13>=b1 && a22+a23<=b2 && a32+a33<=b3){
                        out_x[i]=0;
                        out_x[i+1]=1;
                        out_x[i+2]=1;
                        FLAG=4;
                   }
                   else {FLAG=0;}
               }
               case 5: {
                   if (a11<=b1 && a12>=b2 && a13>=b3){
                        out_x[i]=1;
                        out_x[i+1]=0;
                        out_x[i+2]=0;
                        FLAG=5;
                   }
                   else {FLAG=0;}
               }
               case 6: {
                   if (a11+a13<=b1 && a21+a23>=b2 && a31+a33<=b3){
                        out_x[i]=1;
                        out_x[i+1]=0;
                        out_x[i+2]=1;
                        FLAG=6;
                   }
                   else {FLAG=0;}
               }
               case 7: {
                   if (a11+a12<=b1 && a21+a22<=b2 && a31+a32>=b3){
                        out_x[i]=1;
                        out_x[i+1]=1;
                        out_x[i+2]=0;
                        FLAG=7;
                   }
                   else {FLAG=0;}
               }
               case 8: {
                   if (a11+a12+a13<=b1 && a21+a22+a23<=b2 && a31+a32+a33<=b3){
                        out_x[i]=1;
                        out_x[i+1]=1;
                        out_x[i+2]=1;
                        FLAG=8;
                   }
                   else {FLAG=0;}
               }
               case 9: {
                   x3 = b3/a33;
                   if (x3>=0 && x3<=1 && a13*x3>=b1 && a23*x3>=b2){
                        out_x[i]=0;
                        out_x[i+1]=0;
                        out_x[i+2]=x3;
                        FLAG=9;
                   }
                   else {FLAG=0;}
               }
               case 10: {
                    x2 = b2/a22;
                    if (x2>=0 && x2<=1 && a12*x2>=b1 && a32*x2>=b3){
                       out_x[i]=0;
                       out_x[i+1]=x2;
                       out_x[i+2]=0;
                       FLAG=10;
                   }
                   else {FLAG=0;}
               }
               case 11: {
                   x1 = b1/a11;
                   if (x1>=0 && x1<=1 && a21*x1>=b2 && a31*x1>=b3){
                       out_x[i]=x1;
                       out_x[i+1]=0;
                       out_x[i+2]=0;
                       FLAG=11;
                   }
                   else {FLAG=0;}
               }
               case 12: {
                   x1=(b1-a12-a13)/a11;
                   if (x1>=0 && x1<=1 && a21*x1+a22+a23<=b2 && a31*x1+a32+a33<=b3){
                       out_x[i]=x1;
                       out_x[i+1]=1;
                       out_x[i+2]=1;
                       FLAG=12;
                   }
                   else {FLAG=0;}
               }
               case 13: {
                   x2=(b2-a21-a23)/a22;
                   if (x2>=0 && x2<=1 && a11+a12*x2+a13<=b1 && a31+a32*x2+a33<=b3){
                       out_x[i]=1;
                       out_x[i+1]=x2;
                       out_x[i+2]=1;
                       FLAG=13;
                   }
                   else {FLAG=0;}
               }
               case 14: {
                   x3=(b3-a31-a32)/a33;
                   if (x3>=0 && x3<=1 && a11+a12+a13*x3<=b1 && a21+a22+a23*x3<=b2){
                       out_x[i]=1;
                       out_x[i+1]=1;
                       out_x[i+2]=x3;
                       FLAG=14;
                   }
                   else {FLAG=0;}
               }
               case 15: {
                   x1=(b1-a13)/a11;
                   if (x1>=0 && x1<=1 && a21*x1+a23>=b2 && a31*x1+a33<=b3){
                       out_x[i]=x1;
                       out_x[i+1]=0;
                       out_x[i+2]=1;
                       FLAG=15;
                   }
                   else {FLAG=0;}
               }
               case 16: {
                   x1 = (b1-a12)/a11;
                   if (x1>=0 && x1<=1 && a21*x1+a22<=b2 && a31*x1+a32>=b3){
                       out_x[i]=x1;
                       out_x[i+1]=1;
                       out_x[i+2]=0;
                       FLAG=16;
                   }
                   else {FLAG=0;}
               }
               case 17: {
                   x2 = (b2-a23)/a22;
                   if (x2>=0 && x2<=1 && a12*x2+a13>=b1 && a32*x2+a33<=b3){
                       out_x[i]=0;
                       out_x[i+1]=x2;
                       out_x[i+2]=1;
                       FLAG=17;
                   }
                   else {FLAG=0;}
               }
               case 18: {
                   x2=(b2-a21)/a22;
                   if (x2>=0 && x2<=1 && a11+a12*x2<=b1 && a31+a32*x2>=b3){
                       out_x[i]=1;
                       out_x[i+1]=x2;
                       out_x[i+2]=0;
                       FLAG=18;
                   }
                   else {FLAG=0;}
               }
               case 19: {
                   x3=(b3-a32)/a33;
                   if (x3>=0 && x3<=1 && a12+a13*x3>=b1 && a22+a23*x3<=b2){
                       out_x[i]=0;
                       out_x[i+1]=1;
                       out_x[i+2]=x3;
                       FLAG=19;
                   }
                   else {FLAG=0;}
               }
               case 20: {
                   x3=(b3-a31)/a33;
                   if (x3>=0 && x3<=1 && a11+a13*x3<=b1 && a21+a23*x3>=b2){
                       out_x[i]=1;
                       out_x[i+1]=0;
                       out_x[i+2]=x3;
                       FLAG=20;
                   }
                   else {FLAG=0;}
               }
               case 21: {
                   detA = a22*a33-a23*a32;
                   x2 = (a33*b2-a23*b3)/detA;
                   x3 = (a22*b3-a32*b2)/detA;
                   if (x2>=0 && x2<=1 && x3>=0 && x3<=1 && a12*x2+a13*x3>=b1){
                       out_x[i]=0;
                       out_x[i+1]=x2;
                       out_x[i+2]=x3;
                       FLAG=21;
                   }
                   else {FLAG=0;}
               }
               case 22: {
                   detA = a22*a33-a23*a32;
                   x2 = (a33*(b2-a21)-a23*(b3-a31))/detA;
                   x3 = (a22*(b3-a31)-a32*(b2-a21))/detA;
                   if (x2>=0 && x2<=1 && x3>=0 && x3<=1 && a11+a12*x2+a13*x3<=b1){
                       out_x[i]=1;
                       out_x[i+1]=x2;
                       out_x[i+2]=x3;
                       FLAG=22;
                   }
                   else {FLAG=0;}
               }
               case 23: {
                   detA = a11*a33-a13*a31;
                   x1 = (a33*b1-a13*b3)/detA;
                   x3 = (a11*b3-a31*b1)/detA;
                   if (x1>=0 && x1<=1 && x3>=0 && x3<=1 && a21*x1+a23*x3>=b2){
                       out_x[i]=x1;
                       out_x[i+1]=0;
                       out_x[i+2]=x3;
                       FLAG=23;
                   }
                   else {FLAG=0;}
               }
               case 24: {
                   detA = a11*a33-a13*a31;
                   x1 = (a33*(b1-a12)-a13*(b3-a32))/detA;
                   x3 = (a11*(b3-a32)-a31*(b1-a12))/detA;
                   if (x1>=0 && x1<=1 && x3>=0 && x3<=1 && a21*x1+a22+a23*x3<=b2){
                       out_x[i]=x1;
                       out_x[i+1]=1;
                       out_x[i+2]=x3;
                       FLAG=24;
                   }
                   else {FLAG=0;}
               }
               case 25: {
                   detA = a11*a22-a12*a21;
                   x1 = (a22*b1-a12*b2)/detA;
                   x2 = (a11*b2-a21*b1)/detA;
                   if (x1>=0 && x1<=1 && x2>=0 && x2<=1 && a31*x1+a32*x2>=b3){
                       out_x[i]=x1;
                       out_x[i+1]=x2;
                       out_x[i+2]=0;
                       FLAG=25;
                   }
                   else {FLAG=0;}
               }
               case 26: {
                   detA = a11*a22-a12*a21;
                   x1 = (a22*(b1-a13)-a12*(b2-a23))/detA;
                   x2 = (a11*(b2-a23)-a21*(b1-a13))/detA;
                   if (x1>=0 && x1<=1 && x2>=0 && x2<=1 && a31*x1+a32*x2+a33<=b3){
                       out_x[i]=x1;
                       out_x[i+1]=x2;
                       out_x[i+2]=1;
                       FLAG=26;
                   }
                   else {FLAG=0;}
               }
               case 27: {
                   detA = a11*a22*a33+a21*a32*a13+a31*a12*a23-a11*a32*a23-a22*a13*a31-a33*a12*a21;
                   x1 = ( (a22*a33-a32*a23)*b1-(a12*a33-a32*a13)*b2+(a12*a23-a22*a13)*b3)/detA;
                   x2 = (-(a21*a33-a31*a23)*b1+(a11*a33-a31*a13)*b2-(a11*a23-a21*a13)*b3)/detA;
                   x3 = ( (a21*a32-a31*a22)*b1-(a11*a32-a31*a12)*b2+(a11*a22-a21*a12)*b3)/detA;
                   if (x1>=0 && x1<=1 && x2>=0 && x2<=1 && x3>=0 && x3<=1){
                       out_x[i]=x1;
                       out_x[i+1]=x2;
                       out_x[i+2]=x3;
                       FLAG=27;
                   }
                   else {FLAG=0;}
               }
               
           }
           //if the last choice does not match this time, 9 cases again
           if (FLAG==0){
               // first discuss three 0 or 1, 8 cases
               if (b1<=0 && b2<=0 && b3<=0){
                   out_x[i]=0;
                   out_x[i+1]=0;
                   out_x[i+2]=0;
                   FLAG=1;
               }
               else if (a13>=b1 && a23>=b2 && a33<=b3){
                   out_x[i]=0;
                   out_x[i+1]=0;
                   out_x[i+2]=1;
                   FLAG=2;
               }
               else if (a12>=b1 && a22<=b2 && a32>=b3){
                   out_x[i]=0;
                   out_x[i+1]=1;
                   out_x[i+2]=0;
                   FLAG=3;
               }
               else if (a12+a13>=b1 && a22+a23<=b2 && a32+a33<=b3){
                   out_x[i]=0;
                   out_x[i+1]=1;
                   out_x[i+2]=1;
                   FLAG=4;
               }
               else if (a11<=b1 && a12>=b2 && a13>=b3){
                   out_x[i]=1;
                   out_x[i+1]=0;
                   out_x[i+2]=0;
                   FLAG=5;
               }
               else if (a11+a13<=b1 && a21+a23>=b2 && a31+a33<=b3){
                   out_x[i]=1;
                   out_x[i+1]=0;
                   out_x[i+2]=1;
                   FLAG=6;
               }
               else if (a11+a12<=b1 && a21+a22<=b2 && a31+a32>=b3){
                   out_x[i]=1;
                   out_x[i+1]=1;
                   out_x[i+2]=0;
                   FLAG=7;
               }
               else if (a11+a12+a13<=b1 && a21+a22+a23<=b2 && a31+a32+a33<=b3){
                   out_x[i]=1;
                   out_x[i+1]=1;
                   out_x[i+2]=1;
                   FLAG=8;
               }
               // second discuss two 0s and two 1s, 6 cases
               if (FLAG==0){
                   x3 = b3/a33;
                   if (x3>=0 && x3<=1 && a13*x3>=b1 && a23*x3>=b2){
                       out_x[i]=0;
                       out_x[i+1]=0;
                       out_x[i+2]=x3;
                       FLAG=9;
                   }
                   if (FLAG==0){
                       x2 = b2/a22;
                       if (x2>=0 && x2<=1 && a12*x2>=b1 && a32*x2>=b3){
                           out_x[i]=0;
                           out_x[i+1]=x2;
                           out_x[i+2]=0;
                           FLAG=10;
                       }
                       if (FLAG==0){
                           x1 = b1/a11;
                           if (x1>=0 && x1<=1 && a21*x1>=b2 && a31*x1>=b3){
                               out_x[i]=x1;
                               out_x[i+1]=0;
                               out_x[i+2]=0;
                               FLAG=11;
                           }
                           if (FLAG==0){
                               x1=(b1-a12-a13)/a11;
                               if (x1>=0 && x1<=1 && a21*x1+a22+a23<=b2 && a31*x1+a32+a33<=b3){
                                   out_x[i]=x1;
                                   out_x[i+1]=1;
                                   out_x[i+2]=1;
                                   FLAG=12;
                               }
                               if (FLAG==0){
                                   x2=(b2-a21-a23)/a22;
                                   if (x2>=0 && x2<=1 && a11+a12*x2+a13<=b1 && a31+a32*x2+a33<=b3){
                                       out_x[i]=1;
                                       out_x[i+1]=x2;
                                       out_x[i+2]=1;
                                       FLAG=13;
                                   }
                                   if (FLAG==0){
                                       x3=(b3-a31-a32)/a33;
                                       if (x3>=0 && x3<=1 && a11+a12+a13*x3<=b1 && a21+a22+a23*x3<=b2){
                                           out_x[i]=1;
                                           out_x[i+1]=1;
                                           out_x[i+2]=x3;
                                           FLAG=14;
                                       }
                                   }
                               }
                           }
                       }
                   }
               }
               // third discuss one 0 and one 1
               if (FLAG==0){
                   x1=(b1-a13)/a11;
                   if (x1>=0 && x1<=1 && a21*x1+a23>=b2 && a31*x1+a33<=b3){
                       out_x[i]=x1;
                       out_x[i+1]=0;
                       out_x[i+2]=1;
                       FLAG=15;
                   }
                   if (FLAG==0){
                       x1 = (b1-a12)/a11;
                       if (x1>=0 && x1<=1 && a21*x1+a22<=b2 && a31*x1+a32>=b3){
                           out_x[i]=x1;
                           out_x[i+1]=1;
                           out_x[i+2]=0;
                           FLAG=16;
                       }
                       if (FLAG==0){
                           x2 = (b2-a23)/a22;
                           if (x2>=0 && x2<=1 && a12*x2+a13>=b1 && a32*x2+a33<=b3){
                               out_x[i]=0;
                               out_x[i+1]=x2;
                               out_x[i+2]=1;
                               FLAG=17;
                           }
                           if (FLAG==0){
                               x2=(b2-a21)/a22;
                               if (x2>=0 && x2<=1 && a11+a12*x2<=b1 && a31+a32*x2>=b3){
                                   out_x[i]=1;
                                   out_x[i+1]=x2;
                                   out_x[i+2]=0;
                                   FLAG=18;
                               }
                               if (FLAG==0){
                                   x3=(b3-a32)/a33;
                                   if (x3>=0 && x3<=1 && a12+a13*x3>=b1 && a22+a23*x3<=b2){
                                       out_x[i]=0;
                                       out_x[i+1]=1;
                                       out_x[i+2]=x3;
                                       FLAG=19;
                                   }
                                   if (FLAG==0){
                                       x3=(b3-a31)/a33;
                                       if (x3>=0 && x3<=1 && a11+a13*x3<=b1 && a21+a23*x3>=b2){
                                           out_x[i]=1;
                                           out_x[i+1]=0;
                                           out_x[i+2]=x3;
                                           FLAG=20;
                                       }
                                   }
                               }
                           }
                       }
                   }
               }
               // fourth discuss one 1 or one 0
               if (FLAG==0){
                   detA = a22*a33-a23*a32;
                   x2 = (a33*b2-a23*b3)/detA;
                   x3 = (a22*b3-a32*b2)/detA;
                   if (x2>=0 && x2<=1 && x3>=0 && x3<=1 && a12*x2+a13*x3>=b1){
                       out_x[i]=0;
                       out_x[i+1]=x2;
                       out_x[i+2]=x3;
                       FLAG=21;
                   }
                   if (FLAG==0){
                       x2 = (a33*(b2-a21)-a23*(b3-a31))/detA;
                       x3 = (a22*(b3-a31)-a32*(b2-a21))/detA;
                       if (x2>=0 && x2<=1 && x3>=0 && x3<=1 && a11+a12*x2+a13*x3<=b1){
                           out_x[i]=1;
                           out_x[i+1]=x2;
                           out_x[i+2]=x3;
                           FLAG=22;
                       }
                       if (FLAG==0){
                           detA = a11*a33-a13*a31;
                           x1 = (a33*b1-a13*b3)/detA;
                           x3 = (a11*b3-a31*b1)/detA;
                           if (x1>=0 && x1<=1 && x3>=0 && x3<=1 && a21*x1+a23*x3>=b2){
                               out_x[i]=x1;
                               out_x[i+1]=0;
                               out_x[i+2]=x3;
                               FLAG=23;
                           }
                           if (FLAG==0){
                               x1 = (a33*(b1-a12)-a13*(b3-a32))/detA;
                               x3 = (a11*(b3-a32)-a31*(b1-a12))/detA;
                               if (x1>=0 && x1<=1 && x3>=0 && x3<=1 && a21*x1+a22+a23*x3<=b2){
                                   out_x[i]=x1;
                                   out_x[i+1]=1;
                                   out_x[i+2]=x3;
                                   FLAG=24;
                               }
                               if (FLAG==0){
                                   detA = a11*a22-a12*a21;
                                   x1 = (a22*b1-a12*b2)/detA;
                                   x2 = (a11*b2-a21*b1)/detA;
                                   if (x1>=0 && x1<=1 && x2>=0 && x2<=1 && a31*x1+a32*x2>=b3){
                                       out_x[i]=x1;
                                       out_x[i+1]=x2;
                                       out_x[i+2]=0;
                                       FLAG=25;
                                   }
                                   if (FLAG==0){
                                       x1 = (a22*(b1-a13)-a12*(b2-a23))/detA;
                                       x2 = (a11*(b2-a23)-a21*(b1-a13))/detA;
                                       if (x1>=0 && x1<=1 && x2>=0 && x2<=1 && a31*x1+a32*x2+a33<=b3){
                                           out_x[i]=x1;
                                           out_x[i+1]=x2;
                                           out_x[i+2]=1;
                                           FLAG=26;
                                       }
                                   }
                               }
                           }
                       }
                   }
               }
               // last case, no 0 or 1
               if (FLAG==0){
                   // solve 3 dim linear system
                   detA = a11*a22*a33+a21*a32*a13+a31*a12*a23-a11*a32*a23-a22*a13*a31-a33*a12*a21;
                   x1 = ( (a22*a33-a32*a23)*b1-(a12*a33-a32*a13)*b2+(a12*a23-a22*a13)*b3)/detA;
                   x2 = (-(a21*a33-a31*a23)*b1+(a11*a33-a31*a13)*b2-(a11*a23-a21*a13)*b3)/detA;
                   x3 = ( (a21*a32-a31*a22)*b1-(a11*a32-a31*a12)*b2+(a11*a22-a21*a12)*b3)/detA;
                   if (x1>=0 && x1<=1 && x2>=0 && x2<=1 && x3>=0 && x3<=1){
                       out_x[i]=x1;
                       out_x[i+1]=x2;
                       out_x[i+2]=x3;
                       FLAG=27;
                   }
                   else {
                       mexPrintf("no update, check code\n");
                   }

               }
           }
           labels[i] = FLAG;
           //update temporal grad
           time2 = clock();////////time//////////
           for (j=0;j<in_d;j++){
               grad[j] = grad[j]+in_A[i*in_d+j]*out_x[i]+in_A[(i+1)*in_d+j]*out_x[i+1]+in_A[(i+2)*in_d+j]*out_x[i+2];
           }
           time3 = clock();
           dt1+=(double)(time1-time0);///((clock_t)1000);
           dt2+=(double)(time2-time1);
           dt3+=(double)(time3-time2);
           dt4+=(double)(time3-time0);
       }
       if (in_d%3==1){
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
       else if (in_d%3==2){
           i=in_d-2;
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
   double *r = mxGetPr(plhs[1]);if(r==NULL){mexErrMsgTxt("pointer r is null");  return;}
   for (i=0;i<epoch;i++){
       r[i]=out_r[i];
   }
   delete grad; delete out_r;
   mexPrintf("dt1 = %.5f, dt2 = %.5f, dt3 = %.5f, dt4 = %.5f\n",dt1,dt2,dt3,dt4);
   mexPrintf("End of CBCD size 3\n");
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