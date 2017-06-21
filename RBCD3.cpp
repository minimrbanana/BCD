#include <math.h>
#include "mex.h"
#define EPSILON 2.220446e-16
/**
 * RBCD size 3 with general constraints
 * random choose from d-2 blocks with overlap
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
     * [9] random parameter alpha
     */
    // first check the number of input values
    if(nrhs!=9){mexErrMsgTxt("Input Value Error");  return;}
    // [1]
    double *in_A;
    mwIndex *irs;// for sparse matrix
    mwIndex *jcs;// for sparse matrix
    // [2]
    double *in_b;
    // [3]
    long int in_d;
    // [4]
    long int in_max_iter;
    // [5]
    double in_precision;
    // [6]
    double lower;
    // [7]
    double upper;
    // [8]
    double in_init;
    // [9]
    double in_alpha;
    /** ouput args
     * [1] the minimizer out_x
     * [2] the KKT condition residual
     *     the length will be adjusted in the end,
     *     and final output is *r
     */
    // [1]
    double *out_x;// the minimizer
    // [2]
    double residual;
    //parameters in the function
    long int i,j,epoch;//loop
    double df;
    double RV;//random variable for RBCD
    double maxA3;//max element in the size 1 block
    //get input args
    // [1]
    in_A = mxGetPr(prhs[0]);if(in_A==NULL){mexErrMsgTxt("pointer in_A is null");  return;}
    irs = mxGetIr(prhs[0]);if(irs==NULL){mexErrMsgTxt("pointer irs is null");  return;}
    jcs = mxGetJc(prhs[0]);if(jcs==NULL){mexErrMsgTxt("pointer jcs is null");  return;}
    // [2]
    in_b = mxGetPr(prhs[1]);if(in_b==NULL){mexErrMsgTxt("pointer in_b is null");  return;}
    // [3]
    in_d = mxGetScalar(prhs[2]);if(in_d<3){mexErrMsgTxt("dimension error");  return;}
    // [4]
    in_max_iter = mxGetScalar(prhs[3]);if(in_max_iter<=0){mexErrMsgTxt("max_iter can not be 0");  return;}
    // [5]
    in_precision = mxGetScalar(prhs[4]);if(in_max_iter<=0){mexErrMsgTxt("precisionr can not be 0");  return;} 
    // [6]
    lower = mxGetScalar(prhs[5]);
    // [7]
    upper = mxGetScalar(prhs[6]);
    // [8]
    in_init = mxGetScalar(prhs[7]);
    // [9]
    in_alpha = mxGetScalar(prhs[8]);
    /* make sure that the upper bound is larger than the lower bound
     * if the lower bound is greater, 
     * then we solve for an unconstrained problem
     */
    if (lower>=upper){mexPrintf("Bounds Error, Results without constraints\n");in_init=0;}
    /* if the in_init is out of the bound, set it as lower */
    else if ((in_init<lower)||(in_init>upper)){in_init=lower;}
    /* Non-Zero elements, is the value of last entry of jcs
     * lengths of in_A and irs are both NZmax
     * length of jcs is in_d + 1, and the last entry of jcs has value NZmax
    */
    double sparsity = jcs[in_d]/double(in_d*in_d);
    mexPrintf("RBCD size 3.cpp...Sparsity = %.5f.\n",sparsity);
    // [1] allocate output, and init as all in_init
    plhs[0] = mxCreateDoubleMatrix(in_d,1,mxREAL);
    out_x = mxGetPr(plhs[0]);if(out_x==NULL){mexErrMsgTxt("pointer out_x is null");  return;} 
    for (i=0;i<in_d;i++){
        out_x[i] = in_init;
    }
    // [2] pre-allocate output of residual, length as max_iter
    double* out_r=new double[in_max_iter]; if(out_r==NULL){mexErrMsgTxt("pointer out_r is null");  return;} 
    // in-code parameters
    //allocate gradient, will delete later
    double* grad=new double[in_d];  if(grad==NULL){mexErrMsgTxt("pointer grad is null");  return;}
     /*allocate diagonal, for solving small block problem, 
     * diag_A0 is the main diagonal, diag_A1/2 are the one below(and above, as symmetric)
     * although the length of diag_A1/2 is in_d-1/in_d-2, set them as in_d, 
     * the in_d th element will be 0 and not be used.
     */
    double* diag_A0=new double[in_d];  if(diag_A0==NULL){mexErrMsgTxt("pointer diag_A0 is null");  return;} 
    double* diag_A1=new double[in_d];  if(diag_A1==NULL){mexErrMsgTxt("pointer diag_A1 is null");  return;} 
    double* diag_A2=new double[in_d];  if(diag_A2==NULL){mexErrMsgTxt("pointer diag_A2 is null");  return;} 
    /*grad and residual of init loop
     *out_x is initialized as in_init
     *in the loop the elements from the diagonal are also extracted
    */
    // init gradient as 0s
    // extract i th element in diagonal and below
    for (i=0;i<in_d;i++){
        grad[i]=0;
        // i th element in diagonal and below
        diag_A0[i]=0;
        diag_A1[i]=0;
        diag_A2[i]=0;
        for (j=jcs[i];j<jcs[i+1];j++){
            if (irs[j]==i){
                diag_A0[i]=in_A[j];
                //mexPrintf("diag_A0[%d]=%.5f;\n",i,in_A[j]);
            }
            else if (irs[j]==i+1){
                diag_A1[i]=in_A[j];
                //mexPrintf("diag_A1[%d]=%.5f;\n",i,in_A[j]);
            }
            else if (irs[j]==i+2){
                diag_A2[i]=in_A[j];
                //mexPrintf("diag_A2[%d]=%.5f;\n",i,in_A[j]);
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
    // parameters for size 3 matrix
    double a11,a12,a13,a21,a22,a23,a31,a32,a33;
    double b1,b2,b3,x1,x2,x3,detA;
    // FLAG stores the label of last choice of each block
    // label stores the last choice of all the blocks in last epoch
    long int* labels=new long int[in_d];  if(labels==NULL){mexErrMsgTxt("pointer labels is null");  return;} ;
    for (i=0;i<in_d;i++){
        labels[i]=0;
    }
    int FLAG;
    /* if the bounds are defined, and lower<upper
     * we take them as [lower,upper]
     * then do the following
     */
    if (lower<upper){
        while ((residual>in_precision)&&(epoch<in_max_iter)){
            // KKT condition is calculated every in_d/2 updates, i.e. one epoch
            for (long int loop_number=0;loop_number<in_d/3;loop_number++){
                // get the random index, in the range of [0,in_d-3]
                // i=[0,3,6,9....]
                RV = ((double)rand())/((double)RAND_MAX+1.0);
                i= 3 * (int) (RV*(in_d/3.0));
                if (i<in_d-2){
                    // calc temporal grad
                    // three for loops to reuse current memory
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i];
                    }
                    for (j=jcs[i+1];j<jcs[i+2];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i+1];
                    }
                    for (j=jcs[i+2];j<jcs[i+3];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i+2];
                    }
                    // update x(i)
                    // define size 3 block
                    a11=diag_A0[i  ]; a12=diag_A1[i  ]; a13=diag_A2[i  ];
                    a21=diag_A1[i  ]; a22=diag_A0[i+1]; a23=diag_A1[i+1];
                    a31=diag_A2[i  ]; a32=diag_A1[i+1]; a33=diag_A0[i+2];
                    b1 =in_b[i]  -grad[i];
                    b2 =in_b[i+1]-grad[i+1];
                    b3 =in_b[i+2]-grad[i+2];
                    // decission tree
                    FLAG = labels[i];
                    switch (FLAG){
                        case 1: {
                            if (b1<=(a11+a12+a13)*lower && b2<=(a21+a22+a23)*lower && b3<=(a31+a32+a33)*lower){
                                out_x[i]  =lower;
                                out_x[i+1]=lower;
                                out_x[i+2]=lower;
                                FLAG=1;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 2: {
                            if (b1<=(a11+a12)*lower+a13*upper && b2<=(a21+a22)*lower+a23*upper && b3>=(a31+a32)*lower+a33*upper){
                                out_x[i]  =lower;
                                out_x[i+1]=lower;
                                out_x[i+2]=upper;
                                FLAG=2;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 3: {
                            if (b1<=(a11+a13)*lower+a12*upper && b2>=(a21+a23)*lower+a22*upper && b3<=(a31+a33)*lower+a32*upper){
                                out_x[i]  =lower;
                                out_x[i+1]=upper;
                                out_x[i+2]=lower;
                                FLAG=3;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 4: {
                            if (b1<=a11*lower+(a12+a13)*upper && b2>=a21*lower+(a22+a23)*upper && b3>=a31*lower+(a32+a33)*upper){
                                out_x[i]  =lower;
                                out_x[i+1]=upper;
                                out_x[i+2]=upper;
                                FLAG=4;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 5: {
                            if (b1>=a11*upper+(a12+a13)*lower && b2<=a21*upper+(a22+a23)*lower && b3<=a31*upper+(a32+a33)*lower){
                                out_x[i]  =upper;
                                out_x[i+1]=lower;
                                out_x[i+2]=lower;
                                FLAG=5;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 6: {
                            if (b1>=(a11+a12)*upper+a13*lower && b2>=(a21+a22)*upper+a23*lower && b3<=(a31+a32)*upper+a33*lower){
                                out_x[i]  =upper;
                                out_x[i+1]=upper;
                                out_x[i+2]=lower;
                                FLAG=7;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 7: {
                            if (b1>=(a11+a12)*upper+a13*lower && b2>=(a21+a22)*upper+a23*lower && b3<=(a31+a32)*upper+a33*lower){
                                out_x[i]  =upper;
                                out_x[i+1]=upper;
                                out_x[i+2]=lower;
                                FLAG=7;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 8: {
                            if (b1>=(a11+a12+a13)*upper && b2>=(a21+a22+a23)*upper && b3>=(a31+a32+a33)*upper){
                                out_x[i]  =upper;
                                out_x[i+1]=upper;
                                out_x[i+2]=upper;
                                FLAG=8;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 9: {
                            x3 = (b3-(a31+a32)*lower)/a33;
                            if (x3>=lower && x3<=upper && (a11+a12)*lower+a13*x3>=b1 && (a21+a22)*lower+a23*x3>=b2){
                                out_x[i]  =lower;
                                out_x[i+1]=lower;
                                out_x[i+2]=x3;
                                FLAG=9;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 10: {
                            x2 = (b2-(a21+a23)*lower)/a22;
                            if (x2>=lower && x2<=upper && (a11+a13)*lower+a12*x2>=b1 && (a31+a33)*lower+a32*x2>=b3){
                                out_x[i]  =lower;
                                out_x[i+1]=x2;
                                out_x[i+2]=lower;
                                FLAG=10;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 11: {
                            x1 = (b1-(a12+a13)*lower)/a11;
                            if (x1>=lower && x1<=upper && (a22+a23)*lower+a21*x1>=b2 && (a32+a33)*lower+a31*x1>=b3){
                                out_x[i]  =x1;
                                out_x[i+1]=lower;
                                out_x[i+2]=lower;
                                FLAG=11;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 12: {
                            x1=(b1-(a12+a13)*upper)/a11;
                            if (x1>=lower && x1<=upper && a21*x1+(a22+a23)*upper<=b2 && a31*x1+(a32+a33)*upper<=b3){
                                out_x[i]  =x1;
                                out_x[i+1]=upper;
                                out_x[i+2]=upper;
                                FLAG=12;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 13: {
                            x2=(b2-(a21+a23)*upper)/a22;
                            if (x2>=lower && x2<=upper && (a11+a13)*upper+a12*x2<=b1 && (a31+a33)*upper+a32*x2<=b3){
                                out_x[i]  =upper;
                                out_x[i+1]=x2;
                                out_x[i+2]=upper;
                                FLAG=13;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 14: {
                            x3=(b3-(a31+a32)*upper)/a33;
                            if (x3>=lower && x3<=upper && (a11+a12)*upper+a13*x3<=b1 && (a21+a22)*upper+a23*x3<=b2){
                                out_x[i]  =upper;
                                out_x[i+1]=upper;
                                out_x[i+2]=x3;
                                FLAG=14;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 15: {
                            x1=(b1-a12*lower-a13*upper)/a11;
                            if (x1>=lower && x1<=upper && a21*x1+a22*lower+a23*upper>=b2 && a31*x1+a32*lower+a33*upper<=b3){
                                out_x[i]  =x1;
                                out_x[i+1]=lower;
                                out_x[i+2]=upper;
                                FLAG=15;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 16: {
                            x1 = (b1-a12*upper-a13*lower)/a11;
                            if (x1>=lower && x1<=upper && a21*x1+a22*upper+a23*lower<=b2 && a31*x1+a32*upper+a33*lower>=b3){
                                out_x[i]  =x1;
                                out_x[i+1]=upper;
                                out_x[i+2]=lower;
                                FLAG=16;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 17: {
                            x2 = (b2-a23*upper-a21*lower)/a22;
                            if (x2>=lower && x2<=upper && a12*x2+a13*upper+a11*lower>=b1 && a32*x2+a33*upper+a31*lower<=b3){
                                out_x[i]  =lower;
                                out_x[i+1]=x2;
                                out_x[i+2]=upper;
                                FLAG=17;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 18: {
                            x2=(b2-a21*upper-a23*lower)/a22;
                            if (x2>=lower && x2<=upper && a11*upper+a13*lower+a12*x2<=b1 && a31*upper+a33*lower+a32*x2>=b3){
                                out_x[i]  =upper;
                                out_x[i+1]=x2;
                                out_x[i+2]=lower;
                                FLAG=18;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 19: {
                            x3=(b3-a32*upper-a31*lower)/a33;
                            if (x3>=lower && x3<=upper && a11*lower+a12*upper+a13*x3>=b1 && a21*lower+a22*upper+a23*x3<=b2){
                                out_x[i]  =lower;
                                out_x[i+1]=upper;
                                out_x[i+2]=x3;
                                FLAG=19;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 20: {
                            x3=(b3-a31*upper-a32*lower)/a33;
                            if (x3>=lower && x3<=upper && a11*upper+a12*lower+a13*x3<=b1 && a21*upper+a22*lower+a23*x3>=b2){
                                out_x[i]  =upper;
                                out_x[i+1]=lower;
                                out_x[i+2]=x3;
                                FLAG=20;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 21: {
                            detA = a22*a33-a23*a32;
                            x2 = (a33*(b2-a21*lower)-a23*(b3-a31*lower))/detA;
                            x3 = (a22*(b3-a31*lower)-a32*(b2-a21*lower))/detA;
                            if (x2>=lower && x2<=upper && x3>=lower && x3<=upper && a11*lower+a12*x2+a13*x3>=b1){
                                out_x[i]  =lower;
                                out_x[i+1]=x2;
                                out_x[i+2]=x3;
                                FLAG=21;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 22: {
                            detA = a22*a33-a23*a32;
                            x2 = (a33*(b2-a21*upper)-a23*(b3-a31*upper))/detA;
                            x3 = (a22*(b3-a31*upper)-a32*(b2-a21*upper))/detA;
                            if (x2>=lower && x2<=upper && x3>=lower && x3<=upper && a11*upper+a12*x2+a13*x3<=b1){
                                out_x[i]  =upper;
                                out_x[i+1]=x2;
                                out_x[i+2]=x3;
                                FLAG=22;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 23: {
                            detA = a11*a33-a13*a31;
                            x1 = (a33*(b1-a12*lower)-a13*(b3-a32*lower))/detA;
                            x3 = (a11*(b3-a32*lower)-a31*(b1-a12*lower))/detA;
                            if (x1>=lower && x1<=upper && x3>=lower && x3<=upper && a22*lower+a21*x1+a23*x3>=b2){
                                out_x[i]  =x1;
                                out_x[i+1]=lower;
                                out_x[i+2]=x3;
                                FLAG=23;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 24: {
                            detA = a11*a33-a13*a31;
                            x1 = (a33*(b1-a12*upper)-a13*(b3-a32*upper))/detA;
                            x3 = (a11*(b3-a32*upper)-a31*(b1-a12*upper))/detA;
                            if (x1>=lower && x1<=upper && x3>=lower && x3<=upper && a21*x1+a22*upper+a23*x3<=b2){
                                out_x[i]  =x1;
                                out_x[i+1]=upper;
                                out_x[i+2]=x3;
                                FLAG=24;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 25: {
                            detA = a11*a22-a12*a21;
                            x1 = (a22*(b1-a13*lower)-a12*(b2-a23*lower))/detA;
                            x2 = (a11*(b2-a23*lower)-a21*(b1-a13*lower))/detA;
                            if (x1>=lower && x1<=upper && x2>=lower && x2<=upper && a31*x1+a32*x2+a33*lower>=b3){
                                out_x[i]  =x1;
                                out_x[i+1]=x2;
                                out_x[i+2]=lower;
                                FLAG=25;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 26: {
                            detA = a11*a22-a12*a21;
                            x1 = (a22*(b1-a13*upper)-a12*(b2-a23*upper))/detA;
                            x2 = (a11*(b2-a23*upper)-a21*(b1-a13*upper))/detA;
                            if (x1>=lower && x1<=upper && x2>=lower && x2<=upper && a31*x1+a32*x2+a33*upper<=b3){
                                out_x[i]  =x1;
                                out_x[i+1]=x2;
                                out_x[i+2]=upper;
                                FLAG=26;
                            }
                            else {FLAG=0;}
                            break;
                        }
                        case 27: {
                            detA = a11*a22*a33+a21*a32*a13+a31*a12*a23-a11*a32*a23-a22*a13*a31-a33*a12*a21;
                            if(detA==0){mexErrMsgTxt("Input Matrix is not positive definite");  return;}
                            x1 = ( (a22*a33-a32*a23)*b1-(a12*a33-a32*a13)*b2+(a12*a23-a22*a13)*b3)/detA;
                            x2 = (-(a21*a33-a31*a23)*b1+(a11*a33-a31*a13)*b2-(a11*a23-a21*a13)*b3)/detA;
                            x3 = ( (a21*a32-a31*a22)*b1-(a11*a32-a31*a12)*b2+(a11*a22-a21*a12)*b3)/detA;
                            if (x1>=lower && x1<=upper && x2>=lower && x2<=upper && x3>=lower && x3<=upper){
                                out_x[i]  =x1;
                                out_x[i+1]=x2;
                                out_x[i+2]=x3;
                                FLAG=27;
                            }
                            else {FLAG=0;}
                            break;
                        }

                    }
                    //if the last choice does not match this time, 27 cases again
                    if (FLAG==0){
                        // first discuss three lower or upper, 8 cases
                        if (b1<=(a11+a12+a13)*lower && b2<=(a21+a22+a23)*lower && b3<=(a31+a32+a33)*lower){
                            out_x[i]  =lower;
                            out_x[i+1]=lower;
                            out_x[i+2]=lower;
                            FLAG=1;
                        }
                        else if (b1<=(a11+a12)*lower+a13*upper && b2<=(a21+a22)*lower+a23*upper && b3>=(a31+a32)*lower+a33*upper){
                            out_x[i]  =lower;
                            out_x[i+1]=lower;
                            out_x[i+2]=upper;
                            FLAG=2;
                        }
                        else if (b1<=(a11+a13)*lower+a12*upper && b2>=(a21+a23)*lower+a22*upper && b3<=(a31+a33)*lower+a32*upper){
                            out_x[i]  =lower;
                            out_x[i+1]=upper;
                            out_x[i+2]=lower;
                            FLAG=3;
                        }
                        else if (b1<=a11*lower+(a12+a13)*upper && b2>=a21*lower+(a22+a23)*upper && b3>=a31*lower+(a32+a33)*upper){
                            out_x[i]  =lower;
                            out_x[i+1]=upper;
                            out_x[i+2]=upper;
                            FLAG=4;
                        }
                        else if (b1>=a11*upper+(a12+a13)*lower && b2<=a21*upper+(a22+a23)*lower && b3<=a31*upper+(a32+a33)*lower){
                            out_x[i]  =upper;
                            out_x[i+1]=lower;
                            out_x[i+2]=lower;
                            FLAG=5;
                        }
                        else if (b1>=(a11+a13)*upper+a12*lower && b2<=(a21+a23)*upper+a22*lower && b3>=(a31+a33)*upper+a32*lower){
                            out_x[i]  =upper;
                            out_x[i+1]=lower;
                            out_x[i+2]=upper;
                            FLAG=6;
                        }
                        else if (b1>=(a11+a12)*upper+a13*lower && b2>=(a21+a22)*upper+a23*lower && b3<=(a31+a32)*upper+a33*lower){
                            out_x[i]  =upper;
                            out_x[i+1]=upper;
                            out_x[i+2]=lower;
                            FLAG=7;
                        }
                        else if (b1>=(a11+a12+a13)*upper && b2>=(a21+a22+a23)*upper && b3>=(a31+a32+a33)*upper){
                            out_x[i]  =upper;
                            out_x[i+1]=upper;
                            out_x[i+2]=upper;
                            FLAG=8;
                        }
                        // second discuss two lower and two upper, 6 cases
                        if (FLAG==0){
                            x3 = (b3-(a31+a32)*lower)/a33;
                            if (x3>=lower && x3<=upper && (a11+a12)*lower+a13*x3>=b1 && (a21+a22)*lower+a23*x3>=b2){
                                out_x[i]  =lower;
                                out_x[i+1]=lower;
                                out_x[i+2]=x3;
                                FLAG=9;
                            }
                            if (FLAG==0){
                                x2 = (b2-(a21+a23)*lower)/a22;
                                if (x2>=lower && x2<=upper && (a11+a13)*lower+a12*x2>=b1 && (a31+a33)*lower+a32*x2>=b3){
                                    out_x[i]  =lower;
                                    out_x[i+1]=x2;
                                    out_x[i+2]=lower;
                                    FLAG=10;
                                }
                                if (FLAG==0){
                                    x1 = (b1-(a12+a13)*lower)/a11;
                                    if (x1>=lower && x1<=upper && (a22+a23)*lower+a21*x1>=b2 && (a32+a33)*lower+a31*x1>=b3){
                                        out_x[i]  =x1;
                                        out_x[i+1]=lower;
                                        out_x[i+2]=lower;
                                        FLAG=11;
                                    }
                                    if (FLAG==0){
                                        x1=(b1-(a12+a13)*upper)/a11;
                                        if (x1>=lower && x1<=upper && a21*x1+(a22+a23)*upper<=b2 && a31*x1+(a32+a33)*upper<=b3){
                                            out_x[i]  =x1;
                                            out_x[i+1]=upper;
                                            out_x[i+2]=upper;
                                            FLAG=12;
                                        }
                                        if (FLAG==0){
                                            x2=(b2-(a21+a23)*upper)/a22;
                                            if (x2>=lower && x2<=upper && (a11+a13)*upper+a12*x2<=b1 && (a31+a33)*upper+a32*x2<=b3){
                                                out_x[i]  =upper;
                                                out_x[i+1]=x2;
                                                out_x[i+2]=upper;
                                                FLAG=13;
                                            }
                                            if (FLAG==0){
                                                x3=(b3-(a31+a32)*upper)/a33;
                                                if (x3>=lower && x3<=upper && (a11+a12)*upper+a13*x3<=b1 && (a21+a22)*upper+a23*x3<=b2){
                                                    out_x[i]  =upper;
                                                    out_x[i+1]=upper;
                                                    out_x[i+2]=x3;
                                                    FLAG=14;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        // third discuss one lower and one upper
                        if (FLAG==0){
                            x1=(b1-a12*lower-a13*upper)/a11;
                            if (x1>=lower && x1<=upper && a21*x1+a22*lower+a23*upper>=b2 && a31*x1+a32*lower+a33*upper<=b3){
                                out_x[i]  =x1;
                                out_x[i+1]=lower;
                                out_x[i+2]=upper;
                                FLAG=15;
                            }
                            if (FLAG==0){
                                x1 = (b1-a12*upper-a13*lower)/a11;
                                if (x1>=lower && x1<=upper && a21*x1+a22*upper+a23*lower<=b2 && a31*x1+a32*upper+a33*lower>=b3){
                                    out_x[i]  =x1;
                                    out_x[i+1]=upper;
                                    out_x[i+2]=lower;
                                    FLAG=16;
                                }
                                if (FLAG==0){
                                    x2 = (b2-a23*upper-a21*lower)/a22;
                                    if (x2>=lower && x2<=upper && a12*x2+a13*upper+a11*lower>=b1 && a32*x2+a33*upper+a31*lower<=b3){
                                        out_x[i]  =lower;
                                        out_x[i+1]=x2;
                                        out_x[i+2]=upper;
                                        FLAG=17;
                                    }
                                    if (FLAG==0){
                                        x2=(b2-a21*upper-a23*lower)/a22;
                                        if (x2>=lower && x2<=upper && a11*upper+a13*lower+a12*x2<=b1 && a31*upper+a33*lower+a32*x2>=b3){
                                            out_x[i]  =upper;
                                            out_x[i+1]=x2;
                                            out_x[i+2]=lower;
                                            FLAG=18;
                                        }
                                        if (FLAG==0){
                                            x3=(b3-a32*upper-a31*lower)/a33;
                                            if (x3>=lower && x3<=upper && a11*lower+a12*upper+a13*x3>=b1 && a21*lower+a22*upper+a23*x3<=b2){
                                                out_x[i]  =lower;
                                                out_x[i+1]=upper;
                                                out_x[i+2]=x3;
                                                FLAG=19;
                                            }
                                            if (FLAG==0){
                                                x3=(b3-a31*upper-a32*lower)/a33;
                                                if (x3>=lower && x3<=upper && a11*upper+a12*lower+a13*x3<=b1 && a21*upper+a22*lower+a23*x3>=b2){
                                                    out_x[i]  =upper;
                                                    out_x[i+1]=lower;
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
                            x2 = (a33*(b2-a21*lower)-a23*(b3-a31*lower))/detA;
                            x3 = (a22*(b3-a31*lower)-a32*(b2-a21*lower))/detA;
                            if (x2>=lower && x2<=upper && x3>=lower && x3<=upper && a11*lower+a12*x2+a13*x3>=b1){
                                out_x[i]  =lower;
                                out_x[i+1]=x2;
                                out_x[i+2]=x3;
                                FLAG=21;
                            }
                            if (FLAG==0){
                                x2 = (a33*(b2-a21*upper)-a23*(b3-a31*upper))/detA;
                                x3 = (a22*(b3-a31*upper)-a32*(b2-a21*upper))/detA;
                                if (x2>=lower && x2<=upper && x3>=lower && x3<=upper && a11*upper+a12*x2+a13*x3<=b1){
                                    out_x[i]  =upper;
                                    out_x[i+1]=x2;
                                    out_x[i+2]=x3;
                                    FLAG=22;
                                }
                                if (FLAG==0){
                                    detA = a11*a33-a13*a31;
                                    x1 = (a33*(b1-a12*lower)-a13*(b3-a32*lower))/detA;
                                    x3 = (a11*(b3-a32*lower)-a31*(b1-a12*lower))/detA;
                                    if (x1>=lower && x1<=upper && x3>=lower && x3<=upper && a22*lower+a21*x1+a23*x3>=b2){
                                        out_x[i]  =x1;
                                        out_x[i+1]=lower;
                                        out_x[i+2]=x3;
                                        FLAG=23;
                                    }
                                    if (FLAG==0){
                                        x1 = (a33*(b1-a12*upper)-a13*(b3-a32*upper))/detA;
                                        x3 = (a11*(b3-a32*upper)-a31*(b1-a12*upper))/detA;
                                        if (x1>=lower && x1<=upper && x3>=lower && x3<=upper && a21*x1+a22*upper+a23*x3<=b2){
                                            out_x[i]  =x1;
                                            out_x[i+1]=upper;
                                            out_x[i+2]=x3;
                                            FLAG=24;
                                        }
                                        if (FLAG==0){
                                            detA = a11*a22-a12*a21;
                                            x1 = (a22*(b1-a13*lower)-a12*(b2-a23*lower))/detA;
                                            x2 = (a11*(b2-a23*lower)-a21*(b1-a13*lower))/detA;
                                            if (x1>=lower && x1<=upper && x2>=lower && x2<=upper && a31*x1+a32*x2+a33*lower>=b3){
                                                out_x[i]  =x1;
                                                out_x[i+1]=x2;
                                                out_x[i+2]=lower;
                                                FLAG=25;
                                            }
                                            if (FLAG==0){
                                                x1 = (a22*(b1-a13*upper)-a12*(b2-a23*upper))/detA;
                                                x2 = (a11*(b2-a23*upper)-a21*(b1-a13*upper))/detA;
                                                if (x1>=lower && x1<=upper && x2>=lower && x2<=upper && a31*x1+a32*x2+a33*upper<=b3){
                                                    out_x[i]  =x1;
                                                    out_x[i+1]=x2;
                                                    out_x[i+2]=upper;
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
                            if (x1>=lower && x1<=upper && x2>=lower && x2<=upper && x3>=lower && x3<=upper){
                                out_x[i]  =x1;
                                out_x[i+1]=x2;
                                out_x[i+2]=x3;
                                FLAG=27;
                            }
                            else {
                                mexPrintf("no update, check code\n");
                            }
                        }
                    }
                    //update the labels of blocks
                    labels[i] = FLAG;
                    //update temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] += in_A[j]*out_x[i];
                    }
                    for (j=jcs[i+1];j<jcs[i+2];j++){
                        grad[irs[j]] += in_A[j]*out_x[i+1];
                    }
                    for (j=jcs[i+2];j<jcs[i+3];j++){
                        grad[irs[j]] += in_A[j]*out_x[i+2];
                    }
                }
                else if (i==in_d-1){
                    //calc temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i];
                    }
                    //descent
                    out_x[i] = (in_b[i]-grad[i])/diag_A0[i];
                    //bounds
                    if (out_x[i]>upper){
                        out_x[i] = upper;
                    }
                    if (out_x[i]<lower){
                        out_x[i] = lower;
                    }
                    //update temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] += in_A[j]*out_x[i];
                    }
                }
                else if (i==in_d-2){
                    //calc temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i];
                    }
                    //descent
                    out_x[i] = (in_b[i]-grad[i])/diag_A0[i];
                    //bounds
                    if (out_x[i]>upper){
                        out_x[i] = upper;
                    }
                    if (out_x[i]<lower){
                        out_x[i] = lower;
                    }
                    //update temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] += in_A[j]*out_x[i];
                    }
                    i++;
                    //calc temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i];
                    }
                    //descent
                    out_x[i] = (in_b[i]-grad[i])/diag_A0[i];
                    //bounds
                    if (out_x[i]>upper){
                        out_x[i] = upper;
                    }
                    if (out_x[i]<lower){
                        out_x[i] = lower;
                    }
                    //update temporal grad, actually not needed
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] += in_A[j]*out_x[i];
                    }
                }
            }
            // when finishes (in_d/2) updates, we calculate the true gradient
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
            //get the residual
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
            // KKT condition is calculated every in_d/2 updates, i.e. one epoch
            for (long int loop_number=0;loop_number<in_d/3;loop_number++){
                // get the random index, in the range of [0,in_d-3]
                // i=[0,3,6,9....]
                RV = ((double)rand())/((double)RAND_MAX+1.0);
                i= 3 * (int) (RV*(in_d/3.0));
                if (i<in_d-2){
                    // calc temporal grad
                    // sparse g=g-A(:,i)*x(i) and i+1, i+2
                    for (j=jcs[i  ];j<jcs[i+1];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i];
                    }
                    for (j=jcs[i+1];j<jcs[i+2];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i+1];
                    }
                    for (j=jcs[i+2];j<jcs[i+3];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i+2];
                    }
                    // update x(i)
                    // define size 3 block
                    a11=diag_A0[i  ]; a12=diag_A1[i  ]; a13=diag_A2[i  ];
                    a21=diag_A1[i  ]; a22=diag_A0[i+1]; a23=diag_A1[i+1];
                    a31=diag_A2[i  ]; a32=diag_A1[i+1]; a33=diag_A0[i+2];
                    b1 =in_b[i]  -grad[i];
                    b2 =in_b[i+1]-grad[i+1];
                    b3 =in_b[i+2]-grad[i+2];
                    // solve 3 dim linear system
                    detA = a11*a22*a33+a21*a32*a13+a31*a12*a23-a11*a32*a23-a22*a13*a31-a33*a12*a21;
                    if(detA==0){mexErrMsgTxt("Input Matrix is not positive definite");  return;}
                    out_x[i  ] = ( (a22*a33-a32*a23)*b1-(a12*a33-a32*a13)*b2+(a12*a23-a22*a13)*b3)/detA;
                    out_x[i+1] = (-(a21*a33-a31*a23)*b1+(a11*a33-a31*a13)*b2-(a11*a23-a21*a13)*b3)/detA;
                    out_x[i+2] = ( (a21*a32-a31*a22)*b1-(a11*a32-a31*a12)*b2+(a11*a22-a21*a12)*b3)/detA;
                    // update temporal grad
                    for (j=jcs[i  ];j<jcs[i+1];j++){
                        grad[irs[j]] += in_A[j]*out_x[i];
                    }
                    for (j=jcs[i+1];j<jcs[i+2];j++){
                        grad[irs[j]] += in_A[j]*out_x[i+1];
                    }
                    for (j=jcs[i+2];j<jcs[i+3];j++){
                        grad[irs[j]] += in_A[j]*out_x[i+2];
                    }
                }
                else if (i==in_d-1){
                    // calc temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i];
                    }
                    // descent
                    // without bounds
                    out_x[i] = (in_b[i]-grad[i])/diag_A0[i];
                    // update temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] += in_A[j]*out_x[i];
                    }
                }
                else if (i==in_d-2){
                    // calc temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i];
                    }
                    // descent
                    // without bounds
                    out_x[i] = (in_b[i]-grad[i])/diag_A0[i];
                    // update temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] += in_A[j]*out_x[i];
                    }
                    i++;
                    // calc temporal grad
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] -= in_A[j]*out_x[i];
                    }
                    // descent
                    // without bounds
                    out_x[i] = (in_b[i]-grad[i])/diag_A0[i];
                    // update temporal grad, actually not needed
                    for (j=jcs[i];j<jcs[i+1];j++){
                        grad[irs[j]] += in_A[j]*out_x[i];
                    }
                }
            }
            // when finishes (in_d/2) updates, we calculate the true gradient
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
    double *r = mxGetPr(plhs[1]);if(r==NULL){mexErrMsgTxt("pointer r is null");  return;}
    for (i=0;i<epoch;i++){
        r[i]=out_r[i];
    }
    delete grad; delete out_r; delete diag_A0; delete diag_A1; delete diag_A2;delete labels;
    //mexPrintf("dt1 = %.5f, dt2 = %.5f, dt3 = %.5f, dt4 = %.5f\n",dt1,dt2,dt3,dt4);
    mexPrintf("epoch:%5d, residual=%.15f\nEnd of RBCD size 3.cpp\n",epoch-1,residual);
}