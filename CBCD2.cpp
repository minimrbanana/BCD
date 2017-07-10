#include <math.h>
#include "mex.h"
#define EPSILON 2.220446e-16
/**
 * CBCD size 2 with general constraints
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
    if(nrhs!=8){mexErrMsgTxt("Input Value Error");  return;}
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
    // parameters in the function
    long int i,j,epoch;//loop
    double df;
    // get input args
    in_A = mxGetPr(prhs[0]);if(in_A==NULL){mexErrMsgTxt("pointer in_A is null");  return;}
    irs = mxGetIr(prhs[0]);if(irs==NULL){mexErrMsgTxt("pointer irs is null");  return;}
    jcs = mxGetJc(prhs[0]);if(jcs==NULL){mexErrMsgTxt("pointer jcs is null");  return;}
    in_b = mxGetPr(prhs[1]);if(in_b==NULL){mexErrMsgTxt("pointer in_b is null");  return;}
    in_d = mxGetScalar(prhs[2]);if(in_d<2){mexErrMsgTxt("dimension error");  return;}
    in_max_iter = mxGetScalar(prhs[3]);if(in_max_iter<=0){mexErrMsgTxt("max_iter can not be 0");  return;}
    in_precision = mxGetScalar(prhs[4]);if(in_max_iter<=0){mexErrMsgTxt("precisionr can not be 0");  return;}
    lower = mxGetScalar(prhs[5]);
    upper = mxGetScalar(prhs[6]);
    in_init = mxGetScalar(prhs[7]);
    /* make sure that the upper bound is larger than the lower bound
     * if the lower bound is greater, 
     * then we solve for an unconstrained problem
     */
    if (lower>=upper){mexPrintf("Bounds Error, Results without constraints\n");in_init=0;}
    /* if the in_init is out of the bound, set it as lower */
    else if ((in_init<lower)||(in_init>upper)){in_init=lower;}
    // print begin
    mexPrintf("CBCD size 2.cpp...\n");
    // allocate output, and init as all in_init
    plhs[0] = mxCreateDoubleMatrix(in_d,1,mxREAL);
    out_x = mxGetPr(plhs[0]);if(out_x==NULL){mexErrMsgTxt("pointer out_x is null");  return;} 
    for (i=0;i<in_d;i++){
        out_x[i] = in_init;
    }
    // pre-allocate output of residual, length as max_iter
    double* out_r=new double[in_max_iter]; if(out_r==NULL){mexErrMsgTxt("pointer out_r is null");  return;} 
    // allocate gradient, will delete later
    double* grad=new double[in_d];  if(grad==NULL){mexErrMsgTxt("pointer grad is null");  return;} 
    /* allocate diagonal, for solving small block problem, 
     * diag_A0 is the main diagonal, diag_A1 is the one below(and above, as symmetric)
     * although the length of diag_A1 is in_d-1, set it as in_d, 
     * the in_d th element will be 0 and not be used.
     */
    double* diag_A0=new double[in_d];  if(diag_A0==NULL){mexErrMsgTxt("pointer diag_A0 is null");  return;} 
    double* diag_A1=new double[in_d];  if(diag_A1==NULL){mexErrMsgTxt("pointer diag_A1 is null");  return;} 
    /* grad and residual of init loop
     * out_x is initialized as in_init
     * in the loop the elements from the diagonal are also extracted
     */
    // init gradient as 0s
    // extract i th element in diagonal and below
    for (i=0;i<in_d;i++){
        grad[i]=0;
        // i th in A
        diag_A0[i]=0;
        diag_A1[i]=0;
        for (j=jcs[i];j<jcs[i+1];j++){
            if (irs[j]==i){
                diag_A0[i]=in_A[j];
            }
            else if (irs[j]==i+1){
                diag_A1[i]=in_A[j];
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
    // parameters for size 2 matrix
    double a11,a12,a21,a22;
    double b1,b2,detA2;
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
            for (i=0;i<in_d-1;i=i+2){
                // calc temporal grad
                // sparse g=g-A(:,i)*x(i) and i+1
                for (j=jcs[i];j<jcs[i+1];j++){
                    grad[irs[j]] -= in_A[j]*out_x[i];
                }
                for (j=jcs[i+1];j<jcs[i+2];j++){
                    grad[irs[j]] -= in_A[j]*out_x[i+1];
                }
                // update x(i)
                // define size 2 block
                a11=diag_A0[i];     a12=diag_A1[i  ];
                a21=diag_A1[i];     a22=diag_A0[i+1];
                b1 =in_b[i]  -grad[i];
                b2 =in_b[i+1]-grad[i+1];
                // NINE choices, frequently used variables
                detA2 = a11*a22-a12*a21;
                if(detA2==0){mexErrMsgTxt("Input Matrix is not positive definite");  return;}
                // this switch first check whether the last choice matches or not
                FLAG = labels[i];
                switch (FLAG){
                    case 1: {
                        if (b1<=(a11+a12)*lower && b2<=(a21+a22)*lower){//case 1
                            out_x[i]  =lower;
                            out_x[i+1]=lower;
                            FLAG = 1;//mexPrintf("case1\n");
                        }
                        else {FLAG=0;}
                        break;
                    }
                    case 2: {
                        if (b1>(a11+a12)*lower && b1<a11*upper+a12*lower && a11*b2-a21*b1<=detA2*lower){//case 2
                             out_x[i]  =(b1-a12*lower)/a11;
                             out_x[i+1]=lower;
                             FLAG = 2;//mexPrintf("case2\n");
                        }
                        else {FLAG=0;}
                        break;
                    }
                    case 3: {
                        if(b1>=a11*upper+a12*lower && b2<=a21*upper+a22*lower){//case 3
                             out_x[i]  =upper;
                             out_x[i+1]=lower;
                             FLAG = 3;//mexPrintf("case3\n");
                        }
                        else {FLAG=0;}
                        break;
                    }
                    case 4: {
                        if (b2>(a21+a22)*lower && b2<a21*lower+a22*upper && a22*b1-a12*b2<=detA2*lower){//case 4
                             out_x[i]  =lower;
                             out_x[i+1]=(b2-a21*lower)/a22;
                             FLAG = 4;//mexPrintf("case4\n");
                        }
                        else {FLAG=0;}
                        break;
                    }
                    case 5: {
                        out_x[i  ]=(a22*b1-a12*b2)/detA2;
                        out_x[i+1]=(a11*b2-a21*b1)/detA2;
                        if (out_x[i]>=lower && out_x[i]<=upper && out_x[i+1]>=lower && out_x[i+1]<=upper){
                             FLAG = 5;
                        }
                        else {FLAG=0;}
                        break;
                    }
                    case 6: {
                        if(b2>a21*upper+a22*lower && b2<(a21+a22)*upper && a22*b1-a12*b2>=detA2*upper){//case 6
                             out_x[i]  =upper;
                             out_x[i+1]=(b2-a21*upper)/a22;
                             FLAG = 6;//mexPrintf("case6\n");
                        }
                        else {FLAG=0;}
                        break;
                    }
                    case 7: {
                        if (b1<=a11*lower+a12*upper && b2>=a21*lower+a22*upper){//case 7
                             out_x[i]  =lower;
                             out_x[i+1]=upper;
                             FLAG = 7;//mexPrintf("case7\n");
                        }
                        else {FLAG=0;}
                        break;
                    }
                    case 8: {
                        if (b1>a11*lower+a12*upper && b1<(a11+a12)*upper && a11*b2-a21*b1>=detA2*upper){//case 8
                             out_x[i]  =(b1-a12*upper)/a11;
                             out_x[i+1]=upper;
                             FLAG = 8;//mexPrintf("case8\n");
                        }
                        else {FLAG=0;}
                        break;
                    }
                    case 9: {
                        if(b1>=(a11+a12)*upper && b2>=(a21+a22)*upper){//case 9
                             out_x[i]  =upper;
                             out_x[i+1]=upper;
                             FLAG = 9;//mexPrintf("case9\n");
                        }
                        else {FLAG=0;}
                        break;
                    }
                }
                //if the last choice does not match this time, 9 cases again
                if (FLAG==0){
                    // first assume x2=lower
                    if (b1<=(a11+a12)*lower && b2<=(a21+a22)*lower){//case 1
                            out_x[i]  =lower;
                            out_x[i+1]=lower;
                            FLAG = 1;//mexPrintf("case1\n");
                    }
                    else if (b1>(a11+a12)*lower && b1<a11*upper+a12*lower && a11*b2-a21*b1<=detA2*lower){//case 2
                            out_x[i]  =(b1-a12*lower)/a11;
                            out_x[i+1]=lower;
                            FLAG = 2;//mexPrintf("case2\n");
                    }
                    else if (b1>=a11*upper+a12*lower && b2<=a21*upper+a22*lower){//case 3
                            out_x[i]  =upper;
                            out_x[i+1]=lower;
                            FLAG = 3;//mexPrintf("case3\n");
                    }
                    // x2~=lower, assume x2=upper
                    else if (b1<=a11*lower+a12*upper && b2>=a21*lower+a22*upper){//case 7
                            out_x[i]  =lower;
                            out_x[i+1]=upper;
                            FLAG = 7;//mexPrintf("case7\n");
                    }
                    else if (b1>=(a11+a12)*upper && b2>=(a21+a22)*upper){//case 9
                            out_x[i]  =upper;
                            out_x[i+1]=upper;
                            FLAG = 9;//mexPrintf("case9\n");
                    }
                    else if (b1>a11*lower+a12*upper && b1<(a11+a12)*upper && a11*b2-a21*b1>=detA2*upper){//case 8
                            out_x[i]  =(b1-a12*upper)/a11;;
                            out_x[i+1]=upper;
                            FLAG = 8;//mexPrintf("case8\n");
                    }
                    // x2~=lower & x2~=upper, x2 in (lower,upper)
                    else if (b2>(a21+a22)*lower && b2<a21*lower+a22*upper && a22*b1-a12*b2<=detA2*lower){//case 4
                            out_x[i]  =lower;
                            out_x[i+1]=(b2-a21*lower)/a22;
                            FLAG = 4;//mexPrintf("case4\n");
                    }
                    else if (b2>a21*upper+a22*lower && b2<(a21+a22)*upper && a22*b1-a12*b2>=detA2*upper){//case 6
                            out_x[i]  =upper;
                            out_x[i+1]=(b2-a21*upper)/a22;
                            FLAG = 6;//mexPrintf("case6\n");
                    }
                    else{//case 5
                            out_x[i  ]=(a22*b1-a12*b2)/detA2;
                            out_x[i+1]=(a11*b2-a21*b1)/detA2;
                            FLAG = 5;
                    }
                }
                // update the label of the block
                labels[i] = FLAG;
                // update temporal grad
                for (j=jcs[i];j<jcs[i+1];j++){
                    grad[irs[j]] += in_A[j]*out_x[i];
                }
                for (j=jcs[i+1];j<jcs[i+2];j++){
                    grad[irs[j]] += in_A[j]*out_x[i+1];
                }
            }
            // in the case if the dimension of A is odd
            if (in_d%2==1){
                i=in_d-1;
                // calc temporal grad
                for (j=jcs[i];j<jcs[i+1];j++){
                    grad[irs[j]] -= in_A[j]*out_x[i];
                }
                // descent
                out_x[i] = (in_b[i]-grad[i])/diag_A0[i];
                // bounds
                if (out_x[i]>upper){
                    out_x[i] = upper;
                }
                if (out_x[i]<lower){
                    out_x[i] = lower;
                }
                // update temporal grad
                for (j=jcs[i];j<jcs[i+1];j++){
                    grad[irs[j]] += in_A[j]*out_x[i];
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
            epoch++;
        }
    }
    /* if the bounds are defined, but lower>=upper
     * we take them as unconstrained
     * then do the following
     */
    else if (lower>=upper){
        while ((residual>in_precision)&&(epoch<in_max_iter)){
            for (i=0;i<in_d-1;i=i+2){
                // calc temporal grad
                // sparse g=g-A(:,i)*x(i) and i+1
                for (j=jcs[i];j<jcs[i+1];j++){
                    grad[irs[j]] -= in_A[j]*out_x[i];
                }
                for (j=jcs[i+1];j<jcs[i+2];j++){
                    grad[irs[j]] -= in_A[j]*out_x[i+1];
                }
                // update x(i)
                // define size 2 block
                a11=diag_A0[i];     a12=diag_A1[i  ];
                a21=diag_A1[i];     a22=diag_A0[i+1];
                b1 =in_b[i]  -grad[i];
                b2 =in_b[i+1]-grad[i+1];
                detA2 = a11*a22-a12*a21;
                if(detA2==0){mexErrMsgTxt("Input Matrix is not positive definite");  return;}
                // solve for linear system 2*2
                out_x[i  ]=(a22*b1-a12*b2)/detA2;
                out_x[i+1]=(a11*b2-a21*b1)/detA2;
                // update temporal grad
                for (j=jcs[i];j<jcs[i+1];j++){
                    grad[irs[j]] += in_A[j]*out_x[i];
                }
                for (j=jcs[i+1];j<jcs[i+2];j++){
                    grad[irs[j]] += in_A[j]*out_x[i+1];
                }
            }
            // in the case if the dimension of A is odd
            if (in_d%2==1){
                i=in_d-1;
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
            // get the residual
            residual = 0;
            for (i=0;i<in_d;i++){
                // i th residual
                df = grad[i]-in_b[i];
                residual += df*df;
            }
            residual = sqrt(residual);
            out_r[epoch] = residual;
            epoch++;
        }
    }
    plhs[1] = mxCreateDoubleMatrix(epoch,1,mxREAL);
    double* r = mxGetPr(plhs[1]);if(r==NULL){mexErrMsgTxt("pointer r is null");  return;}
    for (i=0;i<epoch;i++){
        r[i]=out_r[i];
    }
    delete [] grad; delete [] out_r; delete [] diag_A0; delete [] diag_A1; delete [] labels;
    mexPrintf("epoch:%5d, residual=%.15f\nEnd of CBCD size 2.cpp\n",epoch-1,residual);
}