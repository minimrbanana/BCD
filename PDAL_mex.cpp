#include <math.h>
#include "mex.h"
/*
%% PDAL from paper
% 'Malitsky, Pock - A first-order primal-dual algorithm with
% linesearch(2016)'
% implementation of algorithm 3
% input          A, b 
% output         x*, y(fval)
%%
*/
#define EPSILON 2.220446e-16
double residual_mex(double *A,double *b,int d,double *x,double *grad);
double *grad_mex(double *A,double *x,int d, double *grad);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 5 || nlhs != 3) {
        mexWarnMsgTxt("Usage: X=CoordinateDescentQPBox(B,M)");
        return;
    }
    mexPrintf("To compare, the input of corresponding CBCD functions should be A^T*A and A^T*b\n");
    //input args
    //here the lower and upper bounds are set in the program
    //not set by the input parameters
    double *in_A;
    double *in_b;
    double *A2;//A2 = A'*A
    double *b2;//b2 = A'^b
    double *p_tau;//tau  = 1/sqrt(max(eig(A'*A)));input from matlab calculation
    //ouput args
    double *out_x; // the minimizer
    double *out_y; // dual
    double *n_loop;// number of loops(epoch)
    //get inout args
    in_A = mxGetPr(prhs[0]);if(in_A==NULL){mexErrMsgTxt("pointer in_A is null");  return;}
    in_b = mxGetPr(prhs[1]);if(in_b==NULL){mexErrMsgTxt("pointer in_b is null");  return;}
    A2   = mxGetPr(prhs[2]);if(A2  ==NULL){mexErrMsgTxt("pointer A2 is null");  return;}
    b2   = mxGetPr(prhs[3]);if(b2  ==NULL){mexErrMsgTxt("pointer b2 is null");  return;}
    p_tau= mxGetPr(prhs[4]);if(p_tau ==NULL){mexErrMsgTxt("pointer p_tau is null");  return;}
    int dim = (int)mxGetM(prhs[1]);//get dimension
    mexPrintf("CBCD size 1...input args get\n");
    //allocate output, and init
    plhs[0] = mxCreateDoubleMatrix(dim,1,mxREAL);
    out_x = mxGetPr(plhs[0]);if(out_x==NULL){mexErrMsgTxt("pointer out_x is null");  return;} 
    for (int i=0;i<dim;i++){
        out_x[i] = 0;
    }
    plhs[1] = mxCreateDoubleMatrix(dim,1,mxREAL);
    out_y = mxGetPr(plhs[1]);if(out_y==NULL){mexErrMsgTxt("pointer out_y is null");  return;} 
    for (int i=0;i<dim;i++){
        out_y[i] = -in_b[i];
    }
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    n_loop = mxGetPr(plhs[2]);if(n_loop==NULL){mexErrMsgTxt("pointer n_loop is null");  return;} 
    n_loop[0] = 0;
    //parameters in the function
    int i,j,epoch=1,count;//loop
    double A_y,A_x_bar,break_line,norm_d_y;
    double miu=0.5,tau=*p_tau,tau_old,beta=1,theta=1,gamma=0.5,sigma,residual=1;
    double *x_old = new double[dim];if(x_old==NULL){mexErrMsgTxt("pointer x_old is null");  return;}
    double *x_bar = new double[dim];if(x_bar==NULL){mexErrMsgTxt("pointer x_bar is null");  return;} 
    double *y_old = new double[dim];if(y_old==NULL){mexErrMsgTxt("pointer y_old is null");  return;}
    double *grad  = new double[dim];if(grad ==NULL){mexErrMsgTxt("pointer grad is null");  return;} 
    //main iteration
    while ((residual>1E-12) && (epoch<=1000)){
        //1. compute
        for (i=0;i<dim;i++){
            x_old[i]=out_x[i];
        }
        for (i=0;i<dim;i++){
            A_y=0;
            for (j=0;j<dim;j++){
                A_y+=in_A[j*dim+i]*out_y[j];
            }
            out_x[i]=out_x[i]-tau*A_y;
            if (out_x[i]>1){out_x[i]=1;}
            else if (out_x[i]<0){out_x[i]=0;}
        }
        beta = beta/(1+gamma*beta*tau);
        //2. line search
        tau_old=tau;
        tau = tau*((sqrt(1+theta)-1)*((double)rand())/((double)RAND_MAX + 1.0)+1);
        //2.a compute
        theta = tau/tau_old;
        sigma = beta*tau;
        for (i=0;i<dim;i++){
            x_bar[i] = out_x[i] + theta*(out_x[i]-x_old[i]);
            y_old[i] = out_y[i];
        }
        for (i=0;i<dim;i++){
            A_x_bar=0;
            for (j=0;j<dim;j++){
                A_x_bar+=in_A[j*dim+i]*x_bar[j];
            }
            out_y[i] = (out_y[i]+sigma*(A_x_bar-in_b[i]))/(1+sigma);
        }
        //2.b break line if
        break_line = 0;norm_d_y=0;
        for (i=0;i<dim;i++){
            // calculate norm(A*(ycur-yold))
            A_y=0;
            for (j=0;j<dim;j++){
                A_y+=in_A[j*dim+i]*(out_y[j]-y_old[j]);
            }
            break_line += A_y*A_y;
            // calculate norm(ycur-yold)
            norm_d_y += (out_y[i]-y_old[i])*(out_y[i]-y_old[i]);
        }
        break_line = sqrt(break_line)*sqrt(beta)*tau/sqrt(norm_d_y);
        count=1;
        while ((break_line>1) && (count<=50)){
            //2.a compute
            tau = tau*miu;
            theta = tau/tau_old;
            sigma = beta*tau;
            for (i=0;i<dim;i++){
                x_bar[i] = out_x[i]+theta*(out_x[i]-x_old[i]);
                y_old[i] = out_y[i];
            }
            for (i=0;i<dim;i++){
                A_x_bar=0;
                for (j=0;j<dim;j++){
                    A_x_bar+=in_A[j*dim+i]*x_bar[j];
                }
                out_y[i] = (out_y[i]+sigma*(A_x_bar-in_b[i]))/(1+sigma);
            }
            //2.b break line
            break_line = 0;norm_d_y=0;
            for (i=0;i<dim;i++){
                // calculate norm(A*(ycur-yold))
                A_y=0;
                for (j=0;j<dim;j++){
                    A_y+=in_A[j*dim+i]*(out_y[j]-y_old[j]);
                }
                break_line += A_y*A_y;
                // calculate norm(ycur-yold)
                norm_d_y += (out_y[i]-y_old[i])*(out_y[i]-y_old[i]);
            }
            break_line = sqrt(break_line)*sqrt(beta)*tau/sqrt(norm_d_y);
            count++;            
        }
        //grad and residual
        grad = grad_mex(A2,out_x,dim,grad);
        residual = residual_mex(A2, b2, dim, out_x, grad);
        //mexPrintf("epoch:%5d, residual=%.15f\n",epoch,residual);
        epoch++;
    }
    delete x_old; delete y_old; delete x_bar;delete grad;
    mexPrintf("End of PDAL\n");

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