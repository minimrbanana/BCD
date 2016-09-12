#include "mex.h"
/* WELL def part begin*/
#define W 32
#define R 16
#define P 0
#define M1 13
#define M2 9
#define M3 5

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT3NEG(t,v) (v<<(-(t)))
#define MAT4NEG(t,b,v) (v ^ ((v<<(-(t))) & b))

#define V0            STATE[state_i                   ]
#define VM1           STATE[(state_i+M1) & 0x0000000fU]
#define VM2           STATE[(state_i+M2) & 0x0000000fU]
#define VM3           STATE[(state_i+M3) & 0x0000000fU]
#define VRm1          STATE[(state_i+15) & 0x0000000fU]
#define VRm2          STATE[(state_i+14) & 0x0000000fU]
#define newV0         STATE[(state_i+15) & 0x0000000fU]
#define newV1         STATE[state_i                 ]
#define newVRm1       STATE[(state_i+14) & 0x0000000fU]

#define FACT 2.32830643653869628906e-10
// WELL parameters
static unsigned int state_i = 0;
static unsigned int STATE[R];
static unsigned int z0, z1, z2;
/* WELL def part end*/
double fval_mex(double *A,double *b,int d,double *x);//function value
double residual_mex(double *A,double *b,int d,double *x);//residual
void InitWELLRNG512a (unsigned int *init);//WELL
double WELLRNG512a (void);//WELL
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
   int i,j,k,epoch;//loop
   int lower_index,upper_index;//search random index in L_sum
   double Aix,residual,L_random;//matrix calc,residual,random L
   double *L_sum;
   // WELL parameters
   unsigned int well_init[16];
   for(i=0;i<16;i++){
       well_init[i]=i+1;
   }
   InitWELLRNG512a(well_init);
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
   plhs[1] = mxCreateDoubleMatrix(in_max_iter+1,1,mxREAL);
   out_y = mxGetPr(plhs[1]);
   out_y[0] = fval_mex(in_A, in_b, in_d, out_x);
   for (i=0;i<=in_max_iter;i++){
       out_y[i] = 0;
   }
   // calc L_sum
   plhs[2] = mxCreateDoubleMatrix(in_d+1,1,mxREAL);
   L_sum = mxGetPr(plhs[2]);
   L_sum[0] = 0;
   for (i=0;i<in_d;i++){
       L_sum[i+1] = in_A[i*in_d + i]+L_sum[i];
   }
   for (i=0;i<=in_d;i++){
       L_sum[i] = L_sum[i]/L_sum[in_d];
   }
   /*mexPrintf("L_sum:%.4f,",L_sum[0]);
   for (i=1;i<in_d-1;i++){
       mexPrintf("  %.4f,",L_sum[i]);
   }
   mexPrintf("  %.4f\n",L_sum[in_d]);*/
   // residual of init
   residual = residual_mex(in_A, in_b, in_d, out_x);
   mexPrintf("epoch:    0, residual=%.15f, fval:%.8f\n",residual,out_y[0]);
   epoch=1;
   while ((residual!=0)&&(epoch<=in_max_iter)){
       for (k=0;k<in_d;k++){// one epoch
           //get random number in [0,1]
           L_random = WELLRNG512a();
           // get corresponding i
           lower_index = 0;
           upper_index = in_d;
           while(upper_index-lower_index>1){
               i = floor((upper_index+lower_index)/2);
               if (L_random>L_sum[i]){
                   lower_index = i;
               }
               if (L_random<=L_sum[i]){
                   upper_index = i;
               }
           }
           i = lower_index;
           //mexPrintf("i:%d\n",i);//show i when debugging
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
       }
       //residual
       residual = residual_mex(in_A,in_b,in_d,out_x);
       out_y[epoch] = fval_mex(in_A, in_b, in_d, out_x);
       mexPrintf("epoch:%5d, residual=%.15f, fval:%.8f\n",epoch,residual,out_y[epoch]);
       epoch++;
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
double residual_mex(double *A,double *b,int d,double *x){
    double r;//output residual, scalar
    int i,j;
    double df;
    r = 0;
    for (i=0;i<d;i++){
        df = -b[i];
        for (j=0;j<d;j++){
            df = df + A[i*d + j]*x[j];
        }
        if (x[i]==0){
            if (df<0){
                r = r + df*df;
            }
        }
        else if (x[i]==1){    
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
void InitWELLRNG512a (unsigned int *init){
   int j;
   state_i = 0;
   for (j = 0; j < R; j++)
     STATE[j] = init[j];
}
double WELLRNG512a (void){
  z0    = VRm1;
  z1    = MAT0NEG (-16,V0)    ^ MAT0NEG (-15, VM1);
  z2    = MAT0POS (11, VM2)  ;
  newV1 = z1                  ^ z2; 
  newV0 = MAT0NEG (-2,z0)     ^ MAT0NEG(-18,z1)    ^ MAT3NEG(-28,z2) ^ MAT4NEG(-5,0xda442d24U,newV1) ;
  state_i = (state_i + 15) & 0x0000000fU;
  return ((double) STATE[state_i]) * FACT;
}