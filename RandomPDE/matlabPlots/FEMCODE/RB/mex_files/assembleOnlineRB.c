#include <math.h>
#include "mex.h"


#define isDoubleRealArray(P) (!mxIsComplex(P) && !mxIsSparse(P) && mxIsDouble(P))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  #define T  prhs[0]
  #define TF prhs[1]
  #define OM prhs[2]
  #define ON prhs[3]
  #define A  plhs[0]
  #define F  plhs[1]
  
  /* INPUTS
   T  - a real vector of length Q
   TF - a real vector of length QF
   OM - a real 4D matrix of size N x N x Q x Q
        we assume that it is symmetric: OM(n1,n2,q1,q2) = OM(n2,n1,q2,q1)
   ON - a real 3D matrix of size N x Q x QF
   
   OUTPUS
   A  - a real 2D matrix of size N x N
        A(n1,n2) = (sum q1 = 1:Q)(sum q2 = 1:Q) T(q1) * T(q2) *OM(n1,n2,q1,q2)
   F  - a real 2D matrix of size N x 1
        F(n, 1) = (sum q = 1:Q)(sum qf = 1:QF) T(q) * TF(qf) * ON(n,q,qf) 
   */
  
  if(!(nrhs==4)){
    mexErrMsgTxt("Exactly four parameter are needed.");
  }
  /* Check input types*/
  mwSize i;
  for (i=0; i<4; i++){
    if(!isDoubleRealArray(prhs[i])){
      mexErrMsgTxt("All four imputs need to be real double non-sparse arrays.");
    }
  }
  
  /* get input sizes */
  mwSize Q  = mxGetNumberOfElements(T);
  mwSize QF = mxGetNumberOfElements(TF);
  mwSize N  = mxGetM(OM);
  const mwSize *OMdim = mxGetDimensions(OM);
  mwSize OMdims = mxGetNumberOfDimensions(OM);
  const mwSize *ONdim = mxGetDimensions(ON);
  mwSize ONdims = mxGetNumberOfDimensions(ON);
 
  /* check input sizes */
  if (!(mxGetM(T)==Q   || mxGetN(T) == Q))   
    mexErrMsgTxt("T is not a vector");
  if (!(mxGetM(TF)==QF || mxGetN(TF) == QF)) 
    mexErrMsgTxt("TF is not a vector");
  if (!(OMdim[0]==N && OMdim[1]==N)) 
    mexErrMsgTxt("OM does NOT have correct dimensions");
  if (!((OMdims==4 && OMdim[2]==Q && OMdim[3]==Q) || (OMdims==2 && Q==1))) 
    mexErrMsgTxt("OM does not have correct dimensions");
  if (!(ONdim[0]==N && ONdim[1]==Q)) 
    mexErrMsgTxt("ON does NOT have correct dimensions");
  if (!((ONdims==3 && ONdim[2]==QF) || (ONdims==2 && QF==1))) 
    mexErrMsgTxt("ON does not have correct dimensions");

  /* Create outputs */
  A = mxCreateDoubleMatrix(N,N,mxREAL);
  F = mxCreateDoubleMatrix(N,1,mxREAL);
  
  /* create pointers to inputs and outputs*/
  double *Tp, *TFp, *OMp, *ONp, *Ap, *Fp;
  Ap = mxGetPr(A);
  Fp = mxGetPr(F);
  Tp = mxGetPr(T);
  TFp = mxGetPr(TF);
  OMp = mxGetPr(OM);
  ONp = mxGetPr(ON);
  
  /* 1. ASSEMBLE THE MATRIX */
  mwSize n1, n2, q1, q2, n, q, qf;
 
  for (q2=0; q2<Q; q2++){
    for (q1=0; q1<Q; q1++){
      for (n2=0; n2<N; n2++){
        for (n1=0; n1<=n2; n1++){
          Ap[n1+N*n2] += Tp[q2] * Tp[q1] * OMp[n1+N*n2+N*N*q1+N*N*Q*q2];
        }
      }
    }
  }
  
  /* the symmetric counterpart */
  for (n2=0; n2<N; n2++){
    for (n1=n2+1; n1<N; n1++){
      Ap[n1+N*n2] = Ap[n2+N*n1];
    }
  }
  
  /* 1. ASSEMBLE THE RHS */
  for (qf=0; qf<QF; qf++){
    for (q=0; q<Q; q++){
      for (n=0; n<N; n++){
        Fp[n] += Tp[q] * TFp[qf] * ONp[n+N*q+N*Q*qf];
      }
    }
  }
  return;
}