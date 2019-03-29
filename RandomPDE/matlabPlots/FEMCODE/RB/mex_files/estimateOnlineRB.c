#include <math.h>
#include "mex.h"


#define isDoubleRealArray(P) (!mxIsComplex(P) && !mxIsSparse(P) && mxIsDouble(P))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  #define T  prhs[0]
  #define TF prhs[1]
  #define S  prhs[2]
  #define OP prhs[3]
  #define OQ prhs[4]
  #define OR prhs[5]
  #define R  plhs[0]
  
  /* INPUTS
   T  - a real vector of length Q
   TF - a real vector of length QF
   S  - a real vector of length N
   OP - a real 2D matrix of size QF x QF
   OQ - a real 3D matrix of size N x Q x QF
   OR - a real 4D matrix of size N x N x Q x Q
        we assume that it is symmetric: OR(n1,n2,q1,q2) = OR(n2,n1,q2,q1)
   
   OUTPUS
   R  - a real number given by
        R = (sum qf1 = 1:QF)(sum qf2 = 1:QF) TF(q1) * TF(q2) * OP(q1, q2)
          + 2 * (sum qf = 1:QF)(sum q = 1:Q)(sum n = 1:N) TF(qf) * T(q) * S(n) * OQ(n,q,qf);
          + (sum q2 = 1:Q)(sum q1 = 1:Q)(sum n2 = 1:N)(sum n1 = 1:N) T(q1) * T(q2) * S(n1) * S(n2) * OR(n1,n2,q1,q2);
   */
  
  if(!(nrhs==6)){
    mexErrMsgTxt("Exactly six parameter are needed.");
  }
  /* Check input types*/
  mwSize i;
  for (i=0; i<6; i++){
    if(!isDoubleRealArray(prhs[i])){
      mexErrMsgTxt("All four imputs need to be real double non-sparse arrays.");
    }
  }
  
  /* get input sizes */
  mwSize Q  = mxGetNumberOfElements(T);
  mwSize QF = mxGetNumberOfElements(TF);
  mwSize N  = mxGetNumberOfElements(S);
  const mwSize *OPdim = mxGetDimensions(OP);
  mwSize OPdims = mxGetNumberOfDimensions(OP);
  const mwSize *OQdim = mxGetDimensions(OQ);
  mwSize OQdims = mxGetNumberOfDimensions(OQ);
  const mwSize *ORdim = mxGetDimensions(OR);
  mwSize ORdims = mxGetNumberOfDimensions(OR);
  
  /* check input sizes */
  if (!(mxGetM(T)==Q   || mxGetN(T) == Q))   
    mexErrMsgTxt("T is not a vector");
  if (!(mxGetM(TF)==QF || mxGetN(TF) == QF)) 
    mexErrMsgTxt("TF is not a vector");
  if (!(mxGetM(S)==N || mxGetN(S) == N))
    mexErrMsgTxt("S is not a vector");
  if (!(OPdim[0]==QF && OPdim[1]==QF && OPdims == 2))
    mexErrMsgTxt("OP does NOT have correct dimensions");
  if (!(OQdim[0]==N && OQdim[1]==Q))
    mexErrMsgTxt("OQ does NOT have correct dimensions");
  if (!((OQdims==3 && OQdim[2]==QF) || (OQdims==2 && QF==1)))
    mexErrMsgTxt("OQ does not have correct dimensions");
  if (!(ORdim[0]==N && ORdim[1]==N))
    mexErrMsgTxt("OR does NOT have correct dimensions");
  if (!((ORdims==4 && ORdim[2]==Q && ORdim[3]==Q) || (ORdims==2 && Q==1))) 
    mexErrMsgTxt("OR does not have correct dimensions");

  /* Create outputs */
  R = mxCreateDoubleMatrix(1,1,mxREAL);

    /* create pointers to inputs and outputs*/
  double *Tp, *TFp, *Sp, *OPp, *OQp, *ORp, *Rp;
  Rp = mxGetPr(R);
  Tp = mxGetPr(T);
  TFp = mxGetPr(TF);
  Sp  = mxGetPr(S);
  OPp = mxGetPr(OP);
  OQp = mxGetPr(OQ);
  ORp = mxGetPr(OR);
  
  /* ASSEMBLE THE RESIDUAL (F-F part) */
  mwSize n, n1, n2, q, q1, q2, qf, qf1, qf2;
  double resFF, resAF, resAA;
  
  resFF = 0;
  for (qf2=0; qf2<QF; qf2++){
    for (qf1=0; qf1<QF; qf1++){
      resFF += TFp[qf2] * TFp[qf1] * OPp[qf1 + QF*qf2];
    }
  }
  
  /* ASSEMBLE THE RESIDUAL (F-A part) */
  resAF = 0;
  for (qf=0; qf<QF; qf++){
    for (q=0; q<Q; q++){
      for (n=0; n<N; n++){
        resAF += TFp[qf] * Tp[q] * Sp[n] * OQp[n+q*N+qf*N*Q];        
      }
    }
  }
  
  /* ASSEMBLE THE RESIDUAL (A-A part) */
  resAA = 0;
  for (q2=0; q2<Q; q2++){
    for (q1=0; q1<Q; q1++){
      for (n2=0; n2<N; n2++){
        for (n1=0; n1<n2; n1++){
          resAA += 2 * Tp[q2] * Tp[q1] * Sp[n2] * Sp[n1] * ORp[n1+N*n2+N*N*q1+N*N*Q*q2];
        }
        resAA += Tp[q2] * Tp[q1] * Sp[n2] * Sp[n2] * ORp[n2+N*n2+N*N*q1+N*N*Q*q2];
      }
    }
  }
  Rp[0] = resFF + 2*resAF + resAA;
  return;
}