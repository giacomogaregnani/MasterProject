#include <math.h>
#include "mex.h"


#define isDoubleRealArray(P) (!mxIsComplex(P) && !mxIsSparse(P) && mxIsDouble(P))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  #define T  prhs[0]
  #define TF prhs[1]
  #define S  prhs[2]
  #define OU prhs[3]
  #define OS prhs[4]
  #define OT prhs[5]
  #define tensor plhs[0]
  
  /* INPUTS
   T  - a real vector of length Q
   TF - a real vector of length QF
   S  - a cell array of dimension d,
        S{i} contains a real vector of length N(i)
   OU - a cell array of size d x d 
        OU{i,j} contains a real 2D matrix of size QF x N(j)
   OS - a cell array of size d x d 
        OS{i,j} contains a real 2D matrix of size N(i) x QF
   OT - a cell array of size d x d 
        OT{i,j} contains a real 3D matrix of size N(i) x N(j) x Q
   
   OUTPUS
   tensor - a real 2D array of size d x d given by
   tensor(i,j) = (sum qf = 1:QF)(sum nj = 1:N(j)) TF(qf) * S{j}(nj) * OU{i,j}(qf, nj)
          + (sum qf = 1:QF)(sum ni = 1:N(i)) TF(qf) * S{i}(ni) * OS{i,j}(ni,qf);
          + (sum q = 1:Q)(sum ni = 1:N(i))(sum nj = 1:N(j)) T(q) * S{i}(ni) * S{j}(nj) * OT{i,j}(ni,nj,q);
   */
  
  if(!(nrhs==6)){
    mexErrMsgTxt("Exactly six parameter are needed.");
  }
  /* Check first two inputs*/
  if(!(isDoubleRealArray(T) && isDoubleRealArray(TF))){
      mexErrMsgTxt("First two inputs need to be real double non-sparse arrays.");
    }
  mwSize Q  = mxGetNumberOfElements(T);
  mwSize QF = mxGetNumberOfElements(TF);
  if (!(mxGetM(T)==Q   || mxGetN(T) == Q))
    mexErrMsgTxt("T is not a vector");
  if (!(mxGetM(TF)==QF || mxGetN(TF) == QF))
    mexErrMsgTxt("TF is not a vector");
  
  /* check if other inputs are cells */
  mwSize i, j, ind;
  for (i=2; i<6; i++){
    if(!mxIsCell(prhs[i])){
      mexErrMsgTxt("Parameters 2-6 need to be cell arrays");
    }
  }
  /* and if the cell arrays have right dimensions */
  mwSize d = mxGetM(S);
  if (!(mxGetNumberOfElements(S)==d)) {
    mexErrMsgTxt("Third parameter is not a cell with the right dimensions");
  }
  for (i=3; i<6; i++){
    if(!(mxGetM(prhs[i])==d && mxGetN(prhs[i])==d && mxGetNumberOfElements(prhs[i])==d*d)) {
      mexErrMsgTxt("Parameters 2-6 need to be cell arrays of size d times d");
    }
  }

  /* pointers to structures */
  const mxArray *SE[d], *OUE[d][d], *OSE[d][d], *OTE[d][d]; 
  double *Tp, *TFp, *Sp[d], *OUp[d][d], *OSp[d][d], *OTp[d][d], *tensorp;
  for (i=0; i<d; i++){
    SE[i] = mxGetCell(S, i);
    Sp[i] = mxGetPr(SE[i]);
    for (j=0; j<d; j++){
      ind = i+d*j;
      OUE[i][j]  = mxGetCell(OU, ind);
      OSE[i][j]  = mxGetCell(OS, ind);
      OTE[i][j]  = mxGetCell(OT, ind);
      OUp[i][j]  = mxGetPr(OUE[i][j]);
      OSp[i][j]  = mxGetPr(OSE[i][j]);
      OTp[i][j]  = mxGetPr(OTE[i][j]);
    }
  }
  Tp = mxGetPr(T);
  TFp = mxGetPr(TF);
  tensor = mxCreateDoubleMatrix(d,d,mxREAL);
  tensorp = mxGetPr(tensor);
  
  /* now find the values N(i) */
  mwSize N[d];
  for (i=0; i<d; i++) {
    N[i] = mxGetNumberOfElements(SE[i]);
  }
  

  /* check the sizes of OU, OS, OT */  
  for (i=0; i<d; i++) {
    for (j=0; j<d; j++) {
      if (!(mxGetM(OUE[i][j])==QF && mxGetN(OUE[i][j])==N[j] && mxGetNumberOfDimensions(OUE[i][j])==2)) {
        mexErrMsgTxt("Parameter 4 does not have the right size");
      }
      if (!(mxGetM(OSE[i][j])==N[i] && mxGetN(OSE[i][j])==QF && mxGetNumberOfDimensions(OSE[i][j])==2)) {
        mexErrMsgTxt("Parameter 5 does not have the right size");
      }        
      if (!(mxGetM(OTE[i][j])==N[i] && mxGetN(OTE[i][j])==N[j]*Q && \
              mxGetNumberOfDimensions(OTE[i][j])<=3 && \
              mxGetNumberOfElements(OTE[i][j]) == Q*N[i]*N[j])) {
        mexErrMsgTxt("Parameter 6 does not have the right size");
      } 
    }
  }
  
  double resU, resS, resT;
  mwSize ni, nj, n2, q, qf;
  /* compute the desired value(s) */
  for (i=0; i<d; i++) {
    for (j=0; j<d; j++) {
      ind = i+j*d;
      
      /* OU part */
      resU = 0;
      for (nj=0; nj<N[j]; nj++){
        for (qf=0; qf<QF; qf++){
          resU += TFp[qf] * Sp[j][nj] * OUp[i][j][qf+QF*nj];
        }
      }
      
      /* OS part */
      resS = 0;      
      for (qf=0; qf<QF; qf++){
        for (ni=0; ni<N[i]; ni++){
          resS += TFp[qf] * Sp[i][ni] * OSp[i][j][ni+N[i]*qf];
        }
      }
      
      /* OT part */
      resT = 0;
      for (q=0; q<Q; q++){
        for (nj=0; nj<N[j]; nj++){
          for (ni=0; ni<N[i]; ni++){
            resT += Tp[q] * Sp[j][nj] * Sp[i][ni] * OTp[i][j][ni+nj*N[i]+q*N[i]*N[j]];
          }
        }
      }
      
      /* put together */
      tensorp[ind] = resU + resS + resT;      
    }
  }  
  return;
}