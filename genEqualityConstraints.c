#include "mex.h"

void genEqualityConstraints(double *Aeq, size_t *ir, size_t *jc, int n,int numVars){
    int t,j,v1,v2,entry,rowIndex,k;
    k = 0;
    for (t=1;t<n+1;++t){
        for (j=t+1;j<n+1;++j){
            if (j < t){
                
                v1 = j;
                
                v2 = t;
            }
            if (j > t){
                v1 = t;
                v2=j;
            }
            if (j != t){
                rowIndex = t-1;
                entry = n*(v1-1)+v2-((v1*(v1+1))/2)-1;
                Aeq[k]=1;
                ir[k]=t-1;
                Aeq[k+1]=1;
                ir[k+1]=j-1;
                jc[entry]=k;
                k=k+2;
            }
        }
    }
    jc[(n*(n-1))/2]=k;
}
    
    

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *out;
    double *in;

    double *Aeq;
    
    size_t *ir;
    size_t *jc;
    int n;
    int numVars;
    
    in = mxGetPr(prhs[0]);
    n = in[0];
    numVars = (n*(n-1))/2;
    out = plhs[0] = mxCreateSparse(n,numVars,n*(n-1),mxREAL);
    Aeq = mxGetPr(out);
    ir = mxGetIr(out);
    jc = mxGetJc(out);
    genEqualityConstraints(Aeq,ir,jc,n,numVars);
}
