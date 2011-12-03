#include <stdio.h>
#include <math.h>
/*#include <Rmath.h>*/
#include <R.h>

void fastperm3(double *X0,double *X1,double *X2,double *XtildePerm,int *B,int *N,int *L,
               double *maxvec, double *minvec){
  int i=0,j=0,m1=0,m2=0,BB=*B,NN=*N,ind=0,LL=*L;
  double tmax=0,tmin=0;
  for (m1=0;m1<BB;m1++){
    for (i=0;i<NN;i++) {
      ind=(i+LL*m1) % NN;
      ind=((ind == 0) ? NN:ind);
      XtildePerm[i]=X0[i]+X1[ind];
    }
    for (m2=0;m2<BB;m2++){
      for (j=0;j<NN;j++){
        ind=(j+LL*m2) % NN;
        ind=((ind == 0) ? NN:ind);
        XtildePerm[j]+=X2[ind];
        tmin=((XtildePerm[j] < tmin) ? XtildePerm[j]:tmin);
        tmax=((XtildePerm[j] > tmax) ? XtildePerm[j]:tmax);
      }
      maxvec[m1*BB+m2]=tmax;
      minvec[m1*BB+m2]=tmin;
    }
  }
}
