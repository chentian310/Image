## fastperm3.R version
fastperm3.R <- function(X,t.min,t.max,L=3){
  N<-length(X);B<-N/L;ind<-seq(0,N-1,L)
  XX <- matrix(X,ncol=L,byrow=T)
  Xtilde <- matrix(0,nrow=N,ncol=L)
  for (l in 1:L){
    Xtilde[,l] <- as.vector(kermat[,ind1+l]%*%XX[,l])
  }
  Xt0 <- Xtilde[,1];Xt1 <- Xtilde[,2]; Xt2 <- Xtilde[,3] # should change if L is not 3
  ## R version
  for (m1 in 1:B){
    Xt1 <- shift(Xt1,3)
    sums01 <- Xt0 + Xt1
    for (m2 in 1:B){
      Xt2 <- shift(Xt2,3)
      Xtilde.perm <- sums01 + Xt2
      perm.ind <- (m1-1)*B+m2
      t.max[perm.ind] <- max(Xtilde.perm)
      t.min[perm.ind] <- min(Xtilde.perm)
    }
  }
  return(list(t.min,t.max))
}


## fastperm3.C version
fastperm3.C <- function(X,t.min,t.max,L=3){
  N<-length(X);B<-N/3;ind<-seq(0,N-1,3)
  XX <- matrix(X,ncol=3,byrow=T)
  Xtilde <- matrix(0,nrow=N,ncol=3)
  for (l in 1:3){
    Xtilde[,l] <- as.vector(kermat[,ind1+l]%*%XX[,l])
  }
  Xt0 <- Xtilde[,1];Xt1 <- Xtilde[,2]; Xt2 <- Xtilde[,3]
  ## implemented in C
  XtildePerm <- rep(0,N)
#  t.mins=.C("fastperm3",as.double(Xt0),as.double(Xt1),as.double(Xt2),
#     as.double(XtildePerm), as.integer(B),as.integer(N),
#     as.double(t.max),as.double(t.min),PACKAGE="flash")$t.min
#  t.maxs=.C("fastperm3",as.double(Xt0),as.double(Xt1),as.double(Xt2),
#     as.double(XtildePerm), as.integer(B),as.integer(N),
#     as.double(t.max),as.double(t.min),PACKAGE="flash")$t.max
   MaxVec <- .C("fastperm3",as.double(Xt0),as.double(Xt1),as.double(Xt2),
     XtildePerm=as.double(XtildePerm), as.integer(B),as.integer(N),as.integer(L),
     maxvec=as.double(t.max),minvec=as.double(t.min),PACKAGE="flash")$maxvec
   Xtildeperm <- .C("fastperm3",as.double(Xt0),as.double(Xt1),as.double(Xt2),
     XtildePerm=as.double(XtildePerm), as.integer(B),as.integer(N),as.integer(L),
     maxvec=as.double(t.max),minvec=as.double(t.min),PACKAGE="flash")$XtildePerm
   MinVec <- .C("fastperm3",as.double(Xt0),as.double(Xt1),as.double(Xt2),
     XtildePerm=as.double(XtildePerm), as.integer(B),as.integer(N),as.integer(L),
     maxvec=as.double(t.max),minvec=as.double(t.min),PACKAGE="flash")$minvec
  return(list(MinVec,MaxVec,Xt0,Xt1,Xt2,Xtildeperm))
}
