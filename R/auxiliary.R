## fastperm3.R version
fastperm3.R <- function(X){
  N<-length(X);B<-N/3;ind<-seq(0,N-1,3)
  XX <- matrix(X,ncol=3,byrow=T)
  Xtilde <- matrix(0,nrow=N,ncol=3)
  for (l in 1:3){
    Xtilde[,l] <- as.vector(kermat[,ind1+l]%*%XX[,l])
  }
  Xt0 <- Xtilde[,1];Xt1 <- Xtilde[,2]; Xt2 <- Xtilde[,3] # should be changed if L is not 3
  ## R version
  t.min <- rep(0,B^2);t.max <- rep(0,B^2)
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
fastperm3.C <- function(X,L=3){ # t.min, t.max seems not necessary here
  N<-length(X);B<-N/L;ind<-seq(0,N-1,3)
  XX <- matrix(X,ncol=3,byrow=T)
  Xtilde <- matrix(0,nrow=N,ncol=3)
  for (l in 1:3){
    Xtilde[,l] <- as.vector(kermat[,ind1+l]%*%XX[,l])     # uses "kermat",better as a input
  }
  Xt0 <- Xtilde[,1];Xt1 <- Xtilde[,2]; Xt2 <- Xtilde[,3]  # need to be more flexible
  ## implemented in C
  Xtildeperm <- rep(0,N)
  t.min <- rep(0,B^2);t.max <- rep(0,B^2)
  t.min=.C("fastperm3",as.double(Xt0),as.double(Xt1),as.double(Xt2),
     as.double(Xtildeperm), as.integer(B),as.integer(N),as.integer(L),
     maxvec=as.double(t.max),minvec=as.double(t.min),PACKAGE="flash")$minvec
  t.max=.C("fastperm3",as.double(Xt0),as.double(Xt1),as.double(Xt2),
     as.double(Xtildeperm), as.integer(B),as.integer(N),as.integer(L),
     maxvec=as.double(t.max),minvec=as.double(t.min),PACKAGE="flash")$maxvec
  return(list(t.min,t.max))
}


## The Gauss Filter Matrix
GaussFilterMat <- function(N,sigma=5,cut=4){
   sigma <- 5.0; cutoff <- min(cut*sigma,N-1);tt<-0:(N-1)
   kern1 <- rep(0,N)
   kern1[1:(cutoff+1)] <- dnorm(0:cutoff,sd=sigma)
   kern1[(N-cutoff+1):N] <- kern1[(N-cutoff+1):N]+dnorm(-cutoff:-1,sd=sigma)
   kermat <- sapply(tt,shift,x=kern1)
   return(kermat)
 }


## A shift function for permutation
shift <-  function(x,i=1){
      n<-length(x);i<-i%%n
    if (i==0){
         return(x)
       }
    return(x[c((n-i+1):n,1:(n-i))])
  }
