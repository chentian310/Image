## fastperm3.R version
fastperm3.R <- function(X,bandwidth){
  N<-length(X);B<-N/3;ind1<-seq(0,N-1,3)
  XX <- matrix(X,ncol=3,byrow=T)
  Xtilde <- matrix(0,nrow=N,ncol=3)
  sigma <- bandwidth*0.3706506
  kermat <- GaussFilterMat(N,sigma,cut=4) 
  
  for (l in 1:3){
    Xtilde[,l] <- as.vector(kermat[,ind1+l]%*%XX[,l])
  }
  Xt0 <- Xtilde[,1];Xt1 <- Xtilde[,2]; Xt2 <- Xtilde[,3]
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
fastperm3.C <- function(X,bandwidth){ 
  N<-length(X);B<-N/3;ind1<-seq(0,N-1,3);L=3
  XX <- matrix(X,ncol=3,byrow=T)
  Xtilde <- matrix(0,nrow=N,ncol=3)
  sigma <- bandwidth*0.3706506
  kermat <- GaussFilterMat(N,sigma,cut=4)
  
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
GaussFilterMat <- function(N,sigma,cut=4){
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

## Find the rank of x in a vector, using C implementation
rev.rank <- function(x, x.perm) {
  y <- sort(x.perm, decreasing=TRUE)
  o <- order(x, decreasing=TRUE); ro <- order(o)
  rv <- rep(-1,length(x))
  z <- .C("vec_rev_rank", as.double(x[o]), as.double(y),
          as.integer(length(x)), as.integer(length(y)),
          rv=as.integer(rv), PACKAGE = "flash")
  return(z$rv[ro])
}

## Simple version P-value
## P-VALUE: each pixel has 150 time points, each point has a pvalue
## so the p-value matrix is 30000+ * 150

## Input the smoothed temporal matrix(each row indicates a pixel)
## tmin: each row for a pixel, 2500 perm->2500 tmins
## tmax: each row for a pixel, 2500 perm->2500 tmaxs
P.adj <- function(sm.mat,tmin,tmax){
  nrow <- dim(sm.mat)[1];ncol <- dim(sm.mat)[2]
  p <- matrix(0,nrow,ncol);p.adj <- matrix(0,nrow,ncol)
  for(i in 1:nrow){
    rank.min <- rev.rank(sm.mat[i,],tmin[i,])
    rank.max <- rev.rank(sm.mat[i,],tmax[i,])
    diff <- p.max-p.min
    p[i,] <- 2*(rank.min*(diff<0)+rank.max*(diff>0))/ncol
    p.adj[i,]<-p.adjust(p[i,],method="BH")
  }
  return(p.adj)
}















## SCB p-value
##    sm.mat: matrix, each row is the smoothed temporal data for a pixel
##      tmin: the min value in each of B^2 simulation 
##      tmax: the max value in each of B^2 simulation
#p.scb <- function(sm.mat,tmin,tmax){
#  n=length(sm.mat)
#  sm.vec <- as.vector(sm.mat)
#  tmin.vec <- as.vector(tmin);tmax.vec <- as.vector(tmax)
#  p.min <- sapply(sm.vec,function(x) sum(tmax.vec >=  x)/n)
#  p.max <- sapply(sm.vec,function(x) sum(tmin.vec <=  x)/n)
#  pmatrix <- matrix(c(p.min,p.max),ncol=2,byrow=F)
#  p <- apply(pmatrix,1,function(x) 2*min(x[1],x[2]))
#  p.adj <- p.adjust(p,method="BH")
#  return(p.adj)
#}

