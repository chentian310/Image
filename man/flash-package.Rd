\name{flash-package}
\alias{flash-package}
\alias{flash}
\docType{package}
\title{
Preprocess pictures and then do inference to detect flash area
}
\description{
Haven't done yet. Wait for updated
}
\details{
\tabular{ll}{
Package: \tab flash\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-11-27\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\seealso{\code{\link{ksmooth}}, \code{\link{fft}}}
\examples{
## This part must be written in valid R syntax!

## generate sample
N <- 150;  tt<-0:(N-1); flash.region <- floor(N/3):ceiling(2*N/3)
N.flash <- length(flash.region)
ff <-rep(0,N);ff[flash.region]<-sin(1:N.flash/N.flash*pi)
X<-0.3*rnorm(N)+ff

N=length(X);L=3;B<-N/L;ind1<-seq(0,N-1,3)

## A shift function
shift <- function(x,i=1){
    n<-length(x);i<-i%%n
    if (i==0){
         return(x)
       }
    return(x[c((n-i+1):n,1:(n-i))])
  }

## The Gauss filter
sigma <- 5.0; cutoff <- min(4*sigma,N-1)
kern1 <- rep(0,N)
kern1[1:(cutoff+1)] <- dnorm(0:cutoff,sd=sigma)
kern1[(N-cutoff+1):N] <- kern1[(N-cutoff+1):N]+dnorm(-cutoff:-1,sd=sigma)
kermat <- sapply(tt,shift,x=kern1)  ## must-have

## fastperm3.C workable ??
t.max <- rep(0,50^2); t.min <- rep(0,50^2)
fastperm3.C(X,t.min,t.max)

## Test the validation of fastperm3.C
## should obtain the same result as fastper3.R
t.max <- rep(0,50^2); t.min <- rep(0,50^2)
perm.R.min=matrix(0,nrow=100,ncol=2500);perm.R.max=matrix(0,nrow=100,ncol=2500)
for (i in 1:100) {
  container=fastperm3.R(X,t.min,t.max)
  perm.R.min[i,]=container[[1]]
  perm.R.max[i,]=container[[2]]
}

}