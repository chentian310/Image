## fastperm3.R version
fastperm3.R <- function(X,t.min,t.max){
   N<-length(X);B<-N/3;ind<-seq(0,N-1,3)
   XX <- matrix(X,ncol=3,byrow=T)
   Xtilde <- matrix(0,nrow=N,ncol=3)
   for (l in 1:3){
     Xtilde[,l] <- as.vector(kermat[,ind1+l]%*%XX[,l])
   }
   Xt0 <- Xtilde[,1];Xt1 <- Xtilde[,2]; Xt2 <- Xtilde[,3]
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
   return(list(t.max,t.min))
 }


## fastperm3.C version
fastperm3.C <- function(X,t.mim,t.max){
   N<-length(X);B<-B/3;ind<-seq(0,N-1,3)
   XX <- matrix(X,ncol=3,byrow=T)
   Xtilde <- matrix(0,nrow=N,ncol=3)
   for (l in 1:3){
     Xtilde[,l] <- as.vector(kermat[,ind1+l]%*%XX[,l])
   }
   Xt0 <- Xtilde[,1];Xt1 <- Xtilde[,2]; Xt2 <- Xtilde[,3]
   ## implemented in C
    XtildePerm <- rep(0,150)
   .C("fastperm3",as.double(Xt0),as.double(Xt1),as.double(Xt2),
      as.integer(B),as.integer(N),as.double(t.max),as.double(t.min),PACKAGE="flash")
   return(list(t.min,t.max))
 }

