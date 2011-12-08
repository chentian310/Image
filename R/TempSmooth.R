# Temporal Smooth:

# TempSmooth

## 
## mask = nonmito.ind(TRUE-nonmito)
## preproc = 3D array after clustering, normlization and spacial smoothing 
##
##

## think about: for mask, need one more option? nonmito=T, if T,then T=nonmito,else T=mito



TempSmooth <- function(preproc,mask,bandwidth,method=c("interval",default="circular"),...){ 
    Timenames<-dimnames(preproc)[[3]]
    nrows <- dim(preproc)[1];ncols <- dim(preproc)[2];TimeLen <- dim(preproc)[3];Time <- 1:TimeLen
    mitolen <- sum(!mask)                                # number of mito pixel
    mito.loc <- which(!mask); mito.vec<-preproc[!mask]   # numeric: all mito infor contained 
    pixel.temp <- array(0,dim=c(mitolen,TimeLen))        # Each row is a temporal data for a pixel in mito area
    Temp.sm <- array(0,dim=dim(preproc))                 # temproal smoothed data container,3D
    sm.mat <- array(0,dim=c(mitolen,TimeLen))            # matrix, each row is the smoothed temporal data.  
    all.vec=rep(0,nrows*ncols*TimeLen)
    
##   obtain temporal data for each mito pixel
     for (i in 1:mitolen) {
       pixel.temp[i,]<-mito.vec[i+mitolen*(0:(TimeLen-1))]    # temproal data for each mito pixel
  }


 if (method=="interval") {
    for (i in 1:mitolen) {
      sm.mat[i,]<-ksmooth(Time,pixel.temp[i,],kernel="normal",bandwidth=bandwidth,x.points=Time,...)$y
      k=mito.loc[i] 
## smoothed data put back in the full vec (length=256*256*150) 
      all.vec[length(mask)*(0:(TimeLen-1))+k] <- sm.mat[i,]   
   
    }
      Temp.sm <- array(all.vec,dim=dim(preproc))
      return(list(TempSmooth=Temp.sm,sm.mat=sm.mat))
  }  else  if (method=="circular"){
      sigma=bandwidth*0.37606506
      kermat <- GaussFilterMat(TimeLen,sigma,cut=4)      
   for (i in 1:mitolen) {
      sm.mat[i,]<-as.vector(kermat%*%pixel.temp[i,])
      k=mito.loc[i] 
## smoothed data put back in the full vec (length=256*256*150) 
      all.vec[length(mask)*(0:(TimeLen-1))+k] <- sm.mat[i,]   
   
    }
      Temp.sm <- array(all.vec,dim=dim(preproc))
      return(list(TempSmooth=Temp.sm,sm.mat=sm.mat))
    }
  }
  



