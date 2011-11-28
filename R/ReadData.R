ReadData <-
function(DIR,skip.first=FALSE) {



 Time<-as.integer(substr(dir(DIR),1,4)) 

 fname1<-paste(DIR,"/",sprintf("%04d",min(Time)),".txt",sep="")

 image1<-read.table(fname1)

#temp=scan(fname1)

#image1<-matrix(temp,nrow=sqrt(length(temp)),byrow=T)  # wrong if not square matrix 

 nrows<-dim(image1)[1];ncols=dim(image1)[2]

 mydata<-array(0,c(nrows,ncols,length(Time)))

 

 for (tt in 1:length(Time)){

  fname<-paste(DIR,"/",sprintf("%04d",min(Time)+tt-1),".txt",sep="")

  mydata[,,tt]<-matrix(scan(fname),nrow=nrows,byrow=T)

 }

 dimnames(mydata)<-list(1:nrows,1:ncols,paste("T",Time,sep="."))

#dimnames(mydata)<-list(1:nrows,1:ncols,paste("T",sprintf("%04d",Time),sep="")) # another "Time" format

 if(skip.first) {mydata<-mydata[,,-1]

  }else{mydata=mydata}

 return(mydata)



}

