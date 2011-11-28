WriteData <-
function(mydata,DIR,skip.first=FALSE){

  Time<-1:dim(mydata)[3]

  if (skip.first){

   for(tt in Time){

      fname=paste(DIR,"/",sprintf("%04d",tt),".txt",sep="")    

      write.table(mydata[,,tt],fname,row.names=FALSE,col.names=FALSE)

    } 

   }else{ for(tt in Time){

      fname=paste(DIR,"/",sprintf("%04d",tt-1),".txt",sep="")    

      write.table(mydata[,,tt],fname,row.names=FALSE,col.names=FALSE)

   }

 }

}

