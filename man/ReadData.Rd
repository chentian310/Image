\name{ReadData}
\alias{ReadData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ReadData(DIR, skip.first = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DIR}{
%%     ~~Describe \code{DIR} here~~
}
  \item{skip.first}{
%%     ~~Describe \code{skip.first} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(DIR,skip.first=FALSE) {



 Time<-as.integer(substr(dir(DIR),1,4)) 

 fname1<-paste(DIR,"/",sprintf("\%04d",min(Time)),".txt",sep="")

 image1<-read.table(fname1)

#temp=scan(fname1)

#image1<-matrix(temp,nrow=sqrt(length(temp)),byrow=T)  # wrong if not square matrix 

 nrows<-dim(image1)[1];ncols=dim(image1)[2]

 mydata<-array(0,c(nrows,ncols,length(Time)))

 

 for (tt in 1:length(Time)){

  fname<-paste(DIR,"/",sprintf("\%04d",min(Time)+tt-1),".txt",sep="")

  mydata[,,tt]<-matrix(scan(fname),nrow=nrows,byrow=T)

 }

 dimnames(mydata)<-list(1:nrows,1:ncols,paste("T",Time,sep="."))

#dimnames(mydata)<-list(1:nrows,1:ncols,paste("T",sprintf("\%04d",Time),sep="")) # another "Time" format

 if(skip.first) {mydata<-mydata[,,-1]

  }else{mydata=mydata}

 return(mydata)



  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
