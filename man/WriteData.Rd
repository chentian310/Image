%\ A very bad documentation .....
\name{WriteData}
\alias{WriteData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write Data to a Directory
}
\description{
The data (3-d array) \texttt{mydata} are written to directory \texttt{DIR} in the form of .txt files (each one consists of  a matrix).The 3rd dimetion indicates the order of .txt files.
}
\usage{
WriteData(mydata, DIR, skip.first = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mydata} { the 3-D array to be written out\code{mydata}

  \item{DIR}{
          The target directory \code{DIR} 
}
  \item{skip.first}{
       if \texttt{TRUE} the 3-D array \texttt{mydata} \code{skip.first} here~~
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
       ReadData()   See Also as \code{\link{help}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(mydata,DIR,skip.first=FALSE){

  Time<-1:dim(mydata)[3]

  if (skip.first){

   for(tt in Time){

      fname=paste(DIR,"/",sprintf("\%04d",tt),".txt",sep="")    

      write.table(mydata[,,tt],fname,row.names=FALSE,col.names=FALSE)

    } 

   }else{ for(tt in Time){

      fname=paste(DIR,"/",sprintf("\%04d",tt-1),".txt",sep="")    

      write.table(mydata[,,tt],fname,row.names=FALSE,col.names=FALSE)

   }

 }

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
