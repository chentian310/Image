%\ A very bad documentation .....
\name{WriteData}
\alias{WriteData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write Data to a Directory
}
\description{
The data (3-d array) \code{mydata} are written to directory \code{DIR} in the form of .txt files (each one consists of  a matrix).The 3rd dimetion indicates the order of .txt files.
}
\usage{
WriteData(mydata, DIR, skip.first = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mydata}{the 3-D array to be written out
}

  \item{DIR}{
          The target directory \code{DIR} to write the data in 
}
  \item{skip.first}{
       logical: if \code{TRUE} then 3-D array \code{mydata} would not write the first matrix into the directory
}
}
\details{ NA
%%  ~~ If necessary, more details than the description above ~~
}
\value{ NA
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ NA
%% ~put references to the literature/web site here ~
}
\author{ NA
%%  ~~who you are~~
}
\note{ NA
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
      See Also as \code{\link{ReadData}}
}
\examples{ NA
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
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line