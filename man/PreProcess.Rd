\name{PreProcess}
\alias{PreProcess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Preparation Before Analysis 
}
\description{
   Perform background filtration, spacial smoothing on a 3-D array
}
\usage{
PreProcess(mydata, default.val = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mydata}{
       The 3D array that is to be processed.
}
  \item{default.val}{
       NA consider get rid of this
}
  \item{\dots}{
       Other ... consider get rid of this
}
}
\details{
   \code{mydata} is a 3D array with first 2 dimensions indicating the
   intensity for each pixel of a picture and the 3rd dimension indicating the time captured the picture. To separate the back-ground and front-ground, k-means method is applied to do the cluster analysis and then normalize the intensity by doing temporal quantile regression on each pixel and finally use \code{image.smooth} to smooth each picture.  
}
\value{
%%  ~Describe the value returned
  \item{preproc}{3-D array, the spacial smoothed \code{mydata} after pre-process}
  \item{mask}{A logic matrix indicating the background or not. \code{TRUE} denotes background.}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
