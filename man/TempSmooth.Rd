\name{TempSmooth}
\alias{TempSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Doing Temporal Smoothing for 3D array
}
\description{
  Do the  temporal smoothing for each pixel in the front-ground(mito
  region) by applying the Nadaray-Watson kernel regression estimate. 
}
\usage{
TempSmooth(preproc,mask,bandwidth,method=c("interval","circular"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{preproc}{
       The 3D array that is preprocessed: filtered background and spacial
       smoothed.
}
  \item{mask}{
       a logical matrix with \code{TRUE} indicating the background
       (non-mito region).
}
  \item{bandwidth}{
       the bandwidth. The kernels are scaled so that their quartiles are
       at +/- \code{0.25*bandwidth}.  
}
  \item{method}{
       the method to be used. \code{interval} simply applies
       \code{ksmooth}function in R and \code{circular} applies Gauss
       filter in a circular manner. 
}
 
\item{\dots}{
       Other ... consider get rid of this
}
}
%\details{
%     
%}
\value{
  %describe the value to be returned
  \item{TempSmooth}{3-D array, temporal smoothed \code{preproc}}
  \item{sm.mat}{A matrix. Each row is the intensity vector of a mito-region
    pixel after doing temporal smoothing }
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
