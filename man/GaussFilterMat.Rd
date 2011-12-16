\name{GaussFilterMat}
\alias{GaussFilterMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gaussian Filter Matrix}
\description{
   An matrix implemention of Nadaraya-Watson kernel smoothing when
   considering the data is circular. 
}
\usage{GaussFilterMat(N,sigma,cut=4)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{The data to be permuted and smoothed}
  \item{sigma}{sigma of the Gaussian kernel. The same as \code{0.3706506*bandwidth}}
  \item{cut}{default is 4. The Gaussian kernel is zero beyond the value
    of \code{sigma*cut/2}}
}
%\details{
%  
%}
%\value{
%  
%}
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
