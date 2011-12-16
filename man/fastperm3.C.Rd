\name{fastperm3.C}
\alias{fastperm3.C}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fast Permutation Implemented in C for L=3 }
\description{
   For inference purpose. An implementation of random permutation for L=3 and smoothing using
   Gaussian kernel and only return the minimum and maximum value from
   each permutation. The permutation is implemented in C for time efficiency.
}
\usage{fastperm3.C(X,bandwidth)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The data to be permuted and smoothed}
  \item{bandwidth}{The bandwidth.The kernels are scaled so that their
    quartiles are at +/- \code{0.25*bandwidth}}
}
%\details{
%  
%}
\value{
%%  ~Describe the value returned
  \item{t.min}{a vector. The minimum value obtained from each permutation.}
  \item{t.max}{a vector. The maximum value obtained from each permutation.}
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
