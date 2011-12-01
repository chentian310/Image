\name{ReadData}
\alias{ReadData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Text Files Input
}
\description{
Reads in all text files in a directory and creates a 3-D array from it, with the 3rd dimension corresponding the index.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ReadData(DIR, skip.first = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DIR}{
     the directory containing the target text files.
%%     ~~Describe \code{DIR} here~~
}
  \item{skip.first}{
   logical: if \code{TRUE} then the first text file in the directory would not skipped and not inputed. 
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
