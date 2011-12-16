\name{TS.mat}
\alias{TS.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smooth Each Row of A Matrix}

\description{Consider each row of the input matrix is circular. Apply
  Nadaraya-Watson kernel regression estimate to do the smoothing.
}

\usage{
TS.mat(ymat,bandwidth)
}
\arguments{
  \item{ymat}{The matrix to be inputted.}
  \item{bandwidth}{the bandwidth to be used in Nadaraya-Waston kernel
    estimate. The same as \code{sigma*0.3706506}}
  }

%\details{}
%\value{An integer or a vector of the reverse rank of \code{x} in \code{x.perm}.}

\author{Tian Chen}
%\seealso{\code{\link{rank}}}
%\examples{}
%\keyword{utilities}% __ONLY ONE__ keyword per line