\name{is.ipoptr}
\alias{is.ipoptr}
\title{
R interface to Ipopt
}
\description{
is.ipoptr preforms checks to see if a fully specified problem is supplied to ipoptr. Mostly for internal use.
}
\usage{
is.ipoptr( x )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object to be tested.
  }
}
\value{
    Logical. Return TRUE if all tests were passed, otherwise return FALSE or exit with Error.
}
\references{
A. Waechter and L. T. Biegler, On the Implementation of a Primal-Dual Interior Point Filter Line Search Algorithm for Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp. 25-57, 2006
}
\author{
Jelmer Ypma
}
\seealso{
\code{\link[ipoptr:ipoptr]{ipoptr}}
\code{\link[ipoptr:print.sparseness]{print.sparseness}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
\keyword{ interface }
