\name{print.ipoptr}
\alias{print.ipoptr}
\title{
Print results after running ipoptr
}
\description{
This function prints the ipoptr object that holds the results from a minimization using \code{ipoptr}.
}
\usage{
\method{print}{ipoptr}( x, show.controls=TRUE, \dots )
}
\arguments{
  \item{x}{
    object containing result from minimization.
  }
  \item{show.controls}{
    Logical or vector with indices. Should we show the value of the control variables in the solution? If code{show.controls} is a vector with indices, it is used to select which control variables should be shown. This can be useful if the model contains a set of parameters of interest and a set of nuisance parameters that are not of immediate interest.
  }
  \item{...}{
    further arguments passed to or from other methods.
  }
}
\author{
Jelmer Ypma
}
\seealso{
\code{\link[ipoptr:ipoptr]{ipoptr}}
\code{\link[ipoptr:make.sparse]{make.sparse}}
\code{\link[ipoptr:print.sparseness]{print.sparseness}}
}
\keyword{ optimize }
\keyword{ interface }
