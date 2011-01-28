\name{make.sparse}
\alias{make.sparse}
\title{
Create sparseness structure from logical matrix
}
\description{
This function creates the sparseness structure of a logical matrix in the format that is required by ipoptr.
}
\usage{
make.sparse( A )
}
\arguments{
  \item{A}{
    Matrix with logicals. TRUE denotes a non-zero element in the matrix.
  }
}
\value{
    List of vectors with indices. Each element of the list corresponds to a row in 
    the matrix. Each index corresponds to a non-zero element in the matrix.
}
\author{
Jelmer Ypma
}
\seealso{
\code{\link[ipoptr:ipoptr]{ipoptr}}
\code{\link[ipoptr:print.sparseness]{print.sparseness}}
}
\examples{
library('ipoptr')
# print lower-diagonal 5x5 matrix generated with make.sparse
A_lower <- make.sparse( lower.tri( matrix(1, nrow=5, ncol=5), diag=TRUE ) )
print.sparseness( A_lower )

# prnit a diagonal 5x5 matrix without indices counts
A_diag  <- make.sparse( diag(5) > 0 )
print.sparseness( A_diag )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
\keyword{ interface }
