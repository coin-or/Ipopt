\name{print.sparseness}
\alias{print.sparseness}
\title{
Show sparseness structure of matrix
}
\description{
This function shows the sparseness structure of a matrix in the format that is required by ipoptr.
}
\usage{
print.sparseness( x, indices=TRUE, data=NULL, ncol=NULL, ... )
}
\arguments{
  \item{x}{
    list of vectors with indices. Each element of the list corresponds to a row in 
    the matrix. Each index corresponds to a non-zero element in the matrix.
  }
  \item{indices}{
    Logical. Should we show the order of the non-zero elements or just whether an element is non-zero?
  }
  \item{data}{
    vector with non-zero elements of the sparse matrix.
  }
  \item{ncol}{
    integer supplying the number of columns of the sparse matrix. If this is not supplied, 
    we take the number of columns as the largest index in \code{s}.
  }
  \item{...}{
    further arguments passed to or from other methods.
  }
}
\value{
    A matrix showing the sparseness structure is returned.
}
\author{
Jelmer Ypma
}
\seealso{
\code{\link[ipoptr:ipoptr]{ipoptr}}
\code{\link[ipoptr:plot.sparseness]{plot.sparseness}}
\code{\link[ipoptr:make.sparse]{make.sparse}}
}
\examples{
library('ipoptr')

# print lower-diagonal 4x4 matrix
print.sparseness( list( c(1), c(1,2), c(1,2,3), c(1,2,3,4) ) )

# print diagonal 3x3 matrix without indices counts
print.sparseness( list( c(1), c(2), c(3) ), indices=FALSE )

# print a third sparse matrix
print.sparseness( list( c(1,3,6,8), c(2,5), c(3,7,9) ) )

# and a fourth one, where the elements are in a different order
print.sparseness( list( c(3,1,6,8), c(2,5), c(3,9,7) ) )

# print lower-diagonal 5x5 matrix generated with make.sparse
A_lower <- make.sparse( lower.tri( matrix(1, nrow=5, ncol=5), diag=TRUE ) )
print.sparseness( A_lower )

# print a diagonal 5x5 matrix without indices counts
A_diag  <- make.sparse( diag(5) > 0 )
print.sparseness( A_diag )

# example from tests/lasso.R
n <- 100    # number of observations
m <- 5      # number of variables

# define hessian function
hessian <- function( A ) {
    H <- t(A) %*% A
    H <- unlist( lapply( 1:m, function(i) { H[i,1:i] } ) )
    
    return( H )
}

# define the structure
hessian_structure <- c( lapply( 1:m, function(x) { return( c(1:x) ) } ),
                        lapply( 1:m, function(x) { return( c() ) } ) )

# generate data
set.seed( 3141 )
A <- hessian( matrix( rnorm( n*m ), nrow=n, ncol=m ) )                        
print.sparseness( x       = hessian_structure,
                  indices = TRUE,
                  data    = format( A, digits=2, nsmall=2, justify='right'),
                  ncol    = 2*m )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
\keyword{ interface }
