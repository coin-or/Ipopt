\name{plot.sparseness}
\alias{plot.sparseness}
\title{
Plot sparseness structure of matrix
}
\description{
This function plots the sparseness structure of a matrix in the format that is required by ipoptr.
}
\usage{
plot.sparseness( x, pch='.', asp=1, xaxs='i', yaxs='i', ... )
}
\arguments{
  \item{x}{
    list of vectors with indices. Each element of the list corresponds to a row in 
    the matrix. Each index corresponds to a non-zero element in the matrix.
  }
  \item{pch}{
    plotting `character'. See \code{\link[graphics:plot]{points}} for possible values.
  }
  \item{asp}{
    aspect ratio, default = 1.
  }
  \item{xaxs, yaxs}{
    style of axis interval calculation, default = 'i' (do not extend the axis). See \code{\link[graphics:plot]{par}} for more information.
  }
  \item{...}{
    further graphical parameters that will be passed to \code{\link[graphics:plot]{plot}}.
  }
}
\value{
    A list with the non-zero x and y indices is returned.
}
\author{
Jelmer Ypma
}
\seealso{
\code{\link[ipoptr:ipoptr]{ipoptr}}
\code{\link[ipoptr:print.sparseness]{print.sparseness}}
\code{\link[ipoptr:make.sparse]{make.sparse}}
}
\examples{
library('ipoptr')

# use different plotting symbol for small matrices
plot.sparseness( make.sparse(diag(5)), pch='x' )

# plot large matrix example
s <- make.sparse( lower.tri( matrix( 1, 500, 500), diag=TRUE ) )
plot.sparseness( s )

# plot another large matrix
s <- do.call( "cbind", lapply( 1:5, function(i) { diag(5) \%x\% matrix(1, nrow=5, ncol=20) } ) )
s <- do.call( "rbind", lapply( 1:10,  function(i) { s } ) )
s <- cbind( matrix( 1, nrow=nrow(s), ncol=40 ), s )
plot.sparseness( make.sparse( s ) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
\keyword{ interface }
