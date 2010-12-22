# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# File:   sparseness.R
# Author: Jelmer Ypma
# Date:   18 April 2010
#
# Example showing how the functions print.sparseness and 
# make.sparse work. These show and create the sparseness 
# structure of a matrix as it should be used for input 
# to ipoptr().

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

# make a large sparseness structure and use plot                 
s <- do.call( "cbind", lapply( 1:5, function(i) { diag(5) %x% matrix(1, nrow=5, ncol=20) } ) )
s <- do.call( "rbind", lapply( 1:5,  function(i) { s } ) )
s <- cbind( matrix( 1, nrow=nrow(s), ncol=40 ), s )
plot.sparseness( make.sparse( s ) )
