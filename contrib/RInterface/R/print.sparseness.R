# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# File:   print.sparseness.R
# Author: Jelmer Ypma
# Date:   18 April 2010
#
# Input: sparse matrix structure (as list with non-zero indices)
# Output: print a table with 'x' for non-zero element and '.' for zero element

print.sparseness <- function( x, indices=TRUE, data=NULL, ncol=NULL, ... ) {
    stopifnot( is.list(x) )
    
    # if number of columns is not supplied, take it as the maximum 
    # value of the indices
    if ( is.null(ncol) ) {
        ncol <- max(unlist(x))
    }
    
    # create matrix with dots
    p <- data.frame( matrix( ".", nrow=length(x), ncol ), stringsAsFactors=FALSE )
    names( p ) <- 1:ncol
    
    # change dots by 'x' or count of index
    cnt=1
    for ( row in 1:length(x) ) {
        for ( col in x[[row]] ) {
            if ( indices ) {
                if ( is.null( data ) ) {
                    p[ row, col ] <- cnt
                } else {
                    p[ row, col ] <- paste( cnt, ':', data[cnt] )
                }
            } else {
                if ( is.null( data ) ) {
                    p[ row, col ] <- 'x'
                } else {
                    p[ row, col ] <- data[cnt]
                }
            }
            cnt = cnt+1
        }	
    }
    
    return( p )
}
