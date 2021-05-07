# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# File:   plot.sparseness.R
# Author: Jelmer Ypma
# Date:   23 June 2010
#
# Input: sparse matrix structure (as list with non-zero indices)
# Output: plot a the non-zero elements in the matrix as dots
#         useful for matrices with many elements, for smaller ones
#         use print.sparseness

plot.sparseness <- function( x, pch='.', asp=1, xaxs='i', yaxs='i', ... ) {
    # make a list of y indices corresponding to the non-zero x indices
    structure.y <- lapply( 1:length( x ), function(i) { rep( i, length( x[[i]] ) ) } )

    indices.x <- unlist( x )
    indices.y <- unlist( structure.y )

    # plot non-zero elements, where we revert the y-axis (top-left element is 1,1),
    # fix the aspect ratio (asp=1) and do not extend the x and y axis (x/yaxs='i')
    plot( indices.x, 
          indices.y, 
          xlim=c(min(indices.x), max(indices.x)), 
          ylim=c(max(indices.y), min(indices.y)), 
          type='p', 
          pch=pch, 
          asp=asp,
          xaxs=xaxs,
          yaxs=yaxs,
          ... )
    
    return( list( x=indices.x, y=indices.y ) )
}
