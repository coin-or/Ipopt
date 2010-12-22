# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# File:   make.sparse.R
# Author: Jelmer Ypma
# Date:   18 April 2010
#
# Input:  matrix with logical elements
# Output: list with as elements a vector of indices denoting non-zero (TRUE) elements of the matrix

make.sparse <- function( A ) {
    
    # start with empty list to append to
    S <- list()    
    
    # loop over matrix by row 
    for ( i in 1:nrow(A) ) {
        indices <- c()
        for ( j in 1:ncol(A) ) {
            if (A[i,j]) {
                indices <- c( indices, j )
            }
        }
        S <- c( S, list(indices) )
    }
    
    return( S )
}
