# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# File:   print.ipoptr.R
# Author: Jelmer Ypma
# Date:   18 April 2010
#
# This function prints some basic output of a ipoptr 
# ojbect. The information is only available after it 
# has been solved.

print.ipoptr <- function(x, show.controls=TRUE, ...) {
	cat("\nCall:\n", deparse(x$call), "\n\n", sep = "", fill=TRUE)
	cat( unlist(strsplit(paste( "Ipopt solver status:", x$status, "(", x$message, ")\n" ),' ')), fill=TRUE )
    cat( paste( "Number of Iterations....:", x$iterations, "\n" ) )
	
    # if show.controls is TRUE or FALSE, show all or none of the controls
    if ( is.logical( show.controls ) ) {
        # show all control variables
        if ( show.controls ) {
            controls.indices = 1:length(x$solution)
        }
    }
    
    # if show.controls is a vector with indices, rename this vector
    # and define show.controls as TRUE
    if ( is.numeric( show.controls ) ) {
        controls.indices = show.controls
        show.controls = TRUE
    }
    
	# if solved successfully
	if ( x$status<=0 ) {
		cat( paste( "Optimal value of objective function: ", x$objective, "\n" ) )
		if ( show.controls ) {
            if ( length( controls.indices ) < length(x$solution) ) {
                cat( "Optimal value of user-defined subset of controls: " )
            } else {
                cat( "Optimal value of controls: " )
            }
            cat( x$solution[ controls.indices ], fill=TRUE)
            cat("\n")
        }
	} else {
		cat( paste( "Current value of objective function: ", x$objective, "\n" ) )
		if ( show.controls ) {
            if ( length( controls.indices ) < length(x$solution) ) {
                cat( "Current value of user-defined subset of controls: " )
            } else {
                cat( "Current value of controls: " )
            }
            cat( x$solution[ controls.indices ], fill=TRUE )
            cat("\n")
        }
    }
	cat("\n")
}
