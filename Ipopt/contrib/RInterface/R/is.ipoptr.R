# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# File:   is.ipoptr.R
# Author: Jelmer Ypma
# Date:   18 April 2010
#
# Input: object
# Output: bool telling whether the object is an ipoptr or not

is.ipoptr <- function(x) {
    
    # Check whether the object exists and is a list
    if( is.null(x) ) { return( FALSE ) }
    if( !is.list(x) ) { return( FALSE ) }

    # Check whether a correct environment was specified
    stopifnot( is.environment(x$ipoptr_environment) )
   
    # Define local flag defining whether we approximate the Hessian or not
    flag_hessian_approximation = FALSE
    if ( !is.null( x$options$string$hessian_approximation ) ) {
        flag_hessian_approximation = ( x$options$string$hessian_approximation == "limited-memory" )
    }
    
    # Check whether the needed functions are supplied
    stopifnot( is.function(x$eval_f) )
    stopifnot( is.function(x$eval_grad_f) )
    stopifnot( is.function(x$eval_g) )
    stopifnot( is.function(x$eval_jac_g) )
    if ( !flag_hessian_approximation ) { stopifnot( is.function(x$eval_h) ) }
    
    # Check whether bounds are defined for all controls
    stopifnot( length( x$x0 ) == length( x$lower_bounds ) )
    stopifnot( length( x$x0 ) == length( x$upper_bounds ) )
    
    # Check whether the initial value is within the bounds
    stopifnot( all( x$x0 >= x$lower_bounds ) )
    stopifnot( all( x$x0 <= x$upper_bounds ) )
    
    num.controls <- length( x$x0 )
    num.constraints <- length( x$constraint_lower_bounds )
    
    # Check the length of some return values
    stopifnot( length(x$eval_f( x$x0 ))==1 )
    stopifnot( length(x$eval_grad_f( x$x0 ))==num.controls )
    stopifnot( length(x$eval_g( x$x0 ))==num.constraints )
    stopifnot( length(x$eval_jac_g( x$x0 ))==length(unlist(x$eval_jac_g_structure)) )		# the number of non-zero elements in the Jacobian
    if ( !flag_hessian_approximation ) { 
        stopifnot( length(x$eval_h( x$x0, 1, rep(1,num.constraints) ))==length(unlist(x$eval_h_structure)) )		# the number of non-zero elements in the Hessian
    }
    
    # Check the whether we don't have NA's in initial values
    stopifnot( all(!is.na(x$eval_f( x$x0 ))) )
    stopifnot( all(!is.na(x$eval_grad_f( x$x0 ))) )
    stopifnot( all(!is.na(x$eval_g( x$x0 ))) )
    stopifnot( all(!is.na(x$eval_jac_g( x$x0 ))) )		# the number of non-zero elements in the Jacobian
    if ( !flag_hessian_approximation ) { 
        stopifnot( all(!is.na(x$eval_h( x$x0, 1, rep(1,num.constraints) ))) )		# the number of non-zero elements in the Hessian
    }
    
    # Check whether a correct structure was supplied, and check the size
    stopifnot( is.list(x$eval_jac_g_structure) )
    
    stopifnot( length(x$eval_jac_g_structure)==num.constraints )
    if ( !flag_hessian_approximation ) {
        stopifnot( length(x$eval_h_structure)==num.controls )
        stopifnot( is.list(x$eval_h_structure) )
    }
    
    # Check the number of non-linear constraints
    stopifnot( length(x$constraint_lower_bounds)==length(x$constraint_upper_bounds) )
    
    # Check whether none of the non-zero indices are larger than the number of controls
    # Also, the smallest index should be bigger than 0
    if ( length( x$eval_jac_g_structure ) > 0 ) {
        stopifnot( max(unlist(x$eval_jac_g_structure)) <= num.controls )
        stopifnot( min(unlist(x$eval_jac_g_structure)) > 0 )
    }
    if ( !flag_hessian_approximation ) {
        stopifnot( max(unlist(x$eval_h_structure)) <= num.controls )
        stopifnot( min(unlist(x$eval_h_structure)) > 0 )
    }
    
    # Check whether option to approximate hessian and eval_h are both set
    # If we approximate the hessian, then we don't want to set eval_h
    if ( flag_hessian_approximation ) {
        if( !is.null( x$eval_h ) ) {
            warning("Option supplied to approximate hessian, but eval_h is defined.\nSolution: remove option hessian_approximation=limited-memory to use analytic derivatives.")
        }
        if( !is.null( x$eval_h_structure ) ) {
            warning("Option supplied to approximate hessian, but eval_h_structure is defined.\nSolution: remove option hessian_approximation=limited-memory to use analytic derivatives.")
        }
    }
    

    return( TRUE )
}
