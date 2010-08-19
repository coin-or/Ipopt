# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Common Public License.
#
# File:   ipoptr.R
# Author: Jelmer Ypma
# Date:   18 April 2010
#
# Input: 
#		x0 : vector with initial values
#		eval_f : function to evaluate objective function
#		eval_grad_f : function to evaluate gradient of objective function
#		lb : lower bounds of the control
#		ub : upper bounds of the control
#		eval_g : function to evaluate (non-)linear constraints that should hold in the solution
#		eval_jac_g : function to evaluate the jacobian of the (non-)linear constraints that should hold in the solution
#		eval_jac_g_structure : sparseness structure of the jacobian
#		constraint_lb : lower bounds of the (non-)linear constraints
#		constraint_ub : upper bounds of the (non-)linear constraints
#		eval_h : function to evaluate the hessian
#		eval_h_structure : sparseness structure of the hessian
#		opts : list with options that are passed to Ipopt
#
# Output: structure with inputs and
#		call : the call that was made to solve
#		status : integer value with the status of the optimization (0 is success)
#		message : more informative message with the status of the optimization
#		iterations : number of iterations that were executed
#		objective : value if the objective function in the solution
#		solution : optimal value of the controls

ipoptr <-
function( x0, 
          eval_f, 
          eval_grad_f, 
          lb = NULL, 
          ub = NULL, 
          eval_g = function( x ) { return( numeric(0) ) }, 
          eval_jac_g = function( x ) { return( numeric(0) ) }, 
          eval_jac_g_structure = list(),
          constraint_lb = numeric(0), 
          constraint_ub = numeric(0),
          eval_h = NULL,
          eval_h_structure = NULL,
          opts = list(),
          ipoptr_environment = new.env() ) {
    
    # define 'infinite' lower and upper bounds of the control if they haven't been set
    if ( is.null( lb ) ) { lb <- rep( -Inf, length(x0) ) }
    if ( is.null( ub ) ) { ub <- rep(  Inf, length(x0) ) }
    
    # change the environment of the functions that we're calling
    environment( eval_f ) <- ipoptr_environment
    environment( eval_grad_f ) <- ipoptr_environment
    environment( eval_g ) <- ipoptr_environment
    environment( eval_jac_g ) <- ipoptr_environment
    
    # approximate Hessian
    if ( is.null( eval_h ) && is.null( eval_h_structure ) ) {
        opts$hessian_approximation <- "limited-memory"
    } else {
        environment( eval_h ) <- ipoptr_environment     # change environment
    }

    # build ipoptr object
    ret <- list( "x0"=x0, 
                 "eval_f"=eval_f, 
                 "eval_grad_f"=eval_grad_f, 
                 "lower_bounds"=lb, 
                 "upper_bounds"=ub, 
                 "eval_g"=eval_g, 
                 "eval_jac_g"=eval_jac_g, 
                 "constraint_lower_bounds"=constraint_lb, 
                 "constraint_upper_bounds"=constraint_ub, 
                 "eval_jac_g_structure"=eval_jac_g_structure,
                 "eval_h"=eval_h,
                 "eval_h_structure"=eval_h_structure,
                 "options"=get.option.types(opts),
                 "ipoptr_environment"=ipoptr_environment )
    
    attr(ret, "class") <- "ipoptr"
    
    # add the current call to the list
    ret$call <- match.call()
    
    # check whether we have a correctly formed ipoptr object
    is.ipoptr( ret )
    
    # pass ipoptr object to C code
    solution <- .Call( IpoptRSolve, ret )
    
    # remove the environment from the return object
    ret$environment <- NULL
    
    # add solution variables to object
    ret$status <- solution$status
    ret$message <- solution$message
    ret$iterations <- solution$iterations
    ret$objective <- solution$objective
    ret$solution <- solution$solution
    
    return( ret )
}
