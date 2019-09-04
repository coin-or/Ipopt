# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# File:   mynlp.R
# Author: Jelmer Ypma
# Date:   18 April 2010
#
# Example NLP for interfacing a problem with IPOPT.
# This example is adapted from the C++ example that 
# goes along with the Ipopt tutorial document.
# This example solves the following problem:
#
# min_x f(x) = -(x2-2)^2
#  s.t.
#       0 = x1^2 + x2 - 1
#       -1 <= x1 <= 1

library('ipoptr')

eval_f <- function( x ) { 
    print( paste( "In R::eval_f, x = ", paste( x, collapse=', ' ) ) )

    return( -(x[2] - 2.0)*(x[2] - 2.0) ) 
}

eval_grad_f <- function( x ) {
    return( c(0.0, -2.0*(x[2] - 2.0) ) )
}

eval_g <- function( x ) {
    return( -(x[1]*x[1] + x[2] - 1.0) );
}

# list with indices of non-zero elements
# each element of the list corresponds to the derivative of one constraint
#
# e.g.
#      / 0 x x \
#      \ x 0 x /
# would be
# list( c(2,3), c(1,3) )
eval_jac_g_structure <- list( c(1,2) )


# this should return a vector with all the non-zero elements
# so, no list here, because that is slower I guess
# TODO: make an R-function that shows the structure in matrix form
eval_jac_g <- function( x ) {
    return ( c ( -2.0 * x[1], -1.0 ) )
}


# diagonal matrix, usually only fill the lower triangle
eval_h_structure <- list( c(1), c(2) )

eval_h <- function( x, obj_factor, hessian_lambda ) {
    return ( c( -2.0*hessian_lambda[1], -2.0*obj_factor ) )
}

x0 <- c(0.5,1.5)

lb <- c( -1, -1.0e19 )
ub <- c(  1,  1.0e19 )

constraint_lb <- 0
constraint_ub <- 0
  
opts <- list("print_level"=0,
             "file_print_level"=12,
             "output_file"="ipopttest.out")
  
print( ipoptr( x0=x0, 
               eval_f=eval_f, 
               eval_grad_f=eval_grad_f, 
               lb=lb, 
               ub=ub, 
               eval_g=eval_g, 
               eval_jac_g=eval_jac_g,
               eval_jac_g_structure=eval_jac_g_structure,
               constraint_lb=constraint_lb,
               constraint_ub=constraint_ub,
               eval_h=eval_h,
               eval_h_structure=eval_h_structure,
               opts=opts) )
