# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# File:   approx_banana.R
# Author: Jelmer Ypma
# Date:   5 July 2010
#
# Example showing how to solve the Rosenbrock Banana function
# with an approximated gradient, which doesn't work so well.

library('ipoptr')

# Rosenbrock Banana function
eval_f <- function(x) {   
    return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
}

# Approximate eval_f using finite differences
# http://en.wikipedia.org/wiki/Numerical_differentiation
approx_grad_f <- function( x ) {
    minAbsValue     <- 0
    stepSize        <- sqrt( .Machine$double.eps )

    # if we evaluate at 0, we need a different step size
    stepSizeVec     <- ifelse(abs(x) <= minAbsValue, stepSize^3, x * stepSize)
    
    x_prime <- x
    f       <- eval_f( x )
    grad_f  <- rep( 0, length(x) )
    for (i in 1:length(x)) {
        x_prime[i]      <- x[i] + stepSizeVec[i]
        stepSizeVec[i]  <- x_prime[i] - x[i]

        f_prime         <- eval_f( x_prime )
        grad_f[i]       <- ( f_prime - f )/stepSizeVec[i]
        x_prime[i]      <- x[i]
    }
    
    return( grad_f )
}

# initial values
x0 <- c( -1.2, 1 )

opts <- list("tol"=1.0e-8, "max_iter"=5000)
 
# solve Rosenbrock Banana function with approximated gradient
print( ipoptr( x0=x0, 
               eval_f=eval_f, 
               eval_grad_f=approx_grad_f, 
               opts=opts) )

opts <- list("tol"=1.0e-7)
 
# solve Rosenbrock Banana function with approximated gradient
# and lower tolerance
print( ipoptr( x0=x0, 
               eval_f=eval_f, 
               eval_grad_f=approx_grad_f, 
               opts=opts) )

               