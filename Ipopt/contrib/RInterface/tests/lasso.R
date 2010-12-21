# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# File:   lasso.R
# Author: Jelmer Ypma
# Date:   18 April 2010
#
# Example showing how to estimate a LASSO model. Based on
# the example from the Matlab interface to Ipopt. There are
# other packages in R that can estimate LASSO models, e.g.
# using the package glmnet.

library('ipoptr')

# Experiment parameters.
lambda <- 1                                # Level of L1 regularization.
n      <- 100                              # Number of training examples.
e      <- 1                                # Std. dev. in noise of outputs.
beta   <- c( 0, 0, 2, -4, 0, 0, -1, 3 )    # "True" regression coefficients.

# Set the random number generator seed.
ranseed <- 7
set.seed( ranseed )

# CREATE DATA SET.
# Generate the input vectors from the standard normal, and generate the
# responses from the regression with some additional noise. The variable 
# "beta" is the set of true regression coefficients.
m     <- length(beta)                           # Number of features.
A     <- matrix( rnorm(n*m), nrow=n, ncol=m )   # The n x m matrix of examples.
noise <- rnorm(n, sd=e)                         # Noise in outputs.
y     <- A %*% beta + noise                     # The outputs.


# DEFINE LASSO FUNCTIONS
# m, lambda, y, A are all defined in the ipoptr_environment
eval_f <- function(x) {
    # separate x in two parts
    w <- x[  1:m ]          # parameters
    u <- x[ (m+1):(2*m) ]
    
    return( sum( (y - A %*% w)^2 )/2 + lambda*sum(u) )
}

# ------------------------------------------------------------------
eval_grad_f <- function(x) {
    w <- x[ 1:m ]
    return( c( -t(A) %*% (y - A %*% w),  
               rep(lambda,m) ) )
}

# ------------------------------------------------------------------
eval_g <- function(x) {
    # separate x in two parts
    w <- x[  1:m ]          # parameters
    u <- x[ (m+1):(2*m) ]
   
    return( c( w + u, u - w ) )
}
  
# ------------------------------------------------------------------
#  J = [  I  I
#        -I  I ],
# where I is and identity matrix of size m
eval_jac_g <- function(x) {
    # return a vector of 1 and minus 1, since those are the values of the non-zero elements
    return( c( rep( 1, 2*m ), rep( c(-1,1), m ) ) )
}

# For m=5, The structure looks like this:
#  1  .  .  .  .  2  .  .  .  .
#  .  3  .  .  .  .  4  .  .  .
#  .  .  5  .  .  .  .  6  .  .
#  .  .  .  7  .  .  .  .  8  .
#  .  .  .  .  9  .  .  .  . 10
# 11  .  .  .  . 12  .  .  .  .
#  . 13  .  .  .  . 14  .  .  .
#  .  . 15  .  .  .  . 16  .  .
#  .  .  . 17  .  .  .  . 18  .
#  .  .  .  . 19  .  .  .  . 20
eval_jac_g_structure <- lapply( c(1:m,1:m), function(x) { return( c(x,m+x) ) } )

  
# ------------------------------------------------------------------
# rename lambda so it doesn't cause confusion with lambda in auxdata
eval_h <- function( x, obj_factor, hessian_lambda ) {
    H <- t(A) %*% A
    H <- unlist( lapply( 1:m, function(i) { H[i,1:i] } ) )
    
    return( obj_factor * H )
}

# For m=5, The structure looks like this:
#  1  .  .  .  .  .  .  .  .  .
#  2  3  .  .  .  .  .  .  .  .
#  4  5  6  .  .  .  .  .  .  .
#  7  8  9 10  .  .  .  .  .  .
# 11 12 13 14 15  .  .  .  .  .
#  .  .  .  .  .  .  .  .  .  .
#  .  .  .  .  .  .  .  .  .  .
#  .  .  .  .  .  .  .  .  .  .
#  .  .  .  .  .  .  .  .  .  .
#  .  .  .  .  .  .  .  .  .  .
eval_h_structure <- c( lapply( 1:m, function(x) { return( c(1:x) ) } ),
                       lapply( 1:m, function(x) { return( c() ) } ) )


# ------------------------------------------------------------------


# The starting point.
x0 = c( rep(0, m), 
        rep(1, m) )
       

# The constraint functions are bounded from below by zero.
constraint_lb = rep(   0, 2*m )
constraint_ub = rep( Inf, 2*m )

ipoptr_opts <- list( "jac_d_constant"   = 'yes',
                     "hessian_constant" = 'yes',
                     "mu_strategy"      = 'adaptive',
                     "max_iter"         = 100,
                     "tol"              = 1e-8 )

# Set up the auxiliary data.
auxdata <- new.env()
auxdata$m <- m
auxdata$A <- A
auxdata$y <- y
auxdata$lambda <- lambda
  
# COMPUTE SOLUTION WITH IPOPT.
# Compute the L1-regularized maximum likelihood estimator.
print( ipoptr( x0=x0, 
               eval_f=eval_f, 
               eval_grad_f=eval_grad_f, 
               eval_g=eval_g, 
               eval_jac_g=eval_jac_g,
               eval_jac_g_structure=eval_jac_g_structure,
               constraint_lb=constraint_lb, 
               constraint_ub=constraint_ub,
               eval_h=eval_h,
               eval_h_structure=eval_h_structure,
               opts=ipoptr_opts,
               ipoptr_environment=auxdata ) )

