\name{ipoptr-package}
\alias{ipoptr-package}
\docType{package}
\title{
R interface to Ipopt
}
\description{
ipoptr is an R interface to Ipopt (Interior Point Optimizer), an open source software 
package for large-scale nonlinear optimization. It can be used to solve general nonlinear 
programming problems with nonlinear constraints and lower and upper bounds for the controls. 
Ipopt is written in C++ and is released as open source code under the Eclipse Public License (EPL). 
It is available from the COIN-OR initiative. The code has been written by Carl Laird and 
Andreas Waechter, who is the COIN project leader for Ipopt.  

Ipopt is designed to find (local) solutions of mathematical optimization problems of the from

   min     f(x)
x in R^n

s.t.       g_L <= g(x) <= g_U
           x_L <=  x   <= x_U

where f(x): R^n --> R is the objective function, and g(x): R^n --> R^m are the constraint 
functions. The vectors g_L and g_U denote the lower and upper bounds on the constraints, 
and the vectors x_L and x_U are the bounds on the variables x. The functions f(x) and g(x) 
can be nonlinear and nonconvex, but should be twice continuously differentiable. Note that 
equality constraints can be formulated in the above formulation by setting the corresponding 
components of g_L and g_U to the same value.
}
\author{
Jelmer Ypma
}
\references{
A. Waechter and L. T. Biegler, On the Implementation of a Primal-Dual Interior Point Filter Line Search Algorithm for Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp. 25-57, 2006
}
\keyword{ optimize }
\keyword{ interface }
\seealso{
\code{\link{optim}}
\code{\link{nlm}}
\code{\link{nlminb}}
\code{\link[Rsolnp:Rsolnp-package]{Rsolnp}}
}
\note{See ?ipoptr for more examples.}
\examples{
# Example problem, number 71 from the Hock-Schittkowsky test suite
#
# \min_{x} x1*x4*(x1 + x2 + x3) + x3
# s.t.
#    x1*x2*x3*x4 >= 25
#    x1^2 + x2^2 + x3^2 + x4^2 = 40
#    1 <= x1,x2,x3,x4 <= 5
#
# x0 = (1,5,5,1)
#
# optimal solution = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
#
# Adapted from the Ipopt C++ interface example.


library('ipoptr')

#
# f(x) = x1*x4*(x1 + x2 + x3) + x3
#
eval_f <- function( x ) { 
    return( x[1]*x[4]*(x[1] + x[2] + x[3]) + x[3] ) 
}

eval_grad_f <- function( x ) {
    return( c( x[1] * x[4] + x[4] * (x[1] + x[2] + x[3]),
               x[1] * x[4],
               x[1] * x[4] + 1.0,
               x[1] * (x[1] + x[2] + x[3]) ) )		  
}

# constraint functions
eval_g <- function( x ) {
    return( c( x[1] * x[2] * x[3] * x[4],
               x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 ) )
}

# The Jacobian for this problem is dense
eval_jac_g_structure <- list( c(1,2,3,4), c(1,2,3,4) )

eval_jac_g <- function( x ) {
    return( c ( x[2]*x[3]*x[4],
                x[1]*x[3]*x[4],
                x[1]*x[2]*x[4],
                x[1]*x[2]*x[3],
                2.0*x[1],
                2.0*x[2],
                2.0*x[3],
                2.0*x[4] ) )
}

# The Hessian for this problem is actually dense, 
# This is a symmetric matrix, fill the lower left triangle only.
eval_h_structure <- list( c(1), c(1,2), c(1,2,3), c(1,2,3,4) )

eval_h <- function( x, obj_factor, hessian_lambda ) {

    values <- numeric(10)
    values[1] = obj_factor * (2*x[4]) # 1,1

    values[2] = obj_factor * (x[4])   # 2,1
    values[3] = 0                     # 2,2

    values[4] = obj_factor * (x[4])   # 3,1
    values[5] = 0                     # 4,2
    values[6] = 0                     # 3,3

    values[7] = obj_factor * (2*x[1] + x[2] + x[3]) # 4,1
    values[8] = obj_factor * (x[1])                 # 4,2
    values[9] = obj_factor * (x[1])                 # 4,3
    values[10] = 0                                  # 4,4


    # add the portion for the first constraint
    values[2] = values[2] + hessian_lambda[1] * (x[3] * x[4]) # 2,1
    
    values[4] = values[4] + hessian_lambda[1] * (x[2] * x[4]) # 3,1
    values[5] = values[5] + hessian_lambda[1] * (x[1] * x[4]) # 3,2

    values[7] = values[7] + hessian_lambda[1] * (x[2] * x[3]) # 4,1
    values[8] = values[8] + hessian_lambda[1] * (x[1] * x[3]) # 4,2
    values[9] = values[9] + hessian_lambda[1] * (x[1] * x[2]) # 4,3

    # add the portion for the second constraint
    values[1] = values[1] + hessian_lambda[2] * 2 # 1,1
    values[3] = values[3] + hessian_lambda[2] * 2 # 2,2
    values[6] = values[6] + hessian_lambda[2] * 2 # 3,3
    values[10] = values[10] + hessian_lambda[2] * 2 # 4,4

    return ( values )
}

# initial values
x0 <- c( 1, 5, 5, 1 )

# lower and upper bounds of control
lb <- c( 1, 1, 1, 1 )
ub <- c( 5, 5, 5, 5 )

# lower and upper bounds of constraints
constraint_lb <- c(  25, 40 )
constraint_ub <- c( Inf, 40 )


opts <- list("print_level"=0,
             "file_print_level"=12,
             "output_file"="hs071_nlp.out")
  
print( ipoptr( x0=x0, 
               eval_f=eval_f, 
               eval_grad_f=eval_grad_f, 
               lb=lb, 
               ub=ub, 
               eval_g=eval_g, 
               eval_jac_g=eval_jac_g, 
               constraint_lb=constraint_lb, 
               constraint_ub=constraint_ub, 
               eval_jac_g_structure=eval_jac_g_structure,
               eval_h=eval_h,
               eval_h_structure=eval_h_structure,
               opts=opts) )
}
