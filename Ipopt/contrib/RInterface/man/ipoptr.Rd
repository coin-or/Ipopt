\name{ipoptr}
\alias{ipoptr}
\title{
R interface to Ipopt
}
\description{
ipoptr is an R interface to Ipopt (Interior Point Optimizer), an open source 
software package for large-scale nonlinear optimization. It can be used to 
solve general nonlinear programming problems with nonlinear constraints and 
lower and upper bounds for the controls. Ipopt is written in C++ and is released 
as open source code under the Eclipse Public License (EPL). It is available from 
the COIN-OR initiative. The code has been written by Carl Laird and Andreas Waechter, 
who is the COIN project leader for Ipopt.  

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
\usage{
ipoptr( x0, 
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
        ipoptr_environment = new.env(),
        ... )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{
    vector with starting values for the optimization.
}
  \item{eval_f}{
    function that returns the value of the objective function.
}
  \item{eval_grad_f}{
    function that returns the value of the gradient of the objective function.
}
  \item{lb}{
    vector with lower bounds of the controls (use -1.0e19 for controls without lower bound).
}
  \item{ub}{
    vector with upper bounds of the controls (use  1.0e19 for controls without upper bound).
}
  \item{eval_g}{
    function to evaluate (non-)linear constraints that should hold in the solution.
}
  \item{eval_jac_g}{
    function to evaluate the jacobian of the (non-)linear constraints that should hold in the solution.
}
  \item{eval_jac_g_structure}{
    list of vectors with indices defining the sparseness structure of the Jacobian. 
    Each element of the list corresponds to a row in the matrix. Each index corresponds 
    to a non-zero element in the matrix (see also \code{\link[ipoptr:print.sparseness]{print.sparseness}}).
}
  \item{constraint_lb}{
    vector with lower bounds of the (non-)linear constraints
}
  \item{constraint_ub}{
    vector with upper bounds of the (non-)linear constraints
}
  \item{eval_h}{
    function to evaluate the hessian.
}
  \item{eval_h_structure}{
    list of vectors with indices defining the sparseness structure of the Hessian. 
    Each element of the list corresponds to a row in the matrix. Each index corresponds 
    to a non-zero element in the matrix (see also \code{\link[ipoptr:print.sparseness]{print.sparseness}}).
}
  \item{opts}{
    list with options, see examples below. For a full list of options use the option 
    "print_options_documentation"='yes', or have a look at the Ipopt documentation at
    \url{http://www.coin-or.org/Ipopt/documentation/}.
  }
  \item{ipoptr_environment}{
    environment that is used to evaluate the functions. Use this to pass 
    additional data or parameters to a function. See the second example in 
    \code{parameters.R} in the \code{tests} directory.
  }
  \item{...}{
    arguments that will be passed to the user-defined objective and constraints functions.
  }
}
\value{
    The return value contains a list with the inputs, and additional elements
    \item{call}{the call that was made to solve}
    \item{status}{integer value with the status of the optimization (0 is success)}
    \item{message}{more informative message with the status of the optimization}
    \item{iterations}{number of iterations that were executed}
    \item{objective}{value if the objective function in the solution}
    \item{solution}{optimal value of the controls}
}
\references{
A. Waechter and L. T. Biegler, On the Implementation of a Primal-Dual Interior Point Filter Line Search Algorithm for Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp. 25-57, 2006
}
\author{
Jelmer Ypma
}
\seealso{
\code{\link{optim}}
\code{\link{nlm}}
\code{\link{nlminb}}
\code{\link[Rsolnp:Rsolnp-package]{Rsolnp}}
\code{\link[Rsolnp:solnp]{ssolnp}}
\code{\link[ipoptr:print.sparseness]{print.sparseness}}
\code{\link[ipoptr:make.sparse]{make.sparse}}
}
\note{See ?`ipoptr-package` for an extended example.}
\examples{
library('ipoptr')

## Rosenbrock Banana function
eval_f <- function(x) {   
    return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
     }

## Gradient of Rosenbrock Banana function
eval_grad_f <- function(x) { 
         c(-400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
            200 * (x[2] - x[1] * x[1]))
     }

# The Hessian for this problem is actually dense, 
# This is a symmetric matrix, fill the lower left triangle only.
eval_h_structure <- list( c(1), c(1,2) )

eval_h <- function( x, obj_factor, hessian_lambda ) {
    
    return( obj_factor*c( 2 - 400*(x[2] - x[1]^2) + 800*x[1]^2,      # 1,1
                          -400*x[1],                                 # 2,1
                          200 ) )                                    # 2,2
}

# initial values
x0 <- c( -1.2, 1 )

opts <- list("print_level"=0,
             "file_print_level"=12,
             "output_file"="banana.out",
             "tol"=1.0e-8)
 
# solve Rosenbrock Banana function with analytic hessian 
print( ipoptr( x0=x0, 
               eval_f=eval_f, 
               eval_grad_f=eval_grad_f, 
               eval_h=eval_h,
               eval_h_structure=eval_h_structure,
               opts=opts) )

# solve Rosenbrock Banana function with approximated hessian			   
print( ipoptr( x0=x0, 
               eval_f=eval_f, 
               eval_grad_f=eval_grad_f, 
               opts=opts) )



##
#
# Solve the example taken from the Ipopt C++
# tutorial document (see Examples/CppTutorial/).
#
# min_x f(x) = -(x2-2)^2
#  s.t.
#       0 = x1^2 + x2 - 1
#       -1 <= x1 <= 1
#
##

eval_f <- function( x ) { 
    print( paste( "In R::eval_f, x = ", paste( c(1,2), collapse=', ' ) ) )

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
               
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
\keyword{ interface }

