% IPOPT Call the IPOPT constrained, nonlinear solver. 
%   The basic function call is
%   
%     [x, info] = IPOPT(x0,funcs,options)
%
%   The first input is either a matrix or a cell array of matrices. It
%   declares the starting point for the solver.
%
%   CALLBACK FUNCTIONS
%
%   The second input must be struct containing function handles for various
%   MATLAB routines. For more information on using functions and function
%   handles in MATLAB, type HELP FUNCTION and HELP FUNCTION_HANDLE in the
%   MATLAB prompt.
%
%     funcs.objective (required)
%
%     Calculates the objective function at the current point. It takes one
%     point, the current iterate x. For example, the definition of the
%     objective function for the Hock & Schittkowski (H&S) test problem #71
%     (with 4 optimization variables) would be
% 
%         function f = objective (x)
%           f = x(1)*x(4)*sum(x(1:3)) + x(3);
%         
%     funcs.gradient (required)
%
%     Computes the gradient of the objective at the current point. It takes
%     one input, the current iterate x. For H&S test problem #71, the
%     definition of the gradient callback would be
%
%         function g = gradient (x)
%           g = [ x(1)*x(4) + x(4)*sum(x(1:3))
%                 x(1)*x(4)
%                 x(1)*x(4) + 1
%                 x(1)*sum(x(1:3)) ]; 
%
%     funcs.constraints (optional)
%
%     This function is only required if there are constraints on your
%     variables. It evaluates the constraint functions at the current
%     point. It takes one input, x. The return value is a vector of length
%     equal to the number of constraints (it must be of the same length as
%     options.cl and options.cu). For H&S test problem #71, the
%     callback definition would be
%
%         function c = constraints (x)
%           c = [ prod(x); sum(x.^2) ];
%
%     funcs.jacobian (optional)
% 
%     This function is only required if there are constraints on your
%     variables. Evaluates the Jacobian of the constraints at the current
%     point. It takes one input, x. The output must always be an M x N
%     sparse matrix, where M is the number of constraints and N is the
%     number of variables. Type HELP SPARSE for more information on
%     constructing sparse matrices in MATLAB. The definition of the
%     callback function for H&S test problem #71 would be
%
%         function J = jacobian (x)
%           sparse([ prod(x)./x; 2*x ]);
%
%     Notice that the return value is a sparse matrix.
%
%     funcs.jacobianstructure (optional)
%
%     This function is only required if there are constraints on your
%     variables. It takes no inputs. The return value is a sparse
%     matrix whereby an entry is nonzero if and only if the Jacobian of
%     the constraints is nonzero at ANY point. The callback function for
%     the H&S test problem #71 simply returns a 2 x 4 matrix of ones in
%     the sparse matrix format:
%
%         function J = jacobianstructure() 
%           J = sparse(ones(2,4));
%
%     funcs.hessian (optional)
%
%     Evaluates the Hessian of the Lagrangian at the current point. It
%     must be specified unless you choose to use the limited-memory
%     quasi-Newton approximation to the Hessian (see below).
% 
%     The callback function has three inputs: the current point (x), a
%     scalar factor on the objective (sigma), and the Lagrange multipliers
%     (lambda), a vector of length equal to the number of constraints. The
%     function should compute
%                  
%        sigma*H + lambda(1)*G1 + ... + lambda(M)*GM
%
%     where M is the number of constraints, H is the Hessian of the
%     objective and the G's are the Hessians of the constraint
%     functions. The output must always be an N x N sparse, lower triangular
%     matrix, where N is the number of variables. In other words, if X is
%     the output value, then X must be the same as TRIL(X).
%
%     Here is an implementation of the Hessian callback routine for the
%     H&S test problem #71:
%
%         function H = hessian (x, sigma, lambda)
%           H = sigma*[ 2*x(4)             0      0   0;
%                       x(4)               0      0   0;
%                       x(4)               0      0   0;
%                       2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ];
%           H = H + lambda(1)*[    0          0         0         0;
%                               x(3)*x(4)     0         0         0;
%                               x(2)*x(4) x(1)*x(4)     0         0;
%                               x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ];
%           H = sparse(H + lambda(2)*diag([2 2 2 2]));
%  
%     funcs.hessianstructure (optional)
% 
%     This function serves the same purpose as funcs.jacobianstructure, but
%     for the Hessian matrix. Again, it is not needed if you are using the
%     limited-memory quasi-Newton approximation to the Hessian. It takes no
%     inputs, and must return a sparse, lower triangular matrix. For H&S
%     test problem #71, the MATLAB callback routine is fairly
%     straightforward:
%
%         function H = hessianstructure() 
%           H = sparse(tril(ones(4)));
%
%     funcs.iterfunc (optional)
%
%     An additional callback routine that is called once per algorithm
%     iteration. It takes two inputs, is the current iteration of the
%     algorithm, and the current value of the objective. This function
%     should always return true unless you want IPOPT to terminate
%     prematurely for whatever reason.
%
%   OPTIONS
%
%   The options are passed through the third input. What follows is a
%   description of the fields you may optionally specify.
%
%     options.lb  
%
%     Specify lower bounds on the variables. It must have the same number
%     of elements as x0. Set an entry to -Inf to specify no lower bound.
%
%     options.ub
%
%     Specify upper bounds on the variables. It must have the same number
%     of elements as x0. Set an entry to Inf to specify no upper bound.
%
%     options.cl, options.cu
%
%     Set lower and upper bounds on the constraints. Each should be a
%     vector of length equal to the number of constraints. As before, a
%     bound is removed by setting the entry to -Inf or +Inf. An equality
%     constraint is achieved by setting cl(i) = cu(i).
%
%     options.auxdata
%
%     Optionally, one may choose to pass additional auxiliary data to the
%     MATLAB callback routines listed above through the function call. For
%     instance, the objective callback function now takes two inputs, x and
%     auxdata. The auxiliary data may not change through the course of the
%     IPOPT optimization. The auxiliary data keep the same values as they
%     possessed in the initial call. If you need variables that change over
%     time, you may want to consider global variables (type HELP
%     GLOBAL). See the lasso.m file in the examples subdirectory for an
%     illustration of how the auxiliary data is passed to the various
%     callback functions.
%
%     options.zl, options.zu, options.lambda
%
%     These fields specify the initial value for the Lagrange multipliers,
%     which is especially useful for "warm starting" the interior-point
%     solver. They specify the Lagrange multipliers corresponding to the
%     lower bounds on the variables, upper bounds on the variables, and
%     constraints, respectively.
%
%     options.ipopt
%
%     Finally, you may also change the settings of IPOPT through this
%     field. For instance, to turn off the IPOPT output, use the
%     limited-memory BFGS approximation to the Hessian, and turn on the
%     derivative checker, do the following:
%
%       options.ipopt.print_level           = 0;
%       options.ipopt.hessian_approximation = 'limited-memory';
%       options.ipopt.derivative_test       = 'first-order';
%
%     For more details, see the documentation on the IPOPT website.
%
%   OUTPUTS
%
%   If the solver successfully converges to a stationary point or terminated
%   without an unrecoverable error, the function IPOPT outputs the candidate
%   solution x. In all other cases, an error is thrown. It also outputs some
%   additional information:
%
%     info.zl, info.zu, info.lambda
%
%     The value of the Lagrange multipliers at the solution. See the
%     "options" for more information on the Lagrange multipliers.
%
%     info.status
%
%     Upon termination, this field will take on one of these following
%     values (for a more up-to-date listing, see the IpReturnCodes.h header
%     file in the IPOPT C++ source directory):
%
%         0  solved
%         1  solved to acceptable level
%         2  infeasible problem detected
%         3  search direction becomes too small
%         4  diverging iterates
%         5  user requested stop
%     
%        -1  maximum number of iterations exceeded
%        -2  restoration phase failed
%        -3  error in step computation
%       -10  not enough degrees of freedom
%       -11  invalid problem definition
%       -12  invalid option
%       -13  invalid number detected
%
%      -100  unrecoverable exception
%      -101  non-IPOPT exception thrown
%      -102  insufficient memory
%      -199  internal error
%
%     info.iter, info.cpu
%
%     Number of iterations and CPU time (in seconds) taken by the Ipopt run
%
%   Finally, for more information, please consult the following webpages:
%
%      http://www.cs.ubc.ca/~pcarbo/ipopt-for-matlab
%      http://projects.coin-or.org/Ipopt
%
%   Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
%   This code is published under the Eclipse Public License.
%
%   Author: Peter Carbonetto
%           Dept. of Computer Science
%           University of British Columbia
%           September 19, 2008
