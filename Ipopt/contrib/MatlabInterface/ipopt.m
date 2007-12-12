% IPOPT Call the IPOPT constrained, nonlinear solver. 
%   The basic function call is
%   
%     x = IPOPT(x0,lb,ub,constraintlb,constraintub,objfunc,gradfunc,...
%               constraintfunc,jacobianfunc,hessianfunc)
%
%   The first input argument x0 is either a matrix or a cell array of
%   matrices. It declares the starting point for the solver.
%
%   The second and third input arguments lb and ub must be of the same
%   structure as x0. They declare the lower and upper bounds on the
%   variables, respectively. Set an entry of lb to Inf to declare no lower
%   bound, and likewise for the upper bounds.
%
%   The fourth and fifth input arguments constraintlb and constraintub set
%   the lower and upper bounds on the constraints, respectively. Each should
%   be a vector of length equal to the number of constraints. As before, a
%   bound is removed by setting the entry to Inf. An equality constraint is
%   achieved by setting constraintlb(i) = constraintub(i).
%
%   If the solver successfully converges to a stationary point or terminated
%   without an unrecoverable error, the function IPOPT outputs the candidate
%   solution. In all other cases, an error is thrown. The number of outputs
%   is equal to the number of cell entries in x0. If x0 is a matrix, IPOPT
%   produces a single output.
%
%   EXAMPLE. Here is an example of a basic call to IPOPT, in which the
%   optimization variables are specified by an array of length n and another
%   array of length m:
%
%   [r s] = ipopt({r s},{ zeros(n,1) zeros(m,1) },...
%                 { repmat(inf,n,1) repmat(inf,m,1) },lb,ub,...
%                 @computeobjective,@computegradient,@computeconstraints',...
%                 @computejacobian,@computeJGHessian,auxdata);
%
%   The optimization variables are bounded below by 0. The bounds on the
%   equality and inequality constraints are specified by the MATLAB
%   variables lb and ub. Upon convergence to a feasible, stationary point,
%   the routine will return the values of the primal variables r and
%   s. Also, auxiliary data ('auxdata') is passed to all the callback
%   routines. For more information, continue reading below.
%
%   The remaining input arguments must be function handles for various 
%   MATLAB routines (M-files or subfunctions):
%
%     objfunc         Calculates the objective function at the current
%                     point. The routine must accept as many inputs as cell
%                     entries in x0 (or one input if x0 is a matrix). The
%                     output must be a scalar representing the objective
%                     evaluated at the current point.
%         
%     gradfunc        Computes the gradient of the objective at the current
%                     point. The input is the same as objfunc, but it must
%                     return as many outputs as there are inputs, and each
%                     of the outputs must have the same matrix structure
%                     as its corresponding input.
%
%     constraintfunc  Evaluates the constraint functions at the current
%                     point. The MATLAB routine must accept the same inputs
%                     as objfunc. The return value is a vector of length
%                     equal to the number of constraints (it must be of the
%                     same length as constraintlb and constraintub). If
%                     there are no constraints, constraintfunc may be set
%                     to the empty string because the MATLAB callback
%                     routine will never be called.
%
%     jacobianfunc    Evaluates the Jacobian of the constraints at the
%                     current point. The inputs are the same as the
%                     previously described callback functions, except that
%                     it has an additional input "returnStructureOnly" which
%                     will either true or false. If returnStructureOnly is
%                     false, return the Jacobian matrix evaluated at the
%                     specified point. If returnStructureOnly is true,
%                     return a sparse matrix whereby an entry is zero if and
%                     only if the Jacobian of the constraints is zero at
%                     EVERY possible (feasible or non-feasible) point. The
%                     output must always be an M x N sparse matrix, where M
%                     is the number of constraints and N is the number of
%                     variables. Type HELP SPARSE for more information on
%                     constructing sparse matrices in MATLAB. If there
%                     are no constraints, jacobianfunc may be set to the
%                     empty string because the MATLAB callback routine
%                     will never be called.
%
%     hessianfunc     Evaluates the Hessian of the Lagrangian at the current
%                     point. Its inputs are the same as the objective and
%                     gradient callback functions, but it also has three
%                     additional inputs. The first additional input is
%                     "sigma". It is a scalar factor on the
%                     objective. The second additional input is "lambda",
%                     which is a vector of length equal to the number of
%                     constraints. The remaining additional input is
%                     returnStructureOnly (see the description of
%                     jacobianfunc above). The function should compute
%                  
%                       sigma*H + lambda(1)*G1 + ... + lambda(M)*GM
%
%                     where M is the number of constraints, H is the Hessian
%                     of the objective and the G's are the Hessians of the
%                     constraint functions. The output must always be an N x
%                     N sparse, lower triangular matrix, where N is the
%                     number of variables. In other words, if X is the
%                     output value, then X must be the same as TRIL(X).
%
%   For more information on using functions and function handles in MATLAB,
%   type HELP FUNCTION and HELP FUNCTION_HANDLE in the MATLAB prompt.
%
%   Optionally, one may choose to pass additional auxiliary data to the
%   MATLAB callback routines listed above through the function call
%   IPOPT(...,auxdata). If auxdata is the empty matrix, no extra information
%   is passed. It is important to observe that the auxiliary data MAY NOT
%   CHANGE through the course of the IPOPT optimization! The auxiliary data
%   keep the same values as they possessed in the initial call. If you need
%   variables that change over time, you may want to consider global
%   variables (type HELP GLOBAL).
%
%   EXAMPLE (continued). Returning to the above example, the 5 callback
%   routines would be of the following form:
%
%     function F        = computeobjective   (r,s,auxdata)
%     function [dr, ds] = computegradient    (r,s,auxdata)
%     function g        = computeconstraints (r,s,auxdata) 
%     function J        = computejacobian    (r,s,returnStructOnly,auxdata)
%     function H        = computehessian     (r,s,sigma,lambda,...
%                                             returnStructOnly,auxdata)
%
%   IPOPT(...,auxdata,iterfunc) specifies a a function handle to an
%   additional callback routine which is called once per algorithm
%   iteration. The callback routine must take the form ITERFUNC(T,F,AUXDATA). 
%   T is the current iteration of the algorithm. F is the current value of
%   the objective. Finally, extra information may be passed through the
%   input AUXDATA. No outputs are expected from iterfunc. Also note that the
%   value of the optimization variables is not passed to this callback
%   routine. If iterfunc is the empty string, no routine is called.
%
%   By default, the print level is set so that IPOPT displays the
%   progress of the algorithm. However, if the user implements an
%   iterative callback routine, the print level is automatically set so
%   that IPOPT does not generate output to the MATLAB console.
%
%   IPOPT(...,auxdata,iterfunc,multipliers) specifies the initial values for
%   the Lagrange multipliers. The input "multipliers" must have three
%   fields. The field "zl" must be an array of length equal to the total
%   number of optimization variables. It contains the initial values for the
%   Lagrange multipliers associated with the lower bounds on the
%   optimization variables. The field "zu" is an array of the same
%   length. It contains the initial values for the upper bounds Lagrange
%   multipliers. The field "lambda" must be an array of length equal to the
%   number of constraints. It specifies the initial values for the Lagrange
%   multipliers associated with the equality and inequality constraints.
%
%   These inputs may be followed by parameter/value pairs to modify the
%   default IPOPT algorithm options. See the IPOPT documentation for more
%   information. If the IPOPT option 'hessian_approximation' is set to
%   'limited-memory', the Hessian callback routine will never be called
%   so it may be set to the empty string.
%
%   [..., multipliers, numiter] = IPOPT(...) returns some additional
%   information to the user: the values of the Lagrange multipliers at the
%   solution (the multipliers corresponding to the lower bounds on the
%   variables, the upper bounds on the variables, and the equality and
%   inequality constraints), and the number of iterations needed to converge
%   to the stationary point (provided it has converged at all).
%
%   For more information, please consult: A. Wachter and L. T. Biegler. "On
%   the Implementation of a Primal-Dual Interior Point Filter Line Search
%   Algorithm for Large-Scale Nonlinear Programming." Mathematical
%   Programming 106(1), pp. 25-57, 2006.
%
%   Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
%   This code is published under the Common Public License.
%
%   Author: Peter Carbonetto
%           Dept. of Computer Science
%           University of British Columbia
%           October 13, 2007
