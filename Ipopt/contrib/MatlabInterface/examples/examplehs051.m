% Test the "ipopt" Matlab interface on the Hock & Schittkowski test problem
% #51. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
% Nonlinear Programming Codes. Lecture Notes in Economics and Mathematical
% Systems Vol. 187, Springer-Verlag.
%
% Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         September 18, 2008
function [x, info] = examplehs051

  x0         = [ 2.5 0.5 2 -1 0.5 ];  % The starting point.
  options.cl = [ 4 0 0 ];             % Lower bounds on constraints.
  options.cu = [ 4 0 0 ];             % Upper bounds on constraints.

  % Set the IPOPT options.
  options.ipopt.jac_c_constant        = 'yes';
  options.ipopt.hessian_approximation = 'limited-memory';
  options.ipopt.mu_strategy           = 'adaptive';
  options.ipopt.tol                   = 1e-7;

  % The callback functions.
  funcs.objective         = @objective;
  funcs.constraints       = @constraints;
  funcs.gradient          = @gradient;
  funcs.jacobian          = @jacobian;
  funcs.jacobianstructure = @jacobian;
  
  % Run IPOPT.
  [x info] = ipopt(x0,funcs,options);
  
% ----------------------------------------------------------------------
function f = objective (x)
  f = (x(1) - x(2))^2 + ...
      (x(2) + x(3) - 2)^2 + ...
      (x(4) - 1)^2 + (x(5) - 1)^2;

% ----------------------------------------------------------------------
function g = gradient (x)
  g = 2*[ x(1) - x(2);
	  x(2) + x(3) - 2 - x(1) + x(2);
	  x(2) + x(3) - 2;
	  x(4) - 1;
	  x(5) - 1 ];

% ----------------------------------------------------------------------
function c = constraints (x)
  c = [ x(1) + 3*x(2);
        x(3) + x(4) - 2*x(5);
        x(2) - x(5) ];

% ----------------------------------------------------------------------
function J = jacobian (x)  
  J = sparse([ 1  3  0  0  0;
	       0  0  1  1 -2;
	       0  1  0  0 -1 ]);
  