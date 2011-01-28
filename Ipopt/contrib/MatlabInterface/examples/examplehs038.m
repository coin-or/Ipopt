% Test the "ipopt" Matlab interface on the Hock & Schittkowski test problem
% #38. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
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
function [x, info] = examplehs038

  x0         = [-3  -1  -3  -1];   % The starting point.
  options.lb = [-10 -10 -10 -10];  % Lower bound on the variables.
  options.ub = [+10 +10 +10 +10];  % Upper bound on the variables.

  % The callback functions.
  funcs.objective        = @objective;
  funcs.gradient         = @gradient;
  funcs.hessian          = @hessian;
  funcs.hessianstructure = @hessianstructure;
  funcs.iterfunc         = @callback;

  % Set the IPOPT options.
  options.ipopt.mu_strategy = 'adaptive';
  options.ipopt.print_level = 0;
  options.ipopt.tol         = 1e-7;
  options.ipopt.max_iter    = 100;

  % Run IPOPT.
  [x info] = ipopt(x0,funcs,options);

% ----------------------------------------------------------------------
function f = objective (x)
  f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);

% ----------------------------------------------------------------------
function g = gradient (x)
  g(1) = -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
  g(2) = 200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1);
  g(3) = -360*x(3)*(x(4)-x(3)^2) -2*(1-x(3));
  g(4) = 180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1);
  
% ----------------------------------------------------------------------
function H = hessianstructure()
  H = sparse([ 1  0  0  0 
               1  1  0  0
               0  0  1  0
               0  1  1  1 ]);

% ----------------------------------------------------------------------
function H = hessian (x, sigma, lambda)
  H = [ 1200*x(1)^2-400*x(2)+2  0       0                          0
        -400*x(1)               220.2   0                          0
         0                      0       1080*x(3)^2- 360*x(4) + 2  0
         0                      19.8   -360*x(3)                   200.2 ];
  H = sparse(sigma*H);
  
% ----------------------------------------------------------------------
function b = callback (t, f, x)
  fprintf('%3d  %0.3g \n',t,f);
  b = true;