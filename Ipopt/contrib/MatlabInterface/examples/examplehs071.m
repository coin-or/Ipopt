% Test the "ipopt" Matlab interface on the Hock & Schittkowski test problem
% #71. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
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
function [x, info] = examplehs071
  
  x0         = [1 5 5 1];  % The starting point.
  options.lb = [1 1 1 1];  % Lower bound on the variables.
  options.ub = [5 5 5 5];  % Upper bound on the variables.
  options.cl = [25  40];   % Lower bounds on the constraint functions.
  options.cu = [inf 40];   % Upper bounds on the constraint functions.
  
  % Initialize the dual point.
  options.zl     = [1 1 1 1];
  options.zu     = [1 1 1 1];
  options.lambda = [1 1];
  
  % Set the IPOPT options.
  options.ipopt.mu_strategy = 'adaptive';
  options.ipopt.tol         = 1e-7;
  
  % The callback functions.
  funcs.objective         = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
  funcs.constraints       = @(x) [ prod(x); sum(x.^2) ];
  funcs.gradient          = @gradient;
  funcs.jacobian          = @(x) sparse([ prod(x)./x; 2*x ]);
  funcs.jacobianstructure = @() sparse(ones(2,4));
  funcs.hessian           = @hessian;
  funcs.hessianstructure  = @() sparse(tril(ones(4)));
  
  % Run IPOPT.
  [x info] = ipopt(x0,funcs,options);

% ----------------------------------------------------------------------
function g = gradient (x)
  g = [ x(1)*x(4) + x(4)*sum(x(1:3))
        x(1)*x(4)
        x(1)*x(4) + 1
        x(1)*sum(x(1:3)) ]; 
  
% ----------------------------------------------------------------------
function H = hessian (x, sigma, lambda)
  
  H = sigma*[ 2*x(4)             0      0   0;
              x(4)               0      0   0;
              x(4)               0      0   0;
              2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ];
  
  H = H + lambda(1)*[    0          0         0         0;
                      x(3)*x(4)     0         0         0;
                      x(2)*x(4) x(1)*x(4)     0         0;
                      x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ];
  H = H + lambda(2)*diag([2 2 2 2]);
  H = sparse(H);
  