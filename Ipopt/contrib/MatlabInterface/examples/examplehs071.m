% Test the "ipopt" Matlab interface on the Hock & Schittkowski test problem
% #71. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
% Nonlinear Programming Codes. Lecture Notes in Economics and Mathematical
% Systems Vol. 187, Springer-Verlag.
%
% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007
function [x, multipliers] = examplehs071
  
% The starting point.
x0  = [1 5 5 1];  % The starting point.
lb  = [1 1 1 1];  % Lower bound on the variables.
ub  = [5 5 5 5];  % Upper bound on the variables.
lbc = [25  40];   % Lower bounds on the constraint functions.
ubc = [inf 40];   % Upper bounds on the constraint functions.

% Initial dual point.
multipliers.zl     = [1 1 1 1];
multipliers.zu     = [1 1 1 1];
multipliers.lambda = [1 1];
  
[x multipliers] = ipopt(x0,lb,ub,lbc,ubc,@computeObjective,...
			@computeGradient,@computeConstraints,...
			@computeJacobian,@computeHessian,[],'',...
			multipliers,'mu_strategy','adaptive','tol',1e-7);

% ----------------------------------------------------------------------
function f = computeObjective (x)

  f = x(1)*x(4)*sum(x(1:3)) + x(3);

% ----------------------------------------------------------------------
function c = computeConstraints (x) 
  
  c(1) = prod(x);
  c(2) = sum(x.^2);

% ----------------------------------------------------------------------
function g = computeGradient (x)
  
  g(1) = x(1)*x(4) + x(4)*sum(x(1:3));
  g(2) = x(1)*x(4);
  g(3) = x(1)*x(4) + 1;
  g(4) = x(1)*sum(x(1:3));

% ----------------------------------------------------------------------
function J = computeJacobian (x, returnStructureOnly)

  if returnStructureOnly
    J = sparse(ones(2,4));
  else
    J = sparse([prod(x) ./ x; 2*x]);
  end
  
% ----------------------------------------------------------------------
function H = computeHessian (x, sigma, lambda, returnStructureOnly)

  if returnStructureOnly
    H = sparse(tril(ones(4)));
  else
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
  end
  