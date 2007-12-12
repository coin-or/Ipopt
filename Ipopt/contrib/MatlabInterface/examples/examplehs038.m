% Test the "ipopt" Matlab interface on the Hock & Schittkowski test problem
% #38. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
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
function x = examplehs038

% The starting point.
x0 = [-3  -1  -3  -1];   % The starting point.
lb = [-10 -10 -10 -10];  % Lower bound on the variables.
ub = [+10 +10 +10 +10];  % Upper bound on the variables.

x = ipopt(x0,lb,ub,[],[],@computeObjective,@computeGradient,'','',...
	  @computeHessian,[],@callback,[],'mu_strategy','adaptive',...
	  'tol',1e-7,'max_iter',100);

% ----------------------------------------------------------------------
function f = computeObjective (x)

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);
  
  f = 100*(x2-x1^2)^2 + (1-x1)^2 + 90*(x4-x3^2)^2 + (1-x3)^2 + ...
      10.1*(x2-1)^2 + 10.1*(x4-1)^2 + 19.8*(x2-1)*(x4-1);

% ----------------------------------------------------------------------
function g = computeGradient (x)

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);

  g(1) = -400*x1*(x2-x1^2) - 2*(1-x1);
  g(2) = 200*(x2-x1^2) + 20.2*(x2-1) + 19.8*(x4-1);
  g(3) = -360*x3*(x4-x3^2) -2*(1-x3);
  g(4) = 180*(x4-x3^2) + 20.2*(x4-1) + 19.8*(x2-1);
  
% ----------------------------------------------------------------------
function H = computeHessian (x, sigma, lambda, returnStructureOnly)

  if returnStructureOnly
    H = sparse([1  0  0  0;
                1  1  0  0;
                0  0  1  0;
                0  1  1  1]);
  else
    
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);

    H = sigma*[ 1200*x1^2 - 400*x2 + 2  0      0                       0;
                -400*x1                 220.2  0                       0;
                0                       0      1080*x3^2 - 360*x4 + 2  0;
                0                       19.8   -360*x3                 200.2];
    H = sparse(H);
  end
  
% ----------------------------------------------------------------------
function callback (t, f, x)
  fprintf('%3d  %0.3g \n', t, f);
  