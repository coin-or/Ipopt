% Test the "ipopt" Matlab interface on the Hock & Schittkowski test problem
% #51. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
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
function x = examplehs051

x0  = [ +2.5 +0.5 +2.0 -1.0 +0.5 ];
lb  = [-inf -inf -inf -inf -inf];
ub  = [+inf +inf +inf +inf +inf];
lbc = [4 0 0]; 
ubc = [4 0 0];

x = ipopt(x0,lb,ub,lbc,ubc,@computeObjective,@computeGradient,...
	  @computeConstraints,@computeJacobian,'',[],'',[],...
	  'jac_c_constant','yes','hessian_approximation','limited-memory',...
	  'mu_strategy','adaptive','tol',1e-7);

% ----------------------------------------------------------------------
function f = computeObjective (x)

  f = (x(1) - x(2))^2 + ...
      (x(2) + x(3) - 2)^2 + ...
      (x(4) - 1)^2 + (x(5) - 1)^2;
  
% ----------------------------------------------------------------------
function g = computeGradient (x)
  
  g = 2*[ x(1) - x(2);
	  x(2) + x(3) - 2 - x(1) + x(2);
	  x(2) + x(3) - 2;
	  x(4) - 1;
	  x(5) - 1 ];

% ----------------------------------------------------------------------
function c = computeConstraints (x)

  c = [ x(1) + 3*x(2);
        x(3) + x(4) - 2*x(5);
        x(2) - x(5) ];

% ----------------------------------------------------------------------
function J = computeJacobian (x, returnStructureOnly)
  
  J = sparse([ 1  3  0  0  0;
	       0  0  1  1 -2;
	       0  1  0  0 -1 ]);
  