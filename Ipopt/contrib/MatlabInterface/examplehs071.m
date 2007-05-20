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

% The starting point.
x0  = [1 5 5 1];  % The starting point.
lb  = [1 1 1 1];  % Lower bound on the variables.
ub  = [5 5 5 5];  % Upper bound on the variables.
lbc = [25  40];   % Lower bounds on the constraint functions.
ubc = [inf 40];   % Upper bounds on the constraint functions.

x = ipopt(x0,lb,ub,lbc,ubc,'computeObjectiveHS071',...
	  'computeGradientHS071','computeConstraintsHS071',...
	  'computeJacobianHS071','computeHessianHS071',[],'',...
	  'mu_strategy','adaptive','tol',1e-7,'derivative_test',...
          'second-order','derivative_test_print_all','yes');
