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


% The starting point.
x0 = [-3  -1  -3  -1];   % The starting point.
lb = [-10 -10 -10 -10];  % Lower bound on the variables.
ub = [+10 +10 +10 +10];  % Upper bound on the variables.

x = ipopt(x0,lb,ub,[],[],'computeObjectiveHS038',...
	  'computeGradientHS038','','','computeHessianHS038',[],...
	  'genericcallback','mu_strategy','adaptive','tol',1e-7,...
	  'max_iter',100);

