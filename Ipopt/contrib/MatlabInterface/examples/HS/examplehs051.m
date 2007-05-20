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

x0  = [ +2.5 +0.5 +2.0 -1.0 +0.5 ];
lb  = [-inf -inf -inf -inf -inf];
ub  = [+inf +inf +inf +inf +inf];
lbc = [4 0 0]; 
ubc = [4 0 0];

x = ipopt(x0,lb,ub,lbc,ubc,'computeObjectiveHS051','computeGradientHS051',...
	  'computeConstraintsHS051','computeJacobianHS051','',[],'',...
	  'hessian_approximation','limited-memory',...
	  'mu_strategy','adaptive','tol',1e-7);
