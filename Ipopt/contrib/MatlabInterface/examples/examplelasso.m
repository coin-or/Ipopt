% In this small MATLAB script, we compute the least squares solution to a
% regression problem subject to L1 regularization, which rewards "sparse"
% models that have regression coefficients of zero. See, for instance,
% the work in the "Lasso" by the statistician Robert Tibshirani.
%
% Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         September 18, 2008

% Experiment parameters.
lambda = 1;                      % Level of L1 regularization.
n      = 100;                    % Number of training examples.
e      = 1;                      % Std. dev. in noise of outputs.
beta   = [ 0 0 2 -4 0 0 -1 3 ]'; % "True" regression coefficients.

% Set the random number generator seed.
seed = 7;
rand('state',seed);
randn('state',seed);

% CREATE DATA SET.
% Generate the input vectors from the standard normal, and generate the
% binary responses from the regression with some additional noise, and then
% transform the results using the logistic function. The variable "beta" is
% the set of true regression coefficients.
m     = length(beta);      % Number of features.
A     = randn(n,m);        % The n x m matrix of examples.
noise = e * randn(n,1);    % Noise in outputs.
y     = A * beta + noise;  % The binary outputs.

% COMPUTE SOLUTION WITH IPOPT.
% Compute the L1-regularized maximum likelihood estimator.
w = lasso(A,y,lambda);
fprintf('Solution:\n');
disp(w);


