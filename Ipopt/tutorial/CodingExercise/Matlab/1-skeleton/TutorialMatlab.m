% Copyright (C) 2009 International Business Machines 
% All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% $Id: hs071_c.c 699 2006-04-05 21:05:18Z andreasw $
%
% Author:  Andreas Waechter               IBM    2009-04-02
%
% This file is part of the Ipopt tutorial.  It is the skeleton for
% the matlab implemention of the coding exercise problem (in AMPL
% formulation):
%
% param n := 4;
%
% var x {1..n} <= 0, >= -1.5, := -0.5;
%
% minimize obj:
%   sum{i in 1..n} (x[i]-1)^2;
%
% subject to constr {i in 2..n-1}:
%   (x[i]^2+1.5*x[i]-i/n)*cos(x[i+1]) - x[i-1] = 0;
%
% The constant term "i/n" in the constraint is supposed to be input data

function [x, info] = TutorialMatlab

  % Size of the problem
  n = 5;

  % Problem data
  a = ((1:n)/n)';

  % Starting point
  x0 = -0.5*ones(n,1);

  % Lower and upper bounds for the variables
  options.lb = -1.5*ones(n,1);
  options.ub = ...

  % Constraint bounds
  options.cl = ...
  options.cu = ...

  % Set the Ipopt options
  % options.ipopt.mu_strategy = 'adaptive';
  options.ipopt.derivative_test = 'first-order';
  % options.ipopt.derivative_test = 'second-order';

  % Set the callback functions
  funcs.objective         = @eval_f;
  funcs.constraints       = @eval_g;
  funcs.gradient          = @eval_grad_f;
  funcs.jacobian          = @eval_jac_g;
  funcs.jacobianstructure = @eval_jac_g_struct;
  funcs.hessian           = @eval_hess;
  funcs.hessianstructure  = @eval_hess_struct;

  [x info] = ipopt(x0, funcs, options);

  % End of main function
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Evaluate value of objective function
  function f = eval_f(x)

    

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Evaluate gradient of objective function
  function df = eval_grad_f(x)

    

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Evaluate value of constraint bodies
  function g = eval_g(x)

    

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Return constraint Jacobian strcture
  function A = eval_jac_g_struct


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Evaluate constraint Jacobian
  function A = eval_jac_g(x)



  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Return Hessian of Lagrangian function structure
  function H = eval_hess_struct

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Evaluate Hessian of Lagrangian function
  function H = eval_hess(x, sigma, lambda)

  end

end
