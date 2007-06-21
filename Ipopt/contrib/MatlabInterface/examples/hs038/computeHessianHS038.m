% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function H = computeHessianHS038 (x, sigma, lambda, returnStructureOnly)

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

    % Compute the Hessian.
    H = sigma*[ 1200*x1^2 - 400*x2 + 2  0      0                       0;
                -400*x1                 220.2  0                       0;
                0                       0      1080*x3^2 - 360*x4 + 2  0;
                0                       19.8   -360*x3                 200.2];
    H = sparse(H);
  end
  