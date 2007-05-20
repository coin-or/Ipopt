% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function H = computeHessianHS071 (x, sigma, lambda, returnStructureOnly)

  if returnStructureOnly
    H = sparse(tril(ones(4)));
  else
    
    % Compute the Hessian.
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
  