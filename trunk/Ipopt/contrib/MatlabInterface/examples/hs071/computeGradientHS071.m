% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function grad = computeGradientHS071 (x)
  
  grad(1) = x(1)*x(4) + x(4)*sum(x(1:3));
  grad(2) = x(1)*x(4);
  grad(3) = x(1)*x(4) + 1;
  grad(4) = x(1)*sum(x(1:3));

  
