% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function grad = computeGradientHS038 (x)
  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);

  grad(1) = -400*x1*(x2-x1^2) - 2*(1-x1);
  grad(2) = 200*(x2-x1^2) + 20.2*(x2-1) + 19.8*(x4-1);
  grad(3) = -360*x3*(x4-x3^2) -2*(1-x3);
  grad(4) = 180*(x4-x3^2) + 20.2*(x4-1) + 19.8*(x2-1);

  