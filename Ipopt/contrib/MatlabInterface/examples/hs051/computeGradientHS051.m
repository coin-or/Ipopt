% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function grad = computeGradientHS051 (x)
  
  grad = 2*[x(1) - x(2);
            x(2) + x(3) - 2 - x(1) + x(2);
            x(2) + x(3) - 2;
            x(4) - 1;
            x(5) - 1 ];
    