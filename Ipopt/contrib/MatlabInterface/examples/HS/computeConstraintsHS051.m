% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function g = computeConstraintsHS051 (x)

  g = [x(1) + 3*x(2);
       x(3) + x(4) - 2*x(5);
       x(2) - x(5)];
  