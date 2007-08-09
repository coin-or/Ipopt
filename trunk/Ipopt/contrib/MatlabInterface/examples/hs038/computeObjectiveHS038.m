% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function f = computeObjectiveHS038 (x)
  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);
  
  f = 100*(x2-x1^2)^2 + (1-x1)^2 + 90*(x4-x3^2)^2 + (1-x3)^2 + ...
      10.1*(x2-1)^2 + 10.1*(x4-1)^2 + 19.8*(x2-1)*(x4-1);

