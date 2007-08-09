% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function f = computeObjectiveHS071 (x)
  f = x(1)*x(4)*sum(x(1:3)) + x(3);
