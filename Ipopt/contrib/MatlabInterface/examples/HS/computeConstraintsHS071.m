% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function g = computeConstraintsHS071 (x) 
  
  g(1) = prod(x);
  g(2) = sum(x.^2);
