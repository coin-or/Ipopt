% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function J = computeJacobianHS051 (x, returnStructureOnly)
  
  J = sparse([ 1  3  0  0  0;
	       0  0  1  1 -2;
	       0  1  0  0 -1 ]);
