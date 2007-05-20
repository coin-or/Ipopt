% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Common Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function J = computeJacobianHS071 (x, returnStructureOnly)

  if returnStructureOnly
    J = sparse(ones(2,4));
  else
    
    % Compute the Jacobian.
    J = sparse([prod(x) ./ x; 2*x]);
  end