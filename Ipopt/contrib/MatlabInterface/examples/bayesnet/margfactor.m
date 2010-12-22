% MARGFACTOR(F,R,S) takes as input a factor F defined on set R, and
% another set S which must be a (non-strict and non-empty) subset of
% R. The return value is a factor defined on the set S.
%
% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function f = margfactor (f, R, S)
  nr = length(R);  % The size of set R.
  ns = length(S);  % The size of set S.
  
  % Get the intersection of sets S and R.
  [ans RI SI] = intersect(R,S);
  
  % Get the set of elements that belong to R but don't belong to S.
  [ans RD] = setdiff(R,S);
  
  % Reorder the factor table so that the elements belonging to the
  % intersection of R and S supersede the elements that are not in the
  % intersection.
  if nr > 1
    f = permute(f,[RI RD]);
  end
  
  % Sum over the dimensions in R - S.
  if nr > ns
    f = ndsum(f,ns+1:nr);
  end
  
  if ns > 1
    
    % This is the reverse mapping to recover a factor table on S.
    SR     = zeros(1,ns);
    SR(SI) = 1:ns;
    
    % Reorder the final result.
    f = permute(f,SR);
  end