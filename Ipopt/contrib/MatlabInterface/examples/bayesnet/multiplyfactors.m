% MULTIPLYFACTORS(FR,R,FS,S) returns the pointwise product of factors FR
% and FS defined on sets R and S, respectively. We require that S be a
% (non-strict and non-empty) subset of R. The resulting factor is defined
% on set R. 
%
% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function fR = multiplyfactors (fR, R, fS, S)
  nr = length(R);  % The size of set R.
  ns = length(S);  % The size of set S.
  K  = size(fR);   % The number of possible discrete assignments to each 
                   % random variable in R.

  if nr == 1
    fR = fR .* fS;
  else
  
    % Get the intersection of sets S and R.
    [ans RI SI] = intersect(R,S);
    
    % Get the set of elements that belong to R but don't belong to C.
    [ans RD] = setdiff(R,S);
    
    % Check whether or not the intersection is the same size as S. If not,
    % it means that R is not a superset of S.
    if length(SI) ~= length(S)
      error('R is not a superset of S');
    end

    % Reorder the factor tables so that the elements belonging to the
    % intersection of R and S supersede the elements that are not in the
    % intersection.
    RT = [RI RD];
    fR = permute(fR,RT);
    if ns > 1
      fS = permute(fS,SI);
    end

    % This is the reverse mapping to recover the original factor table.
    RR     = zeros(1,nr);
    RR(RT) = 1:nr;
    
    % Expand the table fS so that it is of the same size as fR.
    fS = repmat(fS,[ones(1,ns) K(RD)]);
    
    % Multiply fR with the expanded version of fS.
    fR = fR .* fS;
  
    % Finally, we need to return the factor table fR back to its original 
    % ordering.
    fR = permute(fR,RR);
  end
  