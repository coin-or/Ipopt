% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function [dR, dS] = computeJGGradient (qR, qS, auxdata)
  [K C f Rv Rf Sv Sf NS d] = deal(auxdata{:});
  
  nr = length(Rv);   % The number of large regions.
  ns = length(Sv);   % The number of small regions.

  % Reshape the input vectors.
  qR = reshapemarginals(qR,Rv,K);
  qS = reshapemarginals(qS,Sv,K);
  
  % Compute the gradient terms for the large regions. Repeat for each
  % large region.
  dR = cell(1,nr);
  for r = 1:nr
    is    = Rv{r};  % The variable nodes in the large region.
    table = 1 ./ qR{r};
    
    % Multiply the table by all the factors in the large region.
    for j = Rf{r}
      table = multiplyfactors(table,is,f{j},C{j});
    end
    
    dR{r} = -log(table) + 1;
  end
  
  % Compute the gradient terms for the small regions. Repeat for each
  % small region.
  dS = cell(1,ns);
  for s = 1:ns
    is    = Sv{s};  % The variable nodes in the small region.
    table = 1 ./ qS{s};

    % Multiply the table by all the factors in the small region.
    for j = Sf{s}
      table = multiplyfactors(table,is,f{j},C{j});
    end

    dS{s} = (1 - d(s)) * (-log(table) + 1);
  end
  
  % Convert the computed gradients into vectors.
  dR = vectorize(dR);
  dS = vectorize(dS);
  