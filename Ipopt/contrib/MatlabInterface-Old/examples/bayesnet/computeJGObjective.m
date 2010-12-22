% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function F = computeJGObjective (qR, qS, auxdata)
  [K C f Rv Rf Sv Sf NS d] = deal(auxdata{:});
    
  % Reshape the input vectors.
  qR = reshapemarginals(qR,Rv,K);
  qS = reshapemarginals(qS,Sv,K);

  % The objective function is the junction graph approximation to the
  % variational free energy.
  F = computefreeenergy(C,f,Rv,Rf,Sv,Sf,qR,qS,d);
  
% ------------------------------------------------------------------
function F = computefreeenergy (C, f, Rv, Rf, Sv, Sf, qR, qS, d)
  nr = length(Rv);   % The number of large regions.
  ns = length(Sv);   % The number of small regions.
  F  = 0;            % The return value.
    
  % Compute the average energy of each large region. Repeat for each
  % large region, then for each factor within the large region.
  for r = 1:nr
    for j = Rf{r}
      table = multiplyfactors(qR{r},Rv{r},log(f{j}),C{j});
      F     = F - sum(table(:));
    end
  end
  
  % Compute the average energy of each small region. Repeat for each
  % small region, then for each factor within the small region.
  for s = 1:ns
    for j = Sf{s}
      table = multiplyfactors(qS{s},Sv{s},log(f{j}),C{j});
      F     = F - (1 - d(s)) * sum(table(:));
    end
  end
  
  % Compute the entropy of each large region.
  for r = 1:nr
    table = qR{r} .* log(qR{r});
    F     = F + sum(table(:));
  end
  
  % Compute the entropy of each small region.
  for s = 1:ns
    table = qS{s} .* log(qS{s});
    F     = F + (1 - d(s)) * sum(table(:));
  end
