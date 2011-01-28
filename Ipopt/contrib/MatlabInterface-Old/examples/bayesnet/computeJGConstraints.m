% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function g = computeJGConstraints (qR, qS, auxdata) 
  [K C f Rv Rf Sv Sf NS d] = deal(auxdata{:});
  
  nr = length(Rv);  % The number of large regions.
  ns = length(Sv);  % The number of small regions.

  % Reshape the input vectors.
  qR = reshapemarginals(qR,Rv,K);
  qS = reshapemarginals(qS,Sv,K);

  % Compute the normalization constraints for the large regions. Repeat
  % for each large region.
  gR = zeros(nr,1);
  for r = 1:nr
    gR(r) = sum(qR{r}(:));
  end
  
  % Compute the normalization constraints for the small regions. Repeat
  % for each small region.
  gS = zeros(ns,1);
  for s = 1:ns
    gS(s) = sum(qS{s}(:));
  end
  
  % These are the consistency constraints. There is a constraint for
  % every small region S, then again for every neighbouring large
  % region R, and then for every table entry of the marginal qS(xS). 
  % These constraints are stored in the cell array "cRS". It has the 
  % same structure as the collection of sum-product messages "mRS" in
  % the MATLAB function "bp".
  cRS = cell(1,ns);
  for s = 1:ns
    is = Sv{s};  % The variable nodes in the small region.
    rs = NS{s};  % The set of neighbouring large regions.
    for i = 1:length(rs)
      r         = rs(i);
      cRS{s}{i} = margfactor(qR{r},Rv{r},is) - qS{s};
    end
  end
  
  g = [gR; gS; vectorize(cRS)];
  