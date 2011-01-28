% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function H = computeJGHessian (qR, qS, sigma, lambda, ...
                               returnStructureOnly, auxdata)
  [K C f Rv Rf Sv Sf NS d] = deal(auxdata{:});
  
  ns  = length(Sv);  % The number of small regions.
  nqr = length(qR);  % The number of variables associated with the
                     % marginals defined on the large regions.
  nqs = length(qS);  % The number of variables associated with the
                     % marginals defined on the small regions.
  
  if returnStructureOnly
    H = speye(nqr + nqs);
  else

    % Compute the second-order partial derivatives on the subblock HR.
    HR = spdiag(1 ./ qR);
    
    % Compute the second-order partial derivatives on the subblock
    % HS. Repeat for each small region.
    entries = zeros(1,nqs);
    is      = 0;
    for s = 1:ns
      tablesize = prod(K(Sv{s}));  % The number of entries necessary to
                                   % store the entire marginal.
      is          = is(end) + (1:tablesize);
      entries(is) = (1-d(s)) ./ qS(is);
    end
    HS = spdiag(entries);
    
    % Construct the full Hessian from its blocks.
    H = sigma * [ HR                spzeros(nqr,nqs);
                  spzeros(nqs,nqr)  HS             ];
  end
  
% ------------------------------------------------------------------
function A = spdiag (xs)
  n = length(xs);  % The height and width of the sparse matrix.
  
  % Produce the sparse, diagonal matrix.
  A = sparse(1:n,1:n,xs,n,n);
  