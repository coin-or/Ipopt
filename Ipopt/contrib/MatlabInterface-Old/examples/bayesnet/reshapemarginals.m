% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function qR = reshapemarginals (q, R, K)
  nr  = length(R);   % The number of regions.
  qR  = cell(1,nr);  % The cell array to return.
  pos = 0;           % Current entry of the vector q.
  
  % Repeat for each region.
  for r = 1:nr
    is        = R{r};         % The variable nodes in region r.
    tablesize = prod(K(is));  % The size of the table.
    table     = q(pos + (1:tablesize));
    qR{r}     = reshape(table,[K(is) 1]);
    
    % Move to the next location in the table.
    pos = pos + tablesize;
  end
  