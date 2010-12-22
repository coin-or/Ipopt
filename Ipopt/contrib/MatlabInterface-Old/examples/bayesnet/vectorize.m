% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function v = vectorize (cells)
  n = length(cells);  % The number of cells in the cell array.
  v = [];             % The return value.
  
  % Repeat for each cell.
  for i = 1:n
    c = cells{i};
    if iscell(c)
      c = vectorize(c);
    end
    v = [v; c(:)];
  end
  