% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function J = computeJGJacobian (qR, qS, returnStructureOnly, auxdata)
  [K C f Rv Rf Sv Sf NS d] = deal(auxdata{:});
  
  nr  = length(Rv);  % The number of large regions.
  ns  = length(Sv);  % The number of small regions.
  nqr = length(qR);  % The number of variables associated with the
                     % marginals defined on the large regions.
  nqs = length(qS);  % The number of variables associated with the
                     % marginals defined on the small regions.

  % Reshape the input vectors.
  qR = reshapemarginals(qR,Rv,K);
  qS = reshapemarginals(qS,Sv,K);
  
  % This block is the part of the Jacobian that represents the 
  % first-order partial derivatives of the normalization constraints on
  % the large regions. We have one constraint function for every small
  % region. 
  rows = zeros(1,nqr);
  is   = 0;
  for r = 1:nr
    tablesize = numel(qR{r});
    is        = is(end) + (1:tablesize);
    rows(is)  = r;
  end
  JR = sparse(rows,1:nqr,ones(1,nqr),nr,nqr);
  
  % This block is the part of the Jacobian that represents the
  % first-order partial derivatives of the normalization constraints on
  % the small regions. We have one constraint function for every small
  % region. 
  rows = zeros(1,nqs);
  is   = 0;
  for s = 1:ns
    tablesize = numel(qS{s});
    is        = is(end) + (1:tablesize);
    rows(is)  = s;
  end
  JS = sparse(rows,1:nqs,ones(1,nqs),ns,nqs);
  
  % This block is the part of the Jacobian that represents the
  % first-order partial derivatives of the consistency constraints with
  % respect to the variables qS(xS). What we will do is form subblocks of
  % this block. We have one constraint function for each small region,
  % then again for each neighbouring large region, and then for every
  % possible configuration on the small region.
  JCS = cell(1,ns);
  for s = 1:ns
    w      = numel(qS{s});  % The width of the Jacobian subblock.
    h      = d(s) * w;      % The height of the Jacobian subblock.
    cols   = repmat(1:w,1,d(s));
    JCS{s} = sparse(1:h,cols,-ones(1,h),h,w);
  end

  % This block is the part of the Jacobian that represents the
  % first-order partial derivatives of the consistency constraints with
  % respect to the variables qR(xR). What we will do is form subblocks of
  % this block. We have one constraint function for each small region,
  % then again for each neighbouring large region, and then for every
  % possible configuration on the small region.
  JCR = cell(sum(d),1);
  t   = 1;  % The row of the current subblock.
  for s = 1:ns
    for r = NS{s}
      h      = numel(qS{s});  % The height of the Jacobian subblock.
      w      = numel(qR{r});  % The width of the Jacobian subblock.
      rows   = multiplyfactors(ones([K(Rv{r}) 1]),Rv{r},...
                               reshape(1:h,[K(Sv{s}) 1]),Sv{s});
      JCR{t} = sparse(rows(:)',1:w,ones(1,w),h,w);
      t      = t + 1;         % Move to the next row.
    end
  end

  % Construct the full Jacobian from its blocks.
  J = [ JR                              spzeros(nr,nqs);
        spzeros(ns,nqr)                 JS;
        spblockrows(JCR,vectorize(NS))  spblockdiag(JCS) ];
  
% ------------------------------------------------------------------
function A = spblockdiag (blocks)
  n       = length(blocks);  % The number of blocks.
  h       = 0;               % The height of the full matrix.
  w       = 0;               % The width of the full matrix.
  rows    = [];
  cols    = [];
  entries = [];
  
  % Repeat for each block.
  for i = 1:n
    [hb wb]    = size(blocks{i});  % The height and width of the block.
    [is js xs] = find(blocks{i});
    rows       = [rows h + is'];
    cols       = [cols w + js'];
    entries    = [entries xs'];
    h          = h + hb;
    w          = w + wb;
  end
  
  % Construct the full matrix from the diagonal subblocks.
  A = sparse(rows,cols,entries,h,w);
  
% ------------------------------------------------------------------
function A = spblockrows (blocks, c)
  n       = length(blocks);  % The number of blocks.
  nc      = max(c);          % The number of block columns.
  h       = 0;               % The height of the full matrix.
  rows    = [];              % The rows of the individual entries.
  cols    = [];              % The columns of the individual entries.
  entries = [];
  
  % Compute the width of every subblock.
  bw = zeros(1,nc);
  for i = 1:n
    bw(c(i)) = size(blocks{i},2);
  end
  
  % Compute the starting column index (minus 1) for every subblock. Also,
  % compute the width of the full matrix, w.
  bw = cumsum(bw);
  w  = bw(nc);
  bw = [0 bw(1:nc-1)];
  
  % Repeat for each block.
  for i = 1:n
    hb         = size(blocks{i},1);  % The height of the block.
    [is js xs] = find(blocks{i});
    rows       = [rows h + is'];
    cols       = [cols bw(c(i)) + js'];
    entries    = [entries xs'];
    h          = h + hb;
  end
  
  % Construct the full matrix from the subblocks.
  A = sparse(rows,cols,entries,h,w);
  