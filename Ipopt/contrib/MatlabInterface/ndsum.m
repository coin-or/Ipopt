% NDSUM    Multi-dimensional summation.
% NDSUM(X,DIM) sums out the dimensions in DIM, and squeezes the result.
%
% (c) Microsoft Corporation. All rights reserved.
% Authors: Written by Tom Minka and modified by Peter Carbonetto.

function x = ndsum(x,dim)

  
  sz = size(x);
  for i = dim
    x = sum(x,i);
  end
  sz(dim) = [];
  x = reshape(x,[sz 1 1]);
