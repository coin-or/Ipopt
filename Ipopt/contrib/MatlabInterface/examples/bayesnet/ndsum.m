function X = ndsum (X, dims)
  
  % Sum over requested dimensions.
  s = size(X);
  for i = dims
    X = sum(X,i);
  end
  
  % Remove singleton dimensions.
  s(dims) = [];
  X       = reshape(X,[s 1 1]);
