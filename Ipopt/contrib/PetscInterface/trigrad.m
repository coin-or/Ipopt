function [Center,Grads] = trigrad(Nodes,E,Val)
  % Copyright (C) 2010 International Business Machines.
  % All Rights Reserved.
  % This code is published under the Common Public License.
  %
  % $Id: tri.m 29 2010-08-31 18:19:25Z jahuber $
  %
  % Authors:  Johannes Huber     IBM        2010-09-03
  dim = size(Nodes,2);
  % assume linear elements
  Center = zeros(size(E,1),dim);
  Grads = zeros(size(E,1),dim);
  GradRef = [-ones(dim,1) eye(dim)];  % Gradients of the shape functions on reference element
  for iEl=1:size(E,1)
    Center(iEl,:) = mean(Nodes(E(iEl,:)',:));
    J = Nodes(E(iEl,2:end),:)-repmat(Nodes(E(iEl,1),:),dim,1);
    grad = J\GradRef;
    grad = grad*Val(E(iEl,:));
    Grads(iEl,:) = grad';
  end
end
