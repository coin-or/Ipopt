function ShowResult(Nodes,Elem,Phi,n,startpts)
  % Copyright (C) 2010 International Business Machines.
  % All Rights Reserved.
  % This code is published under the Common Public License.
  %
  % $Id: ShowSolution.m 29 2010-08-31 18:19:25Z jahuber $
  %
  % Authors:  Johannes Huber     IBM        2010-09-03
  dim = size(Nodes,2);
  if ~exist('n','var')
    n = 20*ones(1,dim);
  elseif numel(n)==1
    n = n*ones(1,dim);
  end
  
  [Center,dPhi] = trigrad(Nodes,Elem,Phi);
  if dim==3
    quiver3(Center(:,1),Center(:,2),Center(:,3),dPhi(:,1),dPhi(:,2),dPhi(:,3));
  else
    quiver(Center(:,1),Center(:,2),dPhi(:,1),dPhi(:,2));
  end
  
  Min = min(Center);
  Max = max(Center);
  t1 = linspace(Min(1),Max(1),n(1));
  t1 = t1(2:end-1);
  t2 = linspace(Min(2),Max(2),n(2));
  t2 = t2(2:end-1);
  if dim==2
    [X,Y] = meshgrid(t1,t2);
    RegGrid = [X(:) Y(:)];
  else
    t3 = linspace(Min(3),Max(3),n(3));
    t3 = t3(3:end-1);
    [X,Y,Z] = meshgrid(t1,t2,t3);
    RegGrid = [X(:) Y(:) Z(:)];
  end
  
  u = griddatan(Center,dPhi(:,1),RegGrid);
  u = reshape(u,n-2);
  v = griddatan(Center,dPhi(:,2),RegGrid);
  v = reshape(v,n-2);
  if dim==3
    w = griddatan(Center,dPhi(:,3),RegGrid);
    w = reshape(w,n-2);
  end

%   t2 = linspace(0.11,0.29,3);
%   t3 = linspace(0.21,0.39,3);
%   [startx, starty, startz] = meshgrid(0.2,t2,t3);
  if(dim==2)
    verts = stream2(X,Y,u,v,startpts(:,1),startpts(:,2));
    %streamline(X,Y,u,v,startpts(:,1),startpts(:,2));
    iverts = interpstreamspeed(X,Y,u,v,verts,.0025);
  else
    %streamline(X,Y,Z,u,v,w,startpts(:,1),startpts(:,2),startpts(:,3));
    verts = stream2(X,Y,Z,u,v,w,startpts(:,1),startpts(:,2),startpts(:,3));
    iverts = interpstreamspeed(X,Y,Z,u,v,w,verts,.0025);
  end
  streamline(verts);
%  streamparticles(iverts,35,'animate',10,'ParticleAlignment','on')
end
