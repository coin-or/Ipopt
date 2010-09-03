function [Nodes,Elem,Phi,T,Cntrl] = LoadAirflowData(Datafilename)
% Copyright (C) 2010 International Business Machines.
% All Rights Reserved.
% This code is published under the Common Public License.
%
% Authors:  Johannes Huber     IBM        2010-09-03
  % read mesh data
  MeshFileName = 'MeshGen';
  tmp = importdata([MeshFileName '.node'],' ',2);
  Nodes = tmp.data(:,2:end);
  tmp = importdata([MeshFileName '.ele'],' ',2);
  Elem = tmp.data(:,2:end);

  % reade solution data
  SolutionFileName = [Datafilename 'State.dat'];
  tmp = importdata(SolutionFileName,'\t');
  Sol = reshape(tmp,[],size(Nodes,1))';
  Phi = Sol(:,1);
  if(size(Sol,2)>1)
    T = Sol(:,2);
  else
    T = [];
  end
  
  SolutionFileName = [Datafilename 'Cntrl.dat'];
  tmp = importdata(SolutionFileName,'\t');
  Cntrl = tmp.data(:,2);
end
