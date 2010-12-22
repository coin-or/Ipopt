% Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         May 19, 2007

function qR = bopt (K, C, f, Rv, Rf, Sv, Sf, NS, verbose)
  maxiter   = 100;         % The maximum number of iterations.
  tolerance = 1e-6;        % The convergence criterion.
  nr        = length(Rv);  % The number of large regions.
  ns        = length(Sv);  % The number of small regions.

  % Convert the input K to a row vector.
  K = K(:)';
  
  % Get the degrees of the separators.
  d = degrees(NS);

  % Initialize the marginals on the large and small regions to uniform
  % probability tables.
  qR = initmarginals(K,Rv);
  qS = initmarginals(K,Sv);
  
  % Compute the total number consistency constaints. We have one
  % constraint for every pair (R,S), then again for every possible
  % configuration xS.
  nc = numconsistencyconstraints(qS,d);
  
  % Run the IPOPT solver.
  qR             = vectorize(qR);
  qS             = vectorize(qS);
  nqr            = length(qR);
  nqs            = length(qS);
  [status qR qS] = ipopt({qR qS},{ repmat(eps,nqr,1) repmat(eps,nqs,1) },...
                  { repmat(inf,nqr,1) repmat(inf,nqs,1) },...
                  [ ones(1,nr) ones(1,ns) zeros(1,nc) ],...
                  [ ones(1,nr) ones(1,ns) zeros(1,nc) ],...
                  @computeJGObjective,@computeJGGradient,...
                  @computeJGConstraints,@computeJGJacobian,...
                  @computeJGHessian,{ K C f Rv Rf Sv Sf NS d },'',...
                  [],'mu_strategy','adaptive','max_iter',maxiter,...
                  'tol',tolerance,'jac_c_constant','yes',...
		  'jac_d_constant','yes','print_level',verbose*5);

  % Reshape the solution.
  qR = reshapemarginals(qR,Rv,K);
  qS = reshapemarginals(qS,Sv,K);

% ------------------------------------------------------------------
function d = degrees (NS)
  ns = length(NS);  % The number of separators.
  d  = zeros(1,ns);
  
  % Repeat for each separator.
  for s = 1:ns
    d(s) = length(NS{s});
  end

% ----------------------------------------------------------------
function qR = initmarginals (K, Rv)
  nr = length(Rv);  % The number of regions.
  qR = cell(1,nr);  % The return value.
  
  % Repeat for each region.
  for r = 1:nr
    is    = Rv{r};  % The variables nodes in the large region.
    table = ones([K(is) 1]);
    qR{r} = table / sum(table(:));
  end

% ----------------------------------------------------------------
function nc = numconsistencyconstraints (qS, d)
  ns = length(qS);  % The number of small regions.
  nc = 0;           % The return value.
  
  % Repeat for each small region.
  for s = 1:ns
    tablesize = numel(qS{s});
    nc        = nc + d(s) * tablesize;
  end
