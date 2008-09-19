function w = lasso (A, y, lambda)
  
  % Get the number of examples (n) and the number of regression
  % coefficients (m).
  [n m] = size(A);
    
  % The starting point.
  x0 = { zeros(m,1) 
         ones(m,1) };         
  
  % The constraint functions are bounded from below by zero.
  options.cl = zeros(2*m,1);
  options.cu = repmat(Inf,2*m,1);
  
  % Set up the auxiliary data.
  options.auxdata = { n m A y lambda };
  
  % Set the IPOPT options.
  options.ipopt.jac_d_constant   = 'yes';
  options.ipopt.hessian_constant = 'yes';
  options.ipopt.mu_strategy      = 'adaptive';
  options.ipopt.max_iter         = 100;
  options.ipopt.tol              = 1e-8;
  
  % The callback functions.
  funcs.objective         = @objective;
  funcs.constraints       = @constraints;
  funcs.gradient          = @gradient;
  funcs.jacobian          = @jacobian;
  funcs.jacobianstructure = @jacobianstructure;
  funcs.hessian           = @hessian;
  funcs.hessianstructure  = @hessianstructure;
  
  % Run IPOPT.
  [x info] = ipopt(x0,funcs,options);
  w        = x{1};
  
% ------------------------------------------------------------------
function f = objective (x, auxdata)
  [n m A y lambda] = deal(auxdata{:});
  [w u] = deal(x{:});
  f     = norm(y - A*w)^2/2 + lambda*sum(u);
  
% ------------------------------------------------------------------
function c = constraints (x, auxdata)
  [w u] = deal(x{:});
  c     = [ w + u; u - w ];
  
% ------------------------------------------------------------------
function g = gradient (x, auxdata)
  [n m A y lambda] = deal(auxdata{:});
  w = x{1};
  g = { -A'*(y - A*w) 
        repmat(lambda,m,1) };
  
% ------------------------------------------------------------------
function J = jacobianstructure (auxdata)  
  m = auxdata{2};
  I = speye(m);
  J = [ I I
        I I ];
  
% ------------------------------------------------------------------
function J = jacobian (x, auxdata)  
  m = auxdata{2};
  I = speye(m);
  J = [  I  I
        -I  I ];
  
% ------------------------------------------------------------------
function H = hessianstructure (auxdata)
  m = auxdata{2};
  H = [ tril(ones(m))  zeros(m)
          zeros(m)     zeros(m) ];
  H = sparse(H);

% ------------------------------------------------------------------
function H = hessian (x, sigma, lambda, auxdata)  
  [n m A y lambda] = deal(auxdata{:});
  H = [ tril(A'*A)  zeros(m)
         zeros(m)   zeros(m) ];
  H = sparse(sigma * H);
  