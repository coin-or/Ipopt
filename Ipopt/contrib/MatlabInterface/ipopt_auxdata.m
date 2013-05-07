function [x, info] = ipopt_auxdata(x0, funcs, options)
% Wrapper function to implement auxdata functionality using Matlab function handles

if ~isfield(options, 'auxdata')
    % no auxdata given, call ipopt as normal
    [x, info] = ipopt(x0, funcs, options);
else
    % remove auxdata from options structure and modify function handles
    auxdata = options.auxdata;
    options_new = rmfield(options, 'auxdata');
    funcs_new.objective = @(x) funcs.objective(x, auxdata);
    funcs_new.gradient  = @(x) funcs.gradient(x, auxdata);
    if isfield(funcs, 'constraints')
        funcs_new.constraints = @(x) funcs.constraints(x, auxdata);
    end
    if isfield(funcs, 'jacobian')
        funcs_new.jacobian = @(x) funcs.jacobian(x, auxdata);
    end
    if isfield(funcs, 'jacobianstructure')
        funcs_new.jacobianstructure = @() funcs.jacobianstructure(auxdata);
    end
    if isfield(funcs, 'hessian')
        funcs_new.hessian = @(x, sigma, lambda) funcs.hessian(x, sigma, lambda, auxdata);
    end
    if isfield(funcs, 'hessianstructure')
        funcs_new.hessianstructure = @() funcs.hessianstructure(auxdata);
    end
    if isfield(funcs, 'iterfunc')
        % This assumes the pre-3.11 convention for iterfunc. If you would
        % like to use the additional information that is available via the
        % third input argument to iterfunc as of Ipopt version 3.11, you
        % will need to modify this section by uncommenting the line below.
        funcs_new.iterfunc = @(t, f) funcs.iterfunc(t, f, auxdata);
        % funcs_new.iterfunc = @(t, f, varstruct) funcs.iterfunc(t, f, varstruct, auxdata);
    end
    [x, info] = ipopt(x0, funcs_new, options_new);
end
