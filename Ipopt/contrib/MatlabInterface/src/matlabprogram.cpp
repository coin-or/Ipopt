// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 25, 2008

#include "matlabprogram.hpp"
#include "sparsematrix.hpp"
#include "matlabexception.hpp"

// Function definitions for class MatlabProgram.
// ---------------------------------------------------------------
MatlabProgram::MatlabProgram (const Iterate& x0, 
			      const CallbackFunctions& funcs,
			      const Options& options, Iterate& x, 
			      MatlabInfo& info)
  : x0(x0), funcs(funcs), options(options), x(x), 
    info(info), J(0), H(0) { }

MatlabProgram::~MatlabProgram() { 
  if (J) delete J;
  if (H) delete H;
}

bool MatlabProgram::get_nlp_info (int& n, int& m, int& sizeOfJ, int& sizeOfH, 
				  IndexStyleEnum& indexStyle) 
  try {

    // Get the number of variables and constraints.
    n = numvars(options);
    m = numconstraints(options);

    // Get the size of the Jacobian.
    if (m) {
      if (!funcs.jacobianFuncIsAvailable())
	throw MatlabException("You need to specify the callback functions \
for computing the Jacobian and the sparsity structure of the Jacobian");
      SparseMatrix* J = funcs.getJacobianStructure(n,m);
      sizeOfJ = J->numelems();
      delete J;
    }
    else
      sizeOfJ = 0;

    // Get the size of the symmetric Hessian matrix. We don't need to
    // store the actual result, we just need to look at the number of
    // non-zero entries in the lower triangular part of the matrix.
    if (!options.ipoptOptions().useQuasiNewton()) {
      if (!funcs.hessianFuncIsAvailable())
	throw MatlabException("You need to specify the callback functions \
for computing the Hessian and the sparsity structure of the Hessian");
      SparseMatrix* H = funcs.getHessianStructure(n);
      sizeOfH = H->numelems();
      delete H;
    }
    else
      sizeOfH = 0;

    // Use C-style indexing.
    indexStyle = C_STYLE;

    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::get_bounds_info (int n, double* lb, double* ub, 
				     int m, double* cl, double* cu) 
  try {

    // Fill in the structures with the bounds information.
    copymemory(options.lowerbounds(),lb,n);
    copymemory(options.upperbounds(),ub,n);
    copymemory(options.constraintlb(),cl,m);
    copymemory(options.constraintub(),cu,m);
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::get_starting_point (int n, bool initializeVars, 
					double* vars, bool initializez, 
					double* zl, double* zu, int m, 
					bool initializeLambda, double* lambda) 
  try {

    // Check to see whether IPOPT is requesting the initial point for
    // the primal variables.
    if (initializeVars)
      x0.copyto(vars);
    
    // Check to see whether IPOPT is requesting the initial point for
    // the Lagrange multipliers associated with the bounds on the
    // optimization variables.
    if (initializez) {
      if (!options.multlb() || !options.multub())
	throw MatlabException("Initialization of Lagrange multipliers \
for lower and upper bounds requested but initial values are not supplied");
      copymemory(options.multlb(),zl,n);
      copymemory(options.multub(),zu,n);
    }
    
    // Check to see whether IPOPT is requesting the initial point for
    // the Lagrange multipliers corresponding to the equality and
    // inequality constraints.
    if (initializeLambda && m>0) {
      if (!options.multconstr())
	throw MatlabException("Initialization of Lagrange multipliers \
for constraints are requested but initial values are not provided");
      copymemory(options.multconstr(),lambda,m);
    }
    
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_f (int n, const double* vars, bool ignore, double& f) 
  try {
    x.inject(vars);
    f = funcs.computeObjective(x);
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_grad_f (int n, const double* vars, bool ignore, 
				 double* grad) 
  try {
    x.inject(vars);
    funcs.computeGradient(x,grad);
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_g (int n, const double* vars, bool ignore, int m, 
			    double* g) 
  try {
    if (m > 0) {
      if (!funcs.constraintFuncIsAvailable())
	throw MatlabException("You need to specify the callback function \
for computing the constraints");
      x.inject(vars);
      funcs.computeConstraints(x,m,g);
    }
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_jac_g (int n, const double* vars, bool ignore, int m, 
				int sizeOfJ, int* rows, int *cols, double* Jx) 
  try {
    if (m > 0) {
      if (!funcs.jacobianFuncIsAvailable())
	throw MatlabException("You need to specify the callback functions \
for computing the Jacobian and the sparsity structure of the Jacobian");

      // If the input Jx is 0, then return the sparsity structure of
      // the Jacobian. Otherwise, return the values of the nonzero
      // entries.
      if (Jx == 0) {
	
 	// Delete any previous structure information concerning the
 	// Jacobian matrix.
 	if (J) {
 	  delete J;
 	  J = 0;
 	}
	  
 	// Get the sparse matrix structure of the Jacobian.
	J = funcs.getJacobianStructure(n,m);
 	if (J->numelems() != sizeOfJ)
 	  throw MatlabException("The constraint Jacobian passed back from \
the MATLAB routine has an incorrect number of nonzero entries");
	
 	// Copy the sparse matrix structure to the IPOPT sparse matrix
 	// format.
 	J->getColsAndRows(cols,rows);
       } else {
	
 	// Return the value of the Jacobian.
 	x.inject(vars);
	funcs.computeJacobian(m,x,*J);
	J->copyto(Jx);
      }
    }
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_h (int n, const double* vars, bool ignore, 
			    double sigma, int m, const double* lambda, 
			    bool ignoretoo, int sizeOfH, int* rows, 
			    int* cols, double* Hx)
  try {

    // If the input Hx is 0, then return the sparsity structure of the
    // Hessian. Otherwise, return the values of the nonzero entries.
    if (Hx == 0) {
      if (!funcs.hessianFuncIsAvailable())
	throw MatlabException("You need to specify the callback functions \
for computing the Hessian and the sparsity structure of the Hessian");

      // Delete any previous structure information concerning the
      // Hessian matrix.
      if (H) {
	delete H;
	H = 0;
      }
      
      // Return the sparse matrix structure of the symmetric Hessian.
      H = funcs.getHessianStructure(n);
      if (H->numelems() != sizeOfH)
	throw MatlabException("The Hessian passed back from the MATLAB \
routine has an incorrect number of nonzero entries");

      // Copy the sparse matrix structure to the IPOPT sparse matrix
      // format.
      H->getColsAndRows(cols,rows);
    } else {
	
      // Return the value of the lower triangular portion of the Hessian.
      x.inject(vars);
      funcs.computeHessian(x,sigma,m,lambda,*H);
      H->copyto(Hx);
    }
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

void MatlabProgram::finalize_solution (SolverReturn status, int n, 
				       const double* vars, const double* zl, 
				       const double* zu, int m,
				       const double* g, const double* lambda, 
				       double f, const IpoptData* ip_data,
				       IpoptCalculatedQuantities* ip_cq) {

  // Get the current solution.
  x.inject(vars);

  // Store the value of the Lagrange multipliers at the solution.
  info.setmultlb(n,zl);
  info.setmultub(n,zu);
  info.setmultconstr(m,lambda);
}
  
bool MatlabProgram::intermediate_callback (AlgorithmMode mode,
					   int t, double f, double inf_pr, 
					   double inf_du, double mu, 
					   double d_norm,
					   double regularization_size,
					   double alpha_du, double alpha_pr,
					   int ls_trials,
					   const IpoptData* ip_data,
					   IpoptCalculatedQuantities* ip_cq)
  try {  
    if (funcs.iterFuncIsAvailable())
      return funcs.iterCallback(t,f,inf_pr,inf_du,mu,d_norm,regularization_size,alpha_du,alpha_pr,ls_trials,ip_data,ip_cq,numvars(x));
    else
      return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }
