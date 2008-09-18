// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabprogram.h"
#include "sparsematrix.h"
#include "matlabexception.h"

// Function definitions for class MatlabProgram.
// ---------------------------------------------------------------
MatlabProgram::MatlabProgram (const Iterate& x0, 
			      const CallbackFunctions& funcs,
			      const Options& options, Iterate& x, 
			      const mxArray* auxdata, MatlabInfo& info)
  : x0(x0), funcs(funcs), options(options), x(x), auxdata(auxdata), 
    info(info), J(0), H(0) {

  // Allocate memory for the input array to a Matlab callback
  // routine. There could be as many as n + 4 input arguments to a
  // Matlab callback routine, where n is the numer of entries of the
  // cell array.
//   int n   = x0.length();
//   prhs    = new mxArray*[n + 5];
//   prhs[0] = 0;

//   // Allocate memory for the variables container.
//   x = new ArrayOfMatrices(&prhs[1],x0);
  
//   // Allocate memory for the Lagrange multipliers container. We assume
//   // that lambda is a row vector.
//   int m  = constraintlb.length();
//   lambda = new Matrix(lambdarhs,1,m);
//   lambda->setvalue(0);
  
//   // Copy the initial point to x.
//   *x = x0;

//   // Set the MATLAB trap flag so that MATLAB returns control to the
//   // MEX file when an error is thrown by one of the MATLAB callback
//   // routines.
//   mexSetTrapFlag(1);
}

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
    SparseMatrix* J = funcs.getJacobianStructure(n,m,auxdata);
    sizeOfJ         = J->numelems();

    // Get the size of the symmetric Hessian matrix. We don't need to
    // store the actual result, we just need to look at the number of
    // non-zero entries in the lower triangular part of the matrix.
    SparseMatrix* H = funcs.getHessianStructure(n,auxdata);
    sizeOfH         = H->numelems();

    // Use C-style indexing.
    indexStyle = C_STYLE;

    // Free the dynamically allocated memory.
    delete J;
    delete H;

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

    // Initialize the number of iterations.
//     numiter = 0;

//     // Check to see whether IPOPT is requesting the initial point for
//     // the Lagrange multipliers associated with the bounds on the
//     // optimization variables.
//     if (initializez) {
//       if (initialMultipliers == 0)
// 	throw MatlabException("Initialization of Lagrange multipliers \
// requested but initial values are not provided");
//       initialMultipliers->lowerbounds().copyto(zl);
//       initialMultipliers->upperbounds().copyto(zu);
//     }

//     // Check to see whether IPOPT is requesting the initial point for
//     // the Lagrange multipliers corresponding to the equality and
//     // inequality constraints.
//     if (initializeLambda) {
      
//       if (initialMultipliers == 0)
// 	throw MatlabException("Initialization of Lagrange multipliers \
// requested but initial values are not provided");
//       initialMultipliers->constraints().copyto(lambda);
//     }

    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_f (int n, const double* vars, bool ignore, double& f) 
  try {
    x.inject(vars);
    f = funcs.computeObjective(x,auxdata);
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
    funcs.computeGradient(x,grad,auxdata);
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
      x.inject(vars);
      funcs.computeConstraints(x,m,g,auxdata);
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
	J = funcs.getJacobianStructure(n,m,auxdata);
 	if (J->numelems() != sizeOfJ)
 	  throw MatlabException("The constraint Jacobian passed back from \
the MATLAB routine has an incorrect number of nonzero entries");
	
 	// Copy the sparse matrix structure to the IPOPT sparse matrix
 	// format.
 	J->getColsAndRows(cols,rows);
       } else {
	
 	// Return the value of the Jacobian.
 	x.inject(vars);
	funcs.computeJacobian(m,x,*J,auxdata);
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

      // Delete any previous structure information concerning the
      // Hessian matrix.
      if (H) {
	delete H;
	H = 0;
      }
      
      // Return the sparse matrix structure of the symmetric Hessian.
      H = funcs.getHessianStructure(n,auxdata);
      if (H->numelems() != sizeOfH)
	throw MatlabException("The Hessian passed back from the MATLAB \
routine has an incorrect number of nonzero entries");

      // Copy the sparse matrix structure to the IPOPT sparse matrix
      // format.
      H->getColsAndRows(cols,rows);
    } else {
	
      // Return the value of the lower triangular portion of the Hessian.
      x.inject(vars);
      funcs.computeHessian(x,sigma,m,lambda,*H,auxdata);
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

  // If requested, store the value of the Lagrange multipliers at the
  // solution.
//   if (multipliers) {
//     multipliers->lowerbounds().inject(zl);
//     multipliers->upperbounds().inject(zu);
//     multipliers->constraints().inject(lambda);
//   }  
}
  
bool MatlabProgram::intermediate_callback (AlgorithmMode mode,
					   int iteration, double objective,
					   double inf_pr, double inf_du,
					   double mu, double d_norm,
					   double regularization_ize,
					   double alpha_du, double alpha_pr,
					   int ls_trials,
					   const IpoptData* ip_data,
					   IpoptCalculatedQuantities* ip_cq) {
  
//   // Update the number of iterations of IPOPT.
//   numiter++;

//   if (iterFunc) {
//     int      nrhs = 3 + (bool) auxData;  
//     mxArray* prhs[4];  // The inputs to the Matlab routine.
//     mxArray* plhs[1];  // The outputs from the Matlab routine.
    
//     // The first argument is the function to call. The second argument
//     // is the current iteration, t. The argument after that is the
//     // current value of the objective function, f.
//     prhs[0] = iterFunc;
//     MatlabScalar t(prhs[1],iteration);
//     MatlabScalar f(prhs[2],objective);
    
//     // Pass the auxiliary data, if necessary.
//     if (auxData)
//       prhs[3] = const_cast<mxArray*>(auxData);
    
//     // Call the designated Matlab iterative callback routine. It
//     // takes as input the values of the variables, the current
//     // iteration, and the current value of the objective function.
//     int exitstatus = mexCallMATLAB(0,plhs,nrhs,prhs,"feval");

//     // Deallocate the dynamically allocated memory.
//     mxDestroyArray(prhs[1]);
//     mxDestroyArray(prhs[2]);

//     if (exitstatus)
//       throw MatlabException("Call to MATLAB iterative callback routine");
//   }    
  
  return true;
}

// void MatlabProgram::computeJacobian (const ArrayOfMatrices& x, 
// 				     double* Jacobian) {

//   // Get the result passed back from the Matlab routine.
//   mxArray* ptr = callMatlabJacobianRoutine(x,false);

//   // Copy the Matlab output into the IPOPT result.
//   SparseMatrixStructure newJacobianStructure(ptr);
//   copyElems(newJacobianStructure,*JacobianStructure,mxGetPr(ptr),Jacobian);

//   // Get rid of the dynamically allocated memory.
//   mxDestroyArray(ptr);
// }

// void MatlabProgram::computeHessian (const ArrayOfMatrices& x, 
// 				    const Array<double>& lambda, 
// 				    double sigma, double* Hessian) {

//   // Get the result passed back from the Matlab routine.
//   mxArray* ptr = callMatlabHessianRoutine(x,lambda,false,sigma);

//   // Copy the Matlab output into the IPOPT result.
//   SparseMatrixStructure newHessianStructure(ptr);
//   copyElems(newHessianStructure,*HessianStructure,mxGetPr(ptr),Hessian);

//   // Free the dynamically allocated memory.
//   mxDestroyArray(ptr);
// }

// mxArray* MatlabProgram::callMatlabJacobianRoutine (const ArrayOfMatrices& x, 
// 						   bool returnStructureOnly) {
//   int      nrhs = x.length() + 2 + (bool) auxData;
//   int      nlhs = 1;  // The number of outputs from Matlab.
//   mxArray* plhs[1];   // The outputs to the Matlab routine.

//   // The first input argument is the function to call.
//   prhs[0] = jacobianFunc;

//   // Pass the input, "returnStructureOnly" boolean variable.
//   MatlabScalar RSOMatlabInput(prhs[1 + x.length()],returnStructureOnly);

//   // Pass the auxiliary data, if necessary.
//   if (auxData)
//     prhs[2 + x.length()] = const_cast<mxArray*>(auxData);

//   // Call the designated Matlab routine for evaluating the
//   // Jacobian of the constraints at the current point. It takes as
//   // input the values of the variables, and returns a single
//   // matrix containing the value of the Jacobian evaluated at the
//   // current point.
//   if (mexCallMATLAB(nlhs,plhs,nrhs,prhs,"feval"))
//     throw MatlabException("Evaluation of constraint Jacobian failed in \
// call to MATLAB routine");

//   // Check whether the output from the Matlab routine makes sense.
//   int      numConstraints = constraintlb.length();
//   int      numVariables   = x.numelems();
//   mxArray* ptr            = plhs[0];
//   if ((int) mxGetM(ptr) != numConstraints || 
//       (int) mxGetN(ptr) != numVariables)
//     throw MatlabException("Invalid constraint Jacobian passed back from \
// Matlab routine");

//   // Clear the dynamically allocated memory, except for the
//   // information that is returned from this function.
//   mxDestroyArray(prhs[1 + x.length()]);

//   return plhs[0];
// }

// mxArray* MatlabProgram::callMatlabHessianRoutine (const ArrayOfMatrices& x, 
// 						  const Array<double>& lambda, 
// 						  bool returnStructureOnly, 
// 						  double sigma) {
//   int      nrhs = x.length() + 4 + (bool) auxData;  
//   int      nlhs = 1;  // The number of outputs from Matlab.
//   mxArray* plhs[1];   // The outputs to the Matlab routine.

//   // The first argument is the function to call.
//   prhs[0] = hessianFunc;

//   // The argument after the variables is "sigma".
//   MatlabScalar sigmaMatlabInput(prhs[1 + x.length()],sigma);

//   // The next argument is "lambda".
//   prhs[2 + x.length()] = lambdarhs;

//   // The next argument is "returnStructureOnly".
//   MatlabScalar RSOMatlabInput(prhs[3 + x.length()],returnStructureOnly);

//   // Pass the auxiliary data, if necessary.
//   if (auxData)
//     prhs[4 + x.length()] = const_cast<mxArray*>(auxData);

//   // Call the designated Matlab routine for evaluating the Hessian
//   // of the Lagrangian function at the current point. It takes as
//   // input the factor "sigma", the current values of the Lagrange
//   // multipliers "lambda" and the current values of the variables,
//   // and returns a single matrix containing the value of the
//   // Hessian evaluated at the current point.
//   if (mexCallMATLAB(nlhs,plhs,nrhs,prhs,"feval"))
//     throw MatlabException("Evaluation of Hessian of the Lagrangian \
// failed in call to MATLAB routine");

//   // Check whether the output from the Matlab routine makes sense.
//   int numVariables = x.numelems();
//   mxArray* ptr     = plhs[0];
//   if ((int) mxGetM(ptr) != numVariables || (int) mxGetN(ptr) != numVariables)
//     throw MatlabException("Invalid Hessian of the Lagrangian passed \
// back from Matlab routine");
//   if (!isSparseLowerTriangular(ptr))
//     throw MatlabException("Matlab array must be a lower triangular \
// matrix; type HELP TRIL");

//   // Free some of the dynamically allocated memory.
//   mxDestroyArray(prhs[1+x.length()]);
//   mxDestroyArray(prhs[3+x.length()]);

//   return plhs[0];
// }
