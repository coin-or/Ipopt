// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabprogram.h"
#include "matlabscalar.h"
#include "matlabstring.h"
#include "sparsematrix.h"

// Function declarations.
// ---------------------------------------------------------------
// Destroy an array of Matlab arrays.
void destroyMatlabArrays (mxArray* ptrs[], int n);

// Function definitions for class MatlabProgram.
// ---------------------------------------------------------------
MatlabProgram::MatlabProgram (const ArrayOfMatrices& x0, 
			      const ArrayOfMatrices& lb,
			      const ArrayOfMatrices& ub, 
			      const Matrix& constraintlb,
			      const Matrix& constraintub, 
			      const char* objFunc, 
			      const char* gradFunc, 
			      const char* constraintFunc, 
			      const char* jacobianFunc,
			      const char* hessianFunc,
			      const char* iterFunc, 
			      const mxArray* auxData, 
			      ArrayOfMatrices& xsol,
			      bool useQuasiNewton,
			      Multipliers* multipliers)
  : lb(lb), ub(ub), constraintlb(constraintlb), 
    constraintub(constraintub), auxData(auxData), xsol(xsol),
    multipliers(multipliers), useQuasiNewton(useQuasiNewton), 
    objFunc(objFunc), gradFunc(gradFunc), constraintFunc(constraintFunc), 
    jacobianFunc(jacobianFunc), hessianFunc(hessianFunc),
    iterFunc(iterFunc) { 
   x                 = 0;
   lambda            = 0;
   prhs              = 0;
   JacobianStructure = 0;
   HessianStructure  = 0;

  // Allocate memory for the input array to a Matlab callback
  // routine. There could be as many as n + 4 input arguments to a
  // Matlab callback routine, where n is the numer of entries of the
  // cell array.
  int n = x0.length();
  prhs  = new mxArray*[n + 4];

  // Allocate memory for the variables container.
  x = new ArrayOfMatrices(prhs,x0);
  
  // Allocate memory for the Lagrange multipliers container. We assume
  // that lambda is a row vector.
  int m  = constraintlb.length();
  lambda = new Matrix(lambdarhs,1,m);
  lambda->setvalue(0);
  
  // Copy the initial point to x.
  *x = x0;

  // Set the MATLAB trap flag so that MATLAB returns control to the
  // MEX file when an error is thrown by one of the MATLAB callback
  // routines.
  mexSetTrapFlag(1);
}

MatlabProgram::MatlabProgram (const MatlabProgram& source)  
  : lb(source.lb), ub(source.ub), constraintlb(source.constraintlb), 
    constraintub(source.constraintub), xsol(source.xsol) { }

MatlabProgram::~MatlabProgram() { 
  if (prhs)              delete[] prhs;
  if (x)                 delete x;
  if (lambda)            delete lambda;
  if (JacobianStructure) delete JacobianStructure;
  if (HessianStructure)  delete HessianStructure;
}

bool MatlabProgram::get_nlp_info (int& numVariables, int& numConstraints, 
				  int& sizeOfJ, int& sizeOfH, 
				  IndexStyleEnum& indexStyle) 
  try {

    // Get the number of variables.
    numVariables = x->numelems();
  
    // Get the number of constraints.
    numConstraints = constraintlb.length();
    
    // To get the size of the Jacobian, we call the Matlab routine
    // that computes the Jacobian. We don't need to store the actual
    // result, we just need to look at the number of non-zero
    // entries. Note that we only need to call the Matlab routine if
    // there is at least one constraint.
    if (numConstraints == 0)
      sizeOfJ = 0;
    else {
      mxArray* ptr = callMatlabJacobianRoutine(*x);
      sizeOfJ      = getSparseMatrixSize(ptr);
      mxDestroyArray(ptr);
    }
      
    // To get the size of the symmetric Hessian matrix, we call the
    // Matlab routine that computes the Hessian. We don't need to
    // store the actual result, we just need to look at the number of
    // non-zero entries in the lower triangular part of the matrix.
    if (useQuasiNewton)
      sizeOfH = 0;
    else { 
      mxArray* ptr = callMatlabHessianRoutine(*x,*lambda);
      if (!isSparseLowerTriangular(ptr))
	throw MatlabException("Matlab array must be a lower triangular \
matrix; type HELP TRIL");
      sizeOfH = getSparseMatrixSize(ptr);
      mxDestroyArray(ptr);
    } 
      
    // Use C-style indexing.
    indexStyle = C_STYLE;
    
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::get_bounds_info (int numVariables, double* lbptr, 
				     double* ubptr, int numConstraints, 
				     double* lbcptr, double* ubcptr) 
  try {

    // Fill in the structures "lb", "ub", "constraintlb" and
    // "constraintub" with the bounds information.
    lb.copyto(lbptr);
    ub.copyto(ubptr);
    constraintlb.copyto(lbcptr);
    constraintub.copyto(ubcptr);
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::get_starting_point (int numVariables, 
					bool initializeVars, 
					double* variables, 
					bool initializez, double* zl, 
					double* zu, int numConstraints, 
					bool initializeLambda, 
					double* lambda) 
  try {
    if (initializez || initializeLambda)
      throw MatlabException("Initialization of Lagrange multipliers \
requested");
    
    if (initializeVars)
      x->copyto(variables);
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_f (int numVariables, const double* variables, 
			    bool ignoreThis, double& objective) 
  try {
    x->inject(variables);
    objective = computeObjective(*x);
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_grad_f (int numVariables, const double* variables, 
				 bool ignoreThis, double* gradient) 
  try {
    x->inject(variables);
    ArrayOfMatrices grad(gradient,*x);
    computeGradient(*x,grad);
  return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_g (int numVariables, const double* variables, 
			    bool ignoreThis, int numConstraints, 
			    double* constraints) 
  try {
    if (numConstraints > 0) {
      x->inject(variables);
      Array<double> g(constraints,numConstraints);
      computeConstraints(*x,g);
    }
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_jac_g (int numVariables, const double* variables, 
				bool ignoreThis, int numConstraints, 
				int sizeOfJ, int* rows, int *cols, 
				double* Jacobian) 
  try {
    if (numConstraints > 0) {
      if (Jacobian == 0) {
	
	// Delete any previous structure information concerning the
	// Jacobian matrix.
	if (JacobianStructure) {
	  delete JacobianStructure;
	  JacobianStructure = 0;
	}
	  
	// Return the sparse matrix structure of the Jacobian.
	mxArray* ptr = callMatlabJacobianRoutine(*x);
	JacobianStructure = new SparseMatrixStructure(ptr,true);
	if ((int)JacobianStructure->size() != sizeOfJ)
	  throw MatlabException("Invalid constraint Jacobian passed back \
from Matlab routine");
	
	// Copy the sparse matrix structure to the IPOPT sparse matrix
	// format.
	JacobianStructure->getColsAndRows(cols,rows);
	  
	// Free the dynamically allocated memory.
	mxDestroyArray(ptr);
      } else {
	
	// Return the value of the Jacobian.
	x->inject(variables);
	computeJacobian(*x,Jacobian);
      }
    }
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

bool MatlabProgram::eval_h (int numVariables, const double* variables, 
			    bool ignoreThis, double sigma, 
			    int numConstraints, const double* multipliers, 
			    bool ignoreThisToo, int sizeOfH, int* rows, 
			    int* cols, double* Hessian)
  try {
    if (Hessian == 0) {

      // Delete any previous structure information concerning the
      // Jacobian matrix.
      if (HessianStructure) {
	delete HessianStructure;
	HessianStructure = 0;
      }
      
      // Return the sparse matrix structure of the symmetric Hessian.
      mxArray* ptr = callMatlabHessianRoutine(*x,*lambda);
      HessianStructure = new SparseMatrixStructure(ptr,true);
      if ((int)HessianStructure->size() != sizeOfH)
	throw MatlabException("Invalid Hessian of the Lagrangian passed \
back from Matlab routine");

      // Copy the sparse matrix structure to the IPOPT sparse matrix
      // format.
      HessianStructure->getColsAndRows(cols,rows);

      // Free the dynamically allocated memory.
      mxDestroyArray(ptr);
    } else {
	
      // Return the value of the lower triangular part of the Hessian.
      x->inject(variables);
      lambda->inject(multipliers);
      computeHessian(*x,*lambda,sigma,Hessian);
    }
    return true;
  } catch (std::exception& error) {
    mexPrintf(error.what());
    mexPrintf("\n");
    throw;
  }

void MatlabProgram::finalize_solution (SolverReturn status, int numVariables, 
				       const double* variables, 
				       const double* zl, const double* zu, 
				       int numConstraints,
				       const double* constraints, 
				       const double* lambda, 
				       double objective,
				       const IpoptData* ip_data,
				       IpoptCalculatedQuantities* ip_cq) {

  // Get the current solution.
  xsol.inject(variables);

  // If requested, store the value of the Lagrange multipliers at the
  // solution.
  if (multipliers) {
    multipliers->lowerbounds().inject(zl);
    multipliers->upperbounds().inject(zu);
    multipliers->constraints().inject(lambda);
  }  

  switch (status) {
  case Ipopt::INTERNAL_ERROR:
  case Ipopt::INVALID_NUMBER_DETECTED:
  case Ipopt::USER_REQUESTED_STOP:
  case Ipopt::RESTORATION_FAILURE:
  case Ipopt::LOCAL_INFEASIBILITY:
  case Ipopt::DIVERGING_ITERATES:
  case Ipopt::SUCCESS:
  case Ipopt::MAXITER_EXCEEDED:
  case Ipopt::STOP_AT_TINY_STEP:
  case Ipopt::STOP_AT_ACCEPTABLE_POINT:
  case Ipopt::ERROR_IN_STEP_COMPUTATION:
  case Ipopt::TOO_FEW_DEGREES_OF_FREEDOM:

    // Do nothing.
    break;
  }
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
  if (strlen(iterFunc) > 0) {
    int      nrhs = 2 + (bool) auxData;  
    mxArray* prhs[3];
    mxArray* plhs[1];     // The outputs from the Matlab routine.
    
    // The first argument is the current iteration, t. The argument
    // after that is the current value of the objective function, f.
    MatlabScalar t(prhs[0],iteration);
    MatlabScalar f(prhs[1],objective);
    
    // Pass the auxiliary data, if necessary.
    if (auxData)
      prhs[2] = const_cast<mxArray*>(auxData);
    
    // Call the designated Matlab iterative callback routine. It
    // takes as input the values of the variables, the current
    // iteration, and the current value of the objective function.
    int exitstatus = mexCallMATLAB(0,plhs,nrhs,prhs,iterFunc);

    // Deallocate the dynamically allocated memory.
    mxDestroyArray(prhs[0]);
    mxDestroyArray(prhs[1]);

    if (exitstatus)
      throw MatlabException("Call to MATLAB iterative callback routine");
  }    
  
  return true;
}

double MatlabProgram::computeObjective (const ArrayOfMatrices& x) {
  int      nrhs = x.length() + (bool) auxData;  
  int      nlhs = 1;    // The number of outputs from Matlab.
  mxArray* plhs[1];     // The outputs from the Matlab routine.

  // Pass the auxiliary data, if necessary.
  if (auxData)
    prhs[x.length()] = const_cast<mxArray*>(auxData);

  // Call the designated Matlab routine for evaluating the objective
  // function. It takes as input the values of the variables, and
  // returns a single output, the value of the objective function.
  int exitstatus = mexCallMATLAB(nlhs,plhs,nrhs,prhs,objFunc);
  if (exitstatus)
    throw MatlabException("Evaluation of objective function failed in \
call to MATLAB routine");

  // Get the result passed back from the Matlab routine.
  MatlabScalar matlabOutput(plhs[0]);
  double       f = matlabOutput;

  // Free the dynamically allocated memory.
  mxDestroyArray(plhs[0]);

  return f;
}

void MatlabProgram::computeGradient (const ArrayOfMatrices& x,
				     ArrayOfMatrices& grad) {
  int       nrhs = x.length() + (bool) auxData;  
  int       nlhs = x.length();          // The number of outputs from Matlab.
  mxArray** plhs = new mxArray*[nlhs];  // The outputs to the Matlab routine.

  // Pass the auxiliary data, if necessary.
  if (auxData)
    prhs[x.length()] = const_cast<mxArray*>(auxData);

  // Call the designated Matlab routine for computing the gradient
  // of the objective function. It takes as input the values of the
  // variables, and returns as many outputs corresponding to the
  // partial derivatives of the objective function with respect to
  // the variables at their curent values.
  int exitstatus = mexCallMATLAB(nlhs,plhs,nrhs,prhs,gradFunc);
  if (exitstatus) {
    delete [] plhs;
    throw MatlabException("Evaluation of objective gradient failed in \
call to MATLAB routine");
  }

  // Get the result passed back from the Matlab routine.
  ArrayOfMatrices matlabOutput((const mxArray**) plhs,nlhs);
  if (matlabOutput.length() != grad.length())
    throw MatlabException("Invalid gradient passed back from MATLAB \
routine");
  for (int i = 0; i < grad.length(); i++)
    if (matlabOutput[i]->length() != grad[i]->length()) {
      delete [] plhs;
      throw MatlabException("Invalid gradient passed back from MATLAB \
routine");
    }
  grad = matlabOutput;

  // Free the dynamically allocated memory.
  destroyMatlabArrays(plhs,nlhs);
  delete [] plhs;
}

void MatlabProgram::computeConstraints (const ArrayOfMatrices& x,
					Array<double>& g) {
  int      nrhs = x.length() + (bool) auxData;  
  int      nlhs = 1;  // The number of outputs from Matlab.
  mxArray* plhs[1];   // The outputs to the Matlab routine.

  // Pass the auxiliary data, if necessary.
  if (auxData)
    prhs[x.length()] = const_cast<mxArray*>(auxData);

  // Call the designated Matlab routine for evaluating the
  // constraints at the current point. It takes as input the values
  // of the variables, and returns the set of constraints evaluated
  // at the current point.
  if (mexCallMATLAB(nlhs,plhs,nrhs,prhs,constraintFunc))
    throw MatlabException("Evaluation of constraint equations failed in \
call to MATLAB routine");

  // Get the result passed back from the Matlab routine.
  Array<double>* matlabOutput = new Matrix(plhs[0]);
  if (matlabOutput->length() != g.length())
    throw MatlabException("Invalid constraint values passed back from \
Matlab routine");
  g.inject(*matlabOutput);
    
  // Free the dynamically allocated memory.
  mxDestroyArray(plhs[0]);
}

void MatlabProgram::computeJacobian (const ArrayOfMatrices& x, 
				     double* Jacobian) {

  // Get the result passed back from the Matlab routine.
  mxArray* ptr = callMatlabJacobianRoutine(x,false);

  // Copy the Matlab output into the IPOPT result.
  SparseMatrixStructure newJacobianStructure(ptr);
  copyElems(newJacobianStructure,*JacobianStructure,mxGetPr(ptr),Jacobian);

  // Get rid of the dynamically allocated memory.
  mxDestroyArray(ptr);
}

void MatlabProgram::computeHessian (const ArrayOfMatrices& x, 
				    const Array<double>& lambda, 
				    double sigma, double* Hessian) {

  // Get the result passed back from the Matlab routine.
  mxArray* ptr = callMatlabHessianRoutine(x,lambda,false,sigma);

  // Copy the Matlab output into the IPOPT result.
  SparseMatrixStructure newHessianStructure(ptr);
  copyElems(newHessianStructure,*HessianStructure,mxGetPr(ptr),Hessian);

  // Free the dynamically allocated memory.
  mxDestroyArray(ptr);
}

mxArray* MatlabProgram::callMatlabJacobianRoutine (const ArrayOfMatrices& x, 
						   bool returnStructureOnly) {
  int      nrhs = x.length() + 1 + (bool) auxData;
  int      nlhs = 1;  // The number of outputs from Matlab.
  mxArray* plhs[1];   // The outputs to the Matlab routine.

  // Pass the "returnStructureOnly" boolean variable.
  MatlabScalar RSOMatlabInput(prhs[x.length()],returnStructureOnly);

  // Pass the auxiliary data, if necessary.
  if (auxData)
    prhs[x.length() + 1] = const_cast<mxArray*>(auxData);

  // Call the designated Matlab routine for evaluating the
  // Jacobian of the constraints at the current point. It takes as
  // input the values of the variables, and returns a single
  // matrix containing the value of the Jacobian evaluated at the
  // current point.
  if (mexCallMATLAB(nlhs,plhs,nrhs,prhs,jacobianFunc))
    throw MatlabException("Evaluation of constraint Jacobian failed in \
call to MATLAB routine");

  // Check whether the output from the Matlab routine makes sense.
  int      numConstraints = constraintlb.length();
  int      numVariables   = x.numelems();
  mxArray* ptr            = plhs[0];
  if ((int) mxGetM(ptr) != numConstraints || 
      (int) mxGetN(ptr) != numVariables)
    throw MatlabException("Invalid constraint Jacobian passed back from \
Matlab routine");

  // Clear the dynamically allocated memory, except for the
  // information that is returned from this function.
  mxDestroyArray(prhs[x.length()]);

  return plhs[0];
}

mxArray* MatlabProgram::callMatlabHessianRoutine (const ArrayOfMatrices& x, 
						  const Array<double>& lambda, 
						  bool returnStructureOnly, 
						  double sigma) {
  int      nrhs = x.length() + 3 + (bool) auxData;  
  int      nlhs = 1;  // The number of outputs from Matlab.
  mxArray* plhs[1];   // The outputs to the Matlab routine.

  // The argument after the variables is "sigma".
  MatlabScalar sigmaMatlabInput(prhs[x.length()],sigma);

  // The next argument is "lambda".
  prhs[x.length() + 1] = lambdarhs;

  // The next argument is "returnStructureOnly".
  MatlabScalar RSOMatlabInput(prhs[x.length() + 2],returnStructureOnly);

  // Pass the auxiliary data, if necessary.
  if (auxData)
    prhs[3 + x.length()] = const_cast<mxArray*>(auxData);

  // Call the designated Matlab routine for evaluating the Hessian
  // of the Lagrangian function at the current point. It takes as
  // input the factor "sigma", the current values of the Lagrange
  // multipliers "lambda" and the current values of the variables,
  // and returns a single matrix containing the value of the
  // Hessian evaluated at the current point.
  if (mexCallMATLAB(nlhs,plhs,nrhs,prhs,hessianFunc))
    throw MatlabException("Evaluation of Hessian of the Lagrangian \
failed in call to MATLAB routine");

  // Check whether the output from the Matlab routine makes sense.
  int numVariables = x.numelems();
  mxArray* ptr     = plhs[0];
  if ((int) mxGetM(ptr) != numVariables || (int) mxGetN(ptr) != numVariables)
    throw MatlabException("Invalid Hessian of the Lagrangian passed \
back from Matlab routine");
  if (!isSparseLowerTriangular(ptr))
    throw MatlabException("Matlab array must be a lower triangular \
matrix; type HELP TRIL");

  // Free some of the dynamically allocated memory.
  mxDestroyArray(prhs[x.length()]);
  mxDestroyArray(prhs[x.length()+2]);

  return plhs[0];
}

// Function definitions.
// ---------------------------------------------------------------
void destroyMatlabArrays (mxArray* ptrs[], int n) {
  for (int i = 0; i < n; i++)
    mxDestroyArray(ptrs[i]);
}
