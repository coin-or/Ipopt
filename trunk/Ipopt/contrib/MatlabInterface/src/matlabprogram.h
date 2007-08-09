// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_MATLABPROGRAM
#define INCLUDE_MATLABPROGRAM

#include "array.h"
#include "matlabexception.h"
#include "matlabmatrix.h"
#include "sparsematrix.h"
#include "arrayofmatrices.h"
#include "ipopt/IpTNLP.hpp"
#include "mex.h"

using Ipopt::TNLP;
using Ipopt::SolverReturn;
using Ipopt::AlgorithmMode;
using Ipopt::IpoptData;
using Ipopt::IpoptCalculatedQuantities;

// Class MatlabProgram
// -----------------------------------------------------------------
class MatlabProgram : public TNLP {
public:
    
  // The constructor. The input "error" will be set to the appropriate
  // exception object if the IPOPT solver terminates abnormally.
  MatlabProgram (const ArrayOfMatrices& x0, const ArrayOfMatrices& lb,
		 const ArrayOfMatrices& ub, const Matrix& constraintlb,
		 const Matrix& constraintub, const char* objFunc, 
		 const char* gradFunc, const char* constraintFunc, 
		 const char* jacobianFunc, const char* hessianFunc,
		 const char* iterFunc, const mxArray* auxData, 
		 ArrayOfMatrices& xsol, bool useQuasiNewton);
    
  // The destructor.
  virtual ~MatlabProgram();

  // Returns the error generated. If no error was generated, returns a
  // null pointer.
  char* geterrormsg() const;
  
  // Method to return some info about the nonlinear program.
  virtual bool get_nlp_info (int& numVariables, int& numConstraints, 
			     int& sizeOfJ, int& sizeOfH, 
			     IndexStyleEnum& indexStyle);
  
  // Return the bounds for the problem.
  virtual bool get_bounds_info (int numVariables, double* lbptr, 
				double* ubptr, int numConstraints, 
				double* clbptr, double* cubptr);
    
  // Return the starting point for the algorithm.
  virtual bool get_starting_point (int numVariables, bool initializeVars, 
				   double* variables, 
				   bool initializez, double* zl, 
				   double* zu, int numConstraints, 
				   bool initializeLambda, double* lambda);
    
  // Compute the value of the objective.
  virtual bool eval_f (int numVariables, const double* variables, 
		       bool ignoreThis, double& objective);
    
  // Compute the gradient of the objective.
  virtual bool eval_grad_f (int numVariables, const double* variables, 
			    bool ignoreThis, double* gradient);
    
  // Evaluate the constraint residuals.
  virtual bool eval_g (int numVariables, const double* variables, 
		       bool ignoreThis, int numConstraints, 
		       double* constraints);
    
  // This method either returns: 1.) The structure of the Jacobian
  // (if "Jacobian" is zero), or 2.) The values of the Jacobian (if
  // "Jacobian" is not zero).
  virtual bool eval_jac_g (int numVariables, const double* variables, 
			   bool ignoreThis, int numConstraints, 
			   int sizeOfJ, int* rows, int *cols, 
			   double* Jacobian);
    
  // This method either returns: 1.) the structure of the Hessian of
  // the Lagrangian (if "Hessian" is zero), or 2.) the values of the
  // Hessian of the Lagrangian (if "Hesson" is not zero).
  virtual bool eval_h (int numVariables, const double* variables, 
		       bool ignoreThis, double sigma, 
		       int numConstraints, const double* multipliers, 
		       bool ignoreThisToo, int sizeOfH, int* rows, 
		       int* cols, double* Hessian);

  // This method is called when the algorithm is complete.
  virtual void finalize_solution (SolverReturn status, int numVariables, 
				  const double* variables, const double* zl, 
				  const double* zu, int numConstraints, 
				  const double* constraints, 
				  const double* lambda, double objective,
                                  const IpoptData* ip_data,
                                  IpoptCalculatedQuantities* ip_cq);

  // Intermediate callback method. It is called once per iteration
  // of the IPOPT algorithm.
  virtual bool intermediate_callback (AlgorithmMode mode,
				      int iteration, double objective,
				      double inf_pr, double inf_du,
				      double mu, double d_norm,
				      double regularization_ize,
				      double alpha_du, double alpha_pr,
				      int ls_trials,
				      const IpoptData* ip_data,
				      IpoptCalculatedQuantities* ip_cq);

protected:
  const ArrayOfMatrices& lb;           // Lower bounds on the variables.
  const ArrayOfMatrices& ub;           // Upper bounds on the variables.
  const Matrix&          constraintlb; // Lower bounds on the constraints.
  const Matrix&          constraintub; // Upper bounds on the constraints.
  const mxArray*         auxData;      // Auxiliary data passed to the 
                                       // Matlab callback routines.
  ArrayOfMatrices&       xsol;         // Storage for the solution.
  ArrayOfMatrices*       x;            // Current value of the variables 
                                       // that's passed to the Matlab 
                                       // callback routines.
  Array<double>*         lambda;       // Current value of the Lagrange 
                                       // multipliers that's passed to the
                                       // Matlab callback routine for
				       // computing the Hessian.
  
  // Array of inputs to a Matlab callback function.
  mxArray** prhs;
  mxArray*  lambdarhs;

  // If true, we don't need to compute the Hessian of the Lagrangian.
  bool useQuasiNewton; 

  // These next two members store information about the structure of
  // the sparse Matlab matrix for the Jacobian of the constraints
  // and the Hessian of the Lagragian.
  SparseMatrixStructure* JacobianStructure;
  SparseMatrixStructure* HessianStructure;
  
  // The following members specify the Matlab callback routines.
  const char* objFunc;        // Computes the objective function.
  const char* gradFunc;       // Computes the gradient of the objective.
  const char* constraintFunc; // Computes the constraints.
  const char* jacobianFunc;   // Computes the Jacobian of the constraints.
  const char* hessianFunc;    // Computes the Hessian of the Lagrangian.
  const char* iterFunc;       // The iterative callback function.

  // The copy constructor and copy assignment operator are kept
  // private so that they are not used.
  MatlabProgram            (const MatlabProgram& source);
  MatlabProgram& operator= (const MatlabProgram& source) { return *this; };

private:
  
  // Return the value of the objective function at the point "x".
  double computeObjective (const ArrayOfMatrices& x);

  // Compute the first derivatives of the objective function at the
  // point "x" and return their values in "grad".
  void computeGradient (const ArrayOfMatrices& x, ArrayOfMatrices& grad);
  
  // Compute the value of each constraint function at the point "x"
  // and return the result in second input argument, which is an
  // array of length equal to the number of constraints.
  void computeConstraints (const ArrayOfMatrices& x, Array<double>& g);

  // Compute the derivatives of the Jacobian of the constraints at
  // the point "x". For this function to work correctly, the Matlab
  // routine must return the expected sparse structure of the
  // Jacobian matrix.
  void computeJacobian (const ArrayOfMatrices& x, double* Jacobian);
  
  // Compute the derivatives of the Hessian of the Lagrangian at the
  // point "x". 
  void computeHessian (const ArrayOfMatrices& x, 
		       const Array<double>& lambda,
		       double sigma, double* Hessian);

  // Call the Matlab routine for computing the Jacobian. The input
  // arguments to the Matlab routine are the current values of the
  // optimization variables. The result must be a sparse
  // matrix. Additionally, there is a second input argument passed
  // to the Matlab routine. If it is true (1), then we ignore the
  // actual values entered into the sparse matrix, as long as every
  // entry that might possibly be a value other than zero at any
  // time be returned as a non-zero value. It is up to the
  // programmer to properly deallocate the return value at a later
  // time by calling mxDestroyArray().
  mxArray* callMatlabJacobianRoutine (const ArrayOfMatrices& x,
				      bool returnStructureOnly = true);

  // Call the Matlab routine for computing the Hessian. The input
  // arguments to the Matlab routine are the current values of the
  // optimization variables. The result must be a sparse, lower
  // triangular matrix. Additionally, there is a second input
  // argument passed to the Matlab routine. If it is true (1), then
  // we ignore the actual values entered into the sparse matrix, as
  // long as every entry that might possibly be a value other than
  // zero at any time be returned as a non-zero value. It is up to
  // the programmer to properly deallocate the return value at a
  // later time by calling mxDestroyArray().
  mxArray* callMatlabHessianRoutine  (const ArrayOfMatrices& x,
				      const Array<double>& lambda, 
				      bool returnStructureOnly = true,
				      double sigma = 0);
};

#endif
