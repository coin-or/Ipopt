// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 18, 2008

#ifndef INCLUDE_CALLBACKFUNCTIONS
#define INCLUDE_CALLBACKFUNCTIONS

#include "mex.h"
#include "iterate.hpp"
#include "sparsematrix.hpp"
#include "matlabfunctionhandle.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"

// Class CallbackFunctions.
// -----------------------------------------------------------------
// An object of this class does two things. First of all, it stores
// handles to MATLAB functions (type HELP FUNCTION_HANDLE in the
// MATLAB console) for all the necessary and optional callback
// functions for IPOPT. Secondly, this class actually provides the
// routines for calling these functions with the necessary inputs and
// outputs.
class CallbackFunctions {
public:

  // The constructor must be provided with a single MATLAB array, and
  // this MATLAB array must be a structure array in which each field
  // is a function handle.
  explicit CallbackFunctions (const mxArray* ptr);

  // The destructor.
  ~CallbackFunctions();

  // These functions return true if the respective callback functions
  // are available.
  bool constraintFuncIsAvailable() const { return *constraintfunc; };
  bool jacobianFuncIsAvailable  () const { return *jacobianfunc;   };
  bool hessianFuncIsAvailable   () const { return *hessianfunc;    };
  bool iterFuncIsAvailable      () const { return *iterfunc;       };

  // These functions execute the various callback functions with the
  // appropriate inputs and outputs. Here, m is the number of constraints.
  // The first function returns the value of the objective at x.
  double computeObjective (const Iterate& x) const;

  // This function computes the value of the gradient at x, and
  // returns the gradient entries in the array g, which must be of
  // length equal to the number of optimization variables.
  void computeGradient (const Iterate& x, double* g) const;

  // This function computes the response of the vector-valued
  // constraint function at x, and stores the result in the array c
  // which must be of length m.
  void computeConstraints (const Iterate& x, int m, double* c) const;

  // This function gets the structure of the sparse m x n Jacobian matrix.
  SparseMatrix* getJacobianStructure (int n, int m) const;

  // This function gets the structure of the sparse n x n Hessian matrix.
  SparseMatrix* getHessianStructure (int n) const;

  // This function computes the Jacobian of the constraints at x.
  void computeJacobian (int m, const Iterate& x, SparseMatrix& J) const;

  // This function computes the Hessian of the Lagrangian at x.
  void computeHessian (const Iterate& x, double sigma, int m, 
		       const double* lambda, SparseMatrix& H) const;

  // Call the intermediate callback function. A return value of false
  // tells IPOPT to terminate.
  bool iterCallback (int t, double f, 
				      double inf_pr, double inf_du, 
				      double mu, double d_norm,
				      double regularization_size,
				      double alpha_du, double alpha_pr,
				      int ls_trials, const Ipopt::IpoptData* ip_data, 
				      Ipopt::IpoptCalculatedQuantities* ip_cq,
				      int n) const;

protected:
  MatlabFunctionHandle* objfunc;        // Objective callback function.
  MatlabFunctionHandle* gradfunc;       // Gradient callback function.
  MatlabFunctionHandle* constraintfunc; // Constraint callback function.
  MatlabFunctionHandle* jacobianfunc;   // Jacobian callback function.
  MatlabFunctionHandle* jacstrucfunc;   // Jacobian structure function.
  MatlabFunctionHandle* hessianfunc;    // Hessian callback function.
  MatlabFunctionHandle* hesstrucfunc;   // Hessian structure function.
  MatlabFunctionHandle* iterfunc;       // Iterative callback function.
};

#endif
