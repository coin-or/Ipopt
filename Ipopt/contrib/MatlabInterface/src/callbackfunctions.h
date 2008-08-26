// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         August 25, 2008

#ifndef INCLUDE_CALLBACKFUNCTIONS
#define INCLUDE_CALLBACKFUNCTIONS

#include "mex.h"
#include "matlabfunctionhandle.h"

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
  CallbackFunctions (const mxArray* ptr);

  // The destructor.
  ~CallbackFunctions();

  // These functions return true if the respective callback functions
  // are available.
  bool constraintFuncIsAvailable() const { return *constraintfunc; };
  bool jacobianFuncIsAvailable  () const { return *jacobianfunc;   };
  bool hessianFuncIsAvailable   () const { return *hessianfunc;    };
  bool iterFuncIsAvailable      () const { return *iterfunc;       };

  // These functions execute the various callback functions with the
  // appropriate inputs and outputs. The auxiliary data may be altered
  // over the course of executing the callback function. If there is
  // no auxiliary data, then simply pass in a null pointer.
  double computeObjective(const Iterate& x, mxArray*& auxdata) const;
  void   computeGradient (const Iterate& x, double* g, mxArray*& auxdata)const;
  void   computeConstrs  (const Iterate& x, double* c, mxArray*& auxdata)const;
  void   computeJacobian (const Iterate& x, double* J, mxArray*& auxdata)const;
  void   computeHessian  (const Iterate& x, double* H, mxArray*& auxdata)const;
  void   getJacobianStruc(mxArray*& auxdata)            const;
  void   getHessianStruc (mxArray*& auxdata)            const;


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
