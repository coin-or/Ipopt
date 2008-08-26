// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         August 25, 2008

#include "callbackfunctions.h"
#include "matlabexception.h"

// Functions for class CallbackFunctions.
// -----------------------------------------------------------------
CallbackFunctions::CallbackFunctions (const mxArray* ptr) 
  : objfunc(0), gradfunc(0), constraintfunc(0), jacobianfunc(0), 
  jacstrucfunc(0), hessianfunc(0), hesstrucfunc(0), iterfunc(0) {
  const mxArray* p;  // A pointer to a MATLAB array.

  // Check whether we are provided with a structure array.
  if (!mxIsStruct(ptr))
    throw MatlabException("The second input must be a STRUCT");

  // Get the function handle for computing the objective.
  p = mxGetField(ptr,0,"objective");
  if (!p)
    throw MatlabException("You must specify a callback routine for \
computing the value of the objective function");
  if (mxIsEmpty(p) || !isFunctionHandle(p))
    throw MatlabException("You did not provide a valid function handle for \
computing the value of the objective function");
  objfunc = new MatlabFunctionHandle(p);

  // Get the function handle for computing the gradient.
  p = mxGetField(ptr,0,"gradient");
  if (!p)
    throw MatlabException("You must specify a callback routine for \
computing the gradient of the objective");
  if (mxIsEmpty(p) || !isFunctionHandle(p))
    throw MatlabException("You did not provide a valid function handle for \
computing the gradient of the objective");
  gradfunc = new MatlabFunctionHandle(p);  

  // Get the function handle for computing the constraints, if such a
  // function was specified.
  p = mxGetField(ptr,0,"constraints");
  if (p) {
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for computing the response of the constraints");
    constraintfunc = new MatlabFunctionHandle(p);      
  }
  else
    constraintfunc = new MatlabFunctionHandle();

  // Get the function handle for computing the Jacobian. This function
  // is necessary if there are constraints.
  p = mxGetField(ptr,0,"jacobian");
  if (p) {
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for computing the first derivatives (Jacobian) of the constraints");
    jacobianfunc = new MatlabFunctionHandle(p);      
  }
  else {
    if (*constraintfunc)
      throw MatlabException("You must provide a function that returns the \
first derivatives (Jacobian) of the constraints");
    jacobianfunc = new MatlabFunctionHandle();
  }

  // Get the function handle for computing the sparsity structure of
  // the Jacobian. This function is necessary if the Jacobian is being
  // computed.
  p = mxGetField(ptr,0,"jacobianstructure");
  if (p) { 
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for computing the sparsity structure of the Jacobian");
    jacstrucfunc = new MatlabFunctionHandle(p);      
  }
  else {
    if (*jacobianfunc)
      throw MatlabException("You must provide a function that returns the \
sparsity structure of the Jacobian");
    jacstrucfunc = new MatlabFunctionHandle();
  }

  // Get the function handle for computing the Hessian. This function
  // is always optional.
  p = mxGetField(ptr,0,"hessian");
  if (p) {
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for computing the Hessian of the Lagrangian");
    hessianfunc = new MatlabFunctionHandle(p);      
  }
  else 
    hessianfunc = new MatlabFunctionHandle();

  // Get the function handle for computing the sparsity structure of
  // the Hessian of the Lagrangian. This function is necessary if the
  // Hessian is being computed.
  p = mxGetField(ptr,0,"hessianstructure");
  if (p) { 
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for computing the sparsity structure of the Hessian");
    hesstrucfunc = new MatlabFunctionHandle(p);      
  }
  else {
    if (*hessianfunc)
      throw MatlabException("You must provide a function that returns the \
sparsity structure of the Hessian");
    hesstrucfunc = new MatlabFunctionHandle();
  }  

  // Get the iterative callback function handle. This function is
  // always optional.
  p = mxGetField(ptr,0,"iterfunc");
  if (p) { 
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for the iterative callback");
    iterfunc = new MatlabFunctionHandle(p);      
  }
  else
    iterfunc = new MatlabFunctionHandle();
}

CallbackFunctions::~CallbackFunctions() {
  if (objfunc)        delete objfunc;
  if (gradfunc)       delete gradfunc;
  if (constraintfunc) delete constraintfunc;
  if (jacobianfunc)   delete jacobianfunc;
  if (hessianfunc)    delete hessianfunc;
  if (iterfunc)       delete iterfunc;
}
