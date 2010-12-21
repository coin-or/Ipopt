// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         August 25, 2008

#ifndef INCLUDE_MATLABFUNCTIONHANDLE
#define INCLUDE_MATLABFUNCTIONHANDLE

#include "mex.h"

// Function declarations.
// -----------------------------------------------------------------
// This function returns true if and only if the MATLAB array is a
// valid function handle.
bool isFunctionHandle (const mxArray* ptr);

// Class MatlabFunctionHandle.
// -----------------------------------------------------------------
// The purpose of this class is twofold. The first aim is to store
// information about a MATLAB function handle. (For more information
// on function handles in MATLAB, type HELP FUNCTION_HANDLE in the
// MATLAB console). The second purpose is to provide a routine for
// evaluating the response of the function, provided inputs to the
// function.
class MatlabFunctionHandle {
public:

  // The default constructor creates a null function handle.
  MatlabFunctionHandle() : f(0) { };

  // This constructor accepts as input a pointer to a MATLAB array. It
  // is up to the user to ensure that the MATLAB array is a valid
  // function handle.
  explicit MatlabFunctionHandle (const mxArray* ptr);
  
  // The destructor.
  ~MatlabFunctionHandle();

  // This method is used to call the MATLAB function, provided inputs
  // to the function. It is up to the user to make sure that the
  // outputs array is of the appropriate size. The function returns
  // "true" on success, or "false" on failure. It is up to the user to
  // properly deallocate the outputs.
  bool evaluate (int nin, int nout, const mxArray** inputs, mxArray** outputs);

  // Returns true if and only if the function handle is not null.
  operator bool() const { return f != 0; };

protected:
  mxArray* f;  // The MATLAB function handle.
};

#endif
