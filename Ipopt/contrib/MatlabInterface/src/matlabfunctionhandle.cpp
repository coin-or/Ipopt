// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         August 25, 2008

#include "matlabfunctionhandle.hpp"
#include "matlabexception.hpp"

// Function definitions.
// -----------------------------------------------------------------
bool isFunctionHandle (const mxArray* ptr) {
  return mxIsClass(ptr,"function_handle");
}

// Function definitions for class MatlabFunctionHandle.
// -----------------------------------------------------------------
MatlabFunctionHandle::MatlabFunctionHandle (const mxArray* ptr) {
  if (mxIsEmpty(ptr))
    f = 0;
  else
    f = mxDuplicateArray(ptr);
}
  
MatlabFunctionHandle::~MatlabFunctionHandle() {
  if (f) mxDestroyArray(f);
}

bool MatlabFunctionHandle::evaluate (int nin, int nout, const mxArray** inputs,
				     mxArray** outputs) {
  // Construct the inputs to "feval".
  mxArray** finputs = new mxArray*[nin+1]; 
  finputs[0] = f;
  for (int i = 0; i < nin; i++)
    finputs[i+1] = mxDuplicateArray(inputs[i]);

  // Call "feval".
  int exitstatus = mexCallMATLAB(nout,outputs,nin+1,finputs,"feval");

  // Free the dynamically allocated memory.
  delete[] finputs;
  return exitstatus == 0;
}
