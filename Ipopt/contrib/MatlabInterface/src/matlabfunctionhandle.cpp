#include "matlabfunctionhandle.h"
#include "matlabexception.h"

// Function definitions for class MatlabFunctionHandle.
// -----------------------------------------------------------------
MatlabFunctionHandle::MatlabFunctionHandle() : f(0) { }

MatlabFunctionHandle::MatlabFunctionHandle (const mxArray* ptr)
  : f(0) {

  // First, check to see if the input arugment is an empty array. If
  // it is, we do nothing.
  if (!mxIsEmpty(ptr)) {

    // The array is not empty, so we are expecting a function handle.
    if (mxIsClass(ptr,"function_handle"))
      f = mxDuplicateArray(ptr);
    else
      throw MatlabException("Matlab array must be a function handle");
  }
}

MatlabFunctionHandle::MatlabFunctionHandle 
(const MatlabFunctionHandle& source) : f(source.f) {
  if (f)
    f = mxDuplicateArray(f);
}
  
MatlabFunctionHandle::~MatlabFunctionHandle() {
  if (f) 
    mxDestroyArray(f);
}
