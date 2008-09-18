#include "matlabinfo.h"

// Function definitions for class MatlabInfo.
// ------------------------------------------------------------------
MatlabInfo::MatlabInfo (mxArray*& ptr) 
  : ptr(0) {

  // Create the structure array.
  const char* fieldnames[2];
  const char* exitstatusfield  = "status";
  const char* multipliersfield = "multipliers";
  fieldnames[0] = exitstatusfield;
  fieldnames[1] = multipliersfield;
  this->ptr = ptr = mxCreateStructMatrix(1,1,2,fieldnames);

  // Initialize the exit status field.
  mxSetField(ptr,0,"status",mxCreateDoubleScalar(0));
}

ApplicationReturnStatus MatlabInfo::getExitStatus() const {
  const mxArray* p = mxGetField(ptr,0,"status");
  return (ApplicationReturnStatus) *mxGetPr(p);  
}

void MatlabInfo::setExitStatus (ApplicationReturnStatus status) {
  mxArray* p = mxGetField(ptr,0,"status");
  *mxGetPr(p) = (double) status;
}
