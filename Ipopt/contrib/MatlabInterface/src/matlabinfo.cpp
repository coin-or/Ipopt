// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 25, 2008

#include "matlabinfo.hpp"
#include "iterate.hpp"

// Function definitions for class MatlabInfo.
// ------------------------------------------------------------------
MatlabInfo::MatlabInfo (mxArray*& ptr) 
  : ptr(0) {

  // Create the structure array.
  const char* fieldnames[6];
  const char* exitstatusfield = "status";
  const char* multlbfield     = "zl";
  const char* multubfield     = "zu";
  const char* multconstrfield = "lambda";
  const char* iterfield       = "iter";
  const char* cpu             = "cpu";
  fieldnames[0] = exitstatusfield;
  fieldnames[1] = multlbfield;
  fieldnames[2] = multubfield;
  fieldnames[3] = multconstrfield;
  fieldnames[4] = iterfield;
  fieldnames[5] = cpu;
  this->ptr = ptr = mxCreateStructMatrix(1,1,6,fieldnames);

  // Initialize some fields.
  mxSetField(ptr,0,"status",mxCreateDoubleScalar(0));
  mxSetField(ptr,0,"iter",mxCreateDoubleScalar(0));
  mxSetField(ptr,0,"cpu",mxCreateDoubleScalar(0));
}

ApplicationReturnStatus MatlabInfo::getExitStatus() const {
  const mxArray* p = mxGetField(ptr,0,"status");
  return (ApplicationReturnStatus) (int) *mxGetPr(p);  
}

void MatlabInfo::setExitStatus (ApplicationReturnStatus status) {
  mxArray* p = mxGetField(ptr,0,"status");
  *mxGetPr(p) = (double) status;
}

void MatlabInfo::setIterationCount (int iter) {
  mxArray* p = mxGetField(ptr,0,"iter");
  *mxGetPr(p) = (double) iter;
}

void MatlabInfo::setCpuTime (double cpu) {
  mxArray* p = mxGetField(ptr,0,"cpu");
  *mxGetPr(p) = cpu;
}

const double* MatlabInfo::getmultlb() const {
  mxArray* p = mxGetField(ptr,0,"zl");
  return mxGetPr(p);
}
const double* MatlabInfo::getmultub() const {
  mxArray* p = mxGetField(ptr,0,"zu");
  return mxGetPr(p);
}

const double* MatlabInfo::getmultconstr() const {
  mxArray* p = mxGetField(ptr,0,"lambda");
  return mxGetPr(p);
}

void MatlabInfo::setmultlb (int n, const double* zl) {

  // First destroy any previous multiplier values.
  mxArray* p = mxGetField(ptr,0,"zl");
  if (p) mxDestroyArray(p);

  // Set the field to the new multiplier values.
  p = mxCreateDoubleMatrix(n,1,mxREAL);
  copymemory(zl,mxGetPr(p),n);
  mxSetField(ptr,0,"zl",p);
}

void MatlabInfo::setmultub (int n, const double* zu) {

  // First destroy any previous multiplier values.
  mxArray* p = mxGetField(ptr,0,"zu");
  if (p) mxDestroyArray(p);

  // Set the field to the new multiplier values.
  p = mxCreateDoubleMatrix(n,1,mxREAL);
  copymemory(zu,mxGetPr(p),n);
  mxSetField(ptr,0,"zu",p);
}

void MatlabInfo::setmultconstr (int m, const double* lambda) {

  // First destroy any previous multiplier values.
  mxArray* p = mxGetField(ptr,0,"lambda");
  if (p) mxDestroyArray(p);

  // Set the field to the new multiplier values.
  p = mxCreateDoubleMatrix(m,1,mxREAL);
  copymemory(lambda,mxGetPr(p),m);
  mxSetField(ptr,0,"lambda",p);
}
