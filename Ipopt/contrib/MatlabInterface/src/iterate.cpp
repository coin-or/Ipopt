// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 15, 2008

#include "iterate.hpp"
#include "matlabexception.hpp"

// Function definitions for class Iterate.
// -----------------------------------------------------------------
Iterate::Iterate (mxArray* ptr) 
  : nv(0), ptr(ptr) {
  const mxArray* p = 0;  // Pointer to a MATLAB array.

  // Compute the number of optimization variables.
  if (mxIsCell(ptr)) {

    // The MATLAB array is a cell array. Repeat for each cell.
    int n = (int) mxGetNumberOfElements(ptr);
    for (int i = 0; i < n; i++) {
      p = mxGetCell(ptr,i);  // Get the ith cell.
      if (!mxIsDouble(p) || mxIsComplex(p) || mxIsSparse(p))
	throw MatlabException("The initial iterate must be either a full real array \
in DOUBLE precision, or a cell array in which each cell is a full real array in \
DOUBLE precision");      
      nv += (int) mxGetNumberOfElements(p);
    }
  } else {
    
    // The MATLAB array should be a numeric array.
    if (!mxIsDouble(ptr) || mxIsComplex(ptr) || mxIsSparse(ptr))
      throw MatlabException("The initial iterate must be either a full real array \
in DOUBLE precision, or a cell array in which each cell is a full real array in \
DOUBLE precision");
    nv = (int) mxGetNumberOfElements(ptr);
  }
}

void Iterate::inject (const double* x) {
  if (mxIsCell(ptr)) {

    // The MATLAB array is a cell array. Repeat for each cell.
    mxArray* p; 
    int      m;
    int      n = (int) mxGetNumberOfElements(ptr);
    for (int i = 0; i < n; i++) {
      p = mxGetCell(ptr,i);  // Get the ith cell.
      m = (int) mxGetNumberOfElements(p);
      copymemory(x,mxGetPr(p),m);
      x += m;
    }
  } else

    // The MATLAB array is a numeric array.
    copymemory(x,mxGetPr(ptr),nv);
}

void Iterate::copyto (double* x) const {
  if (mxIsCell(ptr)) {

    // The MATLAB array is a cell array. Repeat for each cell.
    const mxArray* p; 
    int            m;
    int            n = (int) mxGetNumberOfElements(ptr);
    for (int i = 0; i < n; i++) {
      p = mxGetCell(ptr,i);  // Get the ith cell.
      m = (int) mxGetNumberOfElements(p);
      copymemory(mxGetPr(p),x,m);
      x += m;
    }
  } else

    // The MATLAB array is a numeric array.
    copymemory(mxGetPr(ptr),x,nv);
}

// Function definitions for static members of class Iterate.
// -----------------------------------------------------------------
int Iterate::getMatlabData (const mxArray* ptr, double*& data) {
  Iterate x(mxDuplicateArray(ptr)); // Create the iterate object.
  int     n = numvars(x);           // The return value.
  data = new double[n];
  x.copyto(data);
  mxDestroyArray(x.ptr);
  return n;
}
