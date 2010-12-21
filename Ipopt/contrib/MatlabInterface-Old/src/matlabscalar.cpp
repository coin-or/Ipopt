// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabscalar.h"
#include "matlabexception.h"

// Function definitions.
// -----------------------------------------------------------------
double& getMatlabScalar (const mxArray* ptr) {
  if (!mxIsDouble(ptr))
    throw MatlabException("Matlab array must be of type double");
  if (mxGetNumberOfElements(ptr) != 1)
    throw MatlabException("The Matlab array must be a scalar");
  return *mxGetPr(ptr);
}

double& createMatlabScalar (mxArray*& ptr) {
  ptr = mxCreateDoubleScalar(0);
  return *mxGetPr(ptr);
}

// Function definitions for class MatlabScalar.
// -----------------------------------------------------------------
MatlabScalar::MatlabScalar (const mxArray* ptr) 
  : x(getMatlabScalar(ptr)) { }

MatlabScalar::MatlabScalar (mxArray*& ptr, double value) 
  : x(createMatlabScalar(ptr)) {
  x = value;
}

MatlabScalar::MatlabScalar (MatlabScalar& source) 
  : x(source.x) { }

MatlabScalar& MatlabScalar::operator= (double value) {
  x = value;
  return *this;
}
