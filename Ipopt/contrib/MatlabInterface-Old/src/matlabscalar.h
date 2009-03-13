// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_MATLABSCALAR
#define INCLUDE_MATLABSCALAR

#include "mex.h"

// Class MatlabScalar
// -----------------------------------------------------------------
// The main appeal of this class is that one can create a scalar
// object that accesses a MATLAB array.
//
// Note that the copy assignment operator is not valid for this class
// because we cannot reassign a reference.
class MatlabScalar {
public:

  // This constructor accepts as input a pointer to a Matlab array
  // which must be a scalar in double precision.
  explicit MatlabScalar (const mxArray* ptr);

  // This constructor creates a new Matlab array which is a scalar
  // in double precision.
  MatlabScalar (mxArray*& ptr, double value);

  // The copy constructor.
  MatlabScalar (MatlabScalar& source);
  
  // The destructor.
  ~MatlabScalar() { };

  // Access the value of the scalar.
  operator const double () const { return x; };
    
  // Assign the value of the scalar.
  MatlabScalar& operator= (double value);
    
protected:
  double& x;

  // The copy assignment operator is kept protected because it is
  // invalid.
  MatlabScalar& operator= (const MatlabScalar& source) { return *this; };
};

#endif
