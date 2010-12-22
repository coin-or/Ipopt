// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_MATLABOPTION
#define INCLUDE_MATLABOPTION

#include "mex.h"
#include <string>

// Class MatlabOption.
// -----------------------------------------------------------------
class MatlabOption {
public:

  // This constructor creates an object starting from a Matlab
  // array. The Matlab array must either be a string or a scalar in
  // double precision.
  MatlabOption (const mxArray* ptr);

  // Return "true" if the option value is a string.
  bool isString() const { return s; };

  // Get the option value.
  operator const char*        () const { return s->c_str(); };
  operator const std::string& () const { return *s;         };
  operator double             () const { return x;          };
  operator int                () const { return (int) x;    };

  // The destructor.
  ~MatlabOption();

protected:
  std::string* s;  // The string value.
  double       x;  // The double value.
};

#endif
