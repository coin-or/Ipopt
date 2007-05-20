// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_MATLABSTRING
#define INCLUDE_MATLABSTRING

#include "mex.h"
#include <string>

// Function declarations.
// -----------------------------------------------------------------
// Copy a C-style string (i.e. a null-terminated character array).
char* copystring (const char* source);

// Class MatlabString.
// -----------------------------------------------------------------
// This class encapsulates read-only access to a MATLAB character
// array.
class MatlabString {
public:

  // The constructor accepts as input a pointer to a Matlab array,
  // which must be a valid string (array of type CHAR).
  explicit MatlabString (const mxArray* ptr);
  
  // This constructor accepts as input a null-terminated string.
  explicit MatlabString (const char* s);

  // The copy constructor makes a full copy of the source string.
  MatlabString (const MatlabString& source);
  
  // The destructor.
  ~MatlabString();

  // Return true if the string is empty.
  bool isempty () const { return strlen(s) == 0; };
  
  // Conversion operator for null-terminated string.
  operator const char* () const { return s; };
  
  // Conversion operator for string object.
  operator std::string () const { return std::string(s); };

protected:
  char* s;  // The null-terminated string.
  
  // The copy assignment operator is not proper, thus remains
  // protected.
  MatlabString& operator= (const MatlabString& source) { return *this; };
};

#endif
