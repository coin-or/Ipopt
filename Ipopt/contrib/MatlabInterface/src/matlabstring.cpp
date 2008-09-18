// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabstring.h"
#include "matlabexception.h"

// Function declarations.
// ----------------------------------------------------------------
// Copy a C-style string (i.e. a null-terminated character array).
char* copystring (const char* source);

// Function definitions for class MatlabString.
// -----------------------------------------------------------------
MatlabString::MatlabString (const mxArray* ptr) 
  : s(0) {

  // Get the string passed as a Matlab array.
  s = mxArrayToString(ptr);
  if (s == 0)
    throw MatlabException("Unable to obtain string from Matlab array");
}

MatlabString::MatlabString (const char* s) {
  this->s = copystring(s);
}

MatlabString::MatlabString (const MatlabString& source) {
  s = copystring(source.s);
}

MatlabString::~MatlabString() { 
  if (s) mxFree(s); 
}

// Function definitions.
// -----------------------------------------------------------------
char* copystring (const char* source) {
  int   n    = strlen(source);
  char* dest = (char*) mxMalloc((n+1)*sizeof(char));
  strcpy(dest,source);
  return dest;
}

