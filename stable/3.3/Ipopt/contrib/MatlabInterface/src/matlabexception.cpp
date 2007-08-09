// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabexception.h"

// Function definitions for class MatlabException
// ---------------------------------------------------------------
MatlabException::MatlabException (const char* message) throw()
  : exception() { 
  this->message = message;
}

MatlabException::MatlabException (const MatlabException& source) throw() 
  : exception() {
  message = source.message;
}

MatlabException& MatlabException::operator= (const MatlabException& source) 
{ 
  message = source.message; 
  return *this;
}
