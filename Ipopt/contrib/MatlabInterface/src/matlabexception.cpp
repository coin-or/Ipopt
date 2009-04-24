// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabexception.hpp"
#include <cstdio>

// Function definitions for class MatlabException
// ---------------------------------------------------------------
MatlabException::MatlabException (const char* message) throw()
  : exception() { 
  snprintf(this->message, ME_BUFLEN, "\n*** Error using Ipopt Matlab interface: ***\n%s.\n", message);
}

MatlabException::MatlabException (const MatlabException& source) throw() 
  : exception() {
  snprintf(message, ME_BUFLEN, "%s.", source.message);
}

MatlabException& MatlabException::operator= (const MatlabException& source) 
{ 
  snprintf(message, ME_BUFLEN, "%s", source.message);
  return *this;
}
