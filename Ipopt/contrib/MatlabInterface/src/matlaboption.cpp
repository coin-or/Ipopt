// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlaboption.h"
#include "matlabscalar.h"
#include "matlabstring.h"

// Function definitions for class MatlabOption.
// -----------------------------------------------------------------
MatlabOption::MatlabOption (const mxArray* ptr) {
  s = 0;
  x = 0;

  if (mxIsChar(ptr)) {
    MatlabString str(ptr);
    s = new std::string((const char*) str);
  } else {
    MatlabScalar scalar(ptr);
    x = scalar;
  }
}

MatlabOption::~MatlabOption() {
  if (s)
    delete s;
}
