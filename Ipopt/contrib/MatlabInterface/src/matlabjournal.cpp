// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabjournal.h"
#include "matlabexception.h"
#include "mex.h"

namespace Ipopt {

  // Function definitions for class MatlabJournal.
  // ---------------------------------------------------------------
  MatlabJournal::MatlabJournal (EJournalLevel default_level)
    : Journal("matlab", default_level) { }

  void MatlabJournal::PrintImpl (EJournalCategory category, EJournalLevel level, const char* str) {
    mexPrintf(str);
  }

  void MatlabJournal::PrintfImpl (EJournalCategory category, EJournalLevel level, const char* pformat, va_list ap) {
    const int maxStrLen = 1024;
    char      s[maxStrLen];
    if (vsnprintf(s,maxStrLen,pformat,ap) >= maxStrLen)
      throw MatlabException("String buffer it too short for all the \
characters to be printed to MATLAB console");
    mexPrintf(s);
  }

  void MatlabJournal::FlushBufferImpl() { }
}
