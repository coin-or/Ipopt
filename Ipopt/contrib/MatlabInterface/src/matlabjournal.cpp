// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabjournal.hpp"
#include "matlabexception.hpp"
#include "mex.h"

namespace Ipopt {

  // Function definitions for class MatlabJournal.
  // ---------------------------------------------------------------
  MatlabJournal::MatlabJournal (EJournalLevel default_level)
    : Journal("matlab", default_level) { }

  void MatlabJournal::PrintImpl (EJournalCategory category, 
				 EJournalLevel level, const char* str) {
    mexPrintf(str);
  }

  void MatlabJournal::PrintfImpl (EJournalCategory category, 
				  EJournalLevel level, const char* pformat, 
				  va_list ap) {
    const int maxStrLen = 1024;
    char      s[maxStrLen];
#ifdef HAVE_VSNPRINTF
# ifdef HAVE_VA_COPY
    va_list apcopy;
    va_copy(apcopy, ap);
    if (vsnprintf(s,maxStrLen,pformat,apcopy) >= maxStrLen)
      throw MatlabException("String buffer it too short for all the \
characters to be printed to MATLAB console");
    va_end(apcopy);
# else
    if (vsnprintf(s,maxStrLen,pformat,ap) >= maxStrLen)
      throw MatlabException("String buffer it too short for all the \
characters to be printed to MATLAB console");
# endif
#else
# ifdef HAVE__VSNPRINTF
#  ifdef HAVE_VA_COPY
    va_list apcopy;
    va_copy(apcopy, ap);
    if (_vsnprintf(s,maxStrLen,pformat,apcopy) >= maxStrLen)
      throw MatlabException("String buffer it too short for all the \
characters to be printed to MATLAB console");
    va_end(apcopy);
#  else
    if (_vsnprintf(s,maxStrLen,pformat,ap) >= maxStrLen)
      throw MatlabException("String buffer it too short for all the \
characters to be printed to MATLAB console");
#  endif
# else
    vsprintf(s,pformat,ap);
# endif
#endif
    mexPrintf(s);
  }

  void MatlabJournal::FlushBufferImpl() { }
}
