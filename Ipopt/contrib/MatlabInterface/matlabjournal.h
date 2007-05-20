// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_MATLABJOURNAL
#define INCLUDE_MATLABJOURNAL

#include "ipopt/IpJournalist.hpp"

namespace Ipopt {

  // Class MatlabJournal.
  // ---------------------------------------------------------------
  // This class encapsulates journal output to the MATLAB console.
  class MatlabJournal : public Journal {
  public:

    // The constructor.
    MatlabJournal (EJournalLevel default_level);

    // The destructor.
    virtual ~MatlabJournal() { };

  protected:
    virtual void PrintImpl      (const char* str);
    virtual void PrintfImpl     (const char* pformat, va_list ap);
    virtual void FlushBufferImpl();
  };
}

#endif
