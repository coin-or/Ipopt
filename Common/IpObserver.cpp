// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpObserver.hpp"

namespace Ipopt
{

  const Index Observer::dbg_verbosity = 0;
  const Index Subject::dbg_verbosity = 0;

  Observer::Observer()
  {
    DBG_START_METH("Observer::Observer", dbg_verbosity);
  }

  Observer::~Observer()
  {
    DBG_START_METH("Observer::~Observer", dbg_verbosity);
    // Detach all subjects
    for (Index i=0; i<(Index)subjects_.size(); i++) {
      DBG_PRINT((1,"subjects_[%d] = 0x%x\n", i, subjects_[i]));
    }
    for (Int i=(Int)(subjects_.size()-1); i>=0; i--) {
      DBG_PRINT((1, "Whatever\n"));
      DBG_ASSERT(i >= 0);
      if (i < 0) {
        DBG_PRINT((1,"I had better get here...\n"));
      }
      DBG_PRINT((1,"About to detach subjects_[%d] = 0x%x\n", i, subjects_[i]));
      RequestDetach(NT_All, subjects_[i]);
    }
  }

  Subject::Subject()
  {
    DBG_START_METH("Subject::Subject", dbg_verbosity);
  }

  Subject::~Subject()
  {
    DBG_START_METH("Subject::~Subject", dbg_verbosity);
    std::vector<Observer*>::iterator iter;
    for (iter = observers_.begin(); iter != observers_.end(); iter++) {
      (*iter)->ProcessNotification(Observer::NT_BeingDestroyed, this);
    }
  }

} // namespace Ipopt
