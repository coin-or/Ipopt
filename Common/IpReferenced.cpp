// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpReferenced.hpp"

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  ReferencedObject::ReferencedObject()
      :
      reference_count_(0)
  {
    //    DBG_START_METH("ReferencedObject::ReferencedObject()", dbg_verbosity);
  }

  ReferencedObject::~ReferencedObject()
  {
    //    DBG_START_METH("ReferencedObject::~ReferencedObject()", dbg_verbosity);
    DBG_ASSERT(reference_count_ == 0);
  }


} // namespace Ipopt
