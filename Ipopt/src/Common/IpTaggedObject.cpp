// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpTaggedObject.hpp"

#include <limits>

namespace Ipopt
{

/** Objects derived from TaggedObject MUST call this
 *  method every time their internal state changes to
 *  update the internal tag for comparison
 */
void TaggedObject::ObjectChanged()
{
   DBG_START_METH("TaggedObject::ObjectChanged()", 0);
   tag_ = unique_tag_;
   unique_tag_++;
   DBG_ASSERT(unique_tag_ < std::numeric_limits<Tag>::max());
   // The Notify method from the Subject base class notifies all
   // registered Observers that this subject has changed.
   Notify(Observer::NT_Changed);
}

TaggedObject::Tag IPOPT_THREAD_LOCAL TaggedObject::unique_tag_ = 1;

} // namespace Ipopt
