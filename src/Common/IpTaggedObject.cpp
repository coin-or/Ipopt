// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpTaggedObject.hpp"

#include <limits>

/* keyword to declare a thread-local variable according to http://en.wikipedia.org/wiki/Thread-local_storage
 * GCC < 4.5 on MacOS X does not support TLS
 * With Intel compiler on MacOS X, problems with __thread were reported. (Hope that C++ thread_local will be OK)
 */
#ifndef IPOPT_THREAD_LOCAL

#if __cplusplus >= 201103L
#define IPOPT_THREAD_LOCAL thread_local
#elif defined(_MSC_VER)
#define IPOPT_THREAD_LOCAL __declspec(thread)
#elif defined(__APPLE__) && ((defined(__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__ < 405)) || defined(__INTEL_COMPILER))
#define IPOPT_THREAD_LOCAL
#else
#define IPOPT_THREAD_LOCAL __thread
#endif

#endif

namespace Ipopt
{

/** Global data that is incremented every time ANY TaggedObject changes.
 *
 * This allows us to obtain a unique Tag when the object changes.
 */
static IPOPT_THREAD_LOCAL TaggedObject::Tag unique_tag =  1;

/** Objects derived from TaggedObject MUST call this
 *  method every time their internal state changes to
 *  update the internal tag for comparison
 */
void TaggedObject::ObjectChanged()
{
   DBG_START_METH("TaggedObject::ObjectChanged()", 0);
   tag_ = unique_tag;
   unique_tag++;
   DBG_ASSERT(unique_tag < std::numeric_limits<Tag>::max());
   // The Notify method from the Subject base class notifies all
   // registered Observers that this subject has changed.
   Notify(Observer::NT_Changed);
}

} // namespace Ipopt
