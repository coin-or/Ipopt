// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPTYPES_HPP__
#define __IPTYPES_HPP__

#include "IpoptConfig.h"
#include "IpTypes.h"

namespace Ipopt
{

/** Type of all numbers */
typedef ipnumber Number;

/** Type of all indices of vectors, matrices etc */
typedef ipindex Index;

/** Type of default integer
 * @deprecated Use int instead.
 */
IPOPT_DEPRECATED
typedef int Int;

} // namespace Ipopt

#endif
