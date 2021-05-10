// Copyright (C) 2020 COIN-OR Foundation
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef __IPTYPES_H__
#define __IPTYPES_H__

#include "IpoptConfig.h"

#ifndef IPOPT_DEPRECATED
#if defined(_MSC_VER)
/** macro to declare symbols as deprecated
 * @since Ipopt 3.14.0
 */
#  define IPOPT_DEPRECATED __declspec(deprecated)
#elif defined(__GNUC__)
/** macro to declare symbols as deprecated
 * @since Ipopt 3.14.0
 */
#  define IPOPT_DEPRECATED __attribute__ ((deprecated))
#else
/** macro to declare symbols as deprecated
 * @since Ipopt 3.14.0
 */
#  define IPOPT_DEPRECATED
#endif
#endif

#ifdef __GNUC__
#  define IPOPT_UNUSED __attribute__((unused))
#else
#  define IPOPT_UNUSED
#endif

#ifndef IPOPT_CALLCONV
#ifdef _MSC_VER
/** @since Ipopt 3.14.0 */
#define IPOPT_CALLCONV __cdecl
#else
/** @since Ipopt 3.14.0 */
#define IPOPT_CALLCONV
#endif
#endif

/** Type for floating-point numbers
 * @since 3.14.0
 */
#ifdef IPOPT_SINGLE
typedef float ipnumber;
#else
typedef double ipnumber;
#endif

#ifdef IPOPT_INT64
#include <inttypes.h>
/** Type of all indices of vectors, matrices etc
 * @since 3.14.0
 */
typedef int64_t ipindex;
/** Format specifier to use for ipindex in printf-format strings
 * @since 3.14.0
 */
#define IPOPT_INDEX_FORMAT PRId64
#else
/** Type of all indices of vectors, matrices etc
 * @since 3.14.0
 */
typedef int ipindex;
/** Format specifier to use for ipindex in printf-format strings
 * @since 3.14.0
 */
#define IPOPT_INDEX_FORMAT "d"
#endif

/** Type of Fortran integer translated into C
 * @deprecated Use ipindex instead.
 */
IPOPT_DEPRECATED
typedef ipindex ipfint;

#endif
