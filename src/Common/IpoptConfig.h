/* Copyright (C) 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 */

/** Include file for the configuration of Ipopt.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header file is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN-OR packages or third party code) are set
 * by the files config_*default.h. The project maintainer needs to remember
 * to update these files and choose reasonable defines.
 * A user can modify the default setting by editing the config_*default.h files.
 */

#ifndef __IPOPTCONFIG_H__
#define __IPOPTCONFIG_H__

#ifdef HAVE_CONFIG_H

#ifdef IPOPTLIB_BUILD
#include "config.h"
#else
#include "config_ipopt.h"
#endif

#else /* HAVE_CONFIG_H */

#ifdef IPOPTLIB_BUILD
#include "config_default.h"
#else
#include "config_ipopt_default.h"
#endif

#endif /* HAVE_CONFIG_H */

/* overwrite XYZ_EXPORT from config.h when building XYZ
 * we want it to be __declspec(dllexport) when building a DLL on Windows
 * we want it to be __attribute__((__visibility__("default"))) when building with GCC,
 *   so user can compile with -fvisibility=hidden
 */
#ifdef IPOPTLIB_BUILD
# ifdef DLL_EXPORT
#  undef IPOPTLIB_EXPORT
#  define IPOPTLIB_EXPORT __declspec(dllexport)
# elif defined(__GNUC__) && __GNUC__ >= 4
#  undef IPOPTLIB_EXPORT
#  define IPOPTLIB_EXPORT __attribute__((__visibility__("default")))
# endif
#endif

#ifdef IPOPTAMPLINTERFACELIB_BUILD
# ifdef DLL_EXPORT
#  undef IPOPTAMPLINTERFACELIB_EXPORT
#  define IPOPTAMPLINTERFACELIB_EXPORT __declspec(dllexport)
# elif defined(__GNUC__) && __GNUC__ >= 4
#  undef IPOPTAMPLINTERFACELIB_EXPORT
#  define IPOPTAMPLINTERFACELIB_EXPORT __attribute__((__visibility__("default")))
# endif
#endif

#ifdef SIPOPTLIB_BUILD
# ifdef DLL_EXPORT
#  undef SIPOPTLIB_EXPORT
#  define SIPOPTLIB_EXPORT __declspec(dllexport)
# elif defined(__GNUC__) && __GNUC__ >= 4
#  undef SIPOPTLIB_EXPORT
#  define SIPOPTLIB_EXPORT __attribute__((__visibility__("default")))
# endif
#endif

#endif /*__IPOPTCONFIG_H__*/
