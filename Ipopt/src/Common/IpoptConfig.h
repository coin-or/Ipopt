/* Copyright (C) 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * $Id$
 *
 * Include file for the configuration of Ipopt.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file, and
 * undefines macros that might configure with other Config.h files.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header files is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN-OR packages or third party code) are set
 * by the files config_*default.h. The project maintainer needs to remember
 * to update these file and choose reasonable defines.
 * A user can modify the default setting by editing the config_*default.h files.
 *
 */

#ifndef __IPOPTCONFIG_H__
#define __IPOPTCONFIG_H__

#ifdef HAVE_CONFIG_H

#ifdef IPOPT_BUILD
#include "config.h"
#else
#include "config_ipopt.h"
#endif

#else /* HAVE_CONFIG_H */

#ifdef IPOPT_BUILD
#include "config_default.h"
#else
#include "config_ipopt_default.h"
#endif

#endif /* HAVE_CONFIG_H */

/* If F77_FUNC has not been defined, then because no Fortran compiler was available.
 * However, Ipopt needs this define to call functions from library that are usually
 * written in Fortran (Blas, Lapack). We make up some definition here, but the user
 * may have to adapt for their system.
 */
#if defined(IPOPT_BUILD) && !defined(F77_FUNC)
#if 1
# define F77_FUNC(name,NAME) NAME
# define F77_FUNC_(name,NAME) NAME
#else
# define F77_FUNC(name,NAME) name ## _
# define F77_FUNC_(name,NAME) name ## __
#endif
#endif

#endif /*__IPOPTCONFIG_H__*/
