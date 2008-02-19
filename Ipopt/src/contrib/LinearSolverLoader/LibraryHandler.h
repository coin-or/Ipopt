/* Copyright (C) 2008   GAMS Development and others
 All Rights Reserved.
 This code is published under the Common Public License.

 $Id: LibraryHandler.h 341 2008-02-14 18:51:25Z stefan $

 Authors:  Stefan Vigerske

 copied from optcc.(h|c) in gams i/o libs
*/

#ifndef LIBRARYHANDLER_H_
#define LIBRARYHANDLER_H_

#include "IpoptConfig.h"

#ifdef HAVE_WINDOWS_H
# include <windows.h>
  typedef HINSTANCE soHandle_t;
#elif defined(HAVE_DLFCN_H)
# include <unistd.h>
# include <dlfcn.h>
  typedef void *soHandle_t;
#else
# define ERROR_LOADLIB
  typedef void *soHandle_t;
#endif

/** Loads shared library.
 * @return Shared library handle, or NULL if failure.
 */
soHandle_t LSL_loadLib(const char* libname, char* msgbuf, int msglen);

/** Unloads shared library.
 * @return Zero on success, nonzero on failure.
 */
int LSL_unloadLib(soHandle_t libhandle);

#endif /*LIBRARYHANDLER_H_*/
