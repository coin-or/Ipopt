/* Copyright (C) 2008 GAMS Development and others
 All Rights Reserved.
 This code is published under the Eclipse Public License.

 $Id$

 Author: Stefan Vigerske

 inspired by optcc.h in gams i/o libs
*/

#ifndef LIBRARYHANDLER_H_
#define LIBRARYHANDLER_H_

#include "IpoptConfig.h"

#ifdef HAVE_WINDOWS_H
# include <windows.h>
  typedef HINSTANCE soHandle_t;
#ifdef small
#undef small
#endif
#else
# ifdef HAVE_DLFCN_H
#  include <unistd.h>
#  include <dlfcn.h>
  typedef void *soHandle_t;
# else
#  define ERROR_LOADLIB
  typedef void *soHandle_t;
# endif
#endif

/** Loads a dynamically linked library.
 * @param libname The name of the library to load.
 * @param msgbuf A buffer to store an error message.
 * @param msglen Length of the message buffer.
 * @return Shared library handle, or NULL if failure.
 */
soHandle_t LSL_loadLib(const char* libname, char* msgbuf, int msglen);

/** Unloads a shared library.
 * @param libhandle Handle of shared library to unload.
 * @return Zero on success, nonzero on failure.
 */
int LSL_unloadLib(soHandle_t libhandle);

#endif /*LIBRARYHANDLER_H_*/
