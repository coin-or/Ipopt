/* Copyright (C) 2008   GAMS Development
 All Rights Reserved.
 This code is published under the Common Public License.

 $Id: LibraryHandler.h 341 2008-02-14 18:51:25Z stefan $

 Authors:  Stefan Vigerske

 copied from optcc.(h|c) in gams i/o libs
*/

#ifndef LIBRARYHANDLER_H_
#define LIBRARYHANDLER_H_

#include "LinearSolverLoaderConfig.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>

#if defined(_WIN32) || defined(BUILD_TYPE_WINDOWS)
  typedef HINSTANCE soHandle_t;
/* no HP support yet
#elif defined(CIA_HP7)
# include <unistd.h>
# include <dl.h>
  typedef shl_t soHandle_t;
*/
#else
# include <unistd.h>
# include <dlfcn.h>
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

void* LSL_loadSym (soHandle_t h, const char *symName, char *msgBuf, int msgLen);

#endif /*LIBRARYHANDLER_H_*/
