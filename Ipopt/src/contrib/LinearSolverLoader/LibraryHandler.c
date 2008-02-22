/* Copyright (C) 2008   GAMS Development and others
 All Rights Reserved.
 This code is published under the Common Public License.

 $Id: LibraryHandler.c 341 2008-02-14 18:51:25Z stefan $

 Authors:  Stefan Vigerske

 copied from optcc.(h|c) in gams i/o libs
*/

#include "LibraryHandler.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#if defined(_WIN32)
#define snprintf _snprintf
#endif

soHandle_t LSL_loadLib(const char *libName, char *msgBuf, int msgLen)
{
  soHandle_t h=NULL;

#ifdef ERROR_LOADLIB
  snprintf(msgBuf, msgLen, "loadLib error: Do not know how to handle shared libraries on this operating system");
  return h;
#else

  if (libName==NULL) {
    snprintf(msgBuf, msgLen, "loadLib error: no library name given (libName is NULL)");
    return NULL;
  }

#ifdef HAVE_WINDOWS_H
  h = LoadLibrary (libName);
  if (NULL == h) {
    snprintf(msgBuf, msgLen, "Windows error while loading dynamic library %s", libName);
  }
#else
  h = dlopen (libName, RTLD_NOW);
  if (NULL == h) {
    strncpy(msgBuf, dlerror(), msgLen);
    msgBuf[msgLen-1]=0;
  }
#endif

  return h;
#endif
} /* LSL_loadLib */

int LSL_unloadLib (soHandle_t h)
{
  int rc=1;

#ifdef HAVE_WINDOWS_H
  rc = FreeLibrary (h);
  rc = ! rc;
#else
# ifdef HAVE_DLFCN_H
  rc = dlclose (h);
# endif
#endif
  return rc;
} /* LSL_unLoadLib */

void* LSL_loadSym (soHandle_t h, const char *symName, char *msgBuf, int msgLen)
{
  void *s;
  const char *from;
  char *to;
  const char *tripSym;
  char* err;
  char lcbuf[257];
  char ucbuf[257];
  char ocbuf[257];
  size_t symLen;
  int trip;

  /* search in this order:
   *  1. original
   *  2. lower_
   *  3. upper_
   *  4. original_
   *  5. lower
   *  6. upper
   */

  symLen = 0;
  for (trip = 1;  trip <= 6;  trip++) {
    switch (trip) {
    case 1:                             /* original */
      tripSym = symName;
      break;
    case 2:                             /* lower_ */
      for (from = symName, to = lcbuf;  *from;  from++, to++) {
        *to = tolower(*from);
      }
      symLen = from - symName;
      *to++ = '_';
      *to = '\0';
      tripSym = lcbuf;
      break;
    case 3:                             /* upper_ */
      for (from = symName, to = ucbuf;  *from;  from++, to++) {
        *to = toupper(*from);
      }
      *to++ = '_';
      *to = '\0';
      tripSym = ucbuf;
      break;
    case 4:                             /* original_ */
      memcpy (ocbuf, symName, symLen);
      ocbuf[symLen] = '_';
      ocbuf[symLen+1] = '\0';
      tripSym = ocbuf;
      break;
    case 5:                             /* lower */
      lcbuf[symLen] = '\0';
      tripSym = lcbuf;
      break;
    case 6:                             /* upper */
      ucbuf[symLen] = '\0';
      tripSym = ucbuf;
      break;
    default:
      tripSym = symName;
    } /* end switch */
#ifdef HAVE_WINDOWS_H
    s = GetProcAddress (h, tripSym);
    if (NULL != s) {
      return s;
    }
#else
# ifdef HAVE_DLFCN_H
    s = dlsym (h, tripSym);
    err = dlerror();  /* we have only one chance; a successive call to dlerror() returns NULL */ 
    if (err) {
    	strncpy(msgBuf, err, msgLen);
    	msgBuf[msgLen-1]=0;
    } else {
      return s;
    }
# else
    s = NULL;
    err = NULL;
# endif
#endif
  } /* end loop over symbol name variations */

#ifdef HAVE_WINDOWS_H
  snprintf(msgBuf, msgLen, "Cannot find symbol %s in dynamic library.", symName);
#endif

  return NULL;
} /* LSL_loadSym */
