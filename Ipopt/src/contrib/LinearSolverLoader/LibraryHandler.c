/* Copyright (C) 2008 GAMS Development and others
 All Rights Reserved.
 This code is published under the Eclipse Public License.

 $Id$

 Author: Stefan Vigerske

 inspired by optcc.c in gams i/o libs
*/

#include "LibraryHandler.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#ifdef HAVE_SNPRINTF
#define mysnprintf snprintf
#else
# ifdef HAVE__SNPRINTF
# define mysnprintf _snprintf
# else
#  define mysnprintf snprintf
#  error "Do not have function for save printing into a C-string (snprintf or _snprintf)"
# endif
#endif

soHandle_t LSL_loadLib(const char *libName, char *msgBuf, int msgLen)
{
  soHandle_t h=NULL;

#ifdef ERROR_LOADLIB
  mysnprintf(msgBuf, msgLen, "loadLib error: Do not know how to handle shared libraries on this operating system");
  return h;
#else

  if (libName==NULL) {
    mysnprintf(msgBuf, msgLen, "loadLib error: no library name given (libName is NULL)");
    return NULL;
  }

# ifdef HAVE_WINDOWS_H
  h = LoadLibrary (libName);
  if (NULL == h) {
    mysnprintf(msgBuf, msgLen, "Windows error while loading dynamic library %s, error = %d.\n(see http://msdn.microsoft.com/en-us/library/ms681381%%28v=vs.85%%29.aspx)\n", libName, GetLastError());
  }
# else
  h = dlopen (libName, RTLD_NOW);
  if (NULL == h) {
    strncpy(msgBuf, dlerror(), msgLen);
    msgBuf[msgLen-1]=0;
  }
# endif

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

#ifdef HAVE_WINDOWS_H
typedef FARPROC symtype;
#else
typedef void* symtype; 
#endif
/** Loads a symbol from a dynamically linked library.
 * This function is not defined in the header to allow a workaround for the problem that dlsym returns an object instead of a function pointer.
 * However, Windows also needs special care.
 * 
 * The method does six attempts to load the symbol. Next to its given name, it also tries variations of lower case and upper case form and with an extra underscore.  
 * @param h Handle of dynamicall linkes library.
 * @param symName Name of the symbol to load.
 * @param msgBuf Buffer for error messages, assumed to be NOT NULL!
 * @param msgLen Length of message buffer.
 * @return A pointer to the symbol, or NULL if not found.  
 */
symtype LSL_loadSym (soHandle_t h, const char *symName, char *msgBuf, int msgLen)
{
  symtype s;
  const char *from;
  char *to;
  const char *tripSym;
  char* err;
  char lcbuf[257];
  char ucbuf[257];
  char ocbuf[257];
  size_t symLen;
  int trip;
  
  s = NULL;
  err = NULL;

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
    } else {
      mysnprintf(msgBuf, msgLen, "Cannot find symbol %s in dynamic library, error = %d.", symName, GetLastError());
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
# endif
#endif
  } /* end loop over symbol name variations */

  return NULL;
} /* LSL_loadSym */
