/* Copyright (C) 2008   GAMS Development
 All Rights Reserved.
 This code is published under the Common Public License.

 $Id: LibraryHandler.c 341 2008-02-14 18:51:25Z stefan $

 Authors:  Stefan Vigerske

 copied from optcc.(h|c) in gams i/o libs
*/

#include "LibraryHandler.h"

#if defined(_WIN32) || defined(BUILD_TYPE_WINDOWS)
#define snprintf _snprintf
#endif


soHandle_t LSL_loadLib(const char *libName, char *msgBuf, int msgLen)
{
	soHandle_t h=NULL;
/* no HP support yet
#if defined(CIA_HP7)
  int flag = 0;
#endif
*/
	
	if (libName==NULL) {
		snprintf(msgBuf, msgLen, "loadLib error: no library name given (libName is NULL)");
		return NULL;
	}

#if defined(_WIN32) || defined(BUILD_TYPE_WINDOWS)
  h = LoadLibrary (libName);
  if (NULL == h) {
  	snprintf(msgBuf, msgLen, "Windows error while loading dynamic library %s", libName);
  }
/* no HP support yet 
#elif defined(CIA_HP7)
  flag = BIND_IMMEDIATE | BIND_VERBOSE | DYNAMIC_PATH;
  h = shl_load (libName, flag, 0L);
  if (NULL == h) {
  	strncpy(msgBuf, strerror(errno), msgLen);
  	msgBuf[msgLen-1]=0;
  }
*/
#else
  h = dlopen (libName, RTLD_NOW);
  if (NULL == h) {
  	strncpy(msgBuf, dlerror(), msgLen);
  	msgBuf[msgLen-1]=0;
  }
#endif

  return h;
} /* LSL_loadLib */

int LSL_unloadLib (soHandle_t h)
{
  int rc;

#if defined(_WIN32) || defined(BUILD_TYPE_WINDOWS)
  rc = FreeLibrary (h);
  rc = ! rc;
/* no HP support yet
#elif defined(CIA_HP7)
  rc = shl_unload (h);
*/
#else
  rc = dlclose (h);
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
/* no HP support yet
#if defined(CIA_HP7)
  int rc;
#endif
*/

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
#if defined(_WIN32) || defined(BUILD_TYPE_WINDOWS)
    s = GetProcAddress (h, tripSym);
    if (NULL != s) {
      return s;
    }
/* no HP support yet
#elif defined(CIA_HP7)
    rc = shl_findsym (&h, tripSym, TYPE_UNDEFINED, &s);
    if (rc) {                     
      strncpy(msgBuf, strerror(errno), msgLen);
      msgBuf[msgLen-1]=0;
    }
    else {                        
      return s;
    }
*/
#else
    s = dlsym (h, tripSym);
    err = dlerror();  /* we have only one chance; a successive call to dlerror() returns NULL */ 
    if (err) {
    	strncpy(msgBuf, err, msgLen);
    	msgBuf[msgLen-1]=0;
    } else {
      return s;
    }
#endif
  } /* end loop over symbol name variations */

#if defined(_WIN32) || defined(BUILD_TYPE_WINDOWS)
	snprintf(msgBuf, msgLen, "Cannot find symbol %s in dynamic library.", symName);
#endif

  return NULL;
} /* LSL_loadSym */
