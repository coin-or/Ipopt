/*
  Copyright (C) 2002, Carnegie Mellon University and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include "f2c_adolc.h"

#ifdef __sgi
#elif __linux__
#elif WIN32
#elif __osf__
#elif __sun
#elif AIX
#define model_init_ model_init
#else
#error Need to specify how characters are passed from C to Fortran
#endif

extern "C" {

#ifdef WIN32    // tell compiler that this function is to be exported
__declspec( dllexport )
#endif 
void model_init_(integer  *nz,
		 integer  *ny,
		 integer  *nu,
		 integer  *np,
		 double   *t,
		 double   *z,
		 double   *dz,
		 double   *y,
		 double   *u,
		 double   *p,
		 double   *f,
		 char     *flibname,
		 int      libnamelen)
{
  char *libname;
  int  len;
#ifdef WIN32
  HINSTANCE handle;
  typedef void (* Func)
    (integer *, integer *, integer *, integer *,
		 double *, double *, double *, double *, double *,
		 double *, double *, char *, int);
  Func MODEL_INIT;
  char EXT[]=".dll";
  LPVOID lpMsgBuf;
#else
  void *handle;
  void (*MODEL_INIT)
    (integer *, integer *, integer *, integer *,
		 double *, double *, double *, double *, double *,
		 double *, double *, char *, int);
  char EXT[]=".so";
  char *error;
#endif

  for( len = libnamelen; flibname[len-1]==' ' && len > 0; len-- );
  if( len == 0 )
    {
      fputs("Error in model_init_add: \nNo name for shared object containing model supplied.", stderr);
      exit(1);
    }

  libname= (char *)malloc(len+1+strlen(EXT));
  memcpy(libname,flibname,len);
  strcpy(&libname[len],EXT);

#ifdef WIN32
  handle = LoadLibrary(libname);
  if( !handle )
  {
    FormatMessage(
      FORMAT_MESSAGE_ALLOCATE_BUFFER | 
      FORMAT_MESSAGE_FROM_SYSTEM, NULL,
      GetLastError(), 0, (LPTSTR) &lpMsgBuf, 0, NULL);
    fprintf(stderr,"Error in model_init_add while opening DLL:\n%s",lpMsgBuf);
    LocalFree( lpMsgBuf );
    exit(1);
  }
  MODEL_INIT = (Func) GetProcAddress(handle,"model_init_");
  if( !MODEL_INIT )
  {
    FormatMessage(
      FORMAT_MESSAGE_ALLOCATE_BUFFER | 
      FORMAT_MESSAGE_FROM_SYSTEM, NULL,
      GetLastError(), 0, (LPTSTR) &lpMsgBuf, 0, NULL);
    fprintf(stderr,"Error in model_init_add while searching for function pointer:\n%s",lpMsgBuf);
    LocalFree( lpMsgBuf );
    exit(1);
  }
#else
  handle = dlopen(libname, RTLD_NOW);
  if (!handle)
    {
      fprintf(stderr, "Error in model_init_add while opening shared object:\n%s\n", dlerror());
      exit(1);
    }
  free(libname);

  MODEL_INIT = (void (*)(integer *, integer *, integer *, integer *,
		     double *, double *, double *, double *, double *,
		     double *, double *, char *, int))dlsym(handle, "model_init_");
  if ((error = dlerror()) != NULL) 
    {
      fprintf(stderr, "%s\n", error);
      exit(1);
    }
#endif

  (*MODEL_INIT)(nz, ny, nu, np, t, z, dz, y, u, p, f, flibname, libnamelen);

}

}
