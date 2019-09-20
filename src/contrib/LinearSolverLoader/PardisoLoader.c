/* Copyright (C) 2008 GAMS Development and others
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Author: Stefan Vigerske
 */

#include "IpoptConfig.h"
#include "LibraryHandler.h"
#include "PardisoLoader.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/* Type of Fortran integer translated into C */
typedef FORTRAN_INTEGER_TYPE ipfint;

static soHandle_t Pardiso_handle = NULL;

void LSL_lateParadisoLibLoad();

typedef void (*voidfun)(void);

voidfun LSL_loadSym(
   soHandle_t  h,
   const char* symName,
   char*       msgBuf,
   int         msgLen
);

/* assuming PARDISO 4.0.0 and above */
typedef void (*pardisoinit_t)(
   void*         PT,
   const ipfint* MTYPE,
   const ipfint* SOLVER,
   ipfint*       IPARM,
   double*       DPARM,
   ipfint*       E
);

typedef void (*pardiso_t)(
   void**        PT,
   const ipfint* MAXFCT,
   const ipfint* MNUM,
   const ipfint* MTYPE,
   const ipfint* PHASE,
   const ipfint* N,
   const double* A,
   const ipfint* IA,
   const ipfint* JA,
   const ipfint* PERM,
   const ipfint* NRHS,
   ipfint*       IPARM,
   const ipfint* MSGLVL,
   double*       B,
   double*       X,
   ipfint*       E,
   double*       DPARM
);

static pardisoinit_t func_pardisoinit = NULL;
static pardiso_t     func_pardiso     = NULL;
static int           pardiso_is_parallel  = 0;

void pardisoinit(
   void*         PT,
   const ipfint* MTYPE,
   const ipfint* SOLVER,
   ipfint*       IPARM,
   double*       DPARM,
   ipfint*       E
)
{
   if( func_pardisoinit == NULL )
   {
      LSL_lateParadisoLibLoad();
   }
   assert(func_pardisoinit != NULL);

   func_pardisoinit(PT, MTYPE, SOLVER, IPARM, DPARM, E);
}

void pardiso(
   void**        PT,
   const ipfint* MAXFCT,
   const ipfint* MNUM,
   const ipfint* MTYPE,
   const ipfint* PHASE,
   const ipfint* N,
   const double* A,
   const ipfint* IA,
   const ipfint* JA,
   const ipfint* PERM,
   const ipfint* NRHS,
   ipfint*       IPARM,
   const ipfint* MSGLVL,
   double*       B,
   double*       X,
   ipfint*       E,
   double*       DPARM
)
{
   if (func_pardiso == NULL)
   {
      LSL_lateParadisoLibLoad();
   }
   assert(func_pardiso != NULL);

   /* if we do not have a parallel version, ensure that IPARM[2] (#threads) is set to 1 */
   if (!pardiso_is_parallel)
   {
      IPARM[2] = 1;
   }

   func_pardiso(PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, E, DPARM);
}

#define PARDISOLIBNAME "libpardiso." SHAREDLIBEXT

int LSL_loadPardisoLib(
   const char* libname,
   char*       msgbuf,
   int         msglen
)
{
   /* load Pardiso library */
   if( libname )
   {
      Pardiso_handle = LSL_loadLib(libname, msgbuf, msglen);
   }
   else
   {
      /* try a default library name */
      Pardiso_handle = LSL_loadLib(PARDISOLIBNAME, msgbuf, msglen);
   }
   if( Pardiso_handle == NULL )
   {
      return 1;
   }

   /* load Pardiso functions, we assume the >= 4.0.0 interface */
   func_pardisoinit = (pardisoinit_t) LSL_loadSym(Pardiso_handle, "pardisoinit", msgbuf, msglen);
   if( func_pardisoinit == NULL )
   {
      return 1;
   }

   func_pardiso = (pardiso_t) LSL_loadSym(Pardiso_handle, "pardiso", msgbuf, msglen);
   if( func_pardiso == NULL )
   {
      return 1;
   }

   /* check if we use a parallel version of pardiso */
   pardiso_is_parallel = LSL_loadSym(Pardiso_handle, "pardiso_exist_parallel", msgbuf, msglen) != NULL;

   return 0;
}

int LSL_unloadPardisoLib()
{
   int rc;

   if( Pardiso_handle == NULL )
   {
      return 0;
   }

   rc = LSL_unloadLib(Pardiso_handle);
   Pardiso_handle = NULL;

   func_pardisoinit = NULL;
   func_pardiso = NULL;

   return rc;
}

int LSL_isPardisoLoaded()
{
   return Pardiso_handle != NULL;
}

void LSL_lateParadisoLibLoad()
{
   char buffer[512];
   int rc;

   sprintf(buffer, "Error unknown.");
   rc = LSL_loadPardisoLib(NULL, buffer, 512);
   if( rc != 0 )
   {
      fprintf(stderr, "Error loading Pardiso dynamic library " PARDISOLIBNAME ": %s\nAbort...\n", buffer);
      exit(EXIT_FAILURE);
   }
}

char* LSL_PardisoLibraryName()
{
   static char name[] = PARDISOLIBNAME;
   return name;
}
