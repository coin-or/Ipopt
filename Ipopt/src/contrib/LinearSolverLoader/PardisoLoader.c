/* Copyright (C) 2008 GAMS Development and others
 All Rights Reserved.
 This code is published under the Eclipse Public License.

 $Id$

 Author: Stefan Vigerske
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

voidfun LSL_loadSym (soHandle_t h, const char *symName, char *msgBuf, int msgLen);

/* Old pre-4.0.0 pardiso interface */
typedef void (*pardisoinit_old_t)(void* PT, const ipfint* MTYPE,
                           ipfint* IPARM);
typedef void (*pardiso_old_t)(void** PT, const ipfint* MAXFCT,
                           const ipfint* MNUM, const ipfint* MTYPE,
                           const ipfint* PHASE, const ipfint* N,
                           const double* A, const ipfint* IA,
                           const ipfint* JA, const ipfint* PERM,
                           const ipfint* NRHS, ipfint* IPARM,
                           const ipfint* MSGLVL, double* B, double* X,
                           ipfint* E);

/* PARDISO 4.0.0 and above */
typedef void (*pardisoinit_new_t)(void* PT, const ipfint* MTYPE,
                              const ipfint * SOLVER,
                              ipfint* IPARM,
                              double* DPARM,
                              ipfint* E);
typedef void (*pardiso_new_t)(void** PT, const ipfint* MAXFCT,
                           const ipfint* MNUM, const ipfint* MTYPE,
                           const ipfint* PHASE, const ipfint* N,
                           const double* A, const ipfint* IA,
                           const ipfint* JA, const ipfint* PERM,
                           const ipfint* NRHS, ipfint* IPARM,
                           const ipfint* MSGLVL, double* B, double* X,
                           ipfint* E, double* DPARM);

static pardisoinit_old_t func_pardisoinit     = NULL;
static pardisoinit_new_t func_new_pardisoinit = NULL;
static pardiso_old_t func_pardiso     = NULL;
static pardiso_new_t func_new_pardiso = NULL;
static int pardiso_is_parallel = 0;

void wrap_old_pardisoinit(void* PT, const ipfint* MTYPE, const ipfint* SOLVER, ipfint* IPARM, double* DPARM, ipfint* E) {
   if (func_pardisoinit == NULL)
      LSL_lateParadisoLibLoad();
   func_pardisoinit(PT, MTYPE, IPARM);
   *E = 0;
}

void wrap_old_pardiso(void** PT, const ipfint* MAXFCT,
                      const ipfint* MNUM, const ipfint* MTYPE,
                      const ipfint* PHASE, const ipfint* N,
                      const double* A, const ipfint* IA,
                      const ipfint* JA, const ipfint* PERM,
                      const ipfint* NRHS, ipfint* IPARM,
                      const ipfint* MSGLVL, double* B, double* X,
                      ipfint* E, double* DPARM) {
   /* Note: we assume dparm is not of importance (only used for indirect solver
    * according to PARDISO 4.1.1 documentation). */
  if (func_pardiso == NULL)
     LSL_lateParadisoLibLoad();
  /* if we do not have a parallel version, ensure that IPARM[2] (#threads) is set to 1 */
  if (!pardiso_is_parallel)
     IPARM[2] = 1;
   func_pardiso(PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS,
      IPARM, MSGLVL, B, X, E);
}

void F77_FUNC(pardisoinit,PARDISOINIT)(void* PT, const ipfint* MTYPE,
                      const ipfint* SOLVER,
                      ipfint* IPARM,
                      double* DPARM,
                      ipfint* E) {
   if (func_new_pardisoinit == NULL)
      LSL_lateParadisoLibLoad();
   assert(func_new_pardisoinit != NULL);
   func_new_pardisoinit(PT, MTYPE, SOLVER, IPARM, DPARM, E);
}

void F77_FUNC(pardiso,PARDISO)(void** PT, const ipfint* MAXFCT,
                  const ipfint* MNUM, const ipfint* MTYPE,
                  const ipfint* PHASE, const ipfint* N,
                  const double* A, const ipfint* IA,
                  const ipfint* JA, const ipfint* PERM,
                  const ipfint* NRHS, ipfint* IPARM,
                  const ipfint* MSGLVL, double* B, double* X,
                  ipfint* E, double* DPARM) {
   if (func_new_pardiso == NULL)
      LSL_lateParadisoLibLoad();
   assert(func_new_pardiso != NULL);
   /* if we do not have a parallel version, ensure that IPARM[2] (#threads) is set to 1 */
   if (!pardiso_is_parallel)
      IPARM[2] = 1;
   func_new_pardiso(PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, E, DPARM);
}

#define PARDISOLIBNAME "libpardiso." SHAREDLIBEXT

int LSL_loadPardisoLib(const char* libname, char* msgbuf, int msglen) {
  /* load Pardiso library */
  if (libname) {
    Pardiso_handle=LSL_loadLib(libname, msgbuf, msglen);
  } else { /* try a default library name */
    Pardiso_handle=LSL_loadLib(PARDISOLIBNAME, msgbuf, msglen); 
  }
  if (Pardiso_handle==NULL)
    return 1;
	
  /* load Pardiso functions */
  /* first check if we are new or old interface */
  if(LSL_loadSym(Pardiso_handle, "pardiso_ipopt_newinterface", msgbuf, msglen) != NULL) {
     func_new_pardisoinit=(pardisoinit_new_t)LSL_loadSym(Pardiso_handle, "pardisoinit", msgbuf, msglen);
     if (func_new_pardisoinit == NULL) return 1;

     func_new_pardiso=(pardiso_new_t)LSL_loadSym(Pardiso_handle, "pardiso", msgbuf, msglen);
     if (func_new_pardiso == NULL) return 1;

  } else {
     func_pardisoinit=(pardisoinit_old_t)LSL_loadSym(Pardiso_handle, "pardisoinit", msgbuf, msglen);
     if (func_pardisoinit == NULL) return 1;

     func_pardiso=(pardiso_old_t)LSL_loadSym(Pardiso_handle, "pardiso", msgbuf, msglen);
     if (func_pardiso == NULL) return 1;

     func_new_pardisoinit = wrap_old_pardisoinit;
     func_new_pardiso = wrap_old_pardiso;
  }

  /* check if we use a parallel version of pardiso */
  pardiso_is_parallel = LSL_loadSym(Pardiso_handle, "pardiso_exist_parallel", msgbuf, msglen) != NULL;

  return 0;
}

int LSL_unloadPardisoLib() {
  int rc;

  if (Pardiso_handle==NULL)
    return 0;

  rc = LSL_unloadLib(Pardiso_handle);
  Pardiso_handle=NULL;

  func_pardisoinit=NULL;
  func_pardiso=NULL;

  return rc;
}

int LSL_isPardisoLoaded() {
  return Pardiso_handle!=NULL;
}

void LSL_lateParadisoLibLoad() {
  char buffer[512];
  int rc;

  sprintf(buffer, "Error unknown.");
  rc = LSL_loadPardisoLib(NULL, buffer, 512);
  if (rc!=0) {
    fprintf(stderr, "Error loading Pardiso dynamic library " PARDISOLIBNAME ": %s\nAbort...\n", buffer);
    exit(EXIT_FAILURE);
  }
}

char* LSL_PardisoLibraryName() {
  static char name[] = PARDISOLIBNAME;
  return name;
}
