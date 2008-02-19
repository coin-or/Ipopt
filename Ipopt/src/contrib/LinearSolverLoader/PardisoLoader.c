/* Copyright (C) 2008   GAMS Development and others
 All Rights Reserved.
 This code is published under the Common Public License.

 $Id: PardisoLoader.c 341 2008-02-14 18:51:25Z stefan $

 Authors:  Stefan Vigerske
*/

#include "IpoptConfig.h"
#include "LibraryHandler.h"
#include "PardisoLoader.h"

#include <stdio.h>
#include <stdlib.h>

soHandle_t Pardiso_handle=NULL;

void LSL_lateParadisoLibLoad();

typedef void (*voidfun)(void);

voidfun LSL_loadSym (soHandle_t h, const char *symName, char *msgBuf, int msgLen);

typedef void (*pardisoinit_t)(void* PT, const ipfint* MTYPE, ipfint* IPARM);
typedef void (*pardiso_t)(void** PT, const ipfint* MAXFCT,
                           const ipfint* MNUM, const ipfint* MTYPE,
                           const ipfint* PHASE, const ipfint* N,
                           const double* A, const ipfint* IA,
                           const ipfint* JA, const ipfint* PERM,
                           const ipfint* NRHS, ipfint* IPARM,
                           const ipfint* MSGLVL, double* B, double* X,
                           ipfint* E);

pardisoinit_t func_pardisoinit=NULL;
pardiso_t func_pardiso=NULL;

void F77_FUNC(pardisoinit,PARDISOINIT)(void* PT, const ipfint* MTYPE, ipfint* IPARM) {
  if (func_pardisoinit==NULL) LSL_lateParadisoLibLoad(); 
  func_pardisoinit(PT, MTYPE, IPARM);
}

void F77_FUNC(pardiso,PARDISO)(void** PT, const ipfint* MAXFCT,
                               const ipfint* MNUM, const ipfint* MTYPE,
                               const ipfint* PHASE, const ipfint* N,
                               const double* A, const ipfint* IA,
                               const ipfint* JA, const ipfint* PERM,
                               const ipfint* NRHS, ipfint* IPARM,
                               const ipfint* MSGLVL, double* B, double* X,
                               ipfint* E) {
  if (func_pardiso==NULL) LSL_lateParadisoLibLoad(); 
  func_pardiso(PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, E);
}

#define PARDISOLIBNAME "libpardiso" SHAREDLIBEXT

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
  func_pardisoinit=(pardisoinit_t)LSL_loadSym(Pardiso_handle, "pardisoinit", msgbuf, msglen);
  if (func_pardisoinit == NULL) return 1;

  func_pardiso=(pardiso_t)LSL_loadSym(Pardiso_handle, "pardiso", msgbuf, msglen);
  if (func_pardiso == NULL) return 1;

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
    fprintf(stderr, "Error loading Pardiso dynamic library " PARDISOLIBNAME ": %s\nAborting...\n", buffer);
    exit(EXIT_FAILURE);
  }
}
