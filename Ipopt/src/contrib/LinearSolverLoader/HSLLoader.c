/* Copyright (C) 2008, 2011 GAMS Development and others
 All Rights Reserved.
 This code is published under the Eclipse Public License.

 $Id$

 Author: Stefan Vigerske
*/

#include "IpoptConfig.h"
#include "LibraryHandler.h"
#include "HSLLoader.h"

#ifdef COIN_HAS_HSL
#include "CoinHslConfig.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#define HSLLIBNAME "libhsl." SHAREDLIBEXT

static soHandle_t HSL_handle=NULL;

void LSL_lateHSLLoad();

typedef void (*voidfun)(void);

voidfun LSL_loadSym (soHandle_t h, const char *symName, char *msgBuf, int msgLen);

#ifndef COINHSL_HAS_MA27
static ma27ad_t func_ma27ad=NULL;
static ma27bd_t func_ma27bd=NULL;
static ma27cd_t func_ma27cd=NULL;
static ma27id_t func_ma27id=NULL;

void F77_FUNC(ma27ad,MA27AD)(ipfint *N, ipfint *NZ, const ipfint *IRN, const ipfint* ICN,
                             ipfint *IW, ipfint* LIW, ipfint* IKEEP, ipfint *IW1,
                             ipfint* NSTEPS, ipfint* IFLAG, ipfint* ICNTL,
                             double* CNTL, ipfint *INFO, double* OPS) {
  if (func_ma27ad==NULL) LSL_lateHSLLoad();
  if (func_ma27ad==NULL) {
    fprintf(stderr, "HSL routine MA27AD not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma27ad(N, NZ, IRN, ICN, IW, LIW, IKEEP, IW1, NSTEPS, IFLAG, ICNTL, CNTL, INFO, OPS);	
}

void F77_FUNC(ma27bd,MA27BD)(ipfint *N, ipfint *NZ, const ipfint *IRN, const ipfint* ICN,
                             double* A, ipfint* LA, ipfint* IW, ipfint* LIW,
                             ipfint* IKEEP, ipfint* NSTEPS, ipfint* MAXFRT,
                             ipfint* IW1, ipfint* ICNTL, double* CNTL,
                             ipfint* INFO) {
  if (func_ma27bd==NULL) LSL_lateHSLLoad();
  if (func_ma27bd==NULL) {
    fprintf(stderr, "HSL routine MA27BD not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma27bd(N, NZ, IRN, ICN, A, LA, IW, LIW, IKEEP, NSTEPS, MAXFRT, IW1, ICNTL, CNTL, INFO);
}

void F77_FUNC(ma27cd,MA27CD)(ipfint *N, double* A, ipfint* LA, ipfint* IW,
                             ipfint* LIW, double* W, ipfint* MAXFRT,
                             double* RHS, ipfint* IW1, ipfint* NSTEPS,
                             ipfint* ICNTL, double* CNTL) {
  if (func_ma27cd==NULL) LSL_lateHSLLoad();
  if (func_ma27cd==NULL) {
    fprintf(stderr, "HSL routine MA27CD not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma27cd(N, A, LA, IW, LIW, W, MAXFRT, RHS, IW1, NSTEPS, ICNTL, CNTL);
}

void F77_FUNC(ma27id,MA27ID)(ipfint* ICNTL, double* CNTL)
{
  if (func_ma27id==NULL) LSL_lateHSLLoad();
  if (func_ma27id==NULL) {
    fprintf(stderr, "HSL routine MA27ID not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma27id(ICNTL, CNTL);
}

#endif

#ifndef COINHSL_HAS_MA28

static ma28ad_t func_ma28ad=NULL;

void F77_FUNC(ma28ad,MA28AD)(void* nsize, void* nz, void* rw, void* licn, void* iw,
		                         void* lirn, void* iw2, void* pivtol, void* iw3, void* iw4, void* rw2, void* iflag) {
  if (func_ma28ad==NULL) LSL_lateHSLLoad();
  if (func_ma28ad==NULL) {
    fprintf(stderr, "HSL routine MA28AD not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma28ad(nsize, nz, rw, licn, iw, lirn, iw2, pivtol, iw3, iw4, rw2, iflag);
}
#endif

#ifndef COINHSL_HAS_MA57

static ma57id_t func_ma57id=NULL;
static ma57ad_t func_ma57ad=NULL;
static ma57bd_t func_ma57bd=NULL;
static ma57cd_t func_ma57cd=NULL;
static ma57ed_t func_ma57ed=NULL;

void  F77_FUNC (ma57id, MA57ID) (double    *cntl,  ipfint    *icntl) {
  if (func_ma57id==NULL) LSL_lateHSLLoad();
  if (func_ma57id==NULL) {
    fprintf(stderr, "HSL routine MA57ID not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma57id(cntl, icntl);
}

void  F77_FUNC (ma57ad, MA57AD) (
    ipfint    *n,     /* Order of matrix. */
    ipfint    *ne,            /* Number of entries. */
    const ipfint    *irn,       /* Matrix nonzero row structure */
    const ipfint    *jcn,       /* Matrix nonzero column structure */
    ipfint    *lkeep,     /* Workspace for the pivot order of lenght 3*n */
    ipfint    *keep,      /* Workspace for the pivot order of lenght 3*n */
    ipfint    *iwork,     /* Integer work space. */
    ipfint    *icntl,     /* Integer Control parameter of length 30*/
    ipfint    *info,      /* Statistical Information; Integer array of length 20 */
    double    *rinfo)    /* Double Control parameter of length 5 */ {
  if (func_ma57ad==NULL) LSL_lateHSLLoad();
  if (func_ma57ad==NULL) {
    fprintf(stderr, "HSL routine MA57AD not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma57ad(n, ne, irn, jcn, lkeep, keep, iwork, icntl, info, rinfo);
}

void  F77_FUNC (ma57bd, MA57BD) (
    ipfint    *n,     /* Order of matrix. */
    ipfint    *ne,            /* Number of entries. */
    double    *a,     /* Numerical values. */
    double    *fact,      /* Entries of factors. */
    ipfint    *lfact,     /* Length of array `fact'. */
    ipfint    *ifact,     /* Indexing info for factors. */
    ipfint    *lifact,    /* Length of array `ifact'. */
    ipfint    *lkeep,     /* Length of array `keep'. */
    ipfint    *keep,      /* Integer array. */
    ipfint    *iwork,     /* Workspace of length `n'. */
    ipfint    *icntl,     /* Integer Control parameter of length 20. */
    double    *cntl,      /* Double Control parameter of length 5. */
    ipfint    *info,      /* Statistical Information; Integer array of length 40. */
    double    *rinfo)    /* Statistical Information; Real array of length 20. */ {
  if (func_ma57bd==NULL) LSL_lateHSLLoad();
  if (func_ma57bd==NULL) {
    fprintf(stderr, "HSL routine MA57BD not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma57bd(n, ne, a, fact, lfact, ifact, lifact, lkeep, keep, iwork, icntl, cntl, info, rinfo);	
}

void  F77_FUNC (ma57cd, MA57CD) (
    ipfint    *job,       /* Solution job.  Solve for... */
    ipfint    *n,         /* Order of matrix. */
    double    *fact,      /* Entries of factors. */
    ipfint    *lfact,     /* Length of array `fact'. */
    ipfint    *ifact,     /* Indexing info for factors. */
    ipfint    *lifact,    /* Length of array `ifact'. */
    ipfint    *nrhs,      /* Number of right hand sides. */
    double    *rhs,       /* Numerical Values. */
    ipfint    *lrhs,      /* Leading dimensions of `rhs'. */
    double    *work,      /* Real workspace. */
    ipfint    *lwork,     /* Length of `work', >= N*NRHS. */
    ipfint    *iwork,     /* Integer array of length `n'. */
    ipfint    *icntl,     /* Integer Control parameter array of length 20. */
    ipfint    *info)     /* Statistical Information; Integer array of length 40. */ {
  if (func_ma57cd==NULL) LSL_lateHSLLoad();
  if (func_ma57cd==NULL) {
    fprintf(stderr, "HSL routine MA57CD not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma57cd(job, n, fact, lfact, ifact, lifact, nrhs, rhs, lrhs, work, lwork, iwork, icntl, info);    
}

void  F77_FUNC (ma57ed, MA57ED) (
    ipfint    *n,
    ipfint    *ic,        /* 0: copy real array.  >=1:  copy integer array. */
    ipfint    *keep,
    double    *fact,
    ipfint    *lfact,
    double    *newfac,
    ipfint    *lnew,
    ipfint    *ifact,
    ipfint    *lifact,
    ipfint    *newifc,
    ipfint    *linew,
    ipfint    *info) {
  if (func_ma57ed==NULL) LSL_lateHSLLoad();
  if (func_ma57ed==NULL) {
    fprintf(stderr, "HSL routine MA57ED not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma57ed(n, ic, keep, fact, lfact, newfac, lnew, ifact, lifact, newifc, linew, info);
}
#endif

#ifndef COINHSL_HAS_MA77

static ma77_default_control_t func_ma77_default_control = NULL;
static ma77_open_nelt_t func_ma77_open_nelt = NULL;
static ma77_open_t func_ma77_open = NULL;
static ma77_input_vars_t func_ma77_input_vars = NULL;
static ma77_input_reals_t func_ma77_input_reals = NULL;
static ma77_analyse_t func_ma77_analyse = NULL;
static ma77_factor_t func_ma77_factor = NULL;
static ma77_factor_solve_t func_ma77_factor_solve = NULL;
static ma77_solve_t func_ma77_solve = NULL;
static ma77_resid_t func_ma77_resid = NULL;
static ma77_scale_t func_ma77_scale = NULL;
static ma77_enquire_posdef_t func_ma77_enquire_posdef = NULL;
static ma77_enquire_indef_t func_ma77_enquire_indef = NULL;
static ma77_alter_t func_ma77_alter = NULL;
static ma77_restart_t func_ma77_restart = NULL;
static ma77_finalise_t func_ma77_finalise = NULL;

/* Initialise control with default values */
void ma77_default_control(struct ma77_control_d *control) {
  if (func_ma77_default_control==NULL) LSL_lateHSLLoad();
  if (func_ma77_default_control==NULL) {
    fprintf(stderr, "HSL routine ma77_default_control not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_default_control(control);
}

void ma77_open_nelt(const int n, const char* fname1, const char* fname2,
   const char *fname3, const char *fname4, void **keep,
   const struct ma77_control_d *control, struct ma77_info_d *info,
   const int nelt) {
  if (func_ma77_open_nelt==NULL) LSL_lateHSLLoad();
  if (func_ma77_open_nelt==NULL) {
    fprintf(stderr, "HSL routine ma77_open_nelt not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_open_nelt(n, fname1, fname2, fname3, fname4, keep, control, info, nelt);
}

void ma77_open(const int n, const char* fname1, const char* fname2,
   const char *fname3, const char *fname4, void **keep,
   const struct ma77_control_d *control, struct ma77_info_d *info) {
  if (func_ma77_open==NULL) LSL_lateHSLLoad();
  if (func_ma77_open==NULL) {
    fprintf(stderr, "HSL routine ma77_open not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_open(n, fname1, fname2, fname3, fname4, keep, control, info);
}

void ma77_input_vars(const int idx, const int nvar, const int list[],
   void **keep, const struct ma77_control_d *control, struct ma77_info_d *info) {
  if (func_ma77_input_vars==NULL) LSL_lateHSLLoad();
  if (func_ma77_input_vars==NULL) {
    fprintf(stderr, "HSL routine ma77_input_vars not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_input_vars(idx, nvar, list, keep, control, info);
}

void ma77_input_reals(const int idx, const int length,
   const ma77pkgtype_d_ reals[], void **keep, const struct ma77_control_d *control,
   struct ma77_info_d *info) {
  if (func_ma77_input_reals==NULL) LSL_lateHSLLoad();
  if (func_ma77_input_reals==NULL) {
    fprintf(stderr, "HSL routine ma77_input_reals not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_input_reals(idx, length, reals, keep, control, info);
}

void ma77_analyse(const int order[], void **keep,
   const struct ma77_control_d *control, struct ma77_info_d *info) {
  if (func_ma77_analyse==NULL) LSL_lateHSLLoad();
  if (func_ma77_analyse==NULL) {
    fprintf(stderr, "HSL routine ma77_analyse not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_analyse(order, keep, control, info);
}

void ma77_factor(const int posdef, void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info,
   const ma77pkgtype_d_ *scale) {
  if (func_ma77_factor==NULL) LSL_lateHSLLoad();
  if (func_ma77_factor==NULL) {
    fprintf(stderr, "HSL routine ma77_factor not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_factor(posdef, keep, control, info, scale);
}

void ma77_factor_solve(const int posdef, void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info,
   const ma77pkgtype_d_ *scale, const int nrhs, const int lx,
   ma77pkgtype_d_ rhs[]) {
  if (func_ma77_factor_solve==NULL) LSL_lateHSLLoad();
  if (func_ma77_factor_solve==NULL) {
    fprintf(stderr, "HSL routine ma77_factor_solve not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_factor_solve(posdef, keep, control, info, scale, nrhs, lx, rhs);
}

void ma77_solve(const int job, const int nrhs, const int lx, ma77pkgtype_d_ x[],
   void **keep, const struct ma77_control_d *control, struct ma77_info_d *info,
   const ma77pkgtype_d_ *scale) {
  if (func_ma77_solve==NULL) LSL_lateHSLLoad();
  if (func_ma77_solve==NULL) {
    fprintf(stderr, "HSL routine ma77_solve not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_solve(job, nrhs, lx, x, keep, control, info, scale);
}

void ma77_resid(const int nrhs, const int lx, const ma77pkgtype_d_ x[],
   const int lresid, ma77pkgtype_d_ resid[], void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info,
   ma77pkgtype_d_ *anorm_bnd) {
  if (func_ma77_resid==NULL) LSL_lateHSLLoad();
  if (func_ma77_resid==NULL) {
    fprintf(stderr, "HSL routine ma77_resid not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_resid(nrhs, lx, x, lresid, resid, keep, control, info, anorm_bnd);
}

void ma77_scale(ma77pkgtype_d_ scale[], void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info,
   ma77pkgtype_d_ *anorm) {
  if (func_ma77_scale==NULL) LSL_lateHSLLoad();
  if (func_ma77_scale==NULL) {
    fprintf(stderr, "HSL routine ma77_scale not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_scale(scale, keep, control, info, anorm);
}

void ma77_enquire_posdef(ma77pkgtype_d_ d[], void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info) {
  if (func_ma77_enquire_posdef==NULL) LSL_lateHSLLoad();
  if (func_ma77_enquire_posdef==NULL) {
    fprintf(stderr, "HSL routine ma77_enquire_posdef not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_enquire_posdef(d, keep, control, info);
}

void ma77_enquire_indef(int piv_order[], ma77pkgtype_d_ d[], void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info) {
  if (func_ma77_enquire_indef==NULL) LSL_lateHSLLoad();
  if (func_ma77_enquire_indef==NULL) {
    fprintf(stderr, "HSL routine ma77_enquire_indef not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_enquire_indef(piv_order, d, keep, control, info);
}

void ma77_alter(const ma77pkgtype_d_ d[], void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info) {
  if (func_ma77_alter==NULL) LSL_lateHSLLoad();
  if (func_ma77_alter==NULL) {
    fprintf(stderr, "HSL routine ma77_alter not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_alter(d, keep, control, info);
}

void ma77_restart(const char *restart_file, const char *fname1, 
   const char *fname2, const char *fname3, const char *fname4, void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info) {
  if (func_ma77_restart==NULL) LSL_lateHSLLoad();
  if (func_ma77_restart==NULL) {
    fprintf(stderr, "HSL routine ma77_restart not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_restart(restart_file, fname1, fname2, fname3, fname4, keep, control, info);
}

void ma77_finalise(void **keep, const struct ma77_control_d *control,
   struct ma77_info_d *info) {
  if (func_ma77_finalise==NULL) LSL_lateHSLLoad();
  if (func_ma77_finalise==NULL) {
    fprintf(stderr, "HSL routine ma77_finalise not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma77_finalise(keep, control, info);
}

#endif


#ifndef COINHSL_HAS_MA86

static ma86_default_control_t func_ma86_default_control=NULL;
static ma86_analyse_t func_ma86_analyse=NULL;
static ma86_factor_t func_ma86_factor=NULL;
static ma86_factor_solve_t func_ma86_factor_solve=NULL;
static ma86_solve_t func_ma86_solve=NULL;
static ma86_finalise_t func_ma86_finalise=NULL;

void ma86_default_control(struct ma86_control *control) {
  if (func_ma86_default_control==NULL) LSL_lateHSLLoad();
  if (func_ma86_default_control==NULL) {
    fprintf(stderr, "HSL routine ma86_default_control not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma86_default_control(control);
}

void ma86_analyse(const int n, const int ptr[], const int row[], int order[],
      void **keep, const struct ma86_control *control, struct ma86_info *info) {
  if (func_ma86_analyse==NULL) LSL_lateHSLLoad();
  if (func_ma86_analyse==NULL) {
    fprintf(stderr, "HSL routine ma86_analyse not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma86_analyse(n, ptr, row, order, keep, control, info);
}

void ma86_factor(const int n, const int ptr[], const int row[],
      const ma86pkgtype_d_ val[], const int order[], void **keep,
      const struct ma86_control *control, struct ma86_info *info,
      const ma86pkgtype_d_ scale[]) {
  if (func_ma86_factor==NULL) LSL_lateHSLLoad();
  if (func_ma86_factor==NULL) {
    fprintf(stderr, "HSL routine ma86_factor not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma86_factor(n, ptr, row, val, order, keep, control, info, scale);
}

void ma86_factor_solve(const int n, const int ptr[], const int row[],
      const ma86pkgtype_d_ val[], const int order[], void **keep,
      const struct ma86_control *control, struct ma86_info *info,
      const int nrhs, const int ldx, ma86pkgtype_d_ x[],
      const ma86pkgtype_d_ scale[]) {
  if (func_ma86_factor_solve==NULL) LSL_lateHSLLoad();
  if (func_ma86_factor_solve==NULL) {
    fprintf(stderr, "HSL routine ma86_factor_solve not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma86_factor_solve(n, ptr, row, val, order, keep, control, info, nrhs,
      ldx, x, scale);
}

void ma86_solve(const int job, const int nrhs, const int ldx, ma86pkgtype_d_ *x,
      const int order[], void **keep, const struct ma86_control *control,
      struct ma86_info *info, const ma86pkgtype_d_ scale[]) {
  if (func_ma86_solve==NULL) LSL_lateHSLLoad();
  if (func_ma86_solve==NULL) {
    fprintf(stderr, "HSL routine ma86_solve not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma86_solve(job, nrhs, ldx, x, order, keep, control, info, scale);
}

void ma86_finalise(void **keep, const struct ma86_control *control) {
  if (func_ma86_finalise==NULL) LSL_lateHSLLoad();
  if (func_ma86_finalise==NULL) {
    fprintf(stderr, "HSL routine ma86_finalise not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma86_finalise(keep, control);
}
#endif

#ifndef COINHSL_HAS_MA97

static ma97_default_control_t func_ma97_default_control=NULL;
static ma97_analyse_t func_ma97_analyse=NULL;
static ma97_factor_t func_ma97_factor=NULL;
static ma97_factor_solve_t func_ma97_factor_solve=NULL;
static ma97_solve_t func_ma97_solve=NULL;
static ma97_finalise_t func_ma97_finalise=NULL;
static ma97_free_akeep_t func_ma97_free_akeep=NULL;

void ma97_default_control(struct ma97_control *control) {
  if (func_ma97_default_control==NULL) LSL_lateHSLLoad();
  if (func_ma97_default_control==NULL) {
    fprintf(stderr, "HSL routine ma97_default_control not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma97_default_control(control);
}

void ma97_analyse(const int check, const int n, const int ptr[],
      const int row[], ma97pkgtype_d_ val[], void **akeep,
      const struct ma97_control *control, struct ma97_info *info, int order[]) {
  if (func_ma97_analyse==NULL) LSL_lateHSLLoad();
  if (func_ma97_analyse==NULL) {
    fprintf(stderr, "HSL routine ma97_analyse not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma97_analyse(check, n, ptr, row, val, akeep, control, info, order);
}

void ma97_factor(const int matrix_type, const int ptr[], const int row[],
      const ma97pkgtype_d_ val[], void **akeep, void **fkeep,
      const struct ma97_control *control, struct ma97_info *info,
      const ma97pkgtype_d_ scale[]) {
  if (func_ma97_factor==NULL) LSL_lateHSLLoad();
  if (func_ma97_factor==NULL) {
    fprintf(stderr, "HSL routine ma97_factor not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma97_factor(matrix_type, ptr, row, val, akeep, fkeep, control, info,
        scale);
}

void ma97_factor_solve(const int matrix_type, const int ptr[], const int row[],
      const ma97pkgtype_d_ val[], const int nrhs, ma97pkgtype_d_ x[],
      const int ldx, void **akeep, void **fkeep,
      const struct ma97_control *control, struct ma97_info *info,
      const ma97pkgtype_d_ scale[]) {
  if (func_ma97_factor_solve==NULL) LSL_lateHSLLoad();
  if (func_ma97_factor_solve==NULL) {
    fprintf(stderr, "HSL routine ma97_factor_solve not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma97_factor_solve(matrix_type, ptr, row, val, nrhs, x, ldx, akeep, fkeep, control, info, scale);
}

void ma97_solve(const int job, const int nrhs, ma97pkgtype_d_ *x, const int ldx,
      void **akeep, void **fkeep, const struct ma97_control *control,
      struct ma97_info *info) {
  if (func_ma97_solve==NULL) LSL_lateHSLLoad();
  if (func_ma97_solve==NULL) {
    fprintf(stderr, "HSL routine ma97_solve not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma97_solve(job, nrhs, x, ldx, akeep, fkeep, control, info);
}

void ma97_finalise(void **akeep, void **fkeep) {
  if (func_ma97_finalise==NULL) LSL_lateHSLLoad();
  if (func_ma97_finalise==NULL) {
    fprintf(stderr, "HSL routine ma97_finalise not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma97_finalise(akeep, fkeep);
}

void ma97_free_akeep(void **akeep) {
  if (func_ma97_free_akeep==NULL) LSL_lateHSLLoad();
  if (func_ma97_free_akeep==NULL) {
    fprintf(stderr, "HSL routine ma97_free_akeep not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_ma97_free_akeep(akeep);
}
#endif

#ifndef COINHSL_HAS_MC19

static mc19ad_t func_mc19ad=NULL;

void F77_FUNC(mc19ad,MC19AD)(ipfint *N, ipfint *NZ, double* A, ipfint *IRN,
                             ipfint* ICN, float* R, float* C, float* W) {
  if (func_mc19ad==NULL) LSL_lateHSLLoad();
  if (func_mc19ad==NULL) {
    fprintf(stderr, "HSL routine MC19AD not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_mc19ad(N, NZ, A, IRN, ICN, R, C, W);
}
#endif

#ifndef COINHSL_HAS_MC68

mc68_default_control_t func_mc68_default_control=NULL;
mc68_order_t func_mc68_order=NULL;

void mc68_default_control_i(struct mc68_control_i *control) {
  if (func_mc68_default_control==NULL) LSL_lateHSLLoad();
  if (func_mc68_default_control==NULL) {
    fprintf(stderr, "HSL routine mc68_default_control not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_mc68_default_control(control);
}

void mc68_order_i(int ord, int n, const int ptr[], const int row[],
                  int perm[], const struct mc68_control_i *control,
                  struct mc68_info_i *info) {
  if (func_mc68_order==NULL) LSL_lateHSLLoad();
  if (func_mc68_order==NULL) {
    fprintf(stderr, "HSL routine mc68_default_control not found in " HSLLIBNAME ".\nAbort...\n");
    exit(EXIT_FAILURE);
  }
  func_mc68_order(ord, n, ptr, row, perm, control, info);
}

#endif

int LSL_loadHSL(const char* libname, char* msgbuf, int msglen) {
  /* load HSL library */
  if (libname) {
    HSL_handle=LSL_loadLib(libname, msgbuf, msglen);
  } else { /* try a default library name */
    HSL_handle=LSL_loadLib(HSLLIBNAME, msgbuf, msglen); 
  }
  if (HSL_handle==NULL)
    return 1;
	
  /* load HSL functions */
#ifndef COINHSL_HAS_MA27
  func_ma27id=(ma27id_t)LSL_loadSym(HSL_handle, "ma27id", msgbuf, msglen);
  func_ma27ad=(ma27ad_t)LSL_loadSym(HSL_handle, "ma27ad", msgbuf, msglen);
  func_ma27bd=(ma27bd_t)LSL_loadSym(HSL_handle, "ma27bd", msgbuf, msglen);
  func_ma27cd=(ma27cd_t)LSL_loadSym(HSL_handle, "ma27cd", msgbuf, msglen);
#endif

#ifndef COINHSL_HAS_MA28
  func_ma28ad=(ma28ad_t)LSL_loadSym(HSL_handle, "ma28ad", msgbuf, msglen);
#endif

#ifndef COINHSL_HAS_MA57
  func_ma57id=(ma57id_t)LSL_loadSym(HSL_handle, "ma57id", msgbuf, msglen);
  func_ma57ad=(ma57ad_t)LSL_loadSym(HSL_handle, "ma57ad", msgbuf, msglen);
  func_ma57bd=(ma57bd_t)LSL_loadSym(HSL_handle, "ma57bd", msgbuf, msglen);
  func_ma57cd=(ma57cd_t)LSL_loadSym(HSL_handle, "ma57cd", msgbuf, msglen);
  func_ma57ed=(ma57ed_t)LSL_loadSym(HSL_handle, "ma57ed", msgbuf, msglen);
#endif

#ifndef COINHSL_HAS_MA77
#define LOADMA77SYM( sym ) \
  func_##sym = (sym##_t)LSL_loadSym(HSL_handle, #sym"_d", msgbuf, msglen);

  LOADMA77SYM(ma77_default_control)
  LOADMA77SYM(ma77_open_nelt)
  LOADMA77SYM(ma77_open)
  LOADMA77SYM(ma77_input_vars)
  LOADMA77SYM(ma77_input_reals)
  LOADMA77SYM(ma77_analyse)
  LOADMA77SYM(ma77_factor)
  LOADMA77SYM(ma77_factor_solve)
  LOADMA77SYM(ma77_solve)
  LOADMA77SYM(ma77_resid)
  LOADMA77SYM(ma77_scale)
  LOADMA77SYM(ma77_enquire_posdef)
  LOADMA77SYM(ma77_enquire_indef)
  LOADMA77SYM(ma77_alter)
  LOADMA77SYM(ma77_restart)
  LOADMA77SYM(ma77_finalise)
#endif

#ifndef COINHSL_HAS_MA86
  func_ma86_default_control=(ma86_default_control_t)LSL_loadSym(HSL_handle, "ma86_default_control_d", msgbuf, msglen);
  func_ma86_analyse=(ma86_analyse_t)LSL_loadSym(HSL_handle, "ma86_analyse_d", msgbuf, msglen);
  func_ma86_factor=(ma86_factor_t)LSL_loadSym(HSL_handle, "ma86_factor_d", msgbuf, msglen);
  func_ma86_factor_solve=(ma86_factor_solve_t)LSL_loadSym(HSL_handle, "ma86_factor_solve_d", msgbuf, msglen);
  func_ma86_solve=(ma86_solve_t)LSL_loadSym(HSL_handle, "ma86_solve_d", msgbuf, msglen);
  func_ma86_finalise=(ma86_finalise_t)LSL_loadSym(HSL_handle, "ma86_finalise_d", msgbuf, msglen);
#endif

#ifndef COINHSL_HAS_MA97
  func_ma97_default_control=(ma97_default_control_t)LSL_loadSym(HSL_handle, "ma97_default_control_d", msgbuf, msglen);
  func_ma97_analyse=(ma97_analyse_t)LSL_loadSym(HSL_handle, "ma97_analyse_d", msgbuf, msglen);
  func_ma97_factor=(ma97_factor_t)LSL_loadSym(HSL_handle, "ma97_factor_d", msgbuf, msglen);
  func_ma97_factor_solve=(ma97_factor_solve_t)LSL_loadSym(HSL_handle, "ma97_factor_solve_d", msgbuf, msglen);
  func_ma97_solve=(ma97_solve_t)LSL_loadSym(HSL_handle, "ma97_solve_d", msgbuf, msglen);
  func_ma97_finalise=(ma97_finalise_t)LSL_loadSym(HSL_handle, "ma97_finalise_d", msgbuf, msglen);
  func_ma97_free_akeep=(ma97_free_akeep_t)LSL_loadSym(HSL_handle, "ma97_free_akeep_d", msgbuf, msglen);
#endif

#ifndef COINHSL_HAS_MC19
  func_mc19ad=(mc19ad_t)LSL_loadSym(HSL_handle, "mc19ad", msgbuf, msglen);
#endif

#ifndef COINHSL_HAS_MC68
  func_mc68_default_control=(mc68_default_control_t)LSL_loadSym(HSL_handle, "mc68_default_control_i", msgbuf, msglen);
  func_mc68_order=(mc68_order_t)LSL_loadSym(HSL_handle, "mc68_order_i", msgbuf, msglen);
#endif

  return 0;
}

int LSL_unloadHSL() {
  int rc;
	
  if (HSL_handle==NULL)
    return 0;

  rc = LSL_unloadLib(HSL_handle);
  HSL_handle=NULL;
	
#ifndef COINHSL_HAS_MA27
  func_ma27id=NULL;
  func_ma27ad=NULL;
  func_ma27bd=NULL;
  func_ma27cd=NULL;
#endif

#ifndef COINHSL_HAS_MA28
  func_ma28ad=NULL;
#endif

#ifndef COINHSL_HAS_MA57
  func_ma57id=NULL;
  func_ma57ad=NULL;
  func_ma57bd=NULL;
  func_ma57cd=NULL;
  func_ma57ed=NULL;
#endif

#ifndef COINHSL_HAS_MA77
  func_ma77_default_control = NULL;
  func_ma77_open_nelt = NULL;
  func_ma77_open = NULL;
  func_ma77_input_vars = NULL;
  func_ma77_input_reals = NULL;
  func_ma77_analyse = NULL;
  func_ma77_factor = NULL;
  func_ma77_factor_solve = NULL;
  func_ma77_solve = NULL;
  func_ma77_resid = NULL;
  func_ma77_scale = NULL;
  func_ma77_enquire_posdef = NULL;
  func_ma77_enquire_indef = NULL;
  func_ma77_alter = NULL;
  func_ma77_restart = NULL;
  func_ma77_finalise = NULL;
#endif
  
#ifndef COINHSL_HAS_MA86
  func_ma86_default_control=NULL;
  func_ma86_analyse=NULL;
  func_ma86_factor=NULL;
  func_ma86_factor_solve=NULL;
  func_ma86_solve=NULL;
  func_ma86_finalise=NULL;
#endif

#ifndef COINHSL_HAS_MA97
  func_ma97_default_control=NULL;
  func_ma97_analyse=NULL;
  func_ma97_factor=NULL;
  func_ma97_factor_solve=NULL;
  func_ma97_solve=NULL;
  func_ma97_finalise=NULL;
#endif

#ifndef COINHSL_HAS_MC19
  func_mc19ad=NULL;
#endif

#ifndef COINHSL_HAS_MC68
  func_mc68_default_control=NULL;
  func_mc68_order=NULL;
#endif

  return rc;
}

int LSL_isHSLLoaded() {
  return HSL_handle!=NULL;
}

int LSL_isMA27available() {
#ifndef COINHSL_HAS_MA27
	return func_ma27id!=NULL && func_ma27ad!=NULL && func_ma27bd!=NULL && func_ma27cd!=NULL;
#else
	return 1;
#endif
}

int LSL_isMA28available() {
#ifndef COINHSL_HAS_MA28
	return func_ma28ad!=NULL;
#else
	return 1;
#endif
}

int LSL_isMA57available() {
#ifndef COINHSL_HAS_MA57
	return func_ma57id!=NULL && func_ma57ad!=NULL && func_ma57bd!=NULL && func_ma57cd!=NULL && func_ma57ed!=NULL;
#else
	return 1;
#endif
}

int LSL_isMA77available() {
#ifndef COINHSL_HAS_MA77
	return func_ma77_default_control!=NULL && func_ma77_open_nelt!=NULL && func_ma77_open!=NULL && func_ma77_input_vars!=NULL && func_ma77_input_reals!=NULL && func_ma77_analyse!=NULL && func_ma77_factor!=NULL && func_ma77_factor_solve!=NULL && func_ma77_solve!=NULL && func_ma77_resid!=NULL && func_ma77_scale!=NULL && func_ma77_enquire_posdef!=NULL && func_ma77_enquire_indef!=NULL && func_ma77_alter!=NULL && func_ma77_restart!=NULL && func_ma77_finalise!=NULL;
#else
	return 1;
#endif
}

int LSL_isMA86available() {
#ifndef COINHSL_HAS_MA86
	return func_ma86_default_control!=NULL && func_ma86_analyse!=NULL && func_ma86_factor!=NULL && func_ma86_factor_solve!=NULL && func_ma86_solve!=NULL && func_ma86_finalise!=NULL;
#else
	return 1;
#endif
}

int LSL_isMA97available() {
#ifndef COINHSL_HAS_MA97
	return func_ma97_default_control!=NULL && func_ma97_analyse!=NULL && func_ma97_factor!=NULL && func_ma97_factor_solve!=NULL && func_ma97_solve!=NULL && func_ma97_finalise!=NULL && func_ma97_free_akeep!=NULL;
#else
	return 1;
#endif
}

int LSL_isMC19available() {
#ifndef COINHSL_HAS_MC19
	return func_mc19ad!=NULL;
#else
	return 1;
#endif
}

int LSL_isMC68available() {
#ifndef COINHSL_HAS_MC68
	return func_mc68_default_control!=NULL && func_mc68_order!=NULL;
#else
	return 1;
#endif
}

void LSL_lateHSLLoad() {
  char buffer[512];
  int rc;

  sprintf(buffer, "Error unknown.");
  rc = LSL_loadHSL(NULL, buffer, 512);
  if (rc!=0) {
    fprintf(stderr, "Error loading HSL dynamic library " HSLLIBNAME ": %s\nThis executable was not compiled with the HSL routine you specified.\nYou need to compile the HSL dynamic library to use deferred loading of the linear solver.\nAbort...\n", buffer);
    exit(EXIT_FAILURE);
  }
}

char* LSL_HSLLibraryName() {
  static char name[] = HSLLIBNAME;
  return name;
}

void LSL_setMA27(ma27ad_t ma27ad, ma27bd_t ma27bd, ma27cd_t ma27cd, ma27id_t ma27id) {
#ifndef COINHSL_HAS_MA27
   func_ma27ad = ma27ad;
   func_ma27bd = ma27bd;
   func_ma27cd = ma27cd;
   func_ma27id = ma27id;
#endif
}

void LSL_setMA28(ma28ad_t ma28ad) {
#ifndef COINHSL_HAS_MA28
   func_ma28ad = ma28ad;
#endif
}

void LSL_setMA57(ma57ad_t ma57ad, ma57bd_t ma57bd, ma57cd_t ma57cd, ma57ed_t ma57ed, ma57id_t ma57id) {
#ifndef COINHSL_HAS_MA57
   func_ma57ad = ma57ad;
   func_ma57bd = ma57bd;
   func_ma57cd = ma57cd;
   func_ma57ed = ma57ed;
   func_ma57id = ma57id;
#endif
}

void LSL_setMA77(ma77_default_control_t ma77_default_control,
   ma77_open_nelt_t ma77_open_nelt,
   ma77_open_t ma77_open,
   ma77_input_vars_t ma77_input_vars,
   ma77_input_reals_t ma77_input_reals,
   ma77_analyse_t ma77_analyse,
   ma77_factor_t ma77_factor,
   ma77_factor_solve_t ma77_factor_solve,
   ma77_solve_t ma77_solve,
   ma77_resid_t ma77_resid,
   ma77_scale_t ma77_scale,
   ma77_enquire_posdef_t ma77_enquire_posdef,
   ma77_enquire_indef_t ma77_enquire_indef,
   ma77_alter_t ma77_alter,
   ma77_restart_t ma77_restart,
   ma77_finalise_t ma77_finalise) {
#ifndef COINHSL_HAS_MA77
   func_ma77_open_nelt = ma77_open_nelt;
   func_ma77_open = ma77_open;
   func_ma77_input_vars = ma77_input_vars;
   func_ma77_input_reals = ma77_input_reals;
   func_ma77_analyse = ma77_analyse;
   func_ma77_factor = ma77_factor;
   func_ma77_factor_solve = ma77_factor_solve;
   func_ma77_solve = ma77_solve;
   func_ma77_resid = ma77_resid;
   func_ma77_scale = ma77_scale;
   func_ma77_enquire_posdef = ma77_enquire_posdef;
   func_ma77_enquire_indef = ma77_enquire_indef;
   func_ma77_alter = ma77_alter;
   func_ma77_restart = ma77_restart;
   func_ma77_finalise = ma77_finalise;
#endif
}

void LSL_setMA86(ma86_default_control_t ma86_default_control,
   ma86_analyse_t ma86_analyse,
   ma86_factor_t ma86_factor,
   ma86_factor_solve_t ma86_factor_solve,
   ma86_solve_t ma86_solve,
   ma86_finalise_t ma86_finalise) {
#ifndef COINHSL_HAS_MA86
   func_ma86_default_control = ma86_default_control;
   func_ma86_analyse = ma86_analyse;
   func_ma86_factor = ma86_factor;
   func_ma86_factor_solve = ma86_factor_solve;
   func_ma86_solve = ma86_solve;
   func_ma86_finalise = ma86_finalise;
#endif
}

void LSL_setMA97(ma97_default_control_t ma97_default_control,
   ma97_analyse_t ma97_analyse,
   ma97_factor_t ma97_factor,
   ma97_factor_solve_t ma97_factor_solve,
   ma97_solve_t ma97_solve,
   ma97_finalise_t ma97_finalise,
   ma97_free_akeep_t ma97_free_akeep) {
#ifndef COINHSL_HAS_MA97
   func_ma97_default_control = ma97_default_control;
   func_ma97_analyse = ma97_analyse;
   func_ma97_factor = ma97_factor;
   func_ma97_factor_solve = ma97_factor_solve;
   func_ma97_solve = ma97_solve;
   func_ma97_finalise = ma97_finalise;
   func_ma97_free_akeep = ma97_free_akeep;
#endif
}

void LSL_setMC19(mc19ad_t mc19ad) {
#ifndef COINHSL_HAS_MC19
   func_mc19ad = mc19ad;
#endif
}

void LSL_setMC68(mc68_default_control_t mc68_default_control, mc68_order_t mc68_order) {
#ifndef COINHSL_HAS_MC68
   func_mc68_default_control = mc68_default_control;
   func_mc68_order = mc68_order;
#endif
}
