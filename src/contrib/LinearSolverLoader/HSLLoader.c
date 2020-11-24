/* Copyright (C) 2008, 2011 GAMS Development and others
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Author: Stefan Vigerske
 */

#include "IpoptConfig.h"
#include "LibraryHandler.h"
#include "HSLLoader.h"

#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#else
/* use normal C-naming style */
#define IPOPT_HSL_FUNC(name,NAME) name
#endif

#include <stdio.h>
#include <stdlib.h>

#define HSLLIBNAME "libhsl." SHAREDLIBEXT

#ifdef IPOPT_SINGLE
#define HSLFUNCNAMESUFFIX ""
#define HSLFUNCNAMESUFFIXUC ""
#define IPOPT_HSL_FUNCP(name,NAME) IPOPT_HSL_FUNC(name,NAME)
#else
#define HSLFUNCNAMESUFFIX "d"
#define HSLFUNCNAMESUFFIXUC "D"
#define IPOPT_HSL_FUNCP(name,NAME) IPOPT_HSL_FUNC(name ## d,NAME ## D)
#endif

/* make MA27 available in HSL loader if MA27(S) not available in HSL */
#if (defined(IPOPT_SINGLE) && !defined(COINHSL_HAS_MA27S)) || (!defined(IPOPT_SINGLE) && !defined(COINHSL_HAS_MA27))
#define HSLLOADER_MA27
#endif

/* make MA57 available in HSL loader if MA57(S) not available in HSL */
#if (defined(IPOPT_SINGLE) && !defined(COINHSL_HAS_MA57S)) || (!defined(IPOPT_SINGLE) && !defined(COINHSL_HAS_MA57))
#define HSLLOADER_MA57
#endif

/* make MC19 available in HSL loader if MC19(S) not available in HSL */
#if (defined(IPOPT_SINGLE) && !defined(COINHSL_HAS_MC19S)) || (!defined(IPOPT_SINGLE) && !defined(COINHSL_HAS_MC19))
#define HSLLOADER_MC19
#endif

static soHandle_t HSL_handle = NULL;

void LSL_lateHSLLoad(void);

typedef void (*voidfun)(void);

voidfun LSL_loadSym(
   soHandle_t  h,
   const char* symName,
   char*       msgBuf,
   int         msgLen
);

#ifdef HSLLOADER_MA27
static ma27a_t func_ma27a = NULL;
static ma27b_t func_ma27b = NULL;
static ma27c_t func_ma27c = NULL;
static ma27i_t func_ma27i = NULL;

void IPOPT_HSL_FUNCP(ma27a, MA27A)(
   ipfint*       N,
   ipfint*       NZ,
   const ipfint* IRN,
   const ipfint* ICN,
   ipfint*       IW,
   ipfint*       LIW,
   ipfint*       IKEEP,
   ipfint*       IW1,
   ipfint*       NSTEPS,
   ipfint*       IFLAG,
   ipfint*       ICNTL,
   ipnumber*     CNTL,
   ipfint*       INFO,
   ipnumber*     OPS
)
{
   if( func_ma27a == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma27a == NULL )
   {
      fprintf(stderr, "HSL routine MA27A" HSLFUNCNAMESUFFIXUC " not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma27a(N, NZ, IRN, ICN, IW, LIW, IKEEP, IW1, NSTEPS, IFLAG, ICNTL, CNTL, INFO, OPS);
}

void IPOPT_HSL_FUNCP(ma27b, MA27B)(
   ipfint*       N,
   ipfint*       NZ,
   const ipfint* IRN,
   const ipfint* ICN,
   ipnumber*     A,
   ipfint*       LA,
   ipfint*       IW,
   ipfint*       LIW,
   ipfint*       IKEEP,
   ipfint*       NSTEPS,
   ipfint*       MAXFRT,
   ipfint*       IW1,
   ipfint*       ICNTL,
   ipnumber*     CNTL,
   ipfint*       INFO
)
{
   if( func_ma27b == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma27b == NULL )
   {
      fprintf(stderr, "HSL routine MA27B" HSLFUNCNAMESUFFIXUC " not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma27b(N, NZ, IRN, ICN, A, LA, IW, LIW, IKEEP, NSTEPS, MAXFRT, IW1, ICNTL, CNTL, INFO);
}

void IPOPT_HSL_FUNCP(ma27c, MA27C)(
   ipfint*   N,
   ipnumber* A,
   ipfint*   LA,
   ipfint*   IW,
   ipfint*   LIW,
   ipnumber* W,
   ipfint*   MAXFRT,
   ipnumber* RHS,
   ipfint*   IW1,
   ipfint*   NSTEPS,
   ipfint*   ICNTL,
   ipnumber* CNTL
)
{
   if( func_ma27c == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma27c == NULL )
   {
      fprintf(stderr, "HSL routine MA27C" HSLFUNCNAMESUFFIXUC " not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma27c(N, A, LA, IW, LIW, W, MAXFRT, RHS, IW1, NSTEPS, ICNTL, CNTL);
}

void IPOPT_HSL_FUNCP(ma27i, MA27I)(
   ipfint*   ICNTL,
   ipnumber* CNTL
)
{
   if( func_ma27i == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma27i == NULL )
   {
      fprintf(stderr, "HSL routine MA27I" HSLFUNCNAMESUFFIXUC " not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma27i(ICNTL, CNTL);
}
#endif

#ifdef HSLLOADER_MA57
static ma57i_t func_ma57i = NULL;
static ma57a_t func_ma57a = NULL;
static ma57b_t func_ma57b = NULL;
static ma57c_t func_ma57c = NULL;
static ma57e_t func_ma57e = NULL;

void IPOPT_HSL_FUNCP(ma57i, MA57I)(
   ipnumber* cntl,
   ipfint*   icntl
)
{
   if( func_ma57i == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma57i == NULL )
   {
      fprintf(stderr, "HSL routine MA57I" HSLFUNCNAMESUFFIXUC " not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma57i(cntl, icntl);
}

void IPOPT_HSL_FUNCP(ma57a, MA57A)(
   ipfint*       n,       /**< Order of matrix. */
   ipfint*       ne,      /**< Number of entries. */
   const ipfint* irn,     /**< Matrix nonzero row structure */
   const ipfint* jcn,     /**< Matrix nonzero column structure */
   ipfint*       lkeep,   /**< Workspace for the pivot order of lenght 3*n */
   ipfint*       keep,    /**< Workspace for the pivot order of lenght 3*n */
   ipfint*       iwork,   /**< Integer work space. */
   ipfint*       icntl,   /**< Integer Control parameter of length 30*/
   ipfint*       info,    /**< Statistical Information; Integer array of length 20 */
   ipnumber*     rinfo    /**< Double Control parameter of length 5 */
)
{
   if( func_ma57a == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma57a == NULL )
   {
      fprintf(stderr, "HSL routine MA57A" HSLFUNCNAMESUFFIXUC " not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma57a(n, ne, irn, jcn, lkeep, keep, iwork, icntl, info, rinfo);
}

void IPOPT_HSL_FUNCP(ma57b, MA57B)(
   ipfint*    n,         /**< Order of matrix. */
   ipfint*    ne,        /**< Number of entries. */
   ipnumber*  a,         /**< Numerical values. */
   ipnumber*  fact,      /**< Entries of factors. */
   ipfint*    lfact,     /**< Length of array `fact'. */
   ipfint*    ifact,     /**< Indexing info for factors. */
   ipfint*    lifact,    /**< Length of array `ifact'. */
   ipfint*    lkeep,     /**< Length of array `keep'. */
   ipfint*    keep,      /**< Integer array. */
   ipfint*    iwork,     /**< Workspace of length `n'. */
   ipfint*    icntl,     /**< Integer Control parameter of length 20. */
   ipnumber*  cntl,      /**< Double Control parameter of length 5. */
   ipfint*    info,      /**< Statistical Information; Integer array of length 40. */
   ipnumber*  rinfo      /**< Statistical Information; Real array of length 20. */
)
{
   if( func_ma57b == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma57b == NULL )
   {
      fprintf(stderr, "HSL routine MA57B" HSLFUNCNAMESUFFIXUC " not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma57b(n, ne, a, fact, lfact, ifact, lifact, lkeep, keep, iwork, icntl, cntl, info, rinfo);
}

void IPOPT_HSL_FUNCP(ma57c, MA57C)(
   ipfint*    job,       /**< Solution job.  Solve for... */
   ipfint*    n,         /**< Order of matrix. */
   ipnumber*  fact,      /**< Entries of factors. */
   ipfint*    lfact,     /**< Length of array `fact'. */
   ipfint*    ifact,     /**< Indexing info for factors. */
   ipfint*    lifact,    /**< Length of array `ifact'. */
   ipfint*    nrhs,      /**< Number of right hand sides. */
   ipnumber*  rhs,       /**< Numerical Values. */
   ipfint*    lrhs,      /**< Leading dimensions of `rhs'. */
   ipnumber*  work,      /**< Real workspace. */
   ipfint*    lwork,     /**< Length of `work', >= N*NRHS. */
   ipfint*    iwork,     /**< Integer array of length `n'. */
   ipfint*    icntl,     /**< Integer Control parameter array of length 20. */
   ipfint*    info       /**< Statistical Information; Integer array of length 40. */
)
{
   if( func_ma57c == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma57c == NULL )
   {
      fprintf(stderr, "HSL routine MA57C" HSLFUNCNAMESUFFIXUC " not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma57c(job, n, fact, lfact, ifact, lifact, nrhs, rhs, lrhs, work, lwork, iwork, icntl, info);
}

void IPOPT_HSL_FUNCP(ma57e, MA57E)(
   ipfint*    n,
   ipfint*    ic,        /**< 0: copy real array.  >=1:  copy integer array. */
   ipfint*    keep,
   ipnumber*  fact,
   ipfint*    lfact,
   ipnumber*  newfac,
   ipfint*    lnew,
   ipfint*    ifact,
   ipfint*    lifact,
   ipfint*    newifc,
   ipfint*    linew,
   ipfint*    info
)
{
   if( func_ma57e == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma57e == NULL )
   {
      fprintf(stderr, "HSL routine MA57E" HSLFUNCNAMESUFFIXUC " not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma57e(n, ic, keep, fact, lfact, newfac, lnew, ifact, lifact, newifc, linew, info);
}
#endif

#ifndef COINHSL_HAS_MA77

static ma77_default_control_t func_ma77_default_control = NULL;
static ma77_open_nelt_t       func_ma77_open_nelt       = NULL;
static ma77_open_t            func_ma77_open            = NULL;
static ma77_input_vars_t      func_ma77_input_vars      = NULL;
static ma77_input_reals_t     func_ma77_input_reals     = NULL;
static ma77_analyse_t         func_ma77_analyse         = NULL;
static ma77_factor_t          func_ma77_factor          = NULL;
static ma77_factor_solve_t    func_ma77_factor_solve    = NULL;
static ma77_solve_t           func_ma77_solve           = NULL;
static ma77_resid_t           func_ma77_resid           = NULL;
static ma77_scale_t           func_ma77_scale           = NULL;
static ma77_enquire_posdef_t  func_ma77_enquire_posdef  = NULL;
static ma77_enquire_indef_t   func_ma77_enquire_indef   = NULL;
static ma77_alter_t           func_ma77_alter           = NULL;
static ma77_restart_t         func_ma77_restart         = NULL;
static ma77_finalise_t        func_ma77_finalise        = NULL;

/** Initialize control with default values */
void ma77_default_control(
   struct ma77_control_d* control
)
{
   if( func_ma77_default_control == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_default_control == NULL )
   {
      fprintf(stderr, "HSL routine ma77_default_control not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_default_control(control);
}

void ma77_open_nelt(
   const int                    n,
   const char*                  fname1,
   const char*                  fname2,
   const char*                  fname3,
   const char*                  fname4,
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info,
   const int                    nelt
)
{
   if( func_ma77_open_nelt == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_open_nelt == NULL )
   {
      fprintf(stderr, "HSL routine ma77_open_nelt not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_open_nelt(n, fname1, fname2, fname3, fname4, keep, control, info, nelt);
}

void ma77_open(
   const int                    n,
   const char*                  fname1,
   const char*                  fname2,
   const char*                  fname3,
   const char*                  fname4,
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info
)
{
   if( func_ma77_open == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_open == NULL )
   {
      fprintf(stderr, "HSL routine ma77_open not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_open(n, fname1, fname2, fname3, fname4, keep, control, info);
}

void ma77_input_vars(
   const int                    idx,
   const int                    nvar,
   const int                    list[],
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info
)
{
   if( func_ma77_input_vars == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_input_vars == NULL )
   {
      fprintf(stderr, "HSL routine ma77_input_vars not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_input_vars(idx, nvar, list, keep, control, info);
}

void ma77_input_reals(
   const int                    idx,
   const int                    length,
   const ma77pkgtype_d_         reals[],
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info
)
{
   if( func_ma77_input_reals == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_input_reals == NULL )
   {
      fprintf(stderr, "HSL routine ma77_input_reals not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_input_reals(idx, length, reals, keep, control, info);
}

void ma77_analyse(
   const int                    order[],
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info
)
{
   if( func_ma77_analyse == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_analyse == NULL )
   {
      fprintf(stderr, "HSL routine ma77_analyse not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_analyse(order, keep, control, info);
}

void ma77_factor(
   const int                    posdef,
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info,
   const ma77pkgtype_d_*        scale
)
{
   if( func_ma77_factor == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_factor == NULL )
   {
      fprintf(stderr, "HSL routine ma77_factor not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_factor(posdef, keep, control, info, scale);
}

void ma77_factor_solve(
   const int                    posdef,
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info,
   const ma77pkgtype_d_*        scale,
   const int                    nrhs,
   const int                    lx,
   ma77pkgtype_d_               rhs[]
)
{
   if( func_ma77_factor_solve == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_factor_solve == NULL )
   {
      fprintf(stderr, "HSL routine ma77_factor_solve not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_factor_solve(posdef, keep, control, info, scale, nrhs, lx, rhs);
}

void ma77_solve(
   const int                    job,
   const int                    nrhs,
   const int                    lx,
   ma77pkgtype_d_               x[],
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info,
   const ma77pkgtype_d_*        scale
)
{
   if( func_ma77_solve == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_solve == NULL )
   {
      fprintf(stderr, "HSL routine ma77_solve not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_solve(job, nrhs, lx, x, keep, control, info, scale);
}

void ma77_resid(
   const int                    nrhs,
   const int                    lx,
   const ma77pkgtype_d_         x[],
   const int                    lresid,
   ma77pkgtype_d_               resid[],
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info,
   ma77pkgtype_d_*              anorm_bnd
)
{
   if( func_ma77_resid == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_resid == NULL )
   {
      fprintf(stderr, "HSL routine ma77_resid not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_resid(nrhs, lx, x, lresid, resid, keep, control, info, anorm_bnd);
}

void ma77_scale(
   ma77pkgtype_d_               scale[],
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info,
   ma77pkgtype_d_*              anorm
)
{
   if( func_ma77_scale == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_scale == NULL )
   {
      fprintf(stderr, "HSL routine ma77_scale not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_scale(scale, keep, control, info, anorm);
}

void ma77_enquire_posdef(
   ma77pkgtype_d_               d[],
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info
)
{
   if( func_ma77_enquire_posdef == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_enquire_posdef == NULL )
   {
      fprintf(stderr, "HSL routine ma77_enquire_posdef not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_enquire_posdef(d, keep, control, info);
}

void ma77_enquire_indef(
   int                          piv_order[],
   ma77pkgtype_d_               d[],
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info
)
{
   if( func_ma77_enquire_indef == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_enquire_indef == NULL )
   {
      fprintf(stderr, "HSL routine ma77_enquire_indef not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_enquire_indef(piv_order, d, keep, control, info);
}

void ma77_alter(
   const ma77pkgtype_d_         d[],
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info
)
{
   if( func_ma77_alter == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_alter == NULL )
   {
      fprintf(stderr, "HSL routine ma77_alter not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_alter(d, keep, control, info);
}

void ma77_restart(
   const char*                  restart_file,
   const char*                  fname1,
   const char*                  fname2,
   const char*                  fname3,
   const char*                  fname4,
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info
)
{
   if( func_ma77_restart == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_restart == NULL )
   {
      fprintf(stderr, "HSL routine ma77_restart not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_restart(restart_file, fname1, fname2, fname3, fname4, keep, control, info);
}

void ma77_finalise(
   void**                       keep,
   const struct ma77_control_d* control,
   struct ma77_info_d*          info
)
{
   if( func_ma77_finalise == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma77_finalise == NULL )
   {
      fprintf(stderr, "HSL routine ma77_finalise not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma77_finalise(keep, control, info);
}

#endif


#ifndef COINHSL_HAS_MA86

static ma86_default_control_t func_ma86_default_control = NULL;
static ma86_analyse_t         func_ma86_analyse         = NULL;
static ma86_factor_t          func_ma86_factor          = NULL;
static ma86_factor_solve_t    func_ma86_factor_solve    = NULL;
static ma86_solve_t           func_ma86_solve           = NULL;
static ma86_finalise_t        func_ma86_finalise        = NULL;

void ma86_default_control(
   struct ma86_control* control
)
{
   if( func_ma86_default_control == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma86_default_control == NULL )
   {
      fprintf(stderr, "HSL routine ma86_default_control not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma86_default_control(control);
}

void ma86_analyse(
   const int                  n,
   const int                  ptr[],
   const int                  row[],
   int                        order[],
   void**                     keep,
   const struct ma86_control* control,
   struct ma86_info*          info
)
{
   if( func_ma86_analyse == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma86_analyse == NULL )
   {
      fprintf(stderr, "HSL routine ma86_analyse not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma86_analyse(n, ptr, row, order, keep, control, info);
}

void ma86_factor(
   const int                  n,
   const int                  ptr[],
   const int                  row[],
   const ma86pkgtype_d_       val[],
   const int                  order[],
   void**                     keep,
   const struct ma86_control* control,
   struct ma86_info*          info,
   const ma86pkgtype_d_       scale[]
)
{
   if( func_ma86_factor == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma86_factor == NULL )
   {
      fprintf(stderr, "HSL routine ma86_factor not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma86_factor(n, ptr, row, val, order, keep, control, info, scale);
}

void ma86_factor_solve(
   const int                  n,
   const int                  ptr[],
   const int                  row[],
   const ma86pkgtype_d_       val[],
   const int                  order[],
   void**                     keep,
   const struct ma86_control* control,
   struct ma86_info*          info,
   const int                  nrhs,
   const int                  ldx,
   ma86pkgtype_d_             x[],
   const ma86pkgtype_d_       scale[]
)
{
   if( func_ma86_factor_solve == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma86_factor_solve == NULL )
   {
      fprintf(stderr, "HSL routine ma86_factor_solve not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma86_factor_solve(n, ptr, row, val, order, keep, control, info, nrhs,
                          ldx, x, scale);
}

void ma86_solve(
   const int                  job,
   const int                  nrhs,
   const int                  ldx,
   ma86pkgtype_d_*            x,
   const int                  order[],
   void**                     keep,
   const struct ma86_control* control,
   struct ma86_info*          info,
   const ma86pkgtype_d_       scale[]
)
{
   if( func_ma86_solve == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma86_solve == NULL )
   {
      fprintf(stderr, "HSL routine ma86_solve not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma86_solve(job, nrhs, ldx, x, order, keep, control, info, scale);
}

void ma86_finalise(
   void**                     keep,
   const struct ma86_control* control
)
{
   if( func_ma86_finalise == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma86_finalise == NULL )
   {
      fprintf(stderr, "HSL routine ma86_finalise not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma86_finalise(keep, control);
}
#endif

#ifndef COINHSL_HAS_MA97

static ma97_default_control_t func_ma97_default_control = NULL;
static ma97_analyse_t         func_ma97_analyse         = NULL;
static ma97_factor_t          func_ma97_factor          = NULL;
static ma97_factor_solve_t    func_ma97_factor_solve    = NULL;
static ma97_solve_t           func_ma97_solve           = NULL;
static ma97_finalise_t        func_ma97_finalise        = NULL;
static ma97_free_akeep_t      func_ma97_free_akeep      = NULL;

void ma97_default_control(
   struct ma97_control* control
)
{
   if( func_ma97_default_control == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma97_default_control == NULL )
   {
      fprintf(stderr, "HSL routine ma97_default_control not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma97_default_control(control);
}

void ma97_analyse(
   const int                  check,
   const int                  n,
   const int                  ptr[],
   const int                  row[],
   ma97pkgtype_d_             val[],
   void**                     akeep,
   const struct ma97_control* control,
   struct ma97_info*          info,
   int                        order[]
)
{
   if( func_ma97_analyse == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma97_analyse == NULL )
   {
      fprintf(stderr, "HSL routine ma97_analyse not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma97_analyse(check, n, ptr, row, val, akeep, control, info, order);
}

void ma97_factor(
   const int                  matrix_type,
   const int                  ptr[],
   const int                  row[],
   const ma97pkgtype_d_       val[],
   void**                     akeep,
   void**                     fkeep,
   const struct ma97_control* control,
   struct ma97_info*          info,
   const ma97pkgtype_d_       scale[]
)
{
   if( func_ma97_factor == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma97_factor == NULL )
   {
      fprintf(stderr, "HSL routine ma97_factor not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma97_factor(matrix_type, ptr, row, val, akeep, fkeep, control, info,
                    scale);
}

void ma97_factor_solve(
   const int                  matrix_type,
   const int                  ptr[],
   const int                  row[],
   const ma97pkgtype_d_       val[],
   const int                  nrhs,
   ma97pkgtype_d_             x[],
   const int                  ldx,
   void**                     akeep,
   void**                     fkeep,
   const struct ma97_control* control,
   struct ma97_info*          info,
   const ma97pkgtype_d_       scale[]
)
{
   if( func_ma97_factor_solve == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma97_factor_solve == NULL )
   {
      fprintf(stderr, "HSL routine ma97_factor_solve not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma97_factor_solve(matrix_type, ptr, row, val, nrhs, x, ldx, akeep, fkeep, control, info, scale);
}

void ma97_solve(
   const int                  job,
   const int                  nrhs,
   ma97pkgtype_d_*            x,
   const int                  ldx,
   void**                     akeep,
   void**                     fkeep,
   const struct ma97_control* control,
   struct ma97_info*          info
)
{
   if( func_ma97_solve == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma97_solve == NULL )
   {
      fprintf(stderr, "HSL routine ma97_solve not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma97_solve(job, nrhs, x, ldx, akeep, fkeep, control, info);
}

void ma97_finalise(
   void** akeep,
   void** fkeep
)
{
   if( func_ma97_finalise == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma97_finalise == NULL )
   {
      fprintf(stderr, "HSL routine ma97_finalise not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma97_finalise(akeep, fkeep);
}

void ma97_free_akeep(
   void** akeep
)
{
   if( func_ma97_free_akeep == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_ma97_free_akeep == NULL )
   {
      fprintf(stderr, "HSL routine ma97_free_akeep not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_ma97_free_akeep(akeep);
}
#endif

#ifdef HSLLOADER_MC19
static mc19a_t func_mc19a = NULL;

void IPOPT_HSL_FUNCP(mc19a, MC19A)(
   ipfint*   N,
   ipfint*   NZ,
   ipnumber* A,
   ipfint*   IRN,
   ipfint*   ICN,
   float*    R,
   float*    C,
   float*    W
)
{
   if( func_mc19a == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_mc19a == NULL )
   {
      fprintf(stderr, "HSL routine MC19A" HSLFUNCNAMESUFFIXUC " not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_mc19a(N, NZ, A, IRN, ICN, R, C, W);
}
#endif

#ifndef COINHSL_HAS_MC68

mc68_default_control_t func_mc68_default_control = NULL;
mc68_order_t           func_mc68_order           = NULL;

void mc68_default_control_i(
   struct mc68_control_i* control
)
{
   if( func_mc68_default_control == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_mc68_default_control == NULL )
   {
      fprintf(stderr, "HSL routine mc68_default_control not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_mc68_default_control(control);
}

void mc68_order_i(
   int                          ord,
   int                          n,
   const int                    ptr[],
   const int                    row[],
   int                          perm[],
   const struct mc68_control_i* control,
   struct mc68_info_i*          info
)
{
   if( func_mc68_order == NULL )
   {
      LSL_lateHSLLoad();
   }
   if( func_mc68_order == NULL )
   {
      fprintf(stderr, "HSL routine mc68_default_control not found in " HSLLIBNAME ".\nAbort...\n");
      exit(EXIT_FAILURE);
   }
   func_mc68_order(ord, n, ptr, row, perm, control, info);
}

#endif

int LSL_loadHSL(
   const char* libname,
   char*       msgbuf,
   int         msglen
)
{
   /* load HSL library */
   if( libname )
   {
      HSL_handle = LSL_loadLib(libname, msgbuf, msglen);
   }
   else /* try a default library name */
   {
      HSL_handle = LSL_loadLib(HSLLIBNAME, msgbuf, msglen);
   }
   if( HSL_handle == NULL )
   {
      return 1;
   }

   /* load HSL functions */
#ifdef HSLLOADER_MA27
   func_ma27i = (ma27i_t)LSL_loadSym(HSL_handle, "ma27i" HSLFUNCNAMESUFFIX, msgbuf, msglen);
   func_ma27a = (ma27a_t)LSL_loadSym(HSL_handle, "ma27a" HSLFUNCNAMESUFFIX, msgbuf, msglen);
   func_ma27b = (ma27b_t)LSL_loadSym(HSL_handle, "ma27b" HSLFUNCNAMESUFFIX, msgbuf, msglen);
   func_ma27c = (ma27c_t)LSL_loadSym(HSL_handle, "ma27c" HSLFUNCNAMESUFFIX, msgbuf, msglen);
#endif

#ifdef HSLLOADER_MA57
   func_ma57i = (ma57i_t)LSL_loadSym(HSL_handle, "ma57i" HSLFUNCNAMESUFFIX, msgbuf, msglen);
   func_ma57a = (ma57a_t)LSL_loadSym(HSL_handle, "ma57a" HSLFUNCNAMESUFFIX, msgbuf, msglen);
   func_ma57b = (ma57b_t)LSL_loadSym(HSL_handle, "ma57b" HSLFUNCNAMESUFFIX, msgbuf, msglen);
   func_ma57c = (ma57c_t)LSL_loadSym(HSL_handle, "ma57c" HSLFUNCNAMESUFFIX, msgbuf, msglen);
   func_ma57e = (ma57e_t)LSL_loadSym(HSL_handle, "ma57e" HSLFUNCNAMESUFFIX, msgbuf, msglen);
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
   func_ma86_default_control = (ma86_default_control_t)LSL_loadSym(HSL_handle, "ma86_default_control_d", msgbuf, msglen);
   func_ma86_analyse = (ma86_analyse_t)LSL_loadSym(HSL_handle, "ma86_analyse_d", msgbuf, msglen);
   func_ma86_factor = (ma86_factor_t)LSL_loadSym(HSL_handle, "ma86_factor_d", msgbuf, msglen);
   func_ma86_factor_solve = (ma86_factor_solve_t)LSL_loadSym(HSL_handle, "ma86_factor_solve_d", msgbuf, msglen);
   func_ma86_solve = (ma86_solve_t)LSL_loadSym(HSL_handle, "ma86_solve_d", msgbuf, msglen);
   func_ma86_finalise = (ma86_finalise_t)LSL_loadSym(HSL_handle, "ma86_finalise_d", msgbuf, msglen);
#endif

#ifndef COINHSL_HAS_MA97
   func_ma97_default_control = (ma97_default_control_t)LSL_loadSym(HSL_handle, "ma97_default_control_d", msgbuf, msglen);
   func_ma97_analyse = (ma97_analyse_t)LSL_loadSym(HSL_handle, "ma97_analyse_d", msgbuf, msglen);
   func_ma97_factor = (ma97_factor_t)LSL_loadSym(HSL_handle, "ma97_factor_d", msgbuf, msglen);
   func_ma97_factor_solve = (ma97_factor_solve_t)LSL_loadSym(HSL_handle, "ma97_factor_solve_d", msgbuf, msglen);
   func_ma97_solve = (ma97_solve_t)LSL_loadSym(HSL_handle, "ma97_solve_d", msgbuf, msglen);
   func_ma97_finalise = (ma97_finalise_t)LSL_loadSym(HSL_handle, "ma97_finalise_d", msgbuf, msglen);
   func_ma97_free_akeep = (ma97_free_akeep_t)LSL_loadSym(HSL_handle, "ma97_free_akeep_d", msgbuf, msglen);
#endif

#ifdef HSLLOADER_MC19
   func_mc19a = (mc19a_t)LSL_loadSym(HSL_handle, "mc19a" HSLFUNCNAMESUFFIX, msgbuf, msglen);
#endif

#ifndef COINHSL_HAS_MC68
   func_mc68_default_control = (mc68_default_control_t)LSL_loadSym(HSL_handle, "mc68_default_control_i", msgbuf, msglen);
   func_mc68_order = (mc68_order_t)LSL_loadSym(HSL_handle, "mc68_order_i", msgbuf, msglen);
#endif

   return 0;
}

int LSL_unloadHSL(void)
{
   int rc;

   if( HSL_handle == NULL )
   {
      return 0;
   }

   rc = LSL_unloadLib(HSL_handle);
   HSL_handle = NULL;

#ifdef HSLLOADER_MA27
   func_ma27i = NULL;
   func_ma27a = NULL;
   func_ma27b = NULL;
   func_ma27c = NULL;
#endif

#ifdef HSLLOADER_MA57
   func_ma57i = NULL;
   func_ma57a = NULL;
   func_ma57b = NULL;
   func_ma57c = NULL;
   func_ma57e = NULL;
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
   func_ma86_default_control = NULL;
   func_ma86_analyse = NULL;
   func_ma86_factor = NULL;
   func_ma86_factor_solve = NULL;
   func_ma86_solve = NULL;
   func_ma86_finalise = NULL;
#endif

#ifndef COINHSL_HAS_MA97
   func_ma97_default_control = NULL;
   func_ma97_analyse = NULL;
   func_ma97_factor = NULL;
   func_ma97_factor_solve = NULL;
   func_ma97_solve = NULL;
   func_ma97_finalise = NULL;
#endif

#ifdef HSLLOADER_MC19
   func_mc19a = NULL;
#endif

#ifndef COINHSL_HAS_MC68
   func_mc68_default_control = NULL;
   func_mc68_order = NULL;
#endif

   return rc;
}

int LSL_isHSLLoaded(void)
{
   return HSL_handle != NULL;
}

int LSL_isMA27available(void)
{
#ifdef HSLLOADER_MA27
   return func_ma27i != NULL && func_ma27a != NULL && func_ma27b != NULL && func_ma27c != NULL;
#else
   return 1;
#endif
}

int LSL_isMA57available(void)
{
#ifdef HSLLOADER_MA57
   return func_ma57i != NULL && func_ma57a != NULL && func_ma57b != NULL && func_ma57c != NULL && func_ma57e != NULL;
#else
   return 1;
#endif
}

int LSL_isMA77available(void)
{
#ifndef COINHSL_HAS_MA77
   return func_ma77_default_control != NULL && func_ma77_open_nelt != NULL && func_ma77_open != NULL && func_ma77_input_vars != NULL && func_ma77_input_reals != NULL && func_ma77_analyse != NULL && func_ma77_factor != NULL && func_ma77_factor_solve != NULL && func_ma77_solve != NULL && func_ma77_resid != NULL && func_ma77_scale != NULL && func_ma77_enquire_posdef != NULL && func_ma77_enquire_indef != NULL && func_ma77_alter != NULL && func_ma77_restart != NULL && func_ma77_finalise != NULL;
#else
   return 1;
#endif
}

int LSL_isMA86available(void)
{
#ifndef COINHSL_HAS_MA86
   return func_ma86_default_control != NULL && func_ma86_analyse != NULL && func_ma86_factor != NULL && func_ma86_factor_solve != NULL && func_ma86_solve != NULL && func_ma86_finalise != NULL;
#else
   return 1;
#endif
}

int LSL_isMA97available(void)
{
#ifndef COINHSL_HAS_MA97
   return func_ma97_default_control != NULL && func_ma97_analyse != NULL && func_ma97_factor != NULL && func_ma97_factor_solve != NULL && func_ma97_solve != NULL && func_ma97_finalise != NULL && func_ma97_free_akeep != NULL;
#else
   return 1;
#endif
}

int LSL_isMC19available(void)
{
#ifdef HSLLOADER_MC19
   return func_mc19a != NULL;
#else
   return 1;
#endif
}

int LSL_isMC68available(void)
{
#ifndef COINHSL_HAS_MC68
   return func_mc68_default_control != NULL && func_mc68_order != NULL;
#else
   return 1;
#endif
}

void LSL_lateHSLLoad(void)
{
   char buffer[512];
   int rc;

   sprintf(buffer, "Error unknown.");
   rc = LSL_loadHSL(NULL, buffer, 512);
   if( rc != 0 )
   {
      fprintf(stderr, "Error loading HSL dynamic library " HSLLIBNAME ": %s\nThis executable was not compiled with the HSL routine you specified.\nYou need to compile the HSL dynamic library to use deferred loading of the linear solver.\nAbort...\n", buffer);
      exit(EXIT_FAILURE);
   }
}

char* LSL_HSLLibraryName(void)
{
   static char name[] = HSLLIBNAME;
   return name;
}

void LSL_setMA27(
   ma27a_t ma27a,
   ma27b_t ma27b,
   ma27c_t ma27c,
   ma27i_t ma27i
)
{
#ifdef HSLLOADER_MA27
   func_ma27a = ma27a;
   func_ma27b = ma27b;
   func_ma27c = ma27c;
   func_ma27i = ma27i;
#else
   (void) ma27a;
   (void) ma27b;
   (void) ma27c;
   (void) ma27i;
#endif // HSLLOADER_MA27
}

void LSL_setMA57(
   ma57a_t ma57a,
   ma57b_t ma57b,
   ma57c_t ma57c,
   ma57e_t ma57e,
   ma57i_t ma57i
)
{
#ifdef HSLLOADER_MA57
   func_ma57a = ma57a;
   func_ma57b = ma57b;
   func_ma57c = ma57c;
   func_ma57e = ma57e;
   func_ma57i = ma57i;
#else
   (void) ma57a;
   (void) ma57b;
   (void) ma57c;
   (void) ma57e;
   (void) ma57i;
#endif // HSLLOADER_MA57
}

void LSL_setMA77(
   ma77_default_control_t ma77_default_control,
   ma77_open_nelt_t       ma77_open_nelt,
   ma77_open_t            ma77_open,
   ma77_input_vars_t      ma77_input_vars,
   ma77_input_reals_t     ma77_input_reals,
   ma77_analyse_t         ma77_analyse,
   ma77_factor_t          ma77_factor,
   ma77_factor_solve_t    ma77_factor_solve,
   ma77_solve_t           ma77_solve,
   ma77_resid_t           ma77_resid,
   ma77_scale_t           ma77_scale,
   ma77_enquire_posdef_t  ma77_enquire_posdef,
   ma77_enquire_indef_t   ma77_enquire_indef,
   ma77_alter_t           ma77_alter,
   ma77_restart_t         ma77_restart,
   ma77_finalise_t        ma77_finalise
)
{
#ifndef COINHSL_HAS_MA77
   func_ma77_default_control = ma77_default_control;
   func_ma77_open_nelt       = ma77_open_nelt;
   func_ma77_open            = ma77_open;
   func_ma77_input_vars      = ma77_input_vars;
   func_ma77_input_reals     = ma77_input_reals;
   func_ma77_analyse         = ma77_analyse;
   func_ma77_factor          = ma77_factor;
   func_ma77_factor_solve    = ma77_factor_solve;
   func_ma77_solve           = ma77_solve;
   func_ma77_resid           = ma77_resid;
   func_ma77_scale           = ma77_scale;
   func_ma77_enquire_posdef  = ma77_enquire_posdef;
   func_ma77_enquire_indef   = ma77_enquire_indef;
   func_ma77_alter           = ma77_alter;
   func_ma77_restart         = ma77_restart;
   func_ma77_finalise        = ma77_finalise;
#else
   (void) ma77_default_control;
   (void) ma77_open_nelt;
   (void) ma77_open;
   (void) ma77_input_vars;
   (void) ma77_input_reals;
   (void) ma77_analyse;
   (void) ma77_factor;
   (void) ma77_factor_solve;
   (void) ma77_solve;
   (void) ma77_resid;
   (void) ma77_scale;
   (void) ma77_enquire_posdef;
   (void) ma77_enquire_indef;
   (void) ma77_alter;
   (void) ma77_restart;
   (void) ma77_finalise;
#endif
}

void LSL_setMA86(
   ma86_default_control_t ma86_default_control,
   ma86_analyse_t         ma86_analyse,
   ma86_factor_t          ma86_factor,
   ma86_factor_solve_t    ma86_factor_solve,
   ma86_solve_t           ma86_solve,
   ma86_finalise_t        ma86_finalise
)
{
#ifndef COINHSL_HAS_MA86
   func_ma86_default_control = ma86_default_control;
   func_ma86_analyse         = ma86_analyse;
   func_ma86_factor          = ma86_factor;
   func_ma86_factor_solve    = ma86_factor_solve;
   func_ma86_solve           = ma86_solve;
   func_ma86_finalise        = ma86_finalise;
#else
   (void) ma86_default_control;
   (void) ma86_analyse;
   (void) ma86_factor;
   (void) ma86_factor_solve;
   (void) ma86_solve;
   (void) ma86_finalise;
#endif
}

void LSL_setMA97(
   ma97_default_control_t ma97_default_control,
   ma97_analyse_t         ma97_analyse,
   ma97_factor_t          ma97_factor,
   ma97_factor_solve_t    ma97_factor_solve,
   ma97_solve_t           ma97_solve,
   ma97_finalise_t        ma97_finalise,
   ma97_free_akeep_t      ma97_free_akeep
)
{
#ifndef COINHSL_HAS_MA97
   func_ma97_default_control = ma97_default_control;
   func_ma97_analyse         = ma97_analyse;
   func_ma97_factor          = ma97_factor;
   func_ma97_factor_solve    = ma97_factor_solve;
   func_ma97_solve           = ma97_solve;
   func_ma97_finalise        = ma97_finalise;
   func_ma97_free_akeep      = ma97_free_akeep;
#else
   (void) ma97_default_control;
   (void) ma97_analyse;
   (void) ma97_factor;
   (void) ma97_factor_solve;
   (void) ma97_solve;
   (void) ma97_finalise;
   (void) ma97_free_akeep;
#endif
}

void LSL_setMC19(
   mc19a_t mc19a
)
{
#ifdef HSLLOADER_MC19
   func_mc19a = mc19a;
#else
   (void) mc19a;
#endif // HSLLOADER_MC19
}

void LSL_setMC68(
   mc68_default_control_t mc68_default_control,
   mc68_order_t           mc68_order
)
{
#ifndef COINHSL_HAS_MC68
   func_mc68_default_control = mc68_default_control;
   func_mc68_order = mc68_order;
#else
   (void) mc68_default_control;
   (void) mc68_order;
#endif
}
