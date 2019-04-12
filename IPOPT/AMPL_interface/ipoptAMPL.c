/*
  Copyright (C) 2002, Carnegie Mellon University and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/
/*****************************************************************************
 $Id: ipoptAMPL.c 674 2005-06-30 20:54:16Z andreasw $
     AMPL interface ipopt
*****************************************************************************/
/****************************************************************************
Authors: 05/01/02  Arvind Raghunathan, Andreas Waechter
                   Release as Version IPOPT 2.0.1

 We thank Marcelo Marazzi and Sven Leyffer for providing example routines.
****************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#include <config.h>

#include <stdio.h>
#if STDC_HEADERS
# include <stdlib.h>
# include <stddef.h>
#else
# if HAVE_STDLIB_H
#  include <stdlib.h>
# endif
#endif
#if HAVE_STRING_H
# if !STDC_HEADERS && HAVE_MEMORY_H
#  include <memory.h>
# endif
# include <string.h>
#endif
#if HAVE_STRINGS_H
# include <strings.h>
#endif
#if HAVE_INTTYPES_H
# include <inttypes.h>
#else
# if HAVE_STDINT_H
#  include <stdint.h>
# endif
#endif
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

#include "asl.h"
#include "getstub.h"
#define asl cur_ASL

#define FTYPES_DEFINED
#include <Ipopt.h>

/* Problem specific variables required in subroutine */
fint numEq, numIneq ;
#ifdef INCLUDE_CC
fint numCcon, numVineq ;
#endif

/* Store the Lower and Upper bounds for constraints and variables */
real *conL, *conU ;

/* Variables used in inequality handling */
fint *ineqnos ;
real *slacks ;

/* Variables used in Complementarity handling */
#ifdef INCLUDE_CC
fint *cvar1nos, *cvar2nos ; /* holds indices of variables bearing
			       complementarity relation */
fint *vineqvar ; /* hold the variable number in the complementarity relation
		    that is unbounded - handle for variational inequalities. */

real *cvar1LBnds , *cvar1UBnds , *cvar2LBnds , *cvar2UBnds ; /* Store the
								bounds */
real *cslacks ; /* slacks introduced for each of non-negative variable that is
                   introduced for the unbounded complementarity variable */
real *rslacks ; /* slacks introduced for each of the complementarity constraint
                   that IPOPT sees. */
#endif

/* Variables used thruout the file */
static fint nerror = 0 ;
static fint nobj = 0 ;
static fint nnzH = 0 ; /* number of nonzeros in Hessian */

static fint intzero = 0, one = 1 ;
static real realzero = 0 ;
static real realminusone  = -1. ;

#define nullintopt -2682238
#define nullrealopt -1.278463e202

enum evalhcall {EVINIT, EVH, EVLAG, EVOBJ, EVCON};
static int lasthcall = EVINIT;

/* Constants pertaining MAXREAL etc. */

static real maxREAL = 1e+20 ; /* largest real number recognised in the RHS */
static real minREAL = -1e+20 ; /* smallest real number recognised in the RHS */
#ifdef INCLUDE_CC
static real myInf = 1e+30 ; /* representation of infinity used in bounds
			       representation - innocuous*/
static real myNegInf = -1e+30 ; /* --------------"------------ */
#endif

/* IPOPT default options */

/* Define the following macro, if only the most important
   IPOPT options should be available */

#define FEWOPTIONS

static fint
  opIFILE = nullintopt,
  opIFULL = nullintopt,
  opIPRINT = nullintopt,
  opIOUTPUT = nullintopt,
  opIMAXITER = nullintopt,
  opIMAXCPUSEC = nullintopt,
  opIMUINIT = nullintopt,
  opIQUASI = nullintopt,
  opIREFINEITER = nullintopt,
  opIRESTO = nullintopt,
  opISCALE = nullintopt,
  opISCALERR = nullintopt,
  opIMERIT = nullintopt,
  opISOC = nullintopt,
#ifdef USE_MALLOC
  opIMAX_IW = 0,
#else
# ifdef LENGTH_IW
  opIMAX_IW = LENGTH_IW,
# else
  opIMAX_IW = 1000000,
# endif
#endif
#ifdef USE_MALLOC
  opIMAX_RW = 0,
#else
# ifdef LENGTH_RW
  opIMAX_RW = LENGTH_RW,
# else
  opIMAX_RW = 1000000,
# endif
#endif
#ifndef FEWOPTIONS
  opICG = nullintopt,
  opIINITB = nullintopt,
  opILMLEN = nullintopt,
  opICORRECT = nullintopt,
  opIDAMP = nullintopt,
  opIALPHA = nullintopt,
  opIITERBLOCKMAX = nullintopt,
  opISELBAS = nullintopt,
  opISYMSOLV = nullintopt,
  opIHESSVECT = nullintopt,
  opITRON2DERIVS = nullintopt,
  opITRONHESS = nullintopt,
  opIKKTSCALE = nullintopt,
#endif
#ifdef INCLUDE_CC
  opIMPEC_CCAGG = 0,
  opIMPEC_TRIGGER = 0,
#endif
  opICNRM = nullintopt;

static real
  opDFSCALE = nullrealopt,
  opDINFMAXTOL = nullrealopt,
  opDMU0 = nullrealopt,
  opDBNDFRAC = nullrealopt,
  opDBNDPUSH = nullrealopt,
  opDPIVTOL = nullrealopt,
  opDPIVTOLMAX = nullrealopt,
  opDTRONCGTOL = nullrealopt,
  opDSCALECUT = nullrealopt,
  opDMOVEBOUNDS = nullrealopt,
  opDLAMINITMAX = nullrealopt,
  opDFILLINFACT = nullrealopt,
#ifndef FEWOPTIONS
  opDNUMIN = nullrealopt,
  opDMULIN = nullrealopt,
  opDMUSUPER = nullrealopt,
  opDWATCHTOL = nullrealopt,
  opDSR1TOL = nullrealopt,
  opDSKIPFACT = nullrealopt,
  opDRHO = nullrealopt,
  opDMAXCOND = nullrealopt,
  opDCGTOL = nullrealopt,
  opDRESTOKKTRED = nullrealopt,
  opDS_F = nullrealopt,
  opDS_THETA = nullrealopt,
  opDGAMMA_F = nullrealopt,
  opDGAMMA_THETA = nullrealopt,
  opDTAU = nullrealopt,
  opDVCORRECTFACT = nullrealopt,
  opDERRSUPER = nullrealopt,
  opDMUERRFAC = nullrealopt,
  opDMAXERR = nullrealopt,
  opDLS_SAFE = nullrealopt,
#endif
#ifdef INCLUDE_CC
  opDMPEC_CCFACT = 1,
  opDMPEC_CCPOW = 1,
  opDMPEC_ETAFACT = 0.1,
  opDMPEC_THRESH = 5e-6,
#endif
  opDCMAXTOL = nullrealopt,
  opDDEPCONDIAG = nullrealopt,
  opDTOL = nullrealopt;

static keyword keywds[] =
{
  KW("dbndfrac"     ,D_val,&opDBNDFRAC,
     "Relative distance of starting point from closest bound"),
  KW("dbndpush"     ,D_val,&opDBNDPUSH,
     "Distance of starting point from closest bound"),
#ifndef FEWOPTIONS
  KW("dcgtol"       ,D_val,&opDCGTOL,
     "Tolerance for conjugate gradient option"),
#endif
  KW("dcmaxtol"      ,D_val,&opDCMAXTOL,
     "Tolerance for unscaled primal infeasibility"),
  KW("ddepcondiag"  ,D_val,&opDDEPCONDIAG,
     "Factor in relaxation of linearly dependent constraints (0: off)"),
#ifndef FEWOPTIONS
  KW("derrsuper"    ,D_val,&opDERRSUPER,
     "Exponent for superlinear decrease of barrier error tolerance"),
#endif
  KW("dfillinfact"  ,D_val,&opDFILLINFACT,
     "Factor for estimating memory requirement"),
  KW("dfscale"      ,D_val,&opDFSCALE,
     "Internal scaling factor for objective function"),
#ifndef FEWOPTIONS
  KW("dgamma_f"     ,D_val,&opDGAMMA_F,
     "Factor in filter margin (obj. func.)"),
  KW("dgamma_theta" ,D_val,&opDGAMMA_THETA,
     "Factor in filter margin (constraints)"),
#endif
  KW("dinfmaxtol"   ,D_val,&opDINFMAXTOL,
     "Tolerance for unscaled dual infeasibility"),
  KW("dlaminitmax"  ,D_val,&opDLAMINITMAX,
     "Upper threshold for initial equality constraint multiplier estimate"),
#ifndef FEWOPTIONS
  KW("dls_safe"     ,D_val,&opDLS_SAFE,
     "Error tolerance under which merit function penalty parameter is\n             only increased"),
  KW("dmaxcond"     ,D_val,&opDMAXCOND,
     "When to repartition"),
  KW("dmaxerr"      ,D_val,&opDMAXERR,
     "Upper bound for barrier error tolerance"),
#endif
  KW("dmovebounds"  ,D_val,&opDMOVEBOUNDS,
     "Initial relative perturbation of bounds"),
#ifdef INCLUDE_CC
  KW("dmpec_ccfact"     ,D_val,&opDMPEC_CCFACT,
     "Coefficient for barrier parameter in the relaxation of complementarity constraints"),
  KW("dmpec_ccpow"      ,D_val,&opDMPEC_CCPOW,
     "Exponent of barrier parameter in the relaxation of complementarity constraint"),
  KW("dmpec_etafact"    ,D_val,&opDMPEC_ETAFACT,
     "Factor multiplying modification"),
  KW("dmpec_thresh"  ,D_val,&opDMPEC_THRESH,
     "Threshold mu below which modification is enforced"),
#endif
  KW("dmu0"         ,D_val,&opDMU0,
     "Initial value of barrier parameter"),
#ifndef FEWOPTIONS
  KW("dmuerrfac"    ,D_val,&opDMUERRFAC,
     "Factor between MU and barrier error tolerance"),
  KW("dmulin"       ,D_val,&opDMULIN,
     "Factor for linear decrease of MU"),
  KW("dmusuper"     ,D_val,&opDMUSUPER,
     "Exponent for superlinear decrease of MU"),
  KW("dnumin"       ,D_val,&opDNUMIN,
     "Minimal (=Starting) value of merit function penalty parameter"),
#endif
  KW("dpivtol"      ,D_val,&opDPIVTOL,
     "Initial pivot tolerance in linear equation solver"),
  KW("dpivtolmax"   ,D_val,&opDPIVTOLMAX,
     "Maximal pivot tolerance in linear equation solver"),
#ifndef FEWOPTIONS
  KW("drestokktred" ,D_val,&opDRESTOKKTRED,
     "factor for sufficient reduction in KKT system based reestoration phase"),
  KW("drho"         ,D_val,&opDRHO,
     "Parameter in update rule for penalty parameter"),
  KW("ds_f"         ,D_val,&opDS_F,
     "Exponent in filter switching rule (obj. func.)"),
  KW("ds_theta"     ,D_val,&opDS_THETA,
     "Exponent in filter switching rule (constraints)"),
#endif
  KW("dscalecut"    ,D_val,&opDSCALECUT,
     "Parameter for automatic scaling option (iscale)"),
#ifndef FEWOPTIONS
  KW("dskipfact"    ,D_val,&opDSKIPFACT,
     "Parameter for BFGS update skip"),
  KW("dsr1tol"      ,D_val,&opDSR1TOL,
     "Parameter for switching from BFGS to SR1"),
  KW("dtau"         ,D_val,&opDTAU,
     "Factor in fraction to the boundary rule"),
#endif
  KW("dtol"         ,D_val,&opDTOL,
     "Overall error tolerance (termination criterion)"),
#ifndef FEWOPTIONS
  KW("dtroncgtol"   ,D_val,&opDTRONCGTOL,
     "Tolerance for CG with TRON"),
  KW("dvcorrectfact",D_val,&opDVCORRECTFACT,
     "Factor for maximum deviation of V from S^{-1}MU"),
  KW("dwatchtol"    ,D_val,&opDWATCHTOL,
     "When to activate watchdog technique"),
  KW("ialpha"       ,L_val,&opIALPHA,
     "Treatment of primal and dual fraction to the boundary rule"),
#endif
#ifndef FEWOPTIONS
  KW("icg"          ,L_val,&opICG,
     "CG option and preconditioner choice"),
#endif
  KW("icnrm"        ,L_val,&opICNRM,
     "Type of norm to measure constraint violation during line search"),
#ifndef FEWOPTIONS
  KW("icorrect"     ,L_val,&opICORRECT,
     "Correction strategy for hessian if indefinite"),
  KW("idamp"        ,L_val,&opIDAMP,
     "Option of damping the cross term"),
#endif
  KW("ifile"         ,L_val,&opIFILE,
     "Set to 1 to create output and debugging file IPOPT.OUT"),
  KW("ifull"        ,L_val,&opIFULL,
     "Set to 0 for reduced space option"),
#ifndef FEWOPTIONS
  KW("ihessvect"    ,L_val,&opIHESSVECT,
     "Choice of computing hessian-vector products"),
  KW("iinitb"       ,L_val,&opIINITB,
     "Initialization of reduced hessian"),
  KW("iiterblockmax",L_val,&opIITERBLOCKMAX,
     "When to delete filter"),
  KW("ikktscale"    ,L_val,&opIKKTSCALE,
     "How to scale the KKT matrix before factorization"),
  KW("ilmlen"       ,L_val,&opILMLEN,
     "number of (s,y) pairs for limited memory BFGS"),
#endif
#ifndef USE_MALLOC
  KW("imax_iw"      ,L_val,&opIMAX_IW,
     "Size of integer workspace allotted to IPOPT"),
  KW("imax_rw"      ,L_val,&opIMAX_RW,
     "Size of double presision workspace allotted to IPOPT"),
#endif
  KW("imaxcpusec"   ,L_val,&opIMAXCPUSEC,
     "Maximum number of CPU seconds allowed for computation"),
  KW("imaxiter"     ,L_val,&opIMAXITER,
     "Maximum number of iterations"),
  KW("imerit"       ,L_val,&opIMERIT,
     "Type of line search"),
#ifdef INCLUDE_CC
  KW("impec_ccagg" ,L_val,&opIMPEC_CCAGG,
     "Option of treating complementarity constraints individually/as aggregated"),
  KW("impec_trigger",L_val,&opIMPEC_TRIGGER,
     "Option of using modification for step"),
#endif
  KW("imuinit"      ,L_val,&opIMUINIT,
     "Determines how barrier parameter and bound multipliers are initialized"),
  KW("ioutput"       ,L_val,&opIOUTPUT,
     "Set to 1 for wider and more detailed output to screen"),
  KW("iprint"       ,L_val,&opIPRINT,
     "Log level (-1 for no output, 10 for most output)"),
  KW("iquasi"       ,L_val,&opIQUASI,
     "Activates Quasi-Newton approximation for Hessian (6 for L-BFGS)"),
  KW("irefineiter"  ,L_val,&opIREFINEITER,
     "Minimal number of iterative refinement steps per linear system"),
  KW("iresto"       ,L_val,&opIRESTO,
     "Choice of filter restoration phase"),
  KW("iscale"       ,L_val,&opISCALE,
     "Determines automatic scaling method of problem functions"),
  KW("iscalerr"       ,L_val,&opISCALERR,
     "Set to 0 to use unscaled error estimates"),
#ifndef FEWOPTIONS
  KW("iselbas"      ,L_val,&opISELBAS,
     "How to select basis partition"),
#endif
  KW("isoc"         ,L_val,&opISOC,
     "Maximal number of second order correction per iteration"),
#ifndef FEWOPTIONS
  KW("isymsolv"     ,L_val,&opISYMSOLV,
     "Solver for symetric system (0: MA47; 1: MA27; 2: MA57)"),
  KW("itron2derivs" ,L_val,&opITRON2DERIVS,
     "Choice of Gauss-Newton approximation in TRON"),
  KW("itronhess"    ,L_val,&opITRONHESS,
     "Determines how Hessian-vectors products are computed in TRON"),
#endif
  KW("outlev"       ,L_val,&opIPRINT,
     "Log level (-1 for no output, 10 for most output) - same as iprint"),
} ;

static Option_Info Oinfo = {PACKAGE_NAME,PACKAGE_STRING,"ipopt_options",keywds,nkeywds} ;

/* IPOPT subroutine */
Cextern real F77_FUNC(ipopt,IPOPT)
     (fint *N, real *X, fint *M, fint *NLB, fint *ILB,
      real *BNDS_L, fint *NUB, fint *IUB, real *BNDS_U,
      real *V_L, real *V_U, real *C, real *LAM, fint *LRW,
      real *RW, fint *LIW, fint *IW, fint *ITER, fint *IERR,
      void *EV_F(), void *EV_C(), void *EV_G(), void *EV_A(),
      void *EV_H(), void *EV_HLV(), void *EV_HOV(), void *EV_HCV(),
      real *DAT, fint *IDAT, fint *NARGS, real *ARGS,
      char *CARGS, int lCARGS) ;

/* Routine for obtaining current value of MU */
#ifdef INCLUDE_CC
Cextern void F77_FUNC_(get_amplmu,GET_AMPLMU)(real *MU);
#endif

/* Some auxiliary routines for IPOPT */
Cextern void F77_FUNC_(ipopt_getdata,IPOPT_GETDATA)
     (fint *FFUNCS, fint *CFUNCS, fint *NOCG, fint *NORES, fint *NONEGCURV);
Cextern void F77_FUNC_(check_flagfile,CHECK_FLAGFILE)(fint *EX);
Cextern void F77_FUNC(timer,TIMER)(real *T);

/* BLAS routines */
Cextern void F77_FUNC(dcopy,DCOPY)
     (fint *Nvar, real *X, fint *incX, real *Y, fint *incY) ;
Cextern void F77_FUNC(dscal,DSCAL)
     (fint *Nvar, real *alpha, real *X,fint *incX) ;
Cextern fint F77_FUNC(idamax,IDAMAX)
     (fint *M, real* X, fint *incX);

/* Standard error message if memory allocation not successful */
#define MALLOCERR(Varname) {fprintf(Stderr,"Cannot allocate memory for " #Varname ".\n"); exit(-1);}

/* function for copying C-string into Fortran character array */
void str2fstr(char *Cstr, char *Fstr, int Flen)
{
  int i, len;

  len = strlen(Cstr);
  if( len > Flen) len = Flen;
  memcpy(Fstr, Cstr, len);
  for ( i = len; i < Flen; Fstr[i++] = ' ' );
}


/* Initialising the slack variables introduced for inequalities */
void initslacks(real *X, real *slacks, fint *ineqnos)
{
  fint i ;

  xknown(X); /* make sure, that ASL doesn't check for each call of conival
		if poitn X has changed... */
  for ( i = 0 ; i < numIneq ; i++ )
    {
      slacks[i] = conival(ineqnos[i]-1,X,&nerror) ;
      /* Error Handling */
      if ( nerror < 0 )
	{
	  fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
	  fprintf(Stderr," nerror = %d\n",nerror) ;
	  fprintf(Stderr,"Exiting from subroutine - initslacks\n") ;
	  exit(1) ;
	}
    }
  xunknown();
}

/* Initialising the variables introduced in variational type constraints to
   convert them into complementarity constratints.
   Implemented a very arbitrary initialisation. Need to look into it later */
#ifdef INCLUDE_CC
void initcompslacks(real *X, real *cslacks, real *slacks, fint * ineqnos)
{
  fint i ;

  for ( i = 0 ; i < (2*numVineq) ; i++ )
    cslacks[i] = 1e-1 ;
}
#endif

/* Initialising the slack variable introduced in each of the relaxed
   complementarity constraint */
#ifdef INCLUDE_CC
void initrslacks(real *X, real *rslacks)
{
  fint i ;

  if ( opIMPEC_CCAGG == 0 )
    for ( i = 0 ; i < (numCcon + 2*numVineq) ; i++ )
      rslacks[i] = 1e-1 ;
  else
    *rslacks = 1e-1 ;
}
#endif

/* Objective function evaluation - objfun */
void eval_f(fint *N, real *X, real *F, real *DAT, fint *IDAT)
{

  /* Check to deal with the no objective AMPL problem. */
  if ( n_obj == 0 )
    {
      *F = 0 ;
    }
  else
    {
      *F = objval(nobj, X, &nerror) ;
      /* Check if this is a maximization problem */
      if (objtype[nobj])
	*F = -*F ;
    }

  /* Error handling */
  if ( nerror < 0 )
    {
      fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
      fprintf(Stderr," nerror = %d\n",nerror) ;
      fprintf(Stderr,"Exiting from subroutine - eval_f_\n") ;
      exit(1) ;
    }
  else if ( nerror > 0 )
    {
      *F = exp(1000000.);
      return;
    }

}

/* Objective gradient evaluation - objgrd */
void eval_g(fint *N, real *X, real *G, real *DAT, fint *IDAT)
{
  fint n_varF , totsize ;
  real max2min = -1 ;

  /* Check to deal with the no objective AMPL problem. */
  if ( n_obj == 0 )
    {
      if ( n_cc == 0 )
	totsize = n_var + numIneq ;
#ifdef INCLUDE_CC
      else if ( opIMPEC_CCAGG == 0 )
	totsize = n_var + numIneq + 2*numVineq + numCcon + 2*numVineq ;
      else if ( opIMPEC_CCAGG == 1 )
	totsize = n_var + numIneq + 2*numVineq + 1 ;
#endif
      F77_FUNC(dcopy,DCOPY)(&totsize,&realzero,&intzero,&G[0],&one) ;
    }
  else
    {
      objgrd(nobj, X, G, &nerror) ;

      /* Check for maximization problems */
      if (objtype[nobj])
	{
	  n_varF = n_var ;
	  F77_FUNC(dscal,DSCAL)(&n_varF,&max2min,&G[0],&one) ;
	}

      /* Error handling */
      if ( nerror < 0 )
	{
	  fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
  fprintf(Stderr," nerror = %d\n",nerror) ;
	  fprintf(Stderr,"Exiting from subroutine - eval_g_\n") ;
	  exit(1) ;
	}
      else if( nerror > 0 )
	{
	  G[0] = exp(1000000.);
	  return;
	}
      if ( numIneq > 0 )
	F77_FUNC(dcopy,DCOPY)(&numIneq,&realzero,&intzero,&G[n_var],&one) ;
#ifdef INCLUDE_CC
      if ( numVineq > 0 )
	{
	  totsize = 2*numVineq ;
	  F77_FUNC(dcopy,DCOPY)
	    (&totsize,&realzero,&intzero,&G[n_var+numIneq],&one) ;
	}
      if ( n_cc > 0 )
	if ( opIMPEC_CCAGG == 0 )
	  {
	    totsize = numCcon + 2*numVineq ;
	    F77_FUNC(dcopy,DCOPY)
	      (&totsize,&realzero,&intzero,&G[n_var+numIneq+2*numVineq],&one) ;
	  }
	else
	  G[n_var + numIneq + 2*numVineq] = 0 ;
#endif
    }

}

/* Constraint value evaluation - conval */
void eval_c(fint *N, real *X, fint *M, real *C, real *DAT, fint *IDAT)
{
  fint i,j ;

#ifdef INCLUDE_CC
  fint presize ;
  real mu, sum ;
#endif

  conval(X, C, &nerror) ;

  /* Error handling */
  if ( nerror < 0 )
    {
      fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
      fprintf(Stderr," nerror = %d\n",nerror) ;
      fprintf(Stderr,"Exiting from subroutine - eval_c_\n") ;
      exit(1) ;
    }
  else if( nerror > 0 )
    {
      C[0] = exp(1000000.);
      return;
    }

  /* Evaluate the constraints given by AMPL */

  j = 0 ; /* inequality nos index */

  for ( i = 0 ; i < n_con ; i++ )
    {
      if ( ( j < numIneq ) && ( i == (ineqnos[j]-1) ) )
	  {
	    C[i] = C[i] - X[n_var+j] ;
	    j++ ;
	  }
      else
	C[i] = C[i] - conL[i] ;
    }

  /* Evaluate the complementarity constraints and append them appropriately */
#ifdef INCLUDE_CC

  /* First add constraints for splitting the unbdd. variables */

  for (i = 0 ; i < numVineq ; i++)
    C[n_con+i] = X[vineqvar[i]-1] - (X[n_var+numIneq+2*i] - X[n_var+numIneq+2*i+1]) ;

  /* Now add the complementarity constraints */

  presize = n_var + numIneq + 2*numVineq ; /* number of variables without
					      taking into account the number
					      of relaxed compl. slacks */
  sum = 0 ; /* sum of the complementarity constraints for the aggregated case
	     */

  F77_FUNC_(get_amplmu,GET_AMPLMU)(&mu);

  for ( i = 0 ; i < (numCcon + 2*numVineq) ; i++ )
    {
      /* Look for different cases and write the complementarity relation */
      
      /* l1 <= s <= Infinity complements l2 <= x <= Infinity */
      if ( (cvar1UBnds[i] > maxREAL) && (cvar2UBnds[i] > maxREAL) )
	{
	  if ( opIMPEC_CCAGG == 0 )
	    C[n_con+numVineq+i] = (X[cvar1nos[i]-1] - cvar1LBnds[i])*(X[cvar2nos[i]-1] - cvar2LBnds[i]) + X[presize+i] - opDMPEC_CCFACT*pow(mu,opDMPEC_CCPOW) ;
	  else
	    sum = sum + (X[cvar1nos[i]-1] - cvar1LBnds[i])*(X[cvar2nos[i]-1] - cvar2LBnds[i]) ;
	}
      /* l1 <= s <= Infinity complements -Infinity <= x <= u2 */
      else if ( (cvar1UBnds[i] > maxREAL) && (cvar2LBnds[i] < minREAL) )
	{
	  if ( opIMPEC_CCAGG == 0 )
	    C[n_con+numVineq+i] = (X[cvar1nos[i]-1] - cvar1LBnds[i])*(cvar2UBnds[i] - X[cvar2nos[i]-1]) + X[presize+i] - opDMPEC_CCFACT*pow(mu,opDMPEC_CCPOW) ;
	  else
	    sum = sum + (X[cvar1nos[i]-1] - cvar1LBnds[i])*(cvar2UBnds[i] - X[cvar2nos[i]-1]) ;
	}
      /* -Infinity <= s <= u1 complements l1 <= x <= Infinity */
      else if ( (cvar1LBnds[i] < minREAL) && (cvar2UBnds[i] > maxREAL) )
	{
	  if ( opIMPEC_CCAGG == 0 )
	    C[n_con+numVineq+i] = (cvar1UBnds[i] - X[cvar1nos[i]-1])*(X[cvar2nos[i]-1] - cvar2LBnds[i]) + X[presize+i] - opDMPEC_CCFACT*pow(mu,opDMPEC_CCPOW) ;
	  else
	    sum = sum + (cvar1UBnds[i] - X[cvar1nos[i]-1])*(X[cvar2nos[i]-1] - cvar2LBnds[i]) ;
	}
      /* -Infinity <= s <= u1 complements -Infinity <= x <= u2 */
      else if ( (cvar1LBnds[i] < minREAL) && (cvar2LBnds[i] < minREAL) )
	{
	  if ( opIMPEC_CCAGG == 0 )
	    C[n_con+numVineq+i] = (cvar1UBnds[i] - X[cvar1nos[i]-1])*(cvar2UBnds[i] - X[cvar2nos[i]-1]) + X[presize+i] - opDMPEC_CCFACT*pow(mu,opDMPEC_CCPOW) ;
	  else
	    sum = sum + (cvar1UBnds[i] - X[cvar1nos[i]-1])*(cvar2UBnds[i] - X[cvar2nos[i]-1]) ;
	}
    }

  if (n_cc > 0 && opIMPEC_CCAGG == 1)
    C[n_con+2*numVineq] = sum + X[presize] - opDMPEC_CCFACT*pow(mu,opDMPEC_CCPOW) ;
#endif

  /* Need to decide if this is really neceesary: */
  lasthcall = EVINIT;

}

/* Constraint gradients evaluation - jacval */
void eval_a(fint *TASK, fint *N, real *X, fint *NZ, real *A, fint *ACON,
	    fint *AVAR, real *DAT, fint *IDAT)
{
  real *internalA ;
  cgrad *cg ;
  fint i , j = 0 ;
#ifdef INCLUDE_CC
  fint presize, m  ;
#endif

  if ( *TASK == 0 )
    {
      if ( n_cc == 0 )
	*NZ = nzc + numIneq ;
#ifdef INCLUDE_CC
      else if ( opIMPEC_CCAGG == 0 )
	*NZ = nzc + numIneq + 2*numCcon + 7*numVineq + numCcon + 2*numVineq ;
      else if ( opIMPEC_CCAGG == 1 )
	*NZ = nzc + numIneq + 2*numCcon + 7*numVineq + 1 ;
#endif
    }
  else if ( *TASK == 1 )
    {
      if( !(internalA = (real *)Malloc(nzc*sizeof(real))) ) MALLOCERR(internalA) ;
      
      jacval(X, internalA, &nerror) ;
      
      /* Error handling */
      if ( nerror < 0 )
	{
	  fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
	  fprintf(Stderr," nerror = %d\n",nerror) ;
	  fprintf(Stderr,"Exiting from subroutine - eval_a_\n") ;
	  exit(1) ;
	}
      else if( nerror > 0 )
	{
	  internalA[0] = exp(1000000.);
	  return;
	}

      for ( i = 0 ; i < n_con ; i++ )
	for ( cg = Cgrad[i] ; cg ; cg = cg->next )
	  {
	    A[j] = internalA[cg->goff] ;
	    ACON[j] = i + 1 ;
	    AVAR[j] = cg->varno + 1 ;
	    j++ ;
	  }
      /* Jacobian of inequalities */
      for ( i = 0 ; i < numIneq ; i++ )
	{
	  A[j] = -1 ;
	  ACON[j] = ineqnos[i] ;
	  AVAR[j] = n_var+i+1 ;
	  j++ ;
	}

      /* Jacobian of complementarity constraints */
#ifdef INCLUDE_CC

      /* First Jacobian of splits of unbndd constraints/variables */

      for (i = 0 ; i < numVineq ; i++)
	{
	  A[j] = 1 ;
	  ACON[j] = n_con + i + 1 ;
	  AVAR[j] = vineqvar[i] ;
	  j++ ;
	  
	  A[j] = -1 ;
	  ACON[j] = n_con + i + 1 ;
	  AVAR[j] = n_var + numIneq + 2*i + 1 ;
	  j++ ;
	  
	  A[j] = 1 ;
	  ACON[j] = n_con + i + 1 ;
	  AVAR[j] = n_var + numIneq + 2*i + 2 ;
	  j++ ; 
	}

      /* Now add Jacobian of the compl. constraints */

      presize = n_var + numIneq + 2*numVineq ; /* size before adding relaxed
						  compl. slacks */

      for ( i = 0 ; i < (numCcon + 2*numVineq) ; i++ )
	{
	  /* Look for different cases and write the complementarity relation */

	  if ( opIMPEC_CCAGG == 0 )
	    m = n_con + numVineq + i + 1 ;
	  else
	    m = n_con + numVineq + 1 ;

	  /* l1 <= s <= Infinity complements l2 <= x <= Infinity */
	  if ( (cvar1UBnds[i] > maxREAL) && (cvar2UBnds[i] > maxREAL) )
	    {
	      A[j] = X[cvar2nos[i]-1] - cvar2LBnds[i] ;
	      ACON[j] = m ;
	      AVAR[j] = cvar1nos[i] ;
  	      j++ ;

	      A[j] = X[cvar1nos[i]-1] - cvar1LBnds[i] ;
	      ACON[j] = m ;
	      AVAR[j] = cvar2nos[i] ;
	      j++ ;
	    }
	  /* l1 <= s <= Infinity complements -Infinity <= x <= u2 */
	  else if ( (cvar1UBnds[i] > maxREAL) && (cvar2LBnds[i] < minREAL) )
	    {
	      A[j] = cvar2UBnds[i] - X[cvar2nos[i]-1] ;
	      ACON[j] = m ;
	      AVAR[j] = cvar1nos[i] ;
	      j++ ;

	      A[j] = -(X[cvar1nos[i]-1] - cvar1LBnds[i]) ;
	      ACON[j] = m ;
	      AVAR[j] = cvar2nos[i] ;
	      j++ ;
	    }
	  /* -Infinity <= s <= u1 complements l1 <= x <= Infinity */
	  else if ( (cvar1LBnds[i] < minREAL) && (cvar2UBnds[i] > maxREAL) )
	    {
	      A[j] = -(X[cvar2nos[i]-1] - cvar2LBnds[i]) ;
	      ACON[j] = m ;
	      AVAR[j] = cvar1nos[i] ;
	      j++ ;

	      A[j] = cvar1UBnds[i] - X[cvar1nos[i]-1] ;
	      ACON[j] = m ;
	      AVAR[j] = cvar2nos[i] ;
	      j++ ;
	    }
	  /* -Infinity <= s <= u1 complements -Infinity <= x <= u2 */
	  else if ( (cvar1LBnds[i] < minREAL) && (cvar2LBnds[i] < minREAL) )
	    {
	      A[j] = -(cvar2UBnds[i] - X[cvar2nos[i]-1]) ;
	      ACON[j] = m ;
	      AVAR[j] = cvar1nos[i] ;	
	      j++ ;

	      A[j] = -(cvar1UBnds[i] - X[cvar1nos[i]-1]) ;
	      ACON[j] = m ;
	      AVAR[j] = cvar2nos[i] ;	
	      j++ ;
	    }
	  if ( opIMPEC_CCAGG == 0 )
	    {
	      A[j] = 1 ;
	      ACON[j] = m ;
	      AVAR[j] = presize + i + 1 ;
	      j++ ;
	    }
	}
      if ( n_cc > 0 && opIMPEC_CCAGG == 1 )
	{
	  A[j] = 1 ;
	  ACON[j] = n_con + numVineq + 1 ;
	  AVAR[j] = presize + 1 ;
	}
#endif

      free(internalA) ;
    }

}

/* Hessian of Lagrangian evaluation */
void eval_h(fint *TASK, fint *N, real *X, fint *M, real *LAM, fint *NNZH,
	    real *HESS, fint *IRNH, fint *ICNH, real *DAT, fint *IDAT)
{
  fint ow , y , i , j , k , uptri ;
  real OW , *H , temp , *cdummy ;

  if ( *TASK == 0 )
    {
      /* Check to deal with the no objective AMPL problem. */
      if ( n_obj == 0 )
	{
	  ow = 0 ; /* coefficient of objective function in the Lagrangian */
	  y = 1 ; /* mulitpliers will be supplied */
	  uptri = 1 ; /* need only the upper triangular part */
	
	  nnzH = sphsetup(nobj,ow,y,uptri) ;
	}
      else
	{
	  ow = 1 ; /* coefficient of objective function in the Lagrangian */
	  y = 1 ; /* multipliers will be supplied */
	  uptri = 1 ; /* need only the upper triangular part */
	
	  nnzH = sphsetup(-1,ow,y,uptri) ;
	}
      *NNZH = nnzH ;
      /* The Hessain does not return the complementarity derivatives. The
	 one below is necessary as AMPL provides complementarity relation
	 between a variable and a constraint which translates to a AMPL
	 variable and an internal slack - so no duplication would occur. */
#ifdef INCLUDE_CC
      if (n_cc > 0)
	*NNZH = *NNZH + numCcon + 2*numVineq ;
#endif
    }
  else if ( *TASK == 1 )
    {

      if( !(H = (real *)Malloc(nnzH*sizeof(real))) ) MALLOCERR(H) ;

      if( 0 ) {
      /* Innocuous call to make sure we are evaluating Hessian at the given
         point */
      if (n_obj != 0)
	temp = objval(nobj,X,&nerror) ;

      /* Error handling */
      if ( nerror < 0 )
	{
	  fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
	  fprintf(Stderr," nerror = %d\n",nerror) ;
	  fprintf(Stderr,"Exiting from subroutine - eval_h_\n") ;
	  exit(1) ;
	}

      /* Calls to evaluate constraint values */

      if( !(cdummy = (real *)Malloc(nzc*sizeof(real))) ) MALLOCERR(cdummy) ;
      conval(X, cdummy, &nerror) ;
      /* Error handling */
      if ( nerror < 0 )
	{
	  fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
	  fprintf(Stderr," nerror = %d\n",nerror) ;
	  fprintf(Stderr,"Exiting from subroutine - eval_h_\n") ;
	  exit(1) ;
	}
      free(cdummy);
      }

      if (n_obj > 0)
	if (objtype[nobj])
	  OW = -1 ;
	else
	  OW = 1 ;

      if ( n_obj > 0 )
	{
	  sphes(H,-1,&OW,LAM) ;
	}
      else
	{
	  sphes(H,nobj,NULL,LAM) ;
	}
      lasthcall = EVH;

      k = 0 ;
      for ( i = 0 ; i < n_var ; i++ )
	for ( j = sputinfo->hcolstarts[i] ; j < sputinfo->hcolstarts[i+1] ; j++ )
	  {
	    HESS[k] = H[j] ;
	    ICNH[k] = i + 1 ;
	    IRNH[k] = sputinfo->hrownos[j] + 1 ;
	    k++ ;
	  }

#ifdef INCLUDE_CC

      /* Add Hessian terms due to complementarity constraints */ 

      for ( i = 0 ; i < (numCcon + 2*numVineq) ; i++ )
	{
	  /* First check for type of bounds on the complementary variables */
	  if ( (cvar1LBnds[i] < minREAL && cvar2LBnds[i] < minREAL) || (cvar1UBnds[i] > maxREAL && cvar2UBnds[i] > maxREAL) )
	    temp = 1 ;
	  else
	    temp = -1 ;
	  
	  if ( opIMPEC_CCAGG == 0 )
	    HESS[k] = temp*LAM[n_con+numVineq+i] ;
	  else
	    HESS[k] = temp*LAM[n_con+numVineq] ;

	  if (cvar1nos[i] > cvar2nos[i])
	    {
	      ICNH[k] = cvar2nos[i] ;
	      IRNH[k] = cvar1nos[i] ;
	    }
	  else
	    {
	      ICNH[k] = cvar1nos[i] ;
	      IRNH[k] = cvar2nos[i] ;
	    }
	  k++ ;
	}
#endif
      free(H) ;
    }

}

/**************Computation of Hessian-Vector products****************/

/***************Hessian of Lagrangian times vector********************/

void eval_hesslag_v(fint *TASK, fint *N, real *X, fint *M, real *LAM,
		    real *VIN, real *VOUT, real *DAT, fint *IDAT)
{
  real OW , temp, *cdummy ;
#ifdef INCLUDE_CC
  fint i, totsize ;
#endif

  /* Dummy call to ensure evaluation of Hessian at X */
  if (n_obj > 0)
    temp = objval(nobj, X, &nerror) ;

  /* Error handling */
  if ( nerror < 0 )
    {
      fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
      fprintf(Stderr," nerror = %d\n",nerror) ;
      fprintf(Stderr,"Exiting from subroutine - eval_hesslag_v_\n") ;
      exit(1) ;
    }

  /* Calls to evaluate constraint values */

  if( !(cdummy = (real *)Malloc(nzc*sizeof(real))) ) MALLOCERR(cdummy) ;
  conval(X, cdummy, &nerror) ;
  /* Error handling */
  if ( nerror < 0 )
    {
      fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
      fprintf(Stderr," nerror = %d\n",nerror) ;
      fprintf(Stderr,"Exiting from subroutine - eval_hesslag_v_\n") ;
      exit(1) ;
    }
  free(cdummy);

  if (n_obj > 0)
    if (objtype[nobj])
      OW = -1 ;
    else
      OW = 1 ;

  if( lasthcall == EVINIT || lasthcall != EVLAG )
    hvinit(-1,&OW,LAM) ;
  hvcomp(VOUT,VIN,-1,&OW,LAM) ;
  lasthcall = EVLAG;

  /* Copy zeros into the remaining part of the vector */
  if ( numIneq > 0 )
    F77_FUNC(dcopy,DCOPY)(&numIneq,&realzero,&intzero,&VOUT[n_var],&one) ;

#ifdef INCLUDE_CC

  /* Add terms due to complementarity constraints */

  if ( numVineq > 0 )
    {
      totsize = 2*numVineq ;
      F77_FUNC(dcopy,DCOPY)
	(&totsize,&realzero,&intzero,&VOUT[n_var+numIneq],&one) ;
    }

  if ( n_cc > 0 )
    if ( opIMPEC_CCAGG == 0 )
      {
	totsize = numCcon + 2*numVineq ;
	F77_FUNC(dcopy,DCOPY)
	  (&totsize,&realzero,&intzero,&VOUT[n_var+numIneq+2*numVineq],&one) ;
      }
    else
      VOUT[n_var+numIneq+2*numVineq] = 0 ;
  
  for ( i = 0 ; i < (numCcon + 2*numVineq) ; i++ )
    {
      /* First check for type of bounds on the complementary variables */
      if ( (cvar1LBnds[i] < minREAL && cvar2LBnds[i] < minREAL) || (cvar1UBnds[i] > maxREAL && cvar2UBnds[i] > maxREAL) )
	temp = 1 ;
      else
	temp = -1 ;

      if ( opIMPEC_CCAGG == 0 )
	{	
	  VOUT[cvar1nos[i]-1] += temp*LAM[n_con+numVineq+i]*VIN[cvar2nos[i]-1] ;
	  VOUT[cvar2nos[i]-1] += temp*LAM[n_con+numVineq+i]*VIN[cvar1nos[i]-1] ;
	}
      else
	{
	  VOUT[cvar1nos[i]-1] += temp*LAM[n_con+numVineq]*VIN[cvar2nos[i]-1] ;
	  VOUT[cvar2nos[i]-1] += temp*LAM[n_con+numVineq]*VIN[cvar1nos[i]-1] ;
	}
    }
#endif
}

/********************Hessian of objective times vector**********************/

void eval_hessobj_v(fint *TASK, fint *N, real *X, fint *M, real *VIN,
		    real *VOUT, real *DAT, fint *IDAT)
{
  real OW , temp , *cdummy ;
#ifdef INCLUDE_CC
  fint totsize ;
#endif

  /* Dummy call to ensure evaluation of Hessian at X */
  if (n_obj > 0)
    temp = objval(nobj, X, &nerror) ;

  /* Error handling */
  if ( nerror < 0 )
    {
      fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
      fprintf(Stderr," nerror = %d\n",nerror) ;
      fprintf(Stderr,"Exiting from subroutine - eval_hessobj_v_\n") ;
      exit(1) ;
    }

  /* Calls to evaluate constraint values */

  if( !(cdummy = (real *)Malloc(nzc*sizeof(real))) ) MALLOCERR(cdummy) ;
  conval(X, cdummy, &nerror) ;
  /* Error handling */
  if ( nerror < 0 )
    {
      fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
      fprintf(Stderr," nerror = %d\n",nerror) ;
      fprintf(Stderr,"Exiting from subroutine - eval_hesslag_v_\n") ;
      exit(1) ;
    }
  free(cdummy);

  if (n_obj > 0)
    if (objtype[nobj])
      OW = -1 ;
    else
      OW = 1 ;

  if( lasthcall == EVINIT || lasthcall != EVOBJ )
    hvinit(-1,&OW,NULL) ;
  hvcomp(VOUT,VIN,-1,&OW,NULL) ;
  lasthcall = EVOBJ;

  /* Copy zeros into the remaining part of the vector */
  if ( numIneq > 0 )
    F77_FUNC(dcopy,DCOPY)(&numIneq,&realzero,&intzero,&VOUT[n_var],&one) ;

#ifdef INCLUDE_CC

  /* Add zeros */

  if ( numVineq > 0 )
    {
      totsize = 2*numVineq ;
      F77_FUNC(dcopy,DCOPY)
	(&totsize,&realzero,&intzero,&VOUT[n_var+numIneq],&one) ;
    }
  if (n_cc > 0)
    if ( opIMPEC_CCAGG == 0 )
      {
	totsize = numCcon + 2*numVineq ;
	F77_FUNC(dcopy,DCOPY)
	  (&totsize,&realzero,&intzero,&VOUT[n_var+numIneq+2*numVineq],&one) ;
      }
    else
      VOUT[n_var+numIneq+2*numVineq] = 0 ;
#endif
}

/**************Hessain of the constraints times vector***************/

void eval_hesscon_v(fint *TASK, fint *N, real *X, fint *M, real *LAM,
		    real *VIN, real *VOUT, real *DAT, fint *IDAT)
{
  real temp , *cdummy ;
#ifdef INCLUDE_CC
  fint i, totsize ;
#endif

  /* Dummy call to ensure evaluation of Hessian at X */
  if (n_obj > 0)
    temp = objval(nobj, X, &nerror) ;

  /* Error handling */
  if ( nerror < 0 )
    {
      fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
      fprintf(Stderr," nerror = %d\n",nerror) ;
      fprintf(Stderr,"Exiting from subroutine - eval_hessobj_v_\n") ;
      exit(1) ;
    }

  /* Calls to evaluate constraint values */

  if( !(cdummy = (real *)Malloc(nzc*sizeof(real))) ) MALLOCERR(cdummy) ;
  conval(X, cdummy, &nerror) ;
  /* Error handling */
  if ( nerror < 0 )
    {
      fprintf(Stderr,"Error detected in AMPL evaluation!\n") ;
      fprintf(Stderr," nerror = %d\n",nerror) ;
      fprintf(Stderr,"Exiting from subroutine - eval_hesslag_v_\n") ;
      exit(1) ;
    }
  free(cdummy);

  if( lasthcall == EVINIT || lasthcall != EVCON )
    hvinit(-1,NULL,LAM) ;
  hvcomp(VOUT,VIN,-1,NULL,LAM) ;
  lasthcall = EVCON;

  /* Copy zeros into the remaining part of the vector */
  if ( numIneq > 0 )
    F77_FUNC(dcopy,DCOPY)(&numIneq,&realzero,&intzero,&VOUT[n_var],&one) ;

#ifdef INCLUDE_CC

  if ( numVineq > 0 )
    {
      totsize = 2*numVineq ;
      F77_FUNC(dcopy,DCOPY)
	(&totsize,&realzero,&intzero,&VOUT[n_var+numIneq],&one) ;
    }
  if ( n_cc > 0 ) 
    if ( opIMPEC_CCAGG == 0 )
      {
	totsize = numCcon + 2*numVineq ;
	F77_FUNC(dcopy,DCOPY)
	  (&totsize,&realzero,&intzero,&VOUT[n_var+numIneq+2*numVineq],&one) ;
      }
    else
      VOUT[n_var+numIneq+2*numVineq] = 0 ;
  
  for ( i = 0 ; i < (numCcon + 2*numVineq) ; i++ )
    {
      /* First check for type of bounds on the complementary variables */
      if ( (cvar1LBnds[i] < minREAL && cvar2LBnds[i] < minREAL) || (cvar1UBnds[i] > maxREAL && cvar2UBnds[i] > maxREAL) )
	temp = 1 ;
      else
	temp = -1 ;
      
      if ( opIMPEC_CCAGG == 0 )
	{	
	  VOUT[cvar1nos[i]-1] += temp*LAM[n_con+numVineq+i]*VIN[cvar2nos[i]-1] ;
	  VOUT[cvar2nos[i]-1] += temp*LAM[n_con+numVineq+i]*VIN[cvar1nos[i]-1] ;
	}
      else
	{
	  VOUT[cvar1nos[i]-1] += temp*LAM[n_con+numVineq]*VIN[cvar2nos[i]-1] ;
	  VOUT[cvar2nos[i]-1] += temp*LAM[n_con+numVineq]*VIN[cvar1nos[i]-1] ;
	}
    }
#endif
}

/************************evaluating constraints residual**********************/
/*** This is to be evaluated whenver mu is changed by IPOPT to obtain the  ***/
/*** correct complementarity constraint residuls for original NLP and      ***/
/*** relaxed barrier problem.                                              ***/
#ifdef INCLUDE_CC
void F77_FUNC_(eval_ccerr,EVAL_CCERR)
     (fint *N,fint *M, real *X,real *C, real *MU)
{
  fint i, presize ;
  real sum ;

  presize = n_var + numIneq + 2*numVineq ; /* number of variables without
					      taking into account the number
					      of relaxed compl. slacks */
  sum = 0 ; /* sum of the complementarity constraints for the aggregated case
	     */

  for ( i = 0 ; i < (numCcon + 2*numVineq) ; i++ )
    {
      /* Look for different cases and write the complementarity relation */
      
      /* l1 <= s <= Infinity complements l2 <= x <= Infinity */
      if ( (cvar1UBnds[i] > maxREAL) && (cvar2UBnds[i] > maxREAL) )
	if ( opIMPEC_CCAGG == 0 )
	  if (*MU > 0)
	    C[n_con+numVineq+i] = (X[cvar1nos[i]-1] - cvar1LBnds[i])*(X[cvar2nos[i]-1] - cvar2LBnds[i]) + X[presize+i] - opDMPEC_CCFACT*pow(*MU,opDMPEC_CCPOW) ;
	  else
	    C[n_con+numVineq+i] = (X[cvar1nos[i]-1] - cvar1LBnds[i])*(X[cvar2nos[i]-1] - cvar2LBnds[i]) ;
	else
	  sum = sum + (X[cvar1nos[i]-1] - cvar1LBnds[i])*(X[cvar2nos[i]-1] - cvar2LBnds[i]) ;
      /* l1 <= s <= Infinity complements -Infinity <= x <= u2 */
      else if ( (cvar1UBnds[i] > maxREAL) && (cvar2LBnds[i] < minREAL) )
	if ( opIMPEC_CCAGG == 0 )
	  if (*MU > 0)
	    C[n_con+numVineq+i] = (X[cvar1nos[i]-1] - cvar1LBnds[i])*(cvar2UBnds[i] - X[cvar2nos[i]-1]) + X[presize+i] - opDMPEC_CCFACT*pow(*MU,opDMPEC_CCPOW) ;
	  else
	    C[n_con+numVineq+i] = (X[cvar1nos[i]-1] - cvar1LBnds[i])*(cvar2UBnds[i] - X[cvar2nos[i]-1]) ;
	else
	  sum = sum + (X[cvar1nos[i]-1] - cvar1LBnds[i])*(cvar2UBnds[i] - X[cvar2nos[i]-1]) ;
      /* -Infinity <= s <= u1 complements l1 <= x <= Infinity */
      else if ( (cvar1LBnds[i] < minREAL) && (cvar2UBnds[i] > maxREAL) )
	if ( opIMPEC_CCAGG == 0 )
	  if (*MU > 0)
	    C[n_con+numVineq+i] = (cvar1UBnds[i] - X[cvar1nos[i]-1])*(X[cvar2nos[i]-1] - cvar2LBnds[i]) + X[presize+i] - opDMPEC_CCFACT*pow(*MU,opDMPEC_CCPOW) ;
	  else
	    C[n_con+numVineq+i] = (cvar1UBnds[i] - X[cvar1nos[i]-1])*(X[cvar2nos[i]-1] - cvar2LBnds[i]) ;
	else
	  sum = sum + (cvar1UBnds[i] - X[cvar1nos[i]-1])*(X[cvar2nos[i]-1] - cvar2LBnds[i]) ;
      /* -Infinity <= s <= u1 complements -Infinity <= x <= u2 */
      else if ( (cvar1LBnds[i] < minREAL) && (cvar2LBnds[i] < minREAL) )
	if ( opIMPEC_CCAGG == 0 )
	  if (*MU > 0)
	    C[n_con+numVineq+i] = (cvar1UBnds[i] - X[cvar1nos[i]-1])*(cvar2UBnds[i] - X[cvar2nos[i]-1]) + X[presize+i] - opDMPEC_CCFACT*pow(*MU,opDMPEC_CCPOW) ;
	  else
	    C[n_con+numVineq+i] = (cvar1UBnds[i] - X[cvar1nos[i]-1])*(cvar2UBnds[i] - X[cvar2nos[i]-1]) ;
	else
	  sum = sum + (cvar1UBnds[i] - X[cvar1nos[i]-1])*(cvar2UBnds[i] - X[cvar2nos[i]-1]) ;      
    }

  if (n_cc > 0 && opIMPEC_CCAGG == 1)
    if (*MU > 0)
      C[n_con+numVineq] = sum + X[presize] - opDMPEC_CCFACT*pow(*MU,opDMPEC_CCPOW) ;
    else
      C[n_con+numVineq] = sum ;
}
#endif

/********************************** main subroutine **************************/
int main(int argc, char **argv)
{

#ifdef sgi
#include <sys/fpu.h>
#endif

  /* Variables used in the call to IPOPT. */
  fint N, M, NLB, *ILB, NUB, *IUB, LIW, LRW , *IW, ITER, IERR, IDAT ;
  real *X, *BNDS_L, *BNDS_U, *V_L, *V_U, *LAM, *C, *RW, DAT ;

  fint i , j , l , clen , vlen , bound, numLBnds ;
  real optObj, conValue, mult;

#ifdef INCLUDE_CC
  real maxlam ;
  fint m, k ;
#endif

  char *stub ;
  char message[100];
  FILE *nl, *soln ;

  /* Array for passing options */
  static fint noptions = 100;  /* This is a generous upper bound of
				  the number of possible options */
  fint poptions ;        /* number of options set thru AMPL */
  real *optvalues ;
  char *optnames ;

#ifdef WRITE_IPOUT
  fint FFUNCS, CFUNCS, NOCG, NORES, NONEGCURV, ex, fe, nb ;
  real cnrm;
#endif

  /* For running using script */
  FILE *script ;
  /*clock_t clkStart , clkEnd ;*/
  real clkStart , clkEnd ;

#ifdef sgi
  union fpc_csr f;
  f.fc_word = get_fpc_csr();
  f.fc_struct.flush = 0;
  set_fpc_csr(f.fc_word);
#endif

  if ( argc < 2 )
    {
      printf("usage: %s stub\n", argv[0]) ;
      exit(1) ;
    }

  ASL_alloc(ASL_read_pfgh) ;

  stub = getstub(&argv,&Oinfo) ;

  if (getopts(argv,&Oinfo))
    return 1 ;

  /* Set problem constants and dimensions in .nl file */
  nl = jac0dim(stub, (fint)strlen(stub)) ;

  /* Allocate memory for initial guess. */
  if( !(X0 = (real *)Malloc(n_var*sizeof(real))) ) MALLOCERR(X0) ;

  /* Read the file */
  pfgh_read(nl,ASL_findgroups) ;
  hesset(1, nobj, 1, 0, nlc);

  /* Check if complementarity constraints are present but INCLUDE_CC is undefined */
#ifndef INCLUDE_CC
  if (n_cc > 0)
    {
      printf("Error : Cannot handle complementarity constraints.\n") ;
      exit(1) ;
    }
#else
  if (n_cc > 0 && opDMOVEBOUNDS == nullrealopt )
    {
      /* Different default for MOVEBOUNDS */
      opDMOVEBOUNDS = 0.;
    }
#endif

  /* Check for incorrect options setting */
#ifndef USE_MALLOC
  if ( opIMAX_IW <= 0 )
    {
      fprintf(Stderr,"Invalid value for imax_iw: %d\n",opIMAX_IW) ;
      exit(1) ;
    }
  if ( opIMAX_RW <= 0 )
    {
      fprintf(Stderr,"Invalid value for imax_rw: %d\n",opIMAX_RW) ;
      exit(1) ;
    }
#endif

  /* Identify the number of equalities, inequalities, number of lower and upper     bounds */

  numEq = 0 ;
  numIneq = 0 ;

  if( !(conL = (real *)Malloc(n_con*sizeof(real))) ) MALLOCERR(conL) ;
  if( !(conU = (real *)Malloc(n_con*sizeof(real))) ) MALLOCERR(conU) ;

  for ( i = 0 ; i < n_con ; i++ )
    {
      conL[i] = LUrhs[2*i] ;
      conU[i] = LUrhs[2*i+1] ;

      if ( conL[i] <= minREAL )
	conL[i] = minREAL ;
      if ( conU[i] >= maxREAL )
	conU[i] = maxREAL ;

      if ( conL[i] == conU[i] )
	numEq = numEq + 1 ;
      else
	numIneq = numIneq + 1 ;
    }		

  /* Obtain the number of given Complementarity constraints involved .
     Do not count complementarity relations involving equations */
#ifdef INCLUDE_CC
  numCcon = 0 ;
  numVineq = 0 ;

  if ( n_cc > 0 )
    for ( i = 0 ; i < n_con ; i++ )
      {
	if ( cvar[i] > 0 && cvar[i] <= n_var )
	  {
	    /* Check if any one of them has equal lower and upper bounds */
	    if ( (conL[i] != conU[i]) && (LUv[2*(cvar[i]-1)] != LUv[2*(cvar[i]-1)+1]) )
	      {
		/* Check if any of them is unconstrained */
		if ( (LUrhs[2*i] < minREAL) && (LUrhs[2*i+1] > maxREAL) || (LUv[2*(cvar[i]-1)] < minREAL) && (LUv[2*(cvar[i]-1)+1] > maxREAL) )
		  numVineq = numVineq + 1 ;
		else
		  numCcon = numCcon + 1 ;
	      }
	  }
      }
#endif

  /* Set the values of IPOPT call parameters. */
  if ( n_cc == 0 )
    {
      N = n_var + numIneq ;
      M = n_con ;
    }
#ifdef INCLUDE_CC
  else if ( opIMPEC_CCAGG == 0 )
    {
      N = n_var + numIneq + 2*numVineq + numCcon + 2*numVineq ;
      M = n_con + numCcon + 3*numVineq ;
    }
  else if ( opIMPEC_CCAGG == 1 )
    {
      N = n_var + numIneq + 2*numVineq + 1 ;
      M = n_con + numVineq + 1 ;
    }
#endif

  /* Identify the inequality constraint numbers */

  if ( numIneq > 0 )
    {
      if( !(slacks = (real *)Malloc(numIneq*sizeof(real))) ) MALLOCERR(slacks) ;
      if( !(ineqnos = (fint *)Malloc(numIneq*sizeof(fint))) ) MALLOCERR(ineqnos) ;
      j = 0 ;

      for ( i = 0 ; i < n_con ; i++ )
	if ( conL[i] != conU[i] )
	  {
	    ineqnos[j] = i+1 ;
	    j++ ;
	  }
    }

  /* Calculate the number of upper and lower bounds of the variables given */
  NLB = 0 ;
  NUB = 0 ;

  for ( i = 0 ; i < n_var ; i++ )
    {
      if ( LUv[2*i] > minREAL )
	NLB = NLB + 1 ;
      if ( LUv[2*i+1] < maxREAL )
	NUB = NUB + 1 ;
    }

  for ( i = 0 ; i < numIneq ; i++ )
    {
      if ( LUrhs[2*(ineqnos[i]-1)] > minREAL )
	NLB = NLB + 1 ;
      if ( LUrhs[2*(ineqnos[i]-1)+1] < maxREAL )
	NUB = NUB + 1 ;
    }

  /* Compensate for the additional slacks added for variational inequalities */
#ifdef INCLUDE_CC
  if ( n_cc > 0 )
    {
      if ( opIMPEC_CCAGG == 0 )
	NLB = NLB + 2*numVineq + numCcon + 2*numVineq ;
      else if ( opIMPEC_CCAGG == 1 )
	NLB = NLB + 2*numVineq + 1 ;
    }
#endif

  /* Store the variable lower and upper bounds appropriately. */

  if( !(BNDS_L = (real *)Malloc(NLB*sizeof(real))) ) MALLOCERR(BNDS_L) ;
  if( !(ILB = (fint *)Malloc(NLB*sizeof(fint))) ) MALLOCERR(ILB) ;
  if( !(BNDS_U = (real *)Malloc(NUB*sizeof(real))) ) MALLOCERR(BNDS_U) ;
  if( !(IUB = (fint *)Malloc(NUB*sizeof(fint))) ) MALLOCERR(IUB) ;

  j = 0 ; /* index the ILB */
  l = 0 ; /* index the IUB */

  for ( i = 0 ; i < n_var ; i++ )
    {
      if ( LUv[2*i] > minREAL )
	{
	  ILB[j] = i+1 ;
	  BNDS_L[j] = LUv[2*i] ;
	  j++ ;
	}
      if ( LUv[2*i+1] < maxREAL )
	{
	  IUB[l] = i+1 ;
	  BNDS_U[l] = LUv[2*i+1] ;
	  l++ ;
	}
    }

  for ( i = 0 ; i < numIneq ; i++ )
    {
      if ( LUrhs[2*(ineqnos[i]-1)] > minREAL )
	{
	  ILB[j] = n_var+i+1 ;
	  BNDS_L[j] = LUrhs[2*(ineqnos[i]-1)] ;
	  j++ ;
	}
      if ( LUrhs[2*(ineqnos[i]-1)+1] < maxREAL )
	{
	  IUB[l] = n_var+i+1 ;
	  BNDS_U[l] = LUrhs[2*(ineqnos[i]-1)+1] ;
	  l++ ;
	}
    }
  numLBnds = j ;

  /* Store the bounds for the complementarity variables */
#ifdef INCLUDE_CC
  if ( n_cc > 0 )
    {

      if( !(vineqvar = (fint *)Malloc(numVineq*sizeof(fint))) ) MALLOCERR(vineqvar) ;
      if( !(cvar1nos = (fint *)Malloc((numCcon + 2*numVineq)*sizeof(fint))) ) MALLOCERR(cvar1nos) ;
      if( !(cvar2nos = (fint *)Malloc((numCcon + 2*numVineq)*sizeof(fint))) ) MALLOCERR(cvar2nos) ;
      if( !(cvar1LBnds = (real *)Malloc((numCcon + 2*numVineq)*sizeof(real))) ) MALLOCERR(cvar1LBnds) ;
      if( !(cvar1UBnds = (real *)Malloc((numCcon + 2*numVineq)*sizeof(real))) ) MALLOCERR(cvar1UBnds) ;
      if( !(cvar2LBnds = (real *)Malloc((numCcon + 2*numVineq)*sizeof(real))) ) MALLOCERR(cvar2LBnds) ;
      if( !(cvar2UBnds = (real *)Malloc((numCcon + 2*numVineq)*sizeof(real))) ) MALLOCERR(cvar2UBnds) ;
      if( !(cslacks = (real *)Malloc(2*numVineq*sizeof(real))) ) MALLOCERR(cslacks) ;
      if( !(rslacks = (real *)Malloc((numCcon + 2*numVineq)*sizeof(real))) ) MALLOCERR(rslacks) ;

      j = 0 ; /* index the inequality slack variables */
      l = 0 ; /* index the variational inequality slack variables */
      k = 0 ; /* index the array varnos and bounds */
      m = 0 ; /* index the bounds to be supplied to IPOPT */
      for ( i = 0 ; i < n_con ; i++ )
	{
	  if ( cvar[i] > 0 & cvar[i] <= n_var )
	    {
	      /* Check for equalities and avoid them */
	      if ( (conL[i] != conU[i]) && (LUv[2*(cvar[i]-1)] != LUv[2*(cvar[i]-1)+1]) )
		{
		  /* Check for unrestricted expression - variational inequality
		   */
		  if ( (LUrhs[2*i] < minREAL) && (LUrhs[2*i+1] > maxREAL) )
		    {
		      /* -Infinity < s < Infinity complements l <= x <= u
			 Reformulate as,
			        s = s1 - s2
				s1 >= 0 complements x >= l    (C1)
				s2 >= 0 complements x <= u    (C2)        */


		      /* Index of slack introduced for unbndd constraint */
		      vineqvar[l] = n_var+j+1 ; /* Inequality slack is
						    unbounded */
		      
		      /* Complementarity relation - (C1) */
		      cvar1nos[k] = n_var+numIneq+2*l+1 ; /* complementarity
							     slack 1 (+) */
		      cvar1LBnds[k] = 0 ;
		      cvar1UBnds[k] = myInf ;
		      cvar2nos[k] = cvar[i] ; /* complentary variable */
		      cvar2LBnds[k] = LUv[2*(cvar[i]-1)] ;
		      cvar2UBnds[k] = myInf ;
		      ILB[numLBnds+m] = cvar1nos[k] ;
		      BNDS_L[numLBnds+m] = 0 ;
		      k++ ;
		      m++ ;

		      /* Complementarity relation - (C2) */
		      cvar1nos[k] = n_var+numIneq+2*l+2 ; /* complementarity
							     slack 2 (-) */
		      cvar1LBnds[k] = 0 ;
		      cvar1UBnds[k] = myInf ;
		      cvar2nos[k] = cvar[i] ; /* complementary variable */
		      cvar2LBnds[k] = myNegInf ;
		      cvar2UBnds[k] = LUv[2*(cvar[i]-1)+1] ;
		      ILB[numLBnds+m] = cvar1nos[k] ;
		      BNDS_L[numLBnds+m] = 0 ;
		      l++ ;
		      k++ ;
		      m++ ;
		      j++ ;

		    }
		  else if ( (LUv[2*(cvar[i]-1)] < minREAL) && (LUv[2*(cvar[i]-1)+1] > maxREAL) )
		    {
		      /* l <= s <= u complements -Infinity < x < Infinity
			 Reformulate as,
			            x = s1 - s2
				    s1 >= 0 complements s >= l  (C1)
				    s2 >= 0 complements s <= u  (C2)       */

		      /* variable index for unbndd variable */
		      vineqvar[l] = cvar[i] ; /* AMPL variable is unbounded */
      
		      /* Complementarity relation - (C1) */
		      cvar1nos[k] = n_var+numIneq+2*l+1 ; /* Complementarity
							     slack 1 (+) */
		      cvar1LBnds[k] = 0 ;
		      cvar1UBnds[k] = myInf ;
		      cvar2nos[k] = n_var+j+1 ; /* Inequality slack */
		      cvar2LBnds[k] = LUrhs[2*i] ;
		      cvar2UBnds[k] = myInf ;
		      ILB[numLBnds+m] = cvar1nos[k] ;
		      BNDS_L[numLBnds+m] = 0 ;
		      k++ ;
		      m++ ;

		      /* Complementarity relation - (C2) */
		      vineqvar[k] = 0 ;
		      cvar1nos[k] = n_var+numIneq+2*l+2 ; /* Complementarity
							     slack 2 (-) */
		      cvar1LBnds[k] = 0 ;
		      cvar1UBnds[k] = myInf ;
		      cvar2nos[k] = n_var+j+1 ; /* Inequality slack */
		      cvar2LBnds[k] = myNegInf ;
		      cvar2UBnds[k] = LUrhs[2*i+1] ;
		      ILB[numLBnds+m] = cvar1nos[k] ;
		      BNDS_L[numLBnds+m] = 0 ;
		      l++ ;
		      k++ ;
		      m++ ;
		      j++ ;
		    }
		  else
		    {
		      /* It is one of,
			   l1 <= s complements l2 <= x
			   l1 <= s complements u2 >= x
			   u1 >= s complements l2 <= x
			   u1 >= s complements u2 >= x  */
		
		      cvar1nos[k] = n_var+j+1 ; /* Inequality slack */
		      cvar1LBnds[k] = LUrhs[2*i] ;
		      cvar1UBnds[k] = LUrhs[2*i+1] ;
		      cvar2nos[k] = cvar[i] ;
		      cvar2LBnds[k] = LUv[2*(cvar[i]-1)] ;
		      cvar2UBnds[k] = LUv[2*(cvar[i]-1)+1] ;
		      k++ ;
		      j++ ;
		    }
		}
	    }
	  /* An inequality not involved in a complementarity relation */
	  else if ( ineqnos[j] == i+1 )
	    j++ ;
	}
    }
  numLBnds = numLBnds + m ;
#endif

  /* Bounds for the variables introduced as part of relaxing the
     complementarity equation */
#ifdef INCLUDE_CC
  if ( n_cc > 0 )
    {
      j = (opIMPEC_CCAGG == 1) ? 1 : (numCcon + 2*numVineq) ;
      for ( i = 0 ; i < j ; i++ )
	{
	  ILB[numLBnds+i] = n_var + numIneq + 2*numVineq + i + 1 ;
	  BNDS_L[numLBnds+i] = 0 ;
	}
    }
#endif

  /* Initialise the slack variables */
  if ( numIneq > 0 )
    initslacks(X0,slacks,ineqnos) ;

  /* Initialise the slacks for variational inequalities */
#ifdef INCLUDE_CC
  if ( n_cc > 0 )
    initcompslacks(X0,cslacks,slacks,ineqnos) ;
#endif

  /* Initialise the slacks for each relaxed complementarity constraint */
#ifdef INCLUDE_CC
  if ( n_cc > 0 )
    initrslacks(X0,rslacks) ;
#endif

  /* Copy initial point into X */
  if( !(X = (real *)Malloc(N*sizeof(real))) ) MALLOCERR(X) ;

  for ( i = 0 ; i < N ; i++ )
    {
      if ( i < n_var )
	X[i] = X0[i] ;
      else if ( i < (n_var + numIneq) )
	X[i] = slacks[i-n_var] ;
#ifdef INCLUDE_CC
      else if ( i < (n_var + numIneq + 2*numVineq) )
	X[i] = cslacks[i-n_var-numIneq] ;
      else
	X[i] = rslacks[i-n_var-numIneq-2*numVineq] ;
#endif
    }

  /* Allocate memory for the Multipliers and constraints */
  if( !(V_L = (real *)Malloc(NLB*sizeof(real))) ) MALLOCERR(V_L) ;
  if( !(V_U = (real *)Malloc(NUB*sizeof(real))) ) MALLOCERR(V_U) ;
  if( !(LAM = (real *)Malloc(M*sizeof(real))) ) MALLOCERR(LAM) ;
  if( !(C = (real *)Malloc(M*sizeof(real))) ) MALLOCERR(C) ;

#ifdef USE_MALLOC
  LRW = 0;
  LIW = 0;
  RW  = NULL;
  IW  = NULL;
#else
  /* Initialize the workspace variables */
  LRW = opIMAX_RW ;
  LIW = opIMAX_IW ;

  if( !(RW = (real *)Malloc(LRW*sizeof(real))) ) MALLOCERR(RW) ;
  if( !(IW = (fint *)Malloc(LIW*sizeof(fint))) ) MALLOCERR(IW) ;
#endif

  /* Read variable names */
  vlen = maxcolnamelen ;
  clen = maxrownamelen ;

#ifdef WRITE_IPOUT
  /* Write ip.out as one-line output for iterlist */
  nb = (NUB>NLB) ? NUB : NLB ;
  F77_FUNC_(check_flagfile,CHECK_FLAGFILE)(&ex) ;
  if ( ex )
    {
      script = fopen("ip.out","w") ;
      fprintf(script,"%18s | %5d %5d %5d | %5d/%5d %5dc%5dr%5dn (%3d) %15.8e %10.3e %9.2e",stub,N,M,nb,0,-1,0,0,0,-1,0e0,0e0,-1e0 ) ;
      fclose(script) ;
    }
  else
    {
      script = fopen("ip.out","w") ;
      fprintf(script," | %5d/%5d %5dc%5dr%5dn (%3d) %15.8e %10.3e %9.2e",0,-1,0,0,0,-1,0e0,0e0,-1e0 )
;
      fclose(script) ;
    }
#endif

  /* Write SLACKS.DAT to tell IPOPT about the slack variables */
  if( numIneq > 0 )
    {
      script = fopen("SLACKS.DAT","w") ;
      for( i = n_var+1 ; i <= (n_var + numIneq) ; ++i)
        {
          fprintf(script,"%6d\n",i) ;
        }
      fclose(script) ;
    }
  else
    {
      remove("SLACKS.DAT") ;
    }


  /* Pass parameters to IPOPT */
#define CARGSLEN 40
  if( !(optnames = (char *)Malloc(noptions*CARGSLEN*sizeof(char))) ) MALLOCERR(optnames) ;
  if( !(optvalues = (real *)Malloc(noptions*sizeof(real))) ) MALLOCERR(optvalues) ;
  j = 0;
  i = 0;

#define SetIARGS(PARAM) { if ( (op ## PARAM) != nullintopt) {optvalues[i++] = op ## PARAM; str2fstr(#PARAM,optnames+j,CARGSLEN);j+=CARGSLEN;} }

#define SetRARGS(PARAM) { if ( (op ## PARAM) != nullrealopt) {optvalues[i++] = op ## PARAM; str2fstr(#PARAM,optnames+j,CARGSLEN);j+=CARGSLEN;} }

  SetIARGS(IFILE);
  SetIARGS(IPRINT);
  SetIARGS(IOUTPUT);
  SetIARGS(IMAXITER);
  SetIARGS(IMAXCPUSEC);
  SetIARGS(IMUINIT);
  SetIARGS(IQUASI);
  SetIARGS(ISCALE);
  SetIARGS(ISCALERR);
  SetIARGS(IMERIT);
  SetIARGS(ISOC);
  SetIARGS(IREFINEITER);
  SetIARGS(IRESTO);
  SetIARGS(ICNRM);
  SetIARGS(IFULL);

  SetRARGS(DTOL);
  SetRARGS(DFSCALE);
  SetRARGS(DINFMAXTOL);
  SetRARGS(DMU0);
  SetRARGS(DBNDFRAC);
  SetRARGS(DBNDPUSH);
  SetRARGS(DPIVTOL);
  SetRARGS(DPIVTOLMAX);
  SetRARGS(DTRONCGTOL);
  SetRARGS(DCMAXTOL);
  SetRARGS(DDEPCONDIAG);
  SetRARGS(DSCALECUT);
  SetRARGS(DMOVEBOUNDS);
  SetRARGS(DLAMINITMAX);
  SetRARGS(DFILLINFACT);

#ifndef FEWOPTIONS
  SetIARGS(ICG);
  SetIARGS(IINITB);
  SetIARGS(ILMLEN);
  SetIARGS(ICORRECT);
  SetIARGS(IDAMP);
  SetIARGS(IALPHA);
  SetIARGS(IITERBLOCKMAX);
  SetIARGS(ISELBAS);
  SetIARGS(ISYMSOLV);
  SetIARGS(IHESSVECT);
  SetIARGS(ITRON2DERIVS);
  SetIARGS(ITRONHESS);
  SetIARGS(IKKTSCALE);
  SetRARGS(DLS_SAFE);
  SetRARGS(DMULIN);
  SetRARGS(DMUSUPER);
  SetRARGS(DNUMIN);
  SetRARGS(DSR1TOL);
  SetRARGS(DSKIPFACT);
  SetRARGS(DRHO);
  SetRARGS(DMAXCOND);
  SetRARGS(DCGTOL);
  SetRARGS(DRESTOKKTRED);
  SetRARGS(DS_F);
  SetRARGS(DS_THETA);
  SetRARGS(DGAMMA_F);
  SetRARGS(DGAMMA_THETA);
  SetRARGS(DTAU);
  SetRARGS(DVCORRECTFACT);
  SetRARGS(DERRSUPER);
  SetRARGS(DMUERRFAC);
  SetRARGS(DMAXERR);
  SetRARGS(DWATCHTOL);
#endif
#ifdef INCLUDE_CC
#ifndef FEWOPTIONS
  SetIARGS(IMPEC_TRIGGER) ;
  SetRARGS(DMPEC_ETAFACT) ;
  SetRARGS(DMPEC_THRESH) ;
#endif
#endif
#undef SetIARGS
#undef SetRARGS

  /* actual options asigned */
  poptions = i ;

  printf("\n");
  /* clkStart = clock() ; */
  F77_FUNC(timer,TIMER)(&clkStart) ;
  /* Call IPOPT to solve the NLP */
  optObj = F77_FUNC(ipopt,IPOPT)
    (&N, X, &M, &NLB, ILB, BNDS_L, &NUB, IUB, BNDS_U, V_L, V_U,
     LAM, C, &LRW, RW, &LIW, IW, &ITER, &IERR, 
     (void *(*)())eval_f, (void *(*)())eval_c,
     (void *(*)())eval_g, (void *(*)())eval_a, (void *(*)())eval_h,
     (void *(*)())eval_hesslag_v, (void *(*)())eval_hessobj_v,
     (void *(*)())eval_hesscon_v, &DAT, &IDAT,
     &poptions, optvalues, optnames, CARGSLEN) ;
  /* clkEnd = clock() ; */
  F77_FUNC(timer,TIMER)(&clkEnd) ;

#ifdef WRITE_IPOUT
  if( M>0 )
    {
#ifdef INCLUDE_CC
      if (n_cc > 0)
	F77_FUNC_(eval_ccerr,EVAL_CCERR)(&N,&M,X,C,&realzero) ;
#endif
      i = F77_FUNC(idamax,IDAMAX)(&M,C,&one);
      cnrm = C[i-1] ;
      cnrm = (cnrm >0) ? cnrm : -cnrm ;
    }
  else
    cnrm = 0;

  F77_FUNC_(ipopt_getdata,IPOPT_GETDATA)
    (&FFUNCS, &CFUNCS, &NOCG, &NORES, &NONEGCURV);
  fe = (FFUNCS>CFUNCS) ? FFUNCS : CFUNCS ;

  /* For running scripts - modification to suit Arvind and Andreas o/p spec */
#ifdef INCLUDE_CC
  if ( n_cc > 0 )
    {
      script = fopen("ipc.out","w") ;
      fprintf(script,"%d\t%d\t%d\t%e\t%d\t%d\t%e\t%d\t%d\t%d\t%s\n",N,M,numCcon,optObj,ITER,fe,cnrm,NORES,NONEGCURV,IERR,stub) ;
      fclose(script) ;
    }
#endif
  F77_FUNC_(check_flagfile,CHECK_FLAGFILE)(&ex) ;
  if ( ex )
    {
      script = fopen("ip.out","w") ;
      fprintf(script,"%18s | %5d %5d %5d | %5d/%5d %5dc%5dr%5dn (%3d) %15.8e %10.3e %9.2e",stub,N,M,nb,fe,ITER,NOCG,NORES,NONEGCURV,IERR,optObj,cnrm,(clkEnd-clkStart) ) ;
      fclose(script) ;
    }
  else
    {
      script = fopen("ip.out","w") ;
      fprintf(script," | %5d/%5d %5dc%5dr%5dn (%3d) %15.8e %10.3e %9.2e",fe,ITER,NOCG,NORES,NONEGCURV,IERR,optObj,cnrm,(clkEnd-clkStart) ) ;
      fclose(script) ;
    }

  /* Clean up a little */
  if(  numIneq > 0 )
    remove("SLACKS.DAT") ;
#endif

#ifdef INCLUDE_CC
  /* Display some information about the problem */

  printf("IPOPT variables                       = %d\n",N) ;
  printf("IPOPT constraints                     = %d\n",M) ;
  printf("Exit code                             = %d\n",IERR) ;
  if (n_cc > 0)
    {
      printf("IPOPT complementarity constraints     = %d\n",(numCcon-2*numVineq)) ;
      if (numVineq > 0)
	printf("IPOPT variational inequalities        = %d\n",numVineq) ;
      maxlam = 0 ;
      if (M > 0)
	{
	  i = F77_FUNC(idamax,IDAMAX)(&M,LAM,&one) ;
	  maxlam = (maxlam >= fabs(LAM[i-1])) ? maxlam : fabs(LAM[i-1]) ;
	}
      if (NLB > 0)
	{
	  i = F77_FUNC(idamax,IDAMAX)(&NLB,V_L,&one) ;
	  maxlam = (maxlam >= V_L[i-1]) ? maxlam : V_L[i-1] ;
	}
      if (NUB > 0)
	{
	  i = F77_FUNC(idamax,IDAMAX)(&NUB,V_U,&one) ;
	  maxlam = (maxlam >= V_U[i-1]) ? maxlam : V_U[i-1] ;
	}
      printf("IPOPT max. multiplier                 = %g\n",maxlam) ;
    }
#endif

  /* Process the output of IPOPT and display appropriate results.
     NLP.solution is written only if variable and constraint names are defined
     in AMPL */
  need_nl = printf(PACKAGE_STRING": ") ;
  need_nl = 0 ;

  if ((clen != 0) && (vlen != 0))
    {
      soln = fopen("NLP.solution", "w") ;
      fprintf(soln,"                  "PACKAGE_STRING"SOLVE SUMMARY \n") ;
      fprintf(soln,"                  ------------------------- \n") ;

      /* Output formatting */
      if ( n_obj == 0 )
	fprintf(soln,"\n                  Finding a Fesible Point\n") ;
      else
	{
	  if (!objtype[nobj])
	    fprintf(soln,"\n                  Minimization Problem\n") ;
	  else
	    fprintf(soln,"\n                  Maximization Problem\n") ;
	  if (!objtype[nobj])
	    fprintf(soln,"Objective value = %-16g\n",optObj) ;
	  else
	    fprintf(soln,"Objective value = %-16g\n",-optObj) ;
	}

      fprintf(soln,"Solve Status : ") ;
      if ( IERR == 0 )
	fprintf(soln,"OPTIMAL SOLUTION FOUND.\n\n") ;
      else
	fprintf(soln,"Solve process aborted. Refer to message on screen.\n\n") ;

      if (vlen != 0)
	{
	  j = 0 ; /* Track of variables with lower bounds */
	  l = 0 ; /* Track of variables with upper bounds */
	  fprintf(soln,"\nVariables\n") ;
	  fprintf(soln,"\nLower          Value         Upper        Multiplier       Name\n") ;

	  for ( i = 0 ; i < n_var ; i++ )
	    {
	      bound = 0 ;
	      if ( ILB[j] == i+1 && IUB[l] == i+1 )
		{
		  mult = (V_L[j] >= V_U[l]) ? V_L[j] : V_U[l] ;
		  j++ ;
		  l++ ;
		  bound = 1 ;
		}
	      else if ( ILB[j] == i+1 )
		{
		  mult = V_L[j] ;
		  j++ ;
		  bound = 1 ;
		}
	      else if ( IUB[l] == i+1 )
		{
		  mult = V_U[l] ;
		  l++ ;
		  bound = 1 ;
		}
	      if (bound == 1)
		fprintf(soln, "\n%+0.2e      %+0.2e     %+0.2e    %+0.2e        %s",LUv[2*i],X[i],LUv[2*i+1],mult,var_name(i)) ;
	      else
		fprintf(soln, "\n%+0.2e      %+0.2e     %+0.2e    .........        %s",LUv[2*i],X[i],LUv[2*i+1],var_name(i)) ;
	    }
	}

      xknown(X);
      if (clen != 0)
	{
	  fprintf(soln,"\n\nConstraints\n") ;
	  fprintf(soln,"\nLower          Value         Upper        Mulitplier       Name\n") ;
 	
	  for ( i = 0 ; i < n_con ; i++ )
	    {
	      conValue = conival(i,X,&nerror) ;
	      fprintf(soln,"\n%+0.2e      %+0.2e     %+0.2e    %+0.2e        %s",LUrhs[2*i],conValue,LUrhs[2*i+1],LAM[i],con_name(i)) ;
	    }
	}
      xunknown();

      fclose(soln) ;
    }

  /* Revert the sign of the multipliers, so that it matches AMPL's
     default definition */
  F77_FUNC(dscal,DSCAL)(&M,&realminusone,LAM,&one) ;
  if ( IERR == 0 )
    write_sol("OPTIMAL SOLUTION FOUND", X, LAM, 0) ;
  else
    {
      sprintf(message, "Exit code %d. Check IPOPT output for details.", IERR);
      write_sol(message, X, LAM, 0) ;
    }

  /* Free up the memory allocations. */
#ifndef USE_MALLOC
  free(IW) ;
  free(RW) ;
#endif
  free(ILB) ;
  free(IUB) ;
  free(X0) ;
  free(X) ;
  free(BNDS_L) ;
  free(BNDS_U) ;
  free(V_L) ;
  free(V_U) ;
  free(LAM) ;
  free(C) ;
  free(conL) ;
  free(conU) ;
  free(optvalues) ;
  free(optnames);

  if ( numIneq > 0 )
    {
      free(ineqnos) ;
      free(slacks) ;
    }
#ifdef INCLUDE_CC
  if ( n_cc > 0 )
    {
      free(vineqvar) ;
      free(cvar1nos) ;
      free(cvar2nos) ;
      free(cvar1LBnds) ;
      free(cvar1UBnds) ;
      free(cvar2LBnds) ;
      free(cvar2UBnds) ;
      free(cslacks) ;
      free(rslacks) ;
    }
#endif

  return 0 ;

}

#ifdef __cplusplus
}
#endif
