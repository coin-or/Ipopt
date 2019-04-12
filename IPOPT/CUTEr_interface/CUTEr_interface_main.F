C Copyright (C) 2002, Carnegie Mellon University, Dominique Orban and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
      PROGRAM           IPOPTMA
C
C     IPOPT CUTEr driver.
C     D. Orban,  adapted from Andreas Wachters CUTE driver.
C
      IMPLICIT NONE
      INTEGER IOUT
      PARAMETER( IOUT = 6 )

#include <config.h>

#include "CUTE.INC"

#ifdef USE_MALLOC
      INTEGER LRW, LIW
      PARAMETER( LRW=0, LIW=0 )
      INTEGER IW
      DOUBLE PRECISION RW
#else
      INTEGER LRW, LIW
CTOY  PARAMETER (LRW = 20000, LIW = 10000)
CMED  PARAMETER (LRW = 10000000, LIW = 10000000)
CBIG  PARAMETER (LRW = 30000000, LIW = 30000000)
CCUS  PARAMETER (LRW = 30000000, LIW = 30000000)
      INTEGER IW( LIW )
      DOUBLE PRECISION RW( LRW )
#endif

      INTEGER N, M, NLB, NUB
      DOUBLE PRECISION X( CUTE_NMAX )
      INTEGER ILB( CUTE_NMAX )
      INTEGER IUB( CUTE_NMAX )
      DOUBLE PRECISION BNDS_L( CUTE_NMAX )
      DOUBLE PRECISION BNDS_U( CUTE_NMAX )
      DOUBLE PRECISION V_L( CUTE_NMAX )
      DOUBLE PRECISION V_U( CUTE_NMAX )
      DOUBLE PRECISION LAM( CUTE_MMAX )
      DOUBLE PRECISION C( CUTE_MMAX )
      DOUBLE PRECISION IPOPT
C
C     Algorithmic Parameters (INITPARAMS)
C
      INTEGER NARGS
      DOUBLE PRECISION ARGS( 50 )
      CHARACTER*20 CARGS( 50 )

      INTEGER ITER
      INTEGER IERR

      EXTERNAL EVAL_F
      EXTERNAL EVAL_C
      EXTERNAL EVAL_G
      EXTERNAL EVAL_A
      EXTERNAL EVAL_H
      EXTERNAL EVAL_HESSLAG_V
      EXTERNAL EVAL_HESSOBJ_V
      EXTERNAL EVAL_HESSCON_V
C
C     The following arrays could be used instead of common blocks to
C     communicate between driver routines and evaluation subroutines
C
      DOUBLE PRECISION DAT(1)
      INTEGER IDAT(1)

      INTEGER fevals, cevals
      COMMON /EVALS/ fevals, cevals

      REAL CALLS( 7 ), CPU( 2 )
      CHARACTER*10 PNAME
      CHARACTER*10 VNAMES( CUTE_NMAX ), GNAMES( CUTE_MMAX )
      DOUBLE PRECISION FinalF, DummyG, cnrm
      INTEGER IDAMAX
C
      fevals = 0
      cevals = 0
C
C     Get problem dimensions and initialize
C
      CALL CUTE_INIT(N, M, CUTE_NMAX, X, NLB, ILB, BNDS_L,
     .     NUB, IUB, BNDS_U)
C
C     Get problem name.
C
      CALL CNAMES( CUTE_N, CUTE_M, PNAME, VNAMES, GNAMES )
C
C     Set algorithmic parameters (None :)
C
      NARGS = 0
C
C     Call IPOPT
C
      FinalF = IPOPT(N, X, M, NLB, ILB, BNDS_L, NUB, IUB, BNDS_U, V_L,
     1     V_U, LAM, C, LRW, RW, LIW, IW, ITER, IERR, EVAL_F, EVAL_C,
     2     EVAL_G, EVAL_A, EVAL_H, EVAL_HESSLAG_V, EVAL_HESSOBJ_V,
     3     EVAL_HESSCON_V, DAT, IDAT, NARGS, ARGS, CARGS)
C
C     Display CUTEr statistics
C
      CALL CREPRT( CALLS, CPU )
      WRITE ( IOUT, 2000 ) PNAME, CUTE_N, CUTE_M, CALLS(1), CALLS(2),
     .     CALLS(3), CALLS(4), CALLS(5), CALLS(6), CALLS(7),
     .     IERR, FinalF, CPU(1), CPU(2)
c
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *     ,/,' Code used               :  IPOPT',    /
     *     ,' Problem                 :  ', A10,    /
     *     ,' # variables             =      ', I10 /
     *     ,' # constraints           =      ', I10 /
     *     ,' # objective functions   =      ', E15.7 /
     *     ,' # objective gradients   =      ', E15.7 /
     *     ,' # objective Hessians    =      ', E15.7 /
     *     ,' # Hessian-vector prdct  =      ', E15.7 /
     *     ,' # constraints functions =      ', E15.7 /
     *     ,' # constraints gradients =      ', E15.7 /
     *     ,' # constraints Hessians  =      ', E15.7 /
     *     ,' Exit code               =      ', I10 /
     *     ,' Final f                 = ', E15.7 /
     *     ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ,' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     ,/,66('*') / )

 9999 CONTINUE
      END



