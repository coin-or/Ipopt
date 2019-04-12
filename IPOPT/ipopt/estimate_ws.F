C Copyright (C) 2003, International Business Machines and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine ESTIMATE_WS(N, X, M, NLB, NUB, FILLINFACT, LRW, LIW,
     1     IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     1     DAT, IDAT)
#ifdef BUILD_VSDLL
CDEC$ ATTRIBUTES DLLEXPORT::ESTIMATE_WS
#endif
C
C*******************************************************************************
C
C    $Id: estimate_ws.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Estimate work space requirement for IPOPT
CT    (INITPARAMS must have been called before!)
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CBC
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      02/213/03  First version
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   N         I    INT    number of variables in problem statement
CP   X         I    DP     starting point
CP   M         I    INT    number of equality constraints
CP   NLB       I    INT    number of lower bounds
CP   NUB       I    INT    number of upper bounds
CP   FILLINFACT I   DP     if compiled without USE_MALLOC:
CP                           your intelligent guess of a fill-in factor for
CP                           full-space KKT matrix...
CP   LRW       O    INT    minimum required double precision work space
CP   LIW       O    INT    minimum required integer work space
CP   IERR      O    INT    error flag (if nonzero, something went wrong)
CP   EV_F      I    EXT    Subroutine for objective function
CP   EV_C      I    EXT    Subroutine for constraints
CP   EV_G      I    EXT    Subroutine for gradient of objective function
CP   EV_A      I    EXT    Subroutine for Jacobian
CP   EV_H      I    EXT    Subroutine for Lagrangian Hessian
CP   EV_HLV    I    EXT    Subroutine for Lagrangian Hessian-vector products
CP   EV_HOV    I    EXT    Subroutine for objective Hessian-vector products
CP   EV_HCV    I    EXT    Subroutine for constraint Hessian-vector products
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CS
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
      include 'IPOPT.INC'
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      double precision X(N)
      integer M
      integer NLB
      integer NUB
      double precision FILLINFACT
      integer LRW
      integer LIW
      integer IERR
      external EV_F
      external EV_C
      external EV_G
      external EV_A
      external EV_H
      external EV_HLV
      external EV_HOV
      external EV_HCV
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer lrs, lis, nza, nzh, kconstr(6)

      double precision dummy(1)
      integer idummy(1)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
C
C     First determine size of store space and obtain number of nonzeros in
C     Jacobian
C
      lrs = 0
      lis = 0
      INMEMCHECK = 1
      call INIT_MEM(N, X, N, N-M, M, NLB, NUB, nza, nzh,
     1     KCONSTR, lrs, lis, 0, dummy, 0, idummy, IERR, EV_F, EV_C,
     5     EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.ne.0 ) return

      call IPOPT_WS(N, X, M, NLB, NUB, nza, nzh, LRW, LIW, FILLINFACT,
     1     EV_H, DAT, IDAT)

      LRW = LRW + lrs + N+NLB+NUB
      LIW = LIW + lis + N+NLB+NUB
      if( QSCALE.ge.3 ) then
         LRW = LRW + N+M
      elseif( QSCALE.ge.2 ) then
         LRW = LRW + M
      else
         LRW = LRW + 1
      endif

      return
      end
