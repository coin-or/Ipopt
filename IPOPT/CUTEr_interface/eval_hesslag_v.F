C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine EVAL_HESSLAG_V(TASK, N, X, M, LAM, VIN, VOUT,
     1     DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: eval_hesslag_v.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute product of Hessian of Lagrangian (individual constraints
CT    weighted by LAM) with vector
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      11/05/00
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
CP   TASK      I    INT     =0: reevaluate Hessian entries (because X or LAM
CP                              changed
CP                         <>0: no need to reevaluate
CP   N         I    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where Hessians are to be evaluated
CP   M         I    INT    number of constraints (including slack equations)
CP   LAM       I    DP     weights for constraints Hessians
CP   VIN       I    DP     vector to be mutliplied with Hessian
CP   VOUT      O    DP     resulting product
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
CS    CPROD
CS    DCOPY
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
#include "CUTE.INC"
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer TASK
      integer N
      integer M
      double precision LAM(M)
      double precision X(N)
      double precision VIN(N)
      double precision VOUT(N)
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      logical goth
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( TASK.eq.0 ) then
         goth = .false.
      else
         goth = .true.
      endif

      call CPROD(CUTE_N, CUTE_M, goth, X, M, LAM, VIN, VOUT)
C
C     slacks only appear linear
C
      if( N.gt.CUTE_N ) then
         call DCOPY(N-CUTE_N, 0d0, 0, VOUT(CUTE_N+1), 1)
      endif

 9999 continue
      return
      end
