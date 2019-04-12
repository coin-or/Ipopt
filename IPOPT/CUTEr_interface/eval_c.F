C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine EVAL_C(N, X, M, C, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: eval_c.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute values of constraints to CUTE problem
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
CA    Andreas Waechter      02/25/99
CA    Andreas Waechter      07/01/99   BUG: problems if ineq not first
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
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where G is to be evaluated
CP   M         I    INT    number of constraints
CP   C         O    DP     values of constraints
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
CS    CCFG
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
      integer N
      double precision X(N)
      integer M
      double precision C(M)
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision f, dummy
      integer i, liq, leq
      logical ineq
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Call CCFG to obtain constraint values, but without slacks
C
      call CCFG(CUTE_N, CUTE_M, X, M, C, .FALSE., 1, 1, dummy, .FALSE.)
C
C     Add entries for slack variables and constant terms
C
      liq = 1
      leq = 1
      do i = 1, M
         ineq = .false.
         if( liq.le.CUTE_NIQ ) then
            if( CUTE_IIQ(liq).eq.i ) then
               C(i) = C(i) - X(CUTE_N+liq)
               liq = liq + 1
               ineq = .true.
            endif
         endif
         if( .not.ineq ) then
            C(i) = C(i) - CUTE_CEQ(leq)
            leq = leq + 1
         endif
      enddo

 9999 continue
      return
      end
