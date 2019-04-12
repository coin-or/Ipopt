C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine EVAL_A(TASK, N, X, NZ, A, ACON, AVAR, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: eval_a.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute Jacobian of constraints to CUTE problem
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
CP   TASK      I    INT     =0: Obtain NZ
CP                         <>0: Compute Jacobian
CP   N         I    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where A is to be evaluated
CP   NZ       I/O   INT    TASK = 0: O: number of nonzero elements
CP                         otherwise: number of nonzero elements
CP                                     (size of A, AVAR, ACON)
CP   A         O    DP     (only TASK<>0) values in Jacobian
CP   ACON      O    INT    (only TASK<>0) row indices
CP   AVAR      O    INT    (only TASK<>0) column indices
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
CS    CDIMSJ
CS    CCFSG
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
      double precision X(N)
      integer NZ
      double precision A(NZ)
      integer ACON(NZ)
      integer AVAR(NZ)
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision c(CUTE_MMAX)
      double precision cjac(CUTE_NZMAX)
      integer indvar(CUTE_NZMAX), indfun(CUTE_NZMAX)
      integer i, nnzj
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( TASK.eq.0 ) then
C
C     Call CDIMSJ to obtain number of nonzero elements
C
         CALL CDIMSJ( NZ )
C
C     Substract contribution of (dense) objective function gradient
C
         NZ = NZ - CUTE_N
         NZ = NZ + CUTE_NIQ
      else
C
C     Call CCFSG to obtain Jacobian for constraints
C
         call CCFSG(CUTE_N, CUTE_M, X, CUTE_MMAX, c, nnzj,
     1        NZ, A, AVAR, ACON, .TRUE.)
C
C     Augment entries for slacks
C
         do i = 1, CUTE_NIQ
            A   (nnzj+i) = -1.d0
            AVAR(nnzj+i) = CUTE_N + i
            ACON(nnzj+i) = CUTE_IIQ(i)
         enddo

      endif

 9999 continue
      return
      end
