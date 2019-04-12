C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_D(N, NIND, M, X, YPY, PZ, ZPZ, NLB, ILB, NUB, IUB,
     1                 S_L, S_U, V_L, V_U, BNDS_L, BNDS_U,
     1                 SIGMA_L, SIGMA_U, MU, DX, DV_L, DV_U, RESTO)
C
C*******************************************************************************
C
C    $Id: get_d.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute primal and dual steps.
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
CA    Andreas Waechter      05/01/02  Release as version IPOPT 2.0
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
CP   N         I    INT    number of (free) variables; first M variables
CP                         are dependent; remaining independent
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of dependent variables
CP   X         I    DP     actual primal iterate
CP   YPY       I    DP     range space step (all variables; ordered as X)
CP   PZ        I    DP     null space step (indepentent variables)
CP   ZPZ       I    DP     null space step (dependent variables)
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. S_L(i) is slack for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. S_U(i) is slack for X(IUB(i)) )
CP   S_L       I    DP     slack variables for lower bounds
CP   S_U       I    DP     slack variables for upper bounds
CP   V_L       I    DP     dual variables for lower bounds
CP   V_U       I    DP     dual variables for upper bounds
CP   BNDS_L    I    DP     values of lower bounds (ordered as S_L)
CP   BNDS_U    I    DP     values of upper bounds (ordered as S_U)
CP   SIGMA_L   I    DP     primal-dual Hessian of lower bound barrier term
CP                            (NLB diagonal elements only)
CP   SIGMA_U   I    DP     primal-dual Hessian of upper bound barrier term
CP                            (NUB diagonal elements only)
CP   MU        I    DP     barrier parameter (=0, if Error should be computed
CP                            for overall NLP)
CP   DX        O    DP     step for X (primal)
CP   DV_L      O    DP     step for V_L (dual variables for lower bounds)
CP   DV_U      O    DP     step for V_U (dual variables for upper bounds)
CP   RESTO     I    INT    <>0: we are in restoration phase
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
CS    DCOPY
CS    DAXPY
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

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      integer NIND
      integer M
      double precision X(N)
      double precision YPY(N)
      double precision PZ(NIND)
      double precision ZPZ(M)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision BNDS_L(NLB)
      double precision BNDS_U(NUB)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision MU
      double precision DX(N)
      double precision DV_L(NLB)
      double precision DV_U(NUB)
      integer RESTO
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Get DX = Y*PY + Z*PZ
C
      if( NIND.gt.0 ) then
         call DCOPY(NIND, PZ, 1, DX(M+1), 1)
      endif
      call DCOPY(M, ZPZ, 1, DX, 1)
      call DAXPY(N, 1.d0, YPY, 1, DX, 1)
C
      if( RESTO.eq.0 ) then
CTODO is this necessary?
C
C     Get DV_L
C
         do i = 1, NLB
            DV_L(i) = MU/S_L(i) - SIGMA_L(i)*DX(ILB(i)) - V_L(i)
         enddo
C
C     Get DV_U
C
         do i = 1, NUB
            DV_U(i) = MU/S_U(i) + SIGMA_U(i)*DX(IUB(i)) - V_U(i)
         enddo
      else
         call DCOPY(NLB, 0d0, 0, DV_L, 1)
         call DCOPY(NUB, 0d0, 0, DV_U, 1)
      endif

      return
      end
