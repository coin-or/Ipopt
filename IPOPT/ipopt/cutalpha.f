C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine CUTALPHA(N, X, DX, DV_L, DV_U, NLB, ILB, NUB, IUB,
     1                    BNDS_L, BNDS_U, S_L, S_U, V_L, V_U, TAU,
     2                    ALPHA, ALPHA_DUAL, CALPHA)
C
C*******************************************************************************
C
C    $Id: cutalpha.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Determine ALPHA, so that overall step still fits into nonnegative
CT    orthant.
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
CP   N         I    INT    number of variables (without fixed)
CP   X         I    DP     actual primal iterate
CP   DX        I    DP     step for X (primal)
CP   DV_L      I    DP     step for V_L (dual variables for lower bounds)
CP   DV_U      I    DP     step for V_U (dual variables for upper bounds)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. S_L(i) is slack for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. S_U(i) is slack for X(IUB(i)) )
CP   BNDS_L    I    DP     values of lower bounds (ordered as S_L)
CP   BNDS_U    I    DP     values of upper bounds (ordered as S_U)
CP   S_L       I    DP     slacks to lower bounds
CP   S_U       I    DP     slacks to upper bounds
CP   V_L       I    DP     dual variables for lower bounds
CP   V_U       I    DP     dual variables for upper bounds
CP   TAU       I    DP     fraction to the boundary parameter
CP   ALPHA     O    DP     maximal step size to ensure that primal variables
CP                           stay sufficiently positive
CP   ALPHA_DUAL  O    DP     maximal step size to ensure that dual variables
CP                           stay sufficiently positive
CP   CALPHA    O    C*1    =' ': ALPHA = 1.d0
CP                         ='p': step restricted due to primal variable
CP                         ='d': step restricted due to dual variable
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
CS    C_OUT
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
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision DX(N)
      double precision DV_L(NLB)
      double precision DV_U(NUB)
      double precision BNDS_L(NLB)
      double precision BNDS_U(NUB)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision TAU
      double precision ALPHA
      double precision ALPHA_DUAL
      character*1 CALPHA
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i
      double precision alphap, alphad, dxi
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Take care of primal variables
C
      alphap = 1.d0
      do i = 1, NLB
         dxi = DX(ILB(i))
         if( dxi.lt.0.d0 ) then
            alphap = dmin1(alphap, TAU*(-S_L(i)/dxi) )
         endif
      enddo
      do i = 1, NUB
         dxi = DX(IUB(i))
         if( dxi.gt.0.d0 ) then
            alphap = dmin1(alphap, TAU*S_U(i)/dxi )
         endif
      enddo
C
C     Take care of dual variables
C
      alphad = 1.d0
      do i = 1, NLB
         if( DV_L(i).lt.0.d0 ) then
            alphad = dmin1(alphad, TAU*(-V_L(i)/DV_L(i)) )
         endif
      enddo
      do i = 1, NUB
         if( DV_U(i).lt.0.d0 ) then
            alphad = dmin1(alphad, TAU*(-V_U(i)/DV_U(i)) )
         endif
      enddo

      if( QALPHA.eq.0 ) then
         ALPHA = dmin1(alphad, alphap)
         ALPHA_DUAL = ALPHA
      elseif( QALPHA.eq.1 .or. QALPHA.eq.2 ) then
         ALPHA = alphap
         ALPHA_DUAL = alphad
      else
         write(line,*) 'cutalpha: INVALID QALPHA = ',QALPHA
      endif

      if( alphad.lt.alphap ) then
         CALPHA = 'd'
      elseif( alphap.lt.1.d0 ) then
         CALPHA = 'p'
      else
         CALPHA = ' '
      endif

      return
      end
