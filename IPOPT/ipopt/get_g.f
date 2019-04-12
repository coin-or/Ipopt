C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_G(N, X, IVAR, NORIG, XORIG, M, CSCALE, NLB, ILB,
     1     S_L, NUB, IUB, S_U, MU, G, LRW, RW, LIW, IW, IERR, EVAL_G,
     1     DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: get_g.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Evaluate gradient of constraint functions at X (after copying
CT    into XORIG) and eliminate fixed variables
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
C
CP   N         I    INT    number of free variables
CP   X         I    DP     values of free variables, where F should be
CP                            evaluated
CP   IVAR      I    INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP                              X(i) corresponds to XORIG(IVAR(i))
CP   NORIG     I    INT    number of all variables including fixed vars
CP   XORIG    I/O   DP     I: values of fixed variables
CP                         O: values of all variables (from X)
CP   M         I    INT    number of equality constraints
CP   CSCALE    I    DP     scaling values (see get_scale)
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. BNDS_L(i) is bound for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. BNDS_U(i) is bound for X(IUB(i)) )
CP   S_L       I    DP     slacks of lower bounds
CP   S_U       I    DP     slacks of upper bounds
CP   MU        I    DP     value of barrier parameter
CP   G         O    DP     gradient of objective function
CP                            (w/o fixed vars, partitioned into ind. and
CP                             dep. varoables)
CP   LRW      I/O   INT    length of RW
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW       I    INT    length of IW
CP   IW       I/O   INT    can be used as INT work space but content will be
CP                            changed between calls
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
CP   EVAL_G    I    EXT    Subroutine for gradient of objective function
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
CS    DSCAL
CS    EVAL_G
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
      integer IVAR(N)
      integer NORIG
      double precision XORIG(NORIG)
      integer M
      double precision CSCALE(*)
      integer NLB
      integer ILB(NLB)
      double precision S_L(NLB)
      integer NUB
      integer IUB(NUB)
      double precision S_U(NUB)
      double precision MU
      double precision G(N)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR
      external EVAL_G
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i, j, p_rwend, p_gorigtmp, p_iwend, p_onebnd
      integer IDAMAX

      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      p_rwend = 0
      p_iwend = 0
      IERR = 0
C
C     Copy new values to XORIG
C
      if( QSCALE.lt.3 ) then
         do i = 1, N
            XORIG(IVAR(i)) = X(i)
         enddo
      else
         do i = 1, N
            XORIG(IVAR(i)) = X(i)*CSCALE(M+i)
         enddo
      endif
C
C     Evaluate gradient of objective function (incl. fixed vars)
C
      p_gorigtmp = p_rwend
      p_rwend    = p_gorigtmp + NORIG
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call EVAL_G(NORIG, XORIG, RW(p_gorigtmp+1), DAT, IDAT)
C
C     reorder accoring to partition
C
      if( QSCALE.lt.3 ) then
         do i = 1, N
            G(i) = RW(p_gorigtmp+IVAR(i))
         enddo
      else
         do i = 1, N
            G(i) = RW(p_gorigtmp+IVAR(i))*CSCALE(M+i)
         enddo
      endif
      p_rwend = p_gorigtmp
C
C     Scale accoring to QFSCALE and scaling factors
C
      if( QFSCALE.ne.1d0 ) then
         call DSCAL(N, QFSCALE, G, 1)
      endif
C
      if( QCNR.gt.0 .and. QPRINT.ge.4 ) then
         i = IDAMAX(N, G, 1)
         write(line,8000) abs(G(i))
 8000    format('|grad F|_max = ',d8.2)
         call C_OUT(1,0,1,line)
      endif
C
C     Add the linear damping term
C
      if( QLINDAMP.gt.0.d0 .and. MU.gt.0.d0 ) then
         p_onebnd = p_iwend
         p_iwend  = p_onebnd + N
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         do i = 1, N
            IW(p_onebnd+i) = 0
         enddo
         do i = 1, NLB
            IW(p_onebnd+ILB(i)) = -i
         enddo
         do i = 1, NUB
            j = IUB(i)
            if( IW(p_onebnd+j).eq.0 ) then
               IW(p_onebnd+j) = i
            else
               IW(p_onebnd+j) = 0
            endif
         enddo
         do i = 1, N
            j = IW(p_onebnd+i)
            if( j.lt.0 ) then
               G(i) = G(i) + QLINDAMP*(MU**QLINDAMPEXP)
            elseif( IW(p_onebnd+i).gt.0 ) then
               G(i) = G(i) - QLINDAMP*(MU**QLINDAMPEXP)
            endif
         enddo
      endif
C
C     I think that's it...
C
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine GET_G_WS(N, M, NLB, NUB, NZA, LRW, LIW)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW

      LRW = N
      if( QLINDAMP.gt.0.d0 ) then
         LIW = N
      else
         LIW = 0
      endif

      return
      end
