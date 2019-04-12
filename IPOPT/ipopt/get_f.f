C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_F(N, X, IVAR, NORIG, XORIG, M, CSCALE, NLB, ILB,
     1     S_L, NUB, IUB, S_U, MU, F, LIW, IW, IERR, EVAL_F, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: get_f.f 574 2004-04-25 22:42:37Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Evaluate objective function at X (after copying into XORIG)
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
CP   F         O    DP     value of objective function
CP   LIW       I    INT    length of IW
CP   IW       I/O   INT    can be used as INT work space but content will be
CP                            changed between calls
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
CP   EVAL_F    I    EXT    Subroutine for objective function
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
      double precision F
      integer LIW
      integer IW(LIW)
      integer IERR
      external EVAL_F
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i, j, p_iwend, p_onebnd
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
      p_iwend = 0

      if( MEMDBG ) then
         write(line,1)'get_f', LIW
 1       format('MEMDBG - ',a20,': LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif
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
C     Evaluate objective function
C
      call EVAL_F(NORIG, XORIG, F, DAT, IDAT)
C
C     Scale according to QFSCALE
C
      F = F*QFSCALE
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
               F = F + QLINDAMP*(MU**QLINDAMPEXP)*S_L(-j)
            elseif( IW(p_onebnd+i).gt.0 ) then
               F = F + QLINDAMP*(MU**QLINDAMPEXP)*S_U(j)
            endif
         enddo
      endif
C
C     I think that's it...
C
      fevals = fevals + 1

 9999 continue
      return
      end

      subroutine UPDATE_FG_MU(TASK, N, NLB, ILB, S_L, NUB, IUB, S_U,
     1     MUOLD, MUNEW, F, G, LIW, IW, IERR)
      implicit none
      include 'IPOPT.INC'
      integer TASK, N, NLB, ILB(NLB)
      double precision S_L(NLB)
      integer NUB, IUB(NUB)
      double precision S_U(NUB)
      double precision MUOLD, MUNEW, F
      double precision G(N)
      integer LIW, IW(LIW)
      integer IERR

      integer i, j, p_iwend, p_onebnd

      IERR = 0
      p_iwend = 0

      if( QLINDAMP.gt.0.d0 .and. MUNEW.ne.MUOLD ) then
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
         if( TASK.eq.0 .or. TASK.eq.1 ) then
            do i = 1, N
               j = IW(p_onebnd+i)
               if( j.lt.0 ) then
                  F = F + QLINDAMP*((MUNEW**QLINDAMPEXP)-
     1                 (MUOLD**QLINDAMPEXP))*S_L(-j)
               elseif( IW(p_onebnd+i).gt.0 ) then
                  F = F + QLINDAMP*((MUNEW**QLINDAMPEXP)-
     1                 (MUOLD**QLINDAMPEXP))*S_U(j)
               endif
            enddo
         endif
         if( TASK.eq.0 .or. TASK.eq.2 ) then
            do i = 1, N
               j = IW(p_onebnd+i)
               if( j.lt.0 ) then
                  G(i) = G(i) + QLINDAMP*((MUNEW**QLINDAMPEXP)-
     1                 (MUOLD**QLINDAMPEXP))
               elseif( IW(p_onebnd+i).gt.0 ) then
                  G(i) = G(i) - QLINDAMP*((MUNEW**QLINDAMPEXP)-
     1                 (MUOLD**QLINDAMPEXP))
               endif
            enddo
         endif
      endif
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine GET_F_WS(N, M, NLB, NUB, NZA, LRW, LIW)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      character*80 line

      LRW = 0
      if( QLINDAMP.ne.0 ) then
         LIW = N
      else
         LIW = 0
      endif

      if( QPRINT.ge.4 ) then
         write(line,1000)'get_f_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end

      subroutine UPDATE_FG_MU_WS(N, M, NLB, NUB, NZA, LRW, LIW)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      character*80 line

      LRW = 0
      if( QLINDAMP.lt.0.d0 ) then
         LIW = N
      else
         LIW = 0
      endif

      if( QPRINT.ge.4 ) then
         write(line,1000)'update_fg_mu_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
