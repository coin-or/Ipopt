C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_ZWZV(TASK, N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, X, CSCALE, NLB, ILB, NUB, IUB, S_L, S_U,
     2     SIGMA_L, SIGMA_U, LAM, VIN, VOUT1, VOUT2,
     1     KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR, EV_F,
     1     EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: get_zwzv.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute VOUT1 = Z^T ( Hessian of Lagrangian + Sigma ) Z * VIN and
CT            VOUT2 = Z^T ( Hessian of Lagrangian ) Z * VIN
CT            VOUT1*= Z^T ( Sigma ) Z * VIN
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
CP   TASK      I    INT    =1: Compute only VOUT1
CP                          2: Compute only VOUT2
CP                          3: Compute VOUT1 and VOUT2
CP                          4: Compute VOUT1* (only with Sigma)
CP   N         I    INT    number of (free) variables
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of equality constraints / dependent variables
CP   ITER      I    INT    iteration counter
CP   IVAR      I    INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP   NORIG     I    INT    total number of all variables (incl. fixed vars)
CP   XORIG     I    DP     (only TASK = 1,2,3): actual iterate
CP                            XORIG is ordered in ORIGINAL order (i.e. not
CP                            partitioned into independent and dependent
CP                            variables)
CP   X         I    DP     actual primal iterate
CP   CSCALE    I    DP     scaling factors for constraints
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. S_L(i) is slack for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. S_U(i) is slack for X(IUB(i)) )
CP   S_L       I    DP     slack variables for lower bounds
CP   S_U       I    DP     slack variables for upper bounds
CP   SIGMA_L   I    DP     primal-dual Hessian of lower bound barrier term
CP                            (NLB diagonal elements only)
CP   SIGMA_U   I    DP     primal-dual Hessian of upper bound barrier term
CP                            (NUB diagonal elements only)
CP   LAM       I    DP     Lagrangian multipliers for constraints
CP   VIN       I    DP     vector to use in pruduct (ordered like X)
CP   VOUT1     O    DP     Z^T (W + Sigma) Z * VIN  (ordered like X)
CP   VOUT2     O    DP     Z^T W Z * VIN (ordered like X)
CP   KCONSTR   I    INT    KCONSTR(1): LRS for CONSTR
CP                         KCONSTR(2): P_LRS for CONSTR
CP                         KCONSTR(3): LIS for CONSTR
CP                         KCONSTR(4): P_LIS for CONSTR
CP                         KCONSTR(5): LRW for CONSTR
CP                         KCONSTR(6): LIW for CONSTR
CP   LRS       I    INT    total length of RS
CP   RS       I/O   DP     DP storage space (all!)
CP   LIS       I    INT    total length of IS
CP   IS       I/O   INT    INT storage space (all!)
CP   LRW      I/O   INT    length of RW
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW      I/O   INT    length of IW
CP   IW       I/O   INT    can be used as INT work space but content will be
CP                            changed between calls
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
CP   EV_F      I    EXT    Subroutine for objective function
CP   EV_C      I    EXT    Subroutine for constraints
CP   EV_G      I    EXT    Subroutine for gradient of objective fuction
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
CS    CONSTR
CS    DSCAL
CS    DAXPY
CS    DCOPY
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
      include 'TIMER.INC'
      include 'IPOPT.INC'
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer TASK
      integer N
      integer NIND
      integer M
      integer ITER
      integer IVAR(N)
      integer NFIX
      integer IFIX(NFIX)
      integer NORIG
      double precision XORIG(NORIG)
      double precision X(N)
      double precision CSCALE(*)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision LAM(M)
      double precision VIN(NIND)
      double precision VOUT1(NIND)
      double precision VOUT2(NIND)
      integer KCONSTR(6)
      integer LRS
      double precision RS(LRS)
      integer LIS
      integer IS(LIS)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
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
      integer p_iwend, p_rwend, p_tmp1, p_tmp2, p_tmp3
      integer idummy, i, k
      double precision times, timef
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
      p_rwend = 0

      if( MEMDBG ) then
         write(line,1)'get_zwzv', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      p_tmp1  = p_rwend
      p_tmp2  = p_tmp1  + N
      p_tmp3  = p_tmp2  + N
      p_rwend = p_tmp3  + N
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
C     Compute tmp1 = Z VIN
C
      call TIMER(times)
      if( M.gt.0 ) then
         call CONSTR(5, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, VIN, RW(p_tmp1+1),
     2        idummy, idummy, KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4        IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.lt.0 ) then
            write(line,*) 'get_zwzv: Warning in CONSTR-5, IERR = ',IERR
            call C_OUT(2,0,1,line)
         elseif( IERR.ne.0 ) then
            write(line,*) 'get_zwzv: Error in CONSTR-5, IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         call DSCAL(M, -1d0, RW(p_tmp1+1), 1)
      endif
      call DCOPY(NIND, VIN, 1, RW(p_tmp1+M+1), 1)
      call TIMER(timef)
      TIME_ZWZY_BACKS = TIME_ZWZY_BACKS + timef - times
C
C     Compute tmp2 = W * tmp1
C
      if( TASK.ne.4 ) then
         call TIMER(times)
         call GET_HV(ITER, N, NIND, NFIX, IFIX, X, IVAR, NORIG,
     1        XORIG, NLB, ILB, NUB, IUB, S_L, S_U, CSCALE, M, LAM,
     2        RW(p_tmp1+1), RW(p_tmp2+1),
     2        KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     5        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.lt.0 ) then
            write(line,*) 'get_zwzv: Warning in get_hv, IERR = ',IERR
            call C_OUT(2,0,1,line)
         elseif( IERR.ne.0 ) then
            write(line,*) 'get_zwzv: Error in get_hv, IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         call TIMER(timef)
         TIME_ZWZY_EVALA = TIME_ZWZY_EVALA + timef - times
      endif
C
C     Compute tmp3 = tmp2 + Sigma * tmp1
C
      if( TASK.ne.2 ) then
         if( TASK.eq.4 ) then
            call DCOPY(N, 0.d0, 0, RW(p_tmp3+1), 1)
         else
            call DCOPY(N, RW(p_tmp2+1), 1, RW(p_tmp3+1), 1)
         endif
         do i = 1, NLB
            k = ILB(i)
            RW(p_tmp3+k) = RW(p_tmp3+k) + SIGMA_L(i)*RW(p_tmp1+k)
         enddo
         do i = 1, NUB
            k = IUB(i)
            RW(p_tmp3+k) = RW(p_tmp3+k) + SIGMA_U(i)*RW(p_tmp1+k)
         enddo
      endif
C
C     Compute VOUT2 = Z^T * tmp2
C
      call TIMER(times)
      if( TASK.eq.2 .or. TASK.eq.3 ) then
         if( M.gt.0 ) then
            call CONSTR(4, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, RW(p_tmp2+1), VOUT2,
     2           idummy, idummy,
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.lt.0 ) then
               write(line,*)
     1              'get_zwzv: Warning in CONSTR-4, IERR = ',IERR
               call C_OUT(2,0,1,line)
            elseif( IERR.ne.0 ) then
               write(line,*)
     1              'get_zwzv: Error in CONSTR-4, IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            call DSCAL(NIND, -1d0, VOUT2, 1)
            call DAXPY(NIND, 1d0, RW(p_tmp2+M+1), 1, VOUT2, 1)
         else
            call DCOPY(NIND, RW(p_tmp2+1), 1, VOUT2, 1)
         endif
      endif
C
C     Finally, compute VOUT1 = Z^T * tmp3
C
      if( TASK.eq.1 .or. TASK.ge.3 ) then
         if( M.gt.0 ) then
            call CONSTR(4, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, RW(p_tmp3+1), VOUT1,
     2           idummy, idummy,
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.lt.0 ) then
               write(line,*)
     1              'get_zwzv: Warning in CONSTR-4, IERR = ',IERR
               call C_OUT(2,0,1,line)
            elseif( IERR.ne.0 ) then
               write(line,*)
     1              'get_zwzv: Error in CONSTR-4, IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            call DSCAL(NIND, -1d0, VOUT1, 1)
            call DAXPY(NIND, 1d0, RW(p_tmp3+M+1), 1, VOUT1, 1)
         else
            call DCOPY(NIND, RW(p_tmp3+1), 1, VOUT1, 1)
         endif
      endif
      call TIMER(timef)
      TIME_ZWZY_BACKS = TIME_ZWZY_BACKS + timef - times
C
C     That's it
C
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine GET_ZWZV_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer liw1, lrw1
      character*80 line

      call GET_HV_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      if( M.gt.0 ) then
         call CONSTR_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)
      endif

      LRW = 3*N + LRW

      if( QPRINT.ge.4 ) then
         write(line,1000)'get_zwzv', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end

