C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      SUBROUTINE AUG_LAG(ITER, N, NIND, M, X, IVAR, NLB, ILB, NUB,
     1     IUB, BNDS_L, BNDS_U, DX, DV_L, DV_U, S_L, S_U,
     2     V_L, V_U, NORIG, XORIG, CSCALE, MU, LAM, LAMOLD, G, F, C,
     3     CNRM0, ALPHA, ALPHA_DUAL, LS_COUNT, NU_OUT, KCONSTR,
     5     LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR,
     2     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     6     DAT, IDAT)
C
C*******************************************************************************
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Line search based on the method given in Cuthrell and Biegler.
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
CP   ITER      I    INT    -1: initialize memory
CP                          0: first iteration
CP                         >0: otherwise
CP   N         I    INT    number of variables (without fixed)
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of constraints
CP   X        I/O   DP     actual iterate (reordered without fixed vars:
CP                             first M entries belong to dependent
CP                             variables, remaining to independent variables)
CP                            I: old point
CP                            O: point after line search
CP   IVAR      I    INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP                            X(i) corresponds to XORIG(IVAR(i))
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. S_L(i) is slack for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. S_U(i) is slack for X(IUB(i)) )
CP   BNDS_L    I    DP     values of lower bounds (ordered as S_L)
CP   BNDS_U    I    DP     values of upper bounds (ordered as S_U)
CP   DX        I    DP     step for X (primal)
CP   DV_L      I    DP     step for V_L (dual variables for lower bounds)
CP   DV_U      I    DP     step for V_U (dual variables for upper bounds)
CP   S_L      I/O   DP     slacks to lower bounds
CP                            I: for start of line search
CP                            O: after line search
CP   S_U      I/O   DP     slacks to upper bounds
CP                            I: for start of line search
CP                            O: after line search
CP   V_L      I/O   DP     dual variables for lower bounds
CP                            I: for start of line search
CP                            O: after line search
CP   V_U      I/O   DP     dual variables for upper bounds
CP                            I: for start of line search
CP                            O: after line search
CP   NORIG     I    INT    number of all variables including fixed vars
CP   XORIG    I/O   DP     actual iterate
CP                            XORIG is ordered in ORIGINAL order (i.e. not
CP                            partitioned into independent and dependent
CP                            variables) (on output: as X)
CP   CSCALE    I    DP     scaling factors for constraints
CP   MU        I    DP     barrier parameter
CP   LAM      I/O   DP     multipliers for equality constraints
CP   LAMOLD    I    DP     multipliers for equality constraints from last iter
CP   G         I    DP     gradient of objective function
CP   F        I/O   DP     value of objective function at X
CP                            I: for start of line search
CP                            O: after line search
CP   C        I/O   DP     values of constraints at X
CP                            I: for start of line search
CP                            O: after line search
CP   CNRM0     O    DP     2-norm of constraints at old point
CP   ALPHA    I/O   DP     step size: I: where to start line search
CP                                    O: step size from X to X_NEW
CP   ALPHA_DUAL I/O DP     step size for dual variables
CP   LS_COUNT  O    INT    number of trial steps
CP   NU_OUT    O    DP     actual value of penalty parameter
CP                            (only output; value is stored internally!)
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
CP   LRW       I    INT    length of RW
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW       I    INT    length of IW
CP   IW       I/O   INT    can be used as INT work space but content will be
CP                            changed between calls
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
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
CS    D1MACH
CS    DDOT
CS    DCOPY
CS    DAXPY
CS    DSCAL
CS    C_OUT
CS    CALC_BAR
CS    GET_F
CS    GET_C
CS    FFINITE
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
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
      integer ITER
      integer N
      integer NIND
      integer M
      double precision X(N)
      integer IVAR(N)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision BNDS_L(NLB)
      double precision BNDS_U(NUB)
      double precision DX(N)
      double precision DV_L(NLB)
      double precision DV_U(NUB)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      integer NORIG
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      double precision MU
      double precision LAM(M)
      double precision LAMOLD(M)
      double precision G(N)
      double precision F
      double precision C(M)
      double precision CNRM0
      double precision ALPHA
      double precision ALPHA_DUAL
      integer LS_COUNT
      double precision NU_OUT
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
      integer p_rwend, p_iwend, p_gb, p_xnew, p_slnew, p_sunew
      integer p_cnew, p_lamnew
      logical ls_ok
      integer k, i
      double precision enrm0, u0c0, ugvh0, dgx, etop, etadd, alag0
      double precision delal0, phi0, enrm, uc, ugvh, phi_new, zlhs, zrhs
      double precision etals, alag, f_new, eta, macheps, precfac
      double precision lhs, rhs, machtiny, dummy
      character*150 line

      double precision DDOT, CALC_BAR, D1MACH
      integer FFINITE

      double precision f_0, bar_0, cnrm_0, f_n, bar_n, cnrm_n

      double precision ETA_OLD
      save             ETA_OLD

      double precision ALPHA_MIN, DELT
CORIG      parameter( ALPHA_MIN = 1d-14, DELT = 1d-1 )
      parameter( ALPHA_MIN = 1d-14, DELT = 1d-4 )
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
      p_rwend = 0
      p_iwend = 0

      precfac = 10d0
      macheps = d1mach(4)
      machtiny = d1mach(1)
C
C     The following is taken form rSQP implementation of aug_lagrangian
C
******************************************************************************
*                                                                            *
*      Calculate ENRM0 and UGVH0 for ETADD calculation.                      *
*                                                                            *
*        ENRM0  = || C(X0) ||**2                                             *
*        UGVH0  = [V(X0+D) - 2.0*V(X0)]*C(X0)                                *
*               = [ UNEXT  - 2.0*UOLD ]*C0                                   *
*                                                                            *
******************************************************************************

      enrm0 = DDOT(M, C, 1, C, 1)
      u0c0  = DDOT(M, LAMOLD, 1, C, 1)
      ugvh0 = DDOT(M, LAM, 1, C, 1)
      ugvh0 = ugvh0 - 2.0d0*u0c0

*  Calculate  ETADD = (DelF*d + UGVH0) / ENRM0
*  Then calculate ETA.  NFLAG follows how ETA is set.

C
C     Compute gradient of barrier function
C
      p_gb    = p_rwend
      p_rwend = p_gb + N
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call DCOPY(N, G, 1, RW(p_gb+1), 1)
      do i = 1, NLB
         k = ILB(i)
         RW(p_gb+k) = RW(p_gb+k) - MU/S_L(i)
      enddo
      do i = 1, NUB
         k = IUB(i)
         RW(p_gb+k) = RW(p_gb+k) + MU/S_U(i)
      enddo
C
      dgx = DDOT(N, DX, 1, RW(p_gb+1), 1)
      etop = dgx + ugvh0
      etadd = 0.d0
      if( M.gt.0 ) then
CTODO Think about this constant 1d-15 !
         if( enrm0.gt.1.D-15 ) etadd = etop/enrm0
      endif
CTODO  think about this:
C      if( ITER.eq.0 ) eta = DABS(etadd)
C      if( ITER.gt.0 ) eta = DMAX1(0.0D0,etadd)
      if( ITER.le.0 ) then
         eta = DABS(etadd) + QRHO
         ETA_OLD = 0d0
      else
         eta = DMAX1(0.0D0,etadd) + QRHO
CTODO Uncomment this if monoton ETA !
C         eta = DMAX1(eta,ETA_OLD)
         ETA_OLD = eta
      endif

      p_rwend = p_gb

*  Calculate ALAG0, the augmented Lagrangian at the initial point, as
*  well as DELAL0.  DELAL0 is an approxomation to the p-directional
*  derivative of the intial augmented Lagrangian where p is the vector
*  containing the update variables d, u(k+1)-u(k) and v(k+1)-v(k). DELAL0
*  is used for the Armijo chord test

      f_0 = F
      bar_0 = CALC_BAR(NLB, NUB, S_L, S_U, V_L, V_U, MU)
      phi0   = F + bar_0
      cnrm_0 = dsqrt(enrm0)

      alag0  = u0c0 + 0.5d0*eta*enrm0 + phi0
      delal0 = ugvh0 - eta*enrm0 + dgx
C      delal0 = ugvh0 - eta*enrm0 + dgx -
C     1     precfac*macheps*dmax1(dabs(ugvh0),dabs(eta*enrm0),dabs(dgx))

CC
CC     Check if indeed descent direction:
CC
C      if( delal0.ge.0d0 ) then
C         write(line,*) 'Problem in aug_lag: delal0 = ',delal0
C         call C_OUT(2,0,1,line)
C         IERR = 589
C         goto 9999
C      endif

C
C     reserve work space for trial points
C
      p_xnew    = p_rwend
      p_slnew   = p_xnew   + N
      p_sunew   = p_slnew  + NLB
      p_cnew    = p_sunew  + NUB
      p_lamnew  = p_cnew   + M
      p_rwend   = p_lamnew + M
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
C     Now do the actual line search
C
      ls_ok = .false.
      LS_COUNT = 0
C
C     While Alpha not too small and not ok try new point
C
 10   continue
      if( ALPHA.le.ALPHA_MIN .or. ls_ok ) goto 20

         LS_COUNT = LS_COUNT + 1

C
C     Compute new trial point and evaluate objective and constraint
C     function at that point
C
         call DCOPY(N, X, 1, RW(p_xnew+1), 1)
         call DAXPY(N, ALPHA, DX, 1, RW(p_xnew+1), 1)
C
C
C     Compute new trial slack variables
C
         call CHECK_SLACKS(N, RW(p_xnew+1), NLB, ILB, S_L, V_L, BNDS_L,
     1        RW(p_slnew+1), NUB, IUB, S_U, V_U, BNDS_U, RW(p_sunew+1),
     1        MU, MACHEPS)
C
         call GET_F(N, RW(p_xnew+1), IVAR, NORIG, XORIG, M, CSCALE, NLB,
     1        ILB, RW(p_slnew+1), NUB, IUB, RW(p_sunew+1), MU, f_new,
     2        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'aug_lag: get_f returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         if( FFINITE(f_new).eq.0 ) then ! we don't want NaN here; cut step!
            lhs = f_new
            goto 800
         endif
C
         call GET_C(ITER, N, NIND, RW(p_xnew+1), IVAR, NORIG, XORIG, M,
     1        CSCALE, RW(p_cnew+1), KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5        EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'aug_lag: Error: get_c returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'aug_lag: Warning: get_c returns IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
C
C     Compute new multipliers
C
         call DCOPY(M, LAM, 1, RW(p_lamnew+1), 1)
         call DSCAL(M, ALPHA, RW(p_lamnew+1), 1)
         call DAXPY(M, (1.d0-ALPHA), LAMOLD, 1,
     1        RW(p_lamnew+1), 1)

* Calculate ENRM and UGVH for ETALS calculation
*
*   ENRM  = || C(X) ||**2
*   UGVH  = U(X)*C(X) - U(X0)*C(X0)   = U*C - U0*C0

         enrm = DDOT(M, RW(p_cnew+1), 1, RW(p_cnew+1)  , 1)
         uc   = DDOT(M, RW(p_cnew+1), 1, RW(p_lamnew+1), 1)
         ugvh = uc - u0c0

*  ETALS is the bound on ETA for which the Armijo test will be satisfied.
*  It can be either an upper or a lower bound and may be positive or
*  negative.  ETALS is defined as
*    ETALS  =  ZLHS / ZRHS
*  where ZLHS and ZRHS are self-explanatory as calculated below.  If ZRHS
*  is very small, ETALS is not calculated.  Instead, whatever ETA already
*  exists is used for the line search.

         f_n = f_new
         bar_n = CALC_BAR(NLB, NUB, RW(p_slnew+1),
     1        RW(p_sunew+1), dummy, dummy, MU)
         phi_new = f_new + bar_n
         cnrm_n = dsqrt(enrm)

         zlhs = phi_new - phi0 + ugvh - DELT*ALPHA*(dgx+ugvh0)
         zrhs = 0.5d0*(enrm0-enrm) - DELT*ALPHA*enrm0

CTODO  CHECK IF THIS IS GOOD...(?)
C         goto 43
         if( M.gt.0 ) then
CTODO check these small values!
            if( dabs(zrhs).le.1.D-15 .and. zlhs.le.1.D-15 ) then
               ls_ok = .true.
               goto 20
            endif
            if( zrhs.ge.0d0 ) then
               ls_ok = .true.
               goto 20
            endif
            etals = zlhs/zrhs
            if( ITER.gt.0 .and. etals.gt.eta ) then
               ls_ok = .true.
               goto 20
            endif
         endif
 43      continue

*  Calculate the augmented Lagrangian at the trial point

         alag = phi_new + uc + 0.5d0*eta*enrm

*  Armijo test

         lhs = alag-alag0-precfac*macheps*abs(alag0)
         rhs = DELT*ALPHA*delal0

         if( QCNR.gt.0 .and. QPRINT.ge.6) then
            write(line,710) LS_COUNT, ALPHA, lhs, rhs
 710        format('LS_COUNT =',i3,', ALPHA =',d13.6,', lhs = ',
     1           d25.17,', rhs = ',d25.17)
            call C_OUT(1,6,1,line)
            write(line,711) alag0, alag, eta, delal0
 711        format('alag0=',d20.12,' alag=',d20.12,' eta=',d20.12,
     1           ' delal0=',d20.12)
            call C_OUT(1,7,1,line)
            write(line,712) f_0, bar_0,cnrm_0,phi0
 712        format('f_0 = ',d25.17,' bar_0 = ',d25.17,' cnrm_0 = ',
     1              d25.17,' phi0 = ',d25.17)
            call C_OUT(1,7,1,line)
            write(line,713) f_n, bar_n,cnrm_n,phi_new
 713        format('f_n = ',d25.17,' bar_n = ',d25.17,' cnrm_n = ',
     1           d25.17,' phin = ',d25.17)
            call C_OUT(1,7,1,line)
         endif

 800     continue
         if( FFINITE(lhs).ne.0 ) then

            if( lhs.le.rhs ) then
               ls_ok = .true.
            else
CTODO decide for interpolation scheme (now same as for exact )
C            alphap = (ALPHA*delal0) / (ALPHA*delal0 - (alag-alag0))
C            ALPHA  = ALPHA * dmax1(0.1D0, dmin1(0.5D0*alphap, 0.9d0))
C            ALPHA = dmin1(6.d-1*ALPHA,
C     1                     dmax1(4.d-1*ALPHA,
C     2          -5.d-1*delal0*(ALPHA**2)/(alag0 - alag -ALPHA*delal0)))
CORIG               ALPHA = dmin1(6.d-1*ALPHA,
CORIG     1              dmax1(4.d-1*ALPHA,
CORIG     2           -5.d-1*delal0*(ALPHA**2)/(alag0 - alag -ALPHA*delal0)))
c               IF(LS_COUNT.eq.1)WRITE(*,*) 'aug_lag: ALPHA=0.5*ALPHA!'
               ALPHA = 0.5d0*ALPHA
            endif

         else
C
C     NaN encountered: Shorten step length
C
            call C_OUT(2,0,1,' --- WARNING: Evaluation of'//
     1           ' some function (obj/contr) produced Nan/Inf.')
            
            ALPHA = 5.d-1*ALPHA

         endif
C
C     end of while loop
C
         goto 10

 20   continue

C
C     line search failure?
C
      if( .not.ls_ok ) then
         write(line,*) 'aug_lag: line search failure'
         call C_OUT(2,0,1,line)
         IERR = 2
         goto 9999
      endif
C
C     Copy trial point to real point
C
      F = f_new
      call DCOPY(N,   RW(p_xnew +1), 1, X  , 1)
      call DCOPY(M,   RW(p_cnew +1), 1, C  , 1)
      call DCOPY(NLB, RW(p_slnew+1), 1, S_L, 1)
      call DCOPY(NUB, RW(p_sunew+1), 1, S_U, 1)
      call DCOPY(M  , RW(p_lamnew+1),1, LAMOLD, 1)
      call DCOPY(M  , RW(p_lamnew+1),1, LAM, 1)
C
C     Do step in dual variables
C
      if( QALPHA.eq.0 ) then
         ALPHA_DUAL = ALPHA
      elseif( QALPHA.eq.1 ) then
         ALPHA_DUAL = min( ALPHA_DUAL, ALPHA )
      endif

      call DAXPY(NLB, ALPHA_DUAL, DV_L, 1, V_L, 1)
      call DAXPY(NUB, ALPHA_DUAL, DV_U, 1, V_U, 1)
C
C     Make sure that each V_NEW is at least machtiny
C
      do i = 1, NLB
         V_L(i) = dmax1(machtiny, V_L(i))
      enddo
      do i = 1, NUB
         V_U(i) = dmax1(machtiny, V_U(i))
      enddo
C
C     Free work space
C
      p_rwend = p_xnew
C
C     For output:
C
      NU_OUT = eta
      CNRM0 = dsqrt(enrm0)
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

      subroutine AUG_LAG_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1, lrw2, liw2

      LRW = 0
      LIW = 0

      LRW = max(LRW, N)         ! p_gb

      call GET_C_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
      call GET_F_WS(N, M, NLB, NUB, NZA, lrw2, liw2)
      liw1 = max(liw1, liw2)
      lrw1 = max(lrw1, lrw2)
      LRW = max(LRW, N+NLB+NUB+2*M+lrw1) ! p_xnew ... p_lamnew GET_C
      LIW = max(LIW, liw1)

      return
      end
