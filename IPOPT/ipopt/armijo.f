C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine ARMIJO(ITER, N, M, NIND, X, IVAR, NLB, ILB, NUB, IUB,
     1     BNDS_L, BNDS_U, DX, DV_L, DV_U,
     1     S_L, S_U, V_L, V_U, NORIG, XORIG, CSCALE,
     2     TAKEFULL, ALPHA_MIN, ALPHA, PHI0, DPHI, NU, NUS, MU,
     3     X_NEW, S_LNEW, S_UNEW, V_LNEW, V_UNEW, F_NEW,
     4     C_NEW, PHI_NEW, PHI_NEW_B, CNRM_NEW, LS_COUNT,
     1     KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW,
     1     IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     2     EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: armijo.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Find ALPHA satisfying armijo line search condition
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
CP   ITER      I    INT    iteration counter
CP   N         I    INT    number of variables (without fixed)
CP   M         I    INT    number of constraints
CP   NIND      I    INT    number of independent variables
CP   X         I    DP     actual iterate (reordered without fixed vars:
CP                             first M entries belong to dependent
CP                             variables, remaining to independent variables)
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
CP   S_L       I    DP     slacks to lower bounds
CP   S_U       I    DP     slacks to upper bounds
CP   V_L       I    DP     dual variables for lower bounds
CP   V_U       I    DP     dual variables for upper bounds
CP   NORIG     I    INT    number of all variables including fixed vars
CP   XORIG    I/O   DP     actual iterate
CP                            XORIG is ordered in ORIGINAL order (i.e. not
CP                            partitioned into independent and dependent
CP                            variables)
CP   CSCALE    I    DP     scaling factors for constraints
CP   TAKEFULL I/O   LOG    Input:
CP                           .true. : accept full step (e.g. in watchdog)
CP                           .false.: do armijo line search
CP                         Output:
CP                           .true. : step really accepted
CP                           .false.: quit ls, but Armijo cond. not satisfied
CP   ALPHA_MIN I    DP     I: minimal step size (otherwise ls failure)
CP   ALPHA    I/O   DP     I: stepsize from which line search should start
CP                            (after CUT_ALPHA)
CP                         O: total step size to take after line search
CP                            (if <0: line search failure!)
CP   PHI0      I    DP     value of merit function at X
CP   DPHI      I    DP     value of directed derivative of merit function
CP                            at X
CP                         (For KNITRO-DIRECT update: G^D)
CP   NU        I    DP     penalty parameter (if only one)
CP   NUS       I    DP     individual penalty parameters (QMERIT<0)
CP                         (For KNITRO-DIRECT update: DHD from get_step_corr
CP                           in NUS(1) and CNRM0 in NUS(2))
CP   MU        I    DP     barrier parameter
CP   X_NEW     O    DP     accepted new primal iterate
CP   S_LNEW    O    DP     slack for lower bounds at accepted new iterate
CP   S_UNEW    O    DP     slack for upper bounds at accepted new iterate
CP   V_LNEW    O    DP     accepted new dual iterated for lower bounds
CP   V_UNEW    O    DP     accepted new dual iterated for upper bounds
CP   F_NEW     O    DP     value of objective function at X_NEW
CP   C_NEW     O    DP     values of constraints at X_NEW
CP   PHI_NEW   O    DP     value of merit function at X_NEW (incl NU)
CP   PHI_NEW_B O    DP     value of barrier function at X_NEW
CP   CNRM_NEW  O    DP     norm of constraints at X_NEW
CP   LS_COUNT  O    INT    number of trial steps
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
CS    DCOPY
CS    DAXPY
CS    D1MACH
CS    C_OUT
CS    CALC_BAR
CS    CALC_NRM
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
      integer M
      integer NIND
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
      logical TAKEFULL
      double precision ALPHA_MIN
      double precision ALPHA
      double precision PHI0
      double precision DPHI
      double precision NU
      double precision NUS(*)
      double precision MU
      double precision X_NEW(N)
      double precision S_LNEW(NLB)
      double precision S_UNEW(NUB)
      double precision V_LNEW(NLB)
      double precision V_UNEW(NUB)
      double precision F_NEW
      double precision C_NEW(M)
      double precision PHI_NEW
      double precision PHI_NEW_B
      double precision CNRM_NEW
      integer LS_COUNT
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
      double precision CALC_BAR, CALC_NRM
      double precision cnrm, rhs, lhs, sigma
      integer i, p_rwend, p_iwend
      logical ls_ok

      integer FFINITE

      double precision macheps, machtiny, d1mach
      character*150 line(4)

      double precision ETA_F
      parameter( ETA_F = 1.d-4 )
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

      ls_ok = .false.
      LS_COUNT = 0
C
C     Determine MACHEPS and "precision factor" for armijo decrease condition
C     MAKE INPUT!
C
      macheps = d1mach(4)
      machtiny = d1mach(1)
C
      if( QPRINT.ge.7 ) then
         write(line,700) ITER
 700     format(/,'  Protocol of Armijo line search in ITER',
     1        i5,':'//)
         call C_OUT(1,7,4,line)
      endif
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
         call DCOPY(N, X, 1, X_NEW, 1)
         call DAXPY(N, ALPHA, DX, 1, X_NEW, 1)
C
C     Compute new trial slack variables
C
         call CHECK_SLACKS(N, X_NEW, NLB, ILB, S_L, V_L, BNDS_L,
     1        S_LNEW, NUB, IUB, S_U, V_U, BNDS_U, S_UNEW,
     1        MU, MACHEPS)
C
         call GET_F(N, X_NEW, IVAR, NORIG, XORIG, M, CSCALE, NLB,
     1        ILB, S_LNEW, NUB, IUB, S_UNEW, MU, F_NEW,
     2        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'armijo: get_f returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         if( FFINITE(F_NEW).eq.0 ) then ! we don't want NaN here; cut step!
            lhs = F_NEW
            goto 800
         endif
C
         call GET_C(ITER, N, NIND, X_NEW, IVAR, NORIG, XORIG, M,
     1        CSCALE, C_NEW, KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     3        IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     2        EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'armijo: Error: get_c returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'armijo: Warning: get_c returns IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif

         if( abs(QMERIT).eq.1 ) then
C
C     Compute new trial dual variables
C
            call DCOPY(NLB, V_L, 1, V_LNEW, 1)
            call DAXPY(NLB, ALPHA, DV_L, 1, V_LNEW, 1)
            call DCOPY(NUB, V_U, 1, V_UNEW, 1)
            call DAXPY(NUB, ALPHA, DV_U, 1, V_UNEW, 1)
C
C     Make sure that each V_NEW is at least macheps
C
            do i = 1, NLB
               V_LNEW(i) = dmax1(machtiny, V_LNEW(i))
            enddo
            do i = 1, NUB
               V_UNEW(i) = dmax1(machtiny, V_UNEW(i))
            enddo
         endif
C
C     Compute value of merit function
C
         cnrm = CALC_NRM(M, C_NEW)
         PHI_NEW_B = F_NEW +
     1           CALC_BAR(NLB, NUB, S_LNEW, S_UNEW, V_LNEW, V_UNEW, MU)
         CNRM_NEW = cnrm
         if( QMERIT.gt.0 ) then
            PHI_NEW = PHI_NEW_B + NU*CNRM_NEW
         else
            PHI_NEW = PHI_NEW_B
            if( cnrm.gt.0d0 ) then
               do i = 1, M
                  PHI_NEW = PHI_NEW + dabs(NUS(i)*C_NEW(i))
               enddo
            endif
         endif
C
C     check armijo condition
C
         lhs = PHI_NEW - PHI0 - QPRECFACT*macheps*abs(PHI0)

 800     continue
         if( FFINITE(lhs).ne.0 ) then

            if( QKNITROLS.eq.0 ) then
               rhs = ETA_F*ALPHA*DPHI
            else
               if( NUS(1).gt.0 ) then
                  sigma = NUS(1)/2.d0
                  if( QKNITROLIN.eq.1 ) sigma = 0.d0
               else
                  sigma = 0.d0
               endif
               rhs = ETA_F*( ALPHA*(DPHI-NU*NUS(2))
     1              + ALPHA**2*sigma )
            endif
            if( QCNR.gt.0 .and. QPRINT.ge.3) then
               write(line,710) LS_COUNT, ALPHA, lhs, rhs
 710           format('LS_COUNT =',i3,', ALPHA =',d13.6,', lhs = ',
     1              d25.17,', rhs = ',d25.17)
               call C_OUT(1,3,1,line)
               write(line,711) PHI_NEW_B, PHI_NEW, NU, cnrm, DPHI
 711           format('PHI_NEW_B=',d20.12,' PHI_NEW=',d20.12,
     1              ' NU=',d20.12, 'cnrm=',d20.12,' DPHI=',d20.12)
               call C_OUT(1,3,1,line)
            endif

            if( TAKEFULL .or. lhs.le.rhs ) then
               TAKEFULL = lhs.le.rhs
               ls_ok = .true.
            else
C
C     do quadratic interpolation of merit function
C
CTODO Find bounds for ALPHA
C           ALPHA = dmin1(9.d-1*ALPHA,
C    1                     dmax1(4.d-1*ALPHA,
C    2          -5.d-1*DPHI*(ALPHA**2)/(PHI0 - PHI_NEW -ALPHA*DPHI)))
CORIG               ALPHA = dmin1(6.d-1*ALPHA,
CORIG     1              dmax1(4.d-1*ALPHA,
CORIG     2            -5.d-1*DPHI*(ALPHA**2)/(PHI0 - PHI_NEW -ALPHA*DPHI)))
c               IF(LS_COUNT.eq.1)WRITE(*,*) 'armijo: ALPHA=0.5*ALPHA'
               ALPHA = 0.5d0*ALPHA
C            ALPHA = dmin1(6.d-1*ALPHA,
C     1                     dmax1(1.d-1*ALPHA,
C     2          -5.d-1*DPHI*(ALPHA**2)/(PHI0 - PHI_NEW -ALPHA*DPHI)))

            endif

         else
C
C     NaN encountered: Shorten step length
C
            call C_OUT(2,0,1,' --- WARNING: Evaluation of'//
     1           ' some function (obj/contr) produced Nan/Inf.')
            ALPHA = 0.5d0*ALPHA

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
         write(line,*) 'armijo: line search failure'
         call C_OUT(2,0,1,line)
         IERR = 2
      endif

 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine ARMIJO_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1

      call GET_C_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)
      call GET_F_WS(N, M, NLB, NUB, NZA, lrw1, liw1)
      LRW = max(LRW, lrw1)
      LIW = max(LIW, liw1)

      return
      end
