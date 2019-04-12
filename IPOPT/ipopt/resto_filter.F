C Copyright (C) 2003, International Business Machines and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine RESTO_FILTER(ITER, N, NIND, M, NORIG, XORIG,
     1     CSCALE, IVAR, NFIX, IFIX, NLB, ILB, NUB,
     1     IUB, BNDS_L, BNDS_U, MU, MUORIG, X, CNRM, NFILTER,
     1     FILTER_PHI, FILTER_C, CNRM_MAX, S_L, S_U, V_L, V_U, F, C,
     1     DX, DV_L, DV_U, SIGMA_L, SIGMA_U, RHS_RESTO, SIGMA_PN,
     1     ETA_RESTO, PRECFACT, MACHEPS, MACHTINY, LAM, LAMOLD, RESTO,
     1     C_ACCEPT, ALPHA, ALPHA_DUAL, LS_COUNT, DELTA,
     1     ALPHA_MIN_FRAC, ETA_F, BACKTRACK_FAC, B, REGU,
     1     KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR, EV_F,
     1     EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: resto_filter.F 664 2004-11-05 18:46:52Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Restoration phase for filter using an interior point approach
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
CA    Andreas Waechter      11/06/02  First version
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
CP   ITER      I    INT    iteration counter (if ITER=-1: initialize)
CP   N         I    INT    number of variables (without fixed)
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of constraints
CP   NORIG    I/O   INT    total number of variables (incl. fixed vars)
CP   XORIG    I/O   INT    actual iterate
CP                            (original order as in problem statement)
CP   IVAR     I/O   INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP                            X(i) corresponds to XORIG(IVAR(i))
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP                            (assumed to be in increasing order)
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. BNDS_L(i) is bound for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. BNDS_U(i) is bound for X(IUB(i)) )
CP   BNDS_L   I/O   DP     values of lower bounds (ordered as S_L)
CP   BNDS_U   I/O   DP     values of upper bounds (ordered as S_U)
CP   MU       I/O   DP     barrier parameter for restoration phase problem
CP   MUORIG   I/O   DP     barrier parameter for original problem
CP   X        I/O   DP     actual iterate (reordered without fixed vars:
CP                             first M entries belong to dependent
CP                             variables, remaining to independent variables)
CP   CNRM      O    DP     norm of constraint violation
CP   NFILTER   I    INT    number of filter entries (corners)
CP                         for ITER = -1: NZORIG
CP   FILTER_PHI I   DP     filter coordinates; phi-axis
CP   FILTER_C  I    DP     filter coordinates; feasibility axis
CP   CNRM_MAX  I    DP     maximal value of CNRM for filter
CP   S_L      I/O   DP     slacks to lower bounds
CP   S_U      I/O   DP     slacks to upper bounds
CP   V_L      I/O   DP     dual variables for lower bounds
CP   V_U      I/O   DP     dual variables for upper bounds
CP   F        I/O   DP     value of objective function at X
CP   C        I/O   DP     values of constraints at X
CP   DX        I    DP     Newton step for dogleg
CP   DV_L      O    DP     Step for dual variables (lower bounds)
CP   DV_U      O    DP     Step for dual variables (upper bounds)
CP   SIGMA_L   I    DP     primal dual barrier Hessian for lower bounds
CP   SIGMA_U   I    DP     primal dual barrier Hessian for upper bounds
CP   RHS_RESTO O    DP     RHS for get_step_full
CP   SIGMA_PN  O    DP     Sigma for P and N variables
CP   ETA_RESTO O    DP     Factor for regularization term
CP   PRECFACT  I    DP     precision factor for roundoffs
CP   DX        I    DP     (Newton) step for X (primal)
CP   MACHEPS   I    DP     machine precision
CP   MACHTINY  I    DP     smallest DP number
CP   LAM       I    DP     current values of multipliers
CP   LAMOLD   I/O   DP     previous values of multipliers
CP   RESTO    I/O   INT    =0: we are not in restoration phase
CP                         =1: we are regularly in the restoration phase
CP                         =2: we are in iterative refinement
CP   C_ACCEPT  O    C*1    flag indication what kind of step was taken
CP   ALPHA     O    DP     fraction of Newton step taken (before projection)
CP   ALPHA_DUAL O   DP     fraction of the Cauchy step taken (before projection)
CP   LS_COUNT I/O   INT    add the number of trial steps taken
CP   DELTA     I    DP     parameter in filter method
CP   ALPHA_MIN_FRAC I DP   parameter in filter method
CP   ETA_F     I    DP     parameter in filter method
CP   BACKTRACK_FAC I DP    parameter in filter method
CP   B        I/O   DP     Quasi-Newton estimate (is restored after
CP                            restoration phase finished)
CP   REGU      I    DP     Regularization parameter (get_step_full)
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
CS    DNRM2
CS    DCOPY
CS    DAXPY
CS    CONSTR
CS    CALC_BAR
CS    CALC_NRM
CS    GET_C
CS    GET_F
CS    FFINITE
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
      include 'TIMER.INC'

#include <config_f.h>
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer ITER
      integer N
      integer NIND
      integer M
      integer NORIG
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      integer IVAR(N)
      integer NFIX
      integer IFIX(NFIX)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision BNDS_L(NLB)
      double precision BNDS_U(NUB)
      double precision MU
      double precision MUORIG
      double precision X(N)
      double precision CNRM
      integer NFILTER
      double precision FILTER_PHI(NFILTER)
      double precision FILTER_C(NFILTER)
      double precision CNRM_MAX
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision F
      double precision C(M)
      double precision DX(N)
      double precision DV_L(NLB)
      double precision DV_U(NUB)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision RHS_RESTO(3*N+2*M)
      double precision SIGMA_PN(M)
      double precision ETA_RESTO
      double precision PRECFACT
      double precision MACHEPS
      double precision MACHTINY
      double precision LAM(M)
      double precision LAMOLD(M)
      integer RESTO
      character*1 C_ACCEPT
      double precision ALPHA
      double precision ALPHA_DUAL
      integer LS_COUNT
      double precision DELTA
      double precision ALPHA_MIN_FRAC
      double precision ETA_F
      double precision BACKTRACK_FAC
      double precision B(*)
      double precision REGU
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
C                            Local variables
C-------------------------------------------------------------------------------
C
      integer p_rwend, p_iwend, p_tmp, p_xnew, p_cnew, p_slnew, p_sunew
      integer p_dp, p_dn, p_dvp, p_dvn, p_pnew, p_nnew, p_check
      integer i, j, k, idummy, irej, irhs
      double precision tmp, resto_f, resto_cnrm, cmax
      double precision cnrm_new, bar_n, phi_new, fi, ci, lhs, rhs
      double precision xstorenorm, dummy, err_totl, scalmin, csum
      double precision err_dinf, err_pinf, err_cmpl, gd, ptmp, ntmp
      double precision alpha_min, resto_f_new, resto_cnrm_new
      double precision model_alpha, f_add, cnrm_add, nrm_ratio
      double precision nrm_ratio_max, nrm_sol
      double precision scal_dinf, scal_cmpl, fnew, nrm_res
      logical ftype, augment, skip_conv_check, first_mu, decmu
      logical mu_changed

      character*140 line(4)
      character*1 cdummy

      double precision DNRM2, CALC_BAR, CALC_NRM, DDOT, DASUM
      integer FFINITE, IDAMAX, FILTER_CHECK

      double precision CNRM_START, RHO_RESTO, NRM_RHS, NRM_MAT
      save             CNRM_START, RHO_RESTO, NRM_RHS, NRM_MAT

      double precision RERR, RESTO_CNRM_MAX, TAUSTORE
      save             RERR, RESTO_CNRM_MAX, TAUSTORE
C      double precision RESTO_CNRM_MAX_FACT
Cc      parameter( RESTO_CNRM_MAX_FACT = 1.d2 )
C      parameter( RESTO_CNRM_MAX_FACT = 1.d2 )

      integer P_XSTORE, P_RFILTERF, P_RFILTERC, NRFILTER
      save    P_XSTORE, P_RFILTERF, P_RFILTERC, NRFILTER
      integer P_P, P_N, P_VP, P_VN, P_BSTORE, P_RHS, P_RHSCOPY
      save    P_P, P_N, P_VP, P_VN, P_BSTORE, P_RHS, P_RHSCOPY
      integer P_DXSTORE, P_DPSTORE, P_DNSTORE, P_LAMSTORE
      save    P_DXSTORE, P_DPSTORE, P_DNSTORE, P_LAMSTORE
      integer P_DVLSTORE, P_DVUSTORE, P_DVPSTORE, P_DVNSTORE
      save    P_DVLSTORE, P_DVUSTORE, P_DVPSTORE, P_DVNSTORE
      integer P_RSCALE
      save    P_RSCALE

      integer IREJ_LAST, COUNT_ITER_BLOCK, ITER_START
      save    IREJ_LAST, COUNT_ITER_BLOCK, ITER_START

      integer COUNT_ITER_THETA_MAX, REFINE_ITER
      save    COUNT_ITER_THETA_MAX, REFINE_ITER

      logical RESTO_INIT_LAM, RESTO_RLAST
      save    RESTO_INIT_LAM, RESTO_RLAST

      double precision ETARESTOFACTINIT, ETARESTOFACT
c      parameter( ETARESTOFACTINIT = 1d-4 )
C#define OLD
#ifdef OLD
      parameter( ETARESTOFACTINIT = 1d-4 )
#else
c      parameter( ETARESTOFACTINIT = 1d-4 )
      parameter( ETARESTOFACTINIT = 1d0 )
#endif
      save ETARESTOFACT

      double precision RHO_INIT
      parameter( RHO_INIT = 1d3 )

      double precision PNINITDIST
      parameter( PNINITDIST = 1.d0 )

      integer QETAUPDATE
      parameter( QETAUPDATE = 3 )

      integer REFINE_ITER_MAX
      parameter( REFINE_ITER_MAX = 5 )

      double precision RESTO_AMAX
      parameter( RESTO_AMAX = 1.d2 )

#ifdef USE_IP_MYRESTOBARRIERTEST
      logical IP_MYRESTOBARRIERTEST
#endif

      integer QITERTHETAMAXMAX
C      parameter( QITERTHETAMAXMAX = 3 )
      parameter( QITERTHETAMAXMAX = 1000000000 )
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

      if( ITER.eq.-1 ) then
         P_XSTORE   = LRS
         P_RFILTERF = P_XSTORE + N
         P_RFILTERC = P_RFILTERF + QMAXITER
         P_P        = P_RFILTERC + QMAXITER
         P_N        = P_P  + M
         P_VP       = P_N  + M
         P_VN       = P_VP + M
         P_RHS      = P_VN + M
         P_RHSCOPY  = P_RHS + 5*M + N + NLB + NUB
         P_DXSTORE  = P_RHSCOPY + 5*M + N + NLB + NUB
         P_DPSTORE  = P_DXSTORE + N
         P_DNSTORE  = P_DPSTORE + M
         P_LAMSTORE = P_DNSTORE + M
         P_DVLSTORE = P_LAMSTORE + M
         P_DVUSTORE = P_DVLSTORE + NLB
         P_DVPSTORE = P_DVUSTORE + NUB
         P_DVNSTORE = P_DVPSTORE + M
         P_RSCALE   = P_DVNSTORE + M
         LRS        = P_RSCALE + M
         if( abs(QQUASI).ge.6 ) then
            P_BSTORE = LRS
            LRS      = P_BSTORE + 2+2*N+QLMLEN*(1+4*N+2*QLMLEN)
         elseif( QQUASI.ne.0 ) then
            write(line,*) 'This restoration phase is not implemented'//
     1           ' for QQUASI = ',QQUASI
            call C_OUT(2,0,1,line)
            IERR = 4
         endif
         goto 9999
      endif

      if( MEMDBG .and. QPRINT.ge.3 ) then
         write(line,1)'resto_filter', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      if( abs(QMERIT).eq.5 ) then
         call C_OUT(2,0,1,'resto_filter: Assuming QMERIT = 4')
      endif

      skip_conv_check = .false.
C
C     Initialize trust region radius and store current iterate
C
      if( RESTO.eq.0 ) then
         ITER_START = ITER
         RESTO = 1
         RESTO_INIT_LAM = .false. ! If this is set to true then the first
                                  ! step in LAM will be undamped
         CNRM_START = CALC_NRM(M, C)
         if( CNRM_START.lt.QTOL ) then
            write(line,*) 'Restoration started with CNRM = ',CNRM_START
            call C_OUT(2,0,1,line)
            IERR = 18
            goto 9999
         endif
         COUNT_RESTO_CALL = COUNT_RESTO_CALL + 1
         if( QSCALE.ne.0 ) then
C
C     Compute scaling for the constraints for the restoration phase
C     Scale it so that the max norm of all gradients is less than
C     RESTO_AMAX
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
            call CONSTR(14, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, dummy, RS(P_RSCALE+1),
     2           idummy, idummy,
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.lt.0 ) then
               write(line,*)
     1              'resto_filter: Warning in CONSTR-14, IERR = ', IERR
               call C_OUT(2,0,1,line)
            elseif( IERR.ne.0 ) then
               write(line,*)
     1              'resto_filter: Error in CONSTR-14, IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            scalmin = 1.d300
            do i = 1, M
               RS(P_RSCALE+i) =
     1              RESTO_AMAX/max(RESTO_AMAX,RS(P_RSCALE+i))
#ifdef OLD
               RS(P_RSCALE+i) = 1.d0
#endif
cCDELETEME
c               RS(P_RSCALE+i) = 1.d0
               scalmin = min(scalmin, RS(P_RSCALE+i))
            enddo
            if( QCNR.gt.0 .and.QPRINT.ge.1 ) then
               write(line,*) 'resto_filter: smallest RSCALE = ',scalmin
               call C_OUT(1,0,1,line)
            endif
            if( QCNR.gt.0 .and. QPRINT.ge.5 ) then
               write(line,8001)
 8001          format(/,
     1'Information about constraint scaling factors restoration phase:',
     2/)
               call C_OUT(1,0,3,line)
               do i = 1, M
                  write(line,8002) i,RS(P_RSCALE+i)
 8002             format('Scaling for constraint ',i8,': ',d10.4)
                  call C_OUT(1,0,1,line)
               enddo
            endif
         else
            call DCOPY(M, 1.d0, 0, RS(P_RSCALE+1), 1)
         endif

c         RESTO_CNRM_MAX = RESTO_CNRM_MAX_FACT * CALC_NRM(M, C)
         p_tmp   = p_rwend
         p_rwend = p_tmp + M
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         do i = 1, M
            RW(p_tmp+i) = C(i)*RS(P_RSCALE+i)
         enddo
         csum = CALC_NRM(M, RW(p_tmp+1))
         RESTO_CNRM_MAX = QTHETAMAXFACT * max(1.d0, csum)
         i = IDAMAX(M, RW(p_tmp+1), 1)
         cmax = abs(RW(p_tmp+i))
         p_rwend = p_tmp

         xstorenorm = dmax1(1.d0, DNRM2(N, X, 1))
         do i = 1, N
#ifdef OLD
            RHS_RESTO(N+M+i) = 1.d0
#else
            RHS_RESTO(N+M+i) = 1.d0/((max(1.d0,abs(X(i)))**2))
#endif
         enddo
#ifndef OLD
         xstorenorm = 1.d0
#endif

         MUORIG  = MU
         TAUSTORE = QTAU
C
C     Call UPDATE_MU to get error tolerance RERR
C
CTODO
         i = IDAMAX(M, C, 1)
         MU = max(MU,abs(C(i)))
#ifndef OLD
c         MU = max(QMUINIT, cmax)
c         MU = max(QMUINIT, csum/dble(M))
c         MU = max(MU,DASUM(M,C,1)/dble(M))
c         MU = QMUINIT
#endif
c         MU = QMUINIT
         first_mu = .true.
         call UPDATE_MU(first_mu, RESTO, MU, RERR)

         ETARESTOFACT = ETARESTOFACTINIT
         if( QETAUPDATE.eq.1 ) then
            ETA_RESTO = ETARESTOFACT*MU/xstorenorm
         elseif( QETAUPDATE.eq.2 ) then
            ETA_RESTO = ETARESTOFACT*MU/xstorenorm
         elseif( QETAUPDATE.eq.3 ) then
c            ETA_RESTO = ETARESTOFACT*min(1.d0,MU)/xstorenorm/xstorenorm
c            ETA_RESTO = ETARESTOFACT*MU/xstorenorm/xstorenorm
            ETA_RESTO = ETARESTOFACT*dsqrt(MU)/xstorenorm/xstorenorm
         else
            ETA_RESTO = ETARESTOFACT*dsqrt(MU)/xstorenorm
         endif

         if( NIND.eq.0 ) ETA_RESTO = 0.d0

         RHO_RESTO = RHO_INIT

         if( QPRINT.ge.1 ) then
            write(line,111) MU, ETA_RESTO, RHO_RESTO
 111        format('Initialize:  MU = ',d12.4,' ETA_RESTO = ',d12.4,
     1           ' RHO_RESTO = ',d12.4)
            call C_OUT(1,0,1,line)
         endif

         C_ACCEPT = '-'

cCDELETEME
c         dummy = MU
c         MU = 1d-10
         call RESTO_RESTO(M, C, RS(P_RSCALE+1), RS(P_P+1), RS(P_N+1),
     1        RS(P_VP+1), RS(P_VN+1), LAM, MU, RHO_RESTO)
         do i = 1, M
            RS(P_VP+i) = MU/RS(P_P+i)
            RS(P_VN+i) = MU/RS(P_N+i)
         enddo
CDELETEME
c         MU = dummy
C
C     Store current iterate and barrier parameter
C
         call DCOPY(N, X, 1, RS(P_XSTORE+1), 1)
         if( QPRINT.ge.5 ) then
            write(line,100) ITER
 100        format(/,'  Reference point (XSTORE) for restoration ',
     1           'phase (Iter ',i5,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, N
               write(line,101) i,RS(P_XSTORE+i)
 101           format(' XSTORE(',i6,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif


C     Initialize restoration phase filter
         NRFILTER = 0
         COUNT_ITER_BLOCK = 0
         COUNT_ITER_THETA_MAX = 0

C     If we do L-BFGS, store current estimate
         if( abs(QQUASI).ge.6 ) then
            call DCOPY(2+2*N+QLMLEN*(1+4*N+2*QLMLEN), B, 1,
     1           RS(P_BSTORE+1), 1)
            B(1) = 0.d0
            B(2) = 1.d0
         endif

         RESTO_RLAST = .false.

CTODO
CTRY  see if it helps to bring original multipliers into some order...
         do i = 1, NLB
            V_L(i) = dmin1(V_L(i),RHO_RESTO)
         enddo
         do i = 1, NUB
            V_U(i) = dmin1(V_U(i),RHO_RESTO)
         enddo

CTODO
         call DCOPY(M, 0.d0, 0, LAM, 1)
         call DCOPY(M, 0.d0, 0, LAMOLD, 1)

C
C     Make sure V doesn't differ too much form MU S^{-1} for p and n variables
C
         call CORRECT_V(M, RS(P_P+1), RS(P_VP+1), M, RS(P_N+1),
     1        RS(P_VN+1), MU)


         goto 5000              ! Recompute step (so that we really have
                                !                 a Newton step)
      endif

C ==============================================================================
C
C              Compute search directions (incl. iterative refinement)
C
C ==============================================================================

C     Reserve space for search directions
      p_dp    = p_rwend
      p_dn    = p_dp  + M
      p_dvp   = p_dn  + M
      p_dvn   = p_dvp + M
      p_rwend = p_dvn + M
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif

      do i = 1, M
         RW(p_dp+i) = (RS(P_P+i)*(RS(P_RHSCOPY+N+i)+
     1        LAM(i)/RS(P_RSCALE+i))+
     1        RS(P_RHSCOPY+N+3*M+NLB+NUB+i))/RS(P_VP+i)
      enddo
      do i = 1, M
         RW(p_dn+i) = (RS(P_N+i)*(RS(P_RHSCOPY+N+M+i)-
     1        LAM(i)/RS(P_RSCALE+i))+
     1        RS(P_RHSCOPY+N+4*M+NLB+NUB+i))/RS(P_VN+i)
      enddo
      do i = 1, M
         RW(p_dvp+i) = (RS(P_RHSCOPY+N+3*M+NLB+NUB+i)-
     1        RW(p_dp+i)*RS(P_VP+i))/RS(P_P+i)
      enddo
      do i = 1, M
         RW(p_dvn+i) = (RS(P_RHSCOPY+N+4*M+NLB+NUB+i)-
     1        RW(p_dn+i)*RS(P_VN+i))/RS(P_N+i)
      enddo
      do i = 1, NLB
         j = ILB(i)
         DV_L(i) = (RS(P_RHSCOPY+N+3*M+i)-V_L(i)*DX(j))/S_L(i)
      enddo
      do i = 1, NUB
         j = IUB(i)
         DV_U(i) = (RS(P_RHSCOPY+N+3*M+NLB+i)+V_U(i)*DX(j))/S_U(i)
      enddo
C
C     If in iterative refinement, update directions
C
      if( RESTO.eq.2 ) then
         call DAXPY(N, 1.d0, RS(P_DXSTORE+1), 1, DX, 1)
         call DAXPY(M, 1.d0, RS(P_DPSTORE+1), 1, RW(p_dp+1), 1)
         call DAXPY(M, 1.d0, RS(P_DNSTORE+1), 1, RW(p_dn+1), 1)
         call DAXPY(M, 1.d0, RS(P_LAMSTORE+1), 1, LAM, 1)
         call DAXPY(NLB, 1.d0, RS(P_DVLSTORE+1), 1, DV_L, 1)
         call DAXPY(NUB, 1.d0, RS(P_DVUSTORE+1), 1, DV_U, 1)
         call DAXPY(M, 1.d0, RS(P_DVPSTORE+1), 1, RW(p_dvp+1), 1)
         call DAXPY(M, 1.d0, RS(P_DVNSTORE+1), 1, RW(p_dvn+1), 1)
      endif
C
      p_check = p_rwend
      p_rwend = p_check + N+5*M+NLB+NUB
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
C     If this is the first solution (no iterative refinement yet), compute
C     norm of right hand side and matrix
C     RHS_RESTO(2*N+M+1..3*N+2*M) contains already entries for Hessian and
C     Jacobian
C
      if( RESTO.eq.1 ) then
         i =  IDAMAX(N+5*M+NLB+NUB, RS(P_RHS+1), 1)
         NRM_RHS = abs(RS(P_RHS+i))
         irhs = 0
         call DCOPY(N, RHS_RESTO(2*N+M+1), 1, RW(p_check+irhs+1), 1)
         do i = 1, NLB
            k = ILB(i)
            RW(p_check+irhs+k) = RW(p_check+irhs+k) + V_L(i)
         enddo
         do i = 1, NUB
            k = IUB(i)
            RW(p_check+irhs+k) = RW(p_check+irhs+k) + V_U(i)
         enddo
         irhs = irhs + N
         do i = 1, M
            RW(p_check+irhs+i) = 1/RS(P_RSCALE+i) + RS(P_VP+i)
         enddo
         irhs = irhs + M
         do i = 1, M
            RW(p_check+irhs+i) = 1/RS(P_RSCALE+i) + RS(P_VN+i)
         enddo
         irhs = irhs + M
         call DCOPY(M, RHS_RESTO(3*N+M+1), 1, RW(p_check+irhs+1), 1)
         do i = 1, M
            RW(p_check+irhs+i) = RW(p_check+irhs+i) +
     1           2.d0*RS(P_RSCALE+i)
         enddo
         irhs = irhs + M
         call DCOPY(NLB, S_L, 1, RW(p_check+irhs+1), 1)
         call DAXPY(NLB, 1.d0, 1.d0, 0, RW(p_check+irhs+1), 1)
         irhs = irhs + NLB
         call DCOPY(NUB, S_U, 1, RW(p_check+irhs+1), 1)
         call DAXPY(NUB, 1.d0, 1.d0, 0, RW(p_check+irhs+1), 1)
         irhs = irhs + NUB
         call DCOPY(M, RS(P_P+1), 1, RW(p_check+irhs+1), 1)
         call DAXPY(M, 1.d0, 1.d0, 0, RW(p_check+irhs+1), 1)
         irhs = irhs + M
         call DCOPY(M, RS(P_N+1), 1, RW(p_check+irhs+1), 1)
         call DAXPY(M, 1.d0, 1.d0, 0, RW(p_check+irhs+1), 1)
         irhs = irhs + M
         i = IDAMAX(N+5*M+NLB+NUB, RW(p_check+1), 1)
         NRM_MAT = abs(RW(p_check+i))
      endif
C
C     Compute residual
C     RHS_RESTO = W * DX from get_step_full
C

C     1. dual infeasibility for X
      irhs = 0
      call DCOPY(N, RHS_RESTO, 1, RW(p_check+irhs+1), 1)
      do i = 1, N
         RW(p_check+irhs+i) = RW(p_check+irhs+i) +
     1        ETA_RESTO*RHS_RESTO(N+M+i)*DX(i)
      enddo
      do i = 1, NLB
         j = ILB(i)
         RW(p_check+irhs+j) = RW(p_check+irhs+j) - DV_L(i)
      enddo
      do i = 1, NUB
         j = IUB(i)
         RW(p_check+irhs+j) = RW(p_check+irhs+j) + DV_U(i)
      enddo

      p_tmp   = p_rwend
      p_rwend = p_tmp + N
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call CONSTR(8, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, LAM, RW(p_tmp+1), idummy,
     3     idummy, KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4     IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.lt.0 ) then
         write(line,*)
     1        'resto_filter: Warning in CONSTR-8a, IERR = ',IERR
         call C_OUT(2,0,1,line)
      elseif( IERR.ne.0 ) then
         write(line,*) 'resto_filter: Error in CONSTR-8a, IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      call DAXPY(N, 1.d0, RW(p_tmp+1), 1, RW(p_check+irhs+1), 1)
      p_rwend = p_tmp
      irhs = irhs + N

C     2. dual infeasibility for P
      do i = 1, M
         RW(p_check+irhs+i) = -LAM(i)/RS(P_RSCALE+i) - RW(p_dvp+i)
      enddo
      irhs = irhs + M

C     3. dual infeasibility for N
      do i = 1, M
         RW(p_check+irhs+i) = LAM(i)/RS(P_RSCALE+i) - RW(p_dvn+i)
      enddo
      irhs = irhs + M

C     4. Primal infeasibility
      call CONSTR(9, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, DX, RW(p_check+irhs+1), idummy,
     3     idummy, KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4     IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.lt.0 ) then
         write(line,*) 'resto_filter: Warning in CONSTR-9, IERR = ',IERR
         call C_OUT(2,0,1,line)
      elseif( IERR.ne.0 ) then
         write(line,*) 'resto_filter: Error in CONSTR-9, IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      do i = 1, M
         RW(p_check+irhs+i) = RW(p_check+irhs+i) +
     1        (RW(p_dn+i)-RW(p_dp+i))/RS(P_RSCALE+i)
      enddo
      irhs = irhs + M

C     5. complementarity X_L
      do i = 1, NLB
         j = ILB(i)
         RW(p_check+irhs+i) = V_L(i)*DX(j) + S_L(i)*DV_L(i)
      enddo
      irhs = irhs + NLB

C     6. complementarity X_U
      do i = 1, NUB
         j = IUB(i)
         RW(p_check+irhs+i) = -V_U(i)*DX(j) + S_U(i)*DV_U(i)
      enddo
      irhs = irhs + NUB

C     7. complementarity P
      do i = 1, M
         RW(p_check+irhs+i) = RS(P_VP+i)*RW(p_dp+i) +
     1        RS(P_P+i)*RW(p_dvp+i)
      enddo
      irhs = irhs + M

C     8. complementarity N
      do i = 1, M
         RW(p_check+irhs+i) = RS(P_VN+i)*RW(p_dn+i) +
     1        RS(P_N+i)*RW(p_dvn+i)
      enddo
      irhs = irhs + M

      i = IDAMAX(N, DX, 1)
      nrm_sol = abs(DX(i))
      i = IDAMAX(M, RW(p_dp+1), 1)
      nrm_sol = max(nrm_sol, abs(RW(p_dp+i)))
      i = IDAMAX(M, RW(p_dn+1), 1)
      nrm_sol = max(nrm_sol, abs(RW(p_dn+i)))
      i = IDAMAX(M, LAM, 1)
      nrm_sol = max(nrm_sol, abs(LAM(i)))
      if( NLB.gt.0 ) then
         i = IDAMAX(NLB, DV_L, 1)
         nrm_sol = max(nrm_sol, abs(DV_L(i)))
      endif
      if( NUB.gt.0 ) then
         i = IDAMAX(NUB, DV_U, 1)
         nrm_sol = max(nrm_sol, abs(DV_U(i)))
      endif
      i = IDAMAX(M, RW(p_dvp+1), 1)
      nrm_sol = max(nrm_sol, abs(RW(p_dvp+i)))
      i = IDAMAX(M, RW(p_dvn+1), 1)
      nrm_sol = max(nrm_sol, abs(RW(p_dvn+i)))

CCDELETEME
c      do i = 1, N+5*M+NLB+NUB
c         write(line,*) 'sol ',i,RW(p_check+i)
c         call C_OUT(1,0,1,line)
c      enddo
C     Now compute residual
      call DSCAL(N+5*M+NLB+NUB, -1.d0, RW(p_check+1), 1)
      call DAXPY(N+5*M+NLB+NUB, 1.d0, RS(P_RHS+1), 1, RW(p_check+1), 1)
      i = IDAMAX(N+5*M+NLB+NUB, RW(p_check+1), 1)
      nrm_res = abs(RW(p_check+i))
C      i = IDAMAX(N+5*M+NLB+NUB, RS(P_RHS+1), 1)
C      nrm_rhs = abs(RS(P_RHS+i))
C      if( nrm_rhs.gt.0.d0 ) then
C         nrm_ratio = nrm_res/nrm_rhs
CCTODO         nrm_ratio = nrm_res/max(1.d0,nrm_rhs)
C      else
C         nrm_ratio = 0.d0
C      endif
      if( NRM_MAT*nrm_sol+NRM_RHS.ne.0.d0 ) then
         nrm_ratio = nrm_res/(NRM_MAT*nrm_sol+NRM_RHS)
      else
         nrm_ratio = 0.d0
      endif

      if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
         write(line,1105) NRM_RHS, nrm_sol, nrm_res, NRM_MAT, nrm_ratio,
     1        REFINE_ITER
 1105    format('Restoration phase system solution: |rhs| = ',d8.2,
     1        ' |sol| = ',d8.2,' |resid| = ',d8.2,' |mat| = ',d8.2,
     2        ' ratio = ',d8.2, ' refiter = ',i2)
         call C_OUT(1,0,1,line)
      endif
CCDELETEME
c      do i = 1, N+5*M+NLB+NUB
c         write(line,*) 'res ',i,RW(p_check+i), 'rhs ',i,RS(P_RHSCOPY+i)
c         call C_OUT(1,0,1,line)
c      enddo
C
C     If iterative refinement is wanted, compute new RHS for get_step_full
C
CTODO find good tolerance
C      if( RESTO.eq.1 .or. nrm_ratio.gt.1.d-11
c      nrm_ratio_max = min(1.d-10,max(1.d-14,1.d-8 * MU))
      nrm_ratio_max = MACHEPS
      if( (REFINE_ITER.lt.QREFINEITER .or. nrm_ratio.gt.nrm_ratio_max )
     1     .and. REFINE_ITER.lt.REFINE_ITER_MAX ) then
         REFINE_ITER = REFINE_ITER + 1

C     Compute RHS for get_step_full
         call RESFIL_RHS(N, M, NLB, ILB, NUB, IUB, RW(p_check+1),
     1        RHS_RESTO, SIGMA_PN, S_L, S_U, RS(P_P+1), RS(P_VP+1),
     2        RS(P_N+1), RS(P_VN+1), RS(P_RSCALE+1))
         call DCOPY(N+5*M+NLB+NUB, RW(p_check+1), 1, RS(P_RHSCOPY+1), 1)
         p_rwend = p_check

C     Store current steps
         call DCOPY(N, DX, 1, RS(P_DXSTORE+1), 1)
         call DCOPY(M, RW(p_dp+1), 1, RS(P_DPSTORE+1), 1)
         call DCOPY(M, RW(p_dn+1), 1, RS(P_DNSTORE+1), 1)
         call DCOPY(M, LAM, 1, RS(P_LAMSTORE+1), 1)
         call DCOPY(NLB, DV_L, 1, RS(P_DVLSTORE+1), 1)
         call DCOPY(NUB, DV_U, 1, RS(P_DVUSTORE+1), 1)
         call DCOPY(M, RW(p_dvp+1), 1, RS(P_DVPSTORE+1), 1)
         call DCOPY(M, RW(p_dvn+1), 1, RS(P_DVNSTORE+1), 1)
         RESTO = 2
         goto 9999
      endif
      if( nrm_ratio.ge.1.d-10 ) then
         write(line,*) 'nrm_ratio in restoration phase = ',nrm_ratio
         call C_OUT(2,1,1,line)
         RESTO = 3
         goto 7950
      endif
C
C     Done with iterative refinement
C
      if( QLAMBDA.eq.-2 ) then
         call DAXPY(M, 1.d0, LAMOLD, 1, LAM, 1)
      endif
      if( RESTO.eq.2 ) RESTO = 1
      p_rwend = p_check

      COUNT_RESTO_ITER = COUNT_RESTO_ITER + 1

      IREJ_LAST = 0
C
C     Compute maximal step sizes (due to iterative refinement, we cannot
C     use those from get_step_full)
C
      call CUTALPHA(N, X, DX, DV_L, DV_U, NLB, ILB, NUB, IUB,
     1     BNDS_L, BNDS_U, S_L, S_U, V_L, V_U, QTAU,
     2     ALPHA, ALPHA_DUAL, cdummy)

CTODO ALPHA might not be correct after iterative refinement!
C     ALPHA is already shortened to take bounds on X into account
C     Now let's do it for the P and N variables and also get directions
      do i = 1, M
         tmp = RW(p_dp+i)
         if( tmp.lt.0.d0 ) then
            ALPHA = dmin1(ALPHA, -QTAU*RS(P_P+i)/tmp)
         endif

         tmp = RW(p_dn+i)
         if( tmp.lt.0.d0 ) then
            ALPHA = dmin1(ALPHA, -QTAU*RS(P_N+i)/tmp)
         endif
      enddo
      if( QPRINT.ge.6 ) then
         write(line,7000) ITER
 7000    format(/,'  Restoration phase steps in P ',
     1             '(Iter ',i5,'):',/)
         call C_OUT(1,0,3,line)
         do i = 1, M
            write(line,7010) i,RW(p_dp+i)
 7010       format(' DP(',i6,') = ',d25.17)
            call C_OUT(1,0,1,line)
         enddo
         write(line,7020) ITER
 7020    format(/,'  Restoration phase steps in N ',
     1             '(Iter ',i5,'):',/)
         call C_OUT(1,0,3,line)
         do i = 1, M
            write(line,7030) i,RW(p_dn+i)
 7030       format(' DN(',i6,') = ',d25.17)
            call C_OUT(1,0,1,line)
         enddo
      endif

C     Compute objective function for restoration problem at current point
      resto_f = RHO_RESTO*DASUM(M, RS(P_P+1), 1)
      resto_f = resto_f + RHO_RESTO*DASUM(M, RS(P_N+1), 1)
      gd = 0.d0
      do i = 1, N
         tmp = X(i) - RS(P_XSTORE+i)
         resto_f = resto_f + ETA_RESTO*RHS_RESTO(N+M+i)*(tmp*tmp)/2.d0
         gd = gd + ETA_RESTO*RHS_RESTO(N+M+i)*tmp*DX(i)
      enddo
      resto_f = resto_f + CALC_BAR(NLB, NUB, S_L, S_U, V_L, V_U, MU)
      do i = 1, NLB
         j = ILB(i)
         gd = gd - MU/S_L(i)*DX(j)
      enddo
      do i = 1, NUB
         j = IUB(i)
         gd = gd + MU/S_U(i)*DX(j)
      enddo
      do i = 1, M
         ptmp = RS(P_P+i)
         ntmp = RS(P_N+i)
         resto_f = resto_f - MU*dlog(ptmp) - MU*dlog(ntmp)
         gd = gd + (RHO_RESTO - MU/ptmp)*RW(p_dp+i) +
     1        (RHO_RESTO - MU/ntmp)*RW(p_dn+i)
      enddo

C     Compute constraint violation at current point
      p_tmp   = p_rwend
      p_rwend = p_tmp + M
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      do i = 1, M
         RW(p_tmp+i) = C(i)*RS(P_RSCALE+i) - RS(P_P+i) + RS(P_N+i)
      enddo
      resto_cnrm = CALC_NRM(M, RW(p_tmp+1))
      p_rwend = p_tmp

      if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
         write(line,*) 'resto_f = ',resto_f,', resto_cnrm = ',resto_cnrm
         call C_OUT(1,0,1,line)
      endif

C     Compute minimal step size
      alpha_min = QGAMMA_THETA
      if( gd.lt.0.d0 ) then
         alpha_min = min(alpha_min, -QGAMMA_F*resto_cnrm/gd)
         if( resto_cnrm.le.QTHETAMINFACT ) then
            alpha_min = min(alpha_min,
     1           DELTA*(resto_cnrm**QS_THETA)/((-gd)**QS_F))
         endif
      endif
      alpha_min = ALPHA_MIN_FRAC*alpha_min
      alpha_min = max(alpha_min, QPRECFACT*MACHEPS)

      LS_COUNT = 0

C
C     Reserve memory for trail steps, trial points, and new constraint values
C
      p_pnew  = p_rwend
      p_nnew  = p_pnew + M
      p_xnew  = p_nnew + M
      p_cnew  = p_xnew + N
      p_slnew = p_cnew + M
      p_sunew = p_slnew + NLB
      p_rwend = p_sunew + NUB
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif

C ==============================================================================
C
C                       Permform the line search
C
C ==============================================================================
 10   continue
      LS_COUNT = LS_COUNT + 1
C
C     Compute trial point and slacks for current trial step size
C
      call DCOPY(N, X, 1, RW(p_xnew+1), 1)
      call DAXPY(N, ALPHA, DX, 1, RW(p_xnew+1), 1)
C
C     Move bounds a little bit if slacks become too small
C
      call CHECK_SLACKS(N, RW(p_xnew+1), NLB, ILB, S_L, V_L, BNDS_L,
     1     RW(p_slnew+1), NUB, IUB, S_U, V_U, BNDS_U, RW(p_sunew+1),
     1     MU, MACHEPS)
C
C     Compute constraints are trial point
C
      call GET_C(ITER, N, NIND, RW(p_xnew+1), IVAR, NORIG, XORIG, M,
     1     CSCALE, RW(p_cnew+1), KCONSTR, LRS, RS, LIS, IS,
     2     LRW-p_rwend, RW(p_rwend+1),
     4     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1        'resto_filter: Error: get_c returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.ne.0 ) then
         write(line,*)
     1        'resto_filter: Warning: get_c returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         IERR = 0
      endif
C
C     P and N for trial step size
C
      call DCOPY(M, RS(P_P+1), 1, RW(p_pnew+1), 1)
      call DAXPY(M, ALPHA, RW(p_dp+1), 1, RW(p_pnew+1), 1)
      call DCOPY(M, RS(P_N+1), 1, RW(p_nnew+1), 1)
      call DAXPY(M, ALPHA, RW(p_dn+1), 1, RW(p_nnew+1), 1)
C
C     Compute barrier function at trial point and g^T d
C
      resto_f_new = RHO_RESTO*DASUM(M, RW(p_pnew+1), 1)
      resto_f_new = resto_f_new + RHO_RESTO*DASUM(M, RW(p_nnew+1), 1)
      do i = 1, N
         tmp = RW(p_xnew+i) - RS(P_XSTORE+i)
         resto_f_new = resto_f_new +
     1        ETA_RESTO*RHS_RESTO(N+M+i)*(tmp*tmp)/2.d0
      enddo
      resto_f_new = resto_f_new + CALC_BAR(NLB, NUB, RW(p_slnew+1),
     1     RW(p_sunew+1), dummy, dummy, MU)
      do i = 1, M
         ptmp = RW(p_pnew+i)
         ntmp = RW(p_nnew+i)
         resto_f_new = resto_f_new - MU*dlog(ptmp) - MU*dlog(ntmp)
      enddo
      if( FFINITE(resto_f_new).eq.0 ) then
         call C_OUT(2,0,1,' --- WARNING: Evaluation of restoration'//
     1        ' objective function produced Nan/Inf.')
         irej = -2
         goto 800
      endif
C
C     Compute constraint violation at trial point
C
      p_tmp   = p_rwend
      p_rwend = p_tmp + M
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      do i = 1, M
         RW(p_tmp+i) = RW(p_cnew+i)*RS(P_RSCALE+i) - RW(p_pnew+i) +
     1        RW(p_nnew+i)
      enddo
      resto_cnrm_new = CALC_NRM(M, RW(p_tmp+1))
      p_rwend = p_tmp
      if( FFINITE(resto_cnrm_new).eq.0 ) then
         call C_OUT(2,0,1,' --- WARNING: Evaluation of constraint'//
     1        ' functions produced Nan/Inf.')
         irej = -2
         goto 800
      endif
C
C     Compute value of original objective function at new point
C
      call GET_F(N, RW(p_xnew+1), IVAR, NORIG, XORIG, M, CSCALE, NLB,
     1     ILB, RW(p_slnew+1), NUB, IUB, RW(p_sunew+1), MUORIG, fnew,
     1     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'resto_filter: get_f returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( FFINITE(fnew).eq.0 ) then
         call C_OUT(2,0,1,' --- WARNING: Evaluation of'//
     1        ' objective function produced Nan/Inf.')
         irej = -2
         goto 800
      endif
C
C     Check if constraint violation is 'exploding'
C
      if( resto_cnrm_new.gt.RESTO_CNRM_MAX ) then
         if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
            write(line,*)
     1           'resto_filter: cresto_cnrm_new = ',resto_cnrm_new,
     1           ' is getting too large.'
            call C_OUT(1,2,1,line)
         endif
         irej = -3
         goto 800
      endif
C
      if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
         write(line,*) 'resto_f_new = ',resto_f_new,
     1        'resto_cnrm_new = ',resto_cnrm_new
         call C_OUT(1,0,1,line)
      endif
C
C     Check if point is acceptable with respect to current point
C
CTODO also quadratic model here?
      model_alpha = ALPHA*gd
      if( model_alpha.ge.0.d0 ) then
         ftype = .false.
      else
         ftype = ( ((-model_alpha)**QS_F)*(ALPHA**(1.d0-QS_F))
     1        .gt.DELTA*(resto_cnrm**QS_THETA) .and.
     1        resto_cnrm.le.QTHETAMINFACT )
      endif

      if( ftype ) then
C
C     Check if enough decrease in barrier function
C
         lhs = resto_f-resto_f_new + PRECFACT*MACHEPS*dabs(resto_f)
         rhs = -ETA_F*model_alpha
         if( lhs.lt.rhs ) then
            irej = -1
            goto 800
         endif

      else
C
C     Check if "filter margin" satisfied at current iterate
C
         fi = resto_f
         ci = resto_cnrm
         if( resto_cnrm_new-(ci-QGAMMA_THETA*ci).gt.PRECFACT*MACHEPS*ci
     1        .and. resto_f_new-(fi-QGAMMA_F*ci).gt.
     1        PRECFACT*MACHEPS*dabs(fi) ) then
            irej = 0
            goto 800
         endif

      endif
C
C     Check if trial point is acceptable to the filter
C
      irej = FILTER_CHECK(NRFILTER, RS(P_RFILTERC+1), RS(P_RFILTERF+1),
     1     resto_cnrm_new, resto_f_new, PRECFACT, MACHEPS)
      if( irej.gt.0 ) goto 800
C
C     OK, point is accepted!
C
      goto 1000

 800  continue
      if( QCNR.gt.0 .and. QPRINT.gt.1 ) then
         write(line,*) 'LS_COUNT = ',LS_COUNT,', ALPHA = ',ALPHA
         call C_OUT(1,2,1,line)
         if( irej.eq.0 ) then
            call C_OUT(1,2,1,
     1           'No sufficient decrease w.r.t. current point.')
            write(line,*) 'resto_f   = ',resto_f,
     1           ' resto_cnrm = ',resto_cnrm
            call C_OUT(1,2,1,line)
C            write(line,*) 'resto_f_new = ',resto_f_new,
C     1           ' resto_cnrm_new = ',resto_cnrm_new
C            call C_OUT(1,2,1,line)
         elseif( irej.gt.0 ) then
            write(line,*) 'Filter entry that rejected: ',irej
            call C_OUT(1,2,1,line)
            write(line,*) 'F_filter = ',RS(P_RFILTERF+irej),
     1           ' C_filter = ',RS(P_RFILTERC+irej)
            call C_OUT(1,2,1,line)
C            write(line,*) 'resto_f_new = ',resto_f_new,
C     1           ' resto_cnrm_new = ',resto_cnrm_new
C            call C_OUT(1,2,1,line)
         elseif( irej.eq.-1 ) then
            write(line,*) 'f-type iteration test:'
            call C_OUT(1,2,1,line)
            write(line,*) 'resto_f = ',resto_f,' resto_f_new = ',
     1           resto_f_new
            call C_OUT(1,2,1,line)
            write(line,*) 'rhs = ',rhs,'  lhs = ',lhs
            call C_OUT(1,2,1,line)
         elseif( irej.eq.-2 ) then
            write(line,*) 'Nan occured...'
            call C_OUT(1,2,1,line)
         else
c            write(line,*) 'irej = ',irej
c            call C_OUT(1,2,1,line)
         endif
      endif
C
C     Cut back step size
C
      ALPHA = BACKTRACK_FAC*ALPHA
      IREJ_LAST = irej
C
      if( ALPHA.lt.ALPHA_MIN .or.
     1     (irej.eq.-3.and.COUNT_ITER_THETA_MAX.ge.QITERTHETAMAXMAX ))
     2     then
         if( irej.eq.-3.and.COUNT_ITER_THETA_MAX.ge.QITERTHETAMAXMAX )
     1        then
            RESTO_CNRM_MAX = 1.d2 * RESTO_CNRM_MAX
            call C_OUT(1,2,1, 'resto_filter: repeatedly hitting upper
     1bound on contraint violation.')
            write(line,*)
     1           '    setting RESTO_CNRM_MAX to ',RESTO_CNRM_MAX
            call C_OUT(1,2,1,line)
         endif
         COUNT_ITER_THETA_MAX = 0
C
C     Step size becomes too small - switch to feasibility restoration
C
         f_add  = resto_f - QGAMMA_F*resto_cnrm
         cnrm_add = resto_cnrm - QGAMMA_THETA*resto_cnrm
         call FILTER_ADD(NRFILTER, RS(P_RFILTERC+1), RS(P_RFILTERF+1),
     1        cnrm_add, f_add)
         C_ACCEPT = 'R'
         if( RESTO_RLAST ) then
            if( MU.gt.1.d-1*QTOL ) then
               call C_OUT(2,1,1,
     1'resto_filter: Repeated call of restoration phase.  Decrease MU.')
               decmu = .true.
               err_totl = 1.d300
               goto 5100
            endif
            call C_OUT(2,1,1,
     1'resto_filter: Repeated call of restoration phase for small MU.')
            call C_OUT(2,1,1,
     1'              MU is already too small to be decreased.')
            IERR = 18
            goto 9999
         endif
         call RESTO_RESTO(M, C, RS(P_RSCALE+1), RS(P_P+1), RS(P_N+1),
     1        RS(P_VP+1), RS(P_VN+1), LAM, MU, RHO_RESTO)

CTODO Decide whether we want to update the dual variables, too...
         if( QPRINT.ge.2 ) then
            call C_OUT(1,0,1,'resto_filter: do restoration phase.')
         endif
         do i = 1, M
            RS(P_VN+i) = MU/RS(P_N+i)
            RS(P_VP+i) = MU/RS(P_P+i)
         enddo
         RESTO_RLAST = .true.
         p_rwend = p_dp
         goto 5000
      endif

      goto 10

 1000 continue
      RESTO_RLAST = .false.
C
C     Check if filter block progress and needs to be deleted
C
      if( IREJ_LAST.gt.0 ) then
CDELETEME
         COUNT_ITER_BLOCK = COUNT_ITER_BLOCK + 1
         if( QITERBLOCKMAX.gt.0 .and.
     1        COUNT_ITER_BLOCK.gt.QITERBLOCKMAX) then
            call C_OUT(2,0,1,
     1      'resto_filter: It is possible that filter blocks progress.')
            if( RESTO_CNRM_MAX.gt.10d0*resto_cnrm_new ) then
               write(line,*)
     1    '      Max allowed constraint violation is decreased from ',
     2              RESTO_CNRM_MAX
               call C_OUT(2,0,1,line)
               RESTO_CNRM_MAX = RESTO_CNRM_MAX/10d0
               write(line,*)
     1        '       to ',RESTO_CNRM_MAX,' and current filter deleted.'
               call C_OUT(2,0,1,line)
               NRFILTER = 0
            else
               call C_OUT(2,0,1,
     1              '      However, cannot delete filter.')
            endif
         endif
      else
         COUNT_ITER_BLOCK = 0
      endif
C
C     Count how many consecutive time the upper bound on theta cut step
C
      if( IREJ_LAST.eq.-3 ) then
         COUNT_ITER_THETA_MAX = COUNT_ITER_THETA_MAX + 1
      else
         COUNT_ITER_THETA_MAX = 0
      endif
C
C     Check if filter is to be augmented
C
      lhs   = resto_f_new - resto_f - PRECFACT*MACHEPS*dabs(resto_f)
      rhs   = ETA_F*model_alpha
      if( model_alpha.ge.0.d0 ) then
         augment = .true.
      else
         augment = ( ((-model_alpha)**QS_F)*(ALPHA**(1.d0-QS_F))
     1        .le.DELTA*(resto_cnrm_new**QS_THETA) .or. lhs.gt.rhs)
      endif
      if( .not.augment ) then
         C_ACCEPT = 'f'
      else
         C_ACCEPT = 'h'
         f_add  = resto_f - QGAMMA_F*resto_cnrm
         cnrm_add = resto_cnrm - QGAMMA_THETA*resto_cnrm
         call FILTER_ADD(NRFILTER, RS(P_RFILTERC+1), RS(P_RFILTERF+1),
     1        cnrm_add, f_add)
      endif
C
C     Copy trial points to real points
C
      call DCOPY(N, RW(p_xnew+1), 1, X, 1)
      call DCOPY(NLB, RW(p_slnew+1), 1, S_L, 1)
      call DCOPY(NUB, RW(p_sunew+1), 1, S_U, 1)
      call DCOPY(M, RW(p_cnew+1), 1, C, 1)
      p_rwend = p_xnew
C
C     Make ALPHA_DUAL shorter if necessary
C
      do i = 1, M
         tmp = RW(p_dvp+i)
         if( tmp.lt.0.d0 ) then
            ALPHA_DUAL = dmin1(ALPHA_DUAL, -QTAU*RS(P_VP+i)/tmp)
         endif

         tmp = RW(p_dvn+i)
         if( tmp.lt.0.d0 ) then
            ALPHA_DUAL = dmin1(ALPHA_DUAL, -QTAU*RS(P_VN+i)/tmp)
         endif
      enddo
      if( QPRINT.ge.5 ) then
         write(line,7050) ITER
 7050    format(/,'  Restoration phase variables DVP ',
     1             '(Iter ',i5,'):',/)
         call C_OUT(1,0,3,line)
         do i = 1, M
            write(line,7060) i,RW(p_dvp+i)
 7060       format(' DVP(',i6,') = ',d25.17)
            call C_OUT(1,0,1,line)
         enddo
         write(line,7070) ITER
 7070    format(/,'  Restoration phase variables DVN ',
     1             '(Iter ',i5,'):',/)
         call C_OUT(1,0,3,line)
         do i = 1, M
            write(line,7080) i,RW(p_dvn+i)
 7080       format(' DVN(',i6,') = ',d25.17)
            call C_OUT(1,0,1,line)
         enddo
      endif

      if( QALPHA.eq.0 ) then
         ALPHA_DUAL = ALPHA
      elseif( QALPHA.eq.1 ) then
         ALPHA_DUAL = min( ALPHA_DUAL, ALPHA )
      endif

      call DCOPY(M, RW(p_pnew+1), 1, RS(P_P+1), 1)
      call DCOPY(M, RW(p_nnew+1), 1, RS(P_N+1), 1)
      call DAXPY(NLB, ALPHA_DUAL, DV_L, 1, V_L, 1)
      call DAXPY(NUB, ALPHA_DUAL, DV_U, 1, V_U, 1)
      call DAXPY(M, ALPHA_DUAL, RW(p_dvp+1), 1, RS(P_VP+1), 1)
      call DAXPY(M, ALPHA_DUAL, RW(p_dvn+1), 1, RS(P_VN+1), 1)
      p_rwend = p_dp
C
C     Make sure V doesn't differ too much form MU S^{-1} for p and n variables
C
      call CORRECT_V(M, RS(P_P+1), RS(P_VP+1), M, RS(P_N+1),
     1     RS(P_VN+1), MU)

      if( abs(QLAMBDA).eq.2 ) then
         if( RESTO_INIT_LAM ) then
            RESTO_INIT_LAM = .false.
         else
            call DSCAL(M, ALPHA, LAM, 1)
            call DAXPY(M, 1d0-ALPHA, LAMOLD, 1, LAM, 1)
         endif
         call DCOPY(M, LAM, 1, LAMOLD, 1)
      endif
C
C     Compute value of barrier function at new point
C
      F = fnew
      bar_n = CALC_BAR(NLB, NUB, S_L, S_U, V_L, V_U, MUORIG)
      phi_new = F + bar_n
      if( abs(QMERIT).eq.5 ) then
         phi_new = phi_new + DDOT(M, C, 1, LAM, 1)
      endif
      if( FFINITE(phi_new).eq.0 ) then
         call C_OUT(2,0,1,' --- WARNING: Evaluation of barrier'//
     1        ' objective function produced Nan/Inf.')
         goto 5000
      endif
C
C     Check if new point has small enough contraint violation compared to
C     initial constraint violation
C
      cnrm_new = CALC_NRM(M, C)
      if( NIND.eq.0 ) then
         i = IDAMAX(M, C, 1)
         if( abs(C(i)).lt.QTOL ) then
            call C_OUT(2,0,1,
     1           'Maximal constraint violation below tolerance.')
            IERR = 16
            goto 9999           ! Done
         endif
         goto 5000
      else
         if( cnrm_new.gt.QRESTOREDFACT*CNRM_START ) goto 5000
      endif
C
C     Check if new point is acceptable to the original filter
C
      CNRM = CALC_NRM(M, C)
      do i = 1, NFILTER
C     Check f component
         fi = FILTER_PHI(i)
         lhs = phi_new - fi - PRECFACT*MACHEPS*dabs(fi)
         if( lhs.le.0d0 ) goto 4000
C     Check h component
         ci = FILTER_C(i)
         lhs = CNRM - ci - PRECFACT*MACHEPS*ci
         if( lhs.le.0d0 ) goto 4000
C     Point not OK for this entry
         if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
            write(line,*)
     1           'In restoration: new point not compatible with entry',
     2           i,' of filter.'
            call C_OUT(1,2,1,line)
         endif
         goto 5000
C
 4000    continue
      enddo
C
C     Check if constraint violation not above maximal allowed value
C
      if( CNRM.gt.CNRM_MAX ) then
         if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
            write(line,*)
     1           'In restoration: new point larger than CNRM_MAX.'
            call C_OUT(1,2,1,line)
         endif
         goto 5000
      endif
C
C     Filter tests passed
C
      RESTO = 0
C
C     Restore original MU and TAU
C
      MU   = MUORIG
      QTAU = TAUSTORE

C     If we do L-BFGS, restore current estimate
      if( abs(QQUASI).ge.6 ) then
         call DCOPY(2+2*N+QLMLEN*(1+4*N+2*QLMLEN),
     1        RS(P_BSTORE+1), 1, B, 1)
C     Set XOLD to XS so that next update will be skipped
         call DCOPY(N, X, 1, B(2+QLMLEN*(1+4*N+2*QLMLEN)+1), 1)
      endif
C
C     Continue with normal IPOPT
C
      goto 9999

 5000 continue
      mu_changed = .false.
cCDELETEME
c      if( ITER.eq.ITER_START )
c     1     mu_changed = .true.
C
C     Compute right hand side and error of barrier restoration problem
C
 5001 continue
C     1. dual infeasibility for X
      irhs = 0
      call DCOPY(N, RS(P_XSTORE+1), 1, RS(P_RHS+irhs+1), 1)
      call DAXPY(N, -1.d0, X, 1, RS(P_RHS+irhs+1), 1)
      do i = 1, N
         RS(P_RHS+irhs+i) = RS(P_RHS+irhs+i)*ETA_RESTO*RHS_RESTO(N+M+i)
      enddo

      do i = 1, NLB
         j = ILB(i)
         RS(P_RHS+irhs+j) = RS(P_RHS+irhs+j) + V_L(i)
      enddo
      do i = 1, NUB
         j = IUB(i)
         RS(P_RHS+irhs+j) = RS(P_RHS+irhs+j) - V_U(i)
      enddo

      p_tmp   = p_rwend
      p_rwend = p_tmp + N
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C     in the following call we need to specify ITER+1, so that the Jacobian
C     is evaluated at the NEW point
      call CONSTR(8, ITER+1, N, NIND, M, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, LAM, RW(p_tmp+1), idummy,
     3     idummy, KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4     IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.lt.0 ) then
         write(line,*) 'resto_filter: Warning in CONSTR-8, IERR = ',IERR
         call C_OUT(2,0,1,line)
      elseif( IERR.ne.0 ) then
         write(line,*) 'resto_filter: Error in CONSTR-8, IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( QLAMBDA.eq.-2 ) then
         call DAXPY(N, -1.d0, RW(p_tmp+1), 1, RS(P_RHS+irhs+1), 1)
         i = IDAMAX(N, RS(P_RHS+irhs+1), 1)
         err_dinf = dabs(RS(P_RHS+irhs+i))
      else
         call DAXPY(N, -1.d0, RS(P_RHS+irhs+1), 1, RW(p_tmp+1), 1)
         i = IDAMAX(N, RW(p_tmp+1), 1)
         err_dinf = dabs(RW(p_tmp+i))
      endif
      if( QPRINT.ge.6 ) then
         write(line,*) 'err_dinf(x)  = ',err_dinf
         call C_OUT(1,0,1,line)
      endif
      p_rwend = p_tmp

      irhs = irhs + N

C     2. dual infeasinbility for P
      call DCOPY(M, RS(P_VP+1), 1, RS(P_RHS+irhs+1), 1)
      call DAXPY(M, -1.d0, RHO_RESTO, 0, RS(P_RHS+irhs+1), 1)
      p_tmp   = p_rwend
      p_rwend = p_tmp + M
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      do i = 1, M
         RW(p_tmp+i) = LAM(i)/RS(P_RSCALE+i)
      enddo
      if( QLAMBDA.eq.-2 ) then
         call DAXPY(M, 1.d0, RW(p_tmp+1), 1, RS(P_RHS+irhs+1), 1)
         i = IDAMAX(M, RS(P_RHS+irhs+1), 1)
         dummy = abs(RS(P_RHS+irhs+i))
      else
         call DAXPY(M, 1.d0, RS(P_RHS+irhs+1), 1, RW(p_tmp+1), 1)
         i = IDAMAX(M, RW(p_tmp+1), 1)
         dummy = abs(RW(p_tmp+i))
      endif
      if( QPRINT.ge.6 ) then
         write(line,*) 'err_dinf(p)  = ',dummy
         call C_OUT(1,0,1,line)
      endif
      err_dinf = max(err_dinf, dummy)
      irhs = irhs + M

C     3. dual infeasinbility for P
      call DCOPY(M, RS(P_VN+1), 1, RS(P_RHS+irhs+1), 1)
      call DAXPY(M, -1.d0, RHO_RESTO, 0, RS(P_RHS+irhs+1), 1)
      do i = 1, M
         RW(p_tmp+i) = - LAM(i)/RS(P_RSCALE+i)
      enddo
      if( QLAMBDA.eq.-2 ) then
         call DAXPY(M, 1.d0, RW(p_tmp+1), 1, RS(P_RHS+irhs+1), 1)
         i = IDAMAX(M, RS(P_RHS+irhs+1), 1)
         dummy = abs(RS(P_RHS+irhs+i))
      else
         call DAXPY(M, 1.d0, RS(P_RHS+irhs+1), 1, RW(p_tmp+1), 1)
         i = IDAMAX(M, RW(p_tmp+1), 1)
         dummy = abs(RW(p_tmp+i))
      endif
      if( QPRINT.ge.6 ) then
         write(line,*) 'err_dinf(n)  = ',dummy
         call C_OUT(1,0,1,line)
      endif
      err_dinf = max(err_dinf, dummy)
      irhs = irhs + M
      p_tmp = p_rwend

C     4. primal infeasibiliy
      call DCOPY(M, RS(P_P+1), 1, RS(P_RHS+irhs+1), 1)
      call DAXPY(M, -1.d0, RS(P_N+1), 1, RS(P_RHS+irhs+1), 1)
      do i = 1, M
         RS(P_RHS+irhs+i) = RS(P_RHS+irhs+i) - C(i)*RS(P_RSCALE+i)
      enddo
      i = IDAMAX(M, RS(P_RHS+irhs+1), 1)
      err_pinf = dabs(RS(P_RHS+irhs+i))
      if( QPRINT.ge.6 ) then
         write(line,*) 'err_pinf     = ',err_pinf
         call C_OUT(1,0,1,line)
      endif
      do i = 1, M
         RS(P_RHS+irhs+i) = RS(P_RHS+irhs+i)/RS(P_RSCALE+i)
      enddo
      irhs = irhs + M

C     5. complementarity X_L
      do i = 1, NLB
         RS(P_RHS+irhs+i) = MU - S_L(i)*V_L(i)
      enddo
      if( NLB.gt.0 ) then
         i = IDAMAX(NLB, RS(P_RHS+irhs+1), 1)
         err_cmpl = dabs(RS(P_RHS+irhs+i))
         if( QPRINT.ge.6 ) then
            write(line,*) 'err_cmpl(xl) = ',err_cmpl
            call C_OUT(1,0,1,line)
         endif
      else
         err_cmpl = 0.d0
      endif
      irhs = irhs + NLB

C     6. complementarity X_U
      do i = 1, NUB
         RS(P_RHS+irhs+i) = MU - S_U(i)*V_U(i)
      enddo
      if( NUB.gt.0 ) then
         i = IDAMAX(NUB, RS(P_RHS+irhs+1), 1)
         if( QPRINT.ge.6 ) then
            write(line,*) 'err_cmpl(xu) = ',dabs(RS(P_RHS+irhs+i))
            call C_OUT(1,0,1,line)
         endif
         err_cmpl = max(err_cmpl, dabs(RS(P_RHS+irhs+i)))
      endif
      irhs = irhs + NUB

C     7. complementarity P
      do i = 1, M
         RS(P_RHS+irhs+i) = MU - RS(P_P+i)*RS(P_VP+i)
      enddo
      i = IDAMAX(M, RS(P_RHS+irhs+1), 1)
      if( QPRINT.ge.6 ) then
         write(line,*) 'err_cmpl(p)  = ',dabs(RS(P_RHS+irhs+i))
         call C_OUT(1,0,1,line)
      endif
      err_cmpl = max(err_cmpl, dabs(RS(P_RHS+irhs+i)))
      irhs = irhs + M

C     8. complementarity P
      do i = 1, M
         RS(P_RHS+irhs+i) = MU - RS(P_N+i)*RS(P_VN+i)
      enddo
      i = IDAMAX(M, RS(P_RHS+irhs+1), 1)
      if( QPRINT.ge.6 ) then
         write(line,*) 'err_cmpl(n)  = ',dabs(RS(P_RHS+irhs+i))
         call C_OUT(1,0,1,line)
      endif
      err_cmpl = max(err_cmpl, dabs(RS(P_RHS+irhs+i)))
      irhs = irhs + M
      
C     Scale errors if requested
      if( QERROR.gt.0 ) then
         scal_cmpl = DASUM(NLB, V_L, 1) + DASUM(NUB, V_U, 1) +
     1        DASUM(M, RS(P_VP+1), 1) + DASUM(M, RS(P_VN+1), 1)
         scal_dinf = scal_cmpl
         do i = 1, M
            scal_dinf = scal_dinf + abs(LAM(i))/RS(P_RSCALE+i)
         enddo
         scal_dinf = max(QSCALMULTMAX, scal_dinf/dble(NLB+NUB+3*M))/
     1        QSCALMULTMAX
         scal_cmpl = max(QSCALMULTMAX, scal_cmpl/dble(NLB+NUB+2*M))/
     1        QSCALMULTMAX
         err_dinf = err_dinf/scal_dinf
         err_cmpl = err_cmpl/scal_cmpl
      endif

      err_totl = max(err_dinf, err_pinf, err_cmpl)
      if( QPRINT.ge.2 ) then
         write(line(1),*) 'In RESTO: err_dinf = ',err_dinf
         write(line(2),*) '          err_pinf = ',err_pinf
         write(line(3),*) '          err_cmpl = ',err_cmpl
         write(line(4),*) '          err_totl = ',err_totl
         call C_OUT(1,0,4,line)
      endif
      if( err_totl.lt.QTOL .and. (MU.lt.QTOL.and..not.mu_changed) ) then
         write(line,*) 'Restoration phase problem converged.'
         call C_OUT(2,0,1,line)

         i = IDAMAX(M, C, 1)
         if( abs(C(i)).lt.QTOL ) then
            write(line,*) 'Maximum constraint violation ',abs(C(i)),
     1           ' below tolerance.'
            call C_OUT(2,0,1,line)
            IERR = 18
            goto 9999
         elseif( abs(C(i)).lt.max(sqrt(QTOL),1.d2*QTOL) ) then
            write(line,*) 'Maximum constraint violation ',abs(C(i)),
     1           '.'
            call C_OUT(2,0,1,line)
            IERR = 18
            goto 9999
         endif

         IERR = 17
         goto 9999
      endif
C
C     A sophisticated user might want to decided on some other criterion
C     of her/his own to decrease the barrier parameter earlier.  In that
C     case you have to define the preprocessor flag
C     '-DUSE_IP_MYRESTOBARRIERTEST' at compilation time and return .true. as
C     result of the function below.
C
      decmu = .false.
#ifdef USE_IP_MYRESTOBARRIERTEST
      decmu = IP_MYRESTOBARRIERTEST(NORIG, XORIG, N, X, M, C, F, ITER,
     1     ITER_START, MU, ALPHA, RESTO, IERR)
      if( IERR.ne.0 ) then
         write(line,*)
     1        'resto_filter: IP_MYRESTOBARRIERTEST returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
#endif
 5100 continue

      if( ((err_totl.lt.RERR) .and. C_ACCEPT.ne.'R'
CORIG     1     .and. C_ACCEPT.ne.'I' .and. .not.mu_changed)
     1     .and. C_ACCEPT.ne.'I' .and. .not.mu_changed .and.ITER.ne.0)
     2     .or. decmu ) then
C
C     Update MU, the error tolerance, and QTAU
C
CTODO
c         if( ITER_START.ne.ITER ) then
            mu_changed = .true.
c         endif
         first_mu = .false.
         call UPDATE_MU(first_mu, RESTO, MU, RERR)

CDELETEME ?
c         if( ITER_START.eq.ITER ) then
c            call RESTO_RESTO(M, C, RS(P_RSCALE+1), RS(P_P+1), RS(P_N+1),
c     1           RS(P_VP+1), RS(P_VN+1), LAM, MU, RHO_RESTO)
c            do i = 1, M
c               RS(P_VP+i) = MU/RS(P_P+i)
c               RS(P_VN+i) = MU/RS(P_N+i)
c            enddo
c         endif

         xstorenorm = max(1.d0,DNRM2(N, RS(P_XSTORE+1), 1))
#ifndef OLD
         xstorenorm = 1.d0
#endif
         

         if( QETAUPDATE.eq.1 ) then
            ETA_RESTO = ETARESTOFACT*MU/xstorenorm
         elseif( QETAUPDATE.eq.2 ) then
            ETA_RESTO = ETARESTOFACT*MU/xstorenorm
         elseif( QETAUPDATE.eq.3 ) then
c            ETA_RESTO = ETARESTOFACT*min(1.d0,MU)/xstorenorm/xstorenorm
c            ETA_RESTO = ETARESTOFACT*MU/xstorenorm/xstorenorm
            ETA_RESTO = ETARESTOFACT*dsqrt(MU)/xstorenorm/xstorenorm
         else
            ETA_RESTO = ETARESTOFACT*dsqrt(MU)/xstorenorm
         endif

         if( NIND.eq.0 ) ETA_RESTO = 0.d0

         if( QPRINT.ge.1 ) then
            write(line,7900) MU, ETA_RESTO
 7900       format('After update:  MU = ',d12.4,' ETA_RESTO = ',d12.4)
            call C_OUT(1,0,1,line)
         endif

         NRFILTER = 0
C
C     Since MU has just been changed, we need to recompute the RHS
C
         goto 5001

      endif

 7950 continue
C     Compute RHS for get_step_full
      call RESFIL_RHS(N, M, NLB, ILB, NUB, IUB, RS(P_RHS+1), RHS_RESTO,
     1     SIGMA_PN, S_L, S_U, RS(P_P+1), RS(P_VP+1), RS(P_N+1),
     2     RS(P_VN+1), RS(P_RSCALE+1))
      call DCOPY(N+5*M+NLB+NUB, RS(P_RHS+1), 1, RS(P_RHSCOPY+1), 1)
      REFINE_ITER = 0
C
C     Some logging
C
      if( QPRINT.ge.5 .and. RESTO.ne.3 ) then
         write(line,8000) ITER
 8000    format(/,'  Restoration phase variables P ',
     1             '(Iter ',i5,'):',/)
         call C_OUT(1,0,3,line)
         do i = 1, M
            write(line,8010) i,RS(P_P+i)
 8010       format(' P(',i6,') = ',d25.17)
            call C_OUT(1,0,1,line)
         enddo
         write(line,8020) ITER
 8020    format(/,'  Restoration phase variables N ',
     1             '(Iter ',i5,'):',/)
         call C_OUT(1,0,3,line)
         do i = 1, M
            write(line,8030) i,RS(P_N+i)
 8030       format(' N(',i6,') = ',d25.17)
            call C_OUT(1,0,1,line)
         enddo
         write(line,8040) ITER
 8040    format(/,'  Values of constraints ',
     1             '(Iter ',i5,'):',/)
         call C_OUT(1,0,3,line)
         do i = 1, M
            write(line,8050) i,C(i)
 8050       format(' C(',i6,') = ',d25.17)
            call C_OUT(1,0,1,line)
         enddo
         write(line,8060) ITER
 8060    format(/,'  Restoration phase dual variables VP ',
     1             '(Iter ',i5,'):',/)
         call C_OUT(1,0,3,line)
         do i = 1, M
            write(line,8070) i,RS(P_VP+i)
 8070       format(' VP(',i6,') = ',d25.17)
            call C_OUT(1,0,1,line)
         enddo
         write(line,8080) ITER
 8080    format(/,'  Restoration phase dual variables VN ',
     1             '(Iter ',i5,'):',/)
         call C_OUT(1,0,3,line)
         do i = 1, M
            write(line,8090) i,RS(P_VN+i)
 8090       format(' VN(',i6,') = ',d25.17)
            call C_OUT(1,0,1,line)
         enddo
      endif

 9999 continue
      return
      end

C
      subroutine RESTO_RESTO(M, C, RSCALE, P, N, VP, VN, LAM, MU,
     1     RHO_RESTO)
C
      implicit none
C
C     Restoration within restoration phase
C
C     Find feasible point by changing P and N so that the objective function
C     in P and N (with fixed X) is minimized.  Set multipliers VP and VN
C     so that they satisfy complementarity for MU
C
      integer M
      double precision C(M), RSCALE(M), P(M), N(M), VP(M), VN(M), LAM(M)
      double precision MU, RHO_RESTO
C
      integer i
      double precision pp, qq, cc

      do i = 1, M
         cc = C(i)*RSCALE(i)
         pp = (cc - MU/RHO_RESTO)/2.d0
         qq = MU*cc/(2.d0*RHO_RESTO)
         N(i)  = -pp + sqrt(pp*pp + qq)
         P(i)  = N(i) + cc
      enddo
      return
      end

      subroutine RESFIL_RHS(N, M, NLB, ILB, NUB, IUB, RHS, RHS_RESTO,
     1     SIGMA_PN, S_L, S_U, P, VP, NN, VN, RSCALE)
      implicit none
      integer N, M, NLB, ILB(NLB), NUB, IUB(NUB)
      double precision RHS(N+5*M+NLB+NUB)
      double precision RHS_RESTO(N+M)
      double precision SIGMA_PN(M)
      double precision S_L(NLB), S_U(NUB), P(M), VP(M), NN(M), VN(M)
      double precision RSCALE(M)

      integer i, j, irhs
C
C     Compute the right hand side that is actually passed to get_step_full
C
      call DCOPY(N, RHS, 1, RHS_RESTO, 1)
      call DCOPY(M, RHS(N+2*M+1), 1, RHS_RESTO(N+1), 1)
      irhs = N
      do i = 1, M
         RHS_RESTO(N+i) = RHS_RESTO(N+i) +
     1        P(i)/VP(i)*RHS(irhs+i)/RSCALE(i)
      enddo
      irhs = irhs + M
      do i = 1, M
         RHS_RESTO(N+i) = RHS_RESTO(N+i) -
     1        NN(i)/VN(i)*RHS(irhs+i)/RSCALE(i)
      enddo
      irhs = irhs + 2*M
      do i = 1, NLB
         j = ILB(i)
         RHS_RESTO(j) = RHS_RESTO(j) + RHS(irhs+i)/S_L(i)
      enddo
      irhs = irhs + NLB
      do i = 1, NUB
         j = IUB(i)
         RHS_RESTO(j) = RHS_RESTO(j) - RHS(irhs+i)/S_U(i)
      enddo
      irhs = irhs + NUB
      do i = 1, M
         RHS_RESTO(N+i) = RHS_RESTO(N+i) + RHS(irhs+i)/VP(i)/RSCALE(i)
      enddo
      irhs = irhs + M
      do i = 1, M
         RHS_RESTO(N+i) = RHS_RESTO(N+i) - RHS(irhs+i)/VN(i)/RSCALE(i)
      enddo
      irhs = irhs + M
C
C     Sigma for P and N variables
C
      do i = 1, M
         SIGMA_PN(i) = -(P(i)/VP(i)+NN(i)/VN(i))/(RSCALE(i)**2)
      enddo
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine RESTO_FILTER_WS(N, M, NLB, NUB, NZA, LRW, LIW,
     1     DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1, lrwc, liwc, lrwf, liwf, lrwconstr, liwconstr
      character*80 line

      call GET_C_WS(N, M, NLB, NUB, NZA, lrwc, liwc, DAT, IDAT)
      call GET_F_WS(N, M, NLB, NUB, NZA, lrwf, liwf)
      call CONSTR_WS(N, M, NLB, NUB, NZA, lrwconstr, liwconstr,
     1     DAT, IDAT)

      LIW = max( liwc, liwf, liwconstr )
      LRW = 4*M + max( 5*M+N+NLB+NUB + N + lrwconstr,
     1     3*M+N+NLB+NUB + max(lrwc,lrwf,M) )

      if( QPRINT.ge.3 ) then
         write(line,1000)'resto_filter_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
