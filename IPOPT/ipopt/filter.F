C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      SUBROUTINE FILTER(ITER, N, NIND, M, X, IVAR, NFIX, IFIX,
     1     NLB, ILB, NUB,
     1     IUB, BNDS_L, BNDS_U, DX, DV_L, DV_U, S_L, S_U,
     2     V_L, V_U, SIGMA_L, SIGMA_U, PZ, NORIG, XORIG, CSCALE,
     3     MU, MUORIG, LAM, LAMOLD, G, F, C, CNRM0, ALPHA, ALPHA_DUAL,
     3     RESTO, LS_COUNT, NU_OUT, C_ALPHA, C_ACCEPT,
     4     SOC_FLAG, WFLAG, INIT_LAM, ERR, RHS_RESTO, SIGMA_PN,
     5     ETA_RESTO, B, DHD, REGU, KCONSTR, LRS, RS, LIS, IS,
     6     LRW, RW, LIW, IW, IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: filter.F 664 2004-11-05 18:46:52Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Filter method implementation with line search
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
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP                            (assumed to be in increasing order)
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
CP   SIGMA_L   I    DP     primal-dual Hessian of lower bound barrier term
CP                            (NLB diagonal elements only)
CP   SIGMA_U   I    DP     primal-dual Hessian of upper bound barrier term
CP                            (NUB diagonal elements only)
CP   PZ       I/O   DP     kept update-to-date depending on SOC
CP   NORIG     I    INT    number of all variables including fixed vars
CP   XORIG    I/O   DP     actual iterate
CP                            XORIG is ordered in ORIGINAL order (i.e. not
CP                            partitioned into independent and dependent
CP                            variables) (on output: as X)
CP   CSCALE    I    DP     scaling factors for constraints
CP   MU        I    DP     barrier parameter
CP   MUORIG    I    DP     barrier parameter (for original problem)
CP                            (only when in restoration phase) 
CP   LAM      I/O   DP     multipliers for equality constraints
CP   LAMOLD   I/O   DP     multipliers for equality constraints from last iter
CP   G         I    DP     gradient of objective function
CP   F        I/O   DP     value of objective function at X
CP                            I: for start of line search
CP                            O: after line search
CP   C        I/O   DP     values of constraints at X
CP                            I: for start of line search
CP                            O: after line search (or for SOC computation)
CP   CNRM0     O    DP     2-norm of constraints at old point
CP   ALPHA    I/O   DP     step size: I: where to start line search (primal)
CP                                    O: step size from X to X_NEW
CP   ALPHA_DUAL I/O DP     step size for dual variables
CP   RESTO    I/O   INT    flag for restoration phase:
CP                            0: not in restoration phase
CP                           -1: in KKT reduction restoration phase
CP                            1: in l_1 restoration phase
CP                            2: in l_1 restoration phase (iterative refinement)
CP   LS_COUNT  O    INT    number of trial steps
CP   NU_OUT    O    DP     actual value of penalty parameter
CP                            (only output; value is stored internally!)
CP   C_ALPHA  I/O   C*1    This flag is set to 'R' if we are in the restoration
CP                            phase.
CP   C_ACCEPT  O    C*1    ='f' if improvement in barrier function
CP                         ='c' if improvement in constraint violation
CP                         ='b' if improvement in both
CP   SOC_FLAG I/O   INT    Flag for Second order correction:
CP                           =0: no SOC has been tried yet (I)
CP                           =1: try SOC correction now (O)
CP                               input step corresponds to SOC (I)
CP                           =2: SOC has been tried but not successful
CP   WFLAG     O    INT    =0: inactive watchdog
CP                         =1: will now compute trial watchdog step
CP                         =2: watchdog step didn't help
CP   INIT_LAM I/O   LOG    Is changed by this routine, if the equality
CP                           multipliers are to be reinitialized.
CP                           (e.g. after restoration phase)
CP   ERR       I    DP     current value of the KKT error of the barrier problem
CP   RHS_RESTO O    DP     QRESTO.eq.1: RHS if in IP restoration
CP   SIGMA_PN  O    DP     QRESTO.eq.1: Sigma for P and N variables
CP   ETA_RESTO O    DP     QRESTO.eq.1: Factor for regularization term
CP   B        I/O   DP     QRESTO.eq.1: QN estimate (to be restored after
CP                           restoration phase)
CP   DHD       I    DP     d^T Hessian d
CP   REGU      I    DP     regularization parameter (get_step_full)
CP   KCONSTR   I    INT    KCONSTR(1): LRS for CONSTR
CP                         KCONSTR(2): P_LRS for CONSTR
CP                         KCONSTR(3): LIS for CONSTR
CP                         KCONSTR(4): P_LIS for CONSTR
CP                         KCONSTR(5): LRW for CONSTR
CP                         KCONSTR(6): LIW for CONSTR
CP   LRS      I/O   INT    (ITER.ne.-1) I: total length of RS
CP                         (ITER.eq.-1) O: increased by needed storage space
CP   RS       I/O   DP     DP storage space (all!)
CP   LIS      I/O   INT    (ITER.ne.-1) I: total length of IS
CP                         (ITER.eq.-1) O: increased by needed storage space
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
CS    DNRM2
CS    DSCAL
CS    C_OUT
CS    CALC_BAR
CS    CALC_NRM
CS    GET_F
CS    GET_C
CS    FFINITE
CS    RESTO_TRON
CS    RESTO_KKT
CS    RESTO_DOGLEG
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
      double precision X(N)
      integer IVAR(N)
      integer NFIX
      integer IFIX(NFIX)
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
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision PZ(NIND)
      integer NORIG
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      double precision MU
      double precision MUORIG
      double precision LAM(M)
      double precision LAMOLD(M)
      double precision G(N)
      double precision F
      double precision C(M)
      double precision CNRM0
      double precision ALPHA
      double precision ALPHA_DUAL
      integer RESTO
      integer LS_COUNT
      double precision NU_OUT
      character*1 C_ALPHA
      character*1 C_ACCEPT
      integer SOC_FLAG
      integer WFLAG
      logical INIT_LAM
      double precision ERR
      double precision RHS_RESTO(N+M)
      double precision SIGMA_PN(M)
      double precision ETA_RESTO
      double precision B(*)
      double precision DHD
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
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer p_rwend, p_iwend, p_xnew, p_slnew, p_sunew
      integer p_cnew, p_tmp, p_gb, p_lamnew, p_dl, p_du, p_tmp2
      integer i, k, irej, ierr_store
      logical ftype, augment
      double precision cnrm_new, phi_new, f_new, bar_n, bar0
      double precision lhs, rhs, fi, ci, dummy, alphad
      double precision cnrm_add, phi_add, alpha_tmp, dltd
      double precision tmp1, tmp2, tmp3, tmp4, model_alpha, gnrm
      double precision reldxnrm, factor, ctmp
      character*150 line(4)

      double precision CALC_NRM, CALC_BAR, D1MACH, DDOT, FILTER_MODEL
      double precision DNRM2
      integer FFINITE
      integer FILTER_CHECK

      integer NFILTER, P_FILTER_PHI, P_FILTER_C, RESTO_KKT_FLAG
      save    NFILTER, P_FILTER_PHI, P_FILTER_C, RESTO_KKT_FLAG

      integer COUNT_ITER_BLOCK, IREJ_LAST, P_XSTORE
      save    COUNT_ITER_BLOCK, IREJ_LAST, P_XSTORE

      integer P_DXSTORE, P_DVLSTORE, P_DVUSTORE, P_LAMSTORE, SOC_ITER
      save    P_DXSTORE, P_DVLSTORE, P_DVUSTORE, P_LAMSTORE, SOC_ITER

      integer P_WXSTORE, P_WVLSTORE, P_WVUSTORE, P_PZSTORE, P_CSTORE
      save    P_WXSTORE, P_WVLSTORE, P_WVUSTORE, P_PZSTORE, P_CSTORE

      integer P_SLSTORE, P_SUSTORE, P_WLAMOLDSTORE, P_WCSTORE
      save    P_SLSTORE, P_SUSTORE, P_WLAMOLDSTORE, P_WCSTORE

      integer P_WDXSTORE, P_WDVLSTORE, P_WDVUSTORE, P_WLAMSTORE
      save    P_WDXSTORE, P_WDVLSTORE, P_WDVUSTORE, P_WLAMSTORE

      double precision CNRM_NEW_STORE, PHI0, GD, CL, CDL
      save             CNRM_NEW_STORE, PHI0, GD, CL, CDL

      double precision LAST_MU, CNRM_MAX, MACHEPS, MACHTINY, ALPHA_MIN
      save             LAST_MU, CNRM_MAX, MACHEPS, MACHTINY, ALPHA_MIN
      double precision CNRM_SML, ALPHA_PRIMAL, DHD_STORE
      save             CNRM_SML, ALPHA_PRIMAL, DHD_STORE

      double precision ALPHA_STORE, ALPHA_DUAL_STORE
      save             ALPHA_STORE, ALPHA_DUAL_STORE

      double precision WALPHA_STORE, WALPHA_DUAL_STORE
      save             WALPHA_STORE, WALPHA_DUAL_STORE
      double precision WALPHA_MIN_STORE
      save             WALPHA_MIN_STORE

      logical SMALL_STEP_BEFORE, TRY_WATCH
      save    SMALL_STEP_BEFORE, TRY_WATCH

      double precision BACKTRACK_FAC
      parameter( BACKTRACK_FAC = 5.d-1 )

      double precision ETA_F
      double precision ALPHA_MIN_FRAC
c      parameter( DELTA       = 1.d-4 )
      parameter( ETA_F       = 1.d-4 )

#ifdef INCLUDE_CC
CMPECSTART
      double precision MPEC_ETA
      common /MPEC/ MPEC_ETA
      save /MPEC/
CMPECEND
#endif

CTODO out?
C      integer QRESETSLACKS
C      parameter( QRESETSLACKS = 0 )
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Initialize Storage space
C
      if( ITER.eq.-1 ) then

         if( M.eq.0 ) QSOC = 0

         NFILTER = 0
         P_FILTER_PHI = LRS
         P_FILTER_C   = P_FILTER_PHI + QMAXITER ! That's actually generous ;)
         P_XSTORE     = P_FILTER_C   + QMAXITER
         P_SLSTORE    = P_XSTORE     + N
         P_SUSTORE    = P_SLSTORE    + NLB
         LRS          = P_SUSTORE    + NUB

         if( QWATCHTOL.gt.0.d0 .or. QITERBLOCKMAX.gt.0 ) then
C
C     Reserve storage for watchdog
C
            P_WXSTORE   = LRS
            P_WVLSTORE  = P_WXSTORE  + N
            P_WVUSTORE  = P_WVLSTORE  + NLB
            P_WDXSTORE  = P_WVUSTORE  + NUB
            P_WDVLSTORE = P_WDXSTORE  + N
            P_WDVUSTORE = P_WDVLSTORE + NLB
            P_WCSTORE   = P_WDVUSTORE + NUB
            LRS         = P_WCSTORE   + M
            if( abs(QLAMBDA).eq.2 ) then
               P_WLAMSTORE    = LRS
               P_WLAMOLDSTORE = P_WLAMSTORE + M
               LRS            = P_WLAMOLDSTORE + M
            endif
         endif

         if( QSOC.ne.0 ) then
C
C     Reserve storage space for second order correction option
C
            P_DXSTORE  = LRS
            P_DVLSTORE = P_DXSTORE  + N
            P_DVUSTORE = P_DVLSTORE + NLB
            P_CSTORE   = P_DVUSTORE + NUB
            LRS        = P_CSTORE   + M
            if( abs(QLAMBDA).eq.2 ) then
               P_LAMSTORE = LRS
               LRS        = P_LAMSTORE + M
            endif
            if( QQUASI.ne.0 .and. abs(QQUASI).lt.6 ) then
               P_PZSTORE = LRS
               LRS       = P_PZSTORE + NIND
            endif

         endif

         LAST_MU  = -1d0
         MACHEPS  = D1MACH(4)
         MACHTINY = sqrt(D1MACH(1))
CORIG         MACHTINY = sqrt(sqrt(sqrt(D1MACH(1))))

         WFLAG            = 0
         COUNT_ITER_BLOCK = 0
         RESTO_KKT_FLAG   = 0

         SMALL_STEP_BEFORE = .false.
         TRY_WATCH         = .false.

         goto 9999
      endif

      if( MEMDBG .and. QPRINT.ge.3 ) then
         write(line,1)'filter', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif
C
C     Initialization
C
      if( abs(QMERIT).eq.5 .and. abs(QLAMBDA).ne.2 ) then
         IERR = 4
         call C_OUT(2,0,1,
     1        'filter:  Need QLAMBDA = 2 for this filter option.')
         goto 9999
      endif

 5    continue
      IERR = 0
      p_rwend = 0
      p_iwend = 0

      C_ACCEPT = ' '
C
C     If in restoration phase, go directly there
C
      if( RESTO.gt.0 ) then
         CNRM0 = CALC_NRM(M, C)
         LS_COUNT = 0
         COUNT_ITER_BLOCK = 0
         goto 3015
      endif
C
C     If mu changed, delete all entries in filter
C
      if( LAST_MU.ne.MU ) then
         NFILTER = 0
         LAST_MU = MU
         SMALL_STEP_BEFORE = .false.
         WFLAG = 0
      endif

      if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
         write(line,7) WFLAG, SOC_FLAG
 7       format(/,'====> Starting filter line search for WFLAG = ',i1,
     1        ' and SOC_FLAG = ',i1,' <=====',/)
         call C_OUT(1,2,3,line)
      endif
C
C     If in trial step for watchdog, use values for previous point
C
      if( WFLAG.eq.1 ) then
         DHD = DHD_STORE
         goto 6
      elseif( WFLAG.eq.2 ) then
         WFLAG = 0
      endif
C
C     The following calculations only have to be done once per iteration
C
      if( QSOC.ne.0 .and. SOC_FLAG.eq.1 ) goto 6
C
C     Check if step is too small (close to solution already?) so that we can
C     quit with a warning message?
C
      reldxnrm = 0.d0
      do i = 1, N
         reldxnrm = max(reldxnrm, abs(DX(i))/(1.d0+abs(X(i))))
      enddo
CTODO
C      if( reldxnrm.lt.QPRECFACT*MACHEPS .and. ALPHA.gt.0.d0 ) then
      if( reldxnrm.lt.QPRECFACT*MACHEPS .and. ALPHA.gt.0.d0 .and.
     1     REGU.le.1d0) then
         write(line,2) reldxnrm
 2       format(' filter: Step is not making much progress:',/,
     1        '         max(|dx|/|1+x|) = ',d14.6)
         call C_OUT(1,0,2,line)
         LS_COUNT = 0
         NU_OUT   = 0d0
         call DAXPY(N  , ALPHA     , DX  , 1, X  , 1)
         do i = 1, NLB
            k = ILB(i)
            S_L(i) = dmax1(MACHTINY, X(k) - BNDS_L(i))
            X(k) = S_L(i) + BNDS_L(i)
         enddo
         do i = 1, NUB
            k = IUB(i)
            S_U(i) = dmax1(MACHTINY, BNDS_U(i) - X(k))
            X(k) = BNDS_U(i) - S_U(i)
         enddo
         call GET_F(N, X, IVAR, NORIG, XORIG, M, CSCALE, NLB, ILB,
     1        S_L, NUB, IUB, S_U, MU, F, LIW-p_iwend, IW(p_iwend+1),
     1        IERR, EV_F, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'filter: get_f(s) returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         call GET_C(ITER, N, NIND, X, IVAR, NORIG,
     1        XORIG, M, CSCALE, C, KCONSTR, LRS, RS,
     2        LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         call DAXPY(NLB, ALPHA_DUAL, DV_L, 1, V_L, 1)
         call DAXPY(NUB, ALPHA_DUAL, DV_U, 1, V_U, 1)
         if( abs(QLAMBDA).eq.2 ) then
            call DSCAL(M, ALPHA, LAM, 1)
            call DAXPY(M, (1.d0-ALPHA), LAMOLD, 1, LAM, 1)
            call DCOPY(M, LAM, 1, LAMOLD, 1)
         endif
         if( SMALL_STEP_BEFORE ) then
            C_ACCEPT = 's'
            IERR = -2
         else
            C_ACCEPT = 'S'
         endif
         SMALL_STEP_BEFORE = .true.
         goto 9999
      else
         SMALL_STEP_BEFORE = .false.
      endif
C
      IREJ_LAST = 0
      SOC_ITER = 0
      ALPHA_PRIMAL = ALPHA
C
C     Compute current THETA (in CNRM0)
C
      CNRM0 = CALC_NRM(M, C)
C
C     If in KKT restoration phase, go exactly there
C
      if( RESTO_KKT_FLAG.ne.0 ) then
         LS_COUNT = 1
         goto 3010
      endif
C
C     Get value of barrier function at current point
C
      bar0 = CALC_BAR(NLB, NUB, S_L, S_U, V_L, V_U, MU)
      PHI0 = F + bar0
C
C     For Lagrangian add contraint times multiplier part
C
      if( abs(QMERIT).eq.5 ) then
         PHI0 = PHI0 + DDOT(M, C, 1, LAMOLD, 1)
      endif
C
C     Determine maximal allowed constrained violation
C
      if( ITER.eq.0 ) then
         CNRM_MAX = QTHETAMAXFACT*dmax1(1d0, CNRM0)
         CNRM_SML = QTHETAMINFACT*dmax1(1d0, CNRM0)
      endif
C
C     Get linear prediction in barrier function g^T d
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
#ifndef INCLUDE_CC
         RW(p_gb+k) = RW(p_gb+k) - MU/S_L(i)
#else
CMPECSTART
         if( QMPEC_TRIGGER.ge.1 ) then
            RW(p_gb+k) = RW(p_gb+k) - V_L(i) - 
     1           (MU - S_L(i)*V_L(i))/(S_L(i)+MPEC_ETA*V_L(i))
         else
            RW(p_gb+k) = RW(p_gb+k) - MU/S_L(i)
         endif
CMPECEND
#endif
      enddo
      do i = 1, NUB
         k = IUB(i)
#ifndef INCLUDE_CC
         RW(p_gb+k) = RW(p_gb+k) + MU/S_U(i)
#else
CMPECSTART
         if( QMPEC_TRIGGER.ge.1 ) then
            RW(p_gb+k) = RW(p_gb+k) + V_U(i) +
     1           (MU - S_U(i)*V_U(i))/(S_U(i)+MPEC_ETA*V_U(i))
         else
            RW(p_gb+k) = RW(p_gb+k) + MU/S_U(i)
         endif
CMPECEND
#endif    
      enddo
      GD = DDOT(N, DX, 1, RW(p_gb+1), 1)
      p_rwend = p_gb
C
C     For Lagranian filter option computer CL and CDL
C
      if( abs(QMERIT).eq.5 ) then
         CL  = DDOT(M, C, 1, LAMOLD, 1)
         CDL = DDOT(M, C, 1, LAM   , 1) - CL
      endif
C
C     Determine minimal step size ALPHA_MIN
C
      if( CNRM0.gt.0d0 ) then
         ALPHA_MIN = QGAMMA_THETA
      else
         ALPHA_MIN = QGAMMA_THETA ! 1d300 ?
CTRY
         if( GD.gt.0d0 ) then
            if( SOC_FLAG.ne.1 ) then
               gnrm = DNRM2(N, G, 1)
               if( GD.gt.QPRECFACT*MACHEPS*gnrm ) then
                  write(line,*)
     1                 'filter: WARNING:  CNRM0 = 0, but GD = ',GD
                  call C_OUT(2,0,1,line)
c                  IERR = 581
c                  goto 9999
               endif
            endif
            GD = -sqrt(MACHTINY)
            ALPHA_MIN = ALPHA
         endif
      endif
C
      if( abs(QMERIT).eq.4 ) then
         if( GD.lt.0d0 ) then
            ALPHA_MIN = min(ALPHA_MIN, -QGAMMA_F*CNRM0/GD)
            if( CNRM0.le.CNRM_SML ) then
               ALPHA_MIN = min(ALPHA_MIN,
     1              QDELTA*(CNRM0**QS_THETA)/((-GD)**QS_F))
            endif
         endif
      elseif( abs(QMERIT).eq.5 ) then
CTODO I have no clue what rule shoule be here...
c$$$         a = GD - CL + CDL
c$$$         b = CDL
c$$$         if( a.lt.0.d0 ) then
c$$$            ALPHA_MIN = min(ALPHA_MIN, -QGAMMA_F*CNRM0/a)
c$$$         endif
c$$$         if( b.gt.0.d0 ) then
c$$$            ALPHA_MIN = min(ALPHA_MIN,
c$$$     1           (a+(QDELTA*(CNRM0**QS_THETA)/ALPHA_MIN)**(-QS_F))/b)
c$$$         elseif( a.lt.0.d0 ) then
c$$$            ALPHA_MIN = min(ALPHA_MIN,
c$$$     1           QDELTA*(CNRM0**QS_THETA)/((-a)**QS_F))
c$$$         endif
         dltd = GD - CL
         if( CDL.eq.0.d0 ) then
            if( dltd .lt. 0.d0 ) then
               ALPHA_MIN = min(ALPHA_MIN,
     1              min(QGAMMA_F*CNRM0,
     2              QDELTA*(CNRM0**QS_THETA))/((-dltd)**QS_F))
            endif
         else
            tmp1 = (1+dltd/CDL)/2.d0
            tmp2 = QDELTA*(CNRM0**QS_THETA)/(CDL**QS_F)
            tmp3 = QGAMMA_F*CNRM0/CDL
            if( CDL.gt.0.d0 ) then
               tmp4 = tmp1 + dsqrt(tmp1**2+tmp2)
               if( tmp4.gt.0.d0 ) ALPHA_MIN = min(ALPHA_MIN, tmp4)
                                ! in extremely ill-conditioned cases, tmp4
                                ! can become zero numerically
               tmp4 = tmp1 + dsqrt(tmp1**2+tmp3)
               if( tmp4.gt.0.d0 ) ALPHA_MIN = min(ALPHA_MIN, tmp4)
            else
               tmp4 = tmp1 - dsqrt(tmp1**2+tmp2)
               if( tmp4.gt.0.d0 ) then
                  ALPHA_MIN = min(ALPHA_MIN, tmp4)
               endif
               tmp4 = tmp1 - dsqrt(tmp1**2+tmp3)
               if( tmp4.gt.0.d0 ) then
                  ALPHA_MIN = min(ALPHA_MIN, tmp4)
               endif
            endif
         endif
      endif
C
      ALPHA_MIN = QALPHAMINFRAC*ALPHA_MIN
      ALPHA_MIN = max(ALPHA_MIN, QPRECFACT*MACHEPS/reldxnrm)
      if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
         write(line,8) ALPHA_MIN
 8       format('- minimal step size ALPHA_MIN = ',d10.4)
         call C_OUT(1,3,1,line)
      endif
      if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
         write(line,7000)
 7000    format(/,'   #     alpha           phi(x_k)           phi(x_k+a
     1lpha d_k)           |c(x_k)|           |c(x_k+alpha d_k)|        m
     2odel(alpha)      SC Rej')
         call C_OUT(1,2,2,line)
      endif
C
C     reserve work space for trial points
C
 6    p_xnew    = p_rwend
      p_slnew   = p_xnew   + N
      p_sunew   = p_slnew  + NLB
      p_cnew    = p_sunew  + NUB
      p_rwend   = p_cnew   + M
      if( abs(QMERIT).eq.5 ) then
         p_lamnew = p_rwend
         p_rwend  = p_lamnew + M
      endif
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
C     Now do the actual line search
C
      LS_COUNT = 0
      cnrm_new = 0d0
      phi_new  = 0d0
      f_new    = 0d0
      ftype = .false.
      model_alpha = 0.d0
C
C ==============================================================================
C
C     BEGINNING OF BACKTRACKING LINE SEARCH LOOP
C
C ==============================================================================
C
 10   continue

      if( WFLAG.eq.1 .and. ALPHA.eq.0.d0 ) goto 900
CDELETEME
c      if( .true. ) goto 3000 ! jump directly to restoration
c      if( ITER.eq.0 ) goto 3000 ! jump directly to restoration
c      if( ITER.ge.3 ) goto 3000 ! jump directly to restoration
      if( ALPHA.eq.0.d0 ) goto 3000 ! jump directly to restoration
                                ! (ALPHA = 0 set in get_step_full or get_pz)

      LS_COUNT = LS_COUNT + 1
C
C     Compute new trial point and evaluate objective and constraint
C     function at that point
C
      call DCOPY(N, X, 1, RW(p_xnew+1), 1)
      call DAXPY(N, ALPHA, DX, 1, RW(p_xnew+1), 1)
      if( abs(QMERIT).eq.5 ) then
         call DCOPY(M, LAM, 1, RW(p_lamnew+1), 1)
         call DSCAL(M, ALPHA, RW(p_lamnew+1), 1)
         call DAXPY(M, (1.d0-ALPHA), LAMOLD, 1, RW(p_lamnew+1), 1)
      endif
C
C     Compute new trial slack variables
C     (If slacks get too small, the problem probably doesn't have a
C      feasible strictly interior point)
C
C     A sophisticated user might want to use a different way to handle this
C     situation.  In that case, the flag '-DUSE_IP_MYMOVEBNDS' has to be
C     used at compilation time.
C
#ifdef USE_IP_MYMOVEBNDS
      call IP_MYMOVEBNDS(N, RW(p_xnew+1), NLB, RW(p_slnew+1), BNDS_L,
     1     ILB, V_L, DV_L, NUB, RW(p_sunew+1), BNDS_U, IUB, V_U, DV_U, 
     1     MU, ITER, QPRECFACT, MACHEPS)
#else
      call CHECK_SLACKS(N, RW(p_xnew+1), NLB, ILB, S_L, V_L, BNDS_L,
     1     RW(p_slnew+1), NUB, IUB, S_U, V_U, BNDS_U, RW(p_sunew+1),
     1     MU, MACHEPS)
#endif
C
      call GET_F(N, RW(p_xnew+1), IVAR, NORIG, XORIG, M, CSCALE, NLB,
     1     ILB, RW(p_slnew+1), NUB, IUB, RW(p_sunew+1), MU, f_new,
     2     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'filter: get_f returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( FFINITE(f_new).eq.0 ) then ! we don't want NaN here; cut step!
         call C_OUT(2,0,1,' --- WARNING: Evaluation of'//
     1        ' objective function produced Nan/Inf.')
         irej = -2
         goto 800
      endif
C
      call GET_C(ITER, N, NIND, RW(p_xnew+1), IVAR, NORIG, XORIG, M,
     1     CSCALE, RW(p_cnew+1), KCONSTR, LRS, RS, LIS, IS,
     2     LRW-p_rwend, RW(p_rwend+1),
     4     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1        'filter: Error: get_c returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.ne.0 ) then
         write(line,*)
     1        'filter: Warning: get_c returns IERR = ',IERR
         call C_OUT(2,0,1,line)
      endif
C
      cnrm_new = CALC_NRM(M, RW(p_cnew+1))
      if( FFINITE(cnrm_new).eq.0 ) then ! we don't want NaN here; cut step!
         call C_OUT(2,0,1,' --- WARNING: Evaluation of constraint'//
     1        ' functions produced Nan/Inf.')
         irej = -2
         goto 800
      endif
C
C     Check if constraint violation is 'exploding'
C
      if( cnrm_new.gt.CNRM_MAX ) then
         if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
            write(line,7001) cnrm_new
 7001       format('REJECTED: cnrm = ',d10.4,' is getting too large.')
            call C_OUT(1,2,1,line)
         endif
         irej = -3
         goto 800
      endif
C
C     Compute new value of barrier (or Lagrangian) function
C
      bar_n = CALC_BAR(NLB, NUB, RW(p_slnew+1),
     1     RW(p_sunew+1), dummy, dummy, MU)
      phi_new = f_new + bar_n
      if( abs(QMERIT).eq.5 ) then
         phi_new = phi_new + DDOT(M, RW(p_cnew+1), 1, RW(p_lamnew+1), 1)
      endif
      if( FFINITE(phi_new).eq.0 ) then ! we don't want NaN here; cut step!
         call C_OUT(2,0,1,' --- WARNING: Evaluation of barrier'//
     1        ' objective function produced Nan/Inf.')
         irej = -2
         goto 800
      endif
C
C     Check if trial point acceptable to current point
C
CORIG      if( WFLAG.ne.1 ) then
      if( WFLAG.eq.1 ) then
         alpha_tmp = WALPHA_STORE
      elseif( SOC_FLAG.eq.1 ) then
         alpha_tmp = ALPHA_STORE
      else
         alpha_tmp = ALPHA
      endif
      model_alpha = FILTER_MODEL(alpha_tmp, GD, CL, CDL, DHD)
      if( model_alpha.ge.0.d0 ) then
         ftype = .false.
      else
         ftype = ( ((-model_alpha)**QS_F)*(alpha_tmp**(1.d0-QS_F))
     1        .gt.QDELTA*(CNRM0**QS_THETA) .and.
     1        CNRM0.le.CNRM_SML ) ! if true, enforce Armijo cond for obj func
      endif
      if( .false. .and. QCNR.gt.0 .and. QPRINT.ge.2 ) then
         if( ftype ) then
            call C_OUT(1,2,1,
     1  'Switching condition satisfied -> enforce f-type decrease.')
         else
            call C_OUT(1,2,1,
     1  'Switching condition NOT satisfied -> enforce h-type decrease.')
         endif
         write(line,7002) model_alpha, CNRM0
 7002    format('    model(alpha) = ',d10.4,'   ||c(x_k)|| = ',d10.4)
         call C_OUT(1,3,1,line)
      endif

      if( ftype ) then
C
C     Check if enough decrease in barrier function
C
         lhs = PHI0-phi_new + QPRECFACT*MACHEPS*dabs(PHI0)
         rhs = -ETA_F*model_alpha
         if( lhs.lt.rhs ) then
            irej = -1
            goto 800
         endif
C         write(line,*) 'f-type iteration test PASSED:'
C         call C_OUT(1,2,1,line)
C         write(line,*) 'phi0 = ',PHI0,' phi_new = ',phi_new
C         call C_OUT(1,2,1,line)
C         write(line,*) '-ETA_F*model_alpha = ',rhs,
C     1        ' PHI0-phi_new = ',lhs
C         call C_OUT(1,2,1,line)

      else
C
C     Check if "filter margin" satisfied at current iterate
C
         fi = PHI0
         ci = CNRM0
         if( cnrm_new-(ci-QGAMMA_THETA*ci).gt.QPRECFACT*MACHEPS*ci .and.
     1        phi_new-(fi-QGAMMA_F*ci).gt.QPRECFACT*MACHEPS*dabs(fi) )
     2        then
            irej = 0
            goto 800
         endif

      endif
C
C     Check if trial point is acceptable to the filter
C
      irej = FILTER_CHECK(NFILTER, RS(P_FILTER_C+1), RS(P_FILTER_PHI+1),
     1     cnrm_new, phi_new, QPRECFACT, MACHEPS)
      if( irej.gt.0 ) goto 800
C
C     OK, point is accepted!
C
      goto 1000

C ==============================================================================
C
C     Not accepted, maybe do some second order correction.
C     Otherwise, cut step size and go back to beginning of loop.
C
C ==============================================================================

 800  continue
      if( .false. .and. QCNR.gt.0 .and. QPRINT.gt.1 ) then
         call C_OUT(1,2,1,'----- Point rejected -----')
         if( irej.eq.0 ) then
            call C_OUT(1,2,1,
     1           'No sufficient decrease w.r.t. current point.')
            write(line,*) 'PHI0    = ',PHI0,
     1           ' C     = ',CNRM0
            call C_OUT(1,2,1,line)
            write(line,*) 'PHI     = ',phi_new,
     1           ' C     = ',cnrm_new
            call C_OUT(1,2,1,line)
         elseif( irej.gt.0 ) then
            write(line,*) 'Filter entry that rejected: ',irej
            call C_OUT(1,2,1,line)
            write(line,*) 'PHI_filter = ',RS(P_FILTER_PHI+irej),
     1           ' C_filter = ',RS(P_FILTER_C+irej)
            call C_OUT(1,2,1,line)
            write(line,*) 'PHI        = ',phi_new,
     1           ' C        = ',cnrm_new
            call C_OUT(1,2,1,line)
         elseif( irej.eq.-1 ) then
            write(line,*) 'f-type iteration test:'
            call C_OUT(1,2,1,line)
            write(line,*) 'phi0 = ',PHI0,' phi_new = ',phi_new
            call C_OUT(1,2,1,line)
            write(line,*) '-ETA_F*model_alpha = ',rhs,
     1           ' PHI0-phi_new = ',lhs
            call C_OUT(1,2,1,line)
         elseif( irej.eq.-2 ) then
            write(line,*) 'Nan occured...'
            call C_OUT(1,2,1,line)
         else
c            write(line,*) 'irej = ',irej
c            call C_OUT(1,2,1,line)
         endif
      endif
      if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
         write(line,7010) LS_COUNT, ALPHA, PHI0, phi_new, CNRM0,
     1        cnrm_new, model_alpha, ftype, irej
 7010    format(i4,d12.4,5d24.16,l2,i4)
         call C_OUT(1,2,1,line)
      endif
C
C     For Second Order Correction: Check, if SOC needed
C
C      if( WFLAG.eq.1 ) then
c         ctmp = CNRM_NEW_STORE
C      else
         ctmp = CNRM0
C      endif
      if( QSOC.ne.0 .and. SOC_FLAG.eq.0 .and. irej.ne.-2 .and.
     1     (ctmp.le.cnrm_new) .and. WFLAG.ne.2 ) then
C
C     Try SOC
C
         write(line,*) 'Try second order correction with CNRM0 = ',CNRM0
         call C_OUT(1,2,1,line)
         write(line,*) '   and cnrm_new = ', cnrm_new
         call C_OUT(1,2,1,line)
         write(line,*) '   and ALPHA = ', ALPHA
         call C_OUT(1,2,1,line)
C
C     Store current step etc
C
         call DCOPY(N,   DX,   1, RS(P_DXSTORE +1), 1)
         call DCOPY(NLB, DV_L, 1, RS(P_DVLSTORE+1), 1)
         call DCOPY(NUB, DV_U, 1, RS(P_DVUSTORE+1), 1)
         call DCOPY(M  , C   , 1, RS(P_CSTORE  +1), 1)
         if( abs(QLAMBDA).eq.2 ) then
            call DCOPY(M, LAM, 1, RS(P_LAMSTORE+1), 1)
            call DCOPY(M, LAMOLD, 1, LAM, 1)
         endif
         if( QQUASI.ne.0 .and. abs(QQUASI).lt.6 ) then
            call DCOPY(NIND, PZ, 1, RS(P_PZSTORE+1), 1)
         endif
         ALPHA_STORE      = ALPHA
         ALPHA_DUAL_STORE = ALPHA_DUAL
         CNRM_NEW_STORE   = cnrm_new
C
C     Compute new right hand side for constraints
C
         call DSCAL(M, ALPHA, C, 1)
         call DAXPY(M, 1.d0, RW(p_cnew+1), 1, C, 1)
C
         SOC_FLAG = 1
         SOC_ITER = 1

         goto 9999

      elseif( QSOC.ne.0 .and. SOC_FLAG.eq.1 ) then

         if( cnrm_new.ge.0.99d0*CNRM_NEW_STORE .or.
     1        SOC_ITER.ge.QSOC .or. irej.eq.-2 ) then
C
C     SOC didn't help
C
            if( cnrm_new.ge.0.99d0*CNRM_NEW_STORE ) then
               write(line,*)
     1              'Second order correction rejected with cnrm_new = ',
     2              cnrm_new
               call C_OUT(1,2,1,line)
            elseif( SOC_ITER.ge.QSOC ) then
               write(line,*) 'Tried ',SOC_ITER,
     1              ' second order corrections. Quit SOC.'
               call C_OUT(1,2,1,line)
            endif
C
C     Restore old search direction
C
            call DCOPY(N,   RS(P_DXSTORE +1), 1, DX,   1)
            call DCOPY(NLB, RS(P_DVLSTORE+1), 1, DV_L, 1)
            call DCOPY(NUB, RS(P_DVUSTORE+1), 1, DV_U, 1)
            call DCOPY(M  , RS(P_CSTORE  +1), 1, C   , 1)
            if( abs(QLAMBDA).eq.2 ) then
               call DCOPY(M, RS(P_LAMSTORE+1), 1, LAM, 1)
            endif
            if( QQUASI.ne.0 .and. abs(QQUASI).lt.6 ) then
               call DCOPY(NIND, RS(P_PZSTORE+1), 1, PZ, 1)
            endif
            ALPHA      = ALPHA_STORE
            ALPHA_DUAL = ALPHA_DUAL_STORE
C
C     Now do back tracking line search with original direction
C
            SOC_FLAG = 2
C
C     Since we already had tried the full step, let's start this time with
C     the reduced trial step size
C
CBLOCK            if( CNRM0.ge.QWATCHTOL ) then
            if( CNRM0.ge.QWATCHTOL .and. .not.TRY_WATCH ) then
               ALPHA = BACKTRACK_FAC*ALPHA
               goto 10
            else
               goto 10
            endif

         else
C
C     Try additional second order corrections
C
            write(line,*)
     1        'Try additional second order correction with cnrm_new = ',
     2           cnrm_new
            call C_OUT(1,2,1,line)
            write(line,*) '   and ALPHA = ', ALPHA
            call C_OUT(1,2,1,line)

            CNRM_NEW_STORE = cnrm_new
            if( abs(QLAMBDA).eq.2 ) then
               call DCOPY(M, LAMOLD, 1, LAM, 1)
            endif
C
C     Compute new right hand side for constraints
C
            call DSCAL(M, ALPHA, C, 1)
            call DAXPY(M, 1.d0, RW(p_cnew+1), 1, C, 1)
C
            SOC_ITER = SOC_ITER + 1

            goto 9999

         endif

      endif
C
C     Watchdog
C
 900  continue
      if( WFLAG.eq.0 .and. irej.ne.-2 .and.
     1     (CNRM0.lt.QWATCHTOL.and.M.gt.0  .or. TRY_WATCH) ) then
C
C     Watchdog currently not active - wake it up
C
         call C_OUT(1,1,1,'filter: Try watchdog step.')
         WFLAG = 1
         call DCOPY(N,   X,    1, RS(P_WXSTORE +1), 1)
         call DCOPY(NLB, V_L,  1, RS(P_WVLSTORE +1), 1)
         call DCOPY(NUB, V_U,  1, RS(P_WVUSTORE +1), 1)
         call DCOPY(N,   DX,   1, RS(P_WDXSTORE +1), 1)
         call DCOPY(NLB, DV_L, 1, RS(P_WDVLSTORE+1), 1)
         call DCOPY(NUB, DV_U, 1, RS(P_WDVUSTORE+1), 1)
         call DCOPY(M,   C,    1, RS(P_WCSTORE+1),   1)
         if( abs(QLAMBDA).eq.2 ) then
            call DCOPY(M, LAM, 1, RS(P_WLAMSTORE+1), 1)
            call DCOPY(M, LAMOLD, 1, RS(P_WLAMOLDSTORE+1), 1)
         endif
         WALPHA_MIN_STORE  = ALPHA_MIN
         WALPHA_STORE      = ALPHA
         WALPHA_DUAL_STORE = ALPHA_DUAL
         C_ACCEPT         = 'w'
         CNRM_NEW_STORE   = cnrm_new
         DHD_STORE        = DHD
         goto 1010

      elseif( WFLAG.eq.1 ) then
C
C     Trial step didn't help, restore old step
C
         call C_OUT(1,1,1,
     1 'filter: Watchdog trial step didn''t help. Restore old iterate.')
         WFLAG = 2
         call DCOPY(N,   RS(P_WXSTORE +1),  1, X,    1)
         call DCOPY(NLB, RS(P_WVLSTORE+1),  1, V_L,  1)
         call DCOPY(NUB, RS(P_WVUSTORE+1),  1, V_U,  1)
         call DCOPY(N,   RS(P_WDXSTORE +1), 1, DX,   1)
         call DCOPY(NLB, RS(P_WDVLSTORE+1), 1, DV_L, 1)
         call DCOPY(NUB, RS(P_WDVUSTORE+1), 1, DV_U, 1)
         call DCOPY(M,   RS(P_WCSTORE+1),   1, C,    1)
         if( abs(QLAMBDA).eq.2 ) then
            call DCOPY(M, RS(P_WLAMSTORE+1), 1, LAM, 1)
            call DCOPY(M, RS(P_WLAMOLDSTORE+1), 1, LAMOLD, 1)
         endif
         ALPHA      = BACKTRACK_FAC * WALPHA_STORE
         ALPHA_DUAL = WALPHA_DUAL_STORE
         ALPHA_MIN  = WALPHA_MIN_STORE
         DHD        = DHD_STORE
         LS_COUNT = 1
         goto 10
      endif
C
C     Cut back step size
C
      ALPHA = BACKTRACK_FAC*ALPHA
      IREJ_LAST = irej
C
      if( ALPHA.lt.ALPHA_MIN ) then
C
C     Step size becomes too small - switch to feasibility restoration
C
         if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
            write(line,*) 'filter: ALPHA         = ',ALPHA
            call C_OUT(1,2,1,line)
         endif
CBLOCK         if( IREJ_LAST.gt.0 ) then
CBLOCK            COUNT_ITER_BLOCK = COUNT_ITER_BLOCK + 1
CBLOCK            if( QITERBLOCKMAX.gt.0 .and.
CBLOCK     1           COUNT_ITER_BLOCK.gt.QITERBLOCKMAX) then
CBLOCK               call C_OUT(2,0,1,
CBLOCK     1              'filter: Seems that filter blocks progress.')
CBLOCK               IERR = 23452
CBLOCK               goto 9999
CBLOCK               if( CNRM_MAX.gt.10d0*CNRM0 ) then
CBLOCK                  write(line,*)
CBLOCK     1 '      Decrease max allowed constraint violation from ',CNRM_MAX
CBLOCK                  call C_OUT(2,0,1,line)
CBLOCK                  CNRM_MAX = CNRM_MAX/10d0
CBLOCK                  write(line,*)
CBLOCK     1               '      to ',CNRM_MAX,' and delete current filter.'
CBLOCK                  call C_OUT(2,0,1,line)
CBLOCK                  NFILTER = 0
CBLOCKCTODO Do we really want to continue with restoration phase here or rather
CBLOCKCTODO do the normal thing?
CBLOCK               else
CBLOCK                  call C_OUT(2,0,1,
CBLOCK     1                 '       However, cannot delete filter.')
CBLOCK               endif
CBLOCK            endif
CBLOCK         else
            COUNT_ITER_BLOCK = 0
CBLOCK         endif
C
C     Can forget about all previously allocated work space
C
         p_rwend = 0
         p_iwend = 0
         goto 3000
      endif

      goto 10

C ==============================================================================
C
C     Ende of loop
C
C ==============================================================================

C ==============================================================================
C
C     Point accepted
C
C ==============================================================================

 1000 continue
      if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
         write(line,7011) LS_COUNT, ALPHA, PHI0, phi_new, CNRM0,
     1        cnrm_new, model_alpha, ftype, '  Ok'
 7011    format(i4,d12.4,5d24.16,l2,a4)
         call C_OUT(1,2,1,line)
      endif

C
C     Check if filter block progress and needs to be deleted
C
CBLOCK      if( IREJ_LAST.gt.0 ) then
      if( IREJ_LAST.ne.-2 .and. WFLAG.eq.0 .and. LS_COUNT.gt.1 ) then
         COUNT_ITER_BLOCK = COUNT_ITER_BLOCK + 1
         if( QITERBLOCKMAX.gt.0 .and.
     1        COUNT_ITER_BLOCK.gt.QITERBLOCKMAX) then
            call C_OUT(2,2,1,
     1           'filter: It is possible that filter blocks progress.')
CWRING            if( IREJ_LAST.gt.0 .and. CNRM_MAX.gt.10d0*CNRM0 ) then
            if( IREJ_LAST.gt.0 .and. CNRM_MAX.gt.10d0*cnrm_new ) then
               write(line,*)
     1'      Max allowed constraint violation is decreased from ',
     2              CNRM_MAX
               call C_OUT(2,2,1,line)
               CNRM_MAX = CNRM_MAX/10d0
               write(line,*)
     1              '       to ',CNRM_MAX,' and current filter deleted.'
               call C_OUT(2,2,1,line)
               NFILTER = 0
               COUNT_ITER_BLOCK = 0
            else
               TRY_WATCH = .true.
               call C_OUT(2,2,1,'Wake up watchdog for next iteration.')
               COUNT_ITER_BLOCK = 0
            endif
CBLOCK            if( CNRM_MAX.gt.10d0*CNRM0 ) then
C               write(line,*)
C     1'      Decrease max allowed constraint violation from ',CNRM_MAX
C               call C_OUT(2,0,1,line)
C               CNRM_MAX = CNRM_MAX/10d0
C               write(line,*)
C     1              '       to ',CNRM_MAX,' and delete current filter.'
C               call C_OUT(2,0,1,line)
C               NFILTER = 0
C            else
C               call C_OUT(2,0,1,
C     1              '      However, cannot delete filter.')
CBLOCK            endif
         endif
      else
         TRY_WATCH = .false.
         COUNT_ITER_BLOCK = 0
      endif
C
C     Check if filter is to be augmented
C
      model_alpha = FILTER_MODEL(ALPHA, GD, CL, CDL, DHD)
      lhs   = phi_new - PHI0 - QPRECFACT*MACHEPS*dabs(PHI0)
      rhs   = ETA_F*model_alpha
      if( model_alpha.ge.0.d0 ) then
         augment = .true.
      else
         augment = ( ((-model_alpha)**QS_F)*(ALPHA**(1.d0-QS_F))
     1        .le.QDELTA*(CNRM0**QS_THETA) .or. lhs.gt.rhs)
      endif
      if( .not.augment ) then
         C_ACCEPT = 'f'
      else
         C_ACCEPT = 'h'
         phi_add  = PHI0 - QGAMMA_F*CNRM0
         cnrm_add = CNRM0 - QGAMMA_THETA*CNRM0
         call FILTER_ADD(NFILTER, RS(P_FILTER_C+1), RS(P_FILTER_PHI+1),
     1        cnrm_add, phi_add)
      endif
C
      if( SOC_FLAG.eq.1 ) then
         SOC_FLAG = 2
         if( augment ) then
            C_ACCEPT = 'H'
         else
            C_ACCEPT = 'F'
         endif
      endif
      if( WFLAG.eq.2 ) then
         if( augment ) then
            C_ACCEPT = 'J'
         else
            C_ACCEPT = 'j'
         endif
      else
         WFLAG = 0
      endif
C
C     Copy trial point to real point
C
 1010 F = f_new
C
      call DCOPY(N,   RW(p_xnew +1), 1, X  , 1)
      call DCOPY(M,   RW(p_cnew +1), 1, C  , 1)
      call DCOPY(NLB, RW(p_slnew+1), 1, S_L, 1)
      call DCOPY(NUB, RW(p_sunew+1), 1, S_U, 1)
      if( abs(QMERIT).eq.5 ) then
         call DCOPY(M, RW(p_lamnew +1), 1, LAM, 1)
         call DCOPY(M, LAM, 1, LAMOLD, 1)
      endif
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
C     Make sure that each V_NEW is at least MACHTINY
C
      do i = 1, NLB
         V_L(i) = dmax1(MACHTINY, V_L(i))
      enddo
      do i = 1, NUB
         V_U(i) = dmax1(MACHTINY, V_U(i))
      enddo
C
C     Free work space
C
      p_rwend = p_xnew
C
C     Take step in LAM
C
CTODO this also for QLAMBDA = 1 ??? Do we really need this here?
 1400 continue
      if( abs(QLAMBDA).eq.2 .and. abs(QMERIT).ne.5 ) then
CTODO decide which option
CORIG         if( ITER.gt.0 ) then
         if( .true. ) then
            call DSCAL(M, ALPHA, LAM, 1)
            call DAXPY(M, 1d0-ALPHA, LAMOLD, 1, LAM, 1)
         endif
         call DCOPY(M, LAM, 1, LAMOLD, 1)
      endif
      NU_OUT = 0d0
C
C     That's it
C
      goto 9999
C
C ==============================================================================
C
C                      FEASIBILITY RESTORATION PHASE
C
C ==============================================================================
C
 3000 continue
C     Reset a few flags
      TRY_WATCH = .false.
      COUNT_ITER_BLOCK = 0

      call C_OUT(1,2,1,'filter: ENTERING RESTORATION PHASE!')
C
      if( abs(QMERIT).eq.5 .and. QRESTO.ne.2 ) then
CTODO CHECK
         IERR = 4
         call C_OUT(2,0,1,'filter: Only TRON for QMERIT = +- 5')
         goto 9999
      endif
C
C     Augment filter
C
      phi_add  = PHI0 - QGAMMA_F*CNRM0
      cnrm_add = CNRM0 - QGAMMA_THETA*CNRM0
      call FILTER_ADD(NFILTER, RS(P_FILTER_C+1), RS(P_FILTER_PHI+1),
     1     cnrm_add, phi_add)
C
C     First try a restoration phase using the usual steps trying only to
C     reduce the primal and dual infeasibility
C
 3010 continue
      if( M.gt.0 .and. QRESTOKKTRED.gt.0.d0 .and. abs(QLAMBDA).eq.2
     1     .and. ALPHA.ne.0.d0 ) then
         call RESTO_KKT(ITER, N, NIND, M, NORIG, XORIG,
     1        CSCALE, IVAR, NFIX, IFIX, NLB, ILB, NUB,
     1        IUB, BNDS_L, BNDS_U, MU, X, CNRM0, NFILTER,
     1        RS(P_FILTER_PHI+1), RS(P_FILTER_C+1),
     1        CNRM_MAX, S_L, S_U, V_L, V_U, F, G, C,
     1        RESTO_KKT_FLAG, LAMOLD, LAM, DX, DV_L, DV_U,
     1        ALPHA, ALPHA_PRIMAL, ALPHA_DUAL,
     1        QPRECFACT, MACHEPS, MACHTINY, NU_OUT, C_ACCEPT,
     1        KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     1        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'filter: RESTO_KKT returned IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         if( RESTO_KKT_FLAG.eq.0 ) then
            RESTO = 0
            goto 9999           ! restoration iteration successful and
                                ! acceptable to filter
         elseif( RESTO_KKT_FLAG.eq.1 ) then
            RESTO = -1
            goto 9999           ! restoration iteration successful, but not yet
                                ! acceptable to filter
         elseif( RESTO_KKT_FLAG.eq.2 ) then
            RESTO = 0
            RESTO_KKT_FLAG = 0  ! need to do 'normal' restoration phase
         else
            write(line,*) 'filter: Invalid RESTO_KKT_FLAG = ',
     1           RESTO_KKT_FLAG
            call C_OUT(2,0,1,line)
            IERR = 943
            goto 9999
         endif
      endif

C     Store current X (for step in dual variables after successful termination
C     of restoration phase
C
      call DCOPY(N, X, 1, RS(P_XSTORE+1), 1)
      call DCOPY(NLB, S_L, 1, RS(P_SLSTORE+1), 1)
      call DCOPY(NUB, S_U, 1, RS(P_SUSTORE+1), 1)

 3015 continue
      C_ALPHA = 'r'
      if( CNRM0.eq.0d0 ) then
         call C_OUT(2,0,1,
     1'filter: current iterate feasible; can''t do restoration. Abort')
         IERR = 16 !562
         goto 9999
      endif
C
C     Since the proposed step was not good, we don't want to accept the
C     freshly computed mutlipliers
C
      if( abs(QLAMBDA).eq.2 .and. QRESTO.ne.1 ) then
CTODO CHECK IF ALSO FOR RESTO_FILTER
         call DCOPY(M, LAMOLD, 1, LAM, 1)
      endif
      if( QRESTO.eq.1 ) then
         call RESTO_FILTER(ITER, N, NIND, M, NORIG, XORIG,
     1        CSCALE, IVAR, NFIX, IFIX, NLB, ILB, NUB,
     1        IUB, BNDS_L, BNDS_U, MU, MUORIG, X, cnrm_new, NFILTER,
     1        RS(P_FILTER_PHI+1), RS(P_FILTER_C+1), CNRM_MAX,
     1        S_L, S_U, V_L, V_U, F, C, DX, DV_L, DV_U, SIGMA_L,
     1        SIGMA_U, RHS_RESTO, SIGMA_PN, ETA_RESTO, QPRECFACT,
     1        MACHEPS, MACHTINY, LAM, LAMOLD, RESTO, C_ACCEPT, ALPHA,
     1        ALPHA_DUAL, LS_COUNT, QDELTA,
     1        QALPHAMINFRAC, ETA_F, BACKTRACK_FAC, B, REGU,
     1        KCONSTR, LRS, RS, LIS, IS,
     1        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     1        IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*)
     1           'filter: Error: resto_filter returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         if( RESTO.eq.0 .and. abs(QLAMBDA).eq.2 ) then
                                ! tell get_step_full to reinitialize LAMBDA
                                ! at next call
            if( QLAMINIT.eq.1 ) then
               call C_OUT(2,3,1,
     1              'filter: Request reinitialization of LAMBDA.')
               INIT_LAM = .true.
            else
               call C_OUT(2,3,1,
     1        'filter: Setting LAMBDA to zero after restoration phase.')
               call DCOPY(M, 0.d0, 0, LAM, 1)
               call DCOPY(M, 0.d0, 0, LAMOLD, 1)
            endif
         endif
         if( RESTO.ne.0 ) goto 9000
      elseif( QRESTO.eq.2 ) then
         RESTO = 0
         COUNT_RESTO_CALL = COUNT_RESTO_CALL + 1
         p_xnew  = p_rwend
         p_cnew  = p_xnew + N
         p_rwend = p_cnew + M
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         IERR = 0
 3020    ierr_store = IERR
         call RESTO_TRON(ITER, N, NIND, M, NORIG, XORIG,
     1        CSCALE, IVAR, NFIX, IFIX, NLB, ILB, NUB,
     1        IUB, BNDS_L, BNDS_U, MU, X, cnrm_new, NFILTER,
     1        RS(P_FILTER_PHI+1), RS(P_FILTER_C+1), CNRM_MAX,
     1        S_L, S_U, V_L, V_U, F, C, QPRECFACT, MACHEPS, MACHTINY,
     1        RW(p_xnew+1), RW(p_cnew+1), LAM,
     1        KCONSTR, LRS, RS, LIS, IS,
     1        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     1        IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.eq.17 .or. IERR.eq.519 ) then
            if( ierr_store.eq.0 ) then
C
C     now solve feasibility problem with strict bounds!
C
               call C_OUT(2,0,1,
     1              'filter: Solve feasibility problem again with '//
     1              'non-relaxed bounds.')
               goto 3020
            endif
            ierr_store = IERR
C
C     Compute ||Ac||
C
            call GET_C(ITER, N, NIND, X, IVAR, NORIG, XORIG,
     1           M, CSCALE, C, KCONSTR, LRS, RS, LIS, IS,
     2           LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.ne.0 ) then
               write(line,*)
     1              'filter: Error: get_c-tron returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            p_tmp   = p_rwend
            p_tmp2  = p_tmp  + N
            p_rwend = p_tmp2 + N
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            call CONSTR(8, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, C, RW(p_tmp+1), i, i,
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.ne.0 ) then
               write(line,*)
     1              'filter: Error: constr8-tron returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            call DCOPY(N, 1.d0, 0, RW(p_tmp2+1), 1)
            do i = 1, NLB
               k = ILB(i)
               RW(p_tmp2+k) = min(RW(p_tmp2+k), S_L(i))
            enddo
            do i = 1, NUB
               k = IUB(i)
               RW(p_tmp2+k) = min(RW(p_tmp2+k), S_U(i))
            enddo
            do i = 1, N
               RW(p_tmp+i) = RW(p_tmp+i) * RW(p_tmp2+i)
            enddo
            cnrm_new = CALC_NRM(N, RW(p_tmp+1))
            write(line,*) '||XAc|| = ',cnrm_new
            call C_OUT(2,0,1,line)
            cnrm_new = CALC_NRM(M, C)
            if( cnrm_new.le.QTOL ) then
               IERR = 16
            else
               IERR = 17
            endif
            goto 9999
         elseif( ierr_store.ne.0 ) then
            IERR = 520
         endif
         if( IERR.gt.0 ) then
            write(line,*)
     1           'filter: Error: resto_tron returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'filter: Warning: resto_tron returns IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
         if( abs(QLAMBDA).eq.2 ) then
                                ! tell get_step_full to reinitialize LAMBDA
                                ! at next call
            if( QLAMINIT.eq.1 ) then
               call C_OUT(2,3,1,
     1              'filter: Request reinitialization of LAMBDA.')
               INIT_LAM = .true.
            else
               call C_OUT(2,3,1,
     1        'filter: Setting LAMBDA to zero after restoration phase.')
               call DCOPY(M, 0.d0, 0, LAM, 1)
               call DCOPY(M, 0.d0, 0, LAMOLD, 1)
            endif
         endif
      else
         write(line,*) 'filter: Invalid QRESTO = ',QRESTO
         call C_OUT(2,0,1,line)
         IERR = 4
         goto 9999
      endif
C
C     Update dual variables
CTODO     one step
C
      p_dl    = p_rwend
      p_du    = p_dl + NLB
      p_rwend = p_du + NUB
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      alphad = 1.d0
      do i = 1, NLB
         k = ILB(i)
         RW(p_dl+i) = MU/(RS(P_XSTORE+k)-BNDS_L(i)) -
     1        SIGMA_L(i)*(X(k)-RS(P_XSTORE+k)) - V_L(i)
         if( RW(p_dl+i).lt.0.d0 ) then
            alphad = dmin1(alphad, QTAU*(-V_L(i)/RW(p_dl+i)) )
         endif
      enddo
      do i = 1, NUB
         k = IUB(i)
         RW(p_du+i) = MU/(BNDS_U(i)-RS(P_XSTORE+k)) +
     1        SIGMA_U(i)*(X(k)-RS(P_XSTORE+k)) - V_U(i)
         if( RW(p_du+i).lt.0.d0 ) then
            alphad = dmin1(alphad, QTAU*(-V_U(i)/RW(p_du+i)) )
         endif
      enddo
CTODO This seems best option (?)
      call DAXPY( NLB, alphad, RW(p_dl+1), 1, V_L, 1)
      call DAXPY( NUB, alphad, RW(p_du+1), 1, V_U, 1)
      write(line,*) 'After Restoration: alpha(dual) = ', alphad
      call C_OUT(1,2,1,line)
cCDELETEME
c      do i = 1, NLB
c         V_L(i) = MU/S_L(i)
c      enddo
c      do i = 1, NUB
c         V_U(i) = MU/S_U(i)
c      enddo

 9000 continue
      NU_OUT = 0d0
c      C_ACCEPT = 'r'
      goto 9999
C
C     THE END
C
 9999 continue
      return
      end
C
C
C
      double precision function FILTER_MODEL(ALPHA, GD, CL, CDL, DHD)
C
C     Compute value of model for f-measure depending on ALPHA
C
      implicit none
      double precision ALPHA, GD, CL, CDL, DHD
      include 'IPOPT.INC'

      if( abs(QMERIT).eq.4 ) then
         FILTER_MODEL = ALPHA * GD
         if( QMERIT.eq.-4 .and. DHD.gt.0.d0 ) then
            FILTER_MODEL = FILTER_MODEL + ALPHA**2 * DHD / 2.d0
         endif
      elseif( abs(QMERIT).eq.5 ) then
         FILTER_MODEL = ALPHA * (GD-CL) +
     1        ALPHA*(1.d0-ALPHA) * CDL
      else
         call C_OUT(2,0,1,'filter_model: Invalid QMERIT. Abort')
         stop
      endif
      return
      end
C
C     The following subroutine computes new trial values for the
C     slacks, and if some become to small, the corresponding bounds
C     are moved a little bit.
C
      subroutine CHECK_SLACKS(N, XNEW, NLB, ILB, S_L, V_L, BNDS_L,
     1     SLNEW, NUB, IUB, S_U, V_U, BNDS_U, SUNEW, MU, MACHEPS)

      implicit none
      include 'IPOPT.INC'
      integer N, NLB, ILB(NLB), NUB, IUB(NUB)
      double precision XNEW(N)
      double precision S_L(NLB), V_L(NLB), BNDS_L(NLB), SLNEW(NLB)
      double precision S_U(NUB), V_U(NUB), BNDS_U(NUB), SUNEW(NUB)
      double precision MU, MACHEPS

      double precision smin, smove, olds, oldb
      integer i, changed
      character*160 line

      changed = 0

      smin = 1.d-16*MU
      smove = MACHEPS**(0.75d0)
      do i = 1, NLB
         SLNEW(i) = XNEW(ILB(i)) - BNDS_L(i)
         if( SLNEW(i).lt.smin ) then
            changed = changed + 1
            olds = SLNEW(i)
            oldb = BNDS_L(i)
            SLNEW(i) = max(smin, MU/V_L(i))
            BNDS_L(i) = XNEW(ILB(i)) - SLNEW(i)
            if( BNDS_L(i).lt.oldb-smove*max(1d0, abs(oldb)) ) then
               BNDS_L(i) = oldb - smove*max(1d0, abs(oldb))
               SLNEW(i) = XNEW(ILB(i)) - BNDS_L(i)
            endif
            if( QCNR.gt.0 .and. QPRINT.ge.4 ) then
               write(line,8000) ILB(i), olds, V_L(i), oldb, BNDS_L(i)
 8000          format('Slack for lower bound',i7,' was ',d10.4,
     1' with multiplier ',d10.4,'. Decrease lower bound from ',d24.17,
     2' to ',d24.17)
               call C_OUT(1,0,1,line)
            endif
         endif
      enddo
      do i = 1, NUB
         SUNEW(i) = BNDS_U(i) - XNEW(IUB(i))
         if( SUNEW(i).lt.smin ) then
            changed = changed + 1
            olds = SUNEW(i)
            oldb = BNDS_U(i)
            SUNEW(i) = max(smin, MU/V_U(i))
            BNDS_U(i) = XNEW(IUB(i)) + SUNEW(i)
            if( BNDS_U(i).gt.oldb+smove*max(1d0, abs(oldb)) ) then
               BNDS_U(i) = oldb + smove*max(1d0, abs(oldb))
               SUNEW(i) = BNDS_U(i) - XNEW(IUB(i))
            endif
            if( QCNR.gt.0 .and. QPRINT.ge.4 ) then
               write(line,8001) IUB(i), olds, V_U(i), oldb, BNDS_U(i)
 8001          format('Slack for upper bound',i7,' was ',d10.4,
     1' with multiplier ',d10.4,'. Increase upper bound from ',d24.17,
     2' to ',d24.17)
               call C_OUT(1,0,1,line)
            endif
         endif
      enddo

      if( changed.gt.0 ) then
         write(line,8002) changed
 8002    format(i7,' bounds have been changed due to small slacks.')
         call C_OUT(2,2,1,line)
      endif
      return
      end


C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine FILTER_WS(N, X, M, NLB, NUB, NZA, LRW, LIW, EV_H,
     1     DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      integer lrw1, liw1, lrw2, liw2
      integer lrw_getc, liw_getc, lrw_getf, liw_getf, liwres, lrwres
      double precision DAT(*), X(N)
      integer IDAT(*)
      external EV_H
      character*80 line

      LRW = 0
      LIW = 0

      call GET_C_WS(N, M, NLB, NUB, NZA, lrw_getc, liw_getc, DAT, IDAT)
      call GET_F_WS(N, M, NLB, NUB, NZA, lrw_getf, liw_getf)

      LRW = max(LRW, lrw_getc, lrw_getf)
      LIW = max(LIW, liw_getc, liw_getf)

      LRW = max(LRW, N)

      lrw1 = N+NLB+NUB+M
      if( abs(QMERIT).eq.5 ) then
         lrw1 = lrw1 + M
      endif
      liw1 = 0
      lrw1 = lrw1 + max(lrw_getc,lrw_getf)
      liw1 = liw1 + max(liw_getc,liw_getf)
      LRW = max(LRW, lrw1)      ! p_xnew
      LIW = max(LIW, liw1)

      lrwres = 0
      liwres = 0
      if( M.gt.0 .and. QRESTOKKTRED.gt.0.d0 .and.
     1     abs(QLAMBDA).eq.2 ) then
         call RESTO_KKT_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         lrwres = max(lrwres, lrw1)
         liwres = max(liwres, liw1)
      endif

      if( QRESTO.eq.1 ) then
         call RESTO_FILTER_WS(N, M, NLB, NUB, NZA, lrw1, liw1,
     1        DAT, IDAT)
         lrwres = max(lrwres, lrw1)
         liwres = max(liwres, liw1)
      elseif( QRESTO.eq.2 ) then
         call RESTO_TRON_WS(N, X, M, NLB, NUB, NZA, lrw1, liw1,
     1        EV_H, DAT, IDAT)
         call CONSTR_WS(N, M, NLB, NUB, NZA, lrw2, liw2, DAT, IDAT)
         lrwres = max(lrwres, N+M+max(lrw_getc,lrw1,2*N+lrw2))
         liwres = max(liwres, liw1, liw_getc, liw2)
      endif
      lrwres = max(lrwres,NLB+NUB)
      LRW = LRW + lrwres
      LIW = LIW + liwres

      if( QPRINT.ge.4 ) then
         write(line,1000)'filter_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
