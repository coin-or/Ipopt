C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine SOLVE_BARRIER(ITER, N, NIND, M, NORIG, XORIG, CSCALE,
     1     IVAR, NFIX, IFIX, NLB, ILB, NUB, IUB, BNDS_L, BNDS_U, MU,
     1     MUORIG, X, S_L, S_U, V_L, V_U, F, G, C, ERR_BAR, ERR_0,
     1     ERR_DINF, ERR_DINF_0, ERR_DINF_US, ERR_PINF, ERR_PINF_0,
     1     ERR_PINF_US, ERR_CMPL_0, ERR_CMPL_US, RV, SIGMA_L, SIGMA_U,
     1     B, W, WCORR, RG, RGOLD, RGB, YPY, RHS_RESTO, SIGMA_PN,
     1     NEWBAS, CONDC,
     1     PZ, ZPZ, LAM, LAMOLD, DX, DV_L, DV_U, GOLD, RESTO, INIT_LAM,
     1     KCONSTR, LRS, LRS_END, RS, LIS, LIS_END, IS,
     1     LRW, RW, LIW, IW, IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV,
     2     EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: solve_barrier.F 655 2004-10-05 17:23:15Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Solve barrier problem (to relaxed accuracy ERR_BAR)
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
CP   ITER      I    INT    iteration counter (if 0, initialize)
CP   N         I    INT    number of variables (without fixed)
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of constraints
CP   NORIG    I/O   INT    total number of variables (incl. fixed vars)
CP   XORIG    I/O   INT    actual iterate
CP                            (original order as in problem statement)
CP   CSCALE    I    DP     scaling factors for constraints
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
CP   MU        I    DP     barrier parameter
CP   MUORIG    I    DP     barrier parameter for original problem
CP                            (might be different in restoration phase)
CP   X        I/O   DP     actual iterate (reordered without fixed vars:
CP                             first M entries belong to dependent
CP                             variables, remaining to independent variables)
CP   S_L      I/O   DP     slacks to lower bounds
CP   S_U      I/O   DP     slacks to upper bounds
CP   V_L      I/O   DP     dual variables for lower bounds
CP   V_U      I/O   DP     dual variables for upper bounds
CP   F        I/O   DP     value of objective function at X
CP                            (at first iteration assumed to be already
CP                             computed)
CP   G        I/O   DP     gradient of objective function
CP                            (at first iteration assumed to be already
CP                             computed)
CP   C        I/O   DP     values of constraints at X
CP                            (at first iteration assumed to be already
CP                             computed)
CP   ERR_BAR   I    DP     required tolerance for barrier problem
CP   ERR_0     O    DP     overall scaled error for original problem
CP   ERR_DINF  O    DP     scaled dual infeasibility
CP   ERR_DINF_0 O   DP     scaled dual infeasibility (with MU = 0 for linear
CP                            damping term)
CP   ERR_DINF_US O  DP     unscaled dual infeasibility
CP   ERR_PINF  O    DP     scaled primal infeasibility
CP   ERR_PINF_0 O   DP     scaled primal infeasibility (with MU=0 for MPCC)
CP   ERR_PINF_US O  DP     unscaled primal infeasibility
CP   ERR_CMPL_0 O   DP     scaled violation of original complementarity
CP   ERR_CMPL_US O  DP     unscaled violation of original complementarity
CP   RV       I/O   DP     reduced (V_U-V_L)
CP                            (needed for computation of KKT error)
CP                            (at first iteration assumed to be already
CP                             computed)
CP   SIGMA_L   S    DP     primal-dual Hessian of lower bound barrier term
CP                            (NLB diagonal elements only)
CP   SIGMA_U   S    DP     primal-dual Hessian of upper bound barrier term
CP                            (NUB diagonal elements only)
CP   B         S    DP     Quasi-Newton estimate of reduced Hessian of
CP                            original NLP
CP                         for CG: containts preconditioner
CP   W         S    DP     reduced Hessian of barrier problem
CP   WCORR     S    DP     correction term for null space step
CP                         (for IP restoration phase: gradient of QP!)
CP   RG        S    DP     reduced gradient of objective function
CP   RGOLD     S    DP     reduced gradient of objective function of
CP                            last iteration
CP                            (for Quasi-Newton update)
CP   RGB       S    DP     reduced gradient of barrier function
CP   YPY      I/S   DP     range space step
CP                            (at first iteration assumed to be already
CP                             computed)
CP   RHS_RESTO I    DP     QRESTO.eq.1: RHS if in IP restoration
CP   SIGMA_PN  I    DP     QRESTO.eq.1: Sigma for P and N variables
CP   NEWBAS    I    LOG
CP   PZ        S    DP     null space step (independent variables)
CP   ZPZ       S    DP     null space step (dependent variables)
CP   LAM       S    DP     Lagrangian multipliers
CP   LAMOLD    S    DP     Lagrangian multipliers from last iteration
CP   DX        S    DP     step for X (primal)
CP   DV_L      S    DP     step for V_L (dual variables for lower bounds)
CP   DV_U      S    DP     step for V_U (dual variables for upper bounds)
CP   GOLD      S    DP     gradient from previous iteration
CP                            (only for QQUASI < 0 )
CP   RESTO    I/O   INT    <>0: we are in restoration phase
CP   INIT_LAM I/O   LOG    Indicates whether equality multipliers are to be
CP                           initialized (might be set to. false. in line
CP                           search) (for get_step_full)
CP   KCONSTR   I    INT    KCONSTR(1): LRS for CONSTR
CP                         KCONSTR(2): P_LRS for CONSTR
CP                         KCONSTR(3): LIS for CONSTR
CP                         KCONSTR(4): P_LIS for CONSTR
CP                         KCONSTR(5): LRW for CONSTR
CP                         KCONSTR(6): LIW for CONSTR
CP   LRS       I    INT    total length of RS
CP   LRS_END  I/O   INT    last used reserved entry in RS
CP   RS       I/O   DP     DP storage space (all!)
CP   LIS       I    INT    total length of IS
CP   LRS_END  I/O   INT    last used reserved entry in IS
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
CS    DCOPY
CS    DNRM2
CS    OPTERROR
CS    GET_LAMBDA
CS    GET_SIGMA
CS    GET_BB
CS    UPDATE_B
CS    UPDATE_B_LM
CS    UPDATE_MPEC_ETA
CS    GET_STEP_FULL
CS    GET_STEP_RED
CS    GET_WCORR
CS    GET_ZPZ
CS    LINESEARCH
CS    GET_G
CS    GET_YPY
CS    GET_RG
CS    GET_RV
CS    ITER_OUT
CS    NEW_BASIS
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
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision F
      double precision G(N)
      double precision C(M)
      double precision ERR_BAR
      double precision ERR_0
      double precision ERR_DINF
      double precision ERR_DINF_0
      double precision ERR_DINF_US
      double precision ERR_PINF
      double precision ERR_PINF_0
      double precision ERR_PINF_US
      double precision ERR_CMPL_0
      double precision ERR_CMPL_US
      double precision RV(NIND)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision B((NIND*(NIND+1))/2)
      double precision W((NIND*(NIND+1))/2)
      double precision WCORR(*)
      double precision RG(NIND)
      double precision RGOLD(NIND)
      double precision RGB(NIND)
      double precision YPY(N)
      double precision RHS_RESTO(N+M)
      double precision SIGMA_PN(M)
      logical NEWBAS
      double precision CONDC
      double precision PZ(NIND)
      double precision ZPZ(M)
      double precision LAM(M)
      double precision LAMOLD(M)
      double precision DX(N)
      double precision DV_L(NLB)
      double precision DV_U(NUB)
      double precision GOLD(N)
      integer RESTO
      logical INIT_LAM
      integer KCONSTR(6)
      integer LRS
      integer LRS_END
      double precision RS(LRS)
      integer LIS
      integer LIS_END
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
      logical SKIP_UPDATE
      double precision     ALPHA
      save    SKIP_UPDATE, ALPHA

      double precision alpha_cut, dummy, nu_out, cnrm
      double precision alpha_dual, eta_resto, fnomu
      integer ls_count, ieigs
      character*1 c_err, c_alpha, c_watch, c_skip
      integer p_rwend, p_iwend, i, idummy, IDAMAX

      integer soc_flag, iter_start
      integer p_rwfull, p_iwfull, prwkeep, piwkeep
      logical small_x_stop, newbas_store, ex, force_newbas
      logical decmu

      double precision DNRM2

      double precision PYNRM_OLD, CNRM0, REGU
      save             PYNRM_OLD, CNRM0, REGU

      double precision ERR, ERR_CMPL, ERR_CMPL_MU
      save             ERR, ERR_CMPL, ERR_CMPL_MU

      integer CHECK_STOP

      character*128 line(3)

      integer ITER_LASTPRINT, WFLAG
      save    ITER_LASTPRINT, WFLAG

      double precision CHECK_XMAX
      parameter( CHECK_XMAX = 1.d50 )

      include 'TIMER.INC'

#ifdef USE_IP_MYBARRIERTEST
      logical IP_MYBARRIERTEST
#endif
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      p_rwend = 0
      p_iwend = 0

      if( MEMDBG ) then
         write(line,1)'solve_barrier', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

C Store current iteration at start of solution of this barrier problem
      iter_start = ITER

C     This flag is set to true if the algorithm terminates because the
C     search direction is getting too small
      small_x_stop = .false.

      MUORIG = MU         ! This is barrier parameter for original problem
                          ! (not restoration phase problem)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                       C
C         START MAIN LOOP FOR BARRIER PROBLEM                           C
C                                                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 100  continue

      force_newbas = .false.

C Initialize second order correction flag
      soc_flag = 0
C
C     Update eta parameter for MPEC algorithm
C
      if( QMPEC_TRIGGER.ge.1 ) then
         call UPDATE_MPEC_ETA(N, M, X, C, LAM, NLB, V_L, NUB, V_U,
     1        MU, RESTO)
      endif
C
C     Compute primal dual Hessian for barrier part SIGMA_L, SIGMA_U
C
      call GET_SIGMA(NLB, NUB, S_L, S_U, V_L, V_U, SIGMA_L, SIGMA_U)
C
C     Compute Error for barrier problem with barrier parameter MU
C
      call OPTERROR(N, NIND, M, X, G, RG, RV, C, CSCALE, NLB, ILB, NUB,
     1     IUB, BNDS_L, BNDS_U, S_L, S_U, V_L, V_U,
     1     LAM, MU, ERR, c_err, ERR_0, ERR_DINF, ERR_DINF_0,
     1     ERR_DINF_US, ERR_PINF, ERR_PINF_0, ERR_PINF_US, ERR_CMPL,
     1     ERR_CMPL_0, ERR_CMPL_US, ITER, IVAR, NFIX, IFIX, NORIG,
     1     XORIG, DX, KCONSTR, LRS, RS, LIS, IS,
     1     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     1     IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     2     EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1        'solve_barrier: Opterror returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.lt.0 ) then
         write(line,*)
     1        'solve_barrier: Opterror returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         IERR = 0
      endif
C
C     Do the output
C
      if( ITER.eq.0 .or. ITER.ne.ITER_LASTPRINT ) then
         if( ITER.eq.0 ) then
            c_skip = ' '
            c_alpha = ' '
            c_watch = ' '
            ALPHA_DUAL = 0.d0
            alpha_cut = 0.d0
            ALPHA = 0.d0
            ls_count = 0
            nu_out = 0.d0
            ieigs = 0
            condc = 0.d0
            REGU = 0.d0
         endif
         cnrm = DNRM2(M, C, 1)
         ITER_LASTPRINT = ITER
         fnomu = F
         call UPDATE_FG_MU(1, N, NLB, ILB, S_L, NUB, IUB, S_U,
     1        MUORIG, 0.d0, fnomu, dummy, LIW, IW, IERR)
         if( IERR.ne.0 ) then
            write(line,*)
     1           'solve_barrier: update_fg_mu returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
C
C     A user might want to his/her own format for printing information
C     about the current iteration.  You may use the following call for
C     that.  Need to compile with '-DUSE_IP_MYITEROUT' in that case.
C
#ifdef USE_IP_MYITEROUT
         call IP_MYITEROUT(ITER, MU, ERR, c_err, ERR_DINF, ERR_PINF,
     1        ERR_CMPL_MU, fnomu, M, C, cnrm, YPY, NIND, PZ, c_skip,
     1        ALPHA_DUAL, alpha_cut, ALPHA, c_alpha, c_watch, ls_count,
     2        nu_out, ieigs, condc, REGU,
     2        RV, RG, RGB, B, SIGMA_L, SIGMA_U, W, ZPZ, DX,
     3        X, S_L, S_U, V_L, V_U, N, NLB, ILB, NUB, IUB, IVAR, LAM,
     4        XORIG)
#else
         call ITER_OUT(ITER, MU, ERR, c_err, ERR_DINF, ERR_PINF,
     1        ERR_CMPL_MU, fnomu, M, C, cnrm, YPY, NIND, PZ, c_skip,
     1        ALPHA_DUAL, alpha_cut, ALPHA, c_alpha, c_watch, ls_count,
     2        nu_out, ieigs, condc, REGU,
     2        RV, RG, RGB, B, SIGMA_L, SIGMA_U, W, ZPZ, DX,
     3        X, S_L, S_U, V_L, V_U, N, NLB, ILB, NUB, IUB, IVAR, LAM,
     4        XORIG)
#endif
      endif
C
C     Check if done with the overall problem
C
      if( QERROR.ne.3 .and.
     1     (ERR_0.le.QTOL .and. (ERR_DINF_US.le.QDINFMAXTOL
     2     .and. ERR_CMPL_US.le.QDINFMAXTOL .and.
     3     ERR_PINF_US.le.QCMAXTOL)) ) then
         IERR = -1
         goto 9999
      elseif( ERR_0.le.QTOL .and. small_x_stop ) then
         IERR = 9
         goto 9999
      endif
C
C     A sophisticated user might want to decided on some other criterion
C     of her/his own to decrease the barrier parameter earlier.  In that
C     case you have to define the preprocessor flag
C     '-DUSE_IP_MYBARRIERTEST' at compilation time and return .true. as
C     result of the function below.  If you want to IPOPT to stop with
C     IERR = 0, return IERR = -3 here.
C
      decmu = .false.
#ifdef USE_IP_MYBARRIERTEST
      decmu = IP_MYBARRIERTEST(NORIG, XORIG, N, X, M, C, F, ITER,
     1     iter_start, MU, ALPHA, RESTO, IERR)
      if( IERR.ne.0 ) then
         if( IERR.ne.-3 ) then
            write(line,*)
     1           'solve_barrier: IP_MYBARRIERTEST returns IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
         goto 9999
      endif
#endif
C
C     Check if there is a forceful decrease of MU
C
      inquire(file='DECMU', exist=ex)
      if( RESTO.eq.0 .and. ex ) then
         call C_OUT(2,0,1,'File ''DECMU'' detected.  Decreasing MU.')
         open(80,file='DECMU',status='old')
         close(80,status='delete')
         decmu = .true.
      endif
C
C     Check the error; if small enough leave this loop
C        (barrier problem solve well enough)
C
CTODO allow to leave solve_barrier problem even in restoration phase
C     but this assume that all pre-restoration phase values are restored
C     (also LAM, V_L, V_U etc...)!!!
      if( (ERR.lt.ERR_BAR .or. small_x_stop .or. decmu)
     1     .and. (RESTO.eq.0 .or. RESTO.eq.-1.and.small_x_stop )) then
C      if( (ERR.lt.ERR_BAR .or. small_x_stop .or. decmu) .and.
C     1     (RESTO.eq.0.or.NIND.ne.0) ) then
         if( small_x_stop ) then
            IERR = -2
         endif
C
C     Don't decrease MU if it just has been decreased
C
         if( (ITER.eq.0.and.QERROR.ne.3) .or. ITER.gt.iter_start
     1        .or. small_x_stop
     1        .or. decmu) then
            MU = MUORIG
            RESTO = 0
            goto 9000
         endif
      endif

C
      if( QFULL.ne.1 ) then

         if( QCG.eq.0 .and.(QQUASI.ne.0 .or. RESTO.ne.0) ) then
C
C     Compute reduced primal-dual Hessian for barrier part
C        (store it in W)
C
            call GET_BB(N, NIND, M, ITER, IVAR, NFIX, IFIX, NORIG,
     1           XORIG, CSCALE, NLB, ILB, NUB, IUB, SIGMA_L, SIGMA_U,
     1           S_L, S_U, W, RESTO, KCONSTR, LRS, RS, LIS, IS,
     2           LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.gt.0 ) then
               write(line,*)
     1              'solve_barrier: get_bb returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            elseif( IERR.lt.0 ) then
               write(line,*)
     1              'solve_barrier: get_bb returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               IERR = 0
            endif
C
         else
C
C     Compute exact reduced Hessian
C
            if( QCG.eq.0 ) then
               call GET_EXACTW(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1              NORIG, XORIG, X, CSCALE, LAM, NLB, ILB, NUB, IUB,
     2              S_L, S_U, SIGMA_L, SIGMA_U, W, RESTO,
     3              KCONSTR, LRS, RS, LIS, IS,
     2              LRW-p_rwend, RW(p_rwend+1),
     4              LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C,
     6              EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
               if( IERR.gt.0 ) then
                  write(line,*)
     1                 'solve_barrier: get_exactw returns IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  goto 9999
               elseif( IERR.lt.0 ) then
                  write(line,*)
     1                 'solve_barrier: get_exactw returns IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  IERR = 0
               endif
            endif
C
         endif
C
      endif
C
C     Calculate or initialize B (Quasi-Newton matrix)
C     and compute overall reduced Hessian (store it in W)
C     or update limited memory representation
C
c      if( QQUASI.ne.0 .and..not.RESTO ) then
      c_skip = ' '
      if( QQUASI.ne.0 ) then
         if( QFULL.eq.0 .and. QCG.ne.0 ) then
            if( (abs(QCG).ne.1 .and. abs(QCG).ne.2) .or.
     1           (abs(QQUASI).ne.1 .and. abs(QQUASI).ne.5) ) then
               IERR = 4
               goto 9999
            endif
         endif
         if( abs(QQUASI).lt.6 .and. RESTO.eq.0 ) then
            if( QLINDAMP.ne.0.d0 ) then
               call C_OUT(2,0,1,'solve_barrier: This quasi-Newton option
     1 not implemented for nonzero DLINDAMP.')
               IERR = 4
               goto 9999
            endif
            call UPDATE_B(N, NIND, M, X, ITER, ERR, PZ, IVAR, NFIX,
     1           IFIX, NORIG, XORIG, CSCALE, RG, RGOLD, PYNRM_OLD,
     2           dummy, ALPHA, G, GOLD, LAM, SKIP_UPDATE, B, W, c_skip,
     1           NLB, ILB, NUB, IUB, S_L, S_U, SIGMA_L, SIGMA_U,
     3           KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     6           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.gt.0 ) then
               write(line,*)
     1              'solve_barrier: update_b returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            elseif( IERR.lt.0 ) then
               write(line,*)
     1              'solve_barrier: update_b returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               IERR = 0
            endif
         elseif( RESTO.ne.0 .or. QRESTO.eq.1 ) then
            call UPDATE_B_LM(N, NIND, M, X, ITER, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, NLB, ILB, S_L, NUB, IUB, S_U,
     2           MUORIG, G, LAM, SKIP_UPDATE, B, C_SKIP, RESTO,
     3           KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.ne.0 ) then
               write(line,*)
     1              'solve_barrier: update_b_lm returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
         endif
      endif
C
C     Compute primal and dual search direction
C
CTODO Do something about multipliers for L-BFGS
c      if ( RESTO.eq.0 ) INIT_LAM = .true.
 1000 continue  ! jump here returning for SOC or iterative refinement in RESTO
C
C     First make sure, that YPY is update to date for SOC
C
      if( SOC_FLAG.eq.1 .and. QFULL.eq.0 ) then
         call GET_YPY(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, NLB, ILB, NUB, IUB,
     2        S_L, S_U, C, YPY, NEWBAS, CONDC,
     1        KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     6        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1       'solve_barrier: Error: get_ypy(soc) ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1     'solve_barrier: Warning: get_ypy(soc) ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
      endif

      if( QFULL.eq.0 ) then

         newbas_store = NEWBAS
         call GET_STEP_RED(N, NIND, M, X, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, LAM, NLB, ILB, NUB, IUB,
     1        C, S_L, S_U, V_L, V_U,
     1        BNDS_L, BNDS_U, SIGMA_L, SIGMA_U, YPY, WCORR, RG, MU, RGB,
     1        B, W, PZ, ieigs, ZPZ, DX, DV_L, DV_U, ALPHA,
     1        alpha_dual, alpha_cut, c_alpha, soc_flag,
     1        NEWBAS, CONDC, RESTO, ERR, ERR_BAR,
     1        KCONSTR, LRS, RS, LIS, IS,
     1        LRW-p_rwend, RW(p_rwend+1),
     1        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     1        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'solve_barrier: get_step_red returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.lt.0 ) then
            write(line,*)
     1           'solve_barrier: get_step_red returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 0
         endif
         if( .not.newbas_store .and. NEWBAS ) then
            force_newbas = .true.
         endif

      elseif( QFULL.eq.1 ) then

         if( (QSOC.eq.0 .or. soc_flag.eq.0).and.RESTO.lt.2 ) then
            p_rwfull = p_rwend
            p_iwfull = p_iwend
         endif
         call GET_STEP_FULL(ITER, NORIG, N, NIND, M, idummy, X,
     1        XORIG, CSCALE, B, NLB, ILB, NUB, IUB, NFIX, IFIX, IVAR,
     1        G, C, S_L, S_U, V_L, V_U, BNDS_L, BNDS_U,
     1        SIGMA_L, SIGMA_U, MU, DX, DV_L, DV_U, LAM,
     1        LAMOLD, PZ, ALPHA,
     1        alpha_dual, alpha_cut, c_alpha, NEWBAS, ieigs, YPY, REGU,
     1        WCORR, RHS_RESTO, SIGMA_PN, eta_resto, soc_flag, wflag,
     1        ERR_BAR, ERR_PINF, RESTO, prwkeep, piwkeep, INIT_LAM,
     1        KCONSTR, LRS, LRS_END, RS, LIS, LIS_END, IS,
     1        LRW-p_rwfull, RW(p_rwfull+1),
     1        LIW-p_iwfull, IW(p_iwfull+1), IERR, EV_F, EV_C, EV_G,
     1        EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'solve_barrier: get_step_full returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.lt.0 ) then
            write(line,*)
     1           'solve_barrier: get_step_full returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 0
         endif

         if( (QSOC.ne.0 .and. soc_flag.eq.0 .and. QMERIT.ne.3
     1        .and. QMERIT.ne.0) .or.
     2        (abs(QMERIT).ge.4 .and. QRESTO.eq.1 .and. RESTO.gt.0) )
     3        then
C
C     readjust p_rwend and p_iwend if we need to keep entries in RW and IW
C     for later call of get_step_full for SOC
C
            p_rwend = p_rwfull + prwkeep
            p_iwend = p_iwfull + piwkeep
         endif
      else
         call C_OUT(2,0,1,'solve_barrier: invalid choice of QFULL!')
         stop
      endif
C
C     Check if we should stop
C
      IERR = CHECK_STOP()
      if( IERR.ne.0 ) goto 9999
C
C     Do the line search and obtain new point (and F and C at new point)
C
CTODO decide which one!
C      if( abs(QMERIT).eq.1 .or. abs(QMERIT).eq.2 ) then
      if( abs(QMERIT).eq.1 .or. abs(QMERIT).eq.2 .or. QMERIT.eq.0 ) then
         call LINESEARCH(ITER, N, NIND, M, X, IVAR, NLB, ILB, NUB,
     1        IUB, BNDS_L, BNDS_U, DX, DV_L, DV_U, S_L, S_U,
     2        V_L, V_U, SIGMA_L, SIGMA_U, NORIG, XORIG, CSCALE, MU,
     3        ERR, YPY, LAM, LAMOLD, regu,
     2        PZ, ZPZ, G, WCORR, F, C, cnrm0, NEWBAS,
     3        ALPHA, alpha_dual, ls_count, c_watch, nu_out,
     4        SKIP_UPDATE, soc_flag, KCONSTR,
     5        LRS, LRS_END, RS, LIS, LIS_END, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'solve_barrier: linesearch returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.eq.-2 ) then
            small_x_stop = .true.
            IERR = 0
         elseif( IERR.lt.0 ) then
            write(line,*)
     1           'solve_barrier: linesearch returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 0
         endif
C
C     Check if SOC step computation is requested
C
         if( soc_flag.eq.1 ) goto 1000
C
C     Can forget KKT factorization
C
         if( QSOC.ne.0 .and. QFULL.eq.1 ) then
            p_rwend = p_rwfull
            p_iwend = p_iwfull
         endif

      elseif( QMERIT.eq.3 ) then
         c_watch = ' '
         call AUG_LAG(ITER, N, NIND, M, X, IVAR, NLB, ILB, NUB,
     1        IUB, BNDS_L, BNDS_U, DX, DV_L, DV_U, S_L, S_U,
     2        V_L, V_U, NORIG, XORIG, CSCALE, MU, LAM, LAMOLD, G, F, C,
     3        cnrm0, ALPHA, alpha_dual, ls_count, nu_out, KCONSTR,
     5        LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     6        EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*) 'solve_barrier: aug_lag returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.lt.0 ) then
            write(line,*) 'solve_barrier: aug_lag returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 0
         endif
      elseif( abs(QMERIT).ge.4 ) then
         call FILTER(ITER, N, NIND, M, X, IVAR, NFIX, IFIX, NLB, ILB,
     1        NUB, IUB, BNDS_L, BNDS_U, DX, DV_L, DV_U, S_L,
     2        S_U, V_L, V_U, SIGMA_L, SIGMA_U, PZ, NORIG, XORIG, CSCALE,
     2        MU, MUORIG, LAM, LAMOLD, G, F, C, cnrm0, ALPHA,
     3        alpha_dual, RESTO, ls_count, nu_out, c_alpha, c_watch,
     5        soc_flag, wflag, INIT_LAM, ERR, RHS_RESTO, SIGMA_PN,
     1        eta_resto, B, WCORR, REGU, KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5        EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
CTODO maybe jump back to step computation if just switched to restoration phase
         if( IERR.gt.0 ) then
            write(line,*) 'solve_barrier: filter returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.eq.-2 ) then
            small_x_stop = .true.
            IERR = 0
         elseif( IERR.lt.0 ) then
            write(line,*) 'solve_barrier: filter returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 0
         endif
C
C     Check if SOC step computation is requested
C
         if( soc_flag.eq.1 .or. RESTO.ge.2 ) goto 1000
C
C     Can forget KKT factorization
C
         if( QFULL.eq.1 ) then
            p_rwend = p_rwfull
            p_iwend = p_iwfull
         endif

      else
         call C_OUT(2,0,1,'solve_barrier: Invalid QMERIT!')
         IERR = 4
         goto 9999
      endif
C
C     Check if iterates become too large
C
      i = IDAMAX(N, X, 1)
      if( abs(X(i)).gt.CHECK_XMAX ) then
         write(line,*) 'solve_barrier: ||x||_max = ', abs(X(i))
         call C_OUT(2,1,1,line)
         IERR = 13
         goto 9999
      endif
C
C     Make sure that the multipliers do not differ too much from MU/S
C
      call CORRECT_V(NLB, S_L, V_L, NUB, S_U, V_U, MU)
C
C     Store old RG
C
      if( (QQUASI.ne.0.and.abs(QQUASI).lt.6) .and.
     1     .not.small_x_stop ) then
         if( QQUASI.gt.0 ) then
            call DCOPY(NIND, RG, 1, RGOLD, 1)
         else
            call DCOPY(N, G, 1, GOLD, 1)
            if( M.gt.0 ) then
               call CONSTR(11, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1              NORIG, XORIG, CSCALE, dummy, dummy, idummy, idummy,
     2              KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4              IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5              LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     6              EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
               if( IERR.gt.0 ) then
                  write(line,*)
     1                 'solve_barrier: constr-11 returns IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  goto 9999
               elseif( IERR.lt.0 ) then
                  write(line,*)
     1                 'solve_barrier: constr-11 returns IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  IERR = 0
               endif
            endif
         endif
      elseif( abs(QQUASI).ge.6 ) then
         if( M.gt.0 ) then
            call CONSTR(11, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, dummy, dummy, idummy, idummy,
     2           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.gt.0 ) then
               write(line,*)
     1              'solve_barrier: constr-11 returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            elseif( IERR.lt.0 ) then
               write(line,*)
     1              'solve_barrier: constr-11 returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               IERR = 0
            endif
         endif
      endif
C
C     Store old pY norm
C
      if( M.gt.0 ) THEN
CTODO norm of WHOLE YPY (incl indep vars)?
         PYNRM_OLD = DNRM2(M, YPY, 1)
      else
         PYNRM_OLD = 0.d0
      endif
C
C     Increase iteration counter
C
      ITER = ITER + 1
C
C     Check if we should stop
C
      IERR = CHECK_STOP()
      if( IERR.ne.0 ) goto 9999
C
C     Check maximal number of iterations
C
      if( ITER.eq.QMAXITER ) then
         write(line,7000) ITER
 7000    format('Maximal number of iterations (',i6,') exceeded!')
         call C_OUT(2,0,1,line)
C The following lines might be more appropriate for interactive usage.
C         writee(*,*) 'Please enter new MAXITER (0: quit)'
C         read(*,*) QMAXITER
         if( QMAXITER.le.ITER ) then
C            if( RESTO.ne.0 ) then
C               IERR = 12
C            else
               IERR = 1
C            endif
            goto 9999
         endif
      endif
C
C     Compute range space step PY at new point
C
      if( QFULL.eq.0 .or. (QQUASI.ne.0.and.abs(QQUASI).lt.6)
     1     .or. QLAMBDA.eq.0 ) then
         call GET_YPY(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, NLB, ILB, NUB, IUB,
     2        S_L, S_U, C, YPY, NEWBAS, CONDC,
     1        KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5        EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'solve_barrier: Error: get_ypy ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1          'solve_barrier: Warning: get_ypy ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
C
C     Check if we should stop
C
         IERR = CHECK_STOP()
         if( IERR.ne.0 ) goto 9999
C
C     if necessary get new partition and reorder everything
C
         if( NEWBAS .or. force_newbas ) then
C
C     Get new partition and reoder everything
C
            call C_OUT(2,0,1,'Changing Basis!')
CTODO Check if OK...
            if( RESTO.ne.0 ) then
               write(*,*) 'Don''t know yet for RESTO.'
               stop
            endif
            call NEW_BASIS(ITER, N, M, NIND, NORIG, X, XORIG, CSCALE,
     1           NFIX, IFIX, IVAR, NLB, ILB, NUB, IUB,
     2           S_L, S_U, C, LAM, LAMOLD, YPY, ALPHA,
     3           DX, PZ, G, RGOLD, B, SKIP_UPDATE,
     3           KCONSTR, LRS, RS, LIS, IS,
     2           LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.gt.0 ) then
               write(line,*)
     1              'solve_barrier: Error: new_basis return IERR = ',
     2              IERR
               call C_OUT(2,0,1,line)
               goto 9999
            elseif( IERR.ne.0 ) then
               write(line,*)
     1              'solve_barrier: Warning: new_basis returns IERR = ',
     2              IERR
               call C_OUT(2,0,1,line)
            endif
         endif
      endif
C
C     Get G at new point
C
      if( RESTO.le.0 ) then
         MUORIG = MU
      endif
      call GET_G(N, X, IVAR, NORIG, XORIG, M, CSCALE, NLB, ILB,
     1     S_L, NUB, IUB, S_U, MUORIG, G, LRW-p_rwend, RW(p_rwend+1),
     2     LIW-p_iwend, IW(p_iwend+1), IERR, EV_G, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1        'solve_barrier: Error: get_g ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.ne.0 ) then
         write(line,*)
     1        'solve_barrier: Warning: get_g ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
      endif
C
C     Compute reduced gradient of objective function RG
C
CTODO do we update during restoration?
      if( QFULL.eq.0 .or. (QQUASI.ne.0.and.abs(QQUASI).lt.6)
     1     .or. abs(QERROR).eq.1 ) then
         call GET_RG(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, G, RG,
     1        KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'solve_barrier: Error: get_rg ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'solve_barrier: Warning: get_rg ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
      endif
C
C     Compute "reduced" dual variables RV
C
      if( QFULL.eq.0 .or. abs(QERROR).eq.1 ) then
         call GET_RV(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE,
     1        NLB, ILB, NUB, IUB, V_L, V_U, RV,
     1        KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'solve_barrier: Error: get_rv ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'solve_barrier: Warning: get_rv ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
      endif
C
C     For QLAMBDA = 1: get coordinate multipliers
C
      if( QLAMBDA.eq.1 ) then
         call GET_LAMBDA(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, G, NLB, ILB, NUB, IUB, V_L, V_U,
     1        LAM, KCONSTR, LRS, RS, LIS, IS,
     1        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     1        IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     2        EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'solve_barrier: get_lambda returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.lt.0 ) then
            write(line,*)
     1           'solve_barrier: get_lambda returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 0
         endif
      endif
C
C     Now we have everything to start a new round in our loop
C
      goto 100

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                       C
C     END OF WHILE LOOP                                                 C
C                                                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C     Prepare exit
C
 9000 continue

 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine SOLVE_BARRIER_WS(N, X, M, NLB, NUB, NZA, LRW, LIW,
     1     FILLINFACT, EV_H, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      external EV_H
      double precision DAT(*)
      integer IDAT(*)
      double precision FILLINFACT, X(N)
      integer liw1, lrw1, liw2, lrw2, prwkeep, piwkeep
      character*80 line

      call OPTERROR_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      if( QFULL.ne.1 ) then
         if( QCG.eq.0 .and.(QQUASI.ne.0 .or. QRESTO.ne.2) ) then
            call GET_BB_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
            LRW = max(LRW, lrw1)
            LIW = max(LIW, liw1)
         else
            if( QCG.eq.0 ) then
               call GET_EXACTW_WS(N, M, NLB, NUB, NZA, lrw1, liw1,
     1              DAT, IDAT)
               LRW = max(LRW, lrw1)
               LIW = max(LIW, liw1)
            endif
         endif
      endif

      if( QQUASI.ne.0 ) then
         if( abs(QQUASI).lt.6 ) then
            call UPDATE_B_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
            LRW = max(LRW, lrw1)
            LIW = max(LIW, liw1)
         else
            call UPDATE_B_LM_WS(N, M, NLB, NUB, NZA, lrw1, liw1,
     1           DAT, IDAT)
            LRW = max(LRW, lrw1)
            LIW = max(LIW, liw1)
         endif
      endif

      lrw2 = 0
      liw2 = 0
      prwkeep = 0
      piwkeep = 0
      if( QSOC.ne.0 .and. QFULL.eq.0 ) then
         call GET_YPY_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         lrw2 = max(lrw2, lrw1)
         liw2 = max(liw2, liw1)
      endif
      if( QFULL.eq.0 ) then
         call GET_STEP_RED_WS(N, M, NLB, NUB, NZA, lrw1, liw1,
     1        DAT, IDAT)
         lrw2 = max(lrw2, lrw1)
         liw2 = max(liw2, liw1)
      else
         call GET_STEP_FULL_WS(N, X, M, NLB, NUB, NZA, lrw1, liw1,
     1        FILLINFACT, prwkeep, piwkeep, EV_H, DAT, IDAT)
         lrw2 = max(lrw2, lrw1)
         liw2 = max(liw2, liw1)
         if( (QSOC.eq.0 .or. QMERIT.eq.3 .or. QMERIT.eq.0) .and.
     1        (abs(QMERIT).lt.4) .or. QRESTO.ne.1 ) then
            prwkeep = 0
            piwkeep = 0
         endif
      endif
      if( abs(QMERIT).eq.1 .or. abs(QMERIT).eq.2 .or. QMERIT.eq.0 ) then
         call LINESEARCH_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
      elseif( QMERIT.eq.3 ) then
         call AUG_LAG_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
      else
         call FILTER_WS(N, X, M, NLB, NUB, NZA, lrw1, liw1, EV_H,
     1        DAT, IDAT)
      endif
      lrw2 = max(lrw2, prwkeep+lrw1)
      liw2 = max(liw2, piwkeep+liw1)
      LRW = max(LRW, lrw2)
      LIW = max(LIW, liw2)
      if( QQUASI.ne.0 .and. M.gt.0 ) then
         call CONSTR_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)
      endif

      if( QFULL.eq.0 .or. (QQUASI.ne.0.and.abs(QQUASI).lt.6)
     1     .or. QLAMBDA.eq.0 ) then
         call GET_YPY_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)
         call NEW_BASIS_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)
      endif

      call GET_G_WS(N, M, NLB, NUB, NZA, lrw1, liw1)
      LRW = max(LRW, lrw1)
      LIW = max(LIW, liw1)

      if( QFULL.eq.0 .or. (QQUASI.ne.0.and.abs(QQUASI).lt.6)
     1     .or. abs(QERROR).eq.1 ) then
         call GET_RG_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)
      endif

      if( QFULL.eq.0 .or. abs(QERROR).eq.1 ) then
         call GET_RV_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)
      endif

      if( QLAMBDA.eq.1 ) then
         call GET_LAMBDA_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)
      endif

      if( QPRINT.ge.4 ) then
         write(line,1000)'solve_barrier_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
