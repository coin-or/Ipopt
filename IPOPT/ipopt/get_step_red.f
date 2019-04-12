C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_STEP_RED(N, NIND, M, X, ITER, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, LAM, NLB, ILB, NUB, IUB,
     1     C, S_L, S_U, V_L, V_U,
     1     BNDS_L, BNDS_U, SIGMA_L, SIGMA_U, YPY, WCORR, RG, MU, RGB,
     1     B, W, PZ, IEIGS, ZPZ, DX, DV_L, DV_U, ALPHA,
     1     ALPHA_DUAL, ALPHA_CUT, C_ALPHA, SOC_FLAG,
     1     NEWBAS, CONDC, RESTO, ERR, ERR_BAR,
     1     KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR, EV_F,
     1     EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: get_step_red.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute primal and dual steps from reduced system
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
CP   N         I    INT    number of (free) variables; first NIND variables
CP                         are independent; remaining dependent
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of dependent variables
CP   X         I    DP     actual primal iterate
CP   ITER      I    INT    iteration counter
CP   IVAR      I    INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP                            X(i) corresponds to XORIG(IVAR(i))
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP                            (assumed to be in increasing order)
CP   NORIG     I    INT    total number of variables (incl. fixed vars)
CP   XORIG     I    INT    actual iterate
CP                            (original order as in problem statement)
CP   CSCALE    I    DP     scaling factors for constraints
CP   LAM       I    DP     Lagrangian multipliers
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. S_L(i) is slack for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. S_U(i) is slack for X(IUB(i)) )
CP   C         I    DP     values of constraints at X
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
CP   YPY       I    DP     range space step (ordered like X)
CP   WCORR     O    DP     correction term for null space step
CP   RG        I    DP     reduced gradient of objective function
CP   MU        I    DP     barrier parameter
CP   RGB       O    DP     reduced gradient of barrier function
CP   B         I    DP     Quasi-Newton estimate of reduced Hessian of
CP                            original NLP
CP                         for CG: contains preconditioner
CP   W         I    DP     reduced Hessian of barrier problem
CP   PZ        O    DP     null space step (indepentent variables)
CP   IEIGS     O    INT    number of negative eigenvalues in overall
CP                            reduced Hessian
CP   ZPZ       O    DP     null space step (dependent variables)
CP   DX        O    DP     step for X (primal)
CP   DV_L      O    DP     step for V_L (dual variables for lower bounds)
CP   DV_U      O    DP     step for V_U (dual variables for upper bounds)
CP   ALPHA     O    DP     maximal steps size compatible with fraction
CP                            to boudary rule (primal vars)
CP   ALPHA_DUAL O   DP     maximal steps size compatible with fraction
CP                            to boudary rule (dual vars)
CP   ALPHA_CUT O    DP     first alpha on output
CP   C_ALPHA   O    C*1    for output on ALPHA
CP   SOC_FLAG  I    INT    =0: not in SOC computation
CP                          1: in SOC computation (don't need to refactorize)
CP   NEWBAS    O    LOG    =.true.: Basis has become pretty bad; get new one!
CP                            (this is determined during heuristic!)
CP   CONDC     O    DP     estimated condition number of basis matrix C
CP   RESTO     I    INT    <>0: we are in restoration phase
CP   ERR       I    DP     current KKT error for barrier problem
CP   ERR_BAR   I    DP     KKT error tolerance for barrier problem
CP   KCONSTR   I    INT    KCONSTR(1): LRS for CONSTR
CP                         KCONSTR(2): P_LRS for CONSTR
CP                         KCONSTR(3): LIS for CONSTR
CP                         KCONSTR(4): P_LIS for CONSTR
CP                         KCONSTR(5): LRW for CONSTR
CP                         KCONSTR(6): LIW for CONSTR
CP   LRS       I    INT    total length of RS
CP   RS        S    DP     DP storage space (all!)
CP   LIS       I    INT    total length of IS
CP   IS        S    INT    INT storage space (all!)
CP   LRW       I    INT    length of RW
CP   RW        W    DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW       I    INT    length of IW
CP   IW        W    INT    can be used as INT work space but content will be
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
CS    DAXPY
CS    DASUM
CS    DSCAL
CS    GET_WCORR
CS    GET_RGB
CS    GET_PZ
CS    GET_PZ_CG
CS    GET_ZPZ
CS    GET_D
CS    CUTALPHA
CS    GET_YPY
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
      integer NIND
      integer M
      double precision X(N)
      integer ITER
      integer IVAR(N)
      integer NFIX
      integer IFIX(NFIX)
      integer NORIG
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      double precision LAM(M)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision C(M)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision BNDS_L(NLB)
      double precision BNDS_U(NUB)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision YPY(N)
      double precision WCORR(NIND)
      double precision RG(NIND)
      double precision MU
      double precision RGB(NIND)
      double precision B( (NIND*(NIND+1))/2 )
      double precision W( (NIND*(NIND+1))/2 )
      double precision PZ(NIND)
      integer IEIGS
      double precision ZPZ(M)
      double precision DX(N)
      double precision DV_L(NLB)
      double precision DV_U(NUB)
      double precision ALPHA
      double precision ALPHA_DUAL
      double precision ALPHA_CUT
      character*1 C_ALPHA
      integer SOC_FLAG
      logical NEWBAS
      double precision CONDC
      double precision ERR
      double precision ERR_BAR
      integer RESTO
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
      integer p_rwend, p_iwend

      double precision cnrm, tau
      double precision DASUM
      character*80 line

      double precision REGU
      save             REGU

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
C
C     If we are in computation of SOC step, for now skip computation of
C     new PZ and jump directly to computation of overall D
C
      if( SOC_FLAG.eq.1 ) goto 1000
C
C     Jump back to here, if range space step has been recomputed
C     for anti-crash heuristic
C
 120  continue
C
C     Compute correction term for null space step WCORR
C     (if in second order correction step computation, no need to recompute
C      wcorr if we use pcg, since here SOC does not change PZ for simplicity,
C      only the YPY part is changed.)
C
      if( SOC_FLAG.ne.1 .or. QCG.eq.0 ) then
         call GET_WCORR(N, NIND, M, X, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, LAM, NLB, ILB, NUB, IUB,
     2        S_L, S_U, BNDS_L, BNDS_U,
     1        SIGMA_L, SIGMA_U, YPY, 0d0, WCORR, RESTO,
     1        KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.lt.0 ) then
            write(line,*)
     1           'get_step_red: Warning in get_wcorr, IERR = ',IERR
            call C_OUT(2,0,1,line)
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'get_step_red: Error in get_wcorr, IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      endif
C
C     Get reduced gradient of barrier function RGB
C
      if( RESTO.eq.0 .and. .not.SOC_FLAG.eq.1 ) then
         call GET_RGB(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, RG, NLB, ILB, NUB, IUB, S_L, S_U,
     1        MU, RGB, KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     3        IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.lt.0 ) then
            write(line,*)
     1           'get_step_red: Warning in get_rgb, IERR = ',IERR
            call C_OUT(2,0,1,line)
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'get_step_red: Error in get_rgb, IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      endif
C
C     Compute null step PZ
C
      if( SOC_FLAG.ne.1 ) IEIGS = 0
      if( QCG.ne.0 .and. RESTO.eq.0 .and. .not.SOC_FLAG.eq.1 ) then
         if( QDAMP.ne.0 ) then
            call C_OUT(2,0,1,'get_step_red: Can''t do damping for cg.')
            IERR = 4
            goto 9999
         endif
         call GET_PZ_CG(ITER, MU, NIND, RGB, WCORR, B, PZ,
     1        IEIGS, N, M, X, IVAR, NFIX, IFIX, NORIG, XORIG, CSCALE,
     2        LAM, NLB, ILB, NUB, IUB, S_L, S_U, SIGMA_L,
     3        SIGMA_U, ERR, ERR_BAR, NEWBAS, KCONSTR, LRS, RS, LIS, IS,
     4        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*) 'get_step_red: get_pz_cg returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.lt.0 ) then
            write(line,*) 'get_step_red: get_pz_cg returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 0
         endif
      endif

      if( QCG.eq.0 .or. RESTO.ne.0 ) then
         if( RESTO.ne.0 ) then
            call C_OUT(2,0,1,'get_step_red: RESTO not yet implemented.')
            IERR = 4
            goto 9999
         endif
         if( QCG.ne.0 .and.RESTO.eq.0 ) then
c            WRITE(*,*)
c     1           'Hey, need to call GET_EXACTW here in get_step_red!'
            IERR = 4582
            goto 9999
         endif
         cnrm = DASUM(M, C, 1)
         call GET_PZ(ITER, MU, NIND, W, RGB, WCORR, cnrm, B, PZ, IEIGS,
     2        N, M, X, IVAR, NFIX, IFIX, NORIG, XORIG, CSCALE, LAM, NLB,
     2        ILB, NUB, IUB, S_L, S_U, BNDS_L, BNDS_U, SIGMA_L,
     3        SIGMA_U, YPY, REGU, SOC_FLAG, RESTO,
     1        KCONSTR, LRS, RS, LIS, IS,
     4        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*) 'get_step_red: get_pz returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.lt.0 ) then
            write(line,*) 'get_step_red: get_pz returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 0
         endif
      endif
C
C     Compute dependent part of null space step ZPZ
C
      call GET_ZPZ(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, PZ, ZPZ,
     1     KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*) 'get_step_red: get_zpz returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.lt.0 ) then
         write(line,*) 'get_step_red: get_zpz returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         IERR = 0
      endif
C
C     Compute primal and dual search directions
C
 1000 continue
      call GET_D(N, NIND, M, X, YPY, PZ, ZPZ, NLB, ILB, NUB, IUB,
     1           S_L, S_U, V_L, V_U, BNDS_L, BNDS_U,
     1           SIGMA_L, SIGMA_U, MU, DX, DV_L, DV_U, RESTO)
C
      if( RESTO.ne.0.and.QRESTO.eq.3 ) goto 9999
C
C     Compute maximal step size alpha within bounds
C
      tau = QTAU
      call CUTALPHA(N, X, DX, DV_L, DV_U, NLB, ILB, NUB, IUB,
     1              BNDS_L, BNDS_U, S_L, S_U, V_L, V_U, tau,
     2              ALPHA, ALPHA_DUAL, C_ALPHA)
      ALPHA_CUT = ALPHA   ! Show this alpha in normal case
C
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine GET_STEP_RED_WS(N, M, NLB, NUB, NZA, LRW, LIW,
     1     DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1

      call GET_WCORR_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      call GET_RGB_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
      LRW = max(LRW,lrw1)
      LIW = max(LIW,liw1)

      if( QCG.ne.0 ) then
         call GET_PZ_CG_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW,lrw1)
         LIW = max(LIW,liw1)
      endif

      if( QCG.eq.0 .or. (abs(QMERIT).ge.4 .and. abs(QRESTO).ne.2) ) then
         call GET_PZ_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW,lrw1)
         LIW = max(LIW,liw1)
      endif

      call GET_ZPZ_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
      LRW = max(LRW,lrw1)
      LIW = max(LIW,liw1)

      return
      end
