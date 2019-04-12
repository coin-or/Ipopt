C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_STEP_FULL(ITER, NORIG, N, NIND, M, NZORIG, X,
     1     XORIG, CSCALE, B, NLB, ILB, NUB, IUB, NFIX, IFIX, IVAR, G, C,
     1     S_L, S_U, V_L, V_U, BNDS_L, BNDS_U,
     1     SIGMA_L, SIGMA_U, MU, DX, DV_L, DV_U, LAM, LAMOLD, PZ, ALPHA,
     1     ALPHA_DUAL, ALPHA_CUT, C_ALPHA, NEWBAS, IEIGS, YPY, REGU,
     1     DHD, RHS_RESTO, SIGMA_PN, ETA_RESTO, SOC_FLAG, WFLAG,
     1     ERR_BAR, ERR_CNSTR, RESTO, PRWKEEP, PIWKEEP, INIT_LAM,
     1     KCONSTR, LRS, LRS_END, RS, LIS, LIS_END, IS,
     1     LRW, RW, LIW, IW, IERR, EV_F, EV_C, EV_G, EV_A,
     1     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: get_step_full.F 675 2005-07-26 18:47:13Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute primal and dual steps from FULL system
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
CP                            (if -1: initialize pointers for storage space
CP                                    and return number of nonzeros for Hessian
CP                                    as NFIX)
CP   NORIG     I    INT    number of all variables including fixed vars
CP   N         I    INT    number of (free) variables; first M variables
CP                         are dependent; remaining independent
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of dependent variables
CP   NZORIG    I    INT    (only for ITER=-1): number of nonzeros in Jacobian
CP                            including columns for fixed vars
CP   X         I    DP     actual primal iterate
CP   XORIG     I    DP     actual iterate
CP   CSCALE    I    DP     scaling factors for constraints
CP   B         I    DP     Quasi-Newton estimate of reduced Hessian of
CP                            original NLP
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. S_L(i) is slack for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. S_U(i) is slack for X(IUB(i)) )
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP                            (assumed to be in increasing order)
CP   IVAR      I    INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP                            X(i) corresponds to XORIG(IVAR(i))
CP   G         I    DP     gradient of objective function
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
CP   MU        I    DP     barrier parameter (=0, if Error should be computed
CP                            for overall NLP)
CP   DX        O    DP     step for X (primal)
CP   DV_L      O    DP     step for V_L (dual variables for lower bounds)
CP   DV_U      O    DP     step for V_U (dual variables for upper bounds)
CP   LAM      I/O   DP     multipliers (only if QQUASI = 0,2)
CP   LAMOLD   I/O   DP     multipliers (only if QQUASI = 0,2)
CP   PZ        O    DP     null space step (indepentent variables)
CP                            (only if QQUASI.ne.0)
CP   ALPHA     O    DP     maximal steps size compatible with fraction
CP                            to boudary rule (primal vars)
CP   ALPHA_DUAL O   DP     maximal steps size compatible with fraction
CP                            to boudary rule (dual vars)
CP   ALPHA_CUT O    DP     first alpha on output
CP   C_ALPHA   O    C*1    for output on ALPHA
CP   NEWBAS    I    LOG    =.true. if basis changed (need new pivots!)
CP   IEIGS     O    I      number of trials for regularization
CP   YPY       I    DP     range space step
CP                           (only if penalty parameter is based on LAM, i.e.
CP                            QLAMBDA <> 0)
CP   REGU      O    DP     regularization factor (added regu*I to diagonal)
CP   DHD       O    DP     d^T Hessian d
CP   RHS_RESTO I    DP     QRESTO.eq.1: RHS if in IP restoration
CP                                      entries N+M+1..N+M+N 
CP                                        are for regularization term in Hessian
CP                                        (DO NOT CHANGE)
CP             O                        W * DX (+old DX for RESTO=2)
CP                                      last (N+M) entries contain sum of
CP                                        matrix column entries
CP   SIGMA_PN  I    DP     QRESTO.eq.1: Sigma for P and N variables
CP             W                        used as workspace
CP   ETA_RESTO I    DP     QRESTO.eq.1: factor for regularization term
CP   SOC_FLAG  I    INT    =0: not in SOC computation
CP                          1: in SOC computation (do not need to refactorize)
CP   WFLAG     I    INT    =0: inactive watchdog
CP                          1: computation of watchdog trial step
CP                               (store regu etc)
CP                          2: last watchdog step did not help (restore regu)
CP   ERR_BAR   I    DP     error tolerance for the barrier problem
CP   ERR_CNSTR I    DP     current infeasibility error
CP   RESTO     I    INT    =0: we are in regular IPOPT iteration
CP                         =1: we are in regular restoration phase iteration
CP                         =2: we are in iterative refinement for restoration
CP                             phase
CP   PRWKEEP   O    INT    If this routine is called lateron for SOC, the
CP                           calling routine must keep the first PRWKEEP
CP                           entries in RW
CP   PIWKEEP   O    INT    If this routine is called lateron for SOC, the
CP                           calling routine must keep the first PIWKEEP
CP                           entris in IW
CP   INIT_LAM I/O   LOG    Indicates whether equality multipliers are to be
CP                           initialized (might be set to. false. in line
CP                           search)
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
CS    DAXPY
CS    DSCAL
CS    DNRM2
CS    GET_H
CS    GET_YPY
CS    MA27_CALL
CS    MA47_CALL
CS    MA57_CALL
CS    WSSMP_CALL
CS    WGSMP_CALL
CS    MA48SYM_CALL
CS    CONSTR
CS    CUTALPHA
CS    MC19_CALL
CS    MC30AD
CS    C_OUT
CS    GSF_BACKSOLVE
CS    LMBACKSOLVE
CS    TIMER
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
      integer NORIG
      integer N
      integer NIND
      integer M
      integer NZORIG
      double precision X(N)
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      double precision B(*)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      integer NFIX
      integer IFIX(NFIX)
      integer IVAR(N)
      double precision G(N)
      double precision C(M)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision BNDS_L(NLB)
      double precision BNDS_U(NUB)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision MU
      double precision DX(N)
      double precision DV_L(NLB)
      double precision DV_U(NUB)
      double precision LAM(M)
      double precision LAMOLD(M)
      double precision PZ(NIND)
      double precision ALPHA
      double precision ALPHA_DUAL
      double precision ALPHA_CUT
      character*1 C_ALPHA
      logical NEWBAS
      integer IEIGS
      double precision YPY(N)
      double precision REGU
      double precision DHD
      double precision RHS_RESTO(3*N+2*M)
      double precision SIGMA_PN(N)
      double precision ETA_RESTO
      integer SOC_FLAG
      integer WFLAG
      double precision ERR_BAR
      double precision ERR_CNSTR
      integer RESTO
      integer PRWKEEP
      integer PIWKEEP
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
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer p_rwend, p_iwend, p_kkt, p_jcntmp, p_check, p_spalte
      integer p_rhs, p_w, p_hcopy, p_rhscopy, p_iwma, p_rhs_store
      integer p_scale, p_tmp, p_scale_old, p_a, p_realw, pb_uk, pb_vk
      integer inz, i, j, idummy, kr, kc, ifail, k
      integer neigs, nzero, nlm, nlmmax, refinement_iter
      double precision dummy, ele, vin(2), nrm_ratio_old
      double precision nrm_ratio, gzbw, regu_old, nrm_sol, nrm_mat
      double precision nrm_res, nrm_rhs, nrm_ratio_max, nrm_ratio_min
      double precision sigmak, times, timef, constrregu
      character*140 line(2)
      logical ldummy, l_lbfgsreset, is_dep_this_iter, naninf
      logical l_cs              ! for RESTO: if true, do CS and ignore Newton
      logical l_iter_ref, recomp_regu, skip_jac

      double precision REGU_STORE, DEPCONDIAGVAL, MACHEPS
      save             REGU_STORE, DEPCONDIAGVAL, MACHEPS
      double precision REGU_STORE_OLD, DEPCONDIAGVAL_OLD
      save             REGU_STORE_OLD, DEPCONDIAGVAL_OLD
      integer NZKKT, NZA, NNZH, LIWMA, LA, N_SWCORNER
      save    NZKKT, NZA, NNZH, LIWMA, LA, N_SWCORNER
      logical HAVE_PIV, L_INC_LA, L_INC_LIWMA, L_INCPIVTOL
      save    HAVE_PIV, L_INC_LA, L_INC_LIWMA, L_INCPIVTOL

      integer P_JCN, P_IRN
      save    P_JCN, P_IRN

      double precision DNRM2, DASUM, DDOT, D1MACH
      integer FFINITE, IDAMAX, CHECK_STOP

      integer CHECKDEPENCOUNTER, CHECKDEPENCOUNTERMAX
      save CHECKDEPENCOUNTER
      parameter( CHECKDEPENCOUNTERMAX = 3 )

C     Define macros for using MALLOC or not
#ifdef USE_MALLOC
      integer type
#ifdef SIZEOF_INT_P_IS_4
      integer*4 PTR_KKT, PTR_IWMA
#elif SIZEOF_INT_P_IS_8
      integer*8 PTR_KKT, PTR_IWMA
#else
      Size of integer pointers unknown
#endif
      integer IP_MALLOC
      common /IP_ALLOC/ PTR_KKT, PTR_IWMA
      save /IP_ALLOC/
#endif

      double precision MPEC_ETA
      common /MPEC/ MPEC_ETA
      save /MPEC/
C
C     The following variables are for a heuristic to handle structurally
C     singular problems (constraints or Hessian).  It can be switched off
C     by setting ITERDEP_MAX to 0.
C
      integer ITERDEP_MAX
      parameter( ITERDEP_MAX = 3 )
      integer ITERDEP, REFINEMENTHELPLESS
      logical POSSIBLE_DEPCON, IS_DEPCON, POSSIBLE_DEPHES, IS_DEPHES
      save ITERDEP, POSSIBLE_DEPCON, IS_DEPCON, REFINEMENTHELPLESS
      save POSSIBLE_DEPHES, IS_DEPHES

      integer CORRECTTYPE
      parameter( CORRECTTYPE = 0 )
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      call TIMER(times)
      if( PRINTTIMING .and. QPRINT.ge.2 ) then
         write(line,*) 'get_step_full start time:', times
         call C_OUT(2,0,1,line)
      endif
C
C     ITER = -1: initialize data structures
C
      if( ITER.eq.-1 ) then
         NZKKT = N + NZORIG     ! worst estimate
         if( QQUASI.ne.0 .and. abs(QQUASI).lt.6 ) then
            NZKKT = NZKKT + (NIND*(NIND+1))/2
         elseif( QQUASI.eq.0 ) then
C
C     call EVAL_H to obtain number of nonzeros in Hessian
C
            call GET_H(-1, N, NIND, idummy, dummy, idummy, NORIG,
     1           XORIG, dummy, M, NNZH, dummy, dummy, idummy, idummy,
     1           idummy, dummy, idummy, idummy, idummy, EV_H, DAT, IDAT)
            NZKKT = NZKKT + NNZH
            if( INMEMCHECK.eq.0 ) then
               write(line,100) NNZH
 100           format('Number of nonzeros in Hessian : ',i8)
               call C_OUT(2,0,1,line)
            endif
            NFIX = NNZH
         endif
         if( QDEPCONDIAG.gt.0.d0 .or. QSYMSOLV.eq.3 .or. QSYMSOLV.eq.4
     1        .or. (abs(QMERIT).ge.4.and.QRESTO.eq.1)) then
            NZKKT = NZKKT + M   ! space heuristic for dependent constraints
                                ! or for restoration phase step calculation
                                ! OR if we use WSMP
            N_SWCORNER = M
         else
            N_SWCORNER = 0
         endif
C
C     reserve space for solver
C
         if( QSYMSOLV.eq.0 ) then
            call MA47_CALL(0, N+M, NZKKT, dummy, idummy, idummy, dummy,
     1           idummy, ldummy, idummy, idummy, ldummy, idummy, idummy,
     1           LRS_END, dummy, LIS_END, idummy, idummy, dummy, idummy,
     1           idummy, IERR)
         elseif( QSYMSOLV.eq.1 ) then
            call MA27_CALL(0, N+M, NZKKT, dummy, idummy, idummy, dummy,
     1           idummy, ldummy, idummy, idummy, ldummy, idummy, idummy,
     1           LRS_END, dummy, LIS_END, idummy, idummy, dummy, idummy,
     1           idummy, IERR)
         elseif( QSYMSOLV.eq.2 ) then
            call MA57_CALL(0, N+M, NZKKT, dummy, idummy, idummy, dummy,
     1           dummy, idummy, ldummy, idummy, idummy, ldummy,
     1           idummy, idummy,
     1           LRS_END, dummy, LIS_END, idummy, idummy, dummy, idummy,
     1           idummy, IERR)
         elseif( QSYMSOLV.eq.5 .or. QSYMSOLV.eq.6 ) then
            call MA48SYM_CALL(0, N+M, NZKKT, dummy, idummy, idummy,
     1           dummy, idummy, ldummy, idummy, idummy, ldummy, idummy,
     1           idummy, LRS_END, dummy, LIS_END, idummy, idummy, dummy,
     1           idummy, idummy, IERR)
         elseif( QSYMSOLV.eq.3 ) then
            call WSSMP_CALL(0, N+M, NZKKT, dummy, idummy, idummy, dummy,
     1           idummy, ldummy, idummy, idummy, ldummy, idummy, idummy,
     1           LRS_END, dummy, LIS_END, idummy, idummy, dummy, idummy,
     1           idummy, IERR)
         elseif( QSYMSOLV.eq.4 ) then
            call WGSMP_CALL(0, N+M, NZKKT, dummy, idummy, idummy, dummy,
     1           idummy, ldummy, idummy, idummy, ldummy, idummy, idummy,
     1           LRS_END, dummy, LIS_END, idummy, idummy, dummy, idummy,
     1           idummy, IERR)
         else
            call C_OUT(2,0,1,'get_step_full: Invalid QSYMSOLV.')
            IERR = 4
            goto 9999
         endif

         P_JCN   = LIS_END
         P_IRN   = P_JCN + NZKKT
         LIS_END = P_IRN + NZKKT

         HAVE_PIV    = .false.
         L_INC_LA    = .false.
         L_INC_LIWMA = .false.
         LA          = NZKKT
         LIWMA       = 0
         REGU_STORE = 0d0
         CHECKDEPENCOUNTER = 0
#ifdef USE_MALLOC
         PTR_KKT  = 0
         PTR_IWMA = 0
#endif

         if( ITERDEP_MAX.gt.0 ) then
            if( M.gt.0 .and. QDEPCONDIAG.ne.0.d0 ) then
               POSSIBLE_DEPCON = .true.
            else
               POSSIBLE_DEPCON = .false.
            endif
            POSSIBLE_DEPHES = .true.
         else
            POSSIBLE_DEPCON = .false.
            POSSIBLE_DEPHES = .false.
         endif
         IS_DEPCON = .false.
         IS_DEPHES = .false.
         ITERDEP = 0

         MACHEPS  = D1MACH(4)
         goto 9999
      endif
C
      p_rwend = 0
      p_iwend = 0

      if( MEMDBG .and. QPRINT.ge.3 ) then
         write(line,1)'get_step_full', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      if( abs(QQUASI).ge.6 ) then
C
C     For limited memory quasi-Newton: Get the structure of B right
C
         nlm       = int(B(1))
         nlmmax    = QLMLEN
         sigmak    = B(2)
         pb_uk     = 2
         pb_vk     = pb_uk + N*nlmmax
      endif
C
      l_cs         = .false.
      l_iter_ref   = .false.
      skip_jac     = .false.
      l_lbfgsreset = .false.
      refinement_iter = 0
      is_dep_this_iter = .false.
      naninf = .false.
C
C     Initialize
C
      if( SOC_FLAG.eq.0 .and. RESTO.ne.2 ) then
         if( RESTO.eq.3 ) then
            if( QPIVTOL.lt.QPIVTOLMAX .and. .not.L_INCPIVTOL ) then
               QPIVTOL = min((QPIVTOL**0.75d0), QPIVTOLMAX)
               write(line,*) 'get_step_full: increase QPIVTOL to ',
     1              QPIVTOL
               call C_OUT(2,2,1,line)
               L_INCPIVTOL   = .true.
            elseif( REGU.lt.REGU_MAX ) then
               if( REGU.eq.0d0 ) then
                  COUNT_NEG_CURV = COUNT_NEG_CURV + 1
                  if( REGU_STORE.eq.0d0 ) then
                     REGU = REGU_INIT
                  else
                     REGU = REGU_DEC_FACT*REGU_STORE
                  endif
               else
                  if( REGU_STORE.eq.0d0 ) then
                     REGU = REGU_INIT_FACT*REGU
                  else
                     REGU = REGU_INC_FACT*REGU
                  endif
               endif
               IEIGS = IEIGS + 1
            else
               write(line,*)
     1'Regularization parameter getting too large in restoration phase:'
     2, REGU
               call C_OUT(2,1,1,line)
               IERR = 10
               goto 9999
            endif
            RESTO = 1
         else
            REGU   = 0.d0
            L_INCPIVTOL   = .false.
            IEIGS  = 0
         endif

         DHD    = 0.d0
         DEPCONDIAGVAL = 0.d0
C
C     if necessary, increase LA or LIWMA
C
         if( L_INC_LA ) then
            if( QPRINT.ge.2 ) then
               call C_OUT(2,0,1,'get_step_full: Double LA.')
            endif
            LA = 2*LA
#ifdef USE_MALLOC
            type = 1
            i = IP_MALLOC(type, LA, PTR_KKT)
            if( i.ne.0 ) then
               write(line,*) 'get_step_full: IP_MALLOC returned ',
     1              i,' for KKT (L_INC_LA)'
               call C_OUT(2,0,1,line)
               IERR = 96
               goto 9999
            endif
#endif
         endif
         if( L_INC_LIWMA ) then
            if( QPRINT.ge.2 ) then
               call C_OUT(2,0,1,'get_step_full: Double LIWMA.')
            endif
            LIWMA = 2*LIWMA
#ifdef USE_MALLOC
            type = 2
            i = IP_MALLOC(type, LIWMA, PTR_IWMA)
            if( i.ne.0 ) then
               write(line,*) 'get_step_full: IP_MALLOC returned ',
     1              i,' for IWMA (L_INC_LIWMA)'
               call C_OUT(2,0,1,line)
               IERR = 96
               goto 9999
            endif
#endif
         endif
         L_INC_LA    = .false.
         L_INC_LIWMA = .false.

         REFINEMENTHELPLESS = 0
         if( WFLAG.eq.1 ) then
            REGU_STORE_OLD = REGU_STORE
            DEPCONDIAGVAL_OLD = DEPCONDIAGVAL
         elseif( WFLAG.eq.2 ) then
            REGU_STORE = REGU_STORE_OLD
            DEPCONDIAGVAL = DEPCONDIAGVAL_OLD
         endif

CTODO Check if this is a good idea (didn't help on ORTHREGE)
         if( QDEPCONDIAG.gt.1.d-6 ) then
            QDEPCONDIAG = QDEPCONDIAG*1.d-1
            write(line,*) 'QDEPCONDIAG decreased: ',QDEPCONDIAG
            call C_OUT(2,1,1,line)
         endif

      elseif( SOC_FLAG.eq.1 .or. RESTO.eq.2 ) then
C
C     Compute SOC step - reuse old factorization and scaling factors
C     reassigne work space pointers
C
         p_kkt   = p_rwend
#ifdef USE_MALLOC
         p_rwend = p_kkt + NZKKT
#else
         p_rwend = p_kkt + LA
#endif
         p_scale = p_rwend
         if( QKKTSCALE.eq.0 ) then
            p_rwend = p_scale
         else
            p_rwend = p_scale + M+N
         endif
C         if( abs(QMERIT).eq.1 .or. abs(QMERIT).eq.2 ) then
C#ifdef USE_MALLOC
            p_hcopy = p_kkt
C#else
C            p_hcopy = p_rwend
C            p_rwend = p_hcopy + NZKKT
C#endif
C         endif
         if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
            write(line,*) 'Out of DP memory in line ',__LINE__
            call C_OUT(2,0,1,line)
#endif
            IERR = 98
            goto 9999
         endif
         p_iwma  = p_iwend
#ifndef USE_MALLOC
         p_iwend = p_iwma + LIWMA
         if( p_iwend.gt.LIW ) then
#ifdef PRINT__LINE__
            write(line,*) 'Out of INT memory in line ',__LINE__
            call C_OUT(2,0,1,line)
#endif
            IERR = 99
            goto 9999
         endif
#endif
#ifdef USE_MALLOC
         p_rhscopy = p_rwend
#else
         p_hcopy = p_rwend
#endif
         p_rhscopy = p_rwend
         p_rwend   = p_rhscopy + N+M+NLB+NUB
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         goto 1050
      endif
C
C     If basis changed, need to call MA?7AD again
C
      if( NEWBAS ) then
         HAVE_PIV = .false.
      endif
C
 1000 continue
C
C     Fill the KKT matrix
C
      p_kkt   = p_rwend
#ifdef USE_MALLOC
      p_rwend = p_kkt + NZKKT
#else
      p_rwend = p_kkt + LA      ! If HAVE_PIV = .false. we have to correct
                                ! p_rwend lateron!
#endif
      if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
         write(line,*) 'Out of DP memory in line ',__LINE__
         call C_OUT(2,0,1,line)
#endif
         IERR = 98
         goto 9999
      endif
C
C     Build full space KKT matrix
C
      inz = 0
C
C     1. Diagonal in Hessian for regularization and barrier term
C
      if( RESTO.eq.1 ) then
         call DCOPY(N, RHS_RESTO(N+M+1), 1, RW(p_kkt+inz+1), 1)
         call DSCAL(N, ETA_RESTO, RW(p_kkt+inz+1), 1)
         if( abs(QQUASI).ge.6 ) then
            call DAXPY(N, 1.d0, sigmak, 0, RW(p_kkt+inz+1), 1)
         endif
         if( REGU.gt.0d0 ) then
            CALL DAXPY(N, 1d0, REGU, 0, RW(p_kkt+inz+1), 1)
            if( CORRECTTYPE.eq.1 ) then
               do i = 1, NLB
                  k = ILB(i)
                  RW(p_kkt+inz+k) = min(RW(p_kkt+inz+k),
     1                 REGU/(S_L(i)**2))
               enddo
               do i = 1, NUB
                  k = IUB(i)
                  RW(p_kkt+inz+k) = min(RW(p_kkt+inz+k),
     1                 REGU/(S_U(i)**2))
               enddo
            endif
         endif
      elseif( INIT_LAM ) then
         call DCOPY(N, 1.d0, 0, RW(p_kkt+inz+1), 1)
      else
         if( ITERDEP_MAX.gt.0 .and. IS_DEPHES ) then
            if( REGU_STORE.eq.0d0 ) then
CTAKEMEOUT
               REGU = REGU_INIT
               IERR = 32884
               goto 9999
            else
               REGU = REGU_DEC_FACT*REGU_STORE
            endif
         endif
         call DCOPY(N, REGU, 0, RW(p_kkt+inz+1), 1)
         if( CORRECTTYPE.eq.1 ) then
            do i = 1, NLB
               k = ILB(i)
               RW(p_kkt+inz+k) = min(RW(p_kkt+inz+k),
     1             REGU/(S_L(i)**2))
            enddo
            do i = 1, NUB
               k = IUB(i)
               RW(p_kkt+inz+k) = min(RW(p_kkt+inz+k),
     1              REGU/(S_U(i)**2))
            enddo
         endif
         if( abs(QQUASI).ge.6 ) then
            call DAXPY(N, 1.d0, sigmak, 0, RW(p_kkt+inz+1), 1)
         endif
      endif
C
C     Barrier Hessian term
C
      if( (NIND.gt.0.and..not.INIT_LAM )
     1     .or. RESTO.eq.1 ) then
         do i = 1, NLB
            k = ILB(i)
            RW(p_kkt+inz+k) = RW(p_kkt+inz+k) + SIGMA_L(i)
         enddo
         do i = 1, NUB
            k = IUB(i)
            RW(p_kkt+inz+k) = RW(p_kkt+inz+k) + SIGMA_U(i)
         enddo
      endif

      if( .not.HAVE_PIV ) then
         do i = 1, N
            IS(P_JCN+i) = i
            IS(P_IRN+i) = i
         enddo
      endif
      inz = inz + N
C
C     2. Hessian matrix
C
      if( QQUASI.ne.0 .and. abs(QQUASI).lt.6 ) then
         if( RESTO.le.0 .and. .not.INIT_LAM .and. NIND.gt.0 ) then
            call DCOPY( (NIND*(NIND+1))/2, B, 1, RW(p_kkt+inz+1), 1)
         else
            call DCOPY( (NIND*(NIND+1))/2, 0d0, 0,
     1           RW(p_kkt+inz+1), 1)
         endif
         if( .not.HAVE_PIV ) then
            do i = 1, NIND
               do j = 1, i
                  inz = inz + 1
                  IS(P_JCN+inz) = i+M
                  IS(P_IRN+inz) = j+M
               enddo
            enddo
         else
            inz = inz + (NIND*(NIND+1))/2
         endif
      elseif( QQUASI.eq.0 ) then
         if( QLAMBDA.eq.0 ) then
            call C_OUT(2,0,1,
     1           'Need multipliers for exact Hessian! Abort')
            IERR = 4
            goto 9999
         endif
         if( (RESTO.le.0 .and. .not.INIT_LAM) .or. .not.HAVE_PIV ) then
            call GET_H(ITER, N, NIND, NFIX, X, IVAR, NORIG, XORIG,
     1           CSCALE, M, NNZH, LAM, RW(p_kkt+inz+1), IS(P_IRN+inz+1),
     1           IS(P_JCN+inz+1), LRW-p_rwend, RW(p_rwend+1),
     1           LIW-p_iwend, IW(p_iwend+1), IERR, EV_H, DAT, IDAT)
            if( IERR.gt.0 ) then
               write(line,*) 'get_step_full: GET_H returns IERR =',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
         elseif( RESTO.eq.1 ) then
CTODO only need to compute Hessian if LAMBDA is nonzero!
            p_tmp   = p_rwend
            p_w     = p_tmp + NNZH
            p_rwend = p_w + M
            if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
               write(line,*) 'Out of DP memory in line ',__LINE__
               call C_OUT(2,0,1,line)
#endif
               IERR = 98
               goto 9999
            endif
            call GET_H(ITER, N, NIND, NFIX, X, IVAR, NORIG, XORIG,
     1           CSCALE, M, NNZH, LAM, RW(p_kkt+inz+1),
     1           IS(P_IRN+inz+1),
     1           IS(P_JCN+inz+1), LRW-p_rwend, RW(p_rwend+1),
     1           LIW-p_iwend, IW(p_iwend+1), IERR, EV_H, DAT, IDAT)
            if( IERR.ne.0 ) then
               write(line,*)
     1              'get_step_full-r2: GET_H returns IERR =',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            call DCOPY(M, 0.d0, 0, RW(p_w+1), 1)
            call GET_H(ITER, N, NIND, NFIX, X, IVAR, NORIG, XORIG,
     1           CSCALE, M, NNZH, RW(p_w+1), RW(p_tmp+1),
     1           IS(P_IRN+inz+1),
     1           IS(P_JCN+inz+1), LRW-p_rwend, RW(p_rwend+1),
     1           LIW-p_iwend, IW(p_iwend+1), IERR, EV_H, DAT, IDAT)
            if( IERR.ne.0 ) then
               write(line,*)
     1              'get_step_full-r2: GET_H returns IERR =',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            call DAXPY(NNZH, -1.d0, RW(p_tmp+1), 1, RW(p_kkt+inz+1), 1)
            p_rwend = p_tmp
         endif
         if( (INIT_LAM .or. NIND.eq.0) .and. RESTO.le.0 ) then
            call DCOPY(NNZH, 0d0, 0, RW(p_kkt+inz+1), 1)
         endif
         inz = inz + NNZH
      endif
C
C     3. Jacobian of constraints
C
      if( M.gt.0 ) then
         if( .not.HAVE_PIV ) then
            call CONSTR(10, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, vin,
     2           RW(p_kkt+inz+1), IS(P_JCN+inz+1), IS(P_IRN+inz+1),
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.ne.0 ) goto 9999
            NZA = int(vin(1))
c$$$CDELETEME
c$$$            do i = 1, NZA
c$$$               write(QCNR,1442) i, IS(P_JCN+inz+i), IS(P_IRN+inz+i),
c$$$     1              RW(p_kkt+inz+i)
c$$$ 1442          format(i6,i6,i6,d25.16)
c$$$            enddo
C
            do i = 1, NZA       ! Correct row indices
               IS(P_JCN+inz+i) = IS(P_JCN+inz+i) + N
            enddo
         elseif( .not.skip_jac ) then
            p_jcntmp = p_iwend
            p_iwend  = p_jcntmp + NZA
            if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
               write(line,*) 'Out of DP memory in line ',__LINE__
               call C_OUT(2,0,1,line)
#endif
               IERR = 98
               goto 9999
            elseif( p_iwend.gt.LIW ) then
#ifdef PRINT__LINE__
               write(line,*) 'Out of INT memory in line ',__LINE__
               call C_OUT(2,0,1,line)
#endif
               IERR = 99
               goto 9999
            endif
            call CONSTR(10, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, vin,
     2           RW(p_kkt+inz+1), IW(p_jcntmp+1), IS(P_IRN+inz+1),
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.ne.0 ) goto 9999
            p_iwend = p_jcntmp
         endif
      else
         NZA = 0
      endif
C
      inz = inz + NZA
C
C     4. Diagonal in right lower corner to deal with dependent constraints
C
      if( QDEPCONDIAG.gt.0.d0 .or. QSYMSOLV.eq.3 .or. QSYMSOLV.eq.4 .or.
     1     (abs(QMERIT).ge.4.and.QRESTO.eq.1) ) then
         if( RESTO.eq.1 ) then
            call DCOPY(M, SIGMA_PN, 1, RW(p_kkt+inz+1), 1)
         else
            if( ITERDEP_MAX.gt.0 ) then
               if( IS_DEPCON .and. .not.INIT_LAM ) then
                  DEPCONDIAGVAL = -QDEPCONDIAG
               endif
               constrregu = min(-QCONSTRRELAX, DEPCONDIAGVAL)*
     1              MU**QMUDEPEXP
            else
               constrregu = min(-QCONSTRRELAX, DEPCONDIAGVAL)*
     1              MU**QMUDEPEXP
            endif
            call DCOPY(M, constrregu, 0, RW(p_kkt+inz+1), 1)
            if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
               write(line,*) 'Lower left corner in matrix = ',constrregu
               call C_OUT(1,0,1,line)
            endif
         endif
         if( .not.HAVE_PIV ) then
            do i = 1, M
               IS(P_JCN+inz+i) = N + i
            enddo
            do i = 1, M
               IS(P_IRN+inz+i) = N + i
            enddo
         endif
         inz = inz + M
      endif
C
C     Obtain pivot sequence if we do not have it yet
C
      if( .not.HAVE_PIV ) then
C
C     Correct NZKKT (in case of fixed vars...) and reallocate memory
C
         NZKKT   = inz
         p_rwend = p_kkt + NZKKT ! correct memory allocation
C
C     Now get pivots from MA47AD
C
         if( QSYMSOLV.eq.0 ) then
            call MA47_CALL(1, N+M, NZKKT, RW(p_kkt+1), IS(P_IRN+1),
     1           IS(P_JCN+1), dummy, LIWMA, L_INC_LIWMA, idummy, LA,
     2           L_INC_LA, neigs, nzero,
     3           LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR)
         elseif( QSYMSOLV.eq.1 ) then
            call MA27_CALL(1, N+M, NZKKT, RW(p_kkt+1), IS(P_IRN+1),
     1           IS(P_JCN+1), dummy, LIWMA, L_INC_LIWMA, idummy, LA,
     2           L_INC_LA, neigs, nzero,
     3           LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR)
         elseif( QSYMSOLV.eq.2 ) then
            call MA57_CALL(1, N+M, NZKKT, dummy, IS(P_IRN+1),
     1           IS(P_JCN+1), dummy, dummy, LA, L_INC_LA, idummy, LIWMA,
     2           L_INC_LIWMA, neigs, nzero,
     3           LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR)
         elseif( QSYMSOLV.eq.5 .or. QSYMSOLV.eq.6 ) then
            call MA48SYM_CALL(1, N+M, NZKKT, RW(p_kkt+1), IS(P_IRN+1),
     1           IS(P_JCN+1), dummy, LIWMA, L_INC_LIWMA, idummy, LA,
     2           L_INC_LA, neigs, nzero,
     3           LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR)
         elseif( QSYMSOLV.eq.3 ) then
            neigs = N
            call WSSMP_CALL(1, N+M, NZKKT, RW(p_kkt+1), IS(P_IRN+1),
     1           IS(P_JCN+1), dummy, LIWMA, L_INC_LIWMA, idummy, LA,
     2           L_INC_LA, neigs, nzero,
     3           LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR)
         elseif( QSYMSOLV.eq.4 ) then
            call WGSMP_CALL(1, N+M, NZKKT, RW(p_kkt+1), IS(P_IRN+1),
     1           IS(P_JCN+1), dummy, LIWMA, L_INC_LIWMA, idummy, LA,
     2           L_INC_LA, neigs, nzero,
     3           LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR)
         endif
         if( IERR.ne.0 ) then
            write(line,*) 'get_step_full: MA?7_CALL(1) returns IERR =',
     1           IERR
            call C_OUT(2,0,1,line)
            if( IERR.gt.99 ) IERR = 505
            goto 9999
         endif
C
C     Check if we should stop
C
         IERR = CHECK_STOP()
         if( IERR.ne.0 ) goto 9999
C
C     Correct storage space for KKT matrix
C
#ifdef USE_MALLOC
         type = 1
         i = IP_MALLOC(type, LA, PTR_KKT)
         if( i.ne.0 ) then
            write(line,*) 'get_step_full: IP_MALLOC returned ',
     1           i,' for KKT'
            call C_OUT(2,0,1,line)
            IERR = 96
            goto 9999
         endif
         type = 2
         i = IP_MALLOC(type, LIWMA, PTR_IWMA)
         if( i.ne.0 ) then
            write(line,*) 'get_step_full: IP_MALLOC returned ',
     1           i,' for IWMA'
            call C_OUT(2,0,1,line)
            IERR = 96
            goto 9999
         endif

#else
         p_rwend = p_kkt + LA
         if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
            write(line,*) 'Out of DP memory in line ',__LINE__
            call C_OUT(2,0,1,line)
#endif
            IERR = 98
            goto 9999
         endif
#endif
C
         HAVE_PIV = .true.

      endif
C
C     Reserve memory
C
      p_scale   = p_rwend
      if( QKKTSCALE.eq.0 ) then
         p_rwend = p_scale
      else
         p_rwend = p_scale + M+N
      endif
#ifdef USE_MALLOC
      p_hcopy   = p_kkt
      p_rhscopy = p_rwend
#else
      p_hcopy   = p_rwend
      p_rhscopy = p_hcopy   + NZKKT
#endif
      p_rwend   = p_rhscopy + N+M+NLB+NUB
      if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
         write(line,*) 'Out of DP memory in line ',__LINE__
         call C_OUT(2,0,1,line)
#endif
         IERR = 98
         goto 9999
      endif
C
C     Let calling routine know much to remember for another SOC call
CTODO: CAN THIS BE ZERO FOR RESTO ?
C
CTODO: Do not delete hcopy for RESTO if we use iterative refinement
      if( abs(QMERIT).eq.1 .or. abs(QMERIT).eq.2 ) then
         PRWKEEP = p_rhscopy
      else
C#ifdef USE_MALLOC
         PRWKEEP = p_rhscopy
C#else
C         PRWKEEP = p_hcopy
C#endif
      endif
#ifdef USE_MALLOC
      PIWKEEP = 0
#else
      PIWKEEP = LIWMA
#endif
C
C     Keep a copy of the elements of the KKT matrix
C
#ifdef USE_MALLOC
      call DCOPY(NZKKT, RW(p_kkt+1), 1, %VAL(PTR_KKT), 1)
#else
      call DCOPY(NZKKT, RW(p_kkt+1), 1, RW(p_hcopy+1), 1)
#endif
 1500 continue
C
C     Check, whether there is a NaN in the KKT matrix
C     (otherwise MC30 might hang)
C
      dummy = DASUM(NZKKT, RW(p_kkt+1), 1)
      if( FFINITE(dummy).eq.0 ) then
         if( WFLAG.eq.1 ) then
            ALPHA = 0.d0
            goto 3000
         endif
         call C_OUT(2,0,1,
     1        'get_step_full: NaN or Inf found in KKT matrix.')
         IERR = 11
         goto 9999
      endif
C
C     Compute scaling parameters if desired
C
      if( QKKTSCALE.eq.1 ) then
C
C     call MC30, MC29 or MC19 to obtain scaling parameters
C
#ifdef HAVE_MC30
         p_w     = p_rwend
         p_rwend = p_w     + 4*(N+M)
         if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
            write(line,*) 'Out of DP memory in line ',__LINE__,
     1           'p_rwend = ',p_rwend
            call C_OUT(2,0,1,line)
#endif
            IERR = 98
            goto 9999
         endif
         call MC30AD(N+M, NZKKT, RW(p_kkt+1), IS(P_IRN+1),
     1        IS(P_JCN+1), RW(p_scale+1), RW(p_w+1), QCNR, ifail)
         p_rwend = p_w
         if( ifail.ne.0 ) then
            write(line,*) 'MC30AD returns ifail = ',ifail
            call C_OUT(2,0,1,line)
            IERR = 517
            goto 9999
         endif
#else
#ifdef HAVE_MC29
         call MC29_CALL(N+M, NZKKT, RW(p_kkt+1), IS(P_IRN+1),
     1        IS(P_JCN+1), RW(p_scale+1), LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR)
         if( IERR.ne.0 ) then
            write(line,*) 'MC29call returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
#else
#ifdef HAVE_MC19
C
         p_realw  = p_rwend     ! This will be REAL work space in MC19_CALL
         p_rwend  = p_realw + (7*(N+M))/2 + 1
         if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
            write(line,*) 'Out of DP memory in line ',__LINE__
            call C_OUT(2,0,1,line)
#endif
            IERR = 98
            goto 9999
         endif
         call MC19_CALL(N+M, NZKKT, RW(p_kkt+1), IS(P_IRN+1),
     1        IS(P_JCN+1), RW(p_scale+1), RW(p_realw+1),
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR)
         p_rwend = p_realw
         if( IERR.ne.0 ) then
            write(line,*) 'MC19call returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
#else
         IERR = 97
         call C_OUT(2,0,1,
     1        'get_step_full:  Compiled without MC30, MC29 or MC19.')
         call C_OUT(2,0,1,
     1        '                Choose a different value for IKKTSCALE.')
         goto 9999
#endif
#endif
#endif
         do i = 1, N+M
            RW(p_scale+i) = exp(RW(p_scale+i))
         enddo
      elseif( QKKTSCALE.eq.2 ) then
C
C     Use affine scaling based on slacks
C
         call DCOPY(N+M, 1.d0, 0, RW(p_scale+1), 1)
         do i = 1, NLB
            j = ILB(i)
            RW(p_scale+j) = dmin1(RW(p_scale+j), S_L(i))
         enddo
         do i = 1, NUB
            j = IUB(i)
            RW(p_scale+j) = dmin1(RW(p_scale+j), S_U(i))
         enddo
      elseif( QKKTSCALE.ne.0 ) then
         call C_OUT(2,0,1,'get_step_full: Invalid value for QKKTSCALE.')
         IERR = 4
         goto 9999
      endif
C
C     Scale the matrix
C
 1600 continue
      if( QKKTSCALE.ne.0 ) then
         call GSF_SCALE(N+M, NZKKT, IS(P_IRN+1), IS(P_JCN+1),
     1        RW(p_scale+1),
#ifdef USE_MALLOC
     1        %VAL(PTR_KKT))
#else
     1        RW(p_kkt+1))
#endif
      endif
C
C     Now let us do the factorization
C
      p_iwma  = p_iwend
#ifndef USE_MALLOC
      p_iwend = p_iwma + LIWMA
      if( p_iwend.gt.LIW ) then
#ifdef PRINT__LINE__
         write(line,*) 'Out of INT memory in line ',__LINE__
         call C_OUT(2,0,1,line)
#endif
         IERR = 99
         goto 9999
      endif
#endif
      if( PRINTTIMING .and. QPRINT.ge.2 ) then
         call timer(dummy)
         write(line,*) 'Time before factorization = ',dummy
         call C_OUT(2,0,1,line)
      endif
      if( QSYMSOLV.eq.0 ) then
         call MA47_CALL(2, N+M, NZKKT,
#ifdef USE_MALLOC
     1        %VAL(PTR_KKT),
#else
     1        RW(p_kkt+1),
#endif
     1        IS(P_IRN+1), IS(P_JCN+1), dummy, LIWMA, L_INC_LIWMA,
#ifdef USE_MALLOC
     1        %VAL(PTR_IWMA),
#else
     1        IW(p_iwma+1),
#endif
     1        LA, L_INC_LA, neigs, nzero,
     1        LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     1        IW(p_iwend+1), IERR)
      elseif( QSYMSOLV.eq.1 ) then
         call MA27_CALL(2, N+M, NZKKT,
#ifdef USE_MALLOC
     1        %VAL(PTR_KKT),
#else
     1        RW(p_kkt+1),
#endif
     1        IS(P_IRN+1), IS(P_JCN+1), dummy, LIWMA, L_INC_LIWMA,
#ifdef USE_MALLOC
     1        %VAL(PTR_IWMA),
#else
     1        IW(p_iwma+1),
#endif
     1        LA, L_INC_LA, neigs, nzero,
     1        LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     1        IW(p_iwend+1), IERR)
      elseif( QSYMSOLV.eq.2 ) then
         p_a     = p_rwend
         p_rwend = p_a + NZKKT
         if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
            write(line,*) 'Out of DP memory in line ',__LINE__
            call C_OUT(2,0,1,line)
#endif
            IERR = 98
            goto 9999
         endif
#ifdef USE_MALLOC
         call DCOPY(NZKKT, %VAL(PTR_KKT), 1, RW(p_a+1), 1)
#else
         call DCOPY(NZKKT, RW(p_kkt+1), 1, RW(p_a+1), 1)
#endif
         call MA57_CALL(2, N+M, NZKKT, RW(p_a+1), IS(P_IRN+1),
     1        IS(P_JCN+1), dummy,
#ifdef USE_MALLOC
     1        %VAL(PTR_KKT),
#else
     1        RW(p_kkt+1),
#endif
     1        LA, L_INC_LA,
#ifdef USE_MALLOC
     1        %VAL(PTR_IWMA),
#else
     1        IW(p_iwma+1),
#endif
     1        LIWMA, L_INC_LIWMA, neigs, nzero,
     3        LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR)
         p_rwend = p_a
      elseif( QSYMSOLV.eq.5 .or. QSYMSOLV.eq.6 ) then
         call MA48SYM_CALL(2, N+M, NZKKT,
#ifdef USE_MALLOC
     1        %VAL(PTR_KKT),
#else
     1        RW(p_kkt+1),
#endif
     1        IS(P_IRN+1), IS(P_JCN+1), dummy, LIWMA, L_INC_LIWMA,
#ifdef USE_MALLOC
     1        %VAL(PTR_IWMA),
#else
     1        IW(p_iwma+1),
#endif
     1        LA, L_INC_LA, neigs, nzero,
     1        LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     1        IW(p_iwend+1), IERR)
         neigs = M
      elseif( QSYMSOLV.eq.3 ) then
         call WSSMP_CALL(2, N+M, NZKKT,
#ifdef USE_MALLOC
     1        %VAL(PTR_KKT),
#else
     1        RW(p_kkt+1),
#endif
     1        IS(P_IRN+1), IS(P_JCN+1), dummy, LIWMA, L_INC_LIWMA,
#ifdef USE_MALLOC
     1        %VAL(PTR_IWMA),
#else
     1        IW(p_iwma+1),
#endif
     1        LA, L_INC_LA, neigs, nzero,
     1        LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     1        IW(p_iwend+1), IERR)
      elseif( QSYMSOLV.eq.4 ) then
         call WGSMP_CALL(2, N+M, NZKKT,
#ifdef USE_MALLOC
     1        %VAL(PTR_KKT),
#else
     1        RW(p_kkt+1),
#endif
     1        IS(P_IRN+1), IS(P_JCN+1), dummy, LIWMA, L_INC_LIWMA,
#ifdef USE_MALLOC
     1        %VAL(PTR_IWMA),
#else
     1        IW(p_iwma+1),
#endif
     1        LA, L_INC_LA, neigs, nzero,
     1        LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     1        IW(p_iwend+1), IERR)
         neigs = M
      endif
      if( PRINTTIMING .and. QPRINT.ge.2 ) then
         call timer(dummy)
         write(line,*) 'Time after factorization = ',dummy
         call C_OUT(2,0,1,line)
      endif
C
C     Check if we should stop
C
      i = CHECK_STOP()
      if( i.ne.0 ) then
         IERR = i
         goto 9999
      endif

      if( IERR.eq.-1 ) then
C
C     LIWMA or LA too small - increase storage for KKT matrix factorization,
C     rearrange memory, and factorize again
C
#ifdef USE_MALLOC
         type = 1
         i = IP_MALLOC(type, LA, PTR_KKT)
         if( i.ne.0 ) then
            write(line,*) 'get_step_full: IP_MALLOC returned ',
     1           i,' for KKT (realloc)'
            call C_OUT(2,0,1,line)
            IERR = 96
            goto 9999
         endif
         type = 2
         i = IP_MALLOC(type, LIWMA, PTR_IWMA)
         if( i.ne.0 ) then
            write(line,*) 'get_step_full: IP_MALLOC returned ',
     1           i,' for IWMA (realloc)'
            call C_OUT(2,0,1,line)
            IERR = 96
            goto 9999
         endif
#else
         p_iwend     = p_iwma

         p_scale_old = p_scale
         p_rwend     = p_kkt + LA

         inz       = NZKKT
         p_scale   = p_rwend
         if( QKKTSCALE.eq.0 ) then
            p_rwend = p_scale
         else
            p_rwend = p_scale + M+N
            inz     = inz + M+N
         endif
         p_hcopy   = p_rwend
         p_rhscopy = p_hcopy   + NZKKT
         p_rwend   = p_rhscopy + N+M+NLB+NUB
         if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
            write(line,*) 'Out of DP memory in line ',__LINE__
            call C_OUT(2,0,1,line)
#endif
            IERR = 98
            goto 9999
         endif
C
C     Move scale and hcopy
C
         do i = inz, 1, -1      ! I'm not sure if we can use DCOPY here if
                                ! the two arrays overlap...?
            RW(p_scale+i) = RW(p_scale_old+i)
         enddo
#endif
#ifdef USE_MALLOC
         call DCOPY(NZKKT, RW(p_hcopy+1), 1, %VAL(PTR_KKT), 1)
#else
         call DCOPY(NZKKT, RW(p_hcopy+1), 1, RW(p_kkt+1), 1)
#endif
         goto 1600

      elseif( IERR.ne.0 ) then
         write(line,*) 'get_step_full: MA?7_CALL(2) returns IERR =',IERR
         call C_OUT(2,0,1,line)
         if( IERR.gt.99 ) IERR = 506
         goto 9999
      endif
C
C     Sometimes, the linear solver reports too few negative eigenvalues.
C     In this case, we increase the pivtol tolerance, and try another
C     factorization.  If that doesn't help, pretend that the system is
C     singular.
C
      if( neigs.lt.M .and. nzero.eq.0 ) then
         write(line,*) 'get_step_full: Inertia is neigs =',neigs,
     1        ', nzero =',nzero,' M = ',M
         call C_OUT(1,2,1,line)
         if( QPRINT.ge.3 ) then
            call C_OUT(1,0,1,'Detected TOO FEW negative eigenvalues.')
         endif
         if( QPIVTOL.lt.QPIVTOLMAX .and. .not.L_INCPIVTOL ) then
            QPIVTOL = min((QPIVTOL**0.75d0), QPIVTOLMAX)
            write(line,*) 'get_step_full: increase QPIVTOL to ',
     1           QPIVTOL
            call C_OUT(2,2,1,line)
#ifdef USE_MALLOC
               call DCOPY(NZKKT, RW(p_hcopy+1), 1, %VAL(PTR_KKT), 1)
#else
               call DCOPY(NZKKT, RW(p_hcopy+1), 1, RW(p_kkt+1), 1)
#endif
            L_INCPIVTOL = .true.      ! Make sure we increase QPIVTOL only once
                                      ! per iteration
            p_iwend = p_iwma
            p_rwend   = p_rhscopy + N+M+NLB+NUB
            goto 1600
         else
            nzero = 1
         endif
      endif
C
C     The following pretends that the system is singular, if the solver
C     finds TOO FEW negative eigenvalues
C
      if( neigs.lt.M ) nzero = 1
C
C     Check inertia
C
      if( ITERDEP_MAX.gt.0 .and. .not.INIT_LAM ) then
         if( SOC_FLAG.eq.0 .and. RESTO.le.0 .and. nzero.eq.0 ) then
            if( POSSIBLE_DEPCON .and. DEPCONDIAGVAL.eq.0.d0 ) then
               POSSIBLE_DEPCON = .false.
               call C_OUT(2,2,1,'Constraints not dependent.')
            endif
            if( POSSIBLE_DEPHES .and. REGU.eq.0.d0 ) then
               POSSIBLE_DEPHES = .false.
               call C_OUT(2,2,1,'Hessian not dependent.')
            endif
         endif
      endif
C
 1030 continue
      if( (neigs.ne.M .or. nzero.ne.0) .and. SOC_FLAG.eq.0 ) then
C
C     Check if constraints or Hessian seem always dependent
C
         if( ITERDEP_MAX.gt.0 .and. .not.INIT_LAM .and. RESTO.le.0
     `        .and. nzero.ne.0 ) then
            if( REGU.eq.0.d0 .and. DEPCONDIAGVAL.eq.0.d0 .and.
     1           (POSSIBLE_DEPCON.or.POSSIBLE_DEPHES) ) then
               ITERDEP = ITERDEP + 1
               if( ITERDEP.ge.ITERDEP_MAX ) then
                  if( POSSIBLE_DEPCON ) then
                     IS_DEPCON = .true.
                     POSSIBLE_DEPCON = .false.
                     call C_OUT(2,0,1,
     1                    'Constraints seem rank-deficient!')
                  endif
                  if( POSSIBLE_DEPHES ) then
                     IS_DEPHES = .true.
                     POSSIBLE_DEPHES = .false.
                     call C_OUT(2,0,1,'Hessian seems rank-deficient!')
                  endif
               endif
            endif
         endif
         write(line,*) 'get_step_full: Inertia is neigs =',neigs,
     1        ', nzero =',nzero,' M = ',M
         call C_OUT(1,2,1,line)
         goto 5000              ! Modify KKT matrix
      endif
C
      write(line,1010) REGU
 1010 format(/,'get_step_full: Regularization is REGU = ',d12.5)
      call C_OUT(1,2,2,line)
C
CTODO Do not know yet if this is good (seems better without)
c      if( INIT_LAM .and. DEPCONDIAGVAL.eq.0.d0 )
cc     1     write(*,*) 'Decide about INIT_LAM and POSSIBLE_DEPCON'
c     1     POSSIBLE_DEPCON = .false.

 1050 continue
C
C     Compute RHS
C
      p_rhs   = p_rwend
      p_rwend = p_rhs + N + M + NLB + NUB
      if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
         write(line,*) 'Out of DP memory in line ',__LINE__
         call C_OUT(2,0,1,line)
#endif
         IERR = 98
         goto 9999
      endif
      if( RESTO.gt.0 ) then
         call DCOPY(N+M, RHS_RESTO, 1, RW(p_rhs+1), 1)
      else
         call DCOPY(N, G, 1, RW(p_rhscopy+1), 1)
         call DSCAL(N, -1.d0, RW(p_rhscopy+1), 1)
         do i = 1, NLB
            j = ILB(i)
            RW(p_rhscopy+j) = RW(p_rhscopy+j) + V_L(i)
         enddo
         do i = 1, NUB
            j = IUB(i)
            RW(p_rhscopy+j) = RW(p_rhscopy+j) - V_U(i)
         enddo
         if( INIT_LAM ) then
            call DCOPY(M, 0.d0, 0, RW(p_rhscopy+N+1), 1)
            call DCOPY(M+N, RW(p_rhscopy+1), 1, RW(p_rhs+1), 1)
         else
            call DCOPY(M, C, 1, RW(p_rhscopy+N+1), 1)
            call DSCAL(M, -1.d0, RW(p_rhscopy+N+1), 1)
            if( M.gt.0 .and. QLAMBDA.eq.-2 ) then
               inz = N
               if( QQUASI.ne.0 .and. abs(QQUASI).lt.6 ) then
                  inz = inz + (NIND*(NIND+1))/2
               elseif( QQUASI.eq.0 ) then
                  inz = inz + NNZH
               endif
               do i = inz+1, inz+NZA+N_SWCORNER
                  kr  = IS(P_IRN+i)
                  kc  = IS(P_JCN+i) - N
                  ele = RW(p_hcopy+i)
                  RW(p_rhscopy+kr) = RW(p_rhscopy+kr) - ele*LAMOLD(kc)
               enddo
            endif
            if( DEPCONDIAGVAL.ne.0.d0 .or. QCONSTRRELAX.gt.0.d0 ) then
               constrregu = min(-QCONSTRRELAX, DEPCONDIAGVAL)*
     1              MU**QMUDEPEXP
               call DAXPY(M, constrregu, LAMOLD, 1,
     1              RW(p_rhscopy+N+1), 1)
            endif
            do i = 1, NLB
               RW(p_rhscopy+N+M+i) = MU - V_L(i)*S_L(i)
            enddo
            do i = 1, NUB
               RW(p_rhscopy+N+M+NLB+i) = MU - V_U(i)*S_U(i)
            enddo
            call GSF_RHS(N, M, NLB, ILB, S_L, V_L, NUB, IUB, S_U, V_U,
     1           RW(p_rhscopy+1), RW(p_rhs+1))
         endif
      endif
c$$$CDELETEME
c$$$      do i = 1, N+M+NLB+NUB
c$$$         write(QCNR,12345) i, RW(p_rhscopy+i)
c$$$12345    format('RHS(',i6,') = ',d24.16)
c$$$      enddo
C
C     Do the backsolve
C
      refinement_iter = 0
      l_iter_ref = .false.
 1100 continue
CTODO: Is RESTO correct in the following line for L-BFGS? NO WE DO NOT USE LBFGS NOW!!!
C      if( abs(QQUASI).lt.6 .or. RESTO.ne.0 .or. INIT_LAM ) then
      if( abs(QQUASI).lt.6 .or. INIT_LAM ) then
         call GSF_BACKSOLVE(N, M, NZKKT,
#ifdef USE_MALLOC
     1        %VAL(PTR_KKT),
#else
     1        RW(p_kkt+1),
#endif
     1        IS(P_IRN+1), IS(P_JCN+1), RW(p_rhs+1), LIWMA, L_INC_LIWMA,
#ifdef USE_MALLOC
     1        %VAL(PTR_IWMA),
#else
     1        IW(p_iwma+1),
#endif
     1        LA, L_INC_LA, 1, nzero, RW(p_scale+1),
     1        LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     1        LIW-p_iwend, IW(p_iwend+1), IERR)
         if( IERR.ne.0 ) then
            write(line,*)
     1           'get_step_full: GSF_BACKSOLVE returns IERR =',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      else
         call LMBACKSOLVE(N, M, NZKKT,
#ifdef USE_MALLOC
     1        %VAL(PTR_KKT),
#else
     1        RW(p_kkt+1),
#endif
     1        IS(P_IRN+1), IS(P_JCN+1), RW(p_rhs+1), LIWMA, L_INC_LIWMA,
#ifdef USE_MALLOC
     1        %VAL(PTR_IWMA),
#else
     1        IW(p_iwma+1),
#endif
     1        LA, L_INC_LA, 1, nzero,
     1        RW(p_scale+1), nlm, B(pb_uk+1), B(pb_vk+1),
     1        LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     1        LIW-p_iwend, IW(p_iwend+1), IERR)
         if( IERR.ne.0 ) then
            write(line,*)
     1           'get_step_full: LMBACKSOLVE returns IERR =',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      endif
C
C     Improve solution from iterative refinement
C
      if( l_iter_ref ) then
         if( INIT_LAM ) then
            call DAXPY(N+M, 1d0, RW(p_rhs_store+1), 1,
     1           RW(p_rhs+1), 1)
         else
            call DAXPY(N+M+NLB+NUB, 1d0, RW(p_rhs_store+1), 1,
     1           RW(p_rhs+1), 1)
         endif
         p_rwend = p_check
      endif
C
C     Get step for dual variables
C
      if( .not.INIT_LAM .and. RESTO.le.0 ) then
         if( QMPEC_TRIGGER.ge.1 ) then
            do i = 1, NLB
               j = ILB(i)
               RW(p_rhs+N+M+i) =
     1              (RW(p_rhscopy+N+M+i)-V_L(i)*RW(p_rhs+j))/
     2              (S_L(i)+MPEC_ETA*V_L(i))
            enddo
            do i = 1, NUB
               j = IUB(i)
               RW(p_rhs+N+M+NLB+i) =
     1              (RW(p_rhscopy+N+M+NLB+i)+V_U(i)*RW(p_rhs+j))/
     2              (S_U(i)+MPEC_ETA*V_U(i))
            enddo
         else
            do i = 1, NLB
               j = ILB(i)
               RW(p_rhs+N+M+i) =
     1              (RW(p_rhscopy+N+M+i)-V_L(i)*RW(p_rhs+j))/S_L(i)
            enddo
            do i = 1, NUB
               j = IUB(i)
               RW(p_rhs+N+M+NLB+i) =
     1              (RW(p_rhscopy+N+M+NLB+i)+V_U(i)*RW(p_rhs+j))/S_U(i)
            enddo
         endif
      endif
C
C     If in restoration, skip iterative refinement (it is done outside)
C
      if( RESTO.gt.0 ) then
C
C     Compute sum of matrix column entries (note that RHS_RESTO(N+M+1:N+M+N)
C     already has the entries for the regularization term
C
         call DCOPY(N, RHS_RESTO(N+M+1), 1, RHS_RESTO(2*N+M+1), 1)
         call DSCAL(N, ETA_RESTO, RHS_RESTO(2*N+M+1), 1)
         if( REGU.gt.0.d0 ) then
            call DAXPY(N, 1.d0, REGU, 0, RHS_RESTO(2*N+M+1), 1)
         endif
         if( abs(QQUASI).ge.6 ) then
CTODO Good estimate?
            call DAXPY(N, 1.d0, abs(sigmak), 0, RHS_RESTO(2*N+M+1), 1)
         endif
         call DCOPY(M, 0.d0, 0, RHS_RESTO(3*N+M+1), 1)
         do i = N+1, NZKKT-M
            kr  = IS(P_IRN+i)
            kc  = IS(P_JCN+i)
            ele = RW(p_hcopy+i)
            if( kr.eq.kc ) then
               RHS_RESTO(2*N+M+kr) = RHS_RESTO(2*N+M+kr) + abs(ele)
            else
               RHS_RESTO(2*N+M+kr) = RHS_RESTO(2*N+M+kr) + abs(ele)
               RHS_RESTO(2*N+M+kc) = RHS_RESTO(2*N+M+kc) + abs(ele)
            endif
         enddo
         
         goto 2400
      endif
C
C     Check if solution is finite
C
      dummy = DASUM(M+N, RW(p_rhs+1), 1)
      if( FFINITE(dummy).eq.0 ) then
         if( WFLAG.eq.1 ) then
            ALPHA = 0.d0
            goto 3000
         endif
         naninf = .true.
         if( INIT_LAM ) then
            goto 2380
         else
            nrm_ratio = 1.d300
            goto 2360
         endif
      else
         naninf = .false.
      endif
C
C     If this is the first solution of the linear system (no iterative
C     refinement yet) compute norm of right hand side and matrix
C
      if( .not.l_iter_ref ) then
         p_spalte = p_rwend
         p_rwend = p_spalte + M+N+NLB+NUB
         if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
            write(line,*) 'Out of DP memory in line ',__LINE__
            call C_OUT(2,0,1,line)
#endif
            IERR = 98
            goto 9999
         endif
         call DCOPY(M+N+NLB+NUB, 0.d0, 0, RW(p_spalte+1), 1)
         if( INIT_LAM ) then
            call DCOPY(N, 1.d0, 0, RW(p_spalte+1), 1)
         elseif( abs(QQUASI).ge.6 ) then
CTODO Get better estimate!
            call DCOPY(N, max(abs(sigmak),1.d0), 0, RW(p_spalte+1), 1)
         else
            call DCOPY(N, REGU, 0, RW(p_spalte+1), 1)
         endif
         do i = N+1, NZKKT
            kr  = IS(P_IRN+i)
            kc  = IS(P_JCN+i)
            ele = RW(p_hcopy+i)
            if( kr.eq.kc ) then
               RW(p_spalte+kr) = RW(p_spalte+kr) + abs(ele)
            else
               RW(p_spalte+kr) = RW(p_spalte+kr) + abs(ele)
               RW(p_spalte+kc) = RW(p_spalte+kc) + abs(ele)
            endif
         enddo
         if( .not.INIT_LAM ) then
            do i = 1, NLB
               j = ILB(i)
               RW(p_spalte+M+N+i) = RW(p_spalte+M+N+i) + 1.d0
            enddo
            do i = 1, NUB
               j = IUB(i)
               RW(p_spalte+M+N+NLB+i) = RW(p_spalte+M+N+NLB+i) + 1.d0
            enddo
            if( QMPEC_TRIGGER.ge.1 ) then
               do i = 1, NLB
                  j = ILB(i)
                  RW(p_spalte+j) = RW(p_spalte+j) + V_L(i)
                  RW(p_spalte+N+M+i) = RW(p_spalte+N+M+i) +
     1                 S_L(i)+MPEC_ETA*V_L(i)
               enddo
               do i = 1, NUB
                  j = IUB(i)
                  RW(p_spalte+j) = RW(p_spalte+j) + V_U(i)
                  RW(p_spalte+N+M+i) = RW(p_spalte+N+M+i) +
     1                 S_U(i)+MPEC_ETA*V_U(i)
               enddo
            else
               do i = 1, NLB
                  j = ILB(i)
                  RW(p_spalte+j) = RW(p_spalte+j) + V_L(i)
                  RW(p_spalte+N+M+i) = RW(p_spalte+N+M+i) + S_L(i)
               enddo
               do i = 1, NUB
                  j = IUB(i)
                  RW(p_spalte+j) = RW(p_spalte+j) + V_U(i)
                  RW(p_spalte+N+M+i) = RW(p_spalte+N+M+i) + S_U(i)
               enddo
            endif
         endif
         if( INIT_LAM ) then
            i = IDAMAX(N+M, RW(p_rhscopy+1), 1)
            j = IDAMAX(N+M, RW(p_spalte+1) , 1)
         else
            i = IDAMAX(N+M+NLB+NUB, RW(p_rhscopy+1), 1)
            j = IDAMAX(N+M+NLB+NUB, RW(p_spalte+1) , 1)
         endif
         nrm_rhs = abs(RW(p_rhscopy+i))
         nrm_mat = abs(RW(p_spalte+j))
         p_rwend = p_spalte
      endif
C
C     Validate solution (in case of very badly scaled matrix, almost singular)
C
      p_check = p_rwend
      p_rwend = p_check + M+N+NLB+NUB
      if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
         write(line,*) 'Out of DP memory in line ',__LINE__
         call C_OUT(2,0,1,line)
#endif
         IERR = 98
         goto 9999
      endif
C
C     Compute residual
C
      call DCOPY(N, RW(p_rhs+1), 1, RW(p_check+1), 1)
      if( INIT_LAM ) then
         call DSCAL(N, -1.d0, RW(p_check+1), 1)
      elseif( abs(QQUASI).ge.6 ) then
         call DSCAL(N, -sigmak, RW(p_check+1), 1)
      else
         call DSCAL(N, -REGU, RW(p_check+1), 1)
      endif
      call DCOPY(M, 0.d0, 0, RW(p_check+N+1), 1)
      do i = N+1, NZKKT
         kr  = IS(P_IRN+i)
         kc  = IS(P_JCN+i)
         ele = RW(p_hcopy+i)
         if( kr.eq.kc ) then
            RW(p_check+kr) = RW(p_check+kr) - ele*RW(p_rhs+kr)
         else
            RW(p_check+kr) = RW(p_check+kr) - ele*RW(p_rhs+kc)
            RW(p_check+kc) = RW(p_check+kc) - ele*RW(p_rhs+kr)
         endif
      enddo
      if( abs(QQUASI).ge.6 .and. .not.INIT_LAM ) then
CTODO  need to implement estimate of matrix norm
         call LMMULT(N, nlm, B(pb_uk+1), B(pb_vk+1), RW(p_rhs+1),
     1        RW(p_check+1), LRW-p_rwend, RW(p_rwend+1), IERR)
         if( IERR.ne.0 ) then
            write(line,*)
     1           'get_step_full: LMMULT returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      endif
      if( .not.INIT_LAM ) then
         do i = 1, NLB
            j = ILB(i)
            RW(p_check+j) = RW(p_check+j) + RW(p_rhs+M+N+i)
         enddo
         do i = 1, NUB
            j = IUB(i)
            RW(p_check+j) = RW(p_check+j) - RW(p_rhs+M+N+NLB+i)
         enddo
         if( QMPEC_TRIGGER.ge.1 ) then
            do i = 1, NLB
               j = ILB(i)
               RW(p_check+N+M+i) = -V_L(i)*RW(p_rhs+j) -
     1              (S_L(i)+MPEC_ETA*V_L(i))*RW(p_rhs+N+M+i)
            enddo
            do i = 1, NUB
               j = IUB(i)
               RW(p_check+N+M+NLB+i) = V_U(i)*RW(p_rhs+j) -
     1              (S_U(i)+MPEC_ETA*V_U(i))*RW(p_rhs+N+M+NLB+i)
            enddo
         else
            do i = 1, NLB
               j = ILB(i)
               RW(p_check+N+M+i) = -V_L(i)*RW(p_rhs+j) -
     1              S_L(i)*RW(p_rhs+N+M+i)
            enddo
            do i = 1, NUB
               j = IUB(i)
               RW(p_check+N+M+NLB+i) = V_U(i)*RW(p_rhs+j) -
     1              S_U(i)*RW(p_rhs+N+M+NLB+i)
            enddo
         endif
      endif

      if( INIT_LAM ) then
         call DAXPY(N+M, 1d0, RW(p_rhscopy+1), 1, RW(p_check+1), 1)
         i = IDAMAX(N+M, RW(p_check+1)  , 1)
         j = IDAMAX(N+M, RW(p_rhs+1)    , 1)
      else
         call DAXPY(N+M+NLB+NUB, 1d0, RW(p_rhscopy+1), 1,
     1        RW(p_check+1), 1)
         i = IDAMAX(N+M+NLB+NUB, RW(p_check+1)  , 1)
         j = IDAMAX(N+M+NLB+NUB, RW(p_rhs+1)    , 1)
      endif
      nrm_res = abs(RW(p_check+i))
      nrm_sol = abs(RW(p_rhs+j))

C      if( nrm_rhs.gt.0.d0 ) then
C         nrm_ratio = nrm_res/nrm_rhs
CCTODO         nrm_ratio = nrm_res/max(1.d0,nrm_rhs)
C      else
C         nrm_ratio = 0.d0
C      endif

      if( nrm_mat*nrm_sol+nrm_rhs.ne.0.d0 ) then
         nrm_ratio = nrm_res/(nrm_mat*nrm_sol+nrm_rhs)
      else
         nrm_ratio = 0.d0
      endif
      if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
         write(line,1105) nrm_rhs, nrm_sol, nrm_res, nrm_mat, nrm_ratio,
     1        refinement_iter
 1105    format('Linear system solution: |rhs| = ',d8.2,' |sol| = ',
     1        d8.2,' |resid| = ',d8.2,' |mat| = ',d8.2,' ratio = ',
     2        d8.2,' refiter = ',i2)
         call C_OUT(1,0,1,line)
      endif
 2360 continue
cCDELETEME
c      do i = 1, N+M+NLB+NUB
c         write(line,*) i,RW(p_check+i),RW(p_rhs+i),RW(p_rhscopy+i)
c         call C_OUT(1,0,1,line)
c      enddo
c      nrm_ratio_max = min(1.d-10,max(1.d-14,1.d-8 * MU))
      nrm_ratio_max = MACHEPS
      nrm_ratio_min = 1.d-4

      if( refinement_iter.lt.QREFINEITER .or.
     1     nrm_ratio.gt.nrm_ratio_max .or. naninf ) then
         refinement_iter = refinement_iter + 1
C
C     Check iterative refinement really helped
C
         if( (l_iter_ref .and. 0.5d0*nrm_ratio_old.le.nrm_ratio .and.
     1        nrm_ratio.gt.nrm_ratio_max .and.
     1        refinement_iter.ge.QREFINEITER) .or. naninf ) then
            if( .not.naninf ) then
               if( nrm_ratio.lt.nrm_ratio_min .or. SOC_FLAG.ne.0 )
     1              goto 5344
C
C     No, it didn't.  Increase pivot tolerance once and otherwise modify
C     KKT matrix instead
C
               call C_OUT(2,1,1,
     1              'get_step_full: iterative refinement didn''t help.')
            else
               call C_OUT(2,1,1,
     1             'get_step_full: NanInf in solution of linear system')
            endif
            if( QPIVTOL.lt.QPIVTOLMAX .and. .not.L_INCPIVTOL ) then
CORIG            QPIVTOL = min(QPIVTOL*1.d1, QPIVTOLMAX)
               QPIVTOL = min((QPIVTOL**0.75d0), QPIVTOLMAX)
               write(line,*) 'get_step_full: increase QPIVTOL to ',
     1              QPIVTOL
               call C_OUT(2,2,1,line)
#ifdef USE_MALLOC
               call DCOPY(NZKKT, RW(p_hcopy+1), 1, %VAL(PTR_KKT), 1)
#else
               call DCOPY(NZKKT, RW(p_hcopy+1), 1, RW(p_kkt+1), 1)
#endif
               L_INCPIVTOL = .true. ! Make sure we increase QPIVTOL only once
                                ! per iteration
               l_iter_ref = .false.
               p_iwend = p_iwma
               p_rwend = p_rhs
               goto 1500
CTODO            elseif( nrm_ratio.lt.nrm_ratio_min .or. SOC_FLAG.ne.0 ) then
C               goto 5344
            endif
            REFINEMENTHELPLESS = REFINEMENTHELPLESS + 1
            l_iter_ref = .false.
            goto 5000
         else
            l_iter_ref = .true.
         endif
         nrm_ratio_old = nrm_ratio
C
         p_rhs_store = p_rwend
         p_rwend     = p_rhs_store + N+M+NLB+NUB
         if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
            write(line,*) 'Out of DP memory in line ',__LINE__
            call C_OUT(2,0,1,line)
#endif
            write(line,*) 'p_rwend = ',p_rwend,' LRW = ',LRW
            call C_OUT(2,0,1,line)
            IERR = 98
            goto 9999
         endif
         if( INIT_LAM ) then
            call DCOPY(N+M, RW(p_rhs+1), 1, RW(p_rhs_store+1), 1)
            call DCOPY(N+M, RW(p_check+1), 1, RW(p_rhs+1), 1)
         else
            call DCOPY(N+M+NLB+NUB, RW(p_rhs+1), 1,
     1           RW(p_rhs_store+1), 1)
            call GSF_RHS(N, M, NLB, ILB, S_L, V_L, NUB, IUB, S_U, V_U,
     1           RW(p_check+1), RW(p_rhs+1))
         endif
         if( refinement_iter.gt.QREFINEITER .and. QPRINT.ge.4 )
     1     call C_OUT(2,0,1,'get_step_full: Try iterative refinement.')
         goto 1100

      endif
 5344 continue
      p_rwend = p_check
c$$$CDELETEME
c$$$      do i = 1, N+M+NLB+NUB
c$$$         write(QCNR,12346) i, RW(p_rhs+i)
c$$$12346    format('RES(',i6,') = ',d24.16)
c$$$      enddo
C
C     If we just computed least square multipliers, get those and go back
C
 2380 continue
      if( INIT_LAM ) then
         INIT_LAM = .false.
         if( naninf ) then
            call C_OUT(2,2,1,'Least square estimate includes NaN/Inf')
            call C_OUT(2,2,1,' Setting LAMBDA to zero.')
            call DCOPY(M, 0.d0, 0, LAM, 1)
         else
            call C_OUT(1,2,1,
     1           'get_step_full: least square LAMBDA computed!')
            call DCOPY(M, RW(p_rhs+N+1), 1, LAM, 1)
            if( QCNR.gt.0 .and. QPRINT.ge.6 ) then
               call C_OUT(1,0,1,' ')
               call C_OUT(1,0,1,' least square values for LAMBDA:')
               call C_OUT(1,0,1,' ')
               do i = 1, M
                  write(line,2390) i,LAM(i)
 2390             format(' LAM(',i6,') = ',d25.15)
                  call C_OUT(1,0,1,line)
               enddo
               call C_OUT(1,0,1,' ')
            endif
            i = IDAMAX(M, LAM, 1)
            if( abs(LAM(i)).gt.QLAMINITMAX ) then
               write(line,*)
     1              ' Largest entry in LAMBDA is too large:',abs(LAM(i))
               call C_OUT(2,2,1,line)
               call C_OUT(2,2,1,' Setting LAMBDA to zero.')
               call DCOPY(M, 0.d0, 0, LAM, 1)
            endif
         endif
         call DCOPY(M, LAM, 1, LAMOLD, 1) ! need this to make sure filter
                                          ! is not overwriting lambda for SOC
#ifdef USE_MALLOC
         call DCOPY(NZKKT, RW(p_hcopy+1), 1, %VAL(PTR_KKT), 1)
#else
         call DCOPY(NZKKT, RW(p_hcopy+1), 1, RW(p_kkt+1), 1)
#endif
         skip_jac = .true.
         REGU     = 0.d0
         DEPCONDIAGVAL = 0.d0
         p_iwend = p_iwma
         p_rwend = p_kkt
         goto 1000
      endif

 2400 continue
C
C     Store current regularization parameter for next regularization
C
      if( REGU.gt.0.d0 ) then
         REGU_STORE = max( REGU, REGU_MIN )
      endif
C
C     Reset counter for repeated need to increase delta_2
C
      if( DEPCONDIAGVAL.eq.0.d0 ) then
         CHECKDEPENCOUNTER = 0
      endif
C
C     Copy current DX for residual computation in case we are
C     in restoration phase
C
      if( RESTO.eq.2 ) then
         call DCOPY(N, DX, 1, SIGMA_PN, 1)
      endif
C
C     Get primal search direction
C
      call DCOPY(N, RW(p_rhs+1), 1, DX, 1)
C
C     Get multipliers if requested
C
      if( abs(QLAMBDA).eq.2 ) then
         call DCOPY(M, RW(p_rhs+N+1), 1, LAM, 1)
         if( QLAMBDA.eq.-2 .and. RESTO.le.0 ) then
            call DAXPY(M, 1.d0, LAMOLD, 1, LAM, 1)
         endif
      endif
C
C     Get dual search direction
C
      if( RESTO.le.0 ) then
         call DCOPY(NLB, RW(p_rhs+N+M+1), 1, DV_L, 1)
         call DCOPY(NUB, RW(p_rhs+N+M+NLB+1), 1, DV_U, 1)
      endif
C
C     For the KNITRO-DIRECT type of update, we need DX^T (Hessian) DX
C
      if( QKNITROLS.eq.1 .or. QMERIT.eq.-4 ) then
         if( abs(QQUASI).ge.6 ) then
            IERR = 4
            goto 9999
         endif
         p_tmp   = p_rwend
         p_rwend = p_tmp + N
         if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
            write(line,*) 'Out of DP memory in line ',__LINE__,
     1           ' p_rwend = ',p_rwend
            call C_OUT(2,0,1,line)
#endif
            IERR = 98
            goto 9999
         endif
         call DCOPY(N, 0.d0, 0, RW(p_tmp+1), 1)
         do i = 1, N+NNZH
            kr = IS(P_IRN+i)
            kc = IS(P_JCN+i)
            ele = RW(p_hcopy+i)
            if( kr.ne.kc ) then
               RW(p_tmp+kr) =  RW(p_tmp+kr) + ele*DX(kc)
               RW(p_tmp+kc) =  RW(p_tmp+kc) + ele*DX(kr)
            else
               RW(p_tmp+kr) =  RW(p_tmp+kr) + ele*DX(kr)
            endif
         enddo
         if( QQUASI.ne.0 ) then
            IERR = 34343
            goto 9999
         endif
         DHD = DDOT(N, DX, 1, RW(p_tmp+1), 1)
         p_rwend = p_tmp
      endif
C
C     If we are in restoration phase, we need the W * DX product
C     (for computation of residual)
C
      if( RESTO.eq.1 ) then
         if( abs(QQUASI).lt.6 ) then
            call DCOPY(N, 0.d0, 0, RHS_RESTO, 1)
            if( QQUASI.ne.0 ) then
               j = (NIND*(NIND+1))/2
            else
               j = NNZH
            endif
            do i = N+1, N+j
               kr = IS(P_IRN+i)
               kc = IS(P_JCN+i)
               ele = RW(p_hcopy+i)
               if( kr.ne.kc ) then
                  RHS_RESTO(kr) =  RHS_RESTO(kr) + ele*DX(kc)
                  RHS_RESTO(kc) =  RHS_RESTO(kc) + ele*DX(kr)
               else
                  RHS_RESTO(kr) =  RHS_RESTO(kr) + ele*DX(kr)
               endif
            enddo
         else
            call DCOPY(N, DX, 1, RHS_RESTO, 1)
            call DSCAL(N, -sigmak, RHS_RESTO, 1)
            call LMMULT(N, nlm, B(pb_uk+1), B(pb_vk+1), DX,
     1           RHS_RESTO, LRW-p_rwend, RW(p_rwend+1), IERR)
            if( IERR.ne.0 ) then
               write(line,*)
     1              'get_step_full: LMMULT-2 returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            call DSCAL(N, -1.d0, RHS_RESTO, 1)
         endif
         if( REGU.gt.0.d0 )
     1        call DAXPY(N, REGU, DX, 1, RHS_RESTO, 1)
      elseif( RESTO.eq.2 ) then
         call DAXPY(N, 1.d0, DX, 1, SIGMA_PN, 1) !SIGMA_PN is previous DX
         if( abs(QQUASI).lt.6 ) then
            call DCOPY(N, 0.d0, 0, RHS_RESTO, 1)
            if( QQUASI.ne.0 ) then
               j = (NIND*(NIND+1))/2
            else
               j = NNZH
            endif
            do i = N+1, N+j
               kr = IS(P_IRN+i)
               kc = IS(P_JCN+i)
               ele = RW(p_hcopy+i)
               if( kr.ne.kc ) then
                  RHS_RESTO(kr) =  RHS_RESTO(kr) + ele*SIGMA_PN(kc)
                  RHS_RESTO(kc) =  RHS_RESTO(kc) + ele*SIGMA_PN(kr)
               else
                  RHS_RESTO(kr) =  RHS_RESTO(kr) + ele*SIGMA_PN(kr)
               endif
            enddo
         else
            call DCOPY(N, SIGMA_PN, 1, RHS_RESTO, 1)
            call DSCAL(N, -sigmak, RHS_RESTO, 1)
            call LMMULT(N, nlm, B(pb_uk+1), B(pb_vk+1), SIGMA_PN,
     1           RHS_RESTO, LRW-p_rwend, RW(p_rwend+1), IERR)
            if( IERR.ne.0 ) then
               write(line,*)
     1              'get_step_full: LMMULT-3 returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            call DSCAL(N, -1.d0, RHS_RESTO, 1)
         endif
         if( REGU.gt.0.d0 )
     1        call DAXPY(N, REGU, SIGMA_PN, 1, RHS_RESTO, 1)
      endif
C
C     Now can forget everything
C
#ifdef USE_MALLOC
      p_rwend = p_rhscopy
#else
      p_rwend = p_hcopy
#endif
C
C     For iterative refinement in restoration phase, don't do more
C
      if( RESTO.eq.2 ) goto 9999
C
C     Compute maximal step size alpha within bounds
C
      call CUTALPHA(N, X, DX, DV_L, DV_U, NLB, ILB, NUB, IUB,
     1     BNDS_L, BNDS_U, S_L, S_U, V_L, V_U, QTAU,
     2     ALPHA, ALPHA_DUAL, C_ALPHA)
C
      ALPHA_CUT = ALPHA
C
C     Check if need to compute Cauchy step in restoration phase
C
 3000 continue
C
C     In case of Quasi Newton updates, store PZ
C
CTODO CHECK RESTO
CTODO only if QQUASI <> 0 ?
      if( NIND.gt.0 .and.RESTO.le.0 .and. abs(QQUASI).lt.6 ) then
         call DCOPY(NIND, DX(M+1), 1, PZ, 1)
C This is necessary for QORTHO = 1 in get_ypy:
         call DAXPY(NIND, -1.d0, YPY(M+1), 1, PZ, 1)
      endif
C
C     That's it
C
      goto 9999

C
C     Modifying the KKT matrix if inertia is wrong, or iterative refinement
C     didn't help.
C     (If coming from iterative refinement, l_iter_ref is true.
C
 5000 continue
      recomp_regu = .false.
      regu_old = REGU
C
C     If this is about initialization of the multipliers, system
C     seems very ill-conditions, and it is probably better to just
C     the them to zero
C
      if( INIT_LAM ) then
         INIT_LAM = .false.
         call C_OUT(2,1,1,'Least square system singular while '//
     1        'initializing equality multipliers.')
         call C_OUT(2,1,1,'Setting multipliers to zero.')
         call DCOPY(M, 0.d0, 0, LAM, 1)
c         call C_OUT(2,1,1,'Leaving multipliers unchanged.')
         p_iwend = p_iwma
         p_rwend = p_kkt
         goto 1000
      endif
C
C     In first iterations, try to find out if there is structural singularity
C
      if( ITERDEP_MAX.gt.0 .and. nzero.ne.0 .and. RESTO.le.0 ) then
         if( POSSIBLE_DEPCON .and. POSSIBLE_DEPHES ) then
            if( REGU.eq.0.d0 .and. DEPCONDIAGVAL.eq.0.d0 ) then
               DEPCONDIAGVAL = -QDEPCONDIAG
            elseif( DEPCONDIAGVAL.ne.0.d0 .and. REGU.eq.0.d0 ) then
               if( .not.is_dep_this_iter ) DEPCONDIAGVAL = 0.d0
               if( REGU_STORE.eq.0d0 ) then
                  REGU = REGU_INIT
               else
                  REGU = REGU_DEC_FACT*REGU_STORE
               endif
            else
               goto 5010
            endif
         elseif( POSSIBLE_DEPCON ) then
            if( REGU.eq.0.d0 .and. DEPCONDIAGVAL.eq.0.d0 ) then
               DEPCONDIAGVAL = -QDEPCONDIAG
            else
               goto 5010
            endif
         elseif( POSSIBLE_DEPHES ) then
            if( REGU.eq.0.d0 .and. DEPCONDIAGVAL.eq.0.d0 ) then
               if( REGU_STORE.eq.0d0 ) then
                  REGU = REGU_INIT
               else
                  REGU = REGU_DEC_FACT*REGU_STORE
               endif
            else
               goto 5010
            endif
         else
            goto 5010
         endif
         IEIGS = IEIGS + 1
         inz = 0
         if( REGU.ne.regu_old ) then
            if( CORRECTTYPE.eq.1 ) then
               do i = 1, NLB
                  k = ILB(i)
                  RW(p_kkt+inz+k) = min(RW(p_kkt+inz+k),
     1                 REGU/(S_L(i)**2))
               enddo
               do i = 1, NUB
                  k = IUB(i)
                  RW(p_kkt+inz+k) = min(RW(p_kkt+inz+k),
     1                 REGU/(S_U(i)**2))
               enddo
               if( abs(QQUASI).ge.6 ) then
                  call DAXPY(N, 1.d0, sigmak, 0, RW(p_hcopy+inz+1), 1)
               endif
               if( (NIND.gt.0 .and. RESTO.le.0) .or. RESTO.eq.1 ) then
                  do i = 1, NLB
                     k = ILB(i)
                     RW(p_hcopy+inz+k) = RW(p_hcopy+inz+k) + SIGMA_L(i)
                  enddo
                  do i = 1, NUB
                     k = IUB(i)
                     RW(p_hcopy+inz+k) = RW(p_hcopy+inz+k) + SIGMA_U(i)
                  enddo
               endif
            else
               call DAXPY(N, 1.d0, REGU-regu_old, 0,
     1              RW(p_hcopy+inz+1), 1)
            endif
         endif
         inz = N + NZA
         if( QQUASI.eq.0 ) then
            inz = inz + NNZH
         elseif( abs(QQUASI).ne.6 ) then
            inz = inz + (NIND*(NIND+1))/2
         endif
         constrregu = min(-QCONSTRRELAX, DEPCONDIAGVAL)*
     1        MU**QMUDEPEXP
         call DCOPY(M, constrregu, 0, RW(p_hcopy+inz+1), 1)
#ifdef USE_MALLOC
         call DCOPY(NZKKT, RW(p_hcopy+1), 1, %VAL(PTR_KKT), 1)
#else
         call DCOPY(NZKKT, RW(p_hcopy+1), 1, RW(p_kkt+1), 1)
#endif
         goto 1500
      endif
 5010 continue
C
C     First try to regularize by adding
C     a small value for the lower right corner of the KKT matrix
C
      if( QDEPCONDIAG.gt.0.d0 .and. M.gt.0 .and. (RESTO.le.0) ) then
         if( ITERDEP_MAX.gt.0 ) then
            if( DEPCONDIAGVAL.eq.0.d0 .and. nzero.ne.0 ) then
               COUNT_DEPCON = COUNT_DEPCON + 1
               DEPCONDIAGVAL = -QDEPCONDIAG
CTODO            DEPCONDIAGVAL = -QDEPCONDIAG*MU
            endif
         else
            if( DEPCONDIAGVAL.eq.0.d0  ) then
               COUNT_DEPCON = COUNT_DEPCON + 1
               DEPCONDIAGVAL = -QDEPCONDIAG
CTODO            DEPCONDIAGVAL = -QDEPCONDIAG*MU
            endif
         endif
      endif

      if( REGU.eq.0d0 ) then
         COUNT_NEG_CURV = COUNT_NEG_CURV + 1
         if( REGU_STORE.eq.0d0 ) then
            REGU = REGU_INIT
         else
            REGU = REGU_DEC_FACT*REGU_STORE
         endif
      else
         if( REGU_STORE.eq.0d0 .or. 1d5*REGU_STORE.lt.REGU ) then
            REGU = REGU_INIT_FACT*REGU
         else
            REGU = REGU_INC_FACT*REGU
         endif
      endif
      IEIGS = IEIGS + 1
C
      if( REGU.gt.REGU_MAX .or. REFINEMENTHELPLESS.ge.3 ) then
         REFINEMENTHELPLESS = 0
         if( abs(QQUASI).eq.6 ) then
C
C     Reset L-BFGS estimate to identity
C     Need to store current X and G for next update
C
            if( l_lbfgsreset ) then
               call C_OUT(2,0,1, 'get_step_full: REGU getting too '
     1              //'large for L-BFGS after reset')
               IERR = 644
               goto 9999
            endif
            nlm    = 0
            B(1)   = dble(nlm)
            sigmak = 1.d0
            B(2)   = sigmak
            l_lbfgsreset = .true.
            recomp_regu  = .true.
            if( ITERDEP_MAX.gt.0 ) then
               IS_DEPCON = .true.
               call C_OUT(2,1,1,
     1              'get_step_full: REGU becomes too large.  Assume f
     2rom now on that constraint are dependent.')
               DEPCONDIAGVAL = -QDEPCONDIAG
               COUNT_DEPCON = COUNT_DEPCON + 1
               REGU = 0.d0      ! Reset REGU, since problem seems caused by
                                ! constraints
            endif
            goto 5100
         endif
C
C     First try a larger value of DEPCONDIAGVAL
C     AW: I disable this now (2004-03-15) since it doesn't seem to help
C
         if( .false. .and. M.gt.0 .and. QDEPCONDIAG.ne.0.d0 .and.
     1        QDEPCONDIAG.lt.1d3 .and. RESTO.le.0 ) then
            if( WFLAG.eq.1 .or. NIND.eq.0 ) then
C     When in watchdog step, just discard this iteration.  Also, switch
C     to restoration phase if solving square problem
               ALPHA = 0.d0
               goto 3000
            endif
            if( ITERDEP_MAX.gt.0 ) then
               if( IS_DEPCON ) then
                  QDEPCONDIAG   = 1.d2 * QDEPCONDIAG
                  write(line,*)
     1                 'get_step_full: REGU becomes too large,',
     1                 ' try larger QDEPCONDIAG = ',QDEPCONDIAG
                  call C_OUT(2,1,1,line)
               else
                  if( POSSIBLE_DEPCON .or.
c     1                 CHECKDEPENCOUNTER.gt.CHECKDEPENCOUNTERMAX ) then 
     1                 CHECKDEPENCOUNTER.gt.CHECKDEPENCOUNTERMAX .or.
     2                 REGU.gt.REGU_MAX ) then
                     IS_DEPCON = .true.
                     if( REGU.le.REGU_MAX ) then
                        call C_OUT(2,1,1,
     1                       'get_step_full: Too many iterative refineme
     2nt trials.  Assume from now on that constraints are dependent.')
                     else
                        call C_OUT(2,1,1,
     1                       'get_step_full: REGU becomes too large.  As
     2sume from now on that constraints are dependent.')
                     endif
                  else
CTODO
c#define NEW
#ifndef NEW
                     POSSIBLE_DEPCON = .true.
                     CHECKDEPENCOUNTER = CHECKDEPENCOUNTER + 1
                     if( REGU.le.REGU_MAX ) then
                        call C_OUT(2,1,1,
     1                       'get_step_full: Too many iterative refineme
     2nt trials. Check again for dependent constraints.')
                     else
                        call C_OUT(2,1,1,
     1                       'get_step_full: REGU becomes too large.  Ch
     2eck again for dependent constraints.')
                     endif
#else
                     if( is_dep_this_iter ) then
                        ALPHA     = 0d0
                        ALPHA_CUT = 0d0
                        C_ALPHA   = 'R'
                        REGU      = 0d0
                        goto 3000
                     endif
#endif
                     is_dep_this_iter = .true.
                  endif
               endif
               DEPCONDIAGVAL = -QDEPCONDIAG
               COUNT_DEPCON = COUNT_DEPCON + 1
               REGU = 0.d0      ! Reset REGU, since problem seems caused by
                                ! constraints
               recomp_regu = .true.
               goto 5100             
            else
               QDEPCONDIAG   = 1.d2 * QDEPCONDIAG
               DEPCONDIAGVAL = -QDEPCONDIAG
               COUNT_DEPCON = COUNT_DEPCON + 1
               write(line,*) 'get_step_full: REGU becomes too large,',
     1              ' try larger QDEPCONDIAG = ',QDEPCONDIAG
               call C_OUT(2,1,1,line)
               REGU = 0.d0      ! Reset REGU, since problem seems caused by
                                ! constraints
               recomp_regu = .true.
               goto 5100
            endif
         endif
         write(line,*)
     1        'Regularization parameter getting too large (a):', REGU
         call C_OUT(2,1,1,line)
         if( abs(QMERIT).ge.4 .and. ERR_CNSTR.gt.0.d0 .and.
     1        RESTO.le.0 ) then
            ALPHA     = 0d0
            ALPHA_CUT = 0d0
            C_ALPHA   = 'R'
            REGU      = 0d0
            goto 3000
         else
            IERR = 10
            goto 9999
         endif
      endif
C
C     Correct KKT matrix entries
C
 5100 continue
      p_rwend = p_rhscopy + N+M+NLB+NUB
      p_iwend = p_iwma
      inz = 0
CTODO      if( QCORRECT.eq.2 .or. recomp_regu ) then
      if( CORRECTTYPE.eq.1 ) then
         recomp_regu = .true.
      endif
      if( recomp_regu ) then
         if( RESTO.eq.1 ) then
            call DCOPY(N, RHS_RESTO(N+M+1), 1, RW(p_hcopy+inz+1), 1)
            call DSCAL(N, ETA_RESTO, RW(p_hcopy+inz+1), 1)
            call DAXPY(N, 1.d0, REGU, 0, RW(p_hcopy+inz+1), 1)
            if( CORRECTTYPE.eq.1 ) then
               do i = 1, NLB
                  k = ILB(i)
                  RW(p_kkt+inz+k) = min(RW(p_kkt+inz+k),
     1                 REGU/(S_L(i)**2))
               enddo
               do i = 1, NUB
                  k = IUB(i)
                  RW(p_kkt+inz+k) = min(RW(p_kkt+inz+k),
     1                 REGU/(S_U(i)**2))
               enddo
            endif
         else
            call DCOPY(N, REGU, 0, RW(p_hcopy+inz+1), 1)
            if( CORRECTTYPE.eq.1 ) then
               do i = 1, NLB
                  k = ILB(i)
                  RW(p_kkt+inz+k) = min(RW(p_kkt+inz+k),
     1                 REGU/(S_L(i)**2))
               enddo
               do i = 1, NUB
                  k = IUB(i)
                  RW(p_kkt+inz+k) = min(RW(p_kkt+inz+k),
     1                 REGU/(S_U(i)**2))
               enddo
            endif
            if( abs(QQUASI).eq.6 ) then
               call DAXPY(N, 1.d0, sigmak, 0, RW(p_hcopy+inz+1), 1)
            endif
         endif
         if( (NIND.gt.0 .and. RESTO.le.0) .or. RESTO.eq.1 ) then
            do i = 1, NLB
               k = ILB(i)
               RW(p_hcopy+inz+k) = RW(p_hcopy+inz+k) + SIGMA_L(i)
            enddo
            do i = 1, NUB
               k = IUB(i)
               RW(p_hcopy+inz+k) = RW(p_hcopy+inz+k) + SIGMA_U(i)
            enddo
         endif
      else
         call DAXPY(N, 1.d0, REGU-regu_old, 0, RW(p_hcopy+inz+1), 1)
      endif
      inz = inz + N
C
      if( QQUASI.eq.0 ) then
         inz = inz + NNZH
      elseif( abs(QQUASI).ne.6 ) then
         inz = inz + (NIND*(NIND+1))/2
      endif
      inz = inz + NZA
      if( QDEPCONDIAG.ne.0.d0 .and. RESTO.le.0 ) then
         if( ITERDEP_MAX.gt.0 ) then
            constrregu = min(-QCONSTRRELAX, DEPCONDIAGVAL)*
     1           MU**QMUDEPEXP
         else
            constrregu = min(-QCONSTRRELAX, DEPCONDIAGVAL)*
     1           MU**QMUDEPEXP
         endif
         call DCOPY(M, constrregu, 0, RW(p_hcopy+inz+1), 1)
         if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
            write(line,*) 'Lower right corner in matrix = ',constrregu
            call C_OUT(1,0,1,line)
         endif
         inz = inz + M
      endif
#ifdef USE_MALLOC
      call DCOPY(NZKKT, RW(p_hcopy+1), 1, %VAL(PTR_KKT), 1)
#else
      call DCOPY(NZKKT, RW(p_hcopy+1), 1, RW(p_kkt+1), 1)
#endif
C
C     factorize new system
C
      goto 1500


 9999 continue
      call TIMER(timef)
      if( PRINTTIMING .and. QPRINT.ge.2 ) then
         write(line,*) 'get_step_full end time:', timef
         call C_OUT(2,0,1,line)
      endif
      TIME_GET_STEP_FULL = TIME_GET_STEP_FULL + timef - times
      return
      end

      subroutine GSF_BACKSOLVE(N, M, NZKKT, KKT, IRN, JCN, RHS, LIWMA,
     1     L_INC_LIWMA, IWMA, LA, L_INC_LA, NRHS, NZERO, SCALE,
     1     LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR )
C
C     Perform one backsolve with matrix factorize by Harwell routine
C     including scaling
C
      implicit none
      integer N, M, IRN(*), JCN(*), LIWMA, LA, IWMA(*)
      integer NRHS
      integer NZERO, NZKKT
      double precision KKT(*), RHS(N+M,NRHS), SCALE(*)
      logical L_INC_LIWMA, L_INC_LA
      integer LRS
      double precision RS(LRS)
      integer LIS
      integer IS(LIS)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR

      include 'IPOPT.INC'

      integer i, j
      double precision dummy(1)
      character*80 line

      if( QKKTSCALE.ne.0 ) then
         do i = 1, NRHS
            do j = 1, N+M
               RHS(j,i) = RHS(j,i)*SCALE(j)
            enddo
         enddo
      endif
C
C     Solve the system
C
      if( PRINTTIMING .and. QPRINT.ge.2 ) then
         call timer(dummy)
         write(line,*) 'Time before solve = ',dummy
         call C_OUT(2,0,1,line)
      endif
      if( QSYMSOLV.eq.0 ) then
         call MA47_CALL(3, N+M, NZKKT, KKT, IRN, JCN, RHS, LIWMA,
     1        L_INC_LIWMA, IWMA, LA, L_INC_LA, NRHS, NZERO,
     1        LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
      elseif( QSYMSOLV.eq.1 ) then
         call MA27_CALL(3, N+M, NZKKT, KKT, IRN, JCN, RHS, LIWMA,
     1        L_INC_LIWMA, IWMA, LA, L_INC_LA, NRHS, NZERO,
     1        LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
      elseif( QSYMSOLV.eq.2 ) then
         call MA57_CALL(3, N+M, NZKKT, dummy, IRN, JCN, RHS, KKT,
     1        LA, L_INC_LA, IWMA, LIWMA, L_INC_LIWMA, NRHS, NZERO,
     1        LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
      elseif( QSYMSOLV.eq.5 .or. QSYMSOLV.eq.6 ) then
         call MA48SYM_CALL(3, N+M, NZKKT, KKT, IRN, JCN, RHS, LIWMA,
     1        L_INC_LIWMA, IWMA, LA, L_INC_LA, NRHS, NZERO,
     1        LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
      elseif( QSYMSOLV.eq.3 ) then
         call WSSMP_CALL(3, N+M, NZKKT, KKT, IRN, JCN, RHS, LIWMA,
     1        L_INC_LIWMA, IWMA, LA, L_INC_LA, NRHS, NZERO,
     1        LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
      elseif( QSYMSOLV.eq.4 ) then
         call WGSMP_CALL(3, N+M, NZKKT, KKT, IRN, JCN, RHS, LIWMA,
     1        L_INC_LIWMA, IWMA, LA, L_INC_LA, NRHS, NZERO,
     1        LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
      endif
      if( PRINTTIMING .and. QPRINT.ge.2 ) then
         call timer(dummy)
         write(line,*) 'Time after solve = ',dummy
         call C_OUT(2,0,1,line)
      endif
      if( IERR.ne.0 ) then
         write(line,*) 'gsf_backsolve: MA?7_CALL(3) returns IERR =',IERR
         call C_OUT(2,0,1,line)
         if( IERR.gt.99 ) IERR = 507
         goto 9999
      endif
C
C     Undo the scaling
C
      if( QKKTSCALE.ne.0 ) then
         do i = 1, NRHS
            do j = 1, N+M
               RHS(j,i) = RHS(j,i)*SCALE(j)
            enddo
         enddo
      endif

 9999 continue
      return
      end

      subroutine LMBACKSOLVE(N, M, NZKKT, KKT, IRN, JCN, RHS, LIWMA,
     1     L_INC_LIWMA, IWMA, LA, L_INC_LA, NRHS, NZERO, SCALE,
     1     NLM, UK, VK,
     1     LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
C
C     A backsolve for overall KKT system including limited memory approximation
C     using Sherman-Morrison formula
C
      implicit none
      include 'IPOPT.INC'
      integer N, M, IRN(*), JCN(*), LIWMA, LA, IWMA(*), NRHS
      integer NZERO, NZKKT
      double precision KKT(*), RHS(N+M,NRHS), SCALE(*)
      logical L_INC_LIWMA, L_INC_LA
      integer NLM
      double precision UK(N,*), VK(N,*)
      integer LRS
      double precision RS(LRS)
      integer LIS
      integer IS(LIS)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR

      integer i, p_rwend, p_iwend, p_m, p_rhs, p_tmp, p_ipiv, info
      character*80 line

      IERR = 0
      p_rwend = 0
      p_iwend = 0

      if( NRHS.ne.1 ) then
         write(line,*) 'lmbacksolve: NRHS must be one. It is ',NRHS
         call C_OUT(2,0,1,line)
         IERR = 368
         goto 9999
      endif

      if( MEMDBG .and. QPRINT.ge.3 ) then
         write(line,1)'lmbacksolve', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif
C
C     if NLM is still zero, do regular backsolve
C
      if( NLM.eq.0 ) then
         call GSF_BACKSOLVE(N, M, NZKKT, KKT, IRN, JCN, RHS,
     1        LIWMA, L_INC_LIWMA, IWMA, LA, L_INC_LA, NRHS, NZERO,
     1        SCALE, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     5        LIW-p_iwend, IW(p_iwend+1), IERR)
         goto 9999
      endif

      p_m     = p_rwend
      p_rhs   = p_m + 4*NLM*NLM
      p_rwend = p_rhs + NLM*(N+M)
      if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
         write(line,*) 'Out of DP memory in line ',__LINE__
         call C_OUT(2,0,1,line)
#endif
         IERR = 98
         goto 9999
      endif
C
C     First compute the dense matrix in the middle of the Sherman-Morrison
C     formula
C
      call DCOPY(4*NLM*NLM, 0.d0, 0, RW(p_m+1), 1)
      do i = 1, 2*NLM
         RW(p_m+i+2*NLM*(i-1)) = 1.d0
      enddo
C
C     Backsolves for the UK part
C
      do i = 1, NLM
         call DCOPY(N, UK(1,i), 1, RW(p_rhs+1+(N+M)*(i-1)), 1)
         call DCOPY(M, 0.d0, 0, RW(p_rhs+N+1+(N+M)*(i-1)), 1)
      enddo
      call GSF_BACKSOLVE(N, M, NZKKT, KKT, IRN, JCN, RW(p_rhs+1),
     1     LIWMA, L_INC_LIWMA, IWMA, LA, L_INC_LA, NLM, NZERO,
     1     SCALE, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR)
      if( IERR.ne.0 ) then
         write(line,*)
     1        'lmbacksolve: GSF_BACKSOLVE-1 returns IERR =',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      do i = 1, NLM
         call DGEMV('T', N, NLM, 1.d0, UK, N, RW(p_rhs+1+(N+M)*(i-1)),
     1        1, 1.d0, RW(p_m+1+2*NLM*(i-1)), 1)
         call DGEMV('T', N, NLM, -1.d0, VK, N, RW(p_rhs+1+(N+M)*(i-1)),
     1        1, 1.d0, RW(p_m+NLM+1+2*NLM*(i-1)), 1)
      enddo
C
C     Backsolves for the VK part
C
      do i = 1, NLM
         call DCOPY(N, VK(1,i), 1, RW(p_rhs+1+(N+M)*(i-1)), 1)
         call DCOPY(M, 0.d0, 0, RW(p_rhs+N+1+(N+M)*(i-1)), 1)
      enddo
      call GSF_BACKSOLVE(N, M, NZKKT, KKT, IRN, JCN, RW(p_rhs+1),
     1     LIWMA, L_INC_LIWMA, IWMA, LA, L_INC_LA, NLM, NZERO,
     1     SCALE, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR)
      if( IERR.ne.0 ) then
         write(line,*)
     1        'lmbacksolve: GSF_BACKSOLVE-2 returns IERR =',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      do i = 1, NLM
         call DGEMV('T', N, NLM, 1.d0, UK, N, RW(p_rhs+1+(N+M)*(i-1)),
     1        1, 1.d0, RW(p_m+1+2*NLM*(NLM+i-1)), 1)
         call DGEMV('T', N, NLM, -1.d0, VK, N, RW(p_rhs+1+(N+M)*(i-1)),
     1        1, 1.d0, RW(p_m+NLM+1+2*NLM*(NLM+i-1)), 1)
      enddo
C
C     Compute K^{-1} * RHS -> RHS   (K: KKT matrix with only diagonal)
C
      call GSF_BACKSOLVE(N, M, NZKKT, KKT, IRN, JCN, RHS,
     1     LIWMA, L_INC_LIWMA, IWMA, LA, L_INC_LA, NRHS, NZERO,
     1     SCALE, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR)
      if( IERR.ne.0 ) then
         write(line,*)
     1        'lmbacksolve: GSF_BACKSOLVE-3 returns IERR =',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     Compute tmp = K^{-1} * UU * M^{-1} * VV^T RHS
C
      p_tmp   = p_rwend
      p_rwend = p_tmp + 2*NLM
      if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
         write(line,*) 'Out of DP memory in line ',__LINE__
         call C_OUT(2,0,1,line)
#endif
         IERR = 98
         goto 9999
      endif
      call DGEMV('T', N, NLM, 1.d0, UK, N, RHS, 1, 0.d0, RW(p_tmp+1), 1)
      call DGEMV('T', N, NLM, -1.d0, VK, N, RHS, 1, 0.d0,
     1     RW(p_tmp+NLM+1), 1)
      p_ipiv  = p_iwend
      p_iwend = p_ipiv + 2*NLM
      if( p_iwend.gt.LIW ) then
#ifdef PRINT__LINE__
         write(line,*) 'Out of INT memory in line ',__LINE__
         call C_OUT(2,0,1,line)
#endif
         IERR = 99
         goto 9999
      endif
      call DGESV( 2*NLM, 1, RW(p_m+1), 2*NLM, IW(p_ipiv+1), RW(p_tmp+1),
     1     2*NLM, info )
      if( info.ne.0 ) then
         write(line,*) 'lmbacksolve: DBESV returns info = ',info
         call C_OUT(2,0,1,line)
         IERR = 823
         goto 9999
      endif
      p_iwend = p_ipiv
      call DGEMV('N', N, NLM, 1.d0, UK, N, RW(p_tmp+1), 1, 0.d0,
     1     RW(p_rhs+1), 1)
      call DGEMV('N', N, NLM, 1.d0, VK, N, RW(p_tmp+NLM+1), 1, 1.d0,
     1     RW(p_rhs+1), 1)
      p_rwend = p_tmp
      call DCOPY(M, 0.d0, 0, RW(p_rhs+N+1), 1)
      call GSF_BACKSOLVE(N, M, NZKKT, KKT, IRN, JCN, RW(p_rhs+1),
     1     LIWMA, L_INC_LIWMA, IWMA, LA, L_INC_LA, NRHS, NZERO,
     1     SCALE, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR)
      if( IERR.ne.0 ) then
         write(line,*)
     1        'lmbacksolve: GSF_BACKSOLVE-4 returns IERR =',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     Finally put everything together
C
      call DAXPY(N+M, -1.d0, RW(p_rhs+1), 1, RHS, 1)

 9999 continue
      return
      end

      subroutine LMMULT(N, NLM, UK, VK, VIN, VOUT, LRW, RW, IERR)
C
C     Compute VOUT = VOUT - (UK * UK^T - VK * VK^T)*VIN
C     (note the MINUS!)   ^
C
      implicit none
      integer N, NLM
      double precision UK(N, NLM)
      double precision VK(N, NLM)
      double precision VIN(N), VOUT(N)
      integer LRW
      double precision RW(LRW)
      integer IERR

      integer p_rwend, p_tmp

      IERR = 0
      p_rwend = 0

      if( NLM.eq.0 ) goto 9999  ! nothing to do

      p_tmp   = p_rwend
      p_rwend = p_tmp + NLM
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif

      call DGEMV( 'T', N, NLM, 1.d0, UK, N, VIN, 1, 0.d0,
     1     RW(p_tmp+1), 1)
      call DGEMV( 'N', N, NLM, -1.d0, UK, N, RW(p_tmp+1), 1, 1.d0,
     1     VOUT, 1)

      call DGEMV( 'T', N, NLM, 1.d0, VK, N, VIN, 1, 0.d0,
     1     RW(p_tmp+1), 1)
      call DGEMV( 'N', N, NLM, 1.d0, VK, N, RW(p_tmp+1), 1, 1.d0,
     1     VOUT, 1)

 9999 continue
      return
      end

      subroutine GSF_SCALE(NTOT, NZKKT, IRN, JCN, SCALE, KKT)
      implicit none
      integer NTOT, NZKKT, IRN(NZKKT), JCN(NZKKT)
      double precision SCALE(NTOT), KKT(NZKKT)

      integer i, kr, kc

      do i = 1, NZKKT
         kr = IRN(i)
         kc = JCN(i)
         KKT(i) = KKT(i)*SCALE(kr)*SCALE(kc)
      enddo
      return
      end

      subroutine GSF_RHS(N, M, NLB, ILB, S_L, V_L, NUB, IUB, S_U, V_U,
     1           RHSCOPY, RHS)
      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, ILB(NLB), NUB, IUB(NUB)
      double precision S_L(NLB), V_L(NLB), S_U(NUB), V_U(NUB)
      double precision RHSCOPY(N+M+NLB+NUB), RHS(N+M)

      double precision MPEC_ETA
      common /MPEC/ MPEC_ETA
      save /MPEC/

      integer i, j

      call DCOPY(N+M, RHSCOPY, 1, RHS, 1)

      if( QMPEC_TRIGGER.ge.1 ) then
         do i = 1, NLB
            j = ILB(i)
            RHS(j) = RHS(j) +
     1           RHSCOPY(N+M+i)/(S_L(i)+MPEC_ETA*V_L(i))
         enddo
         do i = 1, NUB
            j = IUB(i)
            RHS(j) = RHS(j) -
     1           RHSCOPY(N+M+NLB+i)/(S_U(i)+MPEC_ETA*V_U(i))
         enddo
      else
         do i = 1, NLB
            j = ILB(i)
            RHS(j) = RHS(j) + RHSCOPY(N+M+i)/S_L(i)
         enddo
         do i = 1, NUB
            j = IUB(i)
            RHS(j) = RHS(j) -  RHSCOPY(N+M+NLB+i)/S_U(i)
         enddo
      endif

      return
      end

#ifdef USE_MALLOC
C ==============================================================================
C
C     Free pointers (at end of optimization)
C
C ==============================================================================

      subroutine GSF_FREEPTR(IERR)
      implicit none
      integer IERR
      integer i, type, IP_MALLOC
      character*80 line

#ifdef SIZEOF_INT_P_IS_4
      integer*4 PTR_KKT, PTR_IWMA
#elif SIZEOF_INT_P_IS_8
      integer*8 PTR_KKT, PTR_IWMA
#else
      Size of integer pointers unknown
#endif
      common /IP_ALLOC/ PTR_KKT, PTR_IWMA
      save /IP_ALLOC/
      IERR = 0
      type = 1
      i = IP_MALLOC(type, 0, PTR_KKT)
      if( i.ne.0 ) then
         write(line,*) 'gsf_freeptr: IP_MALLOC returned ',
     1        i,' for KKT'
         call C_OUT(2,0,1,line)
         IERR = 96
         goto 9999
      endif
      type = 2
      i = IP_MALLOC(type, 0, PTR_IWMA)
      if( i.ne.0 ) then
         write(line,*) 'gsf_freeptr: IP_MALLOC returned ',
     1        i,' for IWMA'
         call C_OUT(2,0,1,line)
         IERR = 96
         goto 9999
      endif
 9999 continue
      return
      end
#endif

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine GET_STEP_FULL_WS(N, X, M, NLB, NUB, NZA, LRW, LIW,
     1     FILLINFACT, PRWKEEP, PIWKEEP, EV_H, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, nzh, LRW, LIW, PRWKEEP, PIWKEEP
      external EV_H
      double precision DAT(*)
      integer IDAT(*)
      double precision FILLINFACT, X(N)
      integer nzkkt, la, liwma, lrw1, liw1, lrw2, liw2, lrw3, liw3
      integer idummy(1)
      double precision dummy(1)
      character*80 line

      LIW = 0
      LRW = 0

C     get number of nonzeros in Hessian
      if( QQUASI.eq.0 ) then
         call GET_H(-1, N, N-M, idummy, dummy, idummy, N,
     1        X, dummy, M, nzh, dummy, dummy, idummy, idummy,
     1        idummy, dummy, idummy, idummy, idummy, EV_H, DAT, IDAT)
      else
         nzh = 0
      endif

      nzkkt = N + NZA     ! worst estimate
      if( QQUASI.ne.0 .and. abs(QQUASI).lt.6 ) then
         nzkkt = nzkkt + ((N-M)*((N-M)+1))/2
      elseif( QQUASI.eq.0 ) then
         nzkkt = nzkkt + nzh
      endif
      if( QDEPCONDIAG.gt.0.d0 .or. QSYMSOLV.eq.3 .or. QSYMSOLV.eq.4
CORIG     1     .or. (abs(QMERIT).ge.4.and.QRESTO.eq.1)) then
     1     .or. (abs(QMERIT).ge.4.and.QRESTO.eq.1) .or.
     1     QKNITROLS.ne.0 ) then
         nzkkt = nzkkt + M
      endif
      if( QSYMSOLV.eq.0 ) then
         call MA47_CALL_WS(0, N, M, NLB, NUB, NZA, nzkkt, la,
     1        liwma, FILLINFACT)
      elseif( QSYMSOLV.eq.1 ) then
         call MA27_CALL_WS(0, N, M, NLB, NUB, NZA, nzkkt, la,
     1        liwma, FILLINFACT)
      elseif( QSYMSOLV.eq.2 ) then
         call MA57_CALL_WS(0, N, M, NLB, NUB, NZA, nzkkt, la,
     1        liwma, FILLINFACT)
      elseif( QSYMSOLV.eq.5 .or. QSYMSOLV.eq.6 ) then
         call MA48SYM_CALL_WS(0, N, M, NLB, NUB, NZA, nzkkt, la,
     1        liwma, FILLINFACT)
      elseif( QSYMSOLV.eq.3 ) then
         call WSSMP_CALL_WS(0, N, M, NLB, NUB, NZA, nzkkt, la,
     1        liwma, FILLINFACT)
      elseif( QSYMSOLV.eq.4 ) then
         call WGSMP_CALL_WS(0, N, M, NLB, NUB, NZA, nzkkt, la,
     1        liwma, FILLINFACT)
      endif

      lrw1 = 0
      liw1 = 0
      if( QQUASI.eq.0 ) then
         call GET_H_WS(N, M, NLB, NUB, NZA, nzh, lrw1, liw1)
         if( QRESTO.eq.1 .and. abs(QMERIT).ge.4 ) then
            lrw1 = lrw1 + NZH + M
         endif
      endif
      if( M.gt.0 ) then
         call CONSTR_WS(N, M, NLB, NUB, NZA, lrw2, liw2, DAT, IDAT)
         lrw1 = max(lrw1,lrw2)
         liw1 = max(liw1,liw2+NZA)
      endif
      if( QSYMSOLV.eq.0 ) then
         call MA47_CALL_WS(1, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      elseif( QSYMSOLV.eq.1 ) then
         call MA27_CALL_WS(1, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      elseif( QSYMSOLV.eq.2 ) then
         call MA57_CALL_WS(1, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      elseif( QSYMSOLV.eq.5 .or. QSYMSOLV.eq.6 ) then
         call MA48SYM_CALL_WS(1, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      elseif( QSYMSOLV.eq.3 ) then
         call WSSMP_CALL_WS(1, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      elseif( QSYMSOLV.eq.4 ) then
         call WGSMP_CALL_WS(1, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      endif
      lrw1 = max(lrw1,lrw2)
      liw1 = max(liw1,liw2)

#ifdef USE_MALLOC
      LIW = liw1
      LRW = nzkkt + lrw1
      PIWKEEP = 0
      PRWKEEP = nzkkt
#else
      LIW = liwma + liw1
      LRW = nzkkt + lrw1
      PIWKEEP = liwma
      PRWKEEP = la
      PRWKEEP = PRWKEEP + nzkkt
#endif
      if( QKKTSCALE.ne.0 ) then
         LRW = LRW + N+M
         PRWKEEP = PRWKEEP + N+M
      endif

      lrw1 = 0
      liw1 = 0
      if( QKKTSCALE.eq.1 ) then
#ifdef HAVE_MC30
         lrw1 = 4*(N+M)
#else
#ifdef HAVE_MC29
         call MC29_CALL_WS(N, M, NLB, NUB, NZA, nzkkt, lrw1, liw1)
#else
#ifdef HAVE_MC19
         call MC19_CALL_WS(N, M, NLB, NUB, NZA, nzkkt, lrw1, liw1)
         lrw1 = lrw1 + (7*(N+M))/2 + 1
#endif
#endif
#endif
      endif
      if( QSYMSOLV.eq.0 ) then
         call MA47_CALL_WS(2, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      elseif( QSYMSOLV.eq.1 ) then
         call MA27_CALL_WS(2, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      elseif( QSYMSOLV.eq.2 ) then
         call MA57_CALL_WS(2, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
         lrw2 = lrw2 + nzkkt
      elseif( QSYMSOLV.eq.5 .or. QSYMSOLV.eq.6 ) then
         call MA48SYM_CALL_WS(2, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      elseif( QSYMSOLV.eq.3 ) then
         call WSSMP_CALL_WS(2, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      elseif( QSYMSOLV.eq.4 ) then
         call WGSMP_CALL_WS(2, N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1        liw2, FILLINFACT)
      endif
      lrw1 = max(lrw1, lrw2)
      liw1 = max(liw1, liw2)

      call GSF_BACKSOLVE_WS(N, M, NLB, NUB, NZA, nzkkt, lrw2,
     1     liw2, FILLINFACT)
      if( abs(QQUASI).eq.6 ) then
         call LMBACKSOLVE_WS(N, M, NLB, NUB, NZA, nzkkt, lrw3,
     1        liw3, FILLINFACT)
         lrw2 = max(lrw2, lrw3)
         liw2 = max(liw2, liw3)
         call LMMULT_WS(N, M, NLB, NUB, NZA, lrw3, liw3)
         lrw2 = max(lrw2, lrw3)
         liw2 = max(liw2, liw3)
      endif
      lrw1 = max(lrw1, 3*(M+N+NLB+NUB)+lrw2)
      liw1 = max(liw1, liw2)

      lrw1 = lrw1 + nzkkt + (N+M+NLB+NUB)
      if( QKKTSCALE.ne.0 ) then
         lrw1 = lrw1 + (N+M)
      endif

#ifdef USE_MALLOC
      LRW = max(LRW, lrw1)
      LIW = max(LIW, liw1)
#else
      LRW = max(LRW, la + lrw1)
      LIW = max(LIW, liwma+liw1)
#endif

      if( QPRINT.ge.3 ) then
         write(line,1000)'get_step_full_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end


      subroutine GSF_BACKSOLVE_WS(N, M, NLB, NUB, NZA, NZKKT, LRW,
     1     LIW, FILLINFACT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, NZKKT, LRW, LIW
      double precision FILLINFACT
      character*80 line

      if( QSYMSOLV.eq.0 ) then
         call MA47_CALL_WS(3, N, M, NLB, NUB, NZA, NZKKT, LRW, LIW,
     1        FILLINFACT)
      elseif( QSYMSOLV.eq.1 ) then
         call MA27_CALL_WS(3, N, M, NLB, NUB, NZA, NZKKT, LRW, LIW,
     1        FILLINFACT)
      elseif( QSYMSOLV.eq.2 ) then
         call MA57_CALL_WS(3, N, M, NLB, NUB, NZA, NZKKT, LRW, LIW,
     1        FILLINFACT)
      elseif( QSYMSOLV.eq.5 .or. QSYMSOLV.eq.6 ) then
         call MA48SYM_CALL_WS(3, N, M, NLB, NUB, NZA, NZKKT, LRW,
     1        LIW, FILLINFACT)
      elseif( QSYMSOLV.eq.3 ) then
         call WSSMP_CALL_WS(3, N, M, NLB, NUB, NZA, NZKKT, LRW, LIW,
     1        FILLINFACT)
      elseif( QSYMSOLV.eq.4 ) then
         call WGSMP_CALL_WS(3, N, M, NLB, NUB, NZA, NZKKT, LRW, LIW,
     1        FILLINFACT)
      endif

      if( QPRINT.ge.3 ) then
         write(line,1000)'gsf_backsolve_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end

      subroutine LMBACKSOLVE_WS(N, M, NLB, NUB, NZA, NZKKT, LRW,
     1     LIW, FILLINFACT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, NZKKT, LRW, LIW
      double precision FILLINFACT
      integer lrw_bs, liw_bs
      character*80 line

      call GSF_BACKSOLVE_WS(N, M, NLB, NUB, NZA, NZKKT, lrw_bs,
     1     liw_bs, FILLINFACT)

      LRW = QLMLEN*(N+M)+4*QLMLEN*QLMLEN+max(lrw_bs,2*QLMLEN)
      LIW = max(liw_bs, 2*QLMLEN)

      if( QPRINT.ge.3 ) then
         write(line,1000)'lmbacksolve_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end

      subroutine LMMULT_WS(N, M, NLB, NUB, NZA, LRW, LIW)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      character*80 line

      LIW = 0
      LRW = QLMLEN

      if( QPRINT.ge.3 ) then
         write(line,1000)'lmmult_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
