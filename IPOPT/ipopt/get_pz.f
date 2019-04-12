C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_PZ(ITER, MU, NIND, W, RGB, WCORR, CNRM, B, PZ,
     1     IEIGS, N, M, X, IVAR, NFIX, IFIX, NORIG, XORIG, CSCALE, LAM,
     2     NLB, ILB, NUB, IUB, S_L, S_U, BNDS_L, BNDS_U,
     3     SIGMA_L, SIGMA_U, YPY, REGU, SOC_FLAG, RESTO,
     3     KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR,
     2     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     7     DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: get_pz.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Solve for null space step pZ
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
CP   ITER      I    INT    iteration counter
CP   MU        I    DP     barrier parameter
CP   NIND      I    INT    number of independent variables
CP   W         I    DP     reduced Hessian
CP   RGB       I    DP     reduced gradient of barrier function
CP   WCORR     I    DP     correction term
CP   CNRM      I    DP     infty(???) norm of constraints
CP   B         I    DP     Quasi-Newton matrix
CP   PZ        O    DP     null space step pY (only for independent variables)
CP   IEIGS     O    INT    (only for |QQUASI|=0,2): Number of corrected
CP                             eigenvalues
CP   N         I    INT    number of (free) variables; first M variables
CP                         are dependent; remaining independent
CP   M         I    INT    number of dependent variables
CP   X         I    DP     actual primal iterate
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
CP   S_L       I    DP     slack variables for lower bounds
CP   S_U       I    DP     slack variables for upper bounds
CP   BNDS_L    I    DP     values of lower bounds (ordered as S_L)
CP   BNDS_U    I    DP     values of upper bounds (ordered as S_U)
CP   SIGMA_L   I    DP     primal-dual Hessian of lower bound barrier term
CP                            (NLB diagonal elements only)
CP   SIGMA_U   I    DP     primal-dual Hessian of upper bound barrier term
CP                            (NUB diagonal elements only)
CP   YPY       I    DP     range space step (ordered like X)
CP   REGU      O    DP     regularization factor (added regu*I to diagonal)
CP   SOC_FLAG  I    INT    =0: not in SOC computation
CP                          1: in SOC computation (don't need to refactorize)
CP   RESTO     I    INT    <>0: we are in restoration phase
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
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW       I    INT    length of IW
CP   IW       I/O   INT    can be used as INT work space but content will be
CP                            changed between calls
CP   IERR      O    INT    INFO from DDSV and so on
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
CS    IDAMAX
CS    DDOT
CS    DSCAL
CS    DCOPY
CS    DAXPY
CS    DSPEV
CS    DSPMV
CS    DGESV
CS    DGEMV
CS    DTPSV
CS    C_OUT
CS    CONSTR
CS    GET_WCORR
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
C    |QQUASI|= 1: BFGS is used for B, simple solve symmetric system
C              2: SR1 is used for B, so need to check if W has to be corrected
C                 so that it is positive definite; for this, eigenvalues
C                 as computed...
C              3: First use BFGS, and if error below QSR1TOL, switch to SR1
C              4: As 2, but correct B (the estimate matrix)
C              5: like 1, but used Powell damping in BFGS update
C
C    |QCORRECT| = 1: Do modified Cholesky factorization
C                 2: Do eigenvalue decomposition and make eigenvalues
C                    sufficiently positive
C                 3: Do Cholesky factorization; if indefinite, add multiple
C                    of identity to B (or equivalently to Sigma_ind)
C                 4: Do Cholesky factorization; if indefinite, add multiple
C                    of identity to OVERALL Sigma (like in full-space version)
C      QCORRECT < 0: Solve for PZ frist using LU decomposition; if direction of
C                    negative curvature, do corrections as outlined above.
C
      include 'IPOPT.INC'
      include 'TIMER.INC'
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer ITER
      double precision MU
      integer NIND
      double precision W((NIND*(NIND+1))/2)
      double precision RGB(NIND)
      double precision WCORR(NIND)
      double precision CNRM
      double precision B((NIND*(NIND+1))/2)
      double precision PZ(NIND)
      integer IEIGS
      integer N
      integer M
      double precision X(N)
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
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision BNDS_L(NLB)
      double precision BNDS_U(NUB)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision YPY(N)
      double precision REGU
      integer SOC_FLAG
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
      double precision DDOT
      integer IDAMAX

      integer p_iwend, p_rwend, p_wtmp, p_eigs, p_z, p_work, p_tmp
      integer p_corr, p_icorr, p_a, p_ipiv, p_w, p_iperm
      integer p_ztz, idummy, p_sigma, p_bpz, p_bzg
      integer info, i, k, j
      double precision lam_b, mineig
      double precision eigmin, eigmax, eigold, val1, val2, regu1, diff
      double precision zeta, gzbw, gzbzg
      character*80 line(4)
C
      logical l_regu, indef_bfgs

C     if you want to use modified Cholseky for ICORRECT = 2, set lmodchol to
C     .true. (from some test on HS set it is not clear which version is better,
C     so I left in both implementations)
      logical lmodchol
      parameter( lmodchol = .false. )

      double precision REGU_STORE
      save             REGU_STORE
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
C
C     There is nothing to do if there are no degrees of freedom
C
      if( NIND.eq.0 ) goto 9999

      if( ITER.eq.0 ) then
         REGU_STORE = 0d0
      endif
      indef_bfgs = .false.

      if( MEMDBG ) then
         write(line,1)'get_pz', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

 100  continue
      p_rwend = 0
      p_iwend = 0
      IEIGS   = 0
C
C     Prepare RHS:  PZ = - ( RGB + WCORR )
C
      if( RESTO.eq.0 ) then
         call DCOPY(NIND, RGB, 1, PZ, 1)
         call DAXPY(NIND, 1d0, WCORR, 1, PZ, 1)
      else
         call DCOPY(NIND, WCORR, 1, PZ, 1)
      endif
      call DSCAL(NIND, -1.d0, PZ, 1)

      if( indef_bfgs ) goto 2000

      goto (1000, 2000, 1000, 2000, 1000) abs(QQUASI)
      if( QQUASI.eq.0 ) goto 2000

      write(line,*) 'Invalid Choice of QQUASI = ', QQUASI
      call C_OUT(2,0,1,line)
      IERR = 4
      goto 9999

C ==============================================================================
C
C     BEGIN : BFGS case
C
C ==============================================================================

 1000 continue
C
C     Solve the system to obtain PZ with (modified) Cholesky factorization
C
      p_wtmp  = p_rwend
      p_rwend = p_wtmp + (NIND*(NIND+1))/2
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call DCOPY((NIND*(NIND+1))/2, W, 1, RW(p_wtmp+1), 1)

      p_iperm = p_iwend
      p_iwend = p_iperm + NIND
      p_w     = p_rwend
      p_rwend = p_w + NIND
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      if( QDAMP.eq.0 .or. RESTO.ne.0 ) then
         if( QCORRECT.ne.1 .and. .not.indef_bfgs ) then
            IEIGS = -1
         else
            IEIGS = 0
         endif
         call MOD_CHOL(NIND, RW(p_wtmp+1), IW(p_iperm+1), PZ,
     1        IEIGS, RW(p_w+1))
         if( IEIGS.ne.0 .and. QCORRECT.ne.1 .and. .not.indef_bfgs ) then
C
C     In this case, W is (almost) singular; go to SR1 corrections after
C     recomputing RHS
C
            indef_bfgs = .true.
            goto 100
         endif
      else
         p_bzg   = p_rwend
         p_rwend = p_bzg + NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         call DCOPY(NIND, RGB, 1, RW(p_bzg+1), 1)
         if( QCORRECT.ne.1 .and. .not.indef_bfgs ) then
            IEIGS = -1
         else
            IEIGS = 0
         endif
         call MOD_CHOL(NIND, RW(p_wtmp+1), IW(p_iperm+1),
     1        RW(p_bzg+1), IEIGS, RW(p_w+1))
         if( IEIGS.ne.0 .and. QCORRECT.ne.1 .and. .not.indef_bfgs ) then
C
C     In this case, W is (almost) singular; go to SR1 corrections after
C     recomputing RHS
C
            indef_bfgs = .true.
            goto 100
         endif
         gzbw = DDOT(NIND, WCORR, 1, RW(p_bzg+1), 1)
         if( gzbw.lt.0.d0 ) then
            gzbzg = DDOT(NIND, RGB, 1, RW(p_bzg+1), 1)
            zeta = min(1.d0, -0.1d0*gzbzg/gzbw)
            if( zeta.lt.1.d0 ) then
               write(line,*) 'Damping active: ZETA = ',zeta
               call C_OUT(2,0,1,line)
               call DAXPY(NIND, 1.d0-zeta, WCORR, 1, PZ, 1)
            endif
         endif
         call MOD_CHOL_SOL(NIND, RW(p_wtmp+1), IW(p_iperm+1), PZ,
     1        RW(p_w+1))
      endif
      p_rwend = p_wtmp
      p_iwend = p_iperm
      goto 9999

C ==============================================================================
C
C     END   : BFGS case
C
C ==============================================================================

C ==============================================================================
C
C     BEGIN : SR1 case
C
C ==============================================================================

 2000 continue
C
C     Choose way to deal with indefiniteness
C
      if( QCORRECT.lt.0 ) goto 2300 ! Try LU factorization first

 2050 continue

      goto( 1000, 2100, 2200, 2200 ) abs(QCORRECT)

      write(line,*) 'Invalid QCORRECT = ',QCORRECT
      call C_OUT(2,0,1,line)
      IERR = 4
      goto 9999

 2100 continue
C
C     Eigenvalue correction
C     =====================
C
C     Compute eigenvalues and eigenvectors
C
      p_eigs  = p_rwend
      p_z     = p_eigs + NIND
      p_wtmp  = p_z    + NIND*NIND
      p_work  = p_wtmp + (NIND*(NIND+1))/2
      p_rwend = p_work + 3*NIND
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call DCOPY((NIND*(NIND+1))/2, W, 1, RW(p_wtmp+1), 1)
C
C     After this call, eigs contains the eigenvalues and z contains the
C     eigenvectors
C
      call DSPEV('V', 'U', NIND, RW(p_wtmp+1), RW(p_eigs+1),
     1     RW(p_z+1), NIND, RW(p_work+1), info)
      if( info.ne.0 ) then
         write(line,*) 'get_pz: DSPEV returned INFO = ', info
         call C_OUT(2,0,1,line)
         IERR = 501
         goto 9999
      endif
      p_rwend = p_wtmp
      if( QCNR.gt.0 .and. QPRINT.ge.4 ) then
         write(line,2110) ITER
 2110    format('  Information about eigenvals of W in ITER', i5)
         call C_OUT(1,4,0,line)
         call C_OUT(1,4,1,line)
         call C_OUT(1,4,0,line)
         do i = 1, NIND
            write(line,2111) i, RW(p_eigs+i)
 2111       format(' eig(',i3,') = ',d16.10)
            call C_OUT(1,4,1,line)
         enddo
      endif
C
C     In case, B will be changed according to eigenvalue correction, reserve
C     space for storing corrected eigenvalues.
C
      if( abs(QQUASI).eq.4 .and. RESTO.eq.0 ) then
         p_corr  = p_rwend
         p_rwend = p_corr + NIND
         p_icorr = p_iwend
         p_iwend = p_icorr + NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         elseif( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
      endif
C
C     Perform the eigenvalue correction if necessary:
C
      lam_b = 1.d200    ! This is just for output...
C
C     Get smallest and largest eigenvalue (in absolute values)
C
      eigmin = 1.d200
      do i = 1, NIND
         eigmin = dmin1(dabs(RW(p_eigs+i)),eigmin)
      enddo
      i = IDAMAX(NIND, RW(p_eigs+1), 1)
      eigmax = dabs(RW(p_eigs+i))
CTODO Need to find good rule for mineig
ccc         mineig = 1.d-12*eigmax
ccc         mineig = dmin1(1d-6*MU,0.01d0)
ccc         mineig = dmax1(1d-6*dmin1(MU,0.01d0),1.d-10*mineig)
c         mineig = dmin1(MU,0.01d0)
      mineig = dmin1(MU,1d-6)
CA         mineig = 1.d2*dmin1(MU,1d-6)
ccc         mineig = dmax1(dmin1(MU,0.01d0),1.d-5)
C         mineig = dmin1(MU,0.01d0)/(iter+1)
CCC         mineig = dmin1(dsqrt(MU),0.01d0)
C      write(line,*) 'eigmin = ',eigmin,' eigmax = ',eigmax,
C     1     ' mineig = ',mineig
C      call C_OUT(1,3,1,line)
      write(line,*) 'Currently: mineig = dmin1(MU,1d-6)'
      call C_OUT(1,3,1,line)
C
C     Here is the correction:
C
      do i = 1, NIND
         if( RW(p_eigs+i).lt.mineig ) then
            IEIGS = IEIGS + 1
            if( lam_b.gt.RW(p_eigs+i) ) then
               lam_b = RW(p_eigs+i)
            endif
            eigold = RW(p_eigs+i)
            RW(p_eigs+i) = dmax1(mineig, dabs(RW(p_eigs+i)))
CTODO Need to find good rule for eigenvalue correction
C               RW(p_eigs+i) = dmin1(dmax1(mineig, dabs(RW(p_eigs+i))),
C     1                              1000d0)
            if( abs(QQUASI).eq.4 .and. RESTO.eq.0 ) then
               IW(p_icorr+ IEIGS) = i
               RW(p_corr + IEIGS) = RW(p_eigs+i) - eigold
            endif
         endif
      enddo
      if( QCNR.gt.0 .and. QPRINT.ge.2) then
         write(line,2115) ITER, IEIGS, lam_b
 2115    format(/,'  Information about SR1 correction in ITER',
     1        i5,//,' ',i3,' corrections with worst violation ',
     2        d16.10)
         call C_OUT(1,2,4,line)
      endif
C
C     Perform correction of B as desired
C
      if( abs(QQUASI).eq.4 .and. RESTO.eq.0 ) then
         do i = 1, IEIGS
            k = IW(p_icorr+i)
            call DSPR('U', NIND, RW(p_corr+i), RW(p_z+1+NIND*(k-1)),
     1           1, B)
         enddo
         p_iwend = p_icorr
         p_rwend = p_corr
      endif
C
C     Solve system using eigenvalue decomposition
C
      p_tmp   = p_rwend
      p_rwend = p_tmp + NIND
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      if( QDAMP.ne.0 .and. RESTO.eq.0 ) then
         p_bzg   = p_rwend
         p_rwend = p_bzg + NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         call DGEMV('T', NIND, NIND, 1.d0, RW(p_z+1), NIND, RGB,
     1        1, 0.d0, RW(p_tmp+1), 1)
         do i = 1, NIND
            RW(p_tmp+i) = RW(p_tmp+i)/RW(p_eigs+i)
         enddo
         call DGEMV('N', NIND, NIND, 1.d0, RW(p_z+1), NIND, RW(p_tmp+1),
     1        1, 0.d0, RW(p_bzg+1), 1)
         gzbw = DDOT(NIND, WCORR, 1, RW(p_bzg+1), 1)
         if( gzbw.lt.0.d0 ) then
            gzbzg = DDOT(NIND, RGB, 1, RW(p_bzg+1), 1)
            zeta = min(1.d0, -0.1d0*gzbzg/gzbw)
            if( zeta.lt.1.d0 ) then
               write(line,*) 'Damping active: ZETA = ',zeta
               call C_OUT(2,0,1,line)
               call DAXPY(NIND, 1.d0-zeta, WCORR, 1, PZ, 1)
            endif
         endif
         p_rwend = p_bzg
      endif
CTODO could use different routine here...
      call DGEMV('T', NIND, NIND, 1.d0, RW(p_z+1), NIND, PZ,
     1     1, 0.d0, RW(p_tmp+1), 1)
      do i = 1, NIND
         RW(p_tmp+i) = RW(p_tmp+i)/RW(p_eigs+i)
      enddo
      call DGEMV('N', NIND, NIND, 1.d0, RW(p_z+1), NIND, RW(p_tmp+1),
     1     1, 0.d0, PZ, 1)
      p_rwend = p_eigs
      goto 9999

 2200 continue
C
C     Add multiple of identity
C     ========================
CTODO Think aboput staring with frobenius norm as first trial for correction
C
      l_regu = .false.
      if( SOC_FLAG.eq.0 ) then
         REGU = 0d0
      endif
C
C     First try Cholesky; if not positive definite, add multiple of Z^T Z
C
      p_wtmp  = p_rwend
      p_rwend = p_wtmp + (NIND*(NIND+1))/2
      if( lmodchol ) then
         p_iperm = p_iwend
         p_iwend = p_iperm + NIND
      endif
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
 2210 continue
      call DCOPY((NIND*(NIND+1))/2, W, 1, RW(p_wtmp+1), 1)
C
C     Try to solve system using Cholesky factorization
C     If successful, wtmp will contain Cholesky factors
C
      if( lmodchol ) then
         info = -1
         call MOD_CHOL_FAC(NIND, RW(p_wtmp+1), IW(p_iperm+1), info)
      else
         call DPPTRF('U', NIND, RW(p_wtmp+1), info)
         if( info.lt.0 ) then
            write(line,*) 'get_pz: DPPTRF returned INFO = ', info
            call C_OUT(2,0,1,line)
            IERR = 502
            goto 9999
         endif
      endif
      if( info.ne.0 ) then      ! not positive definite
C
C     Need to regularize
C
         write(line,*) '   INFO from Cholesky = ',info,'; regularize...'
         call C_OUT(1,2,1,line)
C
C     If the first trial in this iteration, get Z^T Z
C
         if( .not.l_regu .and. M.gt.0 .and. abs(QCORRECT).eq.4 ) then
            p_ztz   = p_rwend
            p_sigma = p_ztz   + (NIND*(NIND+1))/2
            p_rwend = p_sigma + N
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
C
C     e to simga
C
            call DCOPY(N, 1d0, 0, RW(p_sigma+1), 1)
C
C     Call CONSTR to get dependent part of BB
C
            call CONSTR(6, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, RW(p_sigma+1), RW(p_ztz+1),
     2           idummy, idummy,
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.lt.0 ) then
               write(line,*) 'get_pz: Warning in CONSTR, IERR = ',
     1              IERR
               call C_OUT(2,0,1,line)
               IERR = 0
            elseif( IERR.ne.0 ) then
               write(line,*) 'get_pz: Error in CONSTR, IERR = ',
     1              IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
C
C     Add independent part to BB
C
            k = 0
            do i = 1, NIND
               k = k + i
               RW(p_ztz+k) = RW(p_ztz+k) + RW(p_sigma+M+i)
            enddo
            p_rwend = p_sigma
C
         endif
C
         regu1 = REGU
         l_regu = .true.
         if( REGU.eq.0d0 ) then
            if( REGU_STORE.eq.0d0 ) then
               REGU = REGU_INIT
            else
               REGU = REGU_DEC_FACT * REGU_STORE
            endif
         else
            if( REGU_STORE.eq.0d0 ) then
               REGU = REGU_INIT_FACT * REGU
            else
               REGU = REGU_INC_FACT  * REGU
            endif
         endif
         IEIGS = IEIGS + 1
         if( REGU.gt.REGU_MAX ) then
            IERR = 10
            write(line,*) 'regu getting too large:',REGU
            call C_OUT(2,0,1,line)
            goto 9999
         endif
C
C     Add multiple of identity to W
C
         diff = REGU - regu1
         if( M.gt.0 .and. abs(QCORRECT).eq.4 ) then
            call DAXPY((NIND*(NIND+1))/2, diff, RW(p_ztz+1), 1, W, 1)
         else
            k = 0
            do i = 1, NIND
               k = k + i
               W(k) = W(k) + diff
            enddo
         endif
C
C     ... and factorize again
C
         goto 2210

      else                      ! positive definite
C
C     If there has been a regularization, we need to recompute WCORR
C     and solve system with that new WCORR
C
         if( l_regu ) then
C
C     If wanted, change estimate B corresponding to correction:
C
            if( abs(QQUASI).eq.4 .and.RESTO.eq.0 ) then
               if( M.gt.0 .and. abs(QCORRECT).eq.4 ) then
                  call DAXPY( (NIND*(NIND+1))/2, REGU, RW(p_ztz+1), 1,
     1                 B, 1)
               else
                  k = 0
                  do i = 1, NIND
                     k = k + i
                     B(k) = B(k) + REGU
                  enddo
               endif
            endif
            if( M.gt.0 ) then
               p_rwend = p_ztz
            endif
C
C     Correct WCORR
C
            if( abs(QCORRECT).eq.4 ) then
               call GET_WCORR(N, NIND, M, X, ITER, IVAR, NFIX, IFIX,
     1              NORIG, XORIG, CSCALE, LAM, NLB, ILB, NUB, IUB,
     2              S_L, S_U, BNDS_L, BNDS_U,
     1              SIGMA_L, SIGMA_U, YPY, regu, WCORR, RESTO,
     1              KCONSTR, LRS, RS, LIS, IS,
     2              LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     3              IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5              EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
               if( IERR.ne.0 ) then
                  write(line,*) 'get_pz: get_wcorr returns IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  goto 9999
               endif
C
C     Recompute RHS:  PZ = - ( RGB + WCORR_new )
C
               if( RESTO.eq.0 ) then
                  call DCOPY(NIND, RGB, 1, PZ, 1)
                  call DAXPY(NIND, 1d0, WCORR, 1, PZ, 1)
               else
                  call DCOPY(NIND, WCORR, 1, PZ, 1)
               endif
            call DSCAL(NIND, -1.d0, PZ, 1)
            endif
C
            write(line,2230) REGU
 2230       format(/,'  Regularization: REGU = ',d12.5)
            call C_OUT(1,2,2,line)
C
            REGU_STORE = REGU
C
         endif
C
C     Solve for new RHS with Cholseky factors stored in wtmp
C
         if( lmodchol ) then
            p_w = p_rwend
            p_rwend = p_w + NIND
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            if( QDAMP.ne.0 .and. RESTO.eq.0 ) then
               p_bzg   = p_rwend
               p_rwend = p_bzg + NIND
               if( p_rwend.gt.LRW ) then
                  IERR = 98
                  goto 9999
               endif
               call DCOPY(NIND, RGB, 1, RW(p_bzg+1), 1)
               call MOD_CHOL_SOL(NIND, RW(p_wtmp+1), IW(p_iperm+1),
     1              RW(p_bzg+1), RW(p_w+1))
               gzbw = DDOT(NIND, WCORR, 1, RW(p_bzg+1), 1)
               if( gzbw.lt.0.d0 ) then
                  gzbzg = DDOT(NIND, RGB, 1, RW(p_bzg+1), 1)
                  zeta = min(1.d0, -0.1d0*gzbzg/gzbw)
                  if( zeta.lt.1.d0 ) then
                     write(line,*) 'Damping active: ZETA = ',zeta
                     call C_OUT(2,0,1,line)
                     call DAXPY(NIND, 1.d0-zeta, WCORR, 1, PZ, 1)
                  endif
               endif
               p_rwend = p_bzg
            endif
            call MOD_CHOL_SOL(NIND, RW(p_wtmp+1), IW(p_iperm+1),
     1           PZ, RW(p_w+1))
            p_iwend = p_iperm
            p_rwend = p_w
         else
            if( QDAMP.ne.0 .and. RESTO.eq.0 ) then
               p_bzg   = p_rwend
               p_rwend = p_bzg + NIND
               if( p_rwend.gt.LRW ) then
                  IERR = 98
                  goto 9999
               endif
               call DCOPY(NIND, RGB, 1, RW(p_bzg+1), 1)
               call DPPTRS( 'U', NIND, 1, RW(p_wtmp+1), RW(p_bzg+1),
     1              NIND, info )
               if( info.ne.0 ) then
                  write(line,*)
     1                 'get_pz: DPPTRS(d) returned INFO = ', info
                  call C_OUT(2,0,1,line)
                  IERR = 503
                  goto 9999
               endif
               gzbw = DDOT(NIND, WCORR, 1, RW(p_bzg+1), 1)
               if( gzbw.lt.0.d0 ) then
                  gzbzg = DDOT(NIND, RGB, 1, RW(p_bzg+1), 1)
                  zeta = min(1.d0, -0.1d0*gzbzg/gzbw)
                  if( zeta.lt.1.d0 ) then
                     write(line,*) 'Damping active: ZETA = ',zeta
                     call C_OUT(2,0,1,line)
                     call DAXPY(NIND, 1.d0-zeta, WCORR, 1, PZ, 1)
                  endif
               endif
               p_rwend = p_bzg
            endif
            call DPPTRS( 'U', NIND, 1, RW(p_wtmp+1), PZ, NIND, info )
            if( info.ne.0 ) then
               write(line,*) 'get_pz: DPPTRS returned INFO = ', info
               call C_OUT(2,0,1,line)
               IERR = 503
               goto 9999
            endif
         endif
         p_rwend = p_wtmp
C
      endif
      goto 9999
C
C     Try LU decomposition, and only bother, if we encounter direction
C     of negative curvature
C     ================================================================
C
 2300 continue
      if( QDAMP.ne.0 ) then
         call C_OUT(2,0,1,
     1        'get_pz: QDAMP not implemented for LU factors')
         IERR = 4
         goto 9999
      endif
C
C     Do solution of system by LU decomposition
C
C     Reserve work space
C
      p_a     = p_rwend
      p_rwend = p_a     + NIND*NIND
      p_ipiv  = p_iwend
      p_iwend = p_ipiv  + NIND
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Copy B into A as square matrix
C
      k = 1
      do j = 1, NIND
         do i = 1, j
            RW(p_a+i+NIND*(j-1)) = W(k)
            RW(p_a+j+NIND*(i-1)) = W(k)
            k = k + 1
         enddo
      enddo
C
C     Now do LU decomposition
C
      call DGESV( NIND, 1, RW(p_a+1), NIND, IW(p_ipiv+1),
     1     PZ, NIND, info)
      p_iwend = p_ipiv
      p_rwend = p_a
      if( info.lt.0 ) then
         write(line,*) 'in get_pz: DGESV returns info = ',info
         call C_OUT(2,0,1,line)
         IERR = 504
         goto 9999
      elseif( info.gt.0 ) then  ! system is singular, do other correction
         goto 2400
      endif
C
C     Check if PZ is direction of negative curvature
C
      p_bpz   = p_rwend
      p_rwend = p_bpz + NIND
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call DSPMV('U', NIND, 1.d0, W, PZ, 1, 0.d0, RW(p_bpz+1), 1)
      val1 = DDOT(NIND, PZ, 1, RW(p_bpz+1), 1)
      val2 = DDOT(NIND, PZ, 1, PZ, 1)
      p_rwend = p_bpz
      if( val1 .le. 1d-8*val2 ) then
         write(line,*) 'val1 = ',val1,' val2 = ',val2,
     1        ' -> neg curv...'
         call C_OUT(1,1,1,line)
         goto 2400
      endif

      goto 9999

 2400 continue
C
C     Restore RHS:  PZ = - ( RGB + WCORR )
C
      if( RESTO.eq.0 ) then
         call DCOPY(NIND, RGB, 1, PZ, 1)
         call DAXPY(NIND, 1d0, WCORR, 1, PZ, 1)
      else
         call DCOPY(NIND, WCORR, 1, PZ, 1)
      endif
      call DSCAL(NIND, -1.d0, PZ, 1)
C
C     Go to other corrections...
C
      goto 2050

C ==============================================================================
C
C     END   : SR1 case
C
C ==============================================================================

C
C     I think that's it...
C
 9999 continue

      if( IEIGS.ne.0 ) then
         COUNT_NEG_CURV = COUNT_NEG_CURV + 1
      endif

      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine GET_PZ_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer nind, lrw1, liw1, lrw2, liw2, lrw3, liw3
      integer lrw_constr, liw_constr
      logical lmodchol          ! set as in code above
      parameter( lmodchol = .false. )
      character*80 line

      nind = N-M

      LIW = 0
      LRW = 0

C     Label 1000
      if( abs(QQUASI).eq.1 .or. abs(QQUASI).eq.3 .or.
     1     abs(QQUASI).eq.5) then
         if( QDAMP.eq.0 ) then
            LRW = max(LRW, (nind*(nind+1))/2+nind)
            LIW = max(LIW, nind)
         else
            LRW = max(LRW, (nind*(nind+1))/2+2*nind)
            LIW = max(LIW, nind)
         endif
      endif

      if( QCORRECT.lt.0 ) then
C     Label 2300
         LRW = max(LRW, nind*nind)
         LIW = max(LIW, nind)
      endif

      if( abs(QCORRECT).eq.2 ) then
C     Label 2100
         LRW = max(LRW, 4*nind +nind*nind+(nind*(nind+1))/2)
         if( abs(QQUASI).eq.4 ) then
            LIW = max(LIW, nind)
         endif
      elseif( abs(QCORRECT).eq.3 ) then
C     Label 2200
         lrw1 = (nind*(nind+1))/2
         if( lmodchol ) then
            liw1 = nind
         else
            liw1 = 0
         endif
         lrw2 = 0
         liw2 = 0
         if( M.gt.0 .and. abs(QCORRECT).eq.4 ) then
            call CONSTR_WS(N, M, NLB, NUB, NZA, lrw_constr,
     1           liw_constr, DAT, IDAT)
            lrw2 = n + (nind*(nind+1))/2 + lrw_constr
            liw2 = liw_constr
         endif
         if( abs(QCORRECT).eq.4 ) then
            call GET_WCORR_WS(N, M, NLB, NUB, NZA, lrw3, liw3,
     1           DAT, IDAT)
            lrw2 = max(lrw2, lrw3)
            liw2 = max(liw2, liw3)
         endif
         if( lmodchol ) then
            if( QDAMP.ne.0 ) then
               lrw2 = max(lrw2, 2*nind)
            else
               lrw2 = max(lrw2, nind)
            endif
         else
            if( QDAMP.ne.0 ) then
               lrw2 = max(lrw2, nind)
            endif
         endif
         LRW = max(LRW, lrw1+lrw2)
         LIW = max(LIW, liw1+liw2)
      endif

      if( QPRINT.ge.4 ) then
         write(line,1000)'get_pz_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end

