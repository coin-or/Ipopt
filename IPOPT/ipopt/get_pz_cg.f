C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_PZ_CG(ITER, MU, NIND, RGB, WCORR, PREC, PZ,
     1     CGITER, N, M, X, IVAR, NFIX, IFIX, NORIG, XORIG, CSCALE, LAM,
     2     NLB, ILB, NUB, IUB, S_L, S_U, SIGMA_L, SIGMA_U,
     1     ERR, ERR_BAR, NEWBAS, KCONSTR, LRS, RS, LIS, IS, LRW, RW,
     1     LIW, IW, IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     1     EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: get_pz_cg.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Solve for null space step pZ using preconditioned CG
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
CP   RGB       I    DP     reduced gradient of barrier function
CP   WCORR     I    DP     correction term
CP   PREC      I    DP     Quasi-Newton approx for preconditioner
CP   PZ        O    DP     null space step pY (only for independent variables)
CP   CGITER    O    INT    number of CG iterations (if negative, then direction
CP                            of negative curvature has been encountered)
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
CP   SIGMA_L   I    DP     primal-dual Hessian of lower bound barrier term
CP                            (NLB diagonal elements only)
CP   SIGMA_U   I    DP     primal-dual Hessian of upper bound barrier term
CP                            (NUB diagonal elements only)
CP   ERR       I    DP     current KKT error for barrier problem
CP   ERR_BAR   I    DP     KKT error tolerance for current barrier problem
CP   NEWBAS   I/O   LOG    If true, basis has just changed, we should therefore
CP                           reset the preconditioner estimate
CP                         Will be set to true, if basis seems bad
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
CS    DDOT
CS    DNRM2
CS    DSCAL
CS    DCOPY
CS    DAXPY
CS    DSPMV
CS    C_OUT
CS    CONSTR
CS    GET_ZWZV
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
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer ITER
      double precision MU
      integer NIND
      double precision RGB(NIND)
      double precision WCORR(NIND)
      double precision PREC((NIND*(NIND+1))/2)
      double precision PZ(NIND)
      integer CGITER
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
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision ERR
      double precision ERR_BAR
      logical NEWBAS
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
      double precision DDOT, DNRM2, DLAMCH

      integer p_iwend, p_rwend, p_prec, p_p, p_y, p_ap, p_bs, p_r
      integer p_apf, p_sc, p_tmp, p_work, p_iwork, p_zwzp, p_w, p_bb
      integer i, k, itermax, info, task, idummy, neigs
      double precision ry, ry_new, sbs, pap, tol, alpha, beta, rnrm
      double precision rnrm0, dummy, rcond, ferr, berr, theta, dot
      double precision times, timef, abstol, mineig, lambda, fact
      double precision snorm, ynorm, timechols, timecholf, increase_diag
      character*120 line
      character*1 equed
      integer damped

C
C     |QCG|     1: QN on inverse of overall reduced Hessian
C               2: QN (damped BFGS) on red. original Hessian;
C                      add barrier part explicitly
C               3: QN (SR1) on red. original Hessian;
C                      add barrier part explicitly
C               4: No preconditioner at all
C              >0: start CG with PZ = 0
C              <0: start CG with PZ = -B^{1}( Z^T g + wcorr )
C
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      call TIMER(times)

      IERR    = 0
      p_rwend = 0
      p_iwend = 0
      damped  = 0
C
      if( MEMDBG ) then
         write(line,1)'get_pz_cg', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif
C
C     There is nothing to do if there are no degrees of freedom
C
      if( NIND.eq.0 ) goto 9999
C
C     If first iteration, initialize preconditioner as identity(?)
C
      if( ITER.eq.0 .or. QCG.eq.4 .or. NEWBAS ) then
         call DCOPY( (NIND*(NIND+1))/2, 0.d0, 0, PREC, 1)
         k = 0
         do i = 1,NIND
            k = k + i
            PREC(k) = 1.d0
         enddo
      endif
C
C     get work space
C
      p_prec  = p_rwend
      p_p     = p_prec + (NIND*(NIND+1))/2
      p_y     = p_p    + NIND
      p_r     = p_y    + NIND
      p_ap    = p_r    + NIND
      p_rwend = p_ap   + NIND
      if( abs(QCG).eq.2 .or. abs(QCG).eq.3 ) then
         p_zwzp  = p_rwend
         p_bb    = p_zwzp  + NIND
         p_rwend = p_bb + (NIND*(NIND+1))/2
      else
         p_zwzp  = 0
      endif
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
C     Determine current preconditioner
C
      increase_diag = 1d-5
 50   continue
      if( abs(QCG).eq.1 ) then
C
C     Use QN estimate of inverse of Hessian
C
         call DCOPY( (NIND*(NIND+1))/2, PREC, 1, RW(p_prec+1), 1)

      elseif( abs(QCG).eq.2 .or. abs(QCG).eq.3 ) then
C
C     Compute Z^T Sigma Z
C
         p_apf   = p_rwend
         p_sc    = p_apf + (NIND*(NIND+1))/2
         p_rwend = p_sc  + NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
C
C     NOTE: IN THIS CALL WE ASSUME THAT S_L & S_U ARE NOT ACCESSED IN GET_BB
C           IF RESTO = .false. AS HERE!
C
         call GET_BB( N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE,
     1        NLB, ILB, NUB, IUB, SIGMA_L, SIGMA_U, dummy, dummy,
     1        RW(p_prec+1), 0, KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     3        IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'get_pz_cg: get_bb returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         call DCOPY( (NIND*(NIND+1))/2, PREC, 1, RW(p_bb+1), 1)
C
C     Add current approximation of original Hessian
C
         call DAXPY( (NIND*(NIND+1))/2, 1.d0, PREC, 1, RW(p_prec+1), 1)
C
C     If SR1 check, if multiple of identity needs to be added in order to
C     make matrix positive definite
C
         if( abs(QCG).eq.3 ) then

            p_w     = p_rwend
            p_work  = p_w    + NIND
            p_rwend = p_work + 8*NIND
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            p_iwork = p_iwend
            p_iwend = p_iwork + 5*NIND
            if( p_iwend.gt.LIW ) then
               IERR = 99
               goto 9999
            endif
            abstol = 2*DLAMCH('S')
            call DCOPY( (NIND*(NIND+1))/2, RW(p_prec+1), 1,
     1           RW(p_apf+1), 1)
            call DSPEVX('N', 'I', 'U', NIND, RW(p_apf+1), 0.d0, 0.d0,
     1           1, 1, abstol, neigs, RW(p_w+1), dummy, 1, RW(p_work+1),
     2           IW(p_iwork+1), idummy, info)
            if( info.ne.0 ) then
               write(line,*) 'get_pz_cg: DSPEVX returns info = ',info
               call C_OUT(2,0,1,line)
               IERR = 498
               goto 9999
            endif
            mineig = RW(p_w+1)
            write(line,*) 'mineig = ',mineig
            call C_OUT(2,0,1,line)
            lambda = dmax1(0.d0, dmax1(MU,
     1           dmin1(1.d0,1.d-1*dabs(mineig))) - mineig)
            if( lambda.gt.0d0 ) then
               write(line,*)
     1              'get_pz_cg: correction SR1 prec with lambda = ',
     2              lambda
               call C_OUT(2,0,1,line)
               k = 0
               do i = 1,NIND
                  k = k + i
                  RW(p_prec+k) = RW(p_prec+k) + lambda
               enddo
            endif
            p_rwend = p_w
            p_iwend = p_iwork
         endif
C
C     Compute Cholesky factorization
C
         p_work  = p_rwend
         p_rwend = p_work + 3*NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         p_iwork = p_iwend
         p_iwend = p_iwork + NIND
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         call timer(timechols)
         call DPPSVX('E', 'U', NIND, 0, RW(p_prec+1), RW(p_apf+1),
     1        equed, RW(p_sc+1), dummy, NIND, dummy, NIND, rcond,
     2        ferr, berr, RW(p_work+1), IW(p_iwork+1), info)
         WRITE(line,*) 'rcond = ',rcond
         call C_OUT(2,0,1,line)
         call timer(timecholf)
         TIME_PZ_CHOL = TIME_PZ_CHOL + timecholf - timechols
         if( info.ne.0 ) then
            write(line,*) 'get_pz_cg: DPPSVX returns info = ',info
            call C_OUT(2,0,1,line)
            write(line,*) 'get_pz_cg: try to increase PREC '
            call C_OUT(2,0,1,line)
            k = 0
            increase_diag = increase_diag * 10d0
            NEWBAS = .true.
            if( increase_diag.gt.1.d100 ) then
               call C_OUT(2,0,1,
     1              'get_pz_cg: increase_diag getting too large.')
               IERR = 673
               goto 9999
            endif
            do i = 1,NIND
               k = k + i
               PREC(k) = PREC(k) + increase_diag
            enddo
            goto 50
         endif
         p_rwend = p_work
         p_iwend = p_iwork

      elseif( QCG.ne.4 ) then
         call C_OUT(2,0,1,'get_pz_cg: Invalid QCG.')
         IERR = 4
         goto 9999
      endif
C
C     Prepare CG (notation like in Nocedal's book, p.118/119)
C
      CGITER = 0
C
C     Compute RHS
C
      call DCOPY(NIND, RGB, 1, RW(p_r+1), 1)
      call DAXPY(NIND, 1.d0, WCORR, 1, RW(p_r+1), 1)
C
C     Starting point for CG
C
      if( QCG.gt.0 .or. ITER.eq.0 .or. NEWBAS ) then
         call DCOPY(NIND, 0.d0, 0, PZ, 1)
      else
         if( QCG.eq.-1) then
            call DSPMV('U', NIND, 1.d0, RW(p_prec+1), RW(p_r+1), 1,
     1           0.d0, PZ, 1)
            call DSCAL(NIND, -1.d0, PZ, 1)
         elseif( QCG.eq.-2 .or. QCG.eq.-3 ) then
            p_tmp   = p_rwend
            p_work  = p_tmp  + NIND
            p_rwend = p_work + 3*NIND
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            p_iwork = p_iwend
            p_iwend = p_iwork + NIND
            if( p_iwend.gt.LIW ) then
               IERR = 99
               goto 9999
            endif
            call DCOPY(NIND, RW(p_r+1), 1, RW(p_tmp+1), 1)
            call DPPSVX('F', 'U', NIND, 1, RW(p_prec+1), RW(p_apf+1),
     1           equed, RW(p_sc+1), RW(p_tmp+1), NIND, PZ, NIND,
     1           rcond, ferr, berr, RW(p_work+1), IW(p_iwork+1), info)
            if( info.ne.0 ) then
               write(line,*) 'get_pz_cg: DPPSVX returns info = ',info
               call C_OUT(2,0,1,line)
               IERR = 498
               goto 9999
            endif
            p_rwend = p_tmp
            p_iwend = p_iwork
            call DSCAL(NIND, -1.d0, PZ, 1)
         endif
         p_tmp   = p_rwend
         p_rwend = p_tmp  + NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         call GET_ZWZV(1, N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, X, CSCALE, NLB, ILB, NUB, IUB, S_L, S_U,
     2        SIGMA_L, SIGMA_U, LAM, PZ, RW(p_tmp+1), dummy,
     1        KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     3        IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.lt.0 ) then
            write(line,*) 'get_pz_cg: Warning in get_zwzv, IERR = ',IERR
            call C_OUT(2,0,1,line)
         elseif( IERR.ne.0 ) then
            write(line,*) 'get_pz_cg: Error in get_zwzv, IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         call DAXPY(NIND, 1.d0, RW(p_tmp+1), 1, RW(p_r+1), 1)
         p_rwend = p_tmp
      endif
C
C     y and p
C
      if( abs(QCG).eq.1 ) then
         call DSPMV('U', NIND, 1.d0, RW(p_prec+1), RW(p_r+1), 1, 0.d0,
     1        RW(p_y+1), 1)
      elseif( abs(QCG).eq.2 .or. abs(QCG).eq.3 ) then
         p_tmp   = p_rwend
         p_work  = p_tmp  + NIND
         p_rwend = p_work + 3*NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         p_iwork = p_iwend
         p_iwend = p_iwork + NIND
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         call DCOPY(NIND, RW(p_r+1), 1, RW(p_tmp+1), 1)
         call DPPSVX('F', 'U', NIND, 1, RW(p_prec+1), RW(p_apf+1),
     1        equed, RW(p_sc+1), RW(p_tmp+1), NIND, RW(p_y+1), NIND,
     1        rcond, ferr, berr, RW(p_work+1), IW(p_iwork+1), info)
         if( info.ne.0 ) then
            write(line,*) 'get_pz_cg: DPPSVX returns info = ',info
            call C_OUT(2,0,1,line)
            IERR = 498
            goto 9999
         endif
         p_rwend = p_tmp
         p_iwend = p_iwork
      elseif( QCG.eq.4 ) then
         call DCOPY(NIND, RW(p_r+1), 1, RW(p_y+1), 1)
      endif

      call DCOPY(NIND, RW(p_y+1), 1, RW(p_p+1), 1)
      call DSCAL(NIND, -1.d0, RW(p_p+1), 1)
      ry = DDOT(NIND, RW(p_y+1), 1, RW(p_r+1), 1)

C
C     decide, when to stop CG
C     if QCGTOL < 0, use preconditioned residual, otherwise normal residual
C
      if( QCGTOL.lt.0.d0 ) then
         rnrm0 = sqrt(ry)
      else
         rnrm0 = DNRM2(NIND, RW(p_r+1), 1)
      endif
      tol     = abs(QCGTOL)*rnrm0
CTODO Decide if this is good...
C      tol     = min(tol, ERR_BAR)
      tol     = min(tol, max(dsqrt(ERR_BAR), 1d2*ERR_BAR))

      if( QMAXCGITER.gt.0 ) then
         itermax = QMAXCGITER
      else
         itermax = NIND
      endif

      write(line,*) 'tol = ',tol, ' rnrm0 = ', rnrm0,' pznrm = ',
     1     DNRM2(NIND, PZ, 1)
      call C_OUT(2,1,1,line)
C
C     START of CG loop
C
 100  continue
      CGITER   = CGITER   + 1
      COUNT_CG = COUNT_CG + 1
C
C     Compute Ap = Z^T (Sigma + W) Z p  and if necessary also
C           zwzp = Z^T W Z p
C
      if( abs(QCG).eq.2 .or. abs(QCG).eq.3 ) then
         task = 3
      else
         task = 1
      endif
      call GET_ZWZV(task, N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, X, CSCALE, NLB, ILB, NUB, IUB, S_L, S_U,
     2     SIGMA_L, SIGMA_U, LAM, RW(p_p+1), RW(p_ap+1), RW(p_zwzp+1),
     1     KCONSTR, LRS, RS, LIS, IS,
     2     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR,
     2     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     1     DAT, IDAT)
      if( IERR.lt.0 ) then
         write(line,*) 'get_pz_cg: Warning in get_zwzv, IERR = ',IERR
         call C_OUT(2,0,1,line)
      elseif( IERR.ne.0 ) then
         write(line,*) 'get_pz_cg: Error in get_zwzv, IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     pap = p^T * A * p
C
      pap = DDOT(NIND, RW(p_p+1), 1, RW(p_ap+1), 1)
C
      if( pap.gt.0.d0 ) then

         alpha = ry / pap
         call DAXPY(NIND, alpha, RW(p_p+1), 1, PZ, 1)
         call DAXPY(NIND, alpha, RW(p_ap+1), 1, RW(p_r+1), 1)
         if( abs(QCG).eq.1 ) then
            call DSPMV('U', NIND, 1.d0, RW(p_prec+1), RW(p_r+1), 1,
     1           0.d0, RW(p_y+1), 1)
         elseif( abs(QCG).eq.2 .or. abs(QCG).eq.3 ) then
            p_tmp   = p_rwend
            p_work  = p_tmp  + NIND
            p_rwend = p_work + 3*NIND
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            p_iwork = p_iwend
            p_iwend = p_iwork + NIND
            if( p_iwend.gt.LIW ) then
               IERR = 99
               goto 9999
            endif
            call DCOPY(NIND, RW(p_r+1), 1, RW(p_tmp+1), 1)
            call DPPSVX('F', 'U', NIND, 1, RW(p_prec+1), RW(p_apf+1),
     1           equed, RW(p_sc+1), RW(p_tmp+1), NIND, RW(p_y+1), NIND,
     1           rcond, ferr, berr, RW(p_work+1), IW(p_iwork+1), info)
            if( info.ne.0 ) then
               write(line,*) 'get_pz_cg: DPPSVX returns info = ',info
               call C_OUT(2,0,1,line)
               IERR = 498
               goto 9999
            endif
            p_rwend = p_tmp
            p_iwend = p_iwork
         elseif( QCG.eq.4 ) then
            call DCOPY(NIND, RW(p_r+1), 1, RW(p_y+1), 1)
         endif

      endif
C
C     Do the QN update for the preconditioner
C
      if( abs(QCG).eq.1 ) then
C
C     Perform BFGS update for inverse of preconditioner
C     (here s = ap, y = p for s,y in BFGS formula)
c$$$C     Perform inverse BFGS update for inverse of preconditioner
c$$$C     (see (8.16) in NocSWriBook)
C
c$$$         snorm = DNRM2(NIND, RW(p_ap+1), 1)
c$$$         ynorm = DNRM2(NIND, RW(p_p+1), 1)
         snorm = DNRM2(NIND, RW(p_p+1), 1)
         ynorm = DNRM2(NIND, RW(p_ap+1), 1)
         if( pap.gt.1.d-8*snorm*ynorm ) then

            if( .false. ) then
            p_bs    = p_rwend
            p_rwend = p_bs + NIND
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
C     bs = B*s (using symmetry of B)
            call DSPMV('U', NIND, 1.d0, PREC, RW(p_ap+1), 1, 0.d0,
     1           RW(p_bs+1), 1)
C     sBs = s'*bs
            sBs = DDOT(NIND, RW(p_ap+1), 1, RW(p_bs+1), 1)
C     B_tmp = B + y*y'/dot
            call DSPR('U', NIND,  1.d0/pap, RW(p_p+1), 1, PREC)
C     B_new = B_tmp -bs*bs'/dot
            call DSPR('U', NIND, -1.d0/sBs, RW(p_bs+1), 1, PREC)
            p_rwend = p_bs
            endif

            if( .true. ) then
            p_bs    = p_rwend
            p_rwend = p_bs + NIND
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
C     bs = PREC * y = PREC * ap
            call DSPMV('U', NIND, 1.d0, PREC, RW(p_ap+1), 1, 0.d0,
     1           RW(p_bs+1), 1)
C     sBs = ap'*bs = y' * PREC * y
            sBs = DDOT(NIND, RW(p_ap+1), 1, RW(p_bs+1), 1)
            fact = (sBs/pap + 1.d0)/pap
C     PREC <- PREC + fact s s^T
            call DSPR('U', NIND, fact, RW(p_p+1), 1, PREC)
C     PREC <- PREC - rho s (bs)^T - rho bs s^T
            fact = -1.d0/pap
            call DSPR2('U', NIND, fact, RW(p_bs+1), 1,
     1           RW(p_p+1), 1, PREC)
            p_rwend = p_bs
            endif
         else
            write(line,*) 'get_pz_cg: skip BFGS update b/c pap = ', pap
            call C_OUT(2,0,1,line)
         endif

      elseif( abs(QCG).eq.2 ) then
C
C     Do damped BFGS update for red Hessian of original Lagrangian
C     here, s = p and y = zwzp
C
         p_bs    = p_rwend
         p_rwend = p_bs + NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
C     bs = B*s (using symmetry of B)
         call DSPMV('U', NIND, 1.d0, PREC, RW(p_p+1), 1, 0.d0,
     1        RW(p_bs+1), 1)
C     sBs = s'*bs
         sBs = DDOT(NIND, RW(p_p+1), 1, RW(p_bs+1), 1)
C     dot = s'*y
         dot = DDOT(NIND, RW(p_p+1), 1, RW(p_zwzp+1), 1)
C
C     Do Powell damping
C
         if( dot.lt.0.2d0*sbs ) then
            damped = damped + 1
            theta = 0.8d0*sBs/(sBs - dot)
            call DSCAL(NIND, theta, RW(p_zwzp+1), 1)
            call DAXPY(NIND, (1d0-theta), RW(p_bs+1), 1,
     1           RW(p_zwzp+1), 1)
            dot = theta*dot + (1d0-theta)*sBs
            if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
               write(line,*) 'get_pz_cg: Powell damping theta = ',
     1              theta
               call C_OUT(1,2,1,line)
               write(line,*) 'get_pz_cg: damped s^T y = ', dot
               call C_OUT(1,2,1,line)
            endif
         endif
C     B_tmp = B + y*y'/dot
         call DSPR('U', NIND,  1.d0/dot, RW(p_zwzp+1), 1, PREC)
C     B_new = B_tmp -bs*bs'/dot
         call DSPR('U', NIND, -1.d0/sBs, RW(p_bs+1), 1, PREC)
 155     continue
         p_rwend = p_bs

      elseif( abs(QCG).eq.3 ) then
C
C     Do SR1 update for red Hessian of original Lagrangian
C     here, s = p and y = zwzp
C
         p_bs    = p_rwend
         p_rwend = p_bs + NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
C     bs = B*s (using symmetry of B)
         call DSPMV('U', NIND, 1.d0, PREC, RW(p_p+1), 1, 0.d0,
     1        RW(p_bs+1), 1)
C     bs = (B*s - y)
         call DAXPY(NIND, -1.d0, RW(p_zwzp+1), 1, RW(p_bs+1), 1)

C     dot = s'*bs
         dot = DDOT(NIND, RW(p_p+1), 1, RW(p_bs+1), 1)

C           snorm = |s|, ynorm = |y-Bs|
         snorm = DNRM2(NIND, RW(p_p+1), 1)
         ynorm = DNRM2(NIND, RW(p_bs+1), 1)

         if( dabs(dot).gt.dmax1(1.d-300,1.d-8*snorm*ynorm) ) then

C     B_new = B_tmp -bs*bs'/dot
            call DSPR('U', NIND, -1.d0/dot, RW(p_bs+1), 1, PREC)
         else
            write(line,710) ITER, dot, snorm, ynorm
 710        format(' In Iter ',i5,
     1           ' no SR1 update since dot = ',
     1           3d11.3)
            call C_OUT(1,2,0,line)
         endif
         p_rwend = p_bs

      endif
C
C     BFGS done
C

C
C     Leave if pap non-positive
C
      if( pap.le.0.d0 ) then
         COUNT_NEG_CURV = COUNT_NEG_CURV + 1
         CGITER = - CGITER
         goto 1000              ! direction of negative curvature!
      endif
C
      ry_new = DDOT(NIND, RW(p_y+1), 1, RW(p_r+1), 1)
C
      beta = ry_new / ry
      ry = ry_new
C
C
C     Check if to quit CG
C
      if( QCGTOL.lt.0.d0 ) then
         rnrm = sqrt(ry)
      else
         rnrm = DNRM2(NIND, RW(p_r+1), 1)
      endif
      write(line,*) 'cgiter = ',CGITER,', rnrm = ', rnrm,' pznrm = ',
     1     DNRM2(NIND, PZ, 1)
      call C_OUT(2,1,1,line)
      if( rnrm.lt.tol .or. CGITER.ge.itermax ) goto 1000
      call DSCAL(NIND, beta, RW(p_p+1), 1)
      call DAXPY(NIND, -1.d0, RW(p_y+1), 1, RW(p_p+1), 1)
      goto 100
C
C     END of CG loop
C
 1000 continue
C
C     Need to do something if negative curvature occoured at first CG iter
C
c      if( CGITER.lt.0 ) then
      if( CGITER.eq.-1 .and. (QCG.gt.0 .or. NEWBAS) ) then
C
C     Compute RHS
C
         call DCOPY(NIND, RGB, 1, RW(p_r+1), 1)
         call DAXPY(NIND, 1.d0, WCORR, 1, RW(p_r+1), 1)
         call DSCAL(NIND, -1.d0, RW(p_r+1), 1)

         if( abs(QCG).eq.1 ) then
C
C     Use QN estimate of inverse of Hessian
C
            call DSPMV('U', NIND, 1.d0, PREC, RW(p_r+1), 1, 0.d0, PZ, 1)

         elseif( abs(QCG).eq.2 ) then
C
C     Add current approximation of original Hessian and store result in bb
C
            call DAXPY( (NIND*(NIND+1))/2, 1.d0, PREC, 1, RW(p_bb+1), 1)
C
C     Compute Cholesky factorization
C
            p_work  = p_rwend
            p_rwend = p_work + 3*NIND
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            p_iwork = p_iwend
            p_iwend = p_iwork + NIND
            if( p_iwend.gt.LIW ) then
               IERR = 99
               goto 9999
            endif
            call timer(timechols)
            call DPPSVX('E', 'U', NIND, 1, RW(p_bb+1), RW(p_apf+1),
     1           equed, RW(p_sc+1), RW(p_r+1), NIND, PZ, NIND, rcond,
     2           ferr, berr, RW(p_work+1), IW(p_iwork+1), info)
            call timer(timecholf)
            TIME_PZ_CHOL = TIME_PZ_CHOL + timecholf - timechols
            if( info.ne.0 ) then
               write(line,*) 'get_pz_cg: DPPSVX returns info = ',info
               call C_OUT(2,0,1,line)
C               IERR = 499
C               goto 9999
               call C_OUT(2,0,1,' Setting PZ to zero and change basis.')
               call DCOPY(NIND, 0.d0, 0, PZ, 1)
               NEWBAS = .true.
            endif
            p_rwend = p_work
            p_iwend = p_iwork

         elseif( abs(QCG).ge.3 ) then
C
            write(line,*)
     1           'Neg. curvature not implemented for QCG = ',QCG
            call C_OUT(2,0,1,line)
            IERR = 4
            goto 9999
         endif

 5544    continue

      endif
C
C     I think that's it...
C
 9999 continue
      call TIMER(timef)
      TIME_CG = TIME_CG + timef - times
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine GET_PZ_CG_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer nind, lrw1, liw1, lrw2, liw2, lrw_zwzv, liw_zwzv
      character*80 line

      nind = N-M

      LIW = 0
      LRW = 0

      if( nind.eq.0 ) return

      call GET_ZWZV_WS(N, M, NLB, NUB, NZA, lrw_zwzv, liw_zwzv,
     1     DAT, IDAT)

      lrw1 = 4*nind + (nind*(nind+1))/2
      if( abs(QCG).eq.2 .or. abs(QCG).eq.3 ) then
         lrw1 = lrw1 + 2*nind + 2*(nind*(nind+1))/2
      endif
      liw1 = 0

      if( abs(QCG).eq.2 .or. abs(QCG).eq.3 ) then
         call GET_BB_WS(N, M, NLB, NUB, NZA, lrw2, liw2, DAT, IDAT)
         if( abs(QCG).eq.3 ) then
            lrw2 = max(lrw2, 9*nind)
            liw2 = max(liw2, 5*nind)
         endif
         lrw2 = max(lrw2, 3*nind)
         liw2 = max(liw2, nind)
      else
         lrw2 = 0
         liw2 = 0
      endif
      if( QCG.le.0 ) then
         if( QCG.eq.-2 .or. QCG.eq.-3 ) then
            lrw2 = max(lrw2, 4*nind)
            liw2 = max(liw2, nind)
         endif
         lrw2 = max(lrw2, nind+lrw_zwzv)
         liw2 = max(liw2, liw_zwzv)
      endif
      if( abs(QCG).eq.2 .or. abs(QCG).eq.3 ) then
         lrw2 = max(lrw2, 4*nind)
         liw2 = max(liw2, nind)
      endif
      lrw2 = max(lrw2, lrw_zwzv)
      liw2 = max(liw2, liw_zwzv)
      if( abs(QCG).eq.2 .or. abs(QCG).eq.3 ) then
         lrw2 = max(lrw2, 4*nind)
         liw2 = max(liw2, nind)
      endif
      if( abs(QCG).eq.1 ) then
         lrw2 = max(lrw2, nind)
      elseif( abs(QCG).eq.2 ) then
         lrw2 = max(lrw2, nind)
      elseif( abs(QCG).eq.3 ) then
         lrw2 = max(lrw2, nind)
      endif
      if( abs(QCG).eq.2 ) then
         lrw2 = max(lrw2, 3*nind)
         liw2 = max(liw2, nind)
      endif

      LRW = lrw1 + lrw2
      LIW = liw1 + liw2

      if( QPRINT.ge.4 ) then
         write(line,1000)'get_pz_cg_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
