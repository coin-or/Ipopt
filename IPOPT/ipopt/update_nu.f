C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine UPDATE_NU(ITER, N, NIND, M, X, NLB, ILB, NUB, IUB,
     1     BNDS_L, BNDS_U, S_L, S_U, SIGMA_L, SIGMA_U, MU, ERR,
     3     YPY, LAM, REGU, ZPZ, PZ, DX, GB, WCORR,
     3     F, C, CNRM0, WFLAG, NU, NUS, LRW, RW, IERR)
C
C*******************************************************************************
C
C    $Id: update_nu.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Update penalty parameter(s) for the line search
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
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of constraints
CP   X         I    DP     actual iterate (reordered without fixed vars:
CP                             first M entries belong to dependent
CP                             variables, remaining to independent variables)
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. S_L(i) is slack for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. S_U(i) is slack for X(IUB(i)) )
CP   BNDS_L    I    DP     values of lower bounds (ordered as S_L)
CP   BNDS_U    I    DP     values of upper bounds (ordered as S_U)
CP   S_L       I    DP     slacks to lower bounds
CP   S_U       I    DP     slacks to upper bounds
CP   SIGMA_L   I    DP     primal-dual Hessian of lower bound barrier term
CP                            (NLB diagonal elements only)
CP   SIGMA_U   I    DP     primal-dual Hessian of upper bound barrier term
CP                            (NUB diagonal elements only)
CP   MU        I    DP     barrier parameter
CP   ERR       I    DP     actual KKT-error (needed for switching on watchdog)
CP   YPY       I    DP     range space step (all variables; ordered like X)
CP   LAM       I    DP     multipliers for equality constraints
CP   REGU      I    DP     regularization factor (added regu*I to diagonal)
CP                            (from get_step_full)
CP   ZPZ       I    DP     null space step (only dependent variables)
CP                            (if NU(s) based on LAM)
CP   PZ        I    DP     null space step (only independent variables)
CP   DX        I    DP     step for X (primal)
CP   GB        I    DP     gradient of barrier function
CP   WCORR     I    DP     correction term for PZ
CP                         (For KNITRO-DIRECT update: d^T Hessian d)
CP   F         I    DP     value of objective function at X
CP   C         I    DP     values of constraints at X
CP   CNRM0     I    DP     2-norm of constraints at old point
CP   WFLAG     I    INT    flag for watchdog (see linesearch)
CP   NU       I/O   DP     actual value of penalty parameter (if only one)
CP                            (for QMERIT<0: largest entry in NUS)
CP   NUS      I/O   DP     actual values of individual penalty parameters
CP                            (if QMERIT<0)
CP   LRW       I    INT    length of RW
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
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
CS    DASUM
CS    DCOPY
CS    DAXPY
CS    IDAMAX
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
      integer ITER
      integer N
      integer NIND
      integer M
      double precision X(N)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision BNDS_L(NLB)
      double precision BNDS_U(NUB)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision MU
      double precision ERR
      double precision YPY(N)
      double precision LAM(M)
      double precision REGU
      double precision ZPZ(M)
      double precision PZ(NIND)
      double precision DX(N)
      double precision GB(N)
      double precision WCORR(NIND)
      double precision F
      double precision C(M)
      double precision CNRM0
      integer WFLAG
      double precision NU
      double precision NUS(M)
      integer LRW
      double precision RW(LRW)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            Local variables
C-------------------------------------------------------------------------------
C
      double precision DDOT, DASUM, CALC_NRM
      integer IDAMAX

      integer i, p_rwend
      double precision gpy, nunew, wpz, delta, cnrm1, lc, cj
      double precision gd, sigma, nu_trial
      character*120 line(4)

      integer QNUUPDATE         ! should later be parameter?!?
C     1: lambda-free without absolute values
C     2: lambda-free with absolute values
C     3: with lambdas
C     4: heuristic for full space?

C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0

C     If there are no constraints, we don't need a penalty parameter
C
      if( M.eq.0 ) then
         NU = 0d0
         goto 9999
      endif

      p_rwend = 0

      if( QMERIT.gt.0 ) goto 1000
      if( QMERIT.lt.0 ) goto 2000
C
C     No line search
C
      NU = 0d0
      goto 9999
C
C     ONLY ONE PENALTY PARAMETER
C
 1000 continue

      if( QFULL.eq.0 ) then
         if( QLAMBDA.eq.0 .or. QDAMP.eq.0 ) then
            QNUUPDATE = 1
         else
            QNUUPDATE = 3
         endif
      else
         if( QKNITROLS.eq.1 ) then
            QNUUPDATE = 4
         else
            QNUUPDATE = 3
         endif
      endif

C
C     Multiplier-free update
C
      if( QNUUPDATE.eq.1 .or. QNUUPDATE.eq.2 ) then
C
C     gpy = YPY'*grad(barrier function)
C
         gpy = DDOT(N, YPY, 1, GB, 1)
         if( QNUUPDATE.eq.2 ) then
            gpy = abs(gpy)
         endif
         if( QDAMP.eq.0 ) then
            wpz = DDOT(NIND, PZ, 1, WCORR, 1)
            wpz = dmin1( 0.d0, wpz )
         else
            wpz = 0.d0
         endif
         if( QPRINT.ge.0 ) then
            write(line,1004) wpz, gpy
 1004       format(/,' In update_nu: wpz = ',d12.5,' gpy = ',d12.5)
            call C_OUT(1,0,2,line)
         endif

         if( (1.d0-QRHO)*NU*(CNRM0**2) .lt.
     1        (gpy-wpz)*CNRM0 ) then
            NU = (gpy-wpz)*CNRM0/((1.d0-QRHO)*(CNRM0**2))
            if( WFLAG.gt.1 ) then
               write(line,*)
     1     'Warning in update_nu, watchdog active and NU increased to',
     2              NU
               call C_OUT(2,0,1,line)
c     WFLAG = 1
            endif
C            NU = dmax1(NU, QNUMIN)
            NU = NU + QNUMIN
CTODO zero or some small number?
         elseif( CNRM0.ne.0.d0 .and. ERR.ge.QLS_SAFE .and.
     1           WFLAG.lt.1 ) then
            nunew = dmax1(0.d0,((gpy-wpz)*CNRM0)/
     1           ((1.d0-QRHO)*(CNRM0**2)))
            NU = dmax1((3.d0*NU+nunew)/4.d0,QNUMIN)
         endif

      elseif( QNUUPDATE.eq.3 ) then

         if( QDAMP.eq.0 .and. QFULL.ne.1 ) then
            call C_OUT(2,0,1,
     1           'update_nu: Need damping for multiplier-based update.')
            IERR = 4
            goto 9999
         endif
         if( QLAMBDA.ne.1 .and. QFULL.ne.1 ) then
            call C_OUT(2,0,1,
     1           'update_nu: Need first order multipliers.')
            IERR = 4
            goto 9999
         endif
C
C     Get dual norm of LAMBDAs
C
         QCNRM = 4 - QCNRM
         nunew = CALC_NRM(M, LAM) + QNUMIN
         if( NU.lt.nunew ) then
            NU = nunew + QNUMIN
         endif
         QCNRM = 4 - QCNRM

      elseif( QNUUPDATE.eq.4 ) then
C
C     KNITRO-DIRECT rule
C
         if( CNRM0.gt.0.d0 ) then
            gd = DDOT(N, GB, 1, DX, 1)
            if( WCORR(1).gt.0.d0 ) then
               sigma = 0.5d0
            else
               sigma = 0.d0
            endif
            nu_trial = (gd+sigma*WCORR(1))/((1.d0-QRHO)*CNRM0)
            if( ERR.ge.QLS_SAFE .or. NU.lt.nu_trial ) then
               NU = max(QNUMIN,nu_trial + 1.d0)
            endif
         endif
      else
         call C_OUT(2,0,1,'update_nu: Invalid choice of QNUUPDATE.')
         IERR = 4
         goto 9999
      endif

      if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
         write(line,700) ITER, NU
 700     format(/,'  Information about penalty parameter in ITER ',i6,
     1        //,' NU = ',d20.10)
         call C_OUT(1,3,4,line)
      endif
      goto 9999
C
C     INDIVIDUAL PENALTY PARAMETERS FOR EACH CONSTRAINT
C
 2000 continue
C
C     check if everything OK
C
      if( QMERIT.lt.-2 .or. QLAMBDA.eq.0 .or. QCNRM.ne.1 ) then
         call C_OUT(2,0,1,'update_nu: some parameters not OK!')
         IERR = 4
         goto 9999
      endif
      cnrm1 = DASUM(M, C, 1)
C
      if( QFULL.eq.0 ) then
C
C     TODO: CLEAN UP THIS MESS!!!
C
C
C     Reduced space approach
C     (only for QLAMBDA = 1)
C
         if( QLAMBDA.ne.1 ) then
            call C_OUT(2,0,1,
     1        'update_nu: for indiv NUS in red space need QLAMBDA=1!')
            IERR = 4
            goto 9999
         endif
C
         if( QDAMP.eq.0 ) then
            delta = DDOT(NIND, WCORR, 1, PZ, 1)
         else
            delta = 0.d0
         endif
C
CTODO zero or some small number?
         if( cnrm1.gt.0d0 ) then
            delta = delta/cnrm1
         endif
CTODO decide if this is better:
         delta = dmin1(0d0, delta) ! this seems to perform better, but in many
                                ! cases, delta = 0 is even better! :(
C      delta = 0d0
C
         if( QPRINT.gt.5 ) then
            write(line,2001) delta
 2001       format(/,' delta in update_nu = ',d12.5)
            call C_OUT(1,6,2,line)
         endif
C
C     update the penatly parameters
C
CTODO update only for cnrm > 0 ?
CTODO zero or some small number?
         if( cnrm1.gt.0d0 ) then
            if( ERR.gt.QLS_SAFE ) then
               do i = 1, M
CTODO could try max(0, lam_j*c_j/|c_j|) instead of |lam_j|
CORIG               NUS(i) = dmax1((dabs(LAM(i))-delta)/(1d0-QRHO),QNUMIN)
                  cj = C(i)
CTODO zero or small number?
                  if ( cj.ne.0d0 ) then
CTODO decide which version
                     lc = LAM(i)*cj/dabs(cj)
C                  lc = dabs(LAM(i))
                  else
                     lc = 0d0
                  endif
                  NUS(i) = dmax1((lc-delta)/(1d0-QRHO),QNUMIN)
               enddo
            else
               do i = 1, M
CORIG               NUS(i) = dmax1((dabs(LAM(i))-delta)/(1d0-QRHO),QNUMIN,NUS(i))
                  cj = C(i)
CTODO zero or small number?
                  if ( cj.ne.0d0 ) then
CTODO decide which version
                     lc = LAM(i)*cj/dabs(cj)
C                  lc = dabs(LAM(i))
                  else
                     lc = 0d0
                  endif
                  NUS(i) = dmax1((lc-delta)/(1d0-QRHO),QNUMIN,NUS(i))
               enddo
            endif
         endif
      else
C
C     Full space approach
C
         if( cnrm1.gt.0d0 ) then
            if( ERR.gt.QLS_SAFE ) then
               do i = 1, M
                  NUS(i) = dmax1(dabs(LAM(i))/(1d0-QRHO),QNUMIN)
               enddo
            else
               do i = 1, M
                  NUS(i) = dmax1(NUS(i),dabs(LAM(i))/(1d0-QRHO),QNUMIN)
               enddo
            endif
         endif
      endif
C
C     return largest parameter in NU
C
      i = IDAMAX(M, NUS, 1)
      if( i.gt.0 ) then
         NU = NUS(i)
      else
         NU = 0d0
      endif
C
      if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
         write(line,1700) ITER, NU
 1700    format(/,'  Information about penalty parameter in ITER ',i6,
     1        'NU = ',d18.8)
         call C_OUT(1,3,2,line)
         call C_OUT(1,3,1,' ')
         do i = 1, M
            write(line,1701) i, NUS(i)
 1701       format(' NU(',i5,') = ',d20.10)
            call C_OUT(1,3,1,line)
         enddo
      endif
      goto 9999

 9999 continue

      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine UPDATE_NU_WS(N, M, NLB, NUB, NZA, LRW, LIW)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      LIW = 0
      LRW = 0
      return
      end

