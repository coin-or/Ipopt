C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine RESTO_KKT(ITER, N, NIND, M, NORIG, XORIG,
     1     CSCALE, IVAR, NFIX, IFIX, NLB, ILB, NUB,
     1     IUB, BNDS_L, BNDS_U, MU, X, CNRM0, NFILTER, FILTER_PHI,
     1     FILTER_C, CNRM_MAX, S_L, S_U, V_L, V_U, F, G, C,
     1     RESTO_KKT_FLAG, LAMOLD, LAM, DX, DV_L, DV_U,
     1     ALPHA, ALPHA_PRIMAL, ALPHA_DUAL,
     1     PRECFACT, MACHEPS, MACHTINY, NU_OUT, C_ACCEPT,
     1     KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR, EV_F,
     1     EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)

C
C*******************************************************************************
C
C    $Id: resto_kkt.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Restoration phase for filter - only trying to use usual steps to
CT    decrease primal and dual infeasibility
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
CP   MU        I    DP     barrier parameter
CP   X        I/O   DP     actual iterate (reordered without fixed vars:
CP                             first M entries belong to dependent
CP                             variables, remaining to independent variables)
CP   CNRM0    I/O   DP     norm of constraint violation
CP   NFILTER  I/O   INT    number of filter entries (corners)
CP   FILTER_PHI I/O DP     filter coordinates; phi-axis
CP   FILTER_C I/O   DP     filter coordinates; feasibility axis
CP   CNRM_MAX  I    DP     maximal value of CNRM for filter
CP   S_L      I/O   DP     slacks to lower bounds
CP   S_U      I/O   DP     slacks to upper bounds
CP   V_L      I/O   DP     dual variables for lower bounds (not really needed)
CP   V_U      I/O   DP     dual variables for upper bounds (not really needed)
CP   F        I/O   DP     value of objective function at X
CP   G         I    DP     gradient of objective function at X (not changed)
CP   C        I/O   DP     values of constraints at X
CP   RESTO_KKT_FLAG I/O INT =0: This is first call for this restoration phase
CP                          =1: We are in KKT restoration phase
CP                          =2: KKT restoration phase not successfull, continue
CP                              with other kind of restoration phase
CP   LAMOLD   I/O   DP     "Old" Value of equality multipliers
CP   LAM      I/O   DP     "New" Value of equality multipliers
CP   DX        I    DP     Step for primal variables
CP   DV_L      I    DP     Step for dual variables (lower bounds)
CP   DV_U      I    DP     Step for dual variables (upper bounds)
CP   ALPHA     O    DP     Step size taken
CP   ALPHA_PRIMAL I DP     Step size for primal variables after
CP                           fraction-to-the-boundary rule
CP   ALPHA_DUAL I   DP     Step size for dual variables after
CP                           fraction-to-the-boundary rule
CP   PRECFACT  I    DP     recision factor for roundoffs
CP   MACHEPS   I    DP     machine precision
CP   MACHTINY  I    DP     smallest DP number
CP   NU_OUT    O    DP     for iter_out output
CP   C_ACCEPT  O    C*1    for iter_out output
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
CS
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
      double precision X(N)
      double precision CNRM0
      integer NFILTER
      double precision FILTER_PHI(NFILTER)
      double precision FILTER_C(NFILTER)
      double precision CNRM_MAX
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision F
      double precision G(N)
      double precision C(M)
      integer RESTO_KKT_FLAG
      double precision LAMOLD(M)
      double precision LAM(M)
      double precision DX(N)
      double precision DV_L(NLB)
      double precision DV_U(NUB)
      double precision ALPHA
      double precision ALPHA_PRIMAL
      double precision ALPHA_DUAL
      double precision PRECFACT
      double precision MACHEPS
      double precision MACHTINY
      double precision NU_OUT
      character*1 C_ACCEPT
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
      integer p_rwend, p_iwend, p_dinf, p_gnew, i, j, idummy, irej
      integer p_xnew, p_cnew, p_slnew, p_sunew, p_vlnew, p_vunew
      integer p_lamnew
      double precision tmp(3), dinf, f_new, cnrmnew, dinfnew
      double precision kktnrmnew, bar_n, phi_new, dummy

      character*120 line
      double precision CALC_BAR, CALC_NRM, DDOT
      integer FFINITE, FILTER_CHECK

      double precision KKTNRM
      save             KKTNRM

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

      if( MEMDBG ) then
         write(line,1)'resto_ktt', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      if( abs(QLAMBDA).ne.2 ) then
         write(line,*)
     1        'resto_kkt: At this point, resto_kkt only works with',
     1        ' ILAMBDA = 2.'
         call C_OUT(2,0,1,line)
         IERR = 4
         goto 9999
      endif

      if( RESTO_KKT_FLAG.eq.0 ) then
C
C     This is the first of the restoration iterations, so we need to compute
C     current dual infeasibility
C
         p_dinf  = p_rwend
         p_rwend = p_dinf + N
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
C     Compute A * LAM
C
         call CONSTR(8, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, LAMOLD, RW(p_dinf+1),
     2        idummy, idummy,
     3        KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4        IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.lt.0 ) then
            write(line,*)
     1           'resto_kkt: Warning in CONSTR-8(1), IERR = ', IERR
            call C_OUT(2,0,1,line)
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'resto_kkt: Error in CONSTR-8(1), IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         call DAXPY(N, 1d0, G, 1, RW(p_dinf+1), 1)
C     Add gradient of barrier function and compute complementarity
         tmp(1) = 0.d0
         do i = 1, NLB
            j = ILB(i)
            RW(p_dinf+j) = RW(p_dinf+j) - V_L(i)
            tmp(1) = tmp(1) + (S_L(i)*V_L(i)-MU)**2
         enddo
         do i = 1, NUB
            j = IUB(i)
            RW(p_dinf+j) = RW(p_dinf+j) + V_U(i)
            tmp(1) = tmp(1) + (S_U(i)*V_U(i)-MU)**2
         enddo
         dinf = CALC_NRM(N, RW(p_dinf+1))
         p_rwend = p_dinf
C
         tmp(1) = dsqrt(tmp(1))
         tmp(2) = CNRM0
         tmp(3) = dinf
         KKTNRM = CALC_NRM(3, tmp)
CTODO    DEcide which norm!
      endif
C
C     Compute trial step size (for now simply same step sizes for primal and
C     dual variables)
C
      ALPHA = dmin1(ALPHA_PRIMAL, ALPHA_DUAL)
C
C     Compute value of variables at trial point
C
      p_xnew   = p_rwend
      p_lamnew = p_xnew   + N
      p_slnew  = p_lamnew + M
      p_sunew  = p_slnew  + NLB
      p_vlnew  = p_sunew  + NUB
      p_vunew  = p_vlnew  + NLB
      p_cnew   = p_vunew  + NUB
      p_rwend  = p_cnew   + M
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
      call DCOPY(N, X, 1, RW(p_xnew+1), 1)
      call DAXPY(N, ALPHA, DX, 1, RW(p_xnew+1), 1)
C
      call DCOPY(M, LAM, 1, RW(p_lamnew+1), 1)
      call DSCAL(M, ALPHA, RW(p_lamnew+1), 1)
      call DAXPY(M, (1.d0-ALPHA), LAMOLD, 1, RW(p_lamnew+1), 1)
C
C     Compute new trial slack variables
C     (If slacks get too small, the problem probably doesn't have a
C      feasible strictly interior point)
C
C     A sophisticated user might want to use a different way to handle this
C     situation.  In that case, the flag '-DUSE_IPMYMOVEBNDS' has to be
C     used at compilation time.
C
#ifdef USE_IPMYMOVEBNDS
      call IPMYMOVEBNDS(N, RW(p_xnew+1), NLB, RW(p_slnew+1), BNDS_L,
     1     ILB, V_L, DV_L, NUB, RW(p_sunew+1), BNDS_U, IUB, V_U, DV_U, 
     1     MU, ITER, PRECFACT, MACHEPS)
#else
      call CHECK_SLACKS(N, RW(p_xnew+1), NLB, ILB, S_L, V_L, BNDS_L,
     1     RW(p_slnew+1), NUB, IUB, S_U, V_U, BNDS_U, RW(p_sunew+1),
     1     MU, MACHEPS)
#endif
C
C     ...and new dual variables
C
      call DCOPY(NLB, V_L, 1, RW(p_vlnew+1), 1)
      call DAXPY(NLB, ALPHA, DV_L, 1, RW(p_vlnew+1), 1)
      call DCOPY(NUB, V_U, 1, RW(p_vunew+1), 1)
      call DAXPY(NUB, ALPHA, DV_U, 1, RW(p_vunew+1), 1)
C
C     Make sure that each V_NEW is at least MACHTINY
C
      do i = 1, NLB
         RW(p_vlnew+i) = dmax1(MACHTINY, RW(p_vlnew+i))
      enddo
      do i = 1, NUB
         RW(p_vunew+i) = dmax1(MACHTINY, RW(p_vunew+i))
      enddo
C
C     Compute value of F and C (to make sure neither is NaN)
C
      call GET_F(N, RW(p_xnew+1), IVAR, NORIG, XORIG, M, CSCALE, NLB,
     1     ILB, RW(p_slnew+1), NUB, IUB, RW(p_sunew+1), MU, f_new,
     1     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'resto_kkt: get_f returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( FFINITE(f_new).eq.0 ) then ! we don't want NaN here; cut step!
         call C_OUT(2,0,1,' --- WARNING: Evaluation of'//
     1        ' objective function produced Nan/Inf.')
         goto 5000
      endif
C
      call GET_C(ITER, N, NIND, RW(p_xnew+1), IVAR, NORIG, XORIG,
     1     M, CSCALE, RW(p_cnew+1), KCONSTR, LRS, RS, LIS, IS,
     2     LRW-p_rwend, RW(p_rwend+1),
     4     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1        'resto_kkt: Error: get_c returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.ne.0 ) then
         write(line,*)
     1        'resto_kkt: Warning: get_c returns IERR = ',IERR
         call C_OUT(2,0,1,line)
      endif
      cnrmnew = CALC_NRM(M, RW(p_cnew+1))
      if( FFINITE(cnrmnew).eq.0 ) then ! we don't want NaN here; cut step!
         call C_OUT(2,0,1,' --- WARNING: Evaluation of constraint'//
     1        ' functions produced Nan/Inf.')
         goto 5000
      endif
C
C     Compute dual infeasibiliy
C
      p_dinf  = p_rwend
      p_gnew  = p_dinf + N
      p_rwend = p_gnew + N
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call GET_G(N, RW(p_xnew+1), IVAR, NORIG, XORIG, M, CSCALE, NLB,
     1     ILB, S_L, NUB, IUB, S_U, MU, RW(p_gnew+1), LRW-p_rwend,
     2     RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR, EV_G, DAT,
     3     IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'resto_kkt: GET_G returned IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C     Compute A * LAM
C
C     (Note: Call CONSTR with ITER+1, so that Jacobian will be reused in next
C      iteration if step is accepted - if not accepted, CONSTR will be called
C      next time with ITER again, and the old Jacobian will be recomputed)
C     NEED TO ENSURE THAT XORIG THEN IS SET ACCORDING TO X (not XNEW) AGAIN!
      call CONSTR(8, ITER+1, N, NIND, M, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, RW(p_lamnew+1), RW(p_dinf+1),
     2     idummy, idummy,
     3     KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4     IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.lt.0 ) then
         write(line,*)
     1        'resto_kkt: Warning in CONSTR-8, IERR = ', IERR
         call C_OUT(2,0,1,line)
      elseif( IERR.ne.0 ) then
         write(line,*)
     1        'resto_kkt: Error in CONSTR-8, IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      call DAXPY(N, 1d0, RW(p_gnew+1), 1, RW(p_dinf+1), 1)
C     Add gradient of barrier function
C     NOTE: Here, we takes the dual variables and not the real gradient
C     or the log barrier terms, since for latter changes for changing MU!
      do i = 1, NLB
         j = ILB(i)
         RW(p_dinf+j) = RW(p_dinf+j) - RW(p_vlnew+i)
      enddo
      do i = 1, NUB
         j = IUB(i)
         RW(p_dinf+j) = RW(p_dinf+j) + RW(p_vunew+i)
      enddo
C
CTODO What kind of norm????
C
      dinfnew = CALC_NRM(N, RW(p_dinf+1))
      p_rwend = p_dinf
C
C     Compute new complementarity
C
      tmp(1) = 0.d0
      do i = 1, NLB
         tmp(1) = tmp(1) + (RW(p_slnew+i)*RW(p_vlnew+i)-MU)**2
      enddo
      do i = 1, NUB
         tmp(1) = tmp(1) + (RW(p_sunew+i)*RW(p_vunew+i)-MU)**2
      enddo
      tmp(1) = dsqrt(tmp(1))
C
C     Check if sufficient progress is obtained
C
      tmp(2)    = cnrmnew
      tmp(3)    = dinfnew
      kktnrmnew = CALC_NRM(3, tmp)

      if( QCNR.gt.0 ) then
         write(line,500) kktnrmnew, KKTNRM
 500     format('resto_kkt: kktnrmnew = ',d15.7,' KKTNRM    = ',d15.7)
         call C_OUT(1,1,1,line)
      endif
      if( kktnrmnew.gt.QRESTOKKTRED*KKTNRM ) goto 5000 ! not suff. progress made
C
C     Ok, sufficient progress is made, now see if we are overall done with
C     the restoration phase
C
C     Check if constraint violation is 'exploding'
C
      if( cnrmnew.gt.CNRM_MAX ) then
         if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
            write(line,*) 'resto_kkt: cnrm = ',cnrmnew,
     1           ' is getting too large.'
            call C_OUT(1,2,1,line)
         endif
         goto 1000
      endif
C
      bar_n   = CALC_BAR(NLB, NUB, RW(p_slnew+i), RW(p_sunew+i),
     1     dummy, dummy, MU)
      phi_new = f_new + bar_n
      if( abs(QMERIT).eq.7 ) then
         phi_new = phi_new + DDOT(M, RW(p_cnew+1), 1, RW(p_lamnew+1), 1)
      endif
      if( FFINITE(phi_new).eq.0 ) then ! we don't want NaN here; cut step!
         call C_OUT(2,0,1,' --- WARNING: Evaluation of barrier'//
     1        ' objective function produced Nan/Inf.')
         goto 1000
      endif
C
C     Check if trial point is acceptable to the filter
C
      irej = FILTER_CHECK(NFILTER, FILTER_C, FILTER_PHI,
     1     cnrmnew, phi_new, PRECFACT, MACHEPS)
      if( irej.gt.0 ) goto 1000
C
C     OK, point is accepted!
C
      RESTO_KKT_FLAG = 0
      goto 2000

 1000 continue
C
C     Point is not acceptable to filter, but we want to continue with this
C     type of restoration phase
C
      RESTO_KKT_FLAG = 1
      KKTNRM = kktnrmnew        ! store for next iteration

 2000 continue
C
C     Accept new point
C
      F = f_new
      call DCOPY(N  , RW(p_xnew+1)  , 1, X     , 1)
      call DCOPY(M  , RW(p_cnew+1)  , 1, C     , 1)
      call DCOPY(NLB, RW(p_slnew+1) , 1, S_L   , 1)
      call DCOPY(NUB, RW(p_sunew+1) , 1, S_U   , 1)
      call DCOPY(NLB, RW(p_vlnew+1) , 1, V_L   , 1)
      call DCOPY(NUB, RW(p_vunew+1) , 1, V_U   , 1)
      call DCOPY(M  , RW(p_lamnew+1), 1, LAM   , 1)
      call DCOPY(M  , LAM           , 1, LAMOLD, 1)
      p_rwend = p_xnew
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
      NU_OUT   = 0d0
      C_ACCEPT = 'K'
C
      goto 9999
C
C     Point has not been accepted - continue with usual restoration phase
C
 5000 continue
      RESTO_KKT_FLAG = 2

 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine RESTO_KKT_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1, lrw2, liw2, lrw3, liw3, lrw_constr, liw_constr
      character*80 line

      call CONSTR_WS(N, M, NLB, NUB, NZA, lrw_constr, liw_constr,
     1     DAT, IDAT)
      LIW = liw_constr
      LRW = lrw_constr + N

      call GET_C_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
      call GET_F_WS(N, M, NLB, NUB, NZA, lrw3, liw3)
      call GET_G_WS(N, M, NLB, NUB, NZA, lrw2, liw2)
      lrw1 = max(lrw1, lrw3, 2*N + max(lrw2,lrw_constr))
      liw1 = max(liw1, liw3, max(liw2,liw_constr))

      LRW = max(LRW, 2*(M+NLB+NUB) + N + lrw1)
      LIW = max(LIW, liw1)

      if( QPRINT.ge.4 ) then
         write(line,1000)'resto_kkt_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
