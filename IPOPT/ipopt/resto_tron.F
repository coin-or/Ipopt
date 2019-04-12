C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine RESTO_TRON(ITER, N, NIND, M, NORIG, XORIG,
     1     CSCALE, IVAR, NFIX, IFIX, NLB, ILB, NUB,
     1     IUB, BNDS_L, BNDS_U, MU, X, CNRM, NFILTER, FILTER_PHI,
     1     FILTER_C, CNRM_MAX, S_L,
     1     S_U, V_L, V_U, F, C, PRECFACT, MACHEPS, MACHTINY,
     1     XHESS, CHESS, LAM, KCONSTR, LRS, RS, LIS, IS,
     1     LRW, RW, LIW, IW, IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV,
     2     EV_HOV, EV_HCV, DAT, IDAT)

C
C*******************************************************************************
C
C    $Id: resto_tron.F 574 2004-04-25 22:42:37Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Restoration phase for filter using TRON
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
CP   CNRM      O    DP     norm of constraint violation
CP   NFILTER   I    INT    number of filter entries (corners)
CP                         for ITER = -1: NZORIG
CP   FILTER_PHI I   DP     filter coordinates; phi-axis
CP   FILTER_C  I    DP     filter coordinates; feasibility axis
CP   CNRM_MAX  I    DP     maximal value of CNRM for filter
CP   S_L      I/O   DP     slacks to lower bounds
CP   S_U      I/O   DP     slacks to upper bounds
CP   V_L       I    DP     dual variables for lower bounds (not really needed)
CP   V_U       I    DP     dual variables for upper bounds (not really needed)
CP   F        I/O   DP     value of objective function at X
CP   C        I/O   DP     values of constraints at X
CP   PRECFACT  I    DP     precision factor for roundoffs
CP   MACHEPS   I    DP     machine precision
CP   MACHTINY  I    DP     smallest DP number
CP   XHESS     W    DP     work space (X at point where Hessian is evaluated)
CP   CHESS     W    DP     work space (C at point where Hessian is evaluated)
CP   LAM       I    DP     current values of multipliers (for Lagrangian)
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
CS    CONSTR
CS    CALC_BAR
CS    CALC_NRM
CS    DTRON
CS    DGPNRM2
CS    GET_C
CS    GET_F
CS    FFINITE
CS    C_OU
CS    MC35AD
CS    MC39AD
CS    MC39BD
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
      double precision PRECFACT
      double precision MACHEPS
      double precision MACHTINY
      double precision XHESS(N)
      double precision CHESS(M)
      double precision LAM(M)
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
      integer p_xl, p_xu, p_g, p_hdiag, p_bdiag, p_ldiag, p_xc, p_s
      integer p_dsave, p_wa, p_h, p_l, p_rwend, p_hcol_ptr, p_iptr
      integer p_bcol_ptr, p_lcol_ptr, p_indfree, p_isave, p_iwa, p_iw
      integer p_hrow_ind, p_brow_ind, p_lrow_ind, p_iwend, p_a, p_acon
      integer p_b, p_zeros, p_hf, p_chc, p_chc_col, p_chcrow_ind
      integer p_chccol_ptr, p_iwork
      integer itermax, i, k, ltron_iter, ncol, noff, nzata_out, idummy
      integer lchc, lh, l_ata, nchccol, natacol, nhcol, rowind_ata
      integer rowind_chc, count_rej, iterstart, cgcount_start

      double precision vin(2), fatol, frtol, fmin, cgtol, cnrm0
      double precision gnorm, gnorm0, lb, ub, bar_n, phi_new, cnrm2
      double precision fi, ci, lhs, dummy

      logical soflag, relax_bounds

      character*120 line(4)
      character*60 task
      character*1 TRflag
      double precision DNRM2, CALC_BAR, CALC_NRM, DGPNRM2, DDOT
      integer FFINITE

      integer NNZA, NNZATA, LATA, NNZH, NNZCHC
      save    NNZA, NNZATA, LATA, NNZH, NNZCHC

      double precision HUGE
      parameter       (HUGE = 1d300)

      double precision DELTA
      save             DELTA

#ifdef HAVE_MC35
      integer LP_MC35
      common /MC35CD/ LP_MC35
#endif
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

      if( ITER.eq.-1 ) then
         if( QTRONHESS.ne.0 ) then
            NNZA = NFILTER
            NNZATA = 3*NNZA+N
            LATA = max(NNZATA, NNZA+N)
            if( QTRON2DERIVS.eq.0 ) then
               NNZCHC = 0
            else
               call GET_H(-1, N, NIND, idummy, dummy, idummy, NORIG,
     1              XORIG, dummy, M, NNZCHC, dummy, dummy, idummy,
     1              idummy, idummy, dummy, idummy, idummy, idummy,
     1              EV_H, DAT, IDAT)
            endif
            NNZH = NNZATA + NNZCHC !This is a very generous estimate!
         endif
         DELTA = -1d0
         goto 9999
      endif

      if( MEMDBG ) then
         write(line,1) 'resto_tron', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif
C
      if( IERR.eq.17 .or. IERR.eq.519 ) then
         relax_bounds = .false.
      else
         relax_bounds = .true.
      endif
      IERR = 0
      p_rwend = 0
      p_iwend = 0
      p_l = 0
      p_b = 0
      p_brow_ind = 0
      p_lrow_ind = 0
C
C     Get work space for TRON
C
      p_xl    = p_rwend
      p_xu    = p_xl + N
      p_g     = p_xu + N
      p_xc    = p_g + N
      p_s     = p_xc + N
      p_dsave = p_s + N
      p_wa    = p_dsave + N
      if( QTRONHESS.eq.0 ) then
         p_rwend = p_wa + 7*N
         p_hdiag = 0
         p_bdiag = 0
         p_ldiag = 0
         p_h     = 0
      else
         p_hdiag = p_wa + 7*N
         p_bdiag = p_hdiag + N
         p_ldiag = p_bdiag + N
         p_h     = p_ldiag + N
         p_rwend = p_h + max(LATA, NNZH)
      endif
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif

      p_indfree  = p_iwend
      p_isave    = p_indfree + N
      p_iwa      = p_isave + N
      if( QTRONHESS.eq.0 ) then
         p_iwend    = p_iwa + 3*N
         p_hcol_ptr = 0
         p_bcol_ptr = 0
         p_lcol_ptr = 0
         p_hrow_ind = 0
      else
         p_hcol_ptr = p_iwa + 3*N
         p_bcol_ptr = p_hcol_ptr + N+1
         p_lcol_ptr = p_bcol_ptr + N+1
         p_hrow_ind = p_lcol_ptr + N+1
         p_iwend    = p_hrow_ind + max(LATA, NNZH)
      endif
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif

c     Set parameters

      itermax = N
c      itermax = 10
      fatol=0.0d0
      frtol=0.0d0
      fmin =0.0d0
      cgtol=QTRONCGTOL

c     We now specify the bounds

      call DCOPY(N, -HUGE, 0, RW(p_xl+1), 1)
      call DCOPY(N,  HUGE, 0, RW(p_xu+1), 1)
C
C     Bounds are set accoring to fraction-to-boundary rule and minimal distance
C     QTHETA
C
      do i = 1, NLB
         k = ILB(i)
         lb = BNDS_L(i)
         if( relax_bounds) then
CORIG         RW(p_xl+k) = dmin1(lb + QTHETA, lb + (1-QTAU)*(X(k)-lb))
            RW(p_xl+k) = dmin1(lb + QTHETA, X(k))
CORIG         RW(p_xl+k) = dmax1(RW(p_xl+k), lb + 1d-10)
         else
            RW(p_xl+k) = lb
         endif
C         X(k) = dmax1(X(k), RW(p_xl+k))
      enddo

      do i = 1, NUB
         k = IUB(i)
         ub = BNDS_U(i)
         if( relax_bounds) then
CORIG         RW(p_xu+k) = dmax1(ub - QTHETA, ub - (1-QTAU)*(ub-X(k)))
            RW(p_xu+k) = dmax1(ub - QTHETA, X(k))
CORIG         RW(p_xu+k) = dmin1(RW(p_xu+k), ub - 1d-10)
         else
            RW(p_xu+k) = ub
         endif
C         X(k) = dmin1(X(k),RW(p_xu+k))
      enddo

CC
CC     The following is for Arvind's new complementarity approach
CC
C      do i = 1, NCOMPL
C         k = ICOMPL(i)
C         RW(p_xl+k) = dmin1(QTHETA, X(k))
C      enddo
C
C     Some output
C
      call C_OUT(2,0,1,'---- Entering restoration phase with TRON ----')

c     We start the iteration by initializing task.
c
      ltron_iter = 0
      task = 'START'
      iterstart = ITER
      cgcount_start = COUNT_TRON_CG

c        ------- the beginning of the loop ----------

 111  continue
      ltron_iter = ltron_iter + 1

C     Function evaluation

      if( task.eq.'F' .or. task.eq.'START' ) then

         call GET_C(-ltron_iter, N, NIND, X, IVAR, NORIG, XORIG, M,
     1        CSCALE, C, KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'resto_tron-f: Error: get_c returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'resto_tron-f: Warning: get_c returns IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
         cnrm2 = DNRM2(M, C, 1)
         if( task.eq.'START' ) cnrm0 = cnrm2

      endif

C     Evaluate grandient and Gauss-Newton Hessian

      if( task.eq.'GH' .or. task.eq.'START' ) then

         if( QTRONHESS.eq.0 ) then
            call DCOPY(N, X, 1, XHESS, 1)
         endif
CTODO figure out if this is needed
         if( task.ne.'START' ) then
            call GET_C(-ltron_iter, N, NIND, X, IVAR, NORIG, XORIG, M,
     1           CSCALE, C, KCONSTR, LRS, RS, LIS, IS,
     2           LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.gt.0 ) then
               write(line,*)
     1              'resto_tron-gh: Error: get_c returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            elseif( IERR.ne.0 ) then
               write(line,*)
     1             'resto_tron-gh: Warning: get_c returns IERR = ',IERR
               call C_OUT(2,0,1,line)
            endif
         endif

c        Compute gradient g = A c

         call CONSTR(8, -ltron_iter, N, NIND, M, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, C, RW(p_g+1), idummy, idummy,
     3        KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4        IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.lt.0 ) then
            write(line,*)
     1           'resto_tron: Warning in CONSTR-8, IERR = ',IERR
            call C_OUT(2,0,1,line)
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'resto_tron: Error in CONSTR-8, IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif

C     Get Jacobian

C
C     If the Hessian is not to be computed explicitely, skip the following part
C
         if( QTRONHESS.eq.0 ) goto 5775
 120     continue
         p_a     = p_rwend
         p_rwend = p_a + NNZA
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         p_acon  = p_iwend
         p_iwend = p_acon + NNZA
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         call CONSTR(10, -ltron_iter, N, NIND, M, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, vin, RW(p_a+1), IW(p_acon+1),
     2        IW(p_hrow_ind+1),
     3        KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4        IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.lt.0 ) then
            write(line,*)
     1           'resto_tron: Warning in CONSTR-10, IERR = ',IERR
            call C_OUT(2,0,1,line)
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'resto_tron: Error in CONSTR-10, IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         NNZA = int(vin(1))
CC
CC     Compute estimate of condition number
CC
C         ljcn    = 10*NNZA
C         lirn    = 10*NNZA
C         p_aa    = p_rwend
C         p_w     = p_aa + ljcn
C         p_rwend = p_w + 4*N
C         if( p_rwend.gt.LRW ) then
C            IERR = 98
C            goto 9999
C         endif
C         p_jcn   = p_iwend
C         p_irn   = p_jcn + ljcn
C         p_iw    = p_irn + lirn
C         p_iwend = p_iw  + 13*N
C         if( p_iwend.gt.LIW ) then
C            IERR = 99
C            goto 9999
C         endif
C         call DCOPY(NNZA, RW(p_a+1), 1, RW(p_aa+1), 1)
C         do i = 1, NNZA
C            IW(p_irn+i) = IW(p_acon+i)
C            IW(p_jcn+i) = IW(p_hrow_ind+i)
C         enddo
C         call MC45AD(N, NNZA, ljcn, IW(p_jcn+1), lirn, IW(p_irn+1),
C     1        RW(p_aa+1), cond, IW(p_iw+1), RW(p_w+1), iflag)
C         write(48,*) 'ltron_iter = ',ltron_iter
C         write(48,*) 'iflag = ', iflag
C         write(48,*) 'irank = ',irank,', M = ',M
C         write(48,*) 'cond  = ',cond
C         p_rwend = p_aa
C         p_iwend = p_jcn

C     Compute A^T A

         p_iptr  = p_iwend
         p_iw    = p_iptr + M+1
         p_iwend = p_iw + M+N
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
#ifndef HAVE_MC35
         IERR = 97
         call C_OUT(2,0,1,
     1        'resto_tron:  Compiled without MC35.')
         call C_OUT(2,0,1,
     1        '             Choose different value for ITRONHESS.')
         goto 9999
#else
         LP_MC35 = QCNR
         call MC35AD(M, N, 3, NNZA, RW(p_a+1), IW(p_acon+1), LATA,
     1        IW(p_hrow_ind+1), RW(p_h+1), IW(p_hcol_ptr+1),
     2        IW(p_iptr+1), IERR, nzata_out, IW(p_iw+1) )
#endif
         if( IERR.eq.3 .or. IERR.eq.4 ) then
C     NNZATA was too small; correct work space sizes and re-evaluate Jacobian
            write(line,*) 'resto_tron: increase LATA from ',LATA,
     1           'to ',nzata_out
            call C_OUT(2,0,1,line)
            LATA    = max(nzata_out, NNZA+N)
            NNZATA  = nzata_out
            NNZH = max(NNZH, NNZATA + NNZA)
            p_rwend = p_h + max(LATA, NNZH)
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            p_iwend = p_hrow_ind + max(LATA, NNZH)
            if( p_iwend.gt.LIW ) then
               IERR = 99
               goto 9999
            endif
            goto 120
         elseif( IERR.ne.0 ) then
            write(line,*) 'resto_tron: MC35AD returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 536
            goto 9999
         endif
         NNZATA = nzata_out
         p_rwend = p_a
         p_iwend = p_acon

         if( QTRON2DERIVS.eq.0 ) then
            NNZH = NNZATA
            goto 160
         endif

C
C     Order rows in A^T A
C
CTODO Find out if really necessary
         p_iw    = p_iwend
         p_iwork = p_iw + N+1
         p_iwend = p_iwork + NNZH
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
#ifndef HAVE_MC39
         IERR = 97
         call C_OUT(2,0,1,
     1        'resto_tron:  Compiled without MC39.')
         call C_OUT(2,0,1,
     1        '             Choose different value for ITRONHESS.')
         goto 9999
#else
         call MC39BD( N, N, NNZATA, IW(p_hrow_ind+1), .true., RW(p_h+1),
     1        IW(p_hcol_ptr+1), IW(p_iwork+1), IW(p_iw+1))
#endif
         p_iwend = p_iw
C
C     Reserve space for Hessian calculations
C
         p_chc   = p_rwend
         p_hf    = p_chc + NNZCHC
         p_rwend = p_hf  + NNZCHC
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         p_chcrow_ind = p_iwend
         p_chccol_ptr = p_chcrow_ind + NNZCHC
         p_chc_col    = p_chccol_ptr + N+1
         p_iwend      = p_chc_col + NNZCHC
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
C
C     Compute sum (c_i * na^2 c_i )
C     Substract hessian of objective function from Hessian of Lagrangian
C     with the values of the constraints as multipliers
C
         p_zeros = p_rwend
         p_rwend = p_zeros + M
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         call DCOPY(M, 0d0, 0, RW(p_zeros+1), 1)
         call GET_H(ITER, N, NIND, NFIX, X, IVAR, NORIG, XORIG, CSCALE,
     1        M, NNZCHC, RW(p_zeros+1), RW(p_hf+1), IW(p_chc_col+1),
     1        IW(p_chcrow_ind+1), LRW-p_rwend, RW(p_rwend+1),
     1        LIW-p_iwend, IW(p_iwend+1), IERR, EV_H, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*) 'resto_tron-f: GET_H returns IERR =',IERR
            call C_OUT(2,0,1,line)
         endif
         p_rwend = p_zeros
         call GET_H(ITER, N, NIND, NFIX, X, IVAR, NORIG, XORIG, CSCALE,
     1        M, NNZCHC, C, RW(p_chc+1), IW(p_chc_col+1),
     1        IW(p_chcrow_ind+1), LRW-p_rwend, RW(p_rwend+1),
     1        LIW-p_iwend, IW(p_iwend+1), IERR, EV_H, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*) 'resto_tron-c: GET_H returns IERR =',IERR
            call C_OUT(2,0,1,line)
         endif
         call DAXPY(NNZCHC, -1d0, RW(p_hf+1), 1, RW(p_chc+1), 1)
         p_rwend = p_hf
C
C     Order Hessian into column format
C
         p_iw          = p_iwend
         p_iwend       = p_iw + N+1
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
#ifndef HAVE_MC39
         IERR = 97
         call C_OUT(2,0,1,
     1        'resto_tron:  Compiled without MC39.')
         call C_OUT(2,0,1,
     1        '             Choose different value for ITRONHESS.')
         goto 9999
#else
         call MC39AD(1, N, N, NNZCHC, IW(p_chcrow_ind+1),
     1        IW(p_chc_col+1),
     1        .true., RW(p_chc+1), IW(p_chccol_ptr+1), IW(p_iw+1))
#endif
         p_iwend = p_chc_col
C
C     Add A^T A and CHC
C
         lchc = NNZCHC
         lh   = NNZH
         l_ata = NNZATA
         nchccol = N
         nhcol   = N
         natacol = N
         IW(p_hcol_ptr+N+1) = NNZH+1
 350     continue
         if( l_ata.gt.0 .or. lchc.gt.0 ) then
            if( lh.lt.1 ) then
               call C_OUT(2,0,1,'resto_tron: NNZH too small! Abort.')
               IERR = 537
               goto 9999
            endif
 360        continue
            if( natacol.gt.1 ) then
               if( IW(p_hcol_ptr+natacol).gt.l_ata ) then
                  natacol = natacol - 1
                  goto 360
               endif
            endif
 365        continue
            if( nchccol.gt.1 ) then
               if( IW(p_chccol_ptr+nchccol).gt.lchc ) then
                  nchccol = nchccol - 1
                  goto 365
               endif
            endif
 370        continue
            if( nhcol.gt.max(natacol,nchccol) ) then
               IW(p_hcol_ptr+nhcol) = lh + 1
               nhcol = nhcol - 1
               if( nhcol.lt.1 ) then
                  WRITE(*,*) 'Shouldn''t get here...'
                  stop
               endif
               goto 370
            endif
            if( l_ata.gt.0 .and. nhcol.eq.natacol ) then
               rowind_ata = IW(p_hrow_ind+l_ata)
            else
               rowind_ata = -1
            endif
            if( lchc.gt.0 .and. nhcol.eq.nchccol ) then
               rowind_chc = IW(p_chcrow_ind+lchc)
            else
               rowind_chc = -1
            endif
            if( rowind_ata.gt.rowind_chc) then
               RW(p_h+lh)        = RW(p_h+l_ata)
               IW(p_hrow_ind+lh) = rowind_ata
               l_ata = l_ata - 1
            elseif( rowind_ata.lt.rowind_chc) then
               RW(p_h+lh)        = RW(p_chc+lchc)
               IW(p_hrow_ind+lh) = rowind_chc
               lchc = lchc - 1
            else
               RW(p_h+lh)        = RW(p_chc+lchc) + RW(p_h+l_ata)
               IW(p_hrow_ind+lh) = rowind_chc
               lchc = lchc - 1
               l_ata = l_ata - 1
            endif
            IF(lchc.lt.-1) then
               write(line,*) 'resto_tron: lchc = ',lchc,'. Abort.'
               stop
            endif
            lh   = lh   - 1
            goto 350
         endif
         IW(p_hcol_ptr+1) = lh + 1
C     If necessary correct NNZH and work space pointers
         if( lh.gt.0 ) then
            if( task.ne.'START' ) then
               WRITE(*,*) ' ERROR in RESTO_TRON!'
               stop
            endif
            NNZH = NNZH - lh
            do i = 1, NNZH
               RW(p_h + i) = RW(p_h + lh + i)
            enddo
            do i = 1, NNZH
               IW(p_hrow_ind + i) = IW(p_hrow_ind + lh + i)
            enddo
            do i = 1, N+1
               IW(p_hcol_ptr+i) = IW(p_hcol_ptr+i) - lh
            enddo
         endif

 160     continue
         if( task.eq.'START' ) then
            p_b     = p_h + max(NNZH, LATA)
            p_l     = p_b + NNZH
            p_rwend = p_l + NNZH + 5*N
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            p_brow_ind = p_hrow_ind + max(NNZH, LATA)
            p_lrow_ind = p_brow_ind + NNZH
            p_iwend    = p_lrow_ind + NNZH + 5*N
            if( p_iwend.gt.LIW ) then
               IERR = 99
               goto 9999
            endif
         elseif( QTRON2DERIVS.eq.1 ) then
            p_rwend = p_chc
            p_iwend = p_chcrow_ind
         endif

C     Extract diagonal matrix

         call DCOPY(N, 0d0, 0, RW(p_hdiag+1), 1)
         noff = 0
         ncol = 1
         do i = 1, NNZH
 130        continue
            if( IW(p_hcol_ptr+ncol+1).le.i ) then
               ncol = ncol + 1
               IW(p_hcol_ptr+ncol) = IW(p_hcol_ptr+ncol) - noff
               goto 130
            endif
            if( IW(p_hrow_ind+i).eq.ncol ) then
               RW(p_hdiag+ncol) = RW(p_h+i)
               noff = noff + 1
            elseif( noff.gt.0 ) then
               RW(p_h+i-noff)        = RW(p_h+i)
               IW(p_hrow_ind+i-noff) = IW(p_hrow_ind+i)
            endif
         enddo
         ncol = ncol + 1
         IW(p_hcol_ptr+ncol) = IW(p_hcol_ptr+ncol) - noff
 5775    continue
      endif

C     initialized the trust region bound

      if( task.eq.'START' .and. DELTA.lt.0d0 ) then
CTODO      if( task.eq.'START' ) then
         gnorm0 = DNRM2(N, RW(p_g+1), 1)
         DELTA  = 1d1*gnorm0
         DELTA  = 1d0
      endif

      if( DELTA.le.1d-30 ) then
         call C_OUT(2,0,1,'resto_tron: DELTA becomming too small.')
         IERR = 18
         goto 9999
      endif

C     stopping criteria

      if( ITER.ge.QMAXITER ) then
         call C_OUT(2,0,1,
     1        'resto_tron: ltron_iter exceeds limit. Abort.')
         IERR = 535
         goto 9999
      endif

C     Check filter stopping criteria
      if( task.eq.'NEWX' ) then
C
C     Compute new trial slack variables
C
         call CHECK_SLACKS(N, X, NLB, ILB, S_L, V_L, BNDS_L,
     1        S_L, NUB, IUB, S_U, V_U, BNDS_U, S_U, MU, MACHEPS)
C
C     Compute value of barrier function at new point
C
         call GET_F(N, X, IVAR, NORIG, XORIG, M, CSCALE, NLB,
     1        ILB, S_L, NUB, IUB, S_U, MU, F,
     2        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'resto_tron: get_f returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         bar_n = CALC_BAR(NLB, NUB, S_L, S_U, V_L, V_U, MU)
         phi_new = F + bar_n
         if( abs(QMERIT).eq.5 ) then
            phi_new = phi_new + DDOT(M, C, 1, LAM, 1)
         endif

         if( FFINITE(phi_new).eq.0 ) then ! we don't want NaN here; cut step!

            call C_OUT(2,0,1,' --- WARNING: Evaluation of barrier'//
     1           ' objective function produced Nan/Inf.')
            if( QCNR.gt.0 .and. QPRINT.ge.5 ) then
               write(line,*) 'F = ',F
               call C_OUT(1,0,1,line)
               do i = 1, NLB
                  k = ILB(i)
                  write(line,*) 'S_L(',k,') = ',S_L(i), X(k)
                  call C_OUT(1,0,1,line)
               enddo
               do i = 1, NUB
                  k = IUB(i)
                  write(line,*) 'S_U(',k,') = ',S_U(i), X(k)
                  call C_OUT(1,0,1,line)
               enddo
            endif
            write(line,*) 'resto_tron: in restoration: phi_new =',
     1           phi_new
            call C_OUT(2,0,1,line)
CTODO decide what to do in that case
C            goto 150
            IERR = 16
            goto 9999
         endif

C
C     Check if this point is acceptable to filter
C
         CNRM = CALC_NRM(M, C)
         do i = 1, NFILTER
C     Check f component
            fi = FILTER_PHI(i)
            lhs = phi_new - fi - PRECFACT*MACHEPS*dabs(fi)
            if( lhs.le.0d0 ) goto 3700
C     Check h component
            ci = FILTER_C(i)
            lhs = CNRM - ci - PRECFACT*MACHEPS*ci
            if( lhs.le.0d0 ) goto 3700
C     Point not OK for this entry
            if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
               write(line,*)
     1           'In restoration: new point not compatible with entry',
     2              i,' of filter.'
               call C_OUT(1,2,1,line)
            endif
            goto 150
C
 3700       continue
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
            goto 150
         endif

C     Filter tests passed

C
C     Correct iteration counter
C
         goto 200

      endif

 150  continue
C
C     Check if problem is locally infeasible
C
      if( task.eq.'NEWX' .or. task.eq.'START' ) then
#ifndef HAVE_TRON
         IERR = 97
         call C_OUT(2,0,1, 'resto_tron:  Compiled without TRON.')
         call C_OUT(2,0,1,
     1        '             Choose a different value for IMERIT.')
         goto 9999
#else
         gnorm = DGPNRM2(N, X, RW(p_xl+1), RW(p_xu+1), RW(p_g+1))
#endif
         if( TASK.eq.'START') then
            count_rej = 0
            TRflag = ' '
         else
            ITER = ITER + 1
         endif
         if(mod(ITER-iterstart,10).eq.0) then
            write(line,114)
            call C_OUT(2,0,2,line)
         endif
 114     format(/,
     1        'ITER    ||C||   ||p-grd||   Delta    #rej   #cg    F')

         COUNT_RESTO_ITER = COUNT_RESTO_ITER + 1
         write(line,115) ITER, cnrm2, gnorm, delta, TRflag, count_rej,
     1        COUNT_TRON_CG-cgcount_start, F
 115     format(i5,d10.3,d10.3,d10.3,a1,i4,i5,d11.3)
         cgcount_start = COUNT_TRON_CG

         count_rej = 0
         call C_OUT(2,0,1,line)
CTODO
         if( gnorm.le.QTRONGTOL ) then
            write(line,116) gnorm
 116        format(/,' resto_tron: GTOL test with gnorm = ',d14.7,
     1           ' satisfied.',/,
     1           '             Problem could be locally infeasible.'/)
            call C_OUT(2,0,4,line)
            IERR = 17
            goto 9999
         endif
      endif

      if( QTRONHESS.eq.0 ) then
         call GET_C(-ltron_iter, N, NIND, XHESS, IVAR, NORIG, XORIG, M,
     1        CSCALE, CHESS, KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     3        IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     2        EV_HCV, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*)
     1           'resto_tron: Error: get_c returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      endif

c     call the optimizer

      soflag = (QTRON2DERIVS.eq.1)
#ifndef HAVE_TRON
      IERR = 97
      call C_OUT(2,0,1, 'resto_tron:  Compiled without TRON.')
      call C_OUT(2,0,1,
     1     '             Choose a different value for IMERIT.')
      goto 9999
#else
      call dtron(N, X, RW(p_xl+1), RW(p_xu+1), (cnrm2**2)/2d0,
     1     RW(p_g+1), RW(p_h+1), RW(p_hdiag+1), IW(p_hcol_ptr+1),
     2     IW(p_hrow_ind+1), frtol, fatol, fmin, cgtol, itermax, DELTA,
     3     task, RW(p_b+1), RW(p_bdiag+1), IW(p_bcol_ptr+1),
     4     IW(p_brow_ind+1), RW(p_l+1), RW(p_ldiag+1), IW(p_lcol_ptr+1),
     5     IW(p_lrow_ind+1), RW(p_xc+1), RW(p_s+1), IW(p_indfree+1),
     6     IW(p_isave+1), RW(p_dsave+1), RW(p_wa+1), IW(p_iwa+1),
     7     TRflag,
     1     ITER, N, NIND, M, NORIG, XORIG, CSCALE, IVAR,
     1     NFIX, IFIX, soflag, CHESS,
     1     NLB, ILB, NUB, IUB, S_L, S_U,
     1     KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     1     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F,
     1     EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
#endif
      if( IERR.ne.0 ) then
         write(line,*) 'resto_tron: dtron returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( TRflag.eq.'r' ) count_rej = count_rej+1

      if( task(1:4).eq.'CONV' .or.
     1     (task.ne.'F' .and. task.ne.'GH' .and. task.ne.'NEWX') ) then
         call C_OUT(2,0,1,'resto_tron: output from TRON:')
         call C_OUT(2,0,1,task)
         if( index(task,'FRTOL').eq.14 .or.
     1        index(task,'FATOL').eq.14 ) then
            IERR = 18
         elseif( task(1:4).eq.'NAN ' ) then
            IERR = 591
            goto 9999
         endif
         goto 200
      endif

      goto 111

c           ---------- the end of the loop -------------

 200  continue
#ifndef HAVE_TRON
      IERR = 97
      call C_OUT(2,0,1, 'resto_tron:  Compiled without TRON.')
      call C_OUT(2,0,1,
     1     '             Choose a different value for IMERIT.')
      goto 9999
#else
      gnorm = DGPNRM2(N, X, RW(p_xl+1), RW(p_xu+1), RW(p_g+1))
#endif
      write(line,115) ITER+1, cnrm2, gnorm, delta, TRflag, count_rej,
     1     COUNT_TRON_CG-cgcount_start,  F
      call C_OUT(2,0,1,line)

      call C_OUT(2,0,1,'---- Leaving restoration phase ----')

      if( .not.relax_bounds .and. gnorm.le.QTRONGTOL ) IERR = 17

 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine RESTO_TRON_WS(N, X, M, NLB, NUB, NZA, LRW, LIW,
     1     EV_H, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*), X(N)
      integer IDAT(*)
      external EV_H
      integer lrw1, liw1, lrw2, liw2, lata, nnzata, nnzchc, nnzh
      integer idummy
      double precision dummy
      character*80 line

      call GET_C_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)
      call GET_F_WS(N, M, NLB, NUB, NZA, lrw1, liw1)
      LRW = max(LRW, lrw1)
      LIW = max(LIW, liw1)
      call CONSTR_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
      LRW = max(LRW, lrw1)
      LIW = max(LIW, liw1)
#ifdef HAVE_TRON
      call DTRON_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
      LRW = max(LRW, lrw1)
      LIW = max(LIW, liw1)
#endif

      LRW = 13*N + LRW
      LIW = 5*N + LIW
C
C     If QTRONHESS not 0, then add stuff for that.  However, size of
C     A^T A is not predictable, and the estimate below might be too small;
C     for that, dynamic memory allocation would be nice (later :)
C
      if( QTRONHESS.ne.0 ) then
         nnzata = 3*NZA+N
         lata   = max(nnzata, NZA+N)
         if( QTRON2DERIVS.eq.0 ) then
            nnzchc = 0
         else
            call GET_H(-1, N, N-M, idummy, dummy, idummy, N,
     1           X, dummy, M, nnzchc, dummy, dummy, idummy, idummy,
     1           idummy, dummy, idummy, idummy, idummy, EV_H, DAT, IDAT)
         endif
         nnzh = nnzata + nnzchc

         LRW = LRW + 3*N + 2*max(lata, nnzh) + 2*nnzh + 5*N
         LIW = LIW + 3*N + 1 + 2*max(lata, nnzh) + 2*nnzh + 5*N

         liw1 = NZA + 2*M + N + 1
         lrw1 = NZA

         if( QTRON2DERIVS.ne.0 ) then
            liw2 = 2*N+2+nnzh + 2*nnzchc
            lrw2 = 2*nnzchc
            liw1 = max(liw1, liw2)
            lrw1 = max(lrw1, lrw2)
         endif

         LRW = LRW + lrw1
         LIW = LIW + liw1
      endif

      if( QPRINT.ge.4 ) then
         write(line,1000)'resto_tron_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
