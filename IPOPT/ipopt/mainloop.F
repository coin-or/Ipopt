C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine MAINLOOP(ITER, N, NIND, M, NORIG, XORIG, IVAR, X, NFIX,
     1     IFIX, NLB, ILB, NUB, IUB, BNDS_L, BNDS_U, NZORIG, V_L, V_U,
     2     LAM, F, C, CSCALE, KCONSTR, LRS, LRS_END, RS, LIS, LIS_END,
     3     IS, LRW, RW, LIW, IW, IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)

C
C*******************************************************************************
C
C    $Id: mainloop.F 655 2004-10-05 17:23:15Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Main (outer) loop of algorithm
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
CP   ITER     I/O   INT    iteration counter (if -1, initialize)
CP                         I: is expected to be 0 at first call for algorithm
CP   N         I    INT    number of variables (without fixed)
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of constraints
CP   NORIG     I    INT    total number of variables (incl. fixed vars)
CP   XORIG    I/O   INT    I: user provided initial point (incl fixed vars)
CP                         O: solution (hopefully)
CP                            (original order as in problem statement)
CP   IVAR     I/O   INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP                            X(i) corresponds to XORIG(IVAR(i))
CP                         I: no partition done, only sort out fixed vars
CP   X        I/O   INT    as XORIG but without fixed vars
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
CP   BNDS_L   I/O   DP     values of lower bounds
CP   BNDS_U   I/O   DP     values of upper bounds
CP   NZORIG    I    INT    number of nonzeros in Jacobian of constraints
CP                            (including rows for FIXED variables!)
CP   V_L       O    DP     multipliers corresponding to the lower bounds
CP   V_U       O    DP     multipliers corresponding to the upper bounds
CP   LAM       O    DP     multipliers corresponding to the equality
CP                            constraints
CP   F         O    DP     objective function value
CP   C         O    DP     value of equality constraints
CP   CSCALE    I    DP     information on scaling (see get_scale)
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
CS    INITPOINT
CS    PARTITON
CS    REORDER_IB
CS    REORDER_X
CS    GET_F
CS    GET_G
CS    GET_C
CS    GET_YPY
CS    GET_RG
CS    GET_RV
CS    OPTERROR
CS    UPDATE_MU
CS    SOLVE_BARRIER
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
      integer IVAR(N)
      double precision X(N)
      integer NFIX
      integer IFIX(NFIX)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision BNDS_L(NLB)
      double precision BNDS_U(NUB)
      integer NZORIG
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision LAM(M)
      double precision F
      double precision C(M)
      double precision CSCALE(*)
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
      integer P_SL, P_SU, P_G, P_YPY, P_RG, P_RV
      save    P_SL, P_SU, P_G, P_YPY, P_RG, P_RV
      integer P_SIGMAL, P_SIGMAU, P_B, P_W, P_WCORR, P_RGOLD
      save    P_SIGMAL, P_SIGMAU, P_B, P_W, P_WCORR, P_RGOLD
      integer P_RGB, P_PZ, P_ZPZ, P_DX, P_DVL, P_DVU, P_LAMOLD
      save    P_RGB, P_PZ, P_ZPZ, P_DX, P_DVL, P_DVU, P_LAMOLD
      integer P_GOLD, P_RHS_RESTO, P_SIGMA_PN, RESTO
      save    P_GOLD, P_RHS_RESTO, P_SIGMA_PN, RESTO
      double precision MU, MUORIG, ERR_BAR, CONDC
      save             MU, MUORIG, ERR_BAR, CONDC
      logical FIRSTMU, NEWBAS
      save    FIRSTMU, NEWBAS

      logical INIT_LAM
      save    INIT_LAM

      integer p_rwend, p_iwend
      integer p_ivarold, p_ivarnew1, p_ib, p_xtmp
      integer i, icmax, blen, ierr_tmp
      double precision err_dinf_us, err_pinf_us, cmax, dummy
      double precision err_0, err_dinf, err_dinf_0, err_pinf, err_cmpl_0
      double precision err_cmpl_us, err_pinf_0, dummy2, scaled
      logical small_x_stop, ex

      character*1 c_err
      character*128 line(4)

      integer IDAMAX, FFINITE
      double precision DNRM2

#ifdef USE_IP_MYCMAX
      double precision IP_MYCMAX
#endif

      double precision mu_old, BETA, BETA_FACT
      parameter( BETA = 0.25d0, BETA_FACT = 10.d0 )
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Initialize pointer to stored variables
C
      if( ITER.eq.-1 ) then
         P_SL     = LRS_END
         P_SU     = P_SL     + NLB
         P_G      = P_SU     + NUB
         P_YPY    = P_G      + N
         P_RG     = P_YPY    + N
         P_RV     = P_RG     + NIND
         P_SIGMAL = P_RV     + NIND
         P_SIGMAU = P_SIGMAL + NLB
         if( QQUASI.ne.0 .or. QFULL.eq.0 ) then
            if( QQUASI.lt.0.and.QQUASI.gt.-6 ) then
               P_GOLD = P_SIGMAU + NUB
               P_B    = P_GOLD   + N
            else
               P_GOLD = 0
               P_B    = P_SIGMAU + NUB
            endif
C            if( QFULL.eq.0 .and. QCG.ne.0 ) then
C               P_W      = 0     ! don't need W for CG
C               P_WCORR  = P_B      + (NIND*(NIND+1))/2
C            else
            if( abs(QQUASI).ge.6 ) then
               P_W      = P_B + 2+2*N+QLMLEN*(1+4*N+2*QLMLEN)
               P_WCORR  = P_W
            else
               P_W      = P_B      + (NIND*(NIND+1))/2
               P_WCORR  = P_W      + (NIND*(NIND+1))/2
            endif
C            endif
         else
            P_GOLD   = 0
            P_B      = 0
            P_W      = 0
            P_WCORR  = P_SIGMAU + NUB
         endif
         P_RGOLD  = P_WCORR  + max(1, NIND) ! 1 for quadractic ls models
         P_RGB    = P_RGOLD  + NIND
         P_PZ     = P_RGB    + NIND
         P_ZPZ    = P_PZ     + NIND
         P_DX     = P_ZPZ    + M
         P_DVL    = P_DX     + N
         P_DVU    = P_DVL    + NLB
         P_LAMOLD = P_DVU    + NUB
         LRS_END  = P_LAMOLD + M

         if( abs(QMERIT).ge.4 .and. QRESTO.eq.1 ) then
            P_RHS_RESTO = LRS_END
            P_SIGMA_PN  = P_RHS_RESTO + 3*N+2*M
            LRS_END     = P_SIGMA_PN + N
         else
            P_RHS_RESTO = 0
            P_SIGMA_PN  = 0
         endif
C
C     Initialize INIT_LAM (for get_step_full)
C
         if( abs(QLAMBDA).eq.2 .and. QINIT.ne.0 ) then
            if( QLAMINIT.eq.1 .and. M.gt.0 ) then
               INIT_LAM = .true.
            else
               INIT_LAM = .false.
            endif
         else
            INIT_LAM = .false.
         endif

         goto 9999
      endif

      p_rwend = 0
      p_iwend = 0

      if( MEMDBG ) then
         write(line,1)'mainloop', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif
C
C     Make sure, no STOP file exists
C
      inquire(file='STOP', exist=ex)
      if( ex ) then
         open(80,file='STOP',status='old')
         close(80,status='delete')
      endif
C
C     Make sure, no DECMU file exists
C
      inquire(file='DECMU', exist=ex)
      if( ex ) then
         open(80,file='DECMU',status='old')
         close(80,status='delete')
      endif
C
C     initialize flag for barrier parameter update
C
      FIRSTMU = .true.
C
C     For filter: We start not in restoration phase
C
      RESTO = 0
      small_x_stop = .false.
C
C     Compute initial point (primal, slack and dual variables) and
C        determine initial barrier parameter
C
      call INITPOINT(N, X, NORIG, XORIG, IVAR, M, NIND,
     1     NLB, ILB, BNDS_L, V_L, RS(P_SL+1), NUB, IUB, BNDS_U, V_U,
     2     RS(P_SU+1), MU, CSCALE, NFIX, IFIX,
     2     KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'mainloop: INITPOINT returns IERR = ', IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     Determine initial partition IVAR
C     (only if necessary for non-full space approach and qN udpates)
C
C      if( QQUASI.ne.0 .or. QFULL.eq.0 .or. QLAMBDA.eq.1 .or.
C     1     abs(QERROR).eq.1 .or. QLAMBDA.eq.0 ) then
      if( QFULL.eq.0 .or. QLAMBDA.eq.1 .or.
     1     (QQUASI.ne.0.and.abs(QQUASI).lt.6) .or.
     1     abs(QERROR).eq.1 .or. QLAMBDA.eq.0 ) then
C
C     Store old partition IVAR in IVAROLD
C
         p_ivarold = p_iwend
         p_iwend   = p_ivarold + N
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         do i = 1, N
            IW(p_ivarold+i) = IVAR(i)
         enddo

         call PARTITION(QSELBAS, ITER, N, NIND, M, NFIX, IFIX, IVAR,
     1        NORIG, XORIG, CSCALE, KCONSTR, LRS, RS, LIS, IS,
     3        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'mainloop: Error: partition ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'mainloop: Warning: partition ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
C
C     Reorder bound pointers ILB, IUB
C
         p_ivarnew1 = p_iwend
         p_ib       = p_ivarnew1 + NORIG
         p_iwend    = p_ib       + N
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         call REORDER_IB(N, NORIG, NLB, ILB, NUB, IUB,
     1        IW(p_ivarold+1), IVAR, IW(p_ivarnew1+1),
     2        IW(p_ib+1))
         p_iwend = p_ivarnew1
C
C     Reorder X
C
         p_xtmp  = p_rwend
         p_rwend = p_xtmp  + NORIG
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         call REORDER_X(N, NORIG, X, IW(p_ivarold+1), IVAR,
     1        RW(p_xtmp+1))
         if( QSCALE.ge.3 ) then
             call REORDER_X(N, NORIG, CSCALE(M+1), IW(p_ivarold+1),
     1           IVAR, RW(p_xtmp+1))
          endif
         p_rwend = p_xtmp
         p_iwend = p_ivarold
      endif
C
C     Compute value of objective function F
C
      call GET_F(N, X, IVAR, NORIG, XORIG, M, CSCALE, NLB, ILB,
     1     RS(P_SL+1), NUB, IUB, RS(P_SU+1), MU, F, LIW-p_iwend,
     2     IW(p_iwend+1), IERR, EV_F, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'mainloop: get_f returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( FFINITE(F).eq.0 ) then
         call C_OUT(2,0,1,' --- ERROR: Evaluation of'//
     1        ' objective function at x0 produced Nan/Inf.')
         IERR = 11
         goto 9999
      endif
C
C     Compute gradient of objective function G
C
      call GET_G(N, X, IVAR, NORIG, XORIG, M, CSCALE, NLB, ILB,
     1     RS(P_SL+1), NUB, IUB, RS(P_SU+1), MU, RS(P_G+1),
     1     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     1     IERR, EV_G, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1        'mainloop: Error: get_g ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.ne.0 ) then
         write(line,*)
     1        'mainloop: Warning: get_g ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
      endif
      dummy = DNRM2(N, RS(P_G+1), 1)
      if( FFINITE(dummy).eq.0 ) then
         call C_OUT(2,0,1,' --- ERROR: Evaluation of'//
     1        ' objective gradient at x0 produced Nan/Inf.')
         IERR = 11
         goto 9999
      endif
C
C     Compute values of constraints C
C
      call GET_C(ITER, N, NIND, X, IVAR, NORIG, XORIG, M,
     1     CSCALE, C, KCONSTR, LRS, RS, LIS, IS,
     2     LRW-p_rwend, RW(p_rwend+1),
     4     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1        'mainloop: Error: get_c returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.ne.0 ) then
         write(line,*)
     1        'mainloop: Warning: get_c returns IERR = ',IERR
         call C_OUT(2,0,1,line)
      endif
      dummy = DNRM2(M, C, 1)
      if( FFINITE(dummy).eq.0 ) then
         call C_OUT(2,0,1,' --- ERROR: Evaluation of'//
     1        ' constraints at x0 produced Nan/Inf.')
         IERR = 11
         goto 9999
      endif
C
C     Obtain first range space step PY
C
      if( QFULL.eq.0 .or. (QQUASI.ne.0.and.abs(QQUASI).lt.6)
     1     .or. QLAMBDA.eq.0 ) then
         call GET_YPY(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, NLB, ILB, NUB, IUB,
     2        RS(P_SL+1), RS(P_SU+1), C, RS(P_YPY+1),
     4        NEWBAS, CONDC, KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'mainloop: Error: get_ypy returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'mainloop: Warning: get_ypy returns IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
      else
         NEWBAS = .false.
CTODO This is bad, if we don't allocate memory for YPY!
         call DCOPY(N, 0d0, 0, RS(P_YPY+1), 1)
      endif
C
C     Compute reduced gradient of objective function RG
C
      if( QFULL.eq.0 .or. (QQUASI.ne.0.and.abs(QQUASI).lt.6)
     1     .or. abs(QERROR).eq.1 ) then
         call GET_RG(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, RS(P_G+1), RS(P_RG+1),
     1        KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'mainloop: Error: get_rg ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'mainloop: Warning: get_rg ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
      endif
C
C     Compute "reduced" dual variables RV
C
      if( QFULL.eq.0 .or. abs(QERROR).eq.1 ) then
         call GET_RV(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, NLB, ILB, NUB, IUB,
     1        V_L, V_U, RS(P_RV+1), KCONSTR, LRS, RS, LIS, IS,
     2        LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'mainloop: Error: get_rv ends with IERR = ',IERR
            goto 9999
            call C_OUT(2,0,1,line)
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'mainloop: Warning: get_rv ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
      endif
C
C     For QLAMBDA = 1: get coordinate multipliers
C
      if( QLAMBDA.eq.1 ) then
         call GET_LAMBDA(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, RS(P_G+1), NLB, ILB, NUB,
     1        IUB, V_L, V_U, LAM, KCONSTR, LRS, RS, LIS, IS,
     1        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     1        IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'mainloop: get_lambda returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.lt.0 ) then
            write(line,*)
     1           'mainloop: get_lambda returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 0
         endif
      elseif( QLAMBDA.ne.0 .and. QINIT.ne.0 ) then
         call DCOPY(M, 0.d0, 0, LAM, 1)
      endif
C
C     Initialize old copy of Lagrangian multipliers
C
      if( QLAMBDA.ne.0 ) then
         call DCOPY(M, 0.d0, 0, RS(P_LAMOLD+1), 1)
      endif
C

C ------------------------------------------------------------------------------
C     Start of main loop !
C ------------------------------------------------------------------------------
 100  continue
      if( err_bar.le.QTOL .and. small_x_stop ) then
         IERR = 9
         goto 9000
      endif
      mu_old = MU
      call UPDATE_MU(FIRSTMU, RESTO, MU, err_bar)
      if( NLB.eq.0 .and. NUB.eq.0 ) then
C     Avoid MU and err_bar becoming too small if err_dinf_us not yet
C     small enough
         err_bar = max(1.d-1*QTOL, err_bar)
         MU = max(1.d-1*QTOL, MU)
      endif
C
C     Adapt F and G's linear damping terms
C
      call UPDATE_FG_MU(0, N, NLB, ILB, RS(P_SL+1), NUB, IUB,
     1     RS(P_SU+1), mu_old, MU, F, RS(P_G+1), LIW-p_iwend,
     2     IW(p_iwend+1), IERR)
      if( IERR.ne.0 ) then
         write(line,*) 'mainloop: update_fg_mu returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9000
      endif
C
C     Solve the barrier problem (uff...)
C
      small_x_stop = .false.
      call SOLVE_BARRIER(ITER, N, NIND, M, NORIG, XORIG, CSCALE,
     1     IVAR, NFIX, IFIX, NLB, ILB, NUB, IUB, BNDS_L, BNDS_U, MU,
     1     MUORIG, X, RS(P_SL+1), RS(P_SU+1), V_L, V_U, F, RS(P_G+1),
     1     C, ERR_BAR, ERR_0, ERR_DINF, ERR_DINF_0, ERR_DINF_US,
     1     ERR_PINF, ERR_PINF_0, ERR_PINF_US, ERR_CMPL_0, ERR_CMPL_US,
     1     RS(P_RV+1), RS(P_SIGMAL+1),
     1     RS(P_SIGMAU+1), RS(P_B+1), RS(P_W+1),
     1     RS(P_WCORR+1), RS(P_RG+1), RS(P_RGOLD+1),
     1     RS(P_RGB+1), RS(P_YPY+1), RS(P_RHS_RESTO+1),
     1     RS(P_SIGMA_PN+1), NEWBAS, CONDC, RS(P_PZ+1), RS(P_ZPZ+1),
     1     LAM, RS(P_LAMOLD+1), RS(P_DX+1),
     1     RS(P_DVL+1), RS(P_DVU+1), RS(P_GOLD+1), RESTO, INIT_LAM,
     1     KCONSTR, LRS, LRS_END, RS, LIS, LIS_END, IS,
     1     LRW-p_rwend, RW(p_rwend+1),
     1     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1      'mainloop: Error: solve_barrier ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
         if( IERR.eq.535 ) then
            IERR = 1
         endif
         goto 9000
      elseif( IERR.eq.-1 ) then
         IERR = 0               ! Done!
         goto 9000
      elseif( IERR.eq.-2 ) then
         small_x_stop = .true.
         IERR = 0
      elseif( IERR.eq.-3 ) then
         IERR = 0
         goto 9000
      elseif( IERR.ne.0 ) then
         write(line,*)
     1      'mainloop: Warning: solve_barrier ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
         IERR = 0
      endif
C ------------------------------------------------------------------------------
C     End of main loop
C ------------------------------------------------------------------------------
      goto 100

 9000 continue
#ifdef USE_MALLOC
C
C     Free allocated memory for factorization
C
      call GSF_FREEPTR(i)
      if( i.ne.0 ) then
         call C_OUT(2,0,1,'mainloop: Warning during freeing memory.')
         write(line,*) '       GSF_FREEPTR returns IERR = ',i
         call C_OUT(2,0,1,line)
      endif
#endif
C
C     Might have to compute error for final report
C
      if( IERR.ne.0 ) then
         call OPTERROR(N, NIND, M, X, RS(P_G+1), RS(P_RG+1), RS(P_RV+1),
     1        C, CSCALE, NLB, ILB, NUB, IUB, BNDS_L,
     1        BNDS_U, RS(P_SL+1), RS(P_SU+1), V_L, V_U,
     1        LAM, 0.d0, dummy, c_err, ERR_0, ERR_DINF,
     1        ERR_DINF_0, ERR_DINF_US, ERR_PINF, ERR_PINF_0,
     1        ERR_PINF_US, dummy2, ERR_CMPL_0, ERR_CMPL_US,
     1        ITER+1, IVAR, NFIX, IFIX, NORIG, XORIG, RS(P_DX+1),
     1        KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     1        LIW-p_iwend, IW(p_iwend+1), ierr_tmp, EV_F, EV_C, EV_G,
     5        EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( ierr_tmp.gt.0 ) then
            write(line,*)
     1           'mainloop: Error: Opterror ends with IERR = ',ierr_tmp
            call C_OUT(2,0,1,line)
            IERR = ierr_tmp
            goto 9999
         elseif( ierr_tmp.ne.0 ) then
            write(line,*)
     1           'mainloop: Warning:Error ends with IERR = ',ierr_tmp
            call C_OUT(2,0,1,line)
         endif
      endif
C
C     If desired, write current quasi-Newton estimate to file BWARM.DAT
C
      if( QBWARMSTART.eq.1 .and. QQUASI.ne.0 ) then
         if( abs(QQUASI).eq.6 ) then
            blen = 2+2*N+QLMLEN*(1+4*N+2*QLMLEN)
         else
            blen = (NIND*(NIND+1))/2
         endif
         call C_OUT(2,1,1,
     1  'mainloop: Writing current quasi-Newton estimate to BWARM.DAT.')
         open(80,file='BWARM.DAT',status='unknown',err=9008)
         write(80,9001,err=9009) blen
 9001    format(i16)
         do i = 1, blen
            write(80,9002,err=9009) RS(P_B+i)
 9002       format(d23.16)
         enddo
         close(80)
      endif
C
C     Undo the scaling of objetcive function
C
      call UPDATE_FG_MU(1, N, NLB, ILB, RS(P_SL+1), NUB, IUB,
     1     RS(P_SU+1), MUORIG, 0.d0, F, dummy, LIW, IW, ierr_tmp)
      if( ierr_tmp.ne.0 ) then
         write(line,*)
     1        'mainloop: UPDATE_FG_MU ends with IERR = ',ierr_tmp
         call C_OUT(2,0,1,line)
         IERR = ierr_tmp
         goto 9999
      endif
      F = F/QFSCALE

      write(line,9010) ITER
 9010 format(/,'Number of iterations taken .............',i23)
      call C_OUT(2,0,2,line)

      write(line,9020) F
 9020 format('Final value of objective function is....',d23.16)
      call C_OUT(2,0,1,line)

      write(line,9015)
 9015 format(/,'Errors at final point                      (scaled)
     1  (unscaled)')
      call C_OUT(2,0,2,line)

      if( M.gt.0 ) then
C
C     A sophisticated user might want to compute the constraint violation
C     differently (e.g. taking into account which variables are slack variables)
C     In that case you can compute 'your' constraint violation with the function
C     below.  At compile time, '-DUSE_IP_MYCMAX'.
C
#ifdef USE_IP_MYCMAX
         cmax = IP_MYCMAX(N, X, M, C, EV_C, DAT, IDAT)
#else
         cmax = err_pinf_us
#endif
         write(line,9021) err_pinf_0, cmax
      else
         write(line,9021) 0.d0, 0.d0
      endif
 9021 format('Final maximal constraint violation is...',d13.6,d16.6)
      call C_OUT(2,0,1,line)

CTODO      write(line,9022) err_dinf_us
      write(line,9022) err_dinf, err_dinf_us
 9022 format('Final value for dual infeasibility is...',d13.6,d16.6)
      call C_OUT(2,0,1,line)
      write(line,9023) err_cmpl_0, err_cmpl_us
 9023 format('Final value of complementarity error is.',d13.6,d16.6)
      call C_OUT(2,0,1,line)     

      write(line,9026) fevals
 9026 format(/,'The objective function was evaluated',i7,' times.')
      call C_OUT(2,0,2,line)

      write(line,9027) cevals
 9027 format('The constraints were evaluated',i13,' times.')
      call C_OUT(2,0,1,line)

      if( QPRINT.ge.2 .and. QCNR.gt.0 ) then
        write(line,9030)
 9030   format(/,'  Variables at Solution(?):',/)
        call C_OUT(1,0,3,line)
        do i = 1, NORIG
           write(line,9040) i,XORIG(i)
 9040      format(' X  (',i7,') = ',d23.15)
           call C_OUT(1,0,1,line)
        enddo
        write(line,9050)
 9050   format(/,'  Multipliers for lower Bounds:',/,/,
     1       24x,'(scaled)',16x,'(unscaled)')
        call C_OUT(1,0,4,line)
        do i = 1, NLB
           scaled = V_L(i)
           V_L(i) = V_L(i)/QFSCALE
           if( QSCALE.ge.3 ) then
              V_L(i) = V_L(i)/CSCALE(M+ILB(i))
           endif
           write(line,9060) IVAR(ILB(i)), scaled, V_L(i)
 9060      format(' V_L(',i7,') = ',d23.15,d25.15)
           call C_OUT(1,0,1,line)
        enddo
        write(line,9070)
 9070   format(/,'  Multipliers for upper Bounds:',/,/,
     1       24x,'(scaled)',16x,'(unscaled)')
        call C_OUT(1,0,3,line)
        do i = 1, NUB
           scaled = V_U(i)
           V_U(i) = V_U(i)/QFSCALE
           if( QSCALE.ge.3 ) then
              V_U(i) = V_U(i)/CSCALE(M+IUB(i))
           endif
           write(line,9080) IVAR(IUB(i)), scaled, V_U(i)
 9080      format(' V_U(',i7,') = ',d23.15,d25.15)
           call C_OUT(1,0,1,line)
        enddo
        write(line,9090)
 9090   format(/,'  Multipliers for equality constraints:',/,/,
     1       24x,'(scaled)',16x,'(unscaled)')
        call C_OUT(1,0,3,line)
        do i = 1, M
           scaled = LAM(i)
           if( QSCALE.lt.2 ) then
              LAM(i) = LAM(i)*CSCALE(1)
           else
              LAM(i) = LAM(i)*CSCALE(i)
           endif
           LAM(i) = LAM(i)/QFSCALE
           write(line,9100) i, scaled, LAM(i)
 9100      format(' LAM(',i7,') = ',d23.15,d25.15)
           call C_OUT(1,0,1,line)
        enddo
      else
        do i = 1, NLB
           V_L(i) = V_L(i)/QFSCALE
           if( QSCALE.ge.3 ) then
              V_L(i) = V_L(i)/CSCALE(M+ILB(i))
           endif
        enddo
        do i = 1, NUB
           V_U(i) = V_U(i)/QFSCALE
           if( QSCALE.ge.3 ) then
              V_U(i) = V_U(i)/CSCALE(M+IUB(i))
           endif
        enddo
        do i = 1, M
           if( QSCALE.lt.2 ) then
              LAM(i) = LAM(i)*CSCALE(1)
           else
              LAM(i) = LAM(i)*CSCALE(i)
           endif
           LAM(i) = LAM(i)/QFSCALE
        enddo
      endif
C
C     A sophisticated user might want to use his/her own routine for printing
C     out information of the solution.  This may be done with the routine below.
C     At compile time, '-DUSE_IP_MYWRITECONSTR'.
C
#ifdef USE_IP_MYWRITESOL
      call IP_MYWRITESOL(N, X, M, C, LAM, NLB, ILB,
     1     V_L, NUB, IUB, V_U)
#else
      if( QPRINT.ge.2 .and. QCNR.gt.0 ) then
        write(line,9110)
 9110   format(/,'  Constraint values:',/,/,
     1       24x,'(scaled)',16x,'(unscaled)')
        call C_OUT(1,0,3,line)
        do i = 1, M
           scaled = C(i)
           if( QSCALE.lt.2 ) then
              C(i) = C(i)/CSCALE(1)
           else
              C(i) = C(i)/CSCALE(i)
           endif
           write(line,9120) i,scaled,C(i)
 9120      format(' C  (',i7,') = ',d23.15,d25.15)
           call C_OUT(1,0,1,line)
        enddo
      endif
#endif

#ifdef OCC
C -------------------------------------------------------------------
C     Report Output to <opt_details.m>
C     (Those lines are contributed by Tobias Jockenhoevel in order
C      to be able to use IPOPT with his OptimizationControlCentre)
C -------------------------------------------------------------------
      open (UNIT = 10, FILE = 'opt_details.m')
      write(10,*) ' % -----------------------------------'
      write(10,*) ' % IPOPT Optimization runtime report        '
      write(10,*) ' % -----------------------------------'
      write(10,*) ' '
      write(10,*) ' % Final value of objective function:'
      write(10,*) ' opt_obj = ',F,';'
      write(10,*) ' '
      write(10,*) ' % Final maximal constraint violation:'
      if( M.gt.0 ) then
         write(10,*) ' c_err = ',dabs(C(icmax)),';'
      else
         write(10,*) ' c_err = 0e0;'
      endif
      write(10,*) ' '
      write(10,*) ' % Final value for overall KKT error:'
      write(10,*) ' kkt_err = ', ERR,';'
      write(10,*) ' '
      write(10,*) ' % Evaluation of the objective function:'
      write(10,*) ' obj_iter = ',fevals,';'
      write(10,*) ' '
      write(10,*) ' % Evaluation of the constraints:'
      write(10,*) ' c_iter = ',cevals,';'

      close(10, STATUS = 'KEEP')
#endif
C
C     END
C
 9999 continue
      return

 9008 IERR = 8
      call C_OUT(2,0,1,'Error while trying to open BWARM.DAT')
      goto 9999
 9009 IERR = 8
      call C_OUT(2,0,1,'Error while writing to BWARM.DAT')
      goto 9999
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine MAINLOOP_WS(N, X, M, NLB, NUB, NZA, LRW, LIW,
     1     FILLINFACT, EV_H, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      external EV_H
      double precision DAT(*)
      integer IDAT(*)
      double precision FILLINFACT, X(N)
      integer lrw1, liw1
      character*80 line

      call INITPOINT_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      if( QFULL.eq.0 .or. QLAMBDA.eq.1 .or.
     1     (QQUASI.ne.0.and.abs(QQUASI).lt.6) .or.
     1     abs(QERROR).eq.1 .or. QLAMBDA.eq.0 ) then
         call PARTITION_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW, max(lrw1,N))
         LIW = max(LIW, max(liw1,3*N))
      endif

      call GET_G_WS(N, M, NLB, NUB, NZA, lrw1, liw1)
      LRW = max(LRW, lrw1)
      LIW = max(LIW, liw1)

      if( QFULL.eq.0 .or. (QQUASI.ne.0.and.abs(QQUASI).lt.6)
     1     .or. QLAMBDA.eq.0 ) then
         call GET_YPY_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)
      endif

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

      call OPTERROR_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
      LRW = max(LRW, lrw1)
      LIW = max(LIW, liw1)

      call SOLVE_BARRIER_WS(N, X, M, NLB, NUB, NZA, lrw1, liw1,
     1     FILLINFACT, EV_H, DAT, IDAT)
      LRW = max(LRW, lrw1)
      LIW = max(LIW, liw1)

      if( QPRINT.ge.4 ) then
         write(line,1000)'mainloop_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end


