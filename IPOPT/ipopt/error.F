C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine OPTERROR(N, NIND, M, X, G, RG, RV, C, CSCALE, NLB, ILB,
     1     NUB, IUB, BNDS_L, BNDS_U, S_L, S_U, V_L, V_U,
     1     LAM, MU, ERR, C_ERR, ERR_0, ERR_DINF, ERR_DINF_0,
     1     ERR_DINF_US, ERR_PINF, ERR_PINF_0, ERR_PINF_US, ERR_CMPL,
     1     ERR_CMPL_0, ERR_CMPL_US, ITER, IVAR, NFIX, IFIX, NORIG,
     1     XORIG, DX, KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR,
     1     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     1     DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: error.F 655 2004-10-05 17:23:15Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    reduced KKT error for barrier problem
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
CP   N         I    INT    number of (free) variables; first M variables
CP                         are dependent; remaining independent
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of equality constraints
CP   X         I    DP     actual primal iterate
CP   G         I    DP     gradient of objective function
CP   RG        I    DP     reduced gradient of objective function
CP   RV        I    DP     reduced dual variables
CP   C         I    DP     values of constraints
CP   CSCALE    I    DP     scaling factors for constraints
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. S_L(i) is slack for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   ILB       I    INT    indices of upper bounds
CP                            (e.g. S_U(i) is slack for X(IUB(i)) )
CP   BNDS_L    I    DP     values of lower bounds (ordered as S_L)
CP   BNDS_U    I    DP     values of upper bounds (ordered as S_U)
CP   S_L       I    DP     slack variables for lower bounds
CP   S_U       I    DP     slack variables for upper bounds
CP   V_L       I    DP     dual variables for lower bounds
CP   V_U       I    DP     dual variables for upper bounds
CP   LAM       I    DP     Lagrangian multipliers
CP   MU        I    DP     barrier parameter (=0, if Error should be computed
CP                            for overall NLP)
CP   ERR       O    DP     overall scaled error for barrier problem
CP   C_ERR     O    C*1    characterizes worst violation
CP                            ='d': dual infeasibility
CP                            ='p': primal infeasibility
CP                            ='c': complementarity
CP   ERR_0     O    DP     overall scaled error for original problem
CP   ERR_DINF  O    DP     scaled dual infeasibility
CP   ERR_DINF_0 O   DP     scaled dual infeasibility (with MU = 0 for linear
CP                            damping term)
CP   ERR_DINF_US O  DP     unscaled dual infeasibility
CP   ERR_PINF  O    DP     scaled primal infeasibility
CP   ERR_PINF_0 O   DP     scaled primal infeasibility (with MU = 0 for MPCC)
CP   ERR_PINF_US O  DP     unscaled primal infeasibility
CP   ERR_CMPL  O    DP     scaled violation of relaxed complementarity
CP   ERR_CMPL_0 O   DP     scaled violation of original complementarity
CP   ERR_CMPL_US O  DP     unscaled violation of original complementarity
CP   ITER      I    INT    iteration counter
CP                         if -1: Initialize LASTITER
CP                         otherwise: It is assumed that for subsequent calls
CP                                    with same ITER, only nrm2 has to be
CP                                    recomputed.
CP   IVAR      I    INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP   NORIG     I    INT    total number of all variables (incl. fixed vars)
CP   XORIG     I    DP     (only TASK = 1,2,3): actual iterate
CP                            XORIG is ordered in ORIGINAL order (i.e. not
CP                            partitioned into independent and dependent
CP                            variables)
CP   DX        I    DP     primal search direction from previous iteration
CP                            (needed for QLAMBDA=2)
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
CS    DASUM
CS    DCOPY
CS    DAXPY
CS    IDAMAX
CS    C_OUT
CS    CONSTR
CS    EVAL_CCERR
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
      integer N
      integer NIND
      integer M
      double precision X(N)
      double precision G(N)
      double precision RG(NIND)
      double precision RV(NIND)
      double precision C(M)
      double precision CSCALE(*)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision BNDS_L(NLB)
      double precision BNDS_U(NUB)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision LAM(M)
      double precision MU
      double precision ERR
      character*1 C_ERR
      double precision ERR_0
      double precision ERR_DINF
      double precision ERR_DINF_0
      double precision ERR_DINF_US
      double precision ERR_PINF
      double precision ERR_PINF_0
      double precision ERR_PINF_US
      double precision ERR_CMPL
      double precision ERR_CMPL_0
      double precision ERR_CMPL_US
      integer ITER
      integer IVAR(N)
      integer NFIX
      integer IFIX(NFIX)
      integer NORIG
      double precision XORIG(NORIG)
      double precision DX(N)
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
      double precision DASUM, DNRM2
      double precision dummy
      integer IDAMAX
      integer i, j, p_iwend, p_rwend, p_tmp, idummy

      double precision SCAL_DINF, SCAL_CMPL
      save SCAL_DINF, SCAL_CMPL
      character*120 line(3)

      integer scaltype          ! = 0: no scaling
                                ! = 1: use mulitpliers

      integer LASTITER
      save LASTITER
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( ITER.eq.-1 ) then
         LASTITER = -1
         goto 9999
      endif

      IERR = 0
      p_rwend = 0
      p_iwend = 0

      if( QERROR.gt.0 ) then
         scaltype = 1
      elseif( QERROR.lt.0 ) then
         scaltype = 0
      else
         call C_OUT(2,0,1,'Error: Invalid value of QERROR!')
         IERR = 4
         goto 9999
      endif
      if( abs(QERROR).eq.2 .and. QLAMBDA.eq.0 ) then
         call C_OUT(2,0,1,'Error: Need multipliers!')
         IERR = 4
         goto 9999
      endif
C
C     First compute scaling factors
C
      if( LASTITER.ne.ITER ) then
         if( scaltype.eq.1 ) then
C     SCAL_DINF = 1+||lam,v||_1/(m+nbnds)
C     SCAL_CMPL = 1 + ||v||_1/(nbnds)
            if( QLAMBDA.eq.0 ) then ! no multipliers given
                                ! Do not know what to do; for now pretend
                                ! lam = e
               SCAL_DINF = dble(M)
            else
               SCAL_DINF = DASUM(M, LAM, 1)
            endif
            SCAL_CMPL = DASUM(NLB, V_L, 1) + DASUM(NUB, V_U, 1)
            SCAL_DINF = SCAL_DINF + SCAL_CMPL
            if( M+NLB+NUB.gt.0 ) then
               SCAL_DINF =
     1            max(QSCALMULTMAX,SCAL_DINF/dble(M+NLB+NUB))/
     2              QSCALMULTMAX
            else
               SCAL_DINF = 1.d0
            endif
            if( NLB+NUB.gt.0 ) then
               SCAL_CMPL =
     1              max(QSCALMULTMAX,SCAL_CMPL/dble(NLB+NUB))/
     2              QSCALMULTMAX
            else
               SCAL_CMPL = 1.d0
            endif
         else
            SCAL_DINF = 1.d0
            SCAL_CMPL = 1.d0
         endif
      endif
C
C     Compute dual infeasibility
C
      if( LASTITER.ne.ITER ) then
         if( NIND.eq.0 ) then
            ERR_DINF = 0.d0
            ERR_DINF_0 = 0.d0
            ERR_DINF_US = 0.d0
         else
            if( abs(QERROR).eq.1 ) then
C     DINF = Red_grad_barrier (ignore ISCALE scaling)
               p_tmp   = p_rwend
               p_rwend = p_tmp + NIND
               if( p_rwend.gt.LRW ) then
                  IERR = 98
                  goto 9999
               endif
               call DCOPY(NIND, RG, 1, RW(p_tmp+1), 1)
               call DAXPY(NIND, 1.d0, RV, 1, RW(p_tmp+1), 1)
               i = IDAMAX(NIND, RW(p_tmp+1), 1)
               ERR_DINF_US = dabs(RW(p_tmp+i))
               ERR_DINF = ERR_DINF/SCAL_DINF
               p_rwend = p_tmp
            else
               p_tmp   = p_rwend
               p_rwend = p_tmp + N
               if( p_rwend.gt.LRW ) then
                  IERR = 98
                  goto 9999
               endif
C
C     Compute A * LAM
C
               if( M.gt.0 ) then
                  call CONSTR(8, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1                 NORIG, XORIG, CSCALE, LAM, RW(p_tmp+1),
     2                 idummy, idummy,
     3                 KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4                 IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5                 LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C,
     5                 EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     7                 DAT, IDAT)
                  if( IERR.lt.0 ) then
                     write(line,*)
     1                    'Error: Warning in CONSTR-8, IERR = ', IERR
                     call C_OUT(2,0,1,line)
                  elseif( IERR.ne.0 ) then
                     write(line,*)
     1                    'Error: Error in CONSTR-8, IERR = ',IERR
                     call C_OUT(2,0,1,line)
                     goto 9999
                  endif
                  call DAXPY(N, 1d0, G, 1, RW(p_tmp+1), 1)
               else
                  call DCOPY(N, G, 1, RW(p_tmp+1), 1)
               endif
C
C     Add dual variables
C
               do i = 1, NLB
                  j = ILB(i)
                  RW(p_tmp+j) = RW(p_tmp+j) - V_L(i)
               enddo
               do i = 1, NUB
                  j = IUB(i)
                  RW(p_tmp+j) = RW(p_tmp+j) + V_U(i)
               enddo
C
C     Compute scaled dual infeasibility
C
               i = IDAMAX(N, RW(p_tmp+1), 1)
               ERR_DINF = dabs(RW(p_tmp+i))/SCAL_DINF
C
C     Compute scaled and unscaled dual infeasibility for original problem
C
               call UPDATE_FG_MU(2, N, NLB, ILB, S_L, NUB, IUB, S_U,
     1              MU, 0.d0, dummy, RW(p_tmp+1), LIW-p_iwend,
     2              IW(p_iwend+1), IERR)
               if( IERR.ne.0 ) then
                  write(line,*)
     1                 'error: update_fg_mu returns IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  goto 9999
               endif
               i = IDAMAX(N, RW(p_tmp+1), 1)
               ERR_DINF_0 = dabs(RW(p_tmp+i))/SCAL_DINF
               call DSCAL(N, 1.d0/abs(QFSCALE), RW(p_tmp+1), 1)
               if( QSCALE.ge.3 ) then
                  do i = 1, N
                     RW(p_tmp+i) = RW(p_tmp+i)/CSCALE(M+i)
                  enddo
               endif
               i = IDAMAX(N, RW(p_tmp+1), 1)
               ERR_DINF_US = dabs(RW(p_tmp+i))
               p_rwend = p_tmp
            endif
         endif
      endif
C
C     Compute primal infeasibility
C
      if( M.eq.0 ) then
         ERR_PINF       = 0.d0
         ERR_PINF_0     = 0.d0
         ERR_PINF_US    = 0.d0
      else
#ifndef INCLUDE_CC
         if( LASTITER.ne.ITER ) then
            i = IDAMAX(M, C, 1)
            ERR_PINF = abs(C(i))
            if( QSCALE.lt.2 ) then
               ERR_PINF_US = ERR_PINF/CSCALE(1)
            else
               ERR_PINF_US = 0.d0
               do i = 1, M
                  ERR_PINF_US = max(ERR_PINF_US, abs(C(i))/CSCALE(i))
               enddo
            endif
            ERR_PINF_0 = ERR_PINF
         endif
#else
         call EVAL_CCERR(N, M, X, C, 0.d0)
         i = IDAMAX(M, C, 1)
         ERR_PINF_0 = abs(C(i))
         if( QSCALE.lt.2 ) then
            i = IDAMAX(M, C, 1)
            ERR_PINF_US = abs(C(i))/CSCALE(1)
         else
            ERR_PINF_US = 0.d0
            do i = 1, M
               ERR_PINF_US = max(ERR_PINF_US, abs(C(i))/CSCALE(i))
            enddo
         endif

         call EVAL_CCERR(N, M, X, C, MU)
         i = IDAMAX(M, C, 1)
         ERR_PINF = abs(C(i))
#endif
      endif
C
C     Compute complementarity
C
      if( NIND.eq.0 ) then
         ERR_CMPL    = 0.d0
         ERR_CMPL_0  = 0.d0
         ERR_CMPL_US = 0.d0
      else
         ERR_CMPL = 0.d0
         ERR_CMPL_0 = 0.d0
         do i = 1, NLB
            ERR_CMPL = dmax1(ERR_CMPL,dabs(V_L(i)*S_L(i)-MU))
            ERR_CMPL_0 = dmax1(ERR_CMPL,dabs(V_L(i)*S_L(i)))
         enddo
         do i = 1, NUB
            ERR_CMPL = dmax1(ERR_CMPL,dabs(V_U(i)*S_U(i)-MU))
            ERR_CMPL_0 = dmax1(ERR_CMPL,dabs(V_U(i)*S_U(i)))
         enddo
         ERR_CMPL = ERR_CMPL/SCAL_CMPL
         ERR_CMPL_0 = ERR_CMPL_0/SCAL_CMPL
         ERR_CMPL_US = 0.d0
         if( QSCALE.lt.3 ) then
            do i = 1, NLB
               ERR_CMPL_US = dmax1(ERR_CMPL_US, dabs(V_L(i)*S_L(i)))
            enddo
            do i = 1, NUB
               ERR_CMPL_US = dmax1(ERR_CMPL_US, dabs(V_U(i)*S_U(i)))
            enddo
         else
            do i = 1, NLB
               ERR_CMPL_US = dmax1(ERR_CMPL_US,
     1              dabs(V_L(i)*S_L(i))/CSCALE(M+ILB(i)))
            enddo
            do i = 1, NUB
               ERR_CMPL_US = dmax1(ERR_CMPL_US,
     1              dabs(V_U(i)*S_U(i))/CSCALE(M+IUB(i)))
            enddo
         endif
      endif
C
C     Compute overall scaled error (for the barrier problem)
C
      ERR = dmax1(ERR_DINF, ERR_PINF, ERR_CMPL)
      ERR_0 = dmax1(ERR_DINF_0, ERR_PINF_0, ERR_CMPL_0)
      if( ERR.eq.ERR_DINF ) then
         C_ERR = 'd'
      elseif( ERR.eq.ERR_PINF ) then
         C_ERR = 'p'
      else
         C_ERR = 'c'
      endif

      if( QCNR.gt.0 .and. QPRINT.ge.4 ) then
         write(line,8000) ERR_DINF, ERR_PINF, ERR_CMPL, ERR
 8000    format(/,'Scaled errors for barrier problem...:',
     1        ' DINF = ',d8.2,' PINF = ',d8.2,' CMPL = ',d8.2,
     2        ' ERR = ',d8.2)
         call C_OUT(1,0,2,line)
         write(line,8001) ERR_DINF_0, ERR_PINF_0, ERR_CMPL_0, ERR_0
 8001    format('Scaled errors for original problem..:',
     1        ' DINF = ',d8.2,' PINF = ',d8.2,' CMPL = ',d8.2,
     2        ' ERR = ',d8.2)
         call C_OUT(1,0,1,line)
         write(line,8002) ERR_DINF_US, ERR_PINF_US, ERR_CMPL_US
 8002    format('Unscaled errors for original problem:',
     1        ' DINF = ',d8.2,' PINF = ',d8.2,' CMPL = ',d8.2)
         call C_OUT(1,0,1,line)
      endif

      LASTITER = ITER
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine OPTERROR_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer scaltype, lrw1, liw1, lrw2, liw2
      character*80 line

      if( QERROR.gt.0 ) then
         scaltype = 1
      elseif( QERROR.lt.0 ) then
         scaltype = 0
      else
         scaltype = 2
      endif

      LRW = 0
      LIW = 0

      if( N-M.gt.0 ) then
         if( abs(QERROR).eq.1 ) then
            LRW = max(LRW, N-M)
         else
            if( M.gt.0 ) then
               call CONSTR_WS(N, M, NLB, NUB, NZA, lrw1, liw1,
     1              DAT, IDAT)
            else
               lrw1 = 0
               liw1 = 0
            endif
            call UPDATE_FG_MU_WS(N, M, NLB, NUB, NZA, lrw2, liw2)
            liw1 = max(liw1, liw2)
            lrw1 = max(lrw1, lrw2)
            LRW = max(LRW, N+lrw1)
            LIW = max(LIW, liw1)
         endif
      endif

      if( M.gt.0 .and. QSCALE.ge.3 ) then
         LRW = max(LRW, M)
      endif

      if( QPRINT.ge.4 ) then
         write(line,1000)'opterror_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
