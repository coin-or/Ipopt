C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_YPY(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, NLB, ILB, NUB, IUB,
     2     S_L, S_U, C, YPY, NEWBAS, CONDC,
     1     KCONSTR, LRS, RS, LIS, IS,
     2     LRW, RW, LIW, IW, IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: get_ypy.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute range space step
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
CP   N         I    INT    number of variables
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of equality constraints / dependent variables
CP   ITER      I    INT    iteration counter
CP   IVAR      I    INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP   NORIG     I    INT    total number of all variables (incl. fixed vars)
CP   XORIG     I    DP     actual iterate
CP                            XORIG is ordered in ORIGINAL order (i.e. not
CP                            partitioned into independent and dependent
CP                            variables)
CP   CSCALE    I    DP     scaling factors for constraints
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. S_L(i) is slack for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. S_U(i) is slack for X(IUB(i)) )
CP   S_L       I    DP     slacks to lower bounds
CP   S_U       I    DP     slacks to upper bounds
CP   C         I    DP     values of constraint function
CP   YPY       O    DP     range space step (ordered like X)
CP   NEWBAS    O    LOG    =.true.: Basis has become pretty bad; get new one!
CP   CONDC     O    DP     estimated condition number of basis matrix C
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
CP   LRW      I/O   INT    length of RW
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW      I/O   INT    length of IW
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
CS    CONSTR
CS    DSCAL
CS    DCOPY
CS    CHECK_BASIS
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
      include 'TIMER.INC'
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      integer NIND
      integer M
      integer ITER
      integer IVAR(N)
      integer NFIX
      integer IFIX(NFIX)
      integer NORIG
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision C(M)
      double precision YPY(N)
      logical NEWBAS
      double precision CONDC
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
      logical CHECK_BASIS

      double precision times, timef
      integer i, idummy

      integer p_iwend, p_rwend

      character*80 line(3)

      double precision vin(2)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      call TIMER(times)

      IERR = 0
C
C     There is nothing to do if there are no constraints
C
      if( M.eq.0 ) then
         call DCOPY(NIND, 0d0, 0, YPY, 1)
         goto 9999
      endif

      p_iwend = 0
      p_rwend = 0

C
C     Call CONSTR to get -PY and dependent part of BB
C
      call CONSTR(3, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, C, YPY, 1, idummy,
     3     KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4     IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.lt.0 ) then
         write(line,*) 'get_ypy: Warning in CONSTR-3, IERR = ',IERR
         call C_OUT(2,0,1,line)
         IERR = 0
      elseif( IERR.ne.0 ) then
         write(line,*) 'get_ypy: Error in CONSTR-3, IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     Check condition number of basis
C
      NEWBAS = CHECK_BASIS(ITER, M, C, YPY, CONDC)
      if( NEWBAS ) then
         goto 9999
      endif
C
C     PY has wrong sign
C
      call DSCAL(M, -1.d0, YPY, 1)
      if( NIND.gt.0 ) then
         call DCOPY(NIND, 0d0, 0, YPY(M+1), 1)
      endif
C
C     That's it
C
 9999 continue
      call TIMER(timef)
      TIME_YPY = TIME_YPY + timef - times
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine GET_YPY_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)

      LRW = 0
      LIW = 0

      if( M.eq.0 ) return

      call CONSTR_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      return
      end

