C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_RV(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, NLB, ILB, NUB, IUB, V_L, V_U, RV,
     1     KCONSTR, LRS, RS, LIS, IS,
     2     LRW, RW, LIW, IW, IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: get_rv.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute reduced dual variables (V_U-V_L)
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
CP   XORIG     I    DP     (only TASK = 1,2,3): actual iterate
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
CP   V_L       I    DP     dual variables for lower bounds
CP   V_U       I    DP     dual variables for upper bounds
CP   RV        O    DP     reduced (V_U-V_L)
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
CS    DASV2F
CS    DAXPY
CS    DCOPY
CS    DSCAL
CS    CONSTR
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
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision RV(NIND)
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
      integer p_iwend, p_rwend, p_v, p_nvu
      integer idummy
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
      if( NIND.eq.0 ) then
         goto 9999
      endif

      p_iwend = 0
      p_rwend = 0

      p_v     = p_rwend
      p_nvu   = p_v     + N
      p_rwend = p_nvu   + NUB
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
C     compute full v = V_L-V_U
C
      call DCOPY(NUB, V_U, 1, RW(p_nvu+1), 1)
      call DSCAL(NUB, -1.d0, RW(p_nvu+1), 1)
      call DASV2F(N, NLB, ILB, NUB, IUB,
     1            V_L, RW(p_nvu+1), RW(p_v+1))
      p_rwend = p_nvu

      if( M.ne.0 ) then
C
C     Call CONSTR to get RV for dependent variables
C
         call CONSTR(4, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, RW(p_v+1), RV, idummy, idummy,
     3        KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4        IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.lt.0 ) then
            write(line,*) 'get_rv: Warning in CONSTR, IERR = ',IERR
            call C_OUT(2,0,1,line)
         elseif( IERR.ne.0 ) then
            write(line,*) 'get_rv: Error in CONSTR, IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      else
         call DCOPY(NIND, 0d0, 0, RV, 1)
      endif
C
C     Add independent part
C
      call DAXPY(NIND, -1.d0, RW(p_v+M+1), 1, RV, 1)
C
C     Free work space
C
      p_rwend = p_v
C
C     That's it
C
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine GET_RV_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1

      LRW = 0
      LIW = 0

      if( N-M.eq.0 ) return

      if( M.ne.0 ) then
         call CONSTR_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         LRW = N + max(NUB,lrw1)
         LIW = liw1
      else
         LRW = N+NUB
         LIW = 0
      endif

      return
      end
