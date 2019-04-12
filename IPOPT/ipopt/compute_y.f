C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine COMPUTE_Y(N, NIND, M, ITER, IVAR, NFIX, IFIX, NORIG,
     1     XORIG, CSCALE,  RG, RGOLD, WBAR, ALPHA, G, GOLD,
     2     LAM, Y, KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR,
     3     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     4     DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: compute_y.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    compute 'y' for BFGS formula (for now no correction!)
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
CP   RG        I    DP     reduced gradient of objective function
CP                            (only for QQUASI > 0)
CP   RGOLD     I    DP     reduced gradient of objective function from
CP                            last iteration
CP                            (only for QQUASI > 0)
CP   WBAR      I    DP     correction term (not for now!)
CP   ALPHA     I    DP     factor for correction term (not for now!)
CP   G         I    DP     actual gradient of objective function
CP                            (only for QQUASI < 0)
CP   GOLD      I    DP     gradient of objective function from last iteration
CP                            (only for QQUASI < 0)
CP   LAM       I    DP     Lagrangian multipliers
CP                            (only for QQUASI < 0)
CP   Y         O    DP     'y' for BFGS formula
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
CS    DCOPY
CS    DAXPY
CS    C_OUT
CS    GET_RG
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
      double precision RG(NIND)
      double precision RGOLD(NIND)
      double precision WBAR(NIND)
      double precision ALPHA
      double precision G(N)
      double precision GOLD(N)
      double precision LAM(M)
      double precision Y(NIND)
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
      integer p_rwend, p_iwend, p_diff, idummy
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      p_rwend = 0
      p_iwend = 0
      IERR = 0
C
      if( QQUASI.gt.0 ) then

         call DCOPY(NIND, RG, 1, Y, 1)
         call DAXPY(NIND, -1.d0, RGOLD, 1, Y, 1)

      elseif( QQUASI.lt.0 ) then

         if( QLAMBDA.eq.0 ) then
            call C_OUT(2,0,1,'Need multipliers for QQUASI<0!')
            IERR = 4
            goto 9999
         endif
C
         p_diff  = p_rwend
         p_rwend = p_diff + N
         if( p_rwend.gt.LRW ) then
            IERR = 99
            goto 9999
         endif

C
C     get A_k * lam_{k+1}
C
         if( M.gt.0 ) then
            call CONSTR(12, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, LAM, RW(p_diff+1), idummy,
     2           idummy, KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.gt.0 ) then
               write(line,*)
     1              'compute_y: constr-12 returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            elseif( IERR.lt.0 ) then
               write(line,*)
     1              'compute_y: constr-12 returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               IERR = 0
            endif
            call DSCAL(N, -1d0, RW(p_diff+1), 1)
            call DAXPY(N, 1d0, G, 1, RW(p_diff+1), 1)
         else
            call DCOPY(N, G, 1, RW(p_diff+1), 1)
         endif
         call DAXPY(N, -1d0, GOLD, 1, RW(p_diff+1), 1)
C
C     Abuse GET_RG to obtain 'reduced' A * lam
C
         call GET_RG(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, RW(p_diff+1), Y,
     1        KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            write(line,*)
     1           'compute_y: Error: get_rg ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'compute_y: Warning: get_rg ends with IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
C
         p_rwend = p_diff
C
C     That should be it
C
      endif
C
C     This has to be included if correction present:
C     call DAXPY(NIND, ALPHA, WBAR, 1, Y, 1)

 9999 continue
      return
      end
C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine COMPUTE_Y_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1, lrw2, liw2, lrw3, liw3

      LRW = 0
      LIW = 0

      if( QQUASI.lt.0 ) then

         lrw1 = N               ! p_diff
         if( M.gt.0 ) then
            call CONSTR_WS(N, M, NLB, NUB, NZA, lrw2, liw2, DAT, IDAT)
         else
            lrw2 = 0
            liw2 = 0
         endif
         call GET_RG_WS(N, M, NLB, NUB, NZA, lrw3, liw3, DAT, IDAT)
         lrw1 = lrw1 + max(lrw2, lrw3)
         liw1 = max(liw2, liw3)

         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)

      endif

      return
      end
