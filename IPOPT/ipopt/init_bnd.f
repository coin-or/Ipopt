C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine INIT_BND(NORIG, N, NIND, M, NFIX, IFIX, XORIG, X, NLBO,
     1     ILBO, BNDS_LO, V_LO, NLB, ILB, BNDS_L, NUBO, IUBO, BNDS_UO,
     1     V_UO, NUB, IUB, BNDS_U, IVAR, CSCALE, LAM, LRW, RW, LIW, IW,
     1     IERR)
C
C*******************************************************************************
C
C    $Id: init_bnd.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Initialize X, ILB, BNDS_L, V_L, IUB, BNDS_U, V_U, IVAR
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
CA    Andreas Waechter      05/03/02  Added copying of V_LO and V_UO to allow
CA                                      warm start option IMUINIT = 0
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
CP   NORIG     I    INT    total number of variables (incl. fixed vars)
CP   N         I    INT    number of variables (without fixed)
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of constraints
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP                            (assumed to be in increasing order)
CP   XORIG     I    INT    user provided initial point (incl fixed vars)
CP   X         O    INT    copy of XORIG without fixed vars
CP   NLBO      I    INT    number of lower bounds (incl. fixed vars)
CP   ILBO      I    INT    indices of lower bounds (incl. fixed vars)
CP                            (e.g. BNDS_LO(i) is bound for XORIG(ILBO(i)) )
CP   BNDS_LO   I    DP     values of lower bounds (incl. fixed vars)
CP   V_LO     I/O   DP     user provided estimates for multipliers for lower
CP                            bounds (only needed for IMUINIT=0)
CP                            Input : (incl. fixed vars.)
CP                            Output: (excl. fixed vars.)
CP   NLB       O    INT    number of lower bounds (excl. fixed vars)
CP   ILB       O    INT    indices of lower bounds (excl. fixed vars)
CP                            (e.g. BNDS_L(i) is bound for X(ILB(i)) )
CP   BNDS_L    O    DP     values of lower bounds (excl. fixed vars)
CP   NUBO      I    INT    number of upper bounds (incl. fixed vars)
CP   IUBO      I    INT    indices of upper bounds (incl. fixed vars)
CP                            (e.g. BNDS_UO(i) is bound for XORIG(IUBO(i)) )
CP   BNDS_UO   I    DP     values of upper bounds (incl. fixed vars)
CP   V_UO     I/O   DP     user provided estimates for multipliers for upper
CP                            bounds (only needed for IMUINIT=0)
CP                            Input : (incl. fixed vars.)
CP                            Output: (excl. fixed vars.)
CP   NUB       O    INT    number of upper bounds (excl. fixed vars)
CP   IUB       O    INT    indices of upper bounds (excl. fixed vars)
CP                            (e.g. BNDS_U(i) is bound for X(IUB(i)) )
CP   BNDS_U    O    DP     values of upper bounds (excl. fixed vars)
CP   IVAR      O    INT    information about partitioning
CP                            in this case there is no partitioning, i.e.
CP                            are only ordered by inceasing index in XORIG
CP                            Note: fixed variables do not occur in IVAR
CP                            XORIG(IVAR(i)) corresponds to X(i)
CP   CSCALE    I    DP     information for scaling of constraints and variables
CP   LAM      I/O   DP     multipliers for equality constraints
CP                            (will be scaled if QINIT=0)
CP   LRW       I    INT    length of RW
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW       I    INT    length of IW
CP   IW       I/O   INT    can be used as INT work space but content will be
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
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer NORIG
      integer N
      integer NIND
      integer M
      integer NFIX
      integer IFIX(NFIX)
      double precision XORIG(NORIG)
      double precision X(N)
      integer NLBO
      integer ILBO(NLBO)
      double precision BNDS_LO(NLBO)
      double precision V_LO(NLBO)
      integer NLB
      integer ILB(NLB)
      double precision BNDS_L(NLB)
      integer NUBO
      integer IUBO(NUBO)
      double precision BNDS_UO(NUBO)
      double precision V_UO(NUBO)
      integer NUB
      integer IUB(NUB)
      double precision BNDS_U(NUB)
      integer IVAR(N)
      double precision CSCALE(*)
      double precision LAM(M)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            Local variables
C-------------------------------------------------------------------------------
C
      integer i, lfix, lx, llb, lub, ib
      logical fixed
      integer p_ilbo1, p_iubo1, p_vcopy, p_iwend, p_rwend
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
      p_iwend = 0
      p_rwend = 0

      p_ilbo1 = p_iwend
      p_iubo1 = p_ilbo1 + NORIG
      p_iwend = p_iubo1 + NORIG
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Compute inverse of ILBO, LUBO
C
      do i = 1, NORIG
         IW(p_ilbo1+i) = 0
         IW(p_iubo1+i) = 0
      enddo
      do i = 1, NLBO
         IW(p_ilbo1+ILBO(i)) = i
      enddo
      do i = 1, NUBO
         IW(p_iubo1+IUBO(i)) = i
      enddo

      lfix = 1
      lx = 1
      llb = 1
      lub = 1
      do i = 1, NORIG
         fixed = .false.
         if( lfix.le.NFIX ) then
            if( IFIX(lfix).eq.i ) then
               lfix = lfix + 1
               fixed = .true.
            endif
         endif
         if( .not.fixed ) then
            ib = IW(p_ilbo1+i)
            if( ib.ne.0 ) then
               ILB   (llb) = lx
               BNDS_L(llb) = BNDS_LO(ib)
               if( QSCALE.ge.3 ) BNDS_L(llb) = BNDS_L(llb)/CSCALE(M+lx)
               llb = llb + 1
            endif
            ib = IW(p_iubo1+i)
            if( ib.ne.0 ) then
               IUB   (lub) = lx
               BNDS_U(lub) = BNDS_UO(ib)
               if( QSCALE.ge.3 ) BNDS_U(lub) = BNDS_U(lub)/CSCALE(M+lx)
               lub = lub + 1
            endif
            X   (lx) = XORIG(i)
            if( QSCALE.ge.3 ) X(lx) = X(lx)/CSCALE(M+lx)
            IVAR(lx) = i
                 lx  = lx + 1
            if( lx.gt.N ) then
               goto 100
            endif
         endif
      enddo
 100  continue
C
C     Reorder the entries in V_LO and V_UO so that they correspond to BNDS_L
C     and BNDS_U (only required for warm start option IMUINIT = 0).  Also,
C     rescale multipliers if necessary
C
      if( QINIT.eq.0 ) then
         p_vcopy = p_rwend
         p_rwend = p_vcopy + max(NLBO,NUBO)
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         call DCOPY(NLBO, V_LO, 1, RW(p_vcopy+1), 1)
         if( QSCALE.le.2 ) then
            do i = 1, NLB
               V_LO(i) = RW(p_vcopy+IW(p_ilbo1+IVAR(ILB(i))))*QFSCALE
            enddo
         else
            do i = 1, NLB
               lx = IVAR(ILB(i))
               V_LO(i) = RW(p_vcopy+IW(p_ilbo1+lx))*CSCALE(M+lx)*QFSCALE
            enddo
         endif
         call DCOPY(NUBO, V_UO, 1, RW(p_vcopy+1), 1)
         if( QSCALE.le.2 ) then
            do i = 1, NUB
               V_UO(i) = RW(p_vcopy+IW(p_iubo1+IVAR(IUB(i))))*QFSCALE
            enddo
         else
            lx = IVAR(IUB(i))
            do i = 1, NUB
               V_UO(i) = RW(p_vcopy+IW(p_iubo1+lx))*CSCALE(M+lx)*QFSCALE
            enddo
         endif
         p_rwend = p_vcopy
         if( QSCALE.lt.2 ) then
            call DSCAL(M, QFSCALE/CSCALE(1), LAM, 1)
         else
            do i = 1, M
               LAM(i) = LAM(i)*QFSCALE/CSCALE(i)
            enddo
         endif
      endif
C
C     END
C
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine INIT_BND_WS(N, M, NLB, NUB, NZA, LRW, LIW)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      character*80 line

      LRW = 0
      LIW = 2*N

      if( QPRINT.ge.4 ) then
         write(line,1000)'init_bnd_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
