C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine NEW_BASIS(ITER, N, M, NIND, NORIG, X, XORIG, CSCALE,
     1     NFIX, IFIX, IVAR, NLB, ILB, NUB, IUB, S_L, S_U, C, LAM,
     2     LAMOLD, YPY, ALPHA, DX, PZ, G, RGOLD, B, SKIP_UPDATE,
     3     KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR, EV_F,
     4     EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: new_basis.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Get new partiton, reorder bound pointers and X and
CT    "update" Quasi-Newton matrix
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
CP
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
CS    PARTIRION
CS    REORDER_IB
CS    REORDER_X
CS    GET_YPY
CS    GET_RG
CS    CONSTR
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
      integer M
      integer NIND
      integer NORIG
      double precision X(N)
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      integer NFIX
      integer IFIX(NFIX)
      integer IVAR(N)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision C(M)
      double precision LAM(M+N)
      double precision LAMOLD(M+N)
      double precision YPY(N)
      double precision ALPHA
      double precision DX(N)
      double precision PZ(NIND)
      double precision G(N)
      double precision RGOLD(NIND)
      double precision B((NIND*(NIND+1))/2)
      logical SKIP_UPDATE
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
      integer p_rwend, p_iwend, p_ivarold, p_ivar1, p_ib, p_xtmp
      integer p_if, p_r, p_rtmp, p_xold, p_ypydummy
      integer p_xorigtmp

      integer i, j, k, lb, idummy
      double precision dummy, ele, condctmp
      logical newbastmp
      character*80 line

      logical LEFTINVERSE       ! if set to true, then we update QN matrix
                                ! according to repartitioning.  otherwise,
                                ! reset QN estimate to I.
      parameter( LEFTINVERSE = .false. )
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

CTODO Think about this!:
      SKIP_UPDATE = .true.

      if( abs(QQUASI).ge.6 ) then
         write(line,*) 'new_basis: QQUASI = ',QQUASI,' not implemented.'
         call C_OUT(2,0,1,line)
         IERR = 4
         goto 9999
      endif

      p_rwend = 0
      p_iwend = 0
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
C
C     Determine new partition IVAR
C
      call PARTITION(0, ITER, N, NIND, M, NFIX, IFIX, IVAR,
     1     NORIG, XORIG, CSCALE, KCONSTR, LRS, RS, LIS, IS,
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR,
     2     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     4     DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1        'new_basis: Error: partition ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.ne.0 ) then
         write(line,*)
     1        'new_basis: Warning: partition ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
      endif
C
C     Reorder bound pointers ILB, IUB
C     Also, compute inverse of IVAR (new)
C
      p_ivar1 = p_iwend
      p_ib    = p_ivar1 + NORIG
      p_iwend = p_ib       + N
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      call REORDER_IB(N, NORIG, NLB, ILB, NUB, IUB,
     1                IW(p_ivarold+1), IVAR, IW(p_ivar1+1),
     2                IW(p_ib+1))
      p_iwend = p_ib
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
     1               RW(p_xtmp+1))
      if( QSCALE.ge.3 ) then
         call REORDER_X(N, NORIG, CSCALE(M+1), IW(p_ivarold+1), IVAR,
     1               RW(p_xtmp+1))
      endif
C
C     The following instructions are only necessary if update will take place:
C
      if( SKIP_UPDATE ) then
         p_rwend = p_xtmp
         goto 1000
      endif
C
C     Reorder DX
C
      call REORDER_X(N, NORIG, DX, IW(p_ivarold+1), IVAR,
     1               RW(p_xtmp+1))
C
      p_rwend = p_xtmp
C
C     Get X_old = X - ALPHA*DX
C
      p_xold  = p_rwend
      p_rwend = p_xold + N
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call DCOPY(N, X, 1, RW(p_xold+1), 1)
      call DAXPY(N, -1.d0*ALPHA, DX, 1, RW(p_xold+1), 1)
C
C     factorize C at old point according to new partition and store CinvN
CAW   MIGHT NOT NEED GET_YPY WITH NEW CONSTR!
C
      p_ypydummy = p_rwend
      p_xorigtmp = p_ypydummy  + N
      p_rwend    = p_xorigtmp + NORIG
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      if( NIND.gt.0 ) then
         call DCOPY(NORIG, XORIG, 1, RW(p_xorigtmp+1), 1)
      endif
      do i = 1, N
         RW(p_xorigtmp+IVAR(i)) = RW(p_xold+i)
      enddo
      call GET_YPY(N, NIND, M, -1, IVAR, NFIX, IFIX,
     1     NORIG, RW(p_xorigtmp+1), CSCALE, NLB, ILB, NUB, IUB,
     2     S_L, S_U, C, RW(p_ypydummy+1), newbastmp,
     1     condctmp, KCONSTR, LRS, RS, LIS, IS,
     2     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     2     IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     6     DAT, IDAT)
      p_rwend = p_xold
C
C     Compute RGOLD (G has still values from old iteration)
C
      call GET_RG(N, NIND, M, -1, IVAR, NFIX, IFIX,
     1     NORIG, RW(p_xorigtmp+1), CSCALE, G, RGOLD,
     1     KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1      'new_basis: Error: get_rg ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.ne.0 ) then
         write(line,*)
     1      'new_basis: Warning: get_rg ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
      endif
C
C     Store PZ in terms of new partition
C
      call DCOPY(NIND, DX(M+1), 1, PZ, 1)
C
 1000 continue

C
C     Factorize C according to new partition, compute CinvN and get
C     range space step according to new partition
C
      call GET_YPY(N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, NLB, ILB, NUB, IUB,
     2     S_L, S_U, C, YPY, newbastmp, condctmp,
     1     KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1      'new_basis: Error: get_py ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.ne.0 ) then
         write(line,*)
     1      'new_basis: Warning: get_py ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
      endif
C
C     Take care of Quasi-Newton matrix (see p. 70 in Claudia's thesis)
C
      if( QQUASI.eq.0 ) goto 2000   ! only need to do this if not exact derivs

      if( LEFTINVERSE ) then
                                ! only need to do this if not exact derivs
C
C     Compute 'filter' for CinvN
C
         p_if    = p_iwend
         p_iwend = p_if    + NIND
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         do i = 1, NIND
            k = IW(p_ivar1+IW(p_ivarold+M+i))
            if( k.le.M ) then
               IW(p_if+i) = k
            else
               IW(p_if+i) = 0
            endif
         enddo
C
C     Obtain R = PNOLD'*PC*CinvN from CONSTR
C
         p_r     = p_rwend
         p_rwend = p_r     + NIND*NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         call DCOPY(NIND*NIND, 0.d0, 0, RW(p_r+1), 1)
         call CONSTR(7, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1        NORIG, XORIG, CSCALE, dummy, RW(p_r+1), IW(p_if+1),
     2        idummy, KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4        IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         p_iwend = p_if
         if( IERR.gt.0 ) then
            write(line,*)
     1           'new_basis: Error: constr returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         elseif( IERR.ne.0 ) then
            write(line,*)
     1           'new_basis: Warning: constr returns IERR = ',IERR
            call C_OUT(2,0,1,line)
         endif
C
C     Copute R = R - PNOLD'*PN
C
         do i = 1, NIND
            k = IW(p_ivar1+IW(p_ivarold+M+i))
            if( k.gt.M ) then
               RW(p_r+i+NIND*(k-M-1)) = RW(p_r+i+NIND*(k-M-1)) - 1.d0
            endif
         enddo
         p_iwend = p_ivarold
C
C     Now get RTMP = B*R
C
         p_rtmp  = p_rwend
         p_rwend = p_rtmp  + NIND*NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif

         do i = 1, NIND
            do j = 1, NIND
               ele = 0.d0
               lb = ((i-1)*i)/2
               do k = 1, i
                  lb = lb + 1
                  ele = ele + B(lb)*RW(p_r+k+NIND*(j-1))
               enddo
               do k = i+1, NIND
                  lb = lb + k - 1
                  ele = ele + B(lb)*RW(p_r+k+NIND*(j-1))
               enddo
               RW(p_rtmp+i+NIND*(j-1)) = ele
            enddo
         enddo
C
C     Finally get B(new) = R'*RTMP (only upper part)
C
         lb = 0
         do i = 1, NIND
            do j = 1, i
               lb = lb + 1
               ele = 0.d0
               do k = 1, NIND
                  ele = ele +
     1                 RW(p_r+k+NIND*(i-1))*RW(p_rtmp+k+NIND*(j-1))
               enddo
               B(lb) = ele
            enddo
         enddo
C
         p_rwend = p_r
      else
C
C     Just reinitialize B to identity
C
         call DCOPY( (NIND*(NIND+1))/2, 0.d0, 0, B, 1)
         k = 0
         do i = 1,NIND
            k = k + i
            B(k) = 1.0d0
         enddo
      endif
C
C     I hope that's it...
C
 2000 continue

 9999 continue

      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine NEW_BASIS_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1, lrw2, liw2

      LRW = 0
      LIW = N

      call PARTITION_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
      liw1 = max(liw1,2*N)
      lrw1 = max(lrw1,N)

      call GET_YPY_WS(N, M, NLB, NUB, NZA, lrw2, liw2, DAT, IDAT)
      liw1 = max(liw1, N+liw2)
      lrw1 = max(lrw1, lrw2)

      return
      end
