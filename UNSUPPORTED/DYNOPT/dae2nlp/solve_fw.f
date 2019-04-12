C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine SOLVE_FW(TRANS, NZ, NY, NCOL, NAC, LA48, NNZ, A,
     1     IRN48, KEEP48, VIN, VOUT, LRW, RW, LIW, IW, IERR)
C
C     $Id: solve_fw.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Solve system involving square part of elemental Jacobian
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      implicit none

      character*(*) TRANS
      integer NZ, NY, NCOL, NAC, LA48, NNZ
      double precision A(*)
      integer IRN48(LA48)
      integer KEEP48(7+10*(NCOL*(NZ+NY)+NAC))
      double precision VIN(NCOL*(NZ+NY)+NAC) !Is changed within this routine!
      double precision VOUT(NCOL*(NZ+NY)+NAC)
      integer LRW               ! need 5*ntot
      double precision RW(LRW)
      integer LIW               ! need ntot
      integer IW(LIW)
      integer IERR

      integer p_rwend, p_iwend, p_iw, p_w
      integer ntot, job
      integer info(12), icntl(9)
      double precision cntl(5), error(3)
      logical ltrans

      p_rwend = 0
      p_iwend = 0
      IERR    = 0
      call MA48ID(cntl, icntl)

      if( TRANS(1:1).eq.'t' .or. TRANS(1:1).eq.'T' ) then
         ltrans = .true.
      else
         ltrans = .false.
      endif

C
C     Change this in case we want iterative refinement
C
      job = 1                   ! It seems that job = 3 often gives nonsense
                                ! if system is ill-conditioned

      ntot = NCOL*(NZ+NY) + NAC

      p_w     = p_rwend
      p_rwend = p_w     + 4*ntot
      if( p_rwend.gt.LRW ) goto 9098
      p_iw    = p_iwend
      p_iwend = p_iw  + ntot
      if( p_iwend.gt.LIW ) goto 9099

      call MA48CD(ntot, ntot, ltrans, job, LA48, A, IRN48,
     1     KEEP48, cntl, icntl, VIN, VOUT, error, RW(p_w+1),
     2     IW(p_iw+1), info)
      if( info(1).ne.0 ) then
         if( info(1).eq.-8 ) then
            IERR = 22
            goto 9999
         endif
c         write(*,*) 'info = ',info(1)
         stop
      endif

      goto 9999
C
C     Error handling
C
 9098 IERR = 98
      goto 9999
 9099 IERR = 99
      goto 9999

 9999 continue
      return
      end
