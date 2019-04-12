C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine FACT_FW(TASK, NZ, NY, NCOL, NAC, LA48, NNZ, A, IRA,
     1     JCA, A48, IRN48, JCN48, KEEP48, LRW, RW, LIW, IW, IERR)
C
C     $Id: fact_fw.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Factorize basis matrix for one finite element
C     
C     TASK   1: call MA48AD to prepare datastructures
C               in:  LA48, NNZ, A(NNZ), IRA(NNZ), JCA(NNZ) (all unchanged)
C               out: IRN48(LA48), JCN48(NNZ), KEEP48(*)
C            2: call MA48BD with job = 1 (i.e. KEEP and IRN can change)
C               in:  LA48, NNZ, A(NNZ), IRN48(LA48), JCN48(NNZ), KEEP48()
C               out: A48(LA48), IRN48(LA48), KEEP48()
C            3: call MA48BD with job = 2 (fast call), doesn't change KEEP, IRN
C               in:  LA48, NNZ, A(NNZ), IRN48(LA48), JCN48(NNZ), KEEP48()
C               out: A48(LA48)
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      implicit none

      integer TASK, NZ, NY, NCOL, NAC, LA48, NNZ
      double precision A(NNZ)
      integer IRA(NNZ)
      integer JCA(NNZ)
      double precision A48(LA48)
      integer IRN48(LA48)
      integer JCN48(NNZ)
      integer KEEP48(7+10*NCOL*(NZ+NY))
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR

      integer p_rwend, p_iwend, p_a, p_iw, p_w, p_jcn
      integer ntot, job, i
      integer info(12), icntl(9)
      double precision rinfo, cntl(5)

      p_rwend = 0
      p_iwend = 0
      IERR    = 0
      call MA48ID(cntl, icntl)
      icntl(2) = 0
C
C     The following statement changes the pivot tolerance
C
c      cntl(2)  = 1.d-2

      goto (1000, 2000, 3000) TASK

c      write(*,*) 'wrong TASK = ',TASK
      STOP
C
C     prepare data structures for MA48AD
C
 1000 continue

      job = 1
      ntot = NCOL*(NZ+NY) + NAC

      p_a     = p_rwend
      p_rwend = p_a + LA48
      if( p_rwend.gt.LRW ) goto 9098
      call DCOPY(NNZ, A, 1, RW(p_a+1), 1)

      p_iw    = p_iwend
      p_jcn   = p_iw  + 9*ntot
      p_iwend = p_jcn + LA48
      if( p_iwend.gt.LIW ) goto 9099

      do i = 1, NNZ
         IRN48(i)    = IRA(i)
         IW(p_jcn+i) = JCA(i)
      enddo

      call MA48AD(ntot, ntot, NNZ, job, LA48, RW(p_a+1), IRN48,
     1     IW(p_jcn+1), KEEP48, cntl, icntl, IW(p_iw+1), info, rinfo)
      if( info(1).ne.0 ) then
         if( info(1).eq.-4 .or. info(1).eq.2 ) then
            IERR = 20           ! system singular
            goto 9999
         endif
         if( info(1).eq.-3 ) then
            IERR = 21           ! LA48 too small
            goto 9999
         endif
         if( info(1).ne.1 .or. info(12).ne.0 ) then
c            write(*,*) 'info = ',info(1)
            IERR = 48
            goto 9999
         endif
      endif

      do i = 1, NNZ
         JCN48(i) = IW(p_jcn+i)
      enddo

      p_rwend = p_a
      p_iwend = p_iw

      goto 9999
C
C     Do normal factorization
C
 2000 continue
      job = 1
      goto 3500
C
C     Do fast factorization
C
 3000 continue
      job = 2

 3500 continue
      ntot = NCOL*(NZ+NY) + NAC

      p_w     = p_rwend
      p_rwend = p_w + ntot
      if( p_rwend.gt.LRW ) goto 9098

      p_iw    = p_iwend
      p_iwend = p_iw + 4*ntot
      if( p_iwend.gt.LIW ) goto 9099

      call DCOPY(NNZ, A, 1, A48, 1)

      call MA48BD(ntot, ntot, NNZ, job, LA48, A48, IRN48, JCN48, KEEP48,
     1     cntl, icntl, RW(p_w+1), IW(p_iw+1), info, rinfo)
      if( info(1).ne.0 ) then
         if( info(1).eq.-3 ) then
            IERR = 21
            goto 9999
         endif
c         write(*,*) 'info = ',info(1)
         IERR = 48
         goto 9999
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
