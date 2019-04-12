C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: push_bnd.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine PUSH_BND(TOL_OPT, ZB, YB, X)
C
C     Change the bounds on the variables, if some entries are at or are
C     violating particular bounds.
C
C     Authors:  Yidong Lang, Andreas Waechter   10-02-01
C
      implicit none
C
      double precision TOL_OPT
      double precision ZB(2,*), YB(2,*)
      double precision X(*)
C
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/
C
      integer i, lx, iele, icol
      double precision tmp, delta
      logical decr, incr

      delta = dsqrt(TOL_OPT)

C     Differential variables
      do i = 1, NZ
         lx = i
         decr = .false.
         incr = .false.
         do iele = 1, NELE + 1
            tmp = X(lx)
            if( tmp.le.ZB(1,i) ) decr = .true.
            if( tmp.ge.ZB(2,i) ) incr = .true.
            lx = lx + NCOL*(NZ+NY)
         enddo
         if( decr ) ZB(1,i) = ZB(1,i) - delta
         if( incr ) ZB(2,i) = ZB(2,i) + delta
      enddo

C     Algebraic variables
      do i = 1, NY
         lx = 2*NZ + i
         decr = .false.
         incr = .false.
         do iele = 1, NELE
            do icol = 1, NCOL
               tmp = X(lx)
               if( tmp.le.YB(1,i) ) decr = .true.
               if( tmp.ge.YB(2,i) ) incr = .true.
               lx = lx + (NZ+NY)
            enddo
            lx = lx + NZ
         enddo
         if( decr ) YB(1,i) = YB(1,i) - delta
         if( incr ) YB(2,i) = YB(2,i) + delta
      enddo

C     That's it!

      return
      end
