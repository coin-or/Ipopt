C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: ipopt_getdata.f 574 2004-04-25 22:42:37Z andreasw $

      subroutine IPOPT_GETDATA(FFUNCS, CFUNCS, NOCG, NORES, NONEGCURV)
      implicit none
      integer FFUNCS, CFUNCS, NOCG, NORES, NONEGCURV
      include 'TIMER.INC'
      include 'IPOPT.INC'
      FFUNCS = FEVALS
      CFUNCS = CEVALS
      NOCG = COUNT_CG
      NORES = COUNT_RESTO_ITER
      NONEGCURV = COUNT_NEG_CURV
C
C     Also close output file if open
C
      if( QCNR.gt.0 ) close(QCNR)
      return
      end

      subroutine CHECK_FLAGFILE( EX )
      integer EX
      logical exist
      inquire(file='FLAG', exist=exist)
      if( exist ) then
         EX = 1
      else
         EX = 0
      endif
      return
      end
