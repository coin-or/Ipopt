C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: get_iele4t.f 531 2004-03-11 01:31:07Z andreasw $
      integer function GET_IELE4T(T)
C
C     Determine which element is responsible for time T
C     (appears literally in appsln.f)
C
C     Author:   Andreas Waechter     10-01-01
C
      implicit none
      double precision T

      include 'DYNAUX.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/

      integer iele, il, iu
C
C     Find element in which T is located (binary search)
C
      if( TI(2).ge.T ) then
         iele = 1
      elseif( TI(NELE).lt.T ) then
         iele = NELE
      else
         il = 2
         iu = NELE
 100     continue
         if( iu-il.gt.1 ) then
            iele = il+(iu-il)/2
            if( TI(iele).lt.T ) then
               il = IELE
            else
               iu = iele
            endif
            goto 100
         endif
         iele = il
      endif

      GET_IELE4T = iele
      return
      end
