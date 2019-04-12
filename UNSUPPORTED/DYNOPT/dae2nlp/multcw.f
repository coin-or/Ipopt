C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine MULTCW(TRANS, PLUS, NZ, NY, NU, NP, NCOL, H,
     1     NNZCW, NNZCWW, IRCW, JCCW, VIN, VOUT)
      implicit none
C
C     $Id: multcw.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Compute product of vector with Jacobian corresponding to continuity
C     equation.
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      character*(*) TRANS       !if = 't' or 'T', compute vout += CW^T*vin
                                !otherwise vout += CW*vin
                                !NOTE: VOUT is NOT INITIALIZED here!!!
      logical PLUS              !if .true., add to VOUT, otherwise substract
      integer NZ, NY, NU, NP, NCOL
      double precision H
      integer NNZCW
      integer NNZCWW(NCOL)
      integer IRCW(NNZCW)
      integer JCCW(NNZCW)
      double precision VIN(*)
      double precision VOUT(*)

      include 'DAE2NLP.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DAENLP/

      integer i, k, ir, jc, lcw
      double precision homega

      lcw = 0

      if( TRANS(1:1).eq.'t' .or. TRANS(1:1).eq.'T' ) then
         if( PLUS )then
            do k = 1, NCOL
               homega = h*OMEGA1(k)
               do i = 1, NNZCWW(k)
                  ir = JCCW(lcw+i)
                  jc = IRCW(lcw+i)
                  VOUT(ir) = VOUT(ir) - homega*VIN(jc)
               enddo
               lcw = lcw + NNZCWW(k)
            enddo
         else
            do k = 1, NCOL
               homega = h*OMEGA1(k)
               do i = 1, NNZCWW(k)
                  ir = JCCW(lcw+i)
                  jc = IRCW(lcw+i)
                  VOUT(ir) = VOUT(ir) + homega*VIN(jc)
               enddo
               lcw = lcw + NNZCWW(k)
            enddo
         endif
      else
         if( PLUS ) then
            do k = 1, NCOL
               homega = h*OMEGA1(k)
               do i = 1, NNZCWW(k)
                  ir = IRCW(lcw+i)
                  jc = JCCW(lcw+i)
                  VOUT(ir) = VOUT(ir) - homega*VIN(jc)
               enddo
               lcw = lcw + NNZCWW(k)
            enddo
         else
            do k = 1, NCOL
               homega = h*OMEGA1(k)
               do i = 1, NNZCWW(k)
                  ir = IRCW(lcw+i)
                  jc = JCCW(lcw+i)
                  VOUT(ir) = VOUT(ir) + homega*VIN(jc)
               enddo
               lcw = lcw + NNZCWW(k)
            enddo
         endif
      endif

      return
      end
