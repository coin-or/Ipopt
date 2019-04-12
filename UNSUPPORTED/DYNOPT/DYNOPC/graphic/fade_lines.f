C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine fade_lines(pre_x, pre_y)
!--------------------------------------------------------------------------------------
!
!        Purpose : fade previous profile Lines
!                  YDL 03-28-2001
!---------------------------------------------------------------------------------------
      USE dflib
      
      implicit none

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/

      double precision pre_x(0:NELEMAX*NCOLMAX)
      double precision pre_y(0:NELEMAX*NCOLMAX,2*NDSPLMAX) ! for previous profiles

      integer dummy
      integer i, j
      TYPE(wxycoord) wxy

      dummy = setcolorRGB(#446644)
C     pre_x(0) = 0 

      do i = 1, nele*ncol+nele
         do j = 1,ndsplz+ndsply+ndsplu 
            call moveto_w(pre_x(i-1),pre_y(i-1,j),wxy)
            dummy = lineto_w(pre_x(i),pre_y(i,j))
         end do
      end do

      return
      end
