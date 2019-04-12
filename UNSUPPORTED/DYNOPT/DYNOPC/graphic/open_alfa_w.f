C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: open_alfa_w.f 531 2004-03-11 01:31:07Z andreasw $
      SUBROUTINE open_alfa_w()

      use dflib

	implicit none
 
	include 'DYNAUX.INC'
	include 'DYNOPC.INC'
      include 'DYNGRA.INC'  
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /GRAPH/

      double precision left_x,down_y,right_x,upper_y
      integer(2)  i2
      INTEGER(4)  i4
      character(7) str
      type (qwinfo) qw
      type(wxycoord)wxy
	integer i

      open(opt_path_w,file='user',title='OPT Iteration Path')
      i4 = INITIALIZEFONTS()
      i4=FOCUSQQ(opt_path_w)
        qw.x = int2(0)
        qw.y = int2(0)
        qw.h = int2(unit_h*0.2)
        qw.w = int2(unit_w/2)-5
        qw.type=QWIN$SET
        i4=setwsizeqq(opt_path_w,qw)
      
      CALL setviewport
     1      (0,0,dfloat(qw.w)*ratiox+5,dfloat(qw.h)*ratioy)
      left_x  =-3.d0
      down_y  =-40.d0
      right_x = 30.d0
      upper_y = 104.d0
      i4=setwindow(.true.,left_x,down_y,right_x,upper_y)
      i4=setcolor(2)
      i4=rectangle_w($GFILLINTERIOR,left_x,
     1      down_y,right_x,upper_y)
      I4=setcolor(0)
      i4=rectangle_w($GFILLINTERIOR,0.D0,0.D0,26.D0,100.D0)
      
      call draw_grids_opt_path()

!     offset=-dlog10(opt_eps)+1
!      jakpot =(dlog10(opt_eps)+offset)*10.d0

!     i4=setcolor(10)
!     CALL moveto_w(0.d0,jakpot,wxy)
!     i4=lineto_w(26.d0,jakpot)

      i2=INITIALIZEFONTS()
      i2=setcolor(14)
      i2=SETFONT('t''Times new Roman''h12w6peb')
      CALL MOVETO_W(9.D0,-20.D0,wxy)
      CALL outGtext('No. of Iterations')
      CALL MOVETO_W(26.3D0,100.D0,wxy)
      i2=setcolor(14)
      CALL outGtext('Alphas')
      i4=SETFONT('t''Times new Roman''h12w6peb')
      do i=0,4
        call moveto_w(26.2d0,2.d0+2*i*10.d0,wxy)
        write(str,'(f4.3)') dfloat(i)*.2
        call outGtext(str)
      end do
!---------------------- For KKT -----------------------------
!     CALL MOVETO_W(26.3D0,100.D0,wxy)
!     i2=setcolor(14)
!     i2=SETFONT('t''Times new Roman''h12w8peb')
!     CALL outGtext('KKT')
!     i4=SETFONT('t''Times new Roman''h12w6peb')
!     str_0 = opt_eps*10
!     do i=1,4
!       call moveto_w(26.2d0,2.d0+2*i*10.d0,wxy)
!       write(str,'(e7.1e1)')str_0
!       call outGtext(str)
!       str_0=str_0*100.d0
!     end do

      return
      END SUBROUTINE open_alfa_w
