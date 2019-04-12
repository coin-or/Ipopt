C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      SUBROUTINE open_report_W()
 
      USE dflib
 
      implicit none
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/
      integer(2)     i4  
      integer marginx, marginy
      integer(2) x0,y0,x1,y1
      real(8) xdelt0,xdelt1,ydelt0,ydelt1 
      real(8) upx, upy, downx, downy
      TYPE (qwinfo)  qw
      TYPE (wxycoord)wxy

      open(report_W,file='user',title='Profiles of OPT Solution')
        i4=FOCUSQQ(Report_W)
       
      marginx = 10
      marginy = 5
      qw.x=marginx
      qw.y=marginy
      qw.h=unit_h-2*marginy
      qw.w=unit_w-2*marginx 
      qw.type=QWIN$set
      i4=SETWSIZEQQ(Report_W,qw)

      !--------------------------------------------------
!-- Set x-y coordinate andDraw coordinate bounder
!-------------------------------------------------- 

      x0 = int2(0)+ 15
      y0 = int2(0)+ 15
      x1 = int2((qw.w)*ratiox) - 30
      y1 = int2((qw.h)*ratioy) - 40

      CALL setviewport(x0,y0,x1,y1) !in unit of pixel

      xdelt0 = 2.5
      ydelt0 = 10
      xdelt1 = 2.5
      ydelt1 = 5
      xx0 = 0.d0
      yy0 = 0.d0
      xx1 = 25.d0
      yy1 = 100.d0
      upx   = xx0 - xdelt0
      downx = xx1 + xdelt1
      upy   = yy0 - ydelt0
      downy = yy1 + ydelt1


      i4 = setwindow(.true.,upx,upy,downx,downy )
      i4 = setcolor(0)
      i4 = rectangle_w($Gborder,upx,upy,downx,downy)  
      i4 = setcolor(12)
      i4 = rectangle_w($Gfillinterior,xx0-.125,yy0-.5,
     &                    xx1+.125,yy1+.5)
      i4 = setcolor(0)
      i4 = rectangle_w($gfillinterior, xx0,yy0,xx1,yy1)
             
!---------------------------------------------------------------------
!   Grid coordinate
!--------------------------------------------------------------------- 
      call gridshape()

      i4 = INITIALIZEFONTS()
      i4 = SETFONT('t''Times New Roman''h20w10peb')
      i4 = setcolor(14)
      CALL MOVETO_W(10.d0,-6.d0,wxy)
      CALL outGtext('Length of Horizon')

  
      call legend_W()
      return
 

      END SUBROUTINE open_report_W
