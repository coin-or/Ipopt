C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: set_unit.f 531 2004-03-11 01:31:07Z andreasw $
!*************************************************************************
! 
!   These subroutine are for display the results from numerical integration
!     and then discretisation. in order to let user have knowledge about 
!    goodness of his/her element placement 
!  
!       YDL
!             third modified : 03-13-2001 
!..........................................................................

      Subroutine Set_Unit()
!----------------------------------------------------------------
!  Purpose: calculate unit of width and height for child windows
!----------------------------------------------------------------
      Use dflib
      implicit none
      include 'DYNAUX.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /GRAPH/
      integer(2) i4
      logical(4) status 
      type (fontinfo) finfo
      type (qwinfo) qw, cqw
      type (windowconfig) wc
!  
!
!-- set highest resolution possible for the system by -1 sets

      i4 = clickqq(qwin$status)  ! turn off status bar

      wc.numxpixels = -1
      wc.numypixels = -1
      wc.numtextcols = -1
      wc.numtextrows = -1
      wc.numcolors = -1
      wc.fontsize = -1
      wc.title = " Dynamic Optimization "C 
      status = setwindowconfig(wc)
      if (.not. status) status = setwindowconfig(wc)
      
C-- Maximize fram window

      qw.TYPE = QWIN$MAX  
       
      i4 = SETWSIZEQQ(QWIN$FRAMEWINDOW,qw)
      i4 = getwsizeqq(qwin$framewindow,QWin$sizemax,qw)
      i4 = ABOUTBOXQQ('DynoPC Version 4.0'C)

C-- in order to measure current window parameters, Open a child window 
C     to get maximum dimensions then delete it

      open(10,file = 'user', title = '')
      cqw.type = QWIN$max
      cqw.w = qw.w
      cqw.h = qw.h
      I4 = SETWSIZEQQ(10,cqw)
      i4 = getwsizeqq(10,QWIN$SIZEcurr,cqw)
      close(10, status='delete')
      
      ratiox = dfloat(qw.w)/dfloat(cqw.w)
      ratioy = dfloat(qw.h)/dfloat(cqw.h)

!  CALL clearscreen($GCLEARSCREEN)
 
      unit_w = cqw.w
      unit_h = cqw.h
      return
      end
