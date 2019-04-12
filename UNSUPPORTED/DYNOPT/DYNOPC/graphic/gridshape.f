C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: gridshape.f 531 2004-03-11 01:31:07Z andreasw $
      SUBROUTINE gridshape()

      USE DFLIB
    
      implicit none
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/
      real(8) xdigit, ydigit
      integer digit, lenstr
      CHARACTER(10) str
      INTEGER * 2 dummy, i, j
      integer * 2 ienhanc, icommon, result
      real(8)   y, tk, tmax
      RECORD/videoconfig/ screen
      RECORD/wxycoord/ wxy
      COMMON screen
      COMMON/MAX/TK, TMAX
C
C    Draw a bordered rectangle around the graph.
c
      dummy = setcolor(12)
      dummy = rectangle_w($Gfillinterior,xx0-.125,yy0-.75,
     &                    xx1+.125,yy1+.75)
      dummy = setcolor(0)
      dummy = rectangle_w($gfillinterior, xx0,yy0,xx1,yy1)   
c
c    Plot the grids.
c
      ienhanc = 7
      icommon = 8
c
c .... Draw vertical lines
c      
      tmax = end_time
      dummy = setcolor(14)
      xdigit=xx0-0.5d0
      Ydigit=yy0-0.5d0
      result = INITIALIZEFONTS()
      result = SETFONT('t''Times New Roman''h14w10peb')
      CALL moveto_w(xdigit,ydigit,wxy)
      CALL outGtext('0')
!
!YDL:  to make the coordinates universal for variety range of end_time 
!
      DO i = 1, 25
        xdigit=dfloat(i)
        IF(MOD(I,5) .eq. 0) THEN
c
c ....  Load digit as characters to str
c
          digit = int(i*(tmax/25))
          WRITE(str,'(i8)') digit
          if(mod(tmax,5.d0) .ne. 0) then
               write(str,'(f8.1)')int(float(i)*tmax/25*100)/100.d0 
          end if
          lenstr = LEN_TRIM(str)
          dummy=setcolor(14)
          call moveto_w(xdigit-0.2d0*lenstr,ydigit,wxy)

          CALL outGtext(trim(str))
          dummy=setcolor(ienhanc)
        ELSE
          dummy = setcolor(icommon)
        END IF
        CALL moveto_w(xdigit, yy0, wxy)
        dummy = lineto_w(xdigit, yy1)
      END DO
C
C ... Draw horizontal lines
c
      y = yy0+2
      DO j = 1,49
        IF(MOD(j,5) .eq. 0) THEN
          dummy=setcolor(14)
          ydigit=y +2.5
          call moveto_w(xx0-1.5,ydigit,wxy)
c
c ....  Load i as a character to str
c
          WRITE(str,'(i2)')int(j/5)
          CALL outGtext(str)
!          call moveto_w(xx1+0.5,ydigit,wxy)
!          write(str,'(i3)') int(j/5)*10
!          call outGtext(str)
          dummy=setcolor(ienhanc)

        ELSE
          dummy = setcolor(icommon)
        END IF
        CALL moveto_w(xx0,y, wxy)
        dummy = lineto_w(xx1,y)
        y = y + 2.
      end do
      return
      END   subroutine gridshape                 
