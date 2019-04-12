C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: draw_grids_opt_path.f 531 2004-03-11 01:31:07Z andreasw $
      SUBROUTINE draw_grids_opt_path()
 
      USE dflib

	implicit none
	include 'DYNAUX.INC'
	include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/

      real(8)             xdigit,ydigit
      CHARACTER(4)           str
      INTEGER(2) dummy, i, j
      integer(2) ienhanc, icommon
      real(8)  x, y, tintvel
      RECORD/videoconfig/ screen
      RECORD/wxycoord/ wxy
      COMMON screen
    

c    Plot the grids.
c
        ienhanc = 7
        icommon = 8
c
c .... Draw vertical lines
c
       dummy=INITIALIZEFONTS()
        tintvel = 1 
        x = tintvel
        dummy = setcolor(14)
        xdigit=0.d0-0.5d0
        Ydigit=-1.5d0
      DUMMY=SETFONT('t''Times new Roman''h12w8pveb')
      CALL moveto_w(xdigit,ydigit,wxy)
      if (i_25 .ge. 1) then
         write(str,'(i4)') i_25*25
      else
         write(str,'(i4)') 0
      end if
      CALL outGtext(str)   
      DO i = 1, 25
        IF(MOD(I,5) .eq. 0) THEN
          dummy=setcolor(14)
          xdigit=x
            call moveto_w(xdigit-.8d0,ydigit,wxy)
c
c ....  Load i as a character to str
c
          WRITE(str,'(i4)')i_25*25+i
          CALL outGtext(str)
          dummy=setcolor(ienhanc)
        ELSE
          dummy = setcolor(icommon)
        END IF
        CALL moveto_w(x, 0.d0, wxy)
        dummy = lineto_w(x, 100.d0)
        x=x+tintvel
      END DO
C
C ... Draw vertical lines
c
      y = 5.d0
      DO j = 1,19
        IF(MOD(j,4) .eq. 0) THEN
          dummy=setcolor(ienhanc)
        ELSE
          dummy = setcolor(icommon)
        END IF
        CALL moveto_w(0.d0,y, wxy)
        dummy = lineto_w(26,y)
        y = y + 5.d0
      end do
      
      END   subroutine draw_grids_opt_path  
