C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: legend_w.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine Legend_W()
 
      USE DFLIB
      implicit none
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/

      real(8)             xtext,ytext,yspace
      character*8         text
      character*4         txt
      INTEGER(2)          result, nvg
      integer(4)          icolor(12)
      integer             i
 
      TYPE (wxycoord)     wxy
      common/color/icolor

C 
c ...  Text discriptions on the graphics
c
 ! set color index for profiles

      icolor(1) = 10  
      icolor(2) = 11  
      icolor(3) = 12 
      icolor(4) = 13     
      icolor(5) = 14  
      icolor(6) = 15   
      icolor(7) = 9      
      icolor(8) = 12
      icolor(9) = 10
      icolor(10) = 11
      icolor(11) = 13
      icolor(12) = 14
      result=INITIALIZEFONTS()
      result=SETFONT('t''Arial''h25w15pvib') 
c
c ... make a color index on the legend
c
      xtext=xx1+0.2d0
      ytext=yy0+1.d0
      yspace=4.d0
      result = INITIALIZEFONTS()
      result=SETFONT('t''Times New Roman''h12w8peb')

      do i = 1,ndsplz
        write(txt,'(I3)') indexz(i)
        text =' Z('//txt//')' 
        result = setcolor(icolor(i))
        CALL MOVETO_W(xtext,(ytext+i*yspace), wxy)
        CALL outGtext(text)
      end do

      do i = 1, ndsplu
        write(txt,'(I2)') indexu(i)
        text =' u('//txt//')' 
        result = setcolor(icolor(ndsplz+i))
        CALL MOVETO_W(xtext,(ytext+(ndsplz+i)*yspace), wxy)
        CALL outGtext(text)
      end do
      if (ndsply .ne. 0) then
        nvg = ndsplz+ndsplu
        do i = 1, ndsply
          write(txt,'(I4)') indexy(i)
             text =' Y('//txt//')' 
          result = setcolor(icolor(nvg+i))
          CALL MOVETO_W(xtext,(ytext+(nvg+i)*yspace), wxy)
          CALL outGtext(text)
        end do
      end if
      return
      END SUBROUTINE 
