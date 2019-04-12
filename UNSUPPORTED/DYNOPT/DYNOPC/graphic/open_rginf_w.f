C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      SUBROUTINE open_rginf_W()
 
      USE dflib
 
      implicit none
	include 'DYNAUX.INC'
	include 'DYNOPC.INC'
      include 'DYNGRA.INC'  
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/
      CHARACTER(6) STR
      DOUBLE PRECISION XCORD,ycord,xdigit,ydigit, xt, intvel
      double precision stepx, stepy, xx,yy,yintval
      integer(2)  i2,i,j
      INTEGER(4)  i4
      type (qwinfo) qw
      type(wxycoord)wxy

      open(ypy_zpz_w,file='user',title='Infeasibility vs. Optimality')
      
      i4=FOCUSQQ(ypy_zpz_w)
        qw.x=int2(unit_w*0.5)-2
        qw.y=int2(0)
        qw.h=int2(unit_h)*0.2
        qw.w=int2(unit_w/2)
        qw.type=QWIN$SET
        i4=setwsizeqq(ypy_zpz_w,qw)
      CALL setviewport(0,0,dfloat(qw.w)*ratiox,dfloat(qw.h)*ratioy)
      i4=setwindow(.true.,-0.7d0,-2d0,6.2d0,6.2d0)
      i2=setbkcolor(2)
      call clearscreen($GCLEARSCREEN)
      I2=INITIALIZEFONTS()
      I2=setcolor(0)
      i4=rectangle_w($GFILLINTERIOR,0.d0,0.d0,6.d0,6.d0)
      i2=setcolor(14)
      CALL MOVETO_W(2.5D0,-0.9D0,wxy)
      I2=SETFONT('t''Times new Roman''h14w6pveb')
      call outGtext('Infeasibility')
      call moveto_w(-0.65d0,6.1d0,wxy)
      call outGtext('RdGrd')

! ----- determine coordinate according to tolerance
!
       origin = 100
       do while(  opt_eps .le. origin) 
         origin = origin/10.d0
       end do
         origin = origin/1.d0

        stepx = 10.d0
        stepy = 10.d0

        spanx = 10**(dlog10(origin)+dlog10(stepx)*12.d0)
       spany = 10**(dlog10(origin)+dlog10(stepy)*10.d0)
       factorx = 6.d0/(dlog10(spanx) - dlog10(origin))
       factory = 6.d0/(dlog10(spany) - dlog10(origin))
         
        ydigit=-0.3
	yintval = 6.d0/10.d0
        xcord = origin
 
        intvel=6.d0/12.d0
        xt=intvel
        
        xcord = xcord*stepx
      I2=setcolorrgb(#f0f000)
!     i4=rectangle_w($GFILLINTERIOR,0.d0,0.d0,intvel,6.d0)
      I2=setcolorrgb(#00f0f0)
!      i4=rectangle_w($GFILLINTERIOR,0.d0,0.d0,6.d0,1.d0)
      I2=setcolorrgb(#00ff00)
      i4=rectangle_w($GFILLINTERIOR,0.d0,0.d0,intvel,yintval)


      i2=setfont('t''Times New Roman''h12w6peb')
      do i = 1,11
        xdigit=xt 
        if ( mod(i,2) .eq. 0 ) then 
          xcord=origin*stepx**i
          call moveto_w(xdigit-0.3,ydigit,wxy)
          write(str,'(e6.1e1)')xcord
          i2=setcolor(14)
          call outGtext(str)
        end if
        i2=setcolor(8)
        call moveto_w(xt,0.d0,wxy)
        i2=lineto_w(xt,6d0)
        xt= xt + intvel
        
      end do
      ycord=origin

      do j=0,9
        i2=setcolor(8)
        call moveto_w(0.d0,yintval*j,wxy)
        i2=lineto_w(6.d0,yintval*j)
	  if ((j .eq.0) .or. (mod(j,2) .eq. 0)) then
           call moveto_w(-0.6d0,yintval*j+0.2d0,wxy)
           i2=setcolor(14)
           write(str,'(e6.1e1)')ycord
           call outGtext(str)
        end if
        ycord=ycord*stepy
      end do
      I2=setcolorrgb(#00ff00)
      xx = ((dlog10(opt_eps)-dlog10(origin))*factorx)
      call moveto_w(xx,0.d0,wxy)
      i2 = lineto_w(xx,6.d0)
      yy = ((dlog10(opt_eps)-dlog10(origin))*factory)
      call moveto_w(0.d0,yy,wxy)
      i2 = lineto_w(6.d0,yy)
      return
      END SUBROUTINE 
