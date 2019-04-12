C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: draw_rg_inf_path.f 531 2004-03-11 01:31:07Z andreasw $
      SUBROUTINE draw_rg_inf_path(rg0,infsbl0, mu0)
 
      USE dflib

      implicit none

	include 'DYNAUX.INC'
	include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/

      save
      INTEGER(2) I2
      integer(4) i4
      INTEGER  FIRST, iclr
      double precision rg, infsbl, mu, mu_old
      double precision rg0,infsbl0, mu0
 
      double precision ypy_0, zpz_0, ypy_1,zpz_1, ypy, zpz
      double precision wx1,wx2,wy1,wy2, deltx, delty
      type(wxycoord)  wxy
c     if( first .le. 1 )goto 10
! make the profile inside of the region

       rg = rg0
       infsbl = infsbl0
       mu = mu0

       if (rg .gt. spany)  rg = spany
       if (rg .lt. origin) rg = origin
       if (infsbl .gt. spanx) infsbl = spanx
       if (infsbl .lt. origin) infsbl = origin

      i4 = focusqq(ypy_zpz_w)
      i2=setcolor(2)
      I2=ELLIPSE_w($GFILLINTERIOR,wx1,wy1,wx2,wy2)
      if(first .gt. 1)  then
        CALL MOVETO_W(zpz_1,ypy_1,wxy)
        i4=setcolor(14)
        i4=LINETO_W(zpz_0,ypy_0)
      end if
10    ypy = (dlog10(rg)-dlog10(origin))*factory
      zpz = (dlog10(infsbl)-dlog10(origin))*factorx
      CALL MOVETO_W(zpz_0,ypy_0,wxy)
      if(first .ne. 0)then
         i4=setcolorrgb(#0000ff)
         i4=LINETO_W(zpz,ypy)
      end if
      deltx = 0.05d0
      delty = 0.1d0
      wx1=zpz-deltx
      wy1=ypy+delty
      wx2=zpz+deltx
      wy2=ypy-delty
      i2=setcolor(11)
      I2=ELLIPSE_w($GFILLINTERIOR,wx1,wy1,wx2,wy2)
      ypy_1 = ypy_0
      zpz_1 = zpz_0
20    ypy_0 = ypy
      zpz_0 = zpz
      FIRST = first+1
!
! ----- display mu in a bar
!
      if (first .eq. 1) then
         iclr = 16
         i2 = setcolor(0)
         i4=rectangle_w($GFILLINTERIOR,0.d0,-0.1d0,6.d0,-0.2d0)
      end if  
      mu = (dlog10(mu)-dlog10(origin))*factorx
      if ( mu .ne. mu_old) then
        i2 = setcolor(iclr)
        i4=rectangle_w($GFILLINTERIOR,0.d0,-0.1d0,mu_old ,-0.2d0)
        i2=setcolor(11)
        i4=rectangle_w($GFILLINTERIOR,0.d0,-0.1d0,mu ,-0.2d0)
        mu_old = mu
	  if (iclr .eq. 13) then
	     iclr = 16
        else
           iclr = iclr - 1
        end if
      end if
      return
      end subroutine draw_rg_inf_path
