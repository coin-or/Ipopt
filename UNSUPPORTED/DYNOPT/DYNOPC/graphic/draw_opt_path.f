C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: draw_opt_path.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine draw_opt_path(alfav,alfatau,alfax, kkt_er, obj0, iter)

      USE dflib
      
      implicit none
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/
      integer(4) iter, iter_n, iter_0, i4
      integer(2) i, i2
      double precision obj, kkt_er
      double precision x_iter, x_iter_0
      double precision delt,wx1,wy1,wx2,wy2
      double precision objmin, objmax, objscale, objfactor, obj_origin 
      double precision alfaV, alfaTau, alfaX 
      double precision obj0, objmult
      logical over
      character*25 str
      type (wxycoord) wxy
      double precision obj_f, obj_f_0

 
      i4=setactiveqq(opt_path_W)
      iter_n = iter 
      if(iter_n .eq. 0) then
!
! --------------- Scaling objective function ----------------
!
      !--- Find the magtitute of the obj and determine its scale------
          obj = dabs(obj0)

        if (obj0 .lt. 0) then
          objmin = dfloat(int4(obj*0.9d0))
          objmax = dfloat(int4(obj*4.0d0))
        else
          objmin = dfloat(int4(obj*0.1d0))
          objmax = dfloat(int4(obj*1.8d0))
        end if
          objmult = 1.d0
 100      if (int4(obj) .gt. 10) then
             obj= obj/10d0
             objmult=objmult*10.d0
             goto 100
          else 
             if(int4(obj) .eq. 0) then
               obj = obj*10.d0
               objmult=objmult/10.d0
             goto 100
             end if
          end if
          objscale = int4((objmax - objmin)/objmult)
          if(objscale .eq. 0) objscale = 1.d0
          objfactor = 100.d0/objscale
          obj_origin = dfloat(int4((objmin/objmult)*10)/10)
          CALL MOVETO_W(-2.7D0,100.D0,wxy)
            i2=setcolor(11)
 !         i2=setcolorRGB(#0011FF)
          i2=SETFONT('t''Times new Roman''h12w8peb')
          CALL outGtext('OBJ')
          i4=SETFONT('t''Times new Roman''h12w6pb')
          do i=0,4
             call moveto_w(-2.d0,2.d0+2*i*10.d0,wxy)
             write(str,'(f5.2)') (objscale/5)*(i) 
     1                           + obj_origin
             call outGtext(str)
          end do
          call moveto_w(-3.2d0, -15.d0, wxy)
          write(str,'(F17.6)')(objmult)
          call outGtext(' Multiplier =')
          call moveto_w(1.d0, -15.d0, wxy)
          call outGtext(str) 
      end if
      if(iter_n .gt. (i_25+1)*25) then
        CALL CLEARSCREEN($GCLEARSCREEN)
        over = .true.
        i_25 = i_25+1
          call open_alfa_w()
            CALL MOVETO_W(-2.7D0,100.D0,wxy)
            i2=setcolor(11)
            i2=SETFONT('t''Times new Roman''h12w8peb')
         CALL outGtext('OBJ')
         i4=SETFONT('t''Times new Roman''h12w6pb')
         do i=0,4
            call moveto_w(-2.d0,2.d0+2*i*10.d0,wxy)
            write(str,'(f5.2)') (objscale/5)*(i) 
     1                           + obj_origin     
            call outGtext(str)
         end do
         call moveto_w(-3.2d0, -15.d0, wxy)
         write(str,'(F17.6)')(objmult)
         call outGtext(' Multiplier =')
         call moveto_w(1.d0, -15.d0, wxy)
         call outGtext(str)  
      end if
c     i2=setcolor(6)
c     I2=ELLIPSE_w($GFILLINTERIOR,wx1,wy1,wx2,wy2)

!     iter      = iter_n
!/ display absolute value of objective function
C     obj_f     = dabs(obj0)
!/ do not exceed the boundary of graphic frame
C     if(obj_f .gt. objmax) obj_f = objmax
C     if(obj_f .lt. objmin) obj_f = objmin
C      obj_f = (obj_f - objmin)/objmult
C      obj_f = obj_f*objfactor
!     kkt_max = opt_eps*1.d+9
!     if(kkt_er .gt. kkt_max) kkt_er = kkt_max
!     max_viol  = kkt_er
      if(over)iter_n=iter - (i_25)*25
!     max_viol=(dlog10(max_viol)+offset)*10.d0
      x_iter=dfloat(iter_n)
      if (iter_n .eq. 0) then 
         iter_0=iter_n
      else
         iter_0 = iter_n - 1
      end if
      x_iter_0=dfloat(iter_0)
 
      if(iter_n .eq. 1) THEN
!         CALL MOVETO_W(x_iter,obj_f,wxy)
!           delt = 0.4d0
!           wx1=x_iter-delt
!           wy1=obj_f+delt*2.d0
!           wx2=x_iter+delt
!           wy2=obj_f-delt*2.d0
!           i2=setcolor(11)
!           I2=ELLIPSE_w($GFILLINTERIOR,wx1,wy1,wx2,wy2)
!         CALL MOVETO_W(x_iter,max_viol,wxy)
!           delt = .4d0
!           wx1=x_iter-delt*2.d0
!           wy1=max_viol+delt*2.d0
!         wx2=x_iter+delt*2.d0
!           wy2=max_viol-delt*2.d0
!         i2=setcolorRGB(#0000ff)
!         I2=ELLIPSE_w($GFILLINTERIOR,wx1,wy1,wx2,wy2)
!     
!          I2=setcolorrgb(#00ff00)
!          i4=rectangle_w($GFILLINTERIOR,
!     1                 iter_n,0.d0,iter_n+.4d0,100.d0*alfaV)
!          I2=setcolorrgb(#00a0ff)
!          i4=rectangle_w($GFILLINTERIOR,
!     1                 iter_n+.35d0,0.d0,iter_n+.8d0,100.d0*alfaTau)
!         I2=setcolorrgb(#00fff0)
!          i4=rectangle_w($GFILLINTERIOR,
!     1                 iter_n+.6d0,0.d0,iter_n+.9d0,100.d0*alfaX)
!          goto 20
      end if

! ------------------------------for KKT--------------
!     delt = 0.4d0
!     wx1=x_iter-delt
!     wy1=obj_f+delt*2.d0
!     wx2=x_iter+delt
!     wy2=obj_f-delt*2.d0
!     i2=setcolor(11)
!     I2=ELLIPSE_w($GFILLINTERIOR,wx1,wy1,wx2,wy2)
  
!      call moveto_w(x_iter_0,max_viol_0,wxy)

!     i4=setcolorRGB(#0000ff)
!     i4=lineto_w(x_iter,max_viol)
!     delt = 0.8d0
!     wx1=x_iter-delt
!     wy1=max_viol+delt*2.d0
!     wx2=x_iter+delt
!     wy2=max_viol-delt*2.d0
!     i2=setcolorRGB(#0000ff)
!     I2=ELLIPSE_w($GFILLINTERIOR,wx1,wy1,wx2,wy2)
!
!-------------For Alfa's ---------------------------------------------------
!
          I2=setcolorrgb(#00ff00)
          i4=rectangle_w($GFILLINTERIOR,
     1                  x_iter,0.d0,x_iter+.4d0,100.d0*alfaV)
          I2=setcolorrgb(#00a0ff)
          i4=rectangle_w($GFILLINTERIOR,
     1                  x_iter+.35d0,0.d0,x_iter+.6d0,100.d0*alfaTau)
          I2=setcolorrgb(#00fff0)
          i4=rectangle_w($GFILLINTERIOR,
     1                  x_iter+.6d0,0.d0,x_iter+.8d0,100.d0*alfaX)
!
! ---------------------- For Obj --------------------------------------------
!
!/ display absolute value of objective function
      obj_f     = dabs(obj0)
!/ do not exceed the boundary of graphic frame
      if(obj_f .gt. objmax) obj_f = objmax
      if(obj_f .lt. objmin) obj_f = objmin
      obj_f = obj_f/objmult - obj_origin
      obj_f = obj_f*objfactor
      call moveto_w(x_iter,obj_f_0,wxy)
 
! Draw a circle at the current point
      delt = 0.3d0
      wx1=x_iter+1-delt
      wy1=obj_f+delt*8.d0
      wx2=x_iter+1+delt
      wy2=obj_f-delt*8.d0
      i4=setcolorRGB(#1111FF)
      I2=ELLIPSE_w($GFILLINTERIOR,wx1,wy1,wx2,wy2)
!     i2=setcolor(11)
      i4=lineto_w(x_iter+1,obj_f)

20    obj_f_0=obj_f
!     max_viol_0 = max_viol

      return
      end subroutine draw_opt_path
