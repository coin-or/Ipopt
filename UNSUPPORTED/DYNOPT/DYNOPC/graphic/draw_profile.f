C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: draw_profile.f 531 2004-03-11 01:31:07Z andreasw $
      SUBROUTINE draw_profile(line_style,X)
 
      USE DFLIB
      implicit none

      include 'DAE2NLP.INC'
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DAENLP/, /DYNAUX/, /DYNOPC/, /GRAPH/

      INTEGER(2) dummy,i
      INTEGER(4) icolor(12)
      double precision TK,TMAX,FACTORTIME
      integer(4) line_style, style
      double precision X(*)    ! discretized  variables
      double precision timex,timepre, varix(12), varixpre(12)
      double precision scaleZ(12),z(nz),y(ny),u(nu)
      double precision factorz(12)
      double precision dt

      integer ix, jy 
      TYPE(windowconfig) myscreen  
      TYPE(wxycoord) wxy

      COMMON myscreen
      COMMON/MAX/TK, TMAX
      common/color/icolor

      double precision pre_x(0:NELEMAX*NCOLMAX)
      double precision pre_y(0:NELEMAX*NCOLMAX,2*NDSPLMAX) ! for previous profiles
      save pre_x, pre_y

      integer nv, nvg, j, k, nvy, iline
      double precision dmz(NZ)

      call INITD2N(ncol)  ! this subroutine is in dae2nlp

      if( OPT_FLAG ) call fade_lines(pre_x, pre_y)

      FACTORTIME = 25D0/TMAX

!/ only display the first 8 of z or specified by user

      nv = ndsplz
      do i = 1,ndsplz
        scalez(i) = zbg(2,indexz(i))-zbg(1,indexz(i))     
        if(zbg(2,indexz(i)) .ge. 1.d+20  .or. 
     1                    zbg(1,indexz(i)) .le.-1.d+20) then
            scalez(i) = 100.d0
            zbg(1,indexz(i)) = -50.d0
        end if
!/ get initial conditions of differential variables taking lower bound as 
!  the origin of Y-axis
        varixpre(i) = z(indexz(i))-zbg(1,indexz(i))
      end do

! This section is not needed any more for user specified number of displayed variables
C     if(nc .le. 10-nv) then 
C           nvu = nc
C      else 
C           nvu = nc+nv-10
C      end if
C
C     nvg = nv+nvu
C
C      if(ny .le. 12-nvu-nv) then
C           nvy = ny
C      else
C           nvy = 12 - nvu - nv
C      end if
!/ get vlaues of control and algebric variables at t=0
   
      timex = 0.d0
      call APPSLN(timex,nz,ny,nu,nele,ncol,ti,x,z,dmz,y,u) !dummy parameters 
                                             ! updated according to Andreas code 
      do i = 1, ndsplu
! Scale control for display 
        scalez(nv+i)=ubg(2,indexu(i))-ubg(1,indexu(i))
        if(ubg(2,indexu(i)) .ge. 1.d+20  
     1             .or. ubg(1,indexu(i)) .le.-1.d+20) then
          scalez(nv+i) = 100.d0
          ubg(1,i) = -50.d0
        end if
!/ graphics takes lower bounds origin of Y-axis 
        if (u(indexu(i)) .gt. ubg(2,indexu(i)) ) then 
          u(indexu(i)) = ubg(2,indexu(i))
        end if   
        varixpre(nv+i) = u(indexu(i))-ubg(1,indexu(i))  
! 
      end do
!
      nvg = ndsplz + ndsplu
      if(ndsply .ge. 1) then
        do i =1,ndsply
          scalez(nvg+i) = ybg(2,indexy(i)) - ybg(1,indexy(i))
          if(ybg(2,indexy(i)) .ge. 1.d+20  
     1             .or. ybg(1,indexy(i)) .le.-1.d+20)then
            scalez(nvg+i) = 100.d0
            ybg(1,indexy(i)) = 0
          end if        
          varixpre(nvg+i) = y(indexy(i))-ybg(1,indexy(i))
        end do
      end if        
!
! To get position in Y-axis of the Window coordinate
      nvg = nvg + ndsply
      do i=1,nvg
        factorz(i) = 100d0/scalez(i)
        varixpre(i) = varixpre(i)*factorz(i)
        pre_y(0,i) = varixpre(i)
      end do
      timepre = TI(1)*FACTORTIME
      pre_x(0) = timepre
      ix = 0
      jy = 0

      do i = 1,NELE
        dt = TI(i+1) - TI(i)
        do j=0,ncol
          if (j .eq. 0) then
            timex = TI(i) + 0.0001d0*dt
          else
            timex = TI(i) + rho(j)*dt
          end if
          call APPSLN(timex,nz,ny,nu,NELE,ncol,TI,x,z,dmz,y,u)
!          call      appsln(timex,z,y,u,p,x,xi,ipar,iflag,iout) 
          do k=1,ndsplz 
            varix(k) = (z(indexz(k))-zbg(1,indexz(k)))*factorz(k)
          end do
          do k = 1,ndsplu 
            if (u(indexu(k)) .gt. ubg(2,indexu(k)) ) then 
              u(indexu(k)) = ubg(2,indexu(k))
            end if   
            varix(nv+k) =
     1          (u(indexu(k))-ubg(1,indexu(k)))*factorz(nv+k)
          end do
          if(ndsply .ge. 1) then
            nvy = ndsplz + ndsplu
            do k =1,ndsply
              varix(nvy+k) = 
     1          (y(indexy(k))-ybg(1,indexy(k)))*factorz(nvy+k)
            end do
          end if
            
          timex = timex*factortime
          ix = ix +1
          pre_x(ix) = timex 
          if(line_style .eq. 0) style=#FF80
          if(line_style .eq. 1) style=#FFFF
          do iline = 1,nvg
            call moveto_w(timepre,varixpre(iline),wxy)
            dummy = setcolor(icolor(iline))
            call setlinestyle(style)
            dummy = lineto_w(timex,varix(iLINE))
            varixpre(iline) = varix(iline)
            pre_y(ix,iline) = varix(iline)
          END DO
          timepre = timex
        end do
      end do
      return
      END
