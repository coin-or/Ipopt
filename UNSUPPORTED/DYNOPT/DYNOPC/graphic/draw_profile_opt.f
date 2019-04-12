C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      SUBROUTINE draw_profile(line_style,iopt)
 
      USE QWgraph
      USE Dyncom


      implicit none

      include 'colcom.inc'
      include 'ip2dae.inc'

      INTEGER(2) dummy,i
      INTEGER(4) icolor(12)
      double precision xx,TK,TMAX,FACTORTIME
      integer    ipar(1)
      integer(2) line_style, style
      double precision timex,timepre, varix(12), varixpre(12)
      double precision scaleZ(12),z(nz),y(ny),u(nu), p(np)
      double precision factorz(12)
      double precision dt

      integer ix, jy 
      TYPE(windowconfig) myscreen  
      TYPE(wxycoord) wxy
      character*12 fname
 

      COMMON myscreen
      COMMON/MAX/TK, TMAX
      common/color/icolor

      integer nv, iopt, nvg, j, iz, iu, iy, k, nvy, iline
      double precision dmz(NZ)

      call INITCOM(ncol)  ! this subroutine is in dae2nlp
      call fade_lines()

      FACTORTIME = 25D0/TMAX


!/ only display the first 8 of z or specified by user

      nv = ndsplz
      do i = 1,ndsplz
        scalez(i) = zb(2,indexz(i))-zb(1,indexz(i))     
        if(zb(2,indexz(i)) .ge. 1.d+20  .or. 
     1                    zb(1,indexz(i)) .le.-1.d+20) then
            scalez(i) = 100.d0
             zb(1,indexz(i)) = -50.d0
        end if
!/ get initial conditions of differential variables taking lower bound as 
!  the origin of Y-axis
        varixpre(i) = x(indexz(i))-zb(1,indexz(i))
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
      call appsln(timex,nz,ny,nu,ne,ncol,Ti,x,z,dmz,y,u) !dummy parameters 
                                             ! updated according to Andreas code 
                                             ! in dae2nlp
      do i = 1, ndsplu
! Scale control for display 
        scalez(nv+i)=ub(2,indexu(i))-ub(1,indexu(i))
        if(ub(2,indexu(i)) .ge. 1.d+20  
     1             .or. ub(1,indexu(i)) .le.-1.d+20) then
             scalez(nv+i) = 100.d0
             ub(1,i) = -50.d0
        end if
!/ graphics takes lower bounds origin of Y-axis 
         if (u(indexu(i)) .gt. ub(2,indexu(i)) ) then 
            u(indexu(i)) = ub(2,indexu(i))
         end if   
         varixpre(nv+i) = u(indexu(i))-ub(1,indexu(i))  
! 
      end do
!
      nvg = ndsplz + ndsplu
      if(ndsply .ge. 1) then
        do i =1,ndsply
           scalez(nvg+i) = yb(2,indexy(i)) - yb(1,indexy(i))
        if(yb(2,indexy(i)) .ge. 1.d+20  
     1              .or. yb(1,indexy(i)) .le.-1.d+20)then
             scalez(nvg+i) = 100.d0
             yb(1,indexy(i)) = 0
        end if        
            varixpre(nvg+i) = y(indexy(i))-yb(1,indexy(i))
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
      timepre = xi(1)*FACTORTIME
      pre_x(0) = timepre
      ix = 0
      jy = 0

      do i = 1,ne
         dt = xi(i+1) - xi(i)
         do j=0,ncol
             if (j .eq. 0) then
                timex = xi(i) + 0.0001d0*dt
             else
                  timex = xi(i) + rho(j)*dt
             end if
             call appsln(timex,nz,ny,nu,ne,ncol,Ti,x,z,dmz,y,u)
                do k=1,ndsplz 
                      varix(k) = (z(indexz(k))-zb(1,indexz(k)))*factorz(k)
              end do
              do k = 1,ndsplu 
                  if (u(indexu(k)) .gt. ub(2,indexu(k)) ) then 
                    u(indexu(k)) = ub(2,indexu(k))
                        end if   
                  varix(nv+k) =
     1                        (u(indexu(k))-ub(1,indexu(k)))*factorz(nv+k)
              end do
                  if(ndsply .ge. 1) then
                  nvy = ndsplz + ndsplu
                  do k =1,ndsply
                         varix(nvy+k) = 
     1                 (y(indexy(k))-yb(1,indexy(k)))*factorz(nvy+k)
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
