C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: write_stp.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine WRITE_STP(WARM, MODNAM, STUBNAM, X )
!DEC$ ATTRIBUTES DLLEXPORT :: WRITE_STP

      implicit none
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
#ifdef OLDSTYLE
#ifdef WIN32
      include 'DYNGRA.INC'
#else
      double precision zbg(2,NZMAX)
      double precision ybg(2,NYMAX)
      double precision ubg(2,NUMAX)
      double precision pbg(2,NPMAX)
#endif
#endif
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/

      logical WARM              ! if .false., use extention ".stp",
                                ! if .true. , use extention ".sol"
      character*(*) MODNAM, STUBNAM
      double precision X(*)

      character*80 fname
      character*40 datestr
      integer i, SLEN, lx, iele, icol
#ifdef OLDSTYLE
      integer CALC_NX
#endif
      if( WARM ) then
         call GET_FILENAME(MODNAM, STUBNAM, 'sol', fname)
      else
         call GET_FILENAME(MODNAM, STUBNAM, 'stp', fname)
      endif

      open(1,file = fname,status = 'unknown')

#ifdef OLDSTYLE
      write(1,*) NZ,NY,NU,NP,NELE,NCOL,0,IADB
      write(1,*) end_time
      write(1,*) eps3
      write(1,*) .false.
      lx = CALC_NX(NZ,NY,NU,NP,NELE,NCOL) + 1
      write(1,*)(x(i), i = 1,lx)
      write(1,*)(TI(i),i=1,NELE+1)
!
!------ store bounds --------------------
!
      write(1,*)(zbg(1,i),zbg(2,i), i = 1,NZ)
      write(1,*)(ubg(1,i),ubg(2,i), i = 1,NU)
      write(1,*)(ybg(1,i),ybg(2,i), i = 1,NY)
      write(1,*)(pbg(1,i),pbg(2,i), i = 1,NP)
!
!--------- store graphic display information ----
!
      write(1,*) ndsplz, ndsplu, ndsply
      write(1,*) (indexz(i), i = 1, ndsplz)
      write(1,*) (indexu(i), i = 1, ndsplu)
      write(1,*) (indexy(i), i = 1, ndsply)

#else
      call GET_TIME(datestr)
      write(1,100) MODNAM(1:SLEN(MODNAM)), STUBNAM(1:SLEN(STUBNAM)),
     1   datestr(1:SLEN(datestr))
 100  format('#',/,'#',3x,'Starting point file for',/,'#',9x,'Model  : '
     1       ,a,/,'#',9x,'Run    : ',a,/,'#',9x,'Created: ',a,/,'#')
      write(1,101) NZ,NY,NU,NP,NELE,NCOL
 101  format('#',/,'#',10x,'Information about problem size:',/,'#',
     1     /,i6,6x,'# NZ   (number of differential variables)',
     1     /,i6,6x,'# NY   (number of algebraic variables)',
     1     /,i6,6x,'# NU   (number of control variables)',
     1     /,i6,6x,'# NP   (number of parameters)',
     1     /,i6,6x,'# NELE (number of finite elements)',
     1     /,i6,6x,'# NCOL (number of collocation points per element)')
C      write(1,101) NZ,NY,NU,NP,NELE,NCOL,iac,IADB
C 101  format('#',/,'#',10x,'Information about problem size:',/,'#',
C     1       /,i6,6x,'# NZ   (number of differential variables)',
C     1       /,i6,6x,'# NY   (number of algebraic variables)',
C     1       /,i6,6x,'# NU   (number of control variables)',
C     1       /,i6,6x,'# NP   (number of parameters)',
C     1       /,i6,6x,'# NELE (number of finite elements)',
C     1       /,i6,6x,'# NCOL (number of collocation points per element)',
C     1       /,i6,6x,'# NAC  (always zero here)',
C     1       /,i6,6x,'# IADB (number of additional bound constraints)')
C      write(1,102) end_time
C 102  format(/,d23.16,6x,'# final time')
C      write(1,103) TOL_OPT
C 103  format(  d23.16,6x,'# tolerance for optimization')
C      write(1,104) duplicat
C 104  format(  l23,6x,'# starting point obtained by copying?')

      write(1,105)
 105  format(/,'#',/,'#   Starting point for optimization:',/,'#')

      lx = 0
      do i = 1, NZ
         write(1,150) X(lx+i),' Z',i,0
      enddo
      lx = lx + NZ
      do iele = 1, NELE
         do icol = 1, NCOL
            do i = 1, NZ
               write(1,155) X(lx+i),'dZ',i,iele,icol
            enddo
            lx = lx + NZ
            do i = 1, NY
               write(1,155) X(lx+i),' Y',i,iele,icol
            enddo
            lx = lx + NY
         enddo
         do i = 1, NZ
            write(1,150) X(lx+i),' Z',i,iele
         enddo
         lx = lx + NZ
      enddo
      do iele = 1, NELE
         do icol = 1, NCOL
            do i = 1, NU
               write(1,155) X(lx+i),' U',i,iele,icol
            enddo
            lx = lx + NU
         enddo
      enddo
      do i = 1, NP
         write(1,160) X(lx+i),' P',i
      enddo

 150  format(d23.15,'  #',a4,'(',i3,')  iele=',i4)
 155  format(d23.15,'  #',a4,'(',i3,')  iele=',i4,' icol=',i2)
 160  format(d23.15,'  #',a4,'(',i3,')')

      write(1,170)
 170  format(/,'#',/,'#   Discretization of time horizon:',/,'#')
      write(1,171)(TI(i),i,i=1,NELE+1)
 171  format(d23.15,'  #  TI(',i3,')')
C
C     I don't think we need this information in .stp or .sol files, since
C     it is already stored in .cmp
C
C!
C!------ store bounds --------------------
C!
C      write(1,180)
C 180  format(/,'#',/,'#   Bounds for variables:',/,'#')
C      write(1,181)(zb(1,i),zb(2,i),'Z',i, i = 1,NZ)
C      write(1,181)(ub(1,i),ub(2,i),'U',i, i = 1,NU)
C      write(1,181)(yb(1,i),yb(2,i),'Y',i, i = 1,NY)
C      write(1,181)(pb(1,i),pb(2,i),'P',i, i = 1,NP)
C 181  format(2d23.15,'  #',a4,'(',i3,')')
C!
C!--------- store graphic display information ----
C!
C      write(1,190)
C 190  format(/,'#',/,'#   Some information for graphic display:',/,'#')
C
C      write(1,191) ndsplz, ndsplu, ndsply
C      write(1,191) (indexz(i), i = 1, ndsplz)
C      write(1,191) (indexu(i), i = 1, ndsplu)
C      write(1,191) (indexy(i), i = 1, ndsply)
C 191  format(3i10)
#endif

      close(1,status='keep')
      return
      end
