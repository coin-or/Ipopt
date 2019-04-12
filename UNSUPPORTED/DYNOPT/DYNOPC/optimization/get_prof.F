C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C =============================================================================
C $Id: get_prof.F 531 2004-03-11 01:31:07Z andreasw $
C =============================================================================
      subroutine FMAIN( MODNAM, STUBNAM )
C =============================================================================
C
C     Get profile of variables from *.stp or *.sol file
C
C     Authors:  Andreas Waechter    01-15-02
C
C =============================================================================

      implicit none

      character*(*) MODNAM
      character*(*) STUBNAM

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
      include 'DAE2NLP.INC'

C     Bounds and initial values for model variables read from .cmp file
      double precision zb(2,NZMAX), z0(NZMAX)
      double precision yb(2,NYMAX), y0(NYMAX)
      double precision ub(2,NUMAX), u0(NUMAX)
      double precision pb(2,NPMAX), p0(NPMAX)

C     Derivatives for differential variables
      double precision dz(NZMAX)

C     Dummy residual for initialization of MODEL (MODEL_INIT)
      double precision f(NZMAX+NYMAX)

C     Time discretization from starting point file
      double precision ti_stp(NELEMAX+1)
      integer nele_stp, ncol_stp

C     Tolerances read from .cmp file
      double precision tol_daspk, tol_discret, tol_opt

C     Vector X for discretization
      double precision x(NXMAX)

      double precision ZVAL(NZMAX)
      double precision DZVAL(NZMAX)
      double precision YVAL(NYMAX)
      double precision UVAL(NUMAX)
      character*1 inchar

      integer ierr, nlines, ndegu, i, lvar, j, lx, lp
      double precision t
      character*80 lines(7)
      character*256 fname
      logical warm
      character*1 answer

C =============================================================================
C
C     Read Data Files
C
C =============================================================================

      write(*,101)
 101  format(/
     1     '>>>>> Starting reading input data.')

C     Read the *.cmp file
      call READ_CMP(MODNAM, STUBNAM, tol_daspk, tol_discret,
     1     tol_opt, zb, z0, yb, y0, ub, u0, pb, p0, ierr, nlines, lines)
      if( ierr.ne.0 ) then
         OutputUnit_COUT = 6
         call C_OUT(0,0,nlines,lines)
         stop
      endif

C     Check, whether there is a *.sol file.  If so, ask if user wants
C     to use this as starting point ('warm start')
      call GET_FILENAME(MODNAM, STUBNAM, 'sol', fname)
      inquire(file=fname, exist=warm)
      if( warm ) then
#ifdef WIN32
         write(*,102)
 102     format(/,'  Do you want to read the *.sol file? (y/n): '\)
#else
         call CPRINTF('\n%s',
     1        '  Do you want to read the *.sol file? (y/n): ')
#endif
         read(*,*) answer
         warm = (answer.eq.'Y' .or. answer.eq.'y')
      endif

C     Read the starting point
      call READ_STP(warm, MODNAM, STUBNAM, x, ti_stp,
     1     nele_stp, ncol_stp, ierr, nlines, lines)
      if( ierr.ne.0 ) then
         write(*,*)
         OutputUnit_COUT = 6
         call C_OUT(0,0,nlines,lines)
         stop
      endif

C     If necessary, fit starting point to a new TI grid
      call UPDATE_STP(x, ti_stp, nele_stp, ncol_stp, ierr)

      open(10,file='PROF.OUT',status='unknown')
C
C     Input
C
 100  write(*,*)
      write(*,*) 'Input choice (for help type ''h'', quit with ''q''):'
      read(*,'(a)') inchar

      if( inchar.eq.'z' .or. inchar.eq.'Z' ) then
C
C     Differential variable
C
         write(*,*) 'Input variable number'
         read(*,*) lvar
         if( lvar.lt.1 .or. lvar.gt.NZ ) then
            write(*,*) 'Invalid choice.  There are only ',NZ,
     1           ' differential variables.'
            goto 100
         endif
         write(10,8100) inchar, lvar
         lx = lvar
         do i = 1, NELE
            T = TI(i)
            write(10,8030) T, X(lx)
            do j = 1, NCOL-1
               T = TI(i) + RHO(j)*(TI(i+1)-TI(i))
               call APPSLN(T, NZ, NY, NU, NELE, NCOL, TI, X,
     1              ZVAL, DZVAL, YVAL, UVAL)
               write(10,8030) T, ZVAL(lvar)
            enddo
            lx = lx + NCOL*(NZ+NY) + NZ
         enddo
         T = TI(NELE+1)
         write(10,8030) T, X(lx)

      elseif( inchar.eq.'y' .or. inchar.eq.'Y' ) then
C
C     Differential variable
C
         write(*,*) 'Input variable number'
         read(*,*) lvar
         if( lvar.lt.1 .or. lvar.gt.NY ) then
            write(*,*) 'Invalid choice.  There are only ',NY,
     1           ' algebraic variables.'
            goto 100
         endif
         write(10,8100) inchar, lvar
         lx = 2*NZ + lvar
         do i = 1, NELE
            do j = 1, NCOL
               T = TI(i) + RHO(j)*(TI(i+1)-TI(i))
               write(10,8030) T, X(lx)
               lx = lx + NZ + NY
            enddo
            lx = lx + NZ
         enddo

      elseif( inchar.eq.'u' .or. inchar.eq.'U' ) then
C
C     Differential variable
C
         write(*,*) 'Input variable number'
         read(*,*) lvar
         if( lvar.lt.1 .or. lvar.gt.NU ) then
            write(*,*) 'Invalid choice.  There are only ',NU,
     1           ' control variables.'
            goto 100
         endif
         write(10,8100) inchar, lvar
         lx = NZ + NELE*(NCOL*(NZ+NY) + NZ) + lvar
         do i = 1, NELE
            do j = 1, NCOL
               T = TI(i) + RHO(j)*(TI(i+1)-TI(i))
               write(10,8030) T, X(lx)
               lx = lx + NU
            enddo
         enddo

      elseif( inchar.eq.'p' .or. inchar.eq.'P' ) then
C
C     Parameter
C
         write(*,*) 'Input variable number'
         read(*,*) lvar
         if( lvar.lt.1 .or. lvar.gt.NP ) then
            write(*,*) 'Invalid choice.  There are only ',NP,
     1           ' parameters.'
            goto 100
         endif
         write(10,8100) inchar, lvar
         lp = NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + lvar
         write(10,8030) X(lp)

      elseif( inchar.ne.'q' .and. inchar.ne.'Q' ) then
C
C     Print user manual
C
         write(*,*)
         write(*,*) 'Usage:'
         write(*,*)
         write(*,*)
         write(*,*) 
     1        '   Write values of variable for file ''PROF.OUT'':'
         write(*,*)
         write(*,*) '   - Input variable type:'
         write(*,*) '      ''z'': differential variable'
         write(*,*) '      ''y'': algebraic variable'
         write(*,*) '      ''u'': control variable'
         write(*,*) '      ''p'': parameter'
         write(*,*) '   - Input variable number'
         write(*,*)

      endif
C
      if( inchar.ne.'q' .and. inchar.ne.'Q' ) goto 100
C
C     That should be it...
C
      close(10)

 8000 format(a)
 8010 format(i10)
 8020 format(3d23.15)
 8030 format(2e23.15)
 8100 format(/,' --------  Variable ',a2,'(',i5,') --------',/)
      

      return
      end
