C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: read_cmp.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine read_cmp(MODNAM, STUBNAM, TOL_DASPK, TOL_DISCRET,
     1     TOL_OPT, ZB, Z0, YB, Y0, UB, U0, PB, P0, IERR, NLINES, LINES)
!DEC$ ATTRIBUTES DLLEXPORT :: READ_CMP
C--------------------------------------------------------------------------
C     Purpose: To read comprehensive data from existing file produced
C              previously by DYNOPC graphical interface
C
C         YDL                03-12-01
C         Andreas Waechter   09-03-01   cleaned up
C--------------------------------------------------------------------------

#ifdef WIN32
      use DFLIB
      use portlib
      use dialogm
#endif

      implicit none

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
#ifdef WIN32
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/
#endif

C     IN:    Name of subdirectory and filebasename for .cmp file
      character*(*) MODNAM, STUBNAM
C     OUT:   Tolerances for DASPK, discretization, and optimization
      double precision TOL_DASPK, TOL_DISCRET, TOL_OPT
C     OUT:   bounds and initial values for differential variables
      double precision ZB(2,*), Z0(*)
C     OUT:   bounds and initial values for algebraic variables
      double precision YB(2,*), Y0(*)
C     OUT:   bounds and initial values for control variables
      double precision UB(2,*), U0(*)
C     OUT:   bounds and initial values for parameters
      double precision PB(2,*), P0(*)
      integer IERR
      integer NLINES
      character*80 LINES(7)       ! for error output

      character*256 fname
      logical new, ex
      integer i, j, k, l_uvari, cnr
      double precision d
      character*1 s
      integer READ_ENTRY, CALC_NX

C
C     What for are the following?!?
C
      integer iac
      integer iEqual, iCnst, iLastU

      call GET_FILENAME(MODNAM, STUBNAM, 'cmp', fname)

      cnr = 1
      open(cnr,file=fname,status='old',err=55)
      new = .true.

C     skip first entry
      if( READ_ENTRY(cnr, new, 's', d, i, s).ne.0 ) goto 50

C     read problem size
      if( READ_ENTRY(cnr, new, 'i', d, NZ, s).ne.0 ) goto 50
      if( NZ.gt.NZMAX .or. NZ.lt.0 ) goto 40

      if( READ_ENTRY(cnr, new, 'i', d, NY, s).ne.0 ) goto 50
      if( NY.gt.NYMAX .or. NY.lt.0 ) goto 41

      if( READ_ENTRY(cnr, new, 'i', d, NU, s).ne.0 ) goto 50
      if( NU.gt.NUMAX .or. NU.lt.0 ) goto 42

      if( READ_ENTRY(cnr, new, 'i', d, NP, s).ne.0 ) goto 50
      if( NP.gt.NPMAX .or. NP.lt.0 ) goto 43

      if( READ_ENTRY(cnr, new, 'i', d, NELE, s).ne.0 ) goto 50
      if( NELE.gt.NELEMAX .or. NELE.le.0 ) goto 44

      if( READ_ENTRY(cnr, new, 'i', d, NCOL, s).ne.0 ) goto 50
      if( NCOL.gt.NCOLMAX .or. NCOL.le.0 ) goto 45

      i = CALC_NX(NZ, NY, NU, NP, NELE, NCOL)
      if( i.gt.NXMAX ) goto 145

C     In this context, number of additinal constraints is always zero
      if( READ_ENTRY(cnr, new, 'i', d, iac, s).ne.0 ) goto 50

C     Number of additional bound constraints
      if( READ_ENTRY(cnr, new, 'i', d, IADB, s).ne.0 ) goto 50

C     Read tolerances for DDASPK, discretization, optimization

      if( READ_ENTRY(cnr, new, 'd', TOL_DASPK, i, s).ne.0 ) goto 50
      if( READ_ENTRY(cnr, new, 'd', TOL_DISCRET, i, s).ne.0 ) goto 50
      if( READ_ENTRY(cnr, new, 'd', TOL_OPT, i, s).ne.0 ) goto 50
C     WANT THIS out for Not OLDSTYLE
      eps3 = TOL_OPT

C     Read iEqual, iCnst, iLastU, end_time  (NEEDED ????)
      if( READ_ENTRY(cnr, new, 'i', d, iEqual, s).ne.0 ) goto 50
      if( READ_ENTRY(cnr, new, 'i', d, iCnst, s).ne.0 ) goto 50
      if( READ_ENTRY(cnr, new, 'i', d, iLastU, s).ne.0 ) goto 50
      if( READ_ENTRY(cnr, new, 'd', end_time, i, s).ne.0 ) goto 50

C     Read discretization of time interval
      do j = 1, NELE+1
         if( READ_ENTRY(cnr, new, 'd', TI(j), i, s).ne.0 ) goto 50
      enddo

C     Read bounds and initial conditions for differential variables
      do j = 1, NZ
         if( READ_ENTRY(cnr, new, 'd', ZB(1,j), i, s).ne.0 ) goto 50
         if( READ_ENTRY(cnr, new, 'd', Z0(  j), i, s).ne.0 ) goto 50
         if( READ_ENTRY(cnr, new, 'd', ZB(2,j), i, s).ne.0 ) goto 50
      enddo
C     copy Z0 into ZINIT from DYNAUX.INC
      call DCOPY(NZ, Z0, 1, ZINIT, 1)

C     Read bounds and initial values for algebraic variables
      do j = 1, NY
         if( READ_ENTRY(cnr, new, 'd', YB(1,j), i, s).ne.0 ) goto 50
         if( READ_ENTRY(cnr, new, 'd', Y0(  j), i, s).ne.0 ) goto 50
         if( READ_ENTRY(cnr, new, 'd', YB(2,j), i, s).ne.0 ) goto 50
      enddo

C     Read bounds and initial values for control variables
      do j = 1, NU
         if( READ_ENTRY(cnr, new, 'd', UB(1,j), i, s).ne.0 ) goto 50
         if( READ_ENTRY(cnr, new, 'd', U0(  j), i, s).ne.0 ) goto 50
         if( READ_ENTRY(cnr, new, 'd', UB(2,j), i, s).ne.0 ) goto 50
      enddo

C     if the control variable has a given profile then read it and store
C     information in IU_PROF and U_PROF (see also DYNAUX.INC)

      if( READ_ENTRY(cnr, new, 'i', d, NU_PROF, s).ne.0 ) goto 50
      do j = 1, NU
         IU_PROF(j) = 0
      end do
      do j = 1, NU_PROF
         if( READ_ENTRY(cnr, new, 'i', d, i, s).ne.0 )
     1        goto 50
         if( i.gt.NU ) goto 46
         if( IU_PROF(i).ne.0 ) goto 47
         IU_PROF(i) = 1
         do k = 1, NELE+1
            if( READ_ENTRY(cnr, new, 'd', U_PROF(i,k), i, s)
     1           .ne.0 ) goto 50
         enddo
      enddo

C     Read bounds and initial point for parameters
      do j = 1, NP
         if( READ_ENTRY(cnr, new, 'd', PB(1,j), i, s).ne.0 ) goto 50
         if( READ_ENTRY(cnr, new, 'd', P0(  j), i, s).ne.0 ) goto 50
         if( READ_ENTRY(cnr, new, 'd', PB(2,j), i, s).ne.0 ) goto 50
      enddo

C     At this point, no information about fixed parameters in *.cmp file
      NP_FIX = 0
      do j = 1, NP
         IP_FIX(j) = 0
      end do

C     Read information for graphical display
      if( READ_ENTRY(cnr, new, 'i', d, ndsplz, s).ne.0 ) goto 50
      if( READ_ENTRY(cnr, new, 'i', d, ndsplu, s).ne.0 ) goto 50
      if( READ_ENTRY(cnr, new, 'i', d, ndsply, s).ne.0 ) goto 50
      do j = 1, ndsplz
         if( READ_ENTRY(cnr, new, 'i', d, indexZ(j), s).ne.0 ) goto 50
      enddo
      do j = 1, ndsplu
         if( READ_ENTRY(cnr, new, 'i', d, indexU(j), s).ne.0 ) goto 50
      enddo
      do j = 1, ndsply
         if( READ_ENTRY(cnr, new, 'i', d, indexY(j), s).ne.0 ) goto 50
      enddo

#ifdef WIN32
C     Copy bound information into /GRAPH/ common block
      call DCOPY(2*NZ, ZB(1,1), 1, ZBG(1,1), 1)
      call DCOPY(2*NY, YB(1,1), 1, YBG(1,1), 1)
      call DCOPY(2*NU, UB(1,1), 1, UBG(1,1), 1)
      call DCOPY(2*NP, PB(1,1), 1, PBG(1,1), 1)
#endif

      inquire(file = 'delete_state_bounds', exist = ex)
      if( ex ) then
         call C_OUT(2,0,1,
     1        '*** File ''delete_state_bounds'' detected'//
     1        '; deleting bounds on states. *** ')
         call DCOPY(NZ, -VLARGE, 0, ZB(1,1), 2)
         call DCOPY(NY, -VLARGE, 0, YB(1,1), 2)
         call DCOPY(NZ,  VLARGE, 0, ZB(2,1), 2)
         call DCOPY(NY,  VLARGE, 0, YB(2,1), 2)
      endif

C     That's it - got everything
      close(cnr,status='keep')
      return
C
C     Error messages
C
 40   write(LINES,90) 'NZ', NZ, 'NZMAX', NZMAX
      NLINES = 5
      IERR = -1
      return
 41   write(LINES,90) 'NY', NY, 'NYMAX', NYMAX
      NLINES = 5
      IERR = -1
      return
 42   write(LINES,90) 'NU', NU, 'NUMAX', NUMAX
      NLINES = 5
      IERR = -1
      return
 43   write(LINES,90) 'NP', NP, 'NPMAX', NPMAX
      NLINES = 5
      IERR = -1
      return
 44   write(LINES,90) 'NELE', NELE, 'NELEMAX', NELEMAX
      NLINES = 5
      IERR = -1
      return
 45   write(LINES,90) 'NCOL', NCOL, 'NCOLMAX', NCOLMAX
      NLINES = 5
      IERR = -1
      return
 145  write(LINES,90) 'NX', i, 'NXMAX', NXMAX
      NLINES = 5
      IERR = -1
      return
 46   write(LINES,100) 'IU_PROF', i, 'NU', NU
      NLINES = 5
      IERR = -1
      return
 47   write(LINES,110) 'U_PROF', i
      NLINES = 4
      IERR = -1
      return
C
 50   write(LINES,70)
      NLINES = 7
      IERR = -2
      return
C
 55   write(LINES,80)
      NLINES = 3
      IERR = -3
      return
C
 70   format(
     1     6x,'******************************************'/
     2     6x,'*     Reading data is not successful     *'/
     3     6x,'*         Please check the file          *'/
     4     6x,'*               "*.cmp"                  *'/
     5     6x,'*         in your directory and          *'/
     6     6x,'*     specified for current problem      *'/
     7     6x,'******************************************')  
 80   format(
     1     6x,'******************************************'/
     2     6x,'*         File *.cmp not found           *'/
     3     6x,'******************************************')
 90   format(
     1     6x,'**********************************************'/
     2     6x,'*         Maximal dimension exceeded:        *'/
     3     6x,'*     ',a4,'=',i7,', but ',a7,'=',i7,'.     *'/
     4     6x,'* You need to edit DYNAUX.INC and recompile! *'/
     7     6x,'**********************************************')  
 100  format(
     1     6x,'******************************************'/
     2     6x,'*              Invalid data:             *'/
     3     6x,'*    ',a7,'=',i7,', but ',a2,'=',i7,'.   *'/
     4     6x,'*          Check the *.cmp file          *'/
     7     6x,'******************************************')  
 110  format(
     1     6x,'**********************************************'/
     2     6x,'*                Invalid data:               *'/
     3     6x,'*  ',a6,'(',i7,') specified more than once. *'/
     4     6x,'*            Check the *.cmp file            *'/
     7     6x,'******************************************')  
      end
