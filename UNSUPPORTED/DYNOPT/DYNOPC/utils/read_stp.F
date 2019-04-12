C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: read_stp.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine READ_STP(WARM, MODNAM, STUBNAM, X_STP, TI_STP,
     1     NELE_STP, NCOL_STP, IERR, NLINES, LINES)
!DEC$ ATTRIBUTES DLLEXPORT :: READ_STP
C--------------------------------------------------------------------------
C     Purpose: To read start point data from existing file produced
C              previously by simulation (.stp) or optimization (.sol)
C
C         Andreas Waechter   09-29-01
C--------------------------------------------------------------------------

      implicit none

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/, /DYNOPC/
#ifdef OLDSTYLE
#ifdef WIN32
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLEXPORT :: /GRAPH/
#else
      double precision zbg(2,NZMAX)
      double precision ybg(2,NYMAX)
      double precision ubg(2,NUMAX)
      double precision pbg(2,NPMAX)
#endif
#endif

C     IN:    If .true., use extention '.sol', otherwise '.stp'
      logical WARM
C     IN:    Name of subdirectory and filebasename for .cmp file
      character*(*) MODNAM, STUBNAM
C     OUT:   starting point (i.e. values of discretized variables)
      double precision X_STP(*)
C     OUT:   discretization used to obtain X_STP
      double precision TI_STP(*)
      integer NELE_STP          ! no. of elements in STP/SOL files
      integer NCOL_STP          ! no. of collocation points in STP/SOL files
      integer IERR
      integer NLINES
      character*80 LINES(7)       ! for error output

      character*256 fname
      logical new
      integer i, j, cnr, nx
      double precision d
      character*1 s
      integer READ_ENTRY, CALC_NX

      if( WARM ) then
         call GET_FILENAME(MODNAM, STUBNAM, 'sol', fname)
      else
         call GET_FILENAME(MODNAM, STUBNAM, 'stp', fname)
      endif

      cnr = 1
      open(cnr,file=fname,status='old',err=55)

#ifdef OLDSTYLE

      read(1,*) nz,ny,nu,np,nele_stp,ncol_stp,i, IADB
      nx = CALC_NX(nz,ny,nu,np,nele_stp,ncol_stp) + 1
      read(1,*) end_time
      read(1,*) eps3
      read(1,*) new

      read(1,*)(x_stp(i), i = 1,nx)
      read(1,*)(Ti_stp(i),i=1,nele+1)

      read(1,*)((zbg(i,j), i = 1,2), j = 1,nz)
      read(1,*)((ubg(i,j), i = 1,2), j = 1,nu)
      read(1,*)((ybg(i,j), i = 1,2), j = 1,ny)
      read(1,*)((pbg(i,j), i = 1,2), j = 1,np)

      read(1,*) ndsplz, ndsplu, ndsply
      read(1,*) (indexz(i), i = 1, ndsplz)
      read(1,*) (indexu(i), i = 1, ndsplu)
      read(1,*) (indexy(i), i = 1, ndsply)

#else

      new = .true.

C     read problem size (check if consistent with .cmp data)
      if( READ_ENTRY(cnr, new, 'i', d, i, s).ne.0 ) goto 50
      if( i.ne.NZ ) goto 40

      if( READ_ENTRY(cnr, new, 'i', d, i, s).ne.0 ) goto 50
      if( i.ne.NY ) goto 41

      if( READ_ENTRY(cnr, new, 'i', d, i, s).ne.0 ) goto 50
      if( i.ne.NU ) goto 42

      if( READ_ENTRY(cnr, new, 'i', d, i, s).ne.0 ) goto 50
      if( i.ne.NP ) goto 43

      if( READ_ENTRY(cnr, new, 'i', d, NELE_STP, s).ne.0 ) goto 50
      if( NELE_STP.gt.NELEMAX .or. NELE_STP.le.0 ) goto 44

      if( READ_ENTRY(cnr, new, 'i', d, NCOL_STP, s).ne.0 ) goto 50
      if( NCOL_STP.gt.NCOLMAX .or. NELE_STP.le.0 ) goto 45

C     Determine size of vector X
      nx = CALC_NX(NZ, NY, NU, NP, NELE_STP, NCOL_STP)
      if( nx.gt.NXMAX ) goto 46

C     Read starting point
      do j = 1, nx
         if( READ_ENTRY(cnr, new, 'd', X_STP(j), i, s).ne.0 ) goto 50
      enddo

C     Read time discretization
      do j = 1, NELE_STP+1
         if( READ_ENTRY(cnr, new, 'd', TI_STP(j), i, s).ne.0 ) goto 50
      enddo

#endif

C     That's it - got everything
      close(cnr,status='keep')
      return
C
C     Error messages
C
 40   write(LINES,85) 'NZ', NZ, 'NZ', i
      NLINES = 6
      IERR = -4
      return
 41   write(LINES,85) 'NY', NY, 'NY', i
      NLINES = 6
      IERR = -4
      return
 42   write(LINES,85) 'NU', NU, 'NU', i
      NLINES = 6
      IERR = -4
      return
 43   write(LINES,85) 'NP', NP, 'NP', i
      NLINES = 6
      IERR = -4
      return
 44   write(LINES,90) 'NELE', NELE, 'NELEMAX', NELEMAX
      NLINES = 5
      IERR = -1
      return
 45   write(LINES,90) 'NCOL', NCOL, 'NCOLMAX', NCOLMAX
      NLINES = 5
      IERR = -1
      return
 46   write(LINES,90) 'NX', nx, 'NXMAX', NXMAX
      NLINES = 5
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
     4     6x,'*          "*.stp" or "*.sol"            *'/
     5     6x,'*         in your directory and          *'/
     6     6x,'*     specified for current problem      *'/
     7     6x,'******************************************')  
 80   format(
     1     6x,'******************************************'/
     2     6x,'*     File *.stp or *.sol not found      *'/
     3     6x,'******************************************')
 85   format(
     1     6x,'**********************************************'/
     2     6x,'*        Inconsisent model dimensions:       *'/
     3     6x,'*      ',a4,' in the CMP file is ', i7 , ',      *'/
     3     6x,'*  but ',a4,' in the STP/SOL file is ', i7 , '.  *'/
     4     6x,'* You need to use consistent starting point. *'/
     7     6x,'**********************************************')  
 90   format(
     1     6x,'**********************************************'/
     2     6x,'*         Maximal dimension exceeded:        *'/
     3     6x,'*          ',a4,'=',i7,', but ',a7,'=',i7,'.        *'/
     4     6x,'* You need to edit DYNAUX.INC and recompile! *'/
     7     6x,'**********************************************')  
      end
