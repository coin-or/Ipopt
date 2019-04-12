C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: read_adb.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine READ_ADB(MODNAM, STUBNAM, IERR, NLINES, LINES)
!DEC$ ATTRIBUTES DLLEXPORT :: READ_ADB
C--------------------------------------------------------------------------
C     Purpose: To read start point data from existing file produced
C              previously by simulation (.stp) or optimization (.sol)
C
C         Andreas Waechter, Yidong Lang   09-29-01
C--------------------------------------------------------------------------

      implicit none

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/, /DYNOPC/

C     IN:    Name of subdirectory and filebasename for .cmp file
      character*(*) MODNAM, STUBNAM
      integer IERR
      integer NLINES
      character*80 LINES(7)       ! for error output

      character*256 fname
      logical new
      integer i, j, cnr
      double precision d
      character*1 s
      integer READ_ENTRY, GET_IELE4T

C     There is nothing to do it IADB is zero
      if( IADB.eq.0 ) return

      call GET_FILENAME(MODNAM, STUBNAM, 'adb', fname)

      cnr = 1
      open(cnr,file=fname,status='old',err=55)

      new = .true.

      NZADB = 0
      NYADB = 0
      NUADB = 0

      if( READ_ENTRY(cnr, new, 'i', d, NZ_IN_ADB, s).ne.0 ) goto 50
      if( NZ_IN_ADB.lt.0 .or. NZ_IN_ADB.gt.NADBMAX ) goto 40

      if( READ_ENTRY(cnr, new, 'i', d, NE_IN_ADB, s).ne.0 ) goto 50
      if( NE_IN_ADB.lt.0 .or. NE_IN_ADB.gt.NADBMAX ) goto 40

      if( READ_ENTRY(cnr, new, 'i', d, NY_IN_ADB, s).ne.0 ) goto 50
      if( NY_IN_ADB.lt.0 .or. NY_IN_ADB.gt.NADBMAX ) goto 40

      if( READ_ENTRY(cnr, new, 'i', d, NU_IN_ADB, s).ne.0 ) goto 50
      if( NU_IN_ADB.lt.0 .or. NU_IN_ADB.gt.NADBMAX ) goto 40

C     Read time points
      do j = 1, NE_IN_ADB
         if( READ_ENTRY(cnr, new, 'd', d, i, s).ne.0 ) goto 50
         if( d.lt.TI(1) .or. d.gt.TI(NELE+1) ) goto 42
         TIME_ADB(j) = d
         IELE_ADB(j) = GET_IELE4T(d)
      enddo

C     Take care of Z bounds
      call READ_ADB2(NZ_IN_ADB, NE_IN_ADB, INDZADB, NZADB, IZADB,
     1     JZADB, ATTRZ, Z_ADB, NADBMAX, CNR, NEW, IERR)
      if( IERR.eq.-1 ) goto 50
      if( IERR.eq.-4 ) goto 41

C     Take care of Y bounds
      call READ_ADB2(NY_IN_ADB, NE_IN_ADB, INDYADB, NYADB, IYADB,
     1     JYADB, ATTRY, Y_ADB, NADBMAX, CNR, NEW, IERR)
      if( IERR.eq.-1 ) goto 50
      if( IERR.eq.-4 ) goto 41

C     Take care of U bounds
      call READ_ADB2(NU_IN_ADB, NE_IN_ADB, INDUADB, NUADB, IUADB,
     1     JUADB, ATTRU, U_ADB, NADBMAX, CNR, NEW, IERR)
      if( IERR.eq.-1 ) goto 50
      if( IERR.eq.-4 ) goto 41

C     That's it - got everything
      close(cnr,status='keep')

      IERR = 0
      return
C
C     Error messages
C
 40   write(LINES,90) NADBMAX
      NLINES = 6
      IERR = -1
      return
C
 41   write(LINES,85)
      NLINES = 5
      IERR = -4
      return
C
 42   write(LINES,86) d
      NLINES = 5
      IERR = -5
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
     2     6x,'*         File *.adb not found           *'/
     3     6x,'******************************************')
 85   format(
     1     6x,'**********************************************'/
     2     6x,'*       There seems to be bad data in        *'/
     3     6x,'*              the .adb file.                *'/
     4     6x,'*            Abort due to error.             *'/
     7     6x,'**********************************************')
 86   format(
     1     6x,'**********************************************'/
     2     6x,'*     There is a time point in .adb file     *'/
     3     6x,'*           outside time horizon!            *'/
     4     6x,'*            Abort due to error.             *'/
     7     6x,'**********************************************')
 90   format(
     1     6x,'**********************************************'/
     2     6x,'*         Maximal dimension exceeded:        *'/
     3     6x,'*       NADBMAX ='  ,i7, ' is too small.     *'/
     4     6x,'* You need to edit DYNOPC.INC and recompile! *'/
     5     6x,'*   Or there is bad data in the .adb file?   *'/
     7     6x,'**********************************************')
      end

C =============================================================================
      subroutine READ_ADB2(N_IN_ADB, NE_IN_ADB, INDADB, NADB, IADB,
     1     JADB, ATTR, D_ADB, NADBMAX, CNR, NEW, IERR)
C
C     "Internal" subroutine for reading information for each kind of varible
C
      implicit none
      integer N_IN_ADB, NE_IN_ADB, NADB
      integer NADBMAX, CNR, IERR
      integer INDADB(NADBMAX), IADB(NADBMAX), JADB(NADBMAX)
      integer ATTR(NADBMAX,NADBMAX)
      double precision D_ADB(NADBMAX,NADBMAX)
      logical NEW

      integer j, k, i
      double precision d
      character*1 s
      integer READ_ENTRY

      do j = 1, N_IN_ADB
         if( READ_ENTRY(CNR, NEW, 'i', d, INDADB(j), s).ne.0 ) goto 50
      enddo
      do j = 1, N_IN_ADB
         do k = 1, NE_IN_ADB
            if( READ_ENTRY(CNR, NEW, 'i', d, i, s).ne.0 ) goto 50
            if( i.ne.0 ) then
               NADB = NADB + 1
C
               if( READ_ENTRY(CNR, NEW, 'i', d, i, s).ne.0 ) goto 50
               if( i.lt.0 .or. i.ge.NADBMAX ) goto 41
               IADB(NADB) = i + 1
C
               if( READ_ENTRY(CNR, NEW, 'i', d, i, s).ne.0 ) goto 50
               if( i.lt.0 .or. i.ge.NADBMAX ) goto 41
               JADB(NADB) = i + 1
C
               if( READ_ENTRY(CNR, NEW, 'i', d, i, s).ne.0 ) goto 50
               ATTR(IADB(NADB),JADB(NADB)) = i
C
               if( READ_ENTRY(CNR, NEW, 'd', d, i, s).ne.0 ) goto 50
               D_ADB(IADB(NADB),JADB(NADB)) = d
            endif
 10      enddo
      enddo

      IERR = 0
      return

 41   IERR = -4
      return

 50   IERR = -2
      return
      end
