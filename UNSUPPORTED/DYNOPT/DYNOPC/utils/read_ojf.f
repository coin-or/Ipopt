C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: read_ojf.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine READ_OJF(MODNAM, STUBNAM, NDEGU, IERR, NLINES, LINES)
!DEC$ ATTRIBUTES DLLEXPORT :: READ_OJF
C--------------------------------------------------------------------------
C     Purpose: To read information about the (linear) objective function
C              from a file .ojf
C
C     Authors: Yidong Lang, Andreas Waechter    09-29-01
C--------------------------------------------------------------------------

      implicit none

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/, /DYNOPC/

C     IN:    Name of subdirectory and filebasename for .cmp file
      character*(*) MODNAM, STUBNAM
      integer NDEGU             ! Degree for control variable polynomial
      integer IERR
      integer NLINES
      character*80 LINES(7)     ! for error output

      character*256 fname
      logical new
      integer i, j, k, cnr
      double precision d
      character*1 s
      integer READ_ENTRY

C
C     What for are the following?!?  (For Visual Basic interface?)
C
      integer obj_mode

      call GET_FILENAME(MODNAM, STUBNAM, 'ojf', fname)

      cnr = 1
      open(cnr,file=fname,status='old',err=55)
      new = .true.

C     read OBJ_MODE   -  IGNORE FOR NOW
      if( READ_ENTRY(cnr, new, 'i', d, obj_mode, s).ne.0 ) goto 50

C     read number of factors
      if( READ_ENTRY(cnr, new, 'i', d, NZ_IN_OBJ, s).ne.0 ) goto 50
      if( NZ_IN_OBJ.gt.NOBJMAX .or. NZ_IN_OBJ.lt.0 ) goto 40
      if( READ_ENTRY(cnr, new, 'i', d, NY_IN_OBJ, s).ne.0 ) goto 50
      if( NY_IN_OBJ.gt.NOBJMAX .or. NY_IN_OBJ.lt.0 ) goto 41
      if( READ_ENTRY(cnr, new, 'i', d, NE_IN_OBJ, s).ne.0 ) goto 50
      if( NE_IN_OBJ.gt.NOBJMAX .or. NE_IN_OBJ.lt.0 ) goto 42
      if( READ_ENTRY(cnr, new, 'i', d, NP_IN_OBJ, s).ne.0 ) goto 50
      if( NP_IN_OBJ.gt.NOBJMAX .or. NP_IN_OBJ.lt.0 ) goto 43

C     Z_INDEX
      do j = 1, NZ_IN_OBJ
         if( READ_ENTRY(cnr, new, 'i', d, Z_INDEX(j), s).ne.0 ) goto 50
         if( Z_INDEX(j).le.0 .or. Z_INDEX(j).gt.NZ ) goto 44
      enddo

C     Y_INDEX
      do j = 1, NY_IN_OBJ
         if( READ_ENTRY(cnr, new, 'i', d, Y_INDEX(j), s).ne.0 ) goto 50
         if( Y_INDEX(j).le.0 .or. Y_INDEX(j).gt.NY ) goto 45
      enddo

C     E_INDEX
      do j = 1, NE_IN_OBJ
         if( READ_ENTRY(cnr, new, 'i', d, E_INDEX(j), s).ne.0 ) goto 50
         if( E_INDEX(j).le.-1 .or. E_INDEX(j).gt.NELE+1 ) goto 46
         if( E_INDEX(j).eq.0 ) E_INDEX(j) = NELE+1
      enddo

C     P_INDEX
      do j = 1, NP_IN_OBJ
         if( READ_ENTRY(cnr, new, 'i', d, P_INDEX(j), s).ne.0 ) goto 50
         if( P_INDEX(j).le.0 .or. P_INDEX(j).gt.NP ) goto 47
      enddo

C     AIJZ
      do j = 1, NE_IN_OBJ
         do k = 1, NZ_IN_OBJ
            if( READ_ENTRY(cnr, new, 'd', AIJZ(k,j), i, s).ne.0 )
     1           goto 50
         enddo
      enddo

C     BIJY
      do j = 1, NE_IN_OBJ
         do k = 1, NY_IN_OBJ
            if( READ_ENTRY(cnr, new, 'd', BIJY(k,j), i, s).ne.0 )
     1           goto 50
         enddo
      enddo

C     CIP
      do k = 1, NP_IN_OBJ
         if( READ_ENTRY(cnr, new, 'd', CIP(k), i, s).ne.0 ) goto 50
      enddo

C     TODO: Put this rather into CMP file?????

      if( READ_ENTRY(cnr, new, 'i', d, NDEGU, s).ne.0 ) goto 50
      if( NDEGU.lt.0 .or. NDEGU.gt.NCOL-1 ) goto 48
CDELETEME
      WRITE(*,*) 'NDEGU = ',ndegu

C     That's it - got everything
      close(cnr,status='keep')
      return
C
C     Error messages
C
 40   write(LINES,85) 'NZ_IN_OBJ', NZ_IN_OBJ
      NLINES = 4
      IERR = -4
      return
 41   write(LINES,85) 'NY_IN_OBJ', NY_IN_OBJ
      NLINES = 4
      IERR = -4
      return
 42   write(LINES,85) 'NE_IN_OBJ', NE_IN_OBJ
      NLINES = 4
      IERR = -4
      return
 43   write(LINES,85) 'NP_IN_OBJ', NP_IN_OBJ
      NLINES = 4
      IERR = -4
      return
C
 44   write(LINES,90) 'Z_INDEX', j, Z_INDEX(j)
      NLINES = 4
      IERR = -1
      return
 45   write(LINES,90) 'Y_INDEX', j, Y_INDEX(j)
      NLINES = 4
      IERR = -1
      return
 46   write(LINES,90) 'E_INDEX', j, E_INDEX(j)
      NLINES = 4
      IERR = -1
      return
 47   write(LINES,90) 'P_INDEX', j, P_INDEX(j)
      NLINES = 4
      IERR = -1
      return
 48   write(LINES,90) 'NDEGU', 0, NDEGU
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
     4     6x,'*                "*.ojf"                 *'/
     5     6x,'*         in your directory and          *'/
     6     6x,'*     specified for current problem      *'/
     7     6x,'******************************************')  
 80   format(
     1     6x,'******************************************'/
     2     6x,'*        "*.ojf" file not found          *'/
     3     6x,'******************************************')
 85   format(
     1     6x,'**********************************************'/
     2     6x,'*         Invalid value in .obf file:        *'/
     3     6x,'*       ' , a9 ,  '   =   ',   i7    , '.    *'/
     7     6x,'**********************************************')  
 90   format(
     1     6x,'******************************************'/
     2     6x,'*  Index out of bounds in "*.ojf:" file  *'/
     3     6x,'*        '  ,a7,'(',i7,')=',i7, ' .      *'/
     7     6x,'******************************************')  
      end
