C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine C_OUT(IWHERE, ILEVEL, NLINES, LINES)
C
C*******************************************************************************
C
C    $Id: c_out.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Output of a line on screen and into file
CT    Also responsible for opening the output file
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      05/01/02  Release as version IPOPT 2.0
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   IWHERE    I    I      0: Output to screen only
CP                         1: Output to file only (to unit QCNR, if nonzero)
CP                         2: Output to screen and file
CP                        -1: Open file for output (other inputs ignored)
CP                        -2: Close output file (other inputs ignored)
CP   ILEVEL    I    I      Output only, if QPRINT >= ILEVEL
CP   NLINES    I    I      Number of lines to be printed
CP                          (if 0, print empty line)
CP   LINES     I    C*(*)  lines of message to be printed
CP                          (spaces at end of a line are cut off)
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CS
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
      include 'IPOPT.INC'
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer IWHERE
      integer ILEVEL
      integer NLINES
      character *(*) LINES(NLINES)
C
C-------------------------------------------------------------------------------
C                             Local Variables
C-------------------------------------------------------------------------------
C
      integer SLEN
      integer sl, i
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( IWHERE.eq.-1 ) then
         open(QCNR, file='IPOPT.OUT', status = 'unknown' )
         return
      elseif( IWHERE.eq.-2 ) then
         close(QCNR)
         return
      endif

      if( ILEVEL.gt.QPRINT ) return

      if( NLINES.le.0 ) then
C
C     Output to screen
C
         if( IWHERE.eq.0 .or. IWHERE.eq.2 ) then
            write(*,8000)
         endif
C
C     Output to file
C
         if( QCNR.gt.0 .and. (IWHERE.eq.1 .or. IWHERE.eq.2 ) ) then
            write(QCNR,8000)
         endif

      else

         do i = 1, NLINES
C
C     Determine actual lenght of string without spaces at end
C
            sl = SLEN(LINES(i))
C
C     Output to screen
C
            if( IWHERE.eq.0 .or. IWHERE.eq.2 ) then
               if( sl.gt.0 ) then
                  write(*,8000) LINES(i)(:sl)
               else
                  write(*,8000)
               endif
            endif
C
C     Output to file
C
            if( QCNR.gt.0 .and. (IWHERE.eq.1 .or. IWHERE.eq.2 ) ) then
               if( sl.gt.0 ) then
                  write(QCNR,8000) LINES(i)(:sl)
               else
                  write(QCNR,8000)
               endif
            endif

         enddo

      endif

      return

 8000 format(a)
      end
