C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C*******************************************************************************
C
      subroutine C_OUT(IWHERE, ILEVEL, NLINES, LINES)
!DEC$ ATTRIBUTES DLLEXPORT::C_OUT
C
C*******************************************************************************
C
C    $Id: c_out_dynopc.f 531 2004-03-11 01:31:07Z andreasw $
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
CA    Andreas Waechter      01/25/00
CA    Andreas Waechter      09/21/01   opening of file here
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
CP   ILEVEL    I    I      Output to file only, if QPRINT >= ILEVEL
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
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLEXPORT :: /PARAMS/, /DYNAUX/, /GRAPH/, /DYNOPC/
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
      integer i
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
      endif

      if( NLINES.le.0 ) then
C
C     Output to screen
C
         if( IWHERE.eq.0 .or. IWHERE.eq.2 ) then
            call mywrite(%VAL(PRINTSUBP),OutputUnit_COUT,' ')
         endif
C
C     Output to file
C
         if( QCNR.gt.0 .and. (QPRINT.ge.ILEVEL .or. ILEVEL.eq.0) .and. 
     1        (IWHERE.eq.1 .or. IWHERE.eq.2 ) ) then
            call mywrite(%VAL(PRINTSUBP),QCNR,' ')
         endif

      else

         do i = 1, NLINES
C
C     Output to screen
C
            if( IWHERE.eq.0 .or. IWHERE.eq.2 ) then
               call mywrite(%VAL(PRINTSUBP),
     1                         OutputUnit_COUT,LINES(i))
            endif
C     
C     Output to file
C     
            if( QCNR.gt.0 .and. QPRINT.ge.ILEVEL .and. 
     1           (IWHERE.eq.1 .or. IWHERE.eq.2 ) ) then
               call mywrite(%VAL(PRINTSUBP),QCNR,LINES(i))
            endif
            
         enddo

      endif

      return
      end
C
C     If this works .... :)  See also printsub.f
C
      subroutine MYWRITE(PRINTSUB, CNR, TEXT)
      implicit none
      external PRINTSUB
      integer CNR
      character*(*) TEXT
      call PRINTSUB(CNR,TEXT)
      return
      end
