C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine FILTER_ADD(NFILTER, FILTER_C, FILTER_PHI,
     1     NEW_CNRM, NEW_PHI)
C
C*******************************************************************************
C
C    $Id: filter_add.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Add antry to filter
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
CP   NFILTER  I/O   INT    number of filter entries
CP   FILTER_C I/O   DP     filter entries - theta part
CP   FILTER_PHI I/O DP     filter entries - phi part
CP   NEW_CNRM  I    DP     theta value to be added (incl. margin)
CP   NEW_PHI   I    DP     phi value to be added (incl. margin)
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
      integer NFILTER
      double precision FILTER_C(*)
      double precision FILTER_PHI(*)
      double precision NEW_CNRM
      double precision NEW_PHI
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i, ndelete
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Delete entries in filter, that are dominated by new entry
C
      ndelete = 0
      do i = 1, NFILTER
         if( FILTER_PHI(i).ge.NEW_PHI .and.
     1        FILTER_C(i) .ge.NEW_CNRM  ) then
            ndelete = ndelete + 1
         elseif( ndelete.gt.0 ) then
            FILTER_PHI(i-ndelete) = FILTER_PHI(i)
            FILTER_C  (i-ndelete) = FILTER_C  (i)
         endif
      enddo
C
C     Add new entry to filter
C
      NFILTER = NFILTER - ndelete + 1
      FILTER_PHI(NFILTER) = NEW_PHI
      FILTER_C  (NFILTER) = NEW_CNRM
      if( QCNR.gt.0 .and. QPRINT.gt.1 ) then
         write(line,*)
     1        'Number of entries deleted from filter: ',ndelete
         call C_OUT(1,2,1,line)
         write(line,*)
     1        'Number of entries now in filter: ',NFILTER
         call C_OUT(1,2,1,line)
         do i = 1, NFILTER
            write(line,900) i, FILTER_PHI(i), FILTER_C(i)
 900        format('Entry',i3,': Phi = ',d20.10,' cnrm = ',d20.10)
            call C_OUT(1,2,1,line)
         enddo
      endif

      return
      end
