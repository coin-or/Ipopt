C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      integer FUNCTION FILTER_CHECK(NFILTER, FILTER_C, FILTER_PHI,
     1     CNRM, PHI, PRECFACT, MACHEPS)
C
C*******************************************************************************
C
C    $Id: filter_check.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Check if point acceptable to filter
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
CP   FILTER_CHECK O INT    =0: point acceptable to filter
CP                         >0: number of filter entry excluding point
C
CP   NFILTER   I    INT    number of filter entries
CP   FILTER_C  I    DP     filter entries - theta part
CP   FILTER_PHI I   DP     filter entries - phi part
CP   CNRM      I    DP     theta value of point to check
CP   PHI       I    DP     objective function value of point to check
CP   PRECFACT  I    DP     precision factor for roundoffs
CP   MACHEPS   I    DP     machine precision
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

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer NFILTER
      double precision FILTER_C(NFILTER)
      double precision FILTER_PHI(NFILTER)
      double precision CNRM
      double precision PHI
      double precision PRECFACT
      double precision MACHEPS
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i
      double precision ci, fi, lhs
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      FILTER_CHECK = 0

      do i = 1, NFILTER
         ci = FILTER_C  (i)
         fi = FILTER_PHI(i)

         lhs = CNRM - ci
         if( lhs.le.PRECFACT*MACHEPS*ci ) goto 100

         lhs = PHI - fi
         if( lhs.le.PRECFACT*MACHEPS*dabs(fi) ) goto 100

         FILTER_CHECK = i
         goto 9999

 100     continue
      enddo

 9999 continue
      return
      end
