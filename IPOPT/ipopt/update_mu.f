C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine UPDATE_MU(FIRST_MU, RESTO, MU, ERR_BAR)
C
C*******************************************************************************
C
C    $Id: update_mu.f 632 2004-07-20 23:52:15Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Update barrier parameter and determine tolerance for barrier problem
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
CP   FIRST_MU I/O   LOG    I: if true, only change ERR_BAR
CP                         O: set to false
CP   RESTO     I    INT    <>0 if in restoration phase (in that case don't
CP                           change AMPLMU)
CP   MU       I/O   DP     barrier parameter
CP   ERR_BAR   O    DP     error for barrier problem
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
      logical FIRST_MU
      integer RESTO
      double precision MU
      double precision ERR_BAR
C
      double precision mumin
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( .not.FIRST_MU ) then
         mumin = QTOL/(2.d0*(QMUERRFAC+1.d0))
         mumin = QTOL/10.d0
         if( MU.gt.mumin ) then
            MU = dmin1(QMULIN*MU, MU**QMUSUPER)
            MU = max(MU,mumin)
         else
            MU = dmin1(QMULIN*MU, MU**QMUSUPER)
         endif
      endif
      ERR_BAR = dmin1(QMAXERR, QMUERRFAC*MU)
C
C     Also increase QTAU to allow full steps at the end
C
      QTAU = max(QTAUMAX, 1.d0-MU)

      FIRST_MU = .false.

      if( RESTO.eq.0 ) AMPLMU = MU

      return
      end
