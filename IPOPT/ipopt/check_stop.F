C Copyright (C) 2004, International Business Machines and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      integer function CHECK_STOP()
C
C*******************************************************************************
C
C    $Id: check_stop.F 574 2004-04-25 22:42:37Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Checkes whether the maximal computation time is exceeded,
CT    a STOP file exists, or a user provided function returns TRUE
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
CA    Andreas Waechter      04/23/04  First version
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
CP   CHECK_STOP O   INT    =20 : Stop file detected
CP                         =21 : User provided USER_REQUESTED_STOP returns TRUE
CP                         = 3 : CPU time limit exceeded
CP                           0 : otherwise
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

#include <config_f.h>
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C

C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision cpusec
      logical ex
#ifdef USE_USER_REQUESTED_STOP
      logical USER_REQUESTED_STOP
#endif
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      CHECK_STOP = 0
C
C     Check if there is a STOP file
C
      inquire(file='STOP', exist=ex)
      if( ex ) then
         call C_OUT(2,0,1,'File ''STOP'' detected.  Abort.')
         CHECK_STOP = 20
         goto 9999
      endif
#ifdef USE_USER_REQUESTED_STOP
C
C     Check if USER_REQUESTED_STOP returns TRUE
C
      if( USER_REQUESTED_STOP() ) then
         call C_OUT(2,0,1,'User requested STOP  Abort.')
         CHECK_STOP = 21
         goto 9999
      endif
#endif

C
C     Check if we are running out of CPU time
C
      call TIMER(cpusec)
      if( cpusec.ge.dble(QMAXCPUSEC) ) then
         call C_OUT(2,0,1,'CPU time limit exceeded  Abort.')
         CHECK_STOP = 3
         goto 9999
      endif

 9999 continue

      return
      end
