C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine TIMER(TOTUSERTIME)
C
C*******************************************************************************
C
C    $Id: timer.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Simple program for measuring user time since program start.
CT    Probably needs to be changed for different operating systems.
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      05/01/02  Release as version IPOPT 2.0
C
C-------------------------------------------------------------------------------
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
CP   TOTUSERTIME O  INT    total user time since program start
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      double precision TOTUSERTIME
C
C-------------------------------------------------------------------------------
C                            Local variables
C-------------------------------------------------------------------------------
C
#include <config_f.h>

      real time
      real tarray(2)
      real ETIME
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
#ifdef HAVE_CPU_TIME
      call CPU_TIME(time)
      TOTUSERTIME = time
#else
# ifdef HAVE_ETIME
      TOTUSERTIME = ETIME(tarray)
      TOTUSERTIME = tarray(1)
# else
      TOTUSERTIME = 0.d0
# endif
#endif
      return
      end
