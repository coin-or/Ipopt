C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      integer function SLEN(STRING)
C
C*******************************************************************************
C
C    $Id: slen.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Simple function that determines the length of a FORTRAN string
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      05/01/02  Release as version IPOPT 2.0
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C     Name     I/O   Type   Meaning
CP    STRING    I    C*(*)  String whos length is to be determined
CP    SLEN      O    I      Position of last non-space character in STRING
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
      character*(*) STRING
C
C-------------------------------------------------------------------------------
C                            Local variables
C-------------------------------------------------------------------------------
C

C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( STRING.eq.' ' ) then
         SLEN = 0
      else
         do SLEN = len(STRING), 1, -1
            if( STRING(slen:slen).ne.' ' ) goto 150
         enddo
      endif
 150  continue
      return
      end
