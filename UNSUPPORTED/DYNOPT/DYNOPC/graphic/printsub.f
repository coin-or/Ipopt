C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: printsub.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine PRINTSUB(CNR, TEXT)
c !DEC$ ATTRIBUTES DLLEXPORT :: PRINTSUB
C
C     A very 'dirty' way to force that 'WRITE' is executed in
C     main program and not one of the DLLs.....   :)
C     (See also C_OUT_DYNOPC.F)
C
C     Author:       09-27-01     Andreas Waechter
C
      implicit none
      integer CNR
      character*(*) TEXT
      integer SLEN, len

      len = SLEN(TEXT)
      if( len.gt.0 ) then
         write(CNR,'(a)') TEXT(:len)
      else
         write(CNR,'(a)')
      endif
      return
      end
C
      subroutine SET_PRINTSUB(PRINTSUB)
C
C     Setting the common block variable
C
      implicit none
      integer PRINTSUB
      include 'DYNAUX.INC'
	include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/
      PRINTSUBP = %LOC(PRINTSUB)
      return
      end
