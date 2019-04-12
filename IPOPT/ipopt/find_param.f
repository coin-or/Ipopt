C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: find_param.f 668 2005-02-22 15:05:12Z andreasw $
      integer function FIND_PARAM(NPARAMS, CPARAMS, NAME)
C
C     Find the entry NAME in the list of names CPARAMS (case insensitive)
C     Returns 0, if NAME was not found, and index in CPARAMS, if NAME was found.
C
C     Author:   Andreas Waechter    09-10-01
C
      implicit none
      integer NPARAMS           ! length of list CPARAMS
      character*(*) CPARAMS(NPARAMS) ! list in which NAME is to be found
      character*(*) NAME        ! entry to be found in CPARAMS

      logical FSTRLCMP

      do FIND_PARAM = 1, NPARAMS
         if( FSTRLCMP(CPARAMS(FIND_PARAM),NAME)) return ! found
      enddo
      FIND_PARAM = 0
      return
      end
C
C
C
      logical function FSTRLCMP(STR1, STR2)
C
C     Compare two character strings (STR1 and STR2) case-insensitive.
C     Returns .true., if STR1 and STR2 are identical; otherwise .false.
C
C     Author:  Andreas Waechter   09-10-01
C
      implicit none
      character*(*) STR1, STR2
C
      integer SLEN, i, length, pos1, pos2
      character*37 uc
      character*26 lc
      save uc, lc
      data uc /'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_'/
      data lc /'abcdefghijklmnopqrstuvwxyz'/
C
      FSTRLCMP = .false.
      length = SLEN(STR1)
      if( length.ne.SLEN(STR2) ) return
      do i = 1, length
         pos1 = index(uc,STR1(i:i))
         if( pos1.eq.0 ) pos1 = index(lc,STR1(i:i))
         pos2 = index(uc,STR2(i:i))
         if( pos2.eq.0 ) pos2 = index(lc,STR2(i:i))
         if( pos1.ne.pos2 .or. pos1.eq.0 ) return
      enddo
      FSTRLCMP = .true.
      return
      end
