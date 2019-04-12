C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      integer function READ_ENTRY(CNR, NEW, TYPE, DP, IN, ST)
C
C*******************************************************************************
C
C    $Id: read_entry.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Read entry from input file
C
C     Returns the next double precision, integer, or string read from a
C     file opened with channel number CNR.  If the file has just been opened
C     and this routine has not been called before for this file, you must set
C     NEW = .true., otherwise NEW = .false. (NEW is always set to .false. on
C     exit.)
C     TYPE specifies, what kind of entry has to be read, and the result - if
C     successful is stored in DB, IN, or ST, respectively.
C         'D' or 'd':  Double precision number  -> DB
C         'I' or 'i':  Integer number           -> IN
C         'S' or 's':  Character string         -> ST
C     Entries in the file are expected to be separated by commas ',' , or
C     new line characters.  In each line, everything after a '#' is ignored,
C     and empty lines (after ignoring '#') are overread.
C     NOTE: It is not possible to read two different files at the same time.
C
C     Return values:  0:  everything Ok
C                     1:  cannot read because end of file has been reached
C                    -1:  error while reading line from file
C                    -2:  error while converting entry into specified format
C                    -3:  invalue of  TYPE
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      05/01/02  Release as version IPOPT 2.0
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
      integer CNR
      logical NEW
      character*1 TYPE
      double precision DP
      integer IN
      character*(*) ST
C
C-------------------------------------------------------------------------------
C                            Local variables
C-------------------------------------------------------------------------------
C
      character*80 BUFFER     ! That's hopefully longer enough...
      integer BUFPOS, LASTPOS
      save BUFFER, BUFPOS, LASTPOS

      integer compos, nextpos, endpos, exppos, dotpos, firstpos
      character*80 entry
      integer SLEN
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C     Reset BUFFER if reading new file

      if( NEW ) BUFPOS = 0

      NEW = .false.

C     If next line has to be read, update buffer while ignoring all comments
C     and empty lines

 10   continue
      if( BUFPOS.eq.0 ) then
 20      continue
         read(CNR,8000,end=101,err=99) BUFFER

         compos = index(BUFFER,'#')

         if( compos.eq.1 ) then
            goto 20
         elseif( compos.gt.0 ) then
            BUFFER = BUFFER(1:compos-1)
         endif

         if( BUFFER.eq.' ' ) goto 20

C     Find last significant entry in BUFFER
         do LASTPOS = len(BUFFER)-1,1,-1
            if( BUFFER(LASTPOS:LASTPOS).ne.' ' ) goto 25
         enddo
 25      LASTPOS = LASTPOS + 1
      endif

C     Copy current entry and update BUFPOS

      nextpos = index(BUFFER(BUFPOS+1:LASTPOS),',') ! is zero, if current entry
                                                    ! is the last one in LINE

      firstpos = BUFPOS
 30   firstpos = firstpos + 1
      if( firstpos.lt.LASTPOS .and.
     1     BUFFER(firstpos:firstpos).eq.' ' ) goto 30

      nextpos = firstpos
 40   nextpos = nextpos + 1
      if( nextpos.lt.LASTPOS .and.
     1     BUFFER(nextpos:nextpos).ne.' ' .and.
     2     BUFFER(nextpos:nextpos).ne.',' ) goto 40

      if( BUFFER(nextpos:nextpos).eq.',' ) then
         endpos = nextpos - 1
      else
         endpos = nextpos
      endif
      entry = BUFFER(firstpos:endpos)

 50   continue
      if( nextpos.lt.LASTPOS .and.
     1     BUFFER(nextpos:nextpos).eq.' ' ) then
         nextpos = nextpos + 1
         goto 50
      endif

      if( BUFFER(nextpos:nextpos).eq.',' .and.
     1     nextpos.lt.LASTPOS ) then
 60      nextpos = nextpos + 1
         if( nextpos.lt.LASTPOS .and.
     1        BUFFER(nextpos:nextpos).eq.' ' ) then
            goto 60
         endif
      endif

C      if( nextpos.eq.LASTPOS .and.
C     1     (BUFFER(nextpos:nextpos).eq.' ' .or.
C     2      BUFFER(nextpos:nextpos).eq.',') ) then
      if( nextpos.eq.LASTPOS ) then
         BUFPOS = 0
      else
         BUFPOS = nextpos - 1
      endif

      if( TYPE.eq.'D' .or. TYPE.eq.'d' ) then
C
C     Read double precision value - since the format statement is very
C     sensitive for some Fortran compilers, see if changes are necessary
C     in order to obtain 123.45D67  format
C
         exppos = max(index(entry,'d'), index(entry,'D'),
     1        index(entry,'e'), index(entry,'E') )
         dotpos = index(entry,'.')

         if( exppos.gt.0 ) then
            entry(exppos:exppos) = 'D'
            if( dotpos.eq.0 ) entry(exppos:) = '.'//entry(exppos:)
         else
            if( dotpos.eq.0 ) then
               entry(slen(entry)+1:) = '.D0'
            else
               entry(slen(entry)+1:) = 'D0'
            endif
         endif
         if( dotpos.eq.1 ) entry = '0'//entry
         read(entry,8001,err=98) DP

      elseif( TYPE.eq.'I' .or. TYPE.eq.'i' ) then
C
C     Read integer value
C
         read(entry,8002,err=98) IN
      elseif( TYPE.eq.'S' .or. TYPE.eq.'s' ) then
         ST = entry
      else
         READ_ENTRY = -3
         return
      endif

      READ_ENTRY = 0
      return

 98   READ_ENTRY = -2
      return

 99   READ_ENTRY = -1
      return

 101  READ_ENTRY = 1
      return

 8000 format(A)
 8001 format(D23.17)
 8002 format(I30)
      end
