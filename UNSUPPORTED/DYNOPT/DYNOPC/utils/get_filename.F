C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: get_filename.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine GET_FILENAME(SUBDIR, NAME, EXT, FILENAME)
!DEC$ ATTRIBUTES DLLEXPORT :: GET_FILENAME
C
C     Determine FILENAME as a comination of
C     ${PROCESS_DATA_DIR}/SUBDIR/NAME.EXT
C     where PROCESS_DATA_DIR is the value of the enrivonment variable with
C     the same name.
C
C     Author:  Andreas Waechter     09-03-01
C
      implicit none

      character*1 slash
#ifdef WIN32
      parameter( slash = '\\' ) ! for Windows
#else
      parameter( slash = '/' )  ! for UNIX
#endif
      character*(*) SUBDIR, NAME, EXT, FILENAME
C
      character*(1024) process_data_dir
      integer SLEN
C
      call GETENV( 'PROCESS_DATA_DIR', FILENAME)

      if( FILENAME.ne.' ' ) then
         FILENAME = FILENAME(1:SLEN(FILENAME))//slash//SUBDIR
      else
         FILENAME = SUBDIR
      endif
      FILENAME = FILENAME(1:SLEN(FILENAME))//slash//NAME
      if( EXT.ne.' ' ) then
         FILENAME = FILENAME(1:SLEN(FILENAME))//'.'//EXT
      endif

      return
      end
