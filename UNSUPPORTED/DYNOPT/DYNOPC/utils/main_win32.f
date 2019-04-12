C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: main_win32.f 531 2004-03-11 01:31:07Z andreasw $
      program main
C
C     Main program for simulation
C
C     Authors:       Andreas Waechter    09-26-01
C
      implicit none
      character*40 MODNAM, STUBNAM
C
C     Read model and run name from file 'MODEL.DAT' in current directory
C
      open(1,file='MODEL.DAT',status='old',err = 999)
      read(1,'(a)',err=999) MODNAM
      read(1,'(a)',err=999) STUBNAM
      close(1)

      call FMAIN(MODNAM, STUBNAM)
      stop

 999  write(*,*) 'Error reading file MODEL.DAT'
      stop
      end
