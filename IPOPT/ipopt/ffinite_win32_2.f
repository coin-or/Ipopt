C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C $Id: ffinite_win32_2.f 531 2004-03-11 01:31:07Z andreasw $
C
C Fortran version of finite.c
C
C Author: Andreas Waechter    10-24-01
C
C (at this point only for WIN32)
C
      integer function FFINITE(X)
C
      implicit none
      double precision X

      double precision zero
      zero = 0d0
      if( (X.eq.1d0/zero) .or. (X.eq.-1d0/zero) .or.
     1     (X.ne.X) ) then
         FFINITE = 0
      else
         FFINITE = 1
      endif

      return
      end
