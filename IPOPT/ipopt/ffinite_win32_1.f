C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C $Id: ffinite_win32_1.f 531 2004-03-11 01:31:07Z andreasw $
C
C Fortran version of finite.c
C
C Author: Andreas Waechter    10-24-01
C
C (at this point only for WIN32)
C
      integer function FFINITE(number)
C
      implicit none
      double precision number
      include 'fordef.for'
      integer class
C
      class = FP_CLASS(number)
      if( class.eq.FOR_K_FP_SNAN .or. class.eq.FOR_K_FP_QNAN
     1   .or. class.eq.FOR_K_FP_POS_INF .or.
     1   class.eq.FOR_K_FP_NEG_INF ) then
         FFINITE = 0
      else
         FFINITE = 1
      endif
      return
      end
