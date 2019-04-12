C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: get_amplmu.f 531 2004-03-11 01:31:07Z andreasw $

      subroutine GET_AMPLMU( MU )
C
C     Obtain the common block variable AMPLMU from /AMPLMU/
C
      implicit none
      include 'IPOPT.INC'
      double precision MU
C
      MU = AMPLMU
      return
      end
