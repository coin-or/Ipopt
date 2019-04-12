C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      SUBROUTINE waiting()
      use msflib
      integer(4) i4,ix,iy,i4k
!     external go_ahead
!  wait forever to allow event-driven action
      do while(.TRUE.)
      i4=waitonmouseevent(MOUSE$RBUTTONDOWN,i4k,ix,iy)
      end do
      end
