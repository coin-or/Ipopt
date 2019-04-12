C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: wait4click.f 531 2004-03-11 01:31:07Z andreasw $
      SUBROUTINE wait4click(flag)
      use dflib
	character*(*) flag
      integer(4) ix,iy,i4k,event

	if( flag.eq.'LR') then
! click mouse left or right button to go ahead
         event = MOUSE$RBUTTONDOWN
         event = IOR(event,MOUSE$LBUTTONDOWN )
	else
         event = MOUSE$RBUTTONDOWN
	endif
      i4=waitonmouseevent(event,i4k,ix,iy)
      end
