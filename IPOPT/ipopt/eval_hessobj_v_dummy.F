C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: eval_hessobj_v_dummy.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine EV_HOV_DUMMY(TASK, N, X, M, VIN, VOUT, DAT, IDAT)
#ifdef BUILD_VSDLL
CDEC$ ATTRIBUTES DLLEXPORT::EV_HOV_DUMMY
#endif
      implicit none
      integer TASK, N, M
      double precision VIN(N), X(N), VOUT(N)
      double precision DAT(*)
      integer IDAT(*)

      write(*,*)
     1  'Wrong choice of IPOPT parameters: Got into EVAL_HESSOBJ_V...'
      stop

      return
      end
