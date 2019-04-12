C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: eval_hesslag_v_dummy.F 531 2004-03-11 01:31:07Z andreasw $

      subroutine EV_HLV_DUMMY(TASK, N, X, M, LAM, VIN, VOUT, DAT, IDAT)
#ifdef BUILD_VSDLL
CDEC$ ATTRIBUTES DLLEXPORT::EV_HLV_DUMMY
#endif
      implicit none
      integer TASK, N, M
      double precision LAM(M), VIN(N), X(N), VOUT(N)
      double precision DAT(*)
      integer IDAT(*)

      write(*,*)
     1  'Wrong choice of IPOPT parameters: Got into EVAL_HESSLAG_V...'
      stop

      return
      end
