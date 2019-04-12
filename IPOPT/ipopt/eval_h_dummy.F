C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: eval_h_dummy.F 531 2004-03-11 01:31:07Z andreasw $

      subroutine EV_H_DUMMY(TASK, N, X, M, LAM, NNZH, HESS, IRNH, ICNH,
     1     DAT, IDAT)
#ifdef BUILD_VSDLL
CDEC$ ATTRIBUTES DLLEXPORT::EV_H_DUMMY
#endif

      implicit none
      integer TASK, N, M, NNZH
      double precision LAM(M), X(N), HESS(NNZH)
      integer IRNH(NNZH), ICNH(NNZH)
      double precision DAT(*)
      integer IDAT(*)

      write(*,*) 'Wrong choice of IPOPT parameters: Got into EVAL_H...'
      stop

      return
      end
