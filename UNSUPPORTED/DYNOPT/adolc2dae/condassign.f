C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine CONDASSIGN(A, B, C, D)
      double precision A, B, C, D
      if( B.gt.0.d0 ) then
         A = C
      else
         A = D
      endif
      return
      end
