C Copyright (C) 2003, International Business Machines and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine MC29_CALL(NTOT, NZKKT, A, IRN, JCN, SCALE,
     1     LRW, RW, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: mc29_call.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface to MC29 to obtain scaling factors for SYMMETRIC matrix
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      05/26/03  Initial implementation
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   NTOT      I    INT    total number of variables
CP   NZKKT     I    INT    number of nonzero elements in A
CP   A         I    DP     elements of matrix
CP   IRN       I    INT    row indices of matrix (MA27 format)
CP   JCN       I    INT    col indices of matrix (MA27 format)
CP   SCALE     O    DP     scaling factors
CP   LIW       I    INT    length of IW
CP   IW        W    INT    integer work space
CP   LRW       I    INT    length of RW
CP   RW        W    DP     double precision work space
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CS    MC29AD
CS    C_OUT
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
      include 'IPOPT.INC'

#include <config_f.h>
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer NTOT
      integer NZKKT
      double precision A(NZKKT)
      integer IRN(NZKKT)
      integer JCN(NZKKT)
      double precision SCALE(NTOT)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer p_iwend, p_rwend, p_kkt2, p_irn2, p_jcn2, p_r, p_c, p_w
      integer i, nzkkt2, ifail
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
#ifndef HAVE_MC29
      IERR = 97
      call C_OUT(2,0,1,
     1 'mc29_call:  Compiled without MC29 - cannot execute this option')
      goto 9999
#else
      IERR = 0
      p_rwend = 0
      p_iwend = 0

      if( MEMDBG ) then
         write(line,1)'mc29_call', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif
C
C     We first need to copy the matrix into non-symmetric format
C
      p_kkt2  = p_rwend
      p_rwend = p_kkt2 + 2*NZKKT ! worst estimate
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      p_irn2  = p_iwend
      p_jcn2  = p_irn2 + 2*NZKKT
      p_iwend = p_jcn2 + 2*NZKKT
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      nzkkt2 = 0
      do i = 1, NZKKT
         if( IRN(i).eq.JCN(i) ) then
            nzkkt2 = nzkkt2 + 1
            RW(p_kkt2+nzkkt2) = A(i)
            IW(p_irn2+nzkkt2) = IRN(i)
            IW(p_jcn2+nzkkt2) = JCN(i)
         else
            nzkkt2 = nzkkt2 + 1
            RW(p_kkt2+nzkkt2) = A(i)
            IW(p_irn2+nzkkt2) = IRN(i)
            IW(p_jcn2+nzkkt2) = JCN(i)
            nzkkt2 = nzkkt2 + 1
            RW(p_kkt2+nzkkt2) = A(i)
            IW(p_irn2+nzkkt2) = JCN(i)
            IW(p_jcn2+nzkkt2) = IRN(i)
         endif
      enddo
C
C     Call MC29 to get scaling factors for both rows and columns
C
      p_r     = p_rwend 
      p_c     = p_r + NTOT
      p_w     = p_c + NTOT
      p_rwend = p_w + 5*NTOT
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call MC29AD(NTOT, NTOT, nzkkt2, RW(p_kkt2+1), IW(p_irn2+1),
     1     IW(p_jcn2+1), RW(p_r+1), RW(p_c+1), RW(p_w+1), QCNR, ifail)
      if( ifail.ne.0 ) then
         write(line,*) 'mc29_call: MC29AD returned IFAIL = ',ifail
         call C_OUT(2,0,1,line)
         IERR = 416
         goto 9999
      endif
C
C     Get overall scaling factors as mean values
C
      do i = 1, NTOT
         SCALE(i) = (RW(p_r+i)+RW(p_c+i))/2.d0
      enddo
C
      p_rwend = p_kkt2
      p_iwend = p_irn2
#endif
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine MC29_CALL_WS(N, M, NLB, NUB, NZA, NZKKT, LRW, LIW)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, NZKKT, LRW, LIW
      character*80 line

      LRW = 2*NZKKT + 7*(N+M)
      LIW = 4*NZKKT

      if( QPRINT.ge.4 ) then
         write(line,1000) 'mc29_call', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end

