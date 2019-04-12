C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine MA28_CALL(TASK, N, M, NZ, A, IROW, ICOL, TRANSP,
     1                     NRHS, RHS, IVAR, LRS, RS, LIS, IS, LRW,
     3                     RW, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: ma28_call.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface to MA28
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
CA    Andreas Waechter      05/01/02  Release as version IPOPT 2.0
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
CP   TASK      I    INT    information about what to do:
CP                          =0: initialize, tell LRS, LIS, LIW, LRW
CP                          =1: factorize nonsquare matrix in order to
CP                              get partition into dependent and independent
CP                              variables
CP                          =2: factorize square matrix
CP                              and solve for multiple rhs;
CP                              if necessary, new structural factorization
CP                              will be done
CP                              (It is assumed that a new pivot sequence
CP                               needs to be determined, if the previous call
CP                               was done with TASK = 0 or TASK = 1, but it
CP                               is assumed that otherwise the old pivot
CP                               sequence exists and can be used, unless
CP                               MA28BD complains...)
CP                          =3: use old factorization in order to solve for
CP                              multiple rhs
CP   N         I    INT    total number of variables (only TASK = 0,1)
CP   M         I    INT    number of constraints = number of depentent vars
CP   NZ        I    INT    number of nonzero elements in A
CP   A         I    DP     TASK =1,2: nonzero elements of matrix
CP                                     (unchanged on exit)
CP                            (note: NZ is different for TASK=1 and others)
CP   IROW      I    INT    TASK =1,2: row indices for A (unchanged on exit)
CP   ICOL      I    INT    TASK =1,2: col indices for A (unchanged on exit)
CP   TRANSP    I    INT    (only for TASK = 2,3)
CP                            = 1: solve system A*x     = rhs
CP                            <>1: solve system (A**T)x = rhs
CP   NRHS      I    INT    TASK =2,3: number of RHS
CP   RHS      I/O   DP     TASK =2,3: I: right hand sides
CP                                    O: solutions
CP   IVAR      O    INT    TASK = 1: IVAR(  1,..,M) containts the indices of
CP                                   a set of linear independent columns of A
CP                                   IVAR(M+1,..,N) containts the indices of the
CP                                   remaining column indices
CP                                   (Note: This does NOT include the fixed
CP                                          variables - those are taken care of
CP                                          in partition.f)
CP   LRS      I/O   INT    length of RS (Output for TASK = 0)
CP   RS       I/O   DP     can be used to store DP variables between calls;
CP                            this array is not touched from the calling
CP                            program
CP   LIS      I/O   INT    length of IS (Output for TASK = 0)
CP   IS       I/O   INT    can be used to store INT variables between calls;
CP                            this array is not touched from the calling
CP                            program
CP   LIW       I    INT    length of IW (Output for TASK = 0)
CP   IW        W    INT    integer work space
CP   LRW       I    INT    length of RW (Output for TASK = 0)
CP   RW        W    DP     double precision work space
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
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
CS    DCOPY
CS    C_OUT
CS    MA28AD
CS    MA28BD
CS    MA28CD
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
      integer TASK
      integer N
      integer M
      integer NZ
      double precision A(NZ)
      integer IROW(NZ)
      integer ICOL(NZ)
      integer TRANSP
      integer NRHS
      double precision RHS(M,NRHS)
      integer IVAR(N)
      integer LRS
      double precision RS(LRS)
      integer LIS
      integer IS(LIS)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            COMMON blocks
C-------------------------------------------------------------------------------
C
      integer LP, MP, IRNCP, ICNCP, MINIRN, MINICN, IRANK
      logical LBLOCK, GROW, ABORT1, ABORT2
      double precision EPS, RMIN, RESID

      COMMON/MA28ED/LP, MP, LBLOCK, GROW
      COMMON/MA28FD/EPS, RMIN, RESID, IRNCP, ICNCP, MINIRN, MINICN,
     *             IRANK, ABORT1, ABORT2
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer LIRN, LICN, P_A, P_ICN, P_IKEEP
      logical HAVE_PIV
      save    LIRN, LICN, P_A, P_ICN, P_IKEEP, HAVE_PIV

      double precision u
      integer i, k, ii, l, j, iflag, nind
      integer p_iwend, p_rwend, p_irn, p_iw, p_w, p_ivect, p_jvect
      integer liw1, liw2, liw3, liw4, lrw1, lrw2, lrw3, lrw4

      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
#ifndef HAVE_MA28
      if( TASK.eq.0 ) then
         LRS = 0
         LIS = 0
         IERR = 0
      else
         IERR = 97
         call C_OUT(2,0,1,
     1'ma28_call:  Compiled without MA28 - cannot execute this option')
      endif
      goto 9999
#else
      p_iwend = 0
      p_rwend = 0
      IERR = 0
C
C     Set MA28 COMMON block parameters
C
      LBLOCK = .false.
      ABORT1 = .true.
      ABORT2 = .true.

      goto (10, 100, 200, 300) TASK+1
      write(line,*) 'ma28_call: wrong TASK = ',TASK
      call C_OUT(2,0,1,line)
      stop

C-------------------------------------------------------------------------------
C     Start: Initialize
C-------------------------------------------------------------------------------
 10   continue

      HAVE_PIV = .false.

C      LICN = 10*NZ
C      LIRN = 4*NZ
      LICN = 40*NZ
      LIRN = 40*NZ
C     Determine storage space
      P_A = 0
      LRS = P_A       + LICN

      P_ICN = 0
      P_IKEEP = P_ICN + LICN
C     IKEEP is reserved a bit generous (5*M is enough for TASK<>1 )
      LIS = P_IKEEP + 5*N

C     Determine work space need
C     TASK = 1
      liw1 = LIRN + 8*N
      lrw1 = N
C     TASK = 2
      liw2 = LIRN + 8*M
      lrw2 = M
C     TASK = 3
      liw3 = 2*NZ + 5*M
      lrw3 = M
C     TASK = 4
      liw4 = M
      lrw4 = 0

      LIW = max0(liw1, liw2, liw3, liw4)
      LRW = max0(lrw1, lrw2, lrw3, lrw4)

C-------------------------------------------------------------------------------
C     End:   Initialize
C-------------------------------------------------------------------------------
      goto 9999

C-------------------------------------------------------------------------------
C     Start: Partitioning
C-------------------------------------------------------------------------------
 100  continue

      HAVE_PIV = .false.
C
C     Get work space pointers
C
      p_irn   = p_iwend
      p_iw    = p_irn   + LIRN
      p_iwend = p_iw    + 8*N

      p_w     = p_rwend
      p_rwend = p_w     + N

      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Copy A, IROW, ICOL into work space (delete old factorization!)
C
      call DCOPY(NZ, A,    1, RS(P_A+1)  , 1)
      do i = 1, NZ
         IW(p_irn+i) = IROW(i)
         IS(P_ICN+i) = ICOL(i)
      enddo
C
C     Do the factorization
C
      ABORT1 = .false.
      ABORT2 = .false.
CTODO
      u = 1.d-1
      u = QPIVTOL
C      u = 5.d-5
      call MA28AD(N, NZ, RS(P_A+1), LICN, IW(p_irn+1), LIRN,
     1            IS(P_ICN+1), u, IS(P_IKEEP+1), IW(p_iw+1),
     1            RW(p_w+1), iflag)
      if( iflag.lt.0 ) then
         write(line,*) 'ma28_call: MA28AD returns iflag = ',iflag
         call C_OUT(2,0,1,line)
         IERR = 514
         goto 9999
      endif
C
C     Get the partitioning out of IKEEP
C
      k = 0
      do i = 1, N
         ii = IS(P_IKEEP+2*N+i)
         if( ii.lt.0 ) then
C           indepentent variable
            k = k + 1
            IS(P_IKEEP+k) = -ii
         endif
      enddo
      nind = N - M
      if( k.gt.nind ) then
         call C_OUT(2,0,1,'Error in partitioning (MA28_CALL):')
         write(line,101) k-nind
 101     format('There seems to be ',i8,' dependent constraint(s).')
         call C_OUT(2,0,1,line)
         if( QPRINT.ge.1 ) then
            call C_OUT(1,0,1,' ')
            call C_OUT(1,0,1,'List of dependent constraints:')
            call C_OUT(1,0,1,' ')
            do i = M-k+nind+1, M
 102           write(line,'(i8)') IS(P_IKEEP+N+i)
               call C_OUT(1,0,1,line)
            enddo
         endif
         IERR = 6
         goto 9999
      endif
      k = M
      l = 0
      do i = 1, N
         do j = 1, nind
            if( i.eq.IS(P_IKEEP+j) ) then
               k = k + 1
               IVAR(k) = i
               goto 110
            endif
         enddo
         l = l + 1
         IVAR(l) = i
 110     continue
      enddo
C-------------------------------------------------------------------------------
C     End:   Partitioning
C-------------------------------------------------------------------------------
      goto 9999

 200  continue
      if( HAVE_PIV ) then
         goto 250
      else
         goto 210
      endif

C-------------------------------------------------------------------------------
C     Start: Structural factorization
C-------------------------------------------------------------------------------
 210  continue

      HAVE_PIV = .true.
C
C     Get work space pointers
C
      p_irn   = p_iwend
      p_iw    = p_irn   + LIRN
      p_iwend = p_iw    + 8*M

      p_w     = p_rwend
      p_rwend = p_w     + M

      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Copy A, IROW, ICOL into work space (delete old factorization!)
C
      call DCOPY(NZ, A,    1, RS(P_A+1)  , 1)
      do i = 1, NZ
         IW(p_irn+i) = IROW(i)
         IS(P_ICN+i) = ICOL(i)
      enddo
C
C     Do the factorization
C
      u = 1.d-1
CCC      u = 2.d-2
      call MA28AD(M, NZ, RS(P_A+1), LICN, IW(p_irn+1), LIRN,
     1            IS(P_ICN+1), u, IS(P_IKEEP+1), IW(p_iw+1),
     1            RW(p_w+1), iflag)
      if( iflag.ne.0 ) then
         write(line,*) 'ma28_call: MA28AD returns iflag = ',iflag
         call C_OUT(2,0,1,line)
         IERR = 512
         goto 9999
      endif
C
C     Free work space and solve for the right hand sides
C
      p_iwend = p_irn
      p_rwend = p_w
      goto 300

C-------------------------------------------------------------------------------
C     End:   Structural factorization
C-------------------------------------------------------------------------------

C-------------------------------------------------------------------------------
C     Start: "Soft" factorization
C-------------------------------------------------------------------------------
 250  continue
C
C     Get work space pointers
C
      p_ivect = p_iwend
      p_jvect = p_ivect + NZ
      p_iw    = p_jvect + NZ
      p_iwend = p_iw    + 5*M

      p_w     = p_rwend
      p_rwend = p_w     + M

      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Copy A, IROW, ICOL into work space (delete old factorization!)
C
      call DCOPY(NZ, A, 1, RS(P_A+1), 1)
      do i = 1, NZ
         IW(p_ivect+i) = IROW(i)
         IW(p_jvect+i) = ICOL(i)
      enddo
C
C     Do the factorization
C
      call MA28BD(M, NZ, RS(P_A+1), LICN, IW(p_ivect+1), IW(p_jvect+1),
     1            IS(P_ICN+1), IS(P_IKEEP+1), IW(p_iw+1), RW(p_w+1),
     1            iflag)
      if( iflag.gt.0 ) then
C
C        Need new pivot sequence; free work space and redo
C        structural factorization
C
         write(line,777) iflag
 777     format('ma28_call: MA28BD iflag = ',i5,
     1          '. Get new pivot sequence from MA28AD.')
         call C_OUT(2,0,1,line)
         p_iwend = p_ivect
         p_rwend = p_w
         goto 210
      elseif( iflag.ne.0 ) then
         write(line,*) 'ma28_call: MA28BD returns iflag = ',iflag
         call C_OUT(2,0,1,line)
         IERR = 513
         goto 9999
      endif
C
C     Free work space and solve for the right hand sides
C
      p_iwend = p_irn
      p_rwend = p_w
      goto 300

C-------------------------------------------------------------------------------
C     End:   "Soft" factorization
C-------------------------------------------------------------------------------

C-------------------------------------------------------------------------------
C     Start: Solve for multiple right hand sides
C-------------------------------------------------------------------------------
 300  continue
C
C     Get work space pointers
C
      p_w     = p_rwend
      p_rwend = p_w     + M

      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
C     Solve for the multiple rhs individually
C
      do i = 1, NRHS
         call MA28CD(M, RS(P_A+1), LICN, IS(P_ICN+1), IS(P_IKEEP+1),
     1               RHS(1,i), RW(p_w+1), TRANSP)
      enddo

      goto 9999

C-------------------------------------------------------------------------------
C     End:   solve for multiple right hand sides
C-------------------------------------------------------------------------------

#endif
 9999 continue
      return
      end
