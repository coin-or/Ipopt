C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine MA48_CALL(TASK, N, M, NZ, A, IROW, ICOL, TRANSP,
     1                     NRHS, RHS, IVAR, LRS, RS, LIS, IS, LRW,
     3                     RW, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: ma48_call.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface to MA48
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
CP                          =1: [factorize nonsquare matrix in order to
CP                               get partition into dependent and independent
CP                               variables]
CP                              THIS IS NOT IMPLEMENTED!
CP                              THE ONLY EFFECT IS TO RESET 'HAVE_PIV=.false.'
CP                              (after a change of basis)
CP                          =2: factorize square matrix
CP                              and solve for multiple rhs;
CP                              if necessary, new structural factorization
CP                              will be done
CP                              (It is assumed that a new pivot sequence
CP                               needs to be determined, if the previous call
CP                               was done with TASK = 0 or TASK = 1.  One
CP                               might also speedup the procedure by
CP                               calling MA48BD with job=2, so that the pivot
CP                               sequence is not changed unless it returns
CP                               an error, and in such a case call it again
CP                               with job=1.  But for now, we always call it
CP                               with job=1 every time.)
CP                          =3: use old factorization in order to solve for
CP                              multiple rhs
CP   N         I    INT    total number of variables
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
CP   IVAR      O    INT    TASK = 1: IVAR(  1,...,M) containts the indices of
CP                                   a set of linear independent columns of A
CP                                   IVAR(M+1,...,N) containts the indices of
CP                                   the remaining column indices
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
CS    MA48AD
CS    MA48BD
CS    MA48CD
CS    MA48ID
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
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer LA, P_A, P_IRN, P_JCN, P_KEEP
      logical HAVE_PIV
      save    LA, P_A, P_IRN, P_JCN, P_KEEP, HAVE_PIV

      integer i, job, ii, j, ifail, rank
      integer p_iwend, p_rwend, p_w, p_iw, p_x, p_ikeep, p_icols
      integer liw1, liw2, liw3, liw4, lrw1, lrw2, lrw3, lrw4
      double precision cntl(5), error(3), rinfo
      integer icntl(9), info(12)
      logical transp_l, copied_a

      character*80 line(2)

      integer JOB_B
      save JOB_B
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
#ifndef HAVE_MA48
      if( TASK.eq.0 ) then
         LRS = 0
         LIS = 0
         IERR = 0
      else
         IERR = 97
         call C_OUT(2,0,1,
     1 'ma48_call:  Compiled without MA48 - cannot execute this option')
      endif
      goto 9999
#else
      p_iwend = 0
      p_rwend = 0
      IERR = 0
      copied_a = .false.
C
C     Initialize parameters for MA48
C
      call MA48ID(cntl,icntl)
      icntl(1) = QCNR   ! redirect output of MA48 into file
      icntl(2) = QCNR   ! redirect output of MA48 into file
C      icntl(3) = 4  ! write much output
      cntl(2) = QPIVTOL
      goto (10, 100, 200, 300) TASK+1
      write(line,*) 'ma48_call: wrong TASK = ', TASK
      call C_OUT(2,0,1,line)
      stop

C-------------------------------------------------------------------------------
C     Start: Initialize
C-------------------------------------------------------------------------------
 10   continue

      HAVE_PIV = .false.

C     If number of garbage collections is large, we should increase this!
      LA = 20*NZ

C     Determine storage space
      P_A     = 0
      LRS     = P_A     + LA

      P_IRN  = 0
      P_JCN  = P_IRN  + LA   ! JCN only needs to be that big for MA48AD
      P_KEEP = P_JCN  + LA
      LIS    = P_KEEP + 10*M+7

C     Determine work space need

      liw1 = 7*N + 12*M + 7
      lrw1 = N

      liw2 = 9*M
      lrw2 = 0

      liw3 = 4*M
      lrw3 = M

      liw4 = M
      lrw4 = 5*M

      LIW = max0(liw1, liw2, liw3, liw4)
      LRW = max0(lrw1, lrw2, lrw3, lrw4)

      JOB_B = 1                 ! Flag for MA48BD

C-------------------------------------------------------------------------------
C     End:   Initialize
C-------------------------------------------------------------------------------
      goto 9999

C-------------------------------------------------------------------------------
C     Start: Partitioning
C-------------------------------------------------------------------------------
 100  continue
C
C     First call MA48AD to perform the structural factorization.
C     It seems necessary to enter the matrix so that there are more
C     rows than columns
C

C
C     Get work space pointers
C
      p_ikeep = p_iwend
      p_iw    = p_ikeep + N + 9*M + 7
      p_iwend = p_iw + 6*N + 3*M
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Copy A, IROW, ICOL into work space (delete old factorization!)
C
      call DCOPY(NZ, A, 1, RS(P_A+1) , 1)
      do i = 1, NZ
         IS(P_IRN+i) = IROW(i)
         IS(P_JCN+i) = ICOL(i)
      enddo
C
C     Do the factorization
C
C      icntl(7) = 0              ! allow structurally deficient matrices
C      cntl(2) = 1d-1
      job = 1
      call MA48AD(N, M, NZ, job, LA, RS(P_A+1), IS(P_JCN+1),
     1     IS(P_IRN+1), IW(p_ikeep+1), cntl, icntl, IW(p_iw+1),
     2     info, rinfo)
      if( info(1).lt.0 ) then
         write(line,113) info(1)
 113     format ('ma48_call: MA48AD-1 returns info(1)=',i3)
         call C_OUT(2,0,1,line)
         if( info(1).eq.-4 ) then
            write(line,114) info(10)
 114        format('   Matrix is structurally rank deficient with rank',
     1           i4)
            call C_OUT(2,0,1,line)
         endif
         IERR = 6
         goto 9999
      elseif( info(1).gt.0 ) then
         write(line,111) info(1)
 111     format('ma48_call: MA48AD-1 returns info(1)=',i2)
         call C_OUT(2,1,1,line)
      endif
      if( info(2).gt.10 ) then
         write(line,112) info(2)
 112     format('ma48_call: Number of garbage collections=',i2,/,
     1          '           LA should be increased!')
         call C_OUT(2,0,2,line)
      endif
C
C     Free work space
C
      p_iwend = p_iw
C
C     Now call MA48BD for the numerical factorizatoin
C

C
C     Get work space pointers
C
      p_w     = p_rwend
      p_rwend = p_w     + N

      p_iw    = p_iwend
      p_iwend = p_iw    + 2*N + 2*M

      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Do the factorization
C
      job = 1
      call MA48BD(N, M, NZ, job, LA, RS(P_A+1), IS(P_JCN+1),
     1     IS(P_IRN+1), IW(p_ikeep+1), cntl, icntl, RW(p_w+1),
     1     IW(p_iw+1), info, rinfo)
      if( info(1).ne.0 .and. info(1).ne.2 ) then
         write(line,151) info(1)
 151     format ('ma48_call: MA48BD-1 returns info(1)=',i3)
         call C_OUT(2,0,1,line)
         if( info(1).eq.-3 ) then
            write(line,152) LA, info(4)
 152        format('ma48_call: LA =',i8,' too small. Make at least',i8)
            call C_OUT(2,0,1,line)
         endif
         IERR = 515
         goto 9999
      endif
C
C     Free work space
C
      p_iwend = p_iw
      p_rwend = p_w
C
C     Finally call MA51AD to get the partition information out of MA48AD's
C     work space.
C
      p_w     = p_iwend
      p_icols = p_w     + N
      p_iwend = p_icols + M

      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif

      call MA51AD(N, M, LA, IS(P_JCN+1), IW(p_ikeep+1), rank,
     1     IVAR, IW(p_icols+1), IW(p_w+1))

      if( rank.ne.M ) then
         call C_OUT(2,0,1,'Error in partitioning (MA48_CALL):')
         write(line,171) M-rank
 171     format('There seems to be ',i8,' dependent constraint(s).')
         call C_OUT(2,0,1,line)
         if( QPRINT.ge.1 ) then
            call C_OUT(1,0,1,' ')
            call C_OUT(1,0,1,'List of dependent constraints:')
            call C_OUT(1,0,1,' ')
            do i = rank+1, M
               write(line,'(i8)') IW(p_icols+i)
               call C_OUT(1,0,1,line)
            enddo
         endif
         IERR = 6
         goto 9999
      endif

      p_iwend = p_ikeep

      HAVE_PIV = .false.

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
C     Start: Structural factorization (using MA48AD)
C-------------------------------------------------------------------------------
 210  continue

      HAVE_PIV = .true.
C
C     Get work space pointers
C
      p_iw    = p_iwend
      p_iwend = p_iw + 9*M

      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Copy A, IROW, ICOL into work space (delete old factorization!)
C
      call DCOPY(NZ, A, 1, RS(P_A+1) , 1)
      copied_a = .true.
      do i = 1, NZ
         IS(P_IRN+i) = IROW(i)
         IS(P_JCN+i) = ICOL(i)
      enddo
C
C     Do the factorization
C
      job = 1
      call MA48AD(M, M, NZ, job, LA, RS(P_A+1), IS(P_IRN+1),
     1     IS(P_JCN+1), IS(P_KEEP+1), cntl, icntl, IW(p_iw+1),
     2     info, rinfo)
      if( info(1).lt.0 ) then
         write(line,213) info(1)
 213     format ('ma48_call: MA48AD returns info(1)=',i3)
         call C_OUT(2,0,1,line)
         if( info(1).eq.-4 ) then
            write(line,214) info(10)
 214        format('   Matrix is structurally rank deficient with rank',
     1           i4)
            call C_OUT(2,0,1,line)
         endif
         IERR = 6
         goto 9999
      elseif( info(1).gt.0 ) then
         write(line,211) info(1)
 211     format('ma48_call: MA48AD returns info(1)=',i2)
         call C_OUT(2,1,1,line)
      endif
      if( info(2).gt.10 ) then
         write(line,212) info(2)
 212     format('ma48_call: Number of garbage collections=',i2,/,
     1          '           LA should be increased!')
         call C_OUT(2,0,2,line)
      endif
      JOB_B = 1
C
C     Free work space and continue factorization
C
      p_iwend = p_iw

C-------------------------------------------------------------------------------
C     End:   Structural factorization
C-------------------------------------------------------------------------------

C-------------------------------------------------------------------------------
C     Start: Factorization
C-------------------------------------------------------------------------------
 250  continue
C
C     Copy A, IROW, ICOL into work space (delete old factorization!)
C
      if( .not.copied_a ) then
         call DCOPY(NZ, A, 1, RS(P_A+1), 1)
      endif
C
C     Get work space pointers
C
      p_w     = p_rwend
      p_rwend = p_w     + M

      p_iw    = p_iwend
      p_iwend = p_iw    + 4*M

      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Do the factorization
C
 260  continue
      call MA48BD(M, M, NZ, JOB_B, LA, RS(P_A+1), IS(P_IRN+1),
     1     IS(P_JCN+1), IS(P_KEEP+1), cntl, icntl, RW(p_w+1),
     1     IW(p_iw+1), info, rinfo)
      if( info(1).eq.-7 ) then
         call DCOPY(NZ, A, 1, RS(P_A+1), 1) ! I'm not sure if this is necessary
         JOB_B = 1
         goto 260
      endif
      JOB_B = 2
      if( info(1).ne.0 ) then
         write(line,251) info(1)
 251     format ('ma48_call: MA48BD returns info(1)=',i3)
         call C_OUT(2,0,1,line)
         if( info(1).eq.-3 ) then
            write(line,252) LA, info(4)
 252        format('ma48_call: LA =',i8,' too small. Make at least',i8)
            call C_OUT(2,0,1,line)
         endif
         IERR = 515
         goto 9999
      endif
C
C     Free work space and solve for the right hand sides
C
      p_iwend = p_iw
      p_rwend = p_w
      goto 300

C-------------------------------------------------------------------------------
C     End:   Factorization
C-------------------------------------------------------------------------------

C-------------------------------------------------------------------------------
C     Start: Solve for multiple right hand sides
C-------------------------------------------------------------------------------
 300  continue
C
C     Get work space pointers
C
      p_x     = p_rwend
      p_w     = p_x     + M
      p_rwend = p_w     + 4*M

      p_iw    = p_iwend
      p_iwend = p_iw    + M

      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Solve for the multiple rhs individually
C
      if( TRANSP.eq.1 ) then
         transp_l = .false.
      else
         transp_l = .true.
      endif
      job = 3   ! for now, with iterative refinement
CTODO?      job = 4   ! for now, with iterative refinement

      do i = 1, NRHS
         call MA48CD(M, M, transp_l, job, LA, RS(P_A+1), IS(P_IRN+1),
     1        IS(P_KEEP+1), cntl, icntl, RHS(1,i), RW(p_x+1), error,
     2        RW(p_w+1), IW(p_iw+1), info)
         if( info(1).ne.0 ) then
            write(line,301) info(1)
 301        format ('ma48_call: MA48CD returns info(1)=',i3)
            call C_OUT(2,0,1,line)
            IERR = 516
            goto 9999
         endif
         call DCOPY(M, RW(p_x+1), 1, RHS(1,i), 1)
      enddo

      goto 9999

C-------------------------------------------------------------------------------
C     End:   solve for multiple right hand sides
C-------------------------------------------------------------------------------

#endif
 9999 continue
      return
      end
