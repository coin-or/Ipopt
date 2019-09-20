C Copyright (C) 2002, 2007 Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Eclipse Public License.
C
C This code is based on the file Ipopt/IPOPT/ipopt/ma28_call.F in the
C Ipopt 2.2.1e distribution.
C
C This code provides an interface to MA28. Since Fortran COMMON blocks
C might be difficult to access from C++, we write this interface in Fortran
C and not C.
C
      subroutine MA28PART(TASK, N, M, NZ, A, IROW, ICOL, PIVTOL,
     1     FILLFACT, IVAR, NDEGEN, IDEGEN, LIW, IW, LRW, RW, IERR)
C
C*******************************************************************************
C
C    $Id$
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface to MA28 for detecting degenerate constraints
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
CA    Andreas Waechter      30/01/07  Based on Ipopt/IPOPT/ipopt/ma28_call.F
CA                                      in the Ipopt 2.2.1e distribution.
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
CP                          =0: initialize, tell LIW, LRW
CP                          =1: factorize nonsquare matrix in order to
CP                              get partition into dependent and independent
CP                              variables and find dependent constraints
CP   N         I    INT    total number of variables (only TASK = 0,1)
CP   M         I    INT    number of constraints = number of depentent vars
CP   NZ        I    INT    number of nonzero elements in A
CP   A         I    DP     TASK =1,2: nonzero elements of matrix
CP                                     (unchanged on exit)
CP                            (note: NZ is different for TASK=1 and others)
CP   IROW      I    INT    TASK =1,2: row indices for A (unchanged on exit)
CP   ICOL      I    INT    TASK =1,2: col indices for A (unchanged on exit)
CP   PIVTOL    I    DP     pivot tolerance  (e.g., 1.d-4 ?)
CP   FILLFACT  I    INT    estimated fillin factor (e.g. 40)
CP   IVAR      O    INT    TASK = 1: IVAR(  1,..,M) containts the indices of
CP                                   a set of linear independent columns of A
CP                                   IVAR(M+1,..,N) containts the indices of the
CP                                   remaining column indices
CP   NDEGEN    O    INT    Number of linearly dependent constraints
CP   IDEGEN    O    INT    List of linearly dependent constraints
CP   LIW       I    INT    length of IW (Output for TASK = 0)
CP   IW        W    INT    integer work space
CP   LRW       I    INT    length of RW (Output for TASK = 0)
CP   RW        W    DP     double precision work space
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
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
      double precision PIVTOL
      integer FILLFACT
      integer IVAR(N)
      integer NDEGEN
      integer IDEGEN(*)
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
      integer LIRN, LICN, p_a, p_icn, p_ikeep
      integer i, k, ii, l, j, iflag, nind, nsize
      integer p_iwend, p_rwend, p_irn, p_iw, p_w
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      LICN = FILLFACT*NZ
      LIRN = FILLFACT*NZ

      p_iwend = 0
      p_rwend = 0
      IERR = 0
C
C     Set MA28 COMMON block parameters
C
      LBLOCK = .false.
      ABORT1 = .true.
      ABORT2 = .true.

C     Allow for more constraints than variables
      nsize = MAX(N, M)

      goto (10, 100) TASK+1

C     Wrong argument for TASK
      IERR = -1
      return

C-------------------------------------------------------------------------------
C     Start: Compute work space requirement
C-------------------------------------------------------------------------------
 10   continue

C     Determine storage space
      LRW = LICN
      LIW = LICN + 5*nsize

C     TASK = 1
      LIW = LIW + LIRN + 8*nsize
      LRW = LRW + nsize

C-------------------------------------------------------------------------------
C     End:   Compute work space requirement
C-------------------------------------------------------------------------------
      goto 9999

C-------------------------------------------------------------------------------
C     Start: Partitioning
C-------------------------------------------------------------------------------
 100  continue

C
C     Get work space pointers
C
      p_icn   = p_iwend
      p_ikeep = p_icn   + LICN
      p_irn   = p_ikeep + 5*nsize
      p_iw    = p_irn   + LIRN
      p_iwend = p_iw    + 8*nsize

      p_a     = p_rwend
      p_w     = p_a     + LICN
      p_rwend = p_w     + nsize

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
      call DCOPY(NZ, A,    1, RW(p_a+1)  , 1)
      do i = 1, NZ
         IW(p_irn+i) = IROW(i)
         IW(p_icn+i) = ICOL(i)
      enddo
C
C     Do the factorization
C
      ABORT1 = .false.
      ABORT2 = .false.

      call MA28AD(nsize, NZ, RW(p_a+1), LICN, IW(p_irn+1), LIRN,
     1            IW(p_icn+1), PIVTOL, IW(p_ikeep+1), IW(p_iw+1),
     1            RW(p_w+1), iflag)
      if( iflag.lt.0 ) then
         IERR = 514
         goto 9999
      endif
C
C     Get the partitioning out of IKEEP
C
      k = 0
      do i = 1, N
         ii = IW(p_ikeep+2*N+i)
         if( ii.lt.0 ) then
C           indepentent variable
            k = k + 1
            IW(p_ikeep+k) = -ii
         endif
      enddo
      nind = N - M
      if( k.gt.nind ) then
C        get the dependent constraints
         NDEGEN = k-nind
         do i = 1, NDEGEN
            IDEGEN(i) = IW(p_ikeep+N+M-NDEGEN+i)
         enddo
      else
         NDEGEN = 0
      endif
      k = M
      l = 0
      do i = 1, N
         do j = 1, nind
            if( i.eq.IW(p_ikeep+j) ) then
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

 9999 continue
      return
      end
