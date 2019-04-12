C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine MOD_CHOL(N, AP, IPERM, B, NMOD, W)
C
C*******************************************************************************
C
C    $Id: mod_chol.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Implements modified Cholesky factorization as described in
CT    Nodeal/Wright's book (p.148) and solves for a vector
CT    Implemented in two subroutines, one for factorization and one for solving
CT    with a factorized matrix, which could be used individually
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
CP   N         I    INT    number of variables
CP   AD       I/O   DB     I: matrix to be factorized
CP                         O: contains L D L^T factorization
CP                            (see p.145-150 in Nocedal's book)
CP                            elements of D are on the diagonal
CP                            Unit-lower triangular L are the off-diagonal
CP                            elements
CP                           (this matrix is stored in upper-triangular packed
CP                            form!)
CP
CP   IPERM     O    INT    permutation matrix to ensure stability
CP                           ( X_reordered(i) = X_orig(IPERM(i)) )
CP   B        I/O   DP     I: right hand side of system
CP                         O: solution
CP   NMOD     I/O   INT    I: if -1 on input, stop factorization if modification
CP                            is necessary with NMOD = -1 on output
CP                         O: number of modifications necessary
CP                            (and see last line)
CP   W         W    DP     work space of size N
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
CS    IDAMAX
CS    DDOT
CS    D1MACH
CS    DSWAP
C
CS    DTPSV
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

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      double precision AP((N*(N+1))/2)
      integer IPERM(N)
      double precision B(N)
      integer NMOD
      double precision W(N)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C

C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      call MOD_CHOL_FAC(N, AP, IPERM, NMOD)
      call MOD_CHOL_SOL(N, AP, IPERM, B, W)
      return
      end
C
      subroutine MOD_CHOL_FAC(N, AP, IPERM, NMOD)
      implicit none
C     input parameters
      integer N
      double precision AP((N*(N+1))/2)
      integer IPERM(N)
      integer NMOD
C     local variables
      double precision gamma, xi, delta, beta, dj, macheps, theta, tmp
      integer i, j, k1, k2, jc, imax
C     functions
      double precision DDOT, D1MACH
      integer IDAMAX
C
      if( NMOD.ne.-1 ) NMOD = 0
C
C     First compute delta and beta (use formulas form Nocedal's book, page 149)
C
      gamma = dabs(AP(1))
      xi = 0d0
      k1 = 1
      do i = 2, N
         k2 = IDAMAX(i-1, AP(k1+1), 1)
         xi = dmax1(xi, AP(k1+k2))
         k1 = k1 + i
         gamma = dmax1(gamma, AP(k1))
      enddo
C
      macheps = D1MACH(4)
      delta = macheps*dmax1(gamma+xi, 1d0)
      beta = dsqrt(dmax1(gamma, xi/dsqrt(N**2-1d0), macheps))
C
C     Initialize permutation
C
      do i = 1, N
         IPERM(i) = i
      enddo
C
C     Now do actual factorization (use only memory of AP)
C
      jc = 0                    ! position in AP before first index in jth row
      do j = 1, N
C
C     Find pivot element
C
         imax = j
         k1 = jc + j            ! position of j-th diagonal element
         tmp = dabs(AP(k1))
         do i = j+1, N
            k1 = k1 + i
            if( tmp.lt.dabs(AP(k1)) ) then
               tmp = dabs(AP(k1))
               imax = i
            endif
         enddo
C
C     Permute rows and columns of necessary
C
         if( imax.ne.j ) then
C
            i = IPERM(j)
            IPERM(j) = IPERM(imax)
            IPERM(imax) = i
C
            i = 1 + ((imax*(imax-1))/2) ! start of imax-th row
            call DSWAP(j-1, AP(jc+1), 1, AP(i), 1)
C
            k1 = jc + j
            k2 = ((imax*(imax+1))/2)
            tmp = AP(k1)
            AP(k1) = AP(k2)
            AP(k2) = tmp
C
            k1 = k1 + j
            k2 = ((imax*(imax-1))/2) + j + 1
            do i = j+1, imax-1
               tmp = AP(k1)
               AP(k1) = AP(k2)
               AP(k2) = tmp
               k1 = k1 + i
               k2 = k2 + 1
            enddo
C
            k1 = (imax*(imax+1))/2
            do i = imax+1, N
               tmp = AP(k1+j)
               AP(k1+j) = AP(k1+imax)
               AP(k1+imax) = tmp
               k1 = k1 + i
            enddo
C
         endif
C
C     Compute j-th row of L
C
         k1 = 0
         do i = 1, j-1
            k1 = k1 + i
            AP(jc+i) = AP(jc+i)/AP(k1)
         enddo
C
C     Compute j-th column of C
C
         theta = 0d0
         k1 = jc + j
         do i = j+1, N
            if( j.gt.1 ) then
               AP(k1+j) = AP(k1+j) - DDOT(j-1, AP(jc+1), 1, AP(k1+1), 1)
            endif
            theta = dmax1(theta,dabs(AP(k1+j)))
            k1 = k1 + i
         enddo
C
C     Get j-th element of D
C
         dj = dmax1(dabs(AP(jc+j)),(theta/beta)**2,delta)
         if( dj.ne.AP(jc+j) ) then
            if( NMOD.eq.-1 ) return
            NMOD = NMOD + 1
         endif
         AP(jc+j) = dj
C
C     Correct remaining diagonal elements in C
C
         k1 = jc + j
         do i = j+1, N
            AP(k1+i) = AP(k1+i) - (AP(k1+j)**2)/dj
            k1 = k1 + i
         enddo
C
         jc = jc + j
C
      enddo
      if( NMOD.eq.-1 ) NMOD = 0
      return
      end
C
      subroutine MOD_CHOL_SOL(N, AP, IPERM, B, W)
      implicit none
C     input parameters
      integer N
      double precision AP((N*(N+1))/2) ! already factorized in LDL^T form
      integer IPERM(N)          ! permutation
      double precision B(N)     ! in: RHS out: solution
      double precision W(N)     ! work space
C     local variables
      integer i, ii
C
      do i = 1, N
         W(i) = B(IPERM(i))
      enddo
C
      call DTPSV('U', 'T', 'U', N, AP, W, 1)
C
      ii = 0
      do i = 1, N
         ii = ii + i
         W(i) = W(i) / AP(ii)
      enddo
C
      call DTPSV('U', 'N', 'U', N, AP, W, 1)
C
      do i = 1, N
         B(IPERM(i)) = W(i)
      enddo
C
      return
      end
