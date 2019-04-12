C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C    $Id: example.f 531 2004-03-11 01:31:07Z andreasw $
C
C =============================================================================
C
C     This is an example for the usage of IPOPT.
C     It implements problem 71 from the Hock-Schittkowsky test suite:
C
C     min   x1*x4*(x1 + x2 + x3)  +  x3
C     s.t.  x1*x2*x3*x4  -  x5             -  25  =  0
C           x1**2 + x2**2 + x3**2 + x4**2  -  40  =  0
C           1 <=  x1,x2,x3,x4  <= 5
C           0 <=  x5
C
C     Starting point:
C        x = (1, 5, 5, 1, -24)
C
C     Optimal solution:
C        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829, 0)
C
C =============================================================================
C
C
C =============================================================================
C
C                            Main driver program
C
C =============================================================================
C
      program example
C
      implicit none
C
C     Define work space:
C
C     If your Fortran compiler understand %VAL and IPOPT can use the C
C     function malloc (i.e. if USE_MALLOC is defined in your config.h),
C     the following variables have to be defined but are ignored.
C
      integer     LRW,     LIW
      parameter  (LRW = 0, LIW = 0)
      integer          IW
      double precision RW
C
C     Otherwise, you need to provide IPOPT with sufficient work space:
C
C      integer     LRW,         LIW
C      parameter  (LRW = 10000, LIW = 10000)
C      integer          IW(LIW)
C      double precision RW(LRW)

C
C     IPOPT is defined as double precision function
C     (returning final value of objective function)
C
      double precision IPOPT
C
C     Size of the problem (number of variables and equality constraints)
C
      integer     N,     M
      parameter  (N = 5, M = 2)
C
C     Space for multipliers and constraints
C
      double precision LAM(M)
      double precision C(M)
C
C     Vector of variables
C
      double precision X(N)
C
C     Number of lower and upper bounds
C
      integer     NLB,     NUB
      parameter  (NLB = 5, NUB = 4)
C
C     Vector of lower and upper bounds
C
      integer             ILB(NLB),    IUB(NUB)
      double precision BNDS_L(NLB), BNDS_U(NUB), V_L(NLB), V_U(NUB)
C
C     Private data for evaluation routines
C     This could be used to pass double precision and integer arrays untouched
C     to the evaluation subroutines EVAL_*
C
      double precision DAT(2)
      integer IDAT(1)
C
C     Algorithmic Parameters (INITPARAMS)
C
      integer NARGS
      double precision ARGS(50)
      character*20 CARGS(50)
C
      integer i, IERR, ITER
      double precision f
      integer lrwe, liwe
C
C     Set initial point and bounds:
C
      data X      / 1d0, 5d0, 5d0, 1d0, -24d0 /
      data ILB    /   1,   2,   3,   4,     5 /
      data BNDS_L / 1d0, 1d0, 1d0, 1d0,   0d0 /
      data IUB    /   1,   2,   3,   4 /
      data BNDS_U / 5d0, 5d0, 5d0, 5d0 /
C
C     The following are the Fortran routines for computing the model
C     functions and their derivatives - their code can be found furhter
C     down in this file.
C
      external EV_F, EV_C, EV_G, EV_A, EV_H
C
C     The following are Fortran routines that IPOPT requires in its
C     argument list, but which are not needed for the chosen option.
C     Those are dummy routines in the IPOPT library.
C
      external EV_HLV_DUMMY, EV_HOV_DUMMY, EV_HCV_DUMMY
C
C     The following is used in the estimate of the required work space
C     to guess the fillin in the linear solver.  This is not used within
C     IPOPT, unless dynamic memory allocation is not used.
C
      double precision fillinfact
C
C     Set algorithmic parameters:
C
      NARGS = 1
      ARGS(1) = 1.d-8
      CARGS(1) = 'dtol'
C
C     As a simple example, we pass the constants in the constraints to
C     the EVAL_C routine via the "private" DAT array.
C
      DAT(1) = 25.d0
      DAT(2) = 40.d0
C
C     Uncomment the follwing lines if you want some estimate of required
C     integer and double precision work space (this could be particularly
C     helpful if USE_MALLOC is not defined, i.e. does not dynamically
C     allocate memory on its own).
C
      fillinfact = 10.d0
      call ESTIMATE_WS(N, X, M, NLB, NUB, fillinfact, lrwe, liwe, IERR,
     1     EV_F, EV_C,EV_G, EV_A, EV_H, EV_HLV_DUMMY, EV_HOV_DUMMY,
     1     EV_HCV_DUMMY, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(*,*) 'Work space estimation return with error IERR = ',
     1        IERR
         stop
      endif
      write(*,*) 'estimated double precision work space requirement = ',
     1     lrwe
      write(*,*) 'estimated integer work space requirement          = ',
     1     liwe
C
C     Call optimization routine:
C
      f = IPOPT(N, X, M, NLB, ILB, BNDS_L, NUB, IUB, BNDS_U, V_L, V_U,
     1     LAM, C, LRW, RW, LIW, IW, ITER, IERR, EV_F, EV_C, EV_G, EV_A,
     1     EV_H, EV_HLV_DUMMY, EV_HOV_DUMMY, EV_HCV_DUMMY, DAT, IDAT,
     1     NARGS, ARGS, CARGS)
C
C     Output:
C
      if( IERR.eq.0 ) then
         write(*,*)
         write(*,*) 'The solution was found after ',ITER,' Iterations.'
         write(*,*)
         write(*,*) 'The final value of the objective function is ',f
         write(*,*)
         write(*,*) 'The optimal values of X are:'
         write(*,*)
         do i = 1, N
            write(*,*) 'X  (',i,') = ',X(i)
         enddo
         write(*,*)
         write(*,*) 'The multipliers for the lower bounds are:'
         write(*,*)
         do i = 1, NLB
            write(*,*) 'V_L(',ILB(i),') = ',V_L(i)
         enddo
         write(*,*)
         write(*,*) 'The multipliers for the upper bounds are:'
         write(*,*)
         do i = 1, NUB
            write(*,*) 'V_U(',IUB(i),') = ',V_U(i)
         enddo
         write(*,*)
         write(*,*) 'The multipliers for the equality constraints are:'
         write(*,*)
         do i = 1, M
            write(*,*) 'LAM(',i,') = ',LAM(i)
         enddo
         write(*,*)
      else
         write(*,*)
         write(*,*) 'An error occoured after ',ITER,' Iterations.'
         write(*,*) 'The error code is ',IERR
         write(*,*)
      endif
C
      end
C
C =============================================================================
C
C                    Computation of objective function
C
C =============================================================================
C
      subroutine EV_F(N, X, F, DAT, IDAT)
      implicit none
      integer N
      double precision F, X(N)
      double precision DAT(*)
      integer IDAT(*)
      F = X(1)*X(4)*(X(1)+X(2)+X(3)) + X(3)
      return
      end
C
C =============================================================================
C
C                Computation of gradient of objective function
C
C =============================================================================
C
      subroutine EV_G(N, X, G, DAT, IDAT)
      implicit none
      integer N
      double precision G(N), X(N)
      double precision DAT(*)
      integer IDAT(*)
      G(1) = X(4)*(2d0*X(1)+X(2)+X(3))
      G(2) = X(1)*X(4)
      G(3) = X(1)*X(4) + 1d0
      G(4) = X(1)*(X(1)+X(2)+X(3))
      G(5) = 0d0
      return
      end
C
C =============================================================================
C
C                     Computation of equality constraints
C
C =============================================================================
C
      subroutine EV_C(N, X, M, C, DAT, IDAT)
      implicit none
      integer N, M
      double precision C(M), X(N)
      double precision DAT(*)
      integer IDAT(*)
      C(1) = X(1)*X(2)*X(3)*X(4) - X(5) - DAT(1)
      C(2) = X(1)**2 + X(2)**2 + X(3)**2 + X(4)**2 - DAT(2)
      return
      end
C
C =============================================================================
C
C                Computation of Jacobian of equality constraints
C
C =============================================================================
C
      subroutine EV_A(TASK, N, X, NZ, A, ACON, AVAR, DAT, IDAT)
      integer TASK, N, NZ
      double precision X(N), A(NZ)
      integer ACON(NZ), AVAR(NZ), I
      double precision DAT(*)
      integer IDAT(*)
C
C     structure of Jacobian:
C
      integer AVAR1(9), ACON1(9)
      data  AVAR1 /1, 2, 3, 4, 5, 1, 2, 3, 4/
      data  ACON1 /1, 1, 1, 1, 1, 2, 2, 2, 2/
      save  AVAR1, ACON1
C
      if( TASK.eq.0 ) then
         NZ = 9
      else
        Do I = 1, 9
          AVAR(I) = AVAR1(I)
          ACON(I) = ACON1(I)
        EndDo
        A(1) = X(2)*X(3)*X(4)
        A(2) = X(1)*X(3)*X(4)
        A(3) = X(1)*X(2)*X(4)
        A(4) = X(1)*X(2)*X(3)
        A(5) = -1d0
        A(6) = 2d0*X(1)
        A(7) = 2d0*X(2)
        A(8) = 2d0*X(3)
        A(9) = 2d0*X(4)
      endif
      return
      end
C
C =============================================================================
C
C                Computation of Hessian of Lagrangian
C
C =============================================================================
C
      subroutine EV_H(TASK, N, X, M, LAM, NNZH, HESS, IRNH, ICNH,
     1     DAT, IDAT)
      implicit none
      integer TASK, N, M, NNZH, i
      double precision X(N), LAM(M), HESS(NNZH)
      integer IRNH(NNZH), ICNH(NNZH)
      double precision DAT(*)
      integer IDAT(*)
C
C     structure of Hessian:
C
      integer IRNH1(10), ICNH1(10)
      data  IRNH1 /1, 2, 2, 3, 3, 3, 4, 4, 4, 4/
      data  ICNH1 /1, 1, 2, 1, 2, 3, 1, 2, 3, 4/
      save  IRNH1, ICNH1

      if( TASK.eq.0 ) then
         NNZH = 10
      else
         do i = 1, 10
            IRNH(i) = IRNH1(i)
            ICNH(i) = ICNH1(i)
            HESS(i) = 0d0
         enddo
C
C     objective function
C
         HESS(1) = 2d0*X(4)
         HESS(2) = X(4)
         HESS(4) = X(4)
         HESS(7) = 2d0*X(1) + X(2) + X(3)
         HESS(8) = X(1)
         HESS(9) = X(1)
C
C     first constraint
C
         HESS(2) = HESS(2) + LAM(1) * X(3)*X(4)
         HESS(4) = HESS(4) + LAM(1) * X(2)*X(4)
         HESS(5) = HESS(5) + LAM(1) * X(1)*X(4)
         HESS(7) = HESS(7) + LAM(1) * X(2)*X(3)
         HESS(8) = HESS(8) + LAM(1) * X(1)*X(3)
         HESS(9) = HESS(9) + LAM(1) * X(1)*X(2)
C
C     second constraint
C
         HESS(1) = HESS(1) + LAM(2) * 2d0
         HESS(3) = HESS(3) + LAM(2) * 2d0
         HESS(6) = HESS(6) + LAM(2) * 2d0
         HESS(10)= HESS(10)+ LAM(2) * 2d0
      endif
      return
      end
