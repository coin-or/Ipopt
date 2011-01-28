C Copyright (C) 2002, 2004, 2005 Carnegie Mellon University,
C                                Dominique Orban and others.
C All Rights Reserved.
C This code is published under the Eclipse Public License.
C*******************************************************************************
      PROGRAM           IPOPTMA
C
C     IPOPT CUTEr driver.
C     D. Orban,  adapted from Andreas Waechter's CUTE driver.
C     Adapted for C++ version by Andreas Waechter, Oct 2004
C
      IMPLICIT NONE
      INTEGER IOUT
      PARAMETER( IOUT = 6 )

C
C     Maximal sizes for CUTEr
C
CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )

C
C     
C
      INTEGER N, M
      DOUBLE PRECISION X( CUTE_NMAX )
      DOUBLE PRECISION X_L( CUTE_NMAX )
      DOUBLE PRECISION X_U( CUTE_NMAX )
      DOUBLE PRECISION Z_L( CUTE_NMAX )
      DOUBLE PRECISION Z_U( CUTE_NMAX )
      DOUBLE PRECISION G( CUTE_MMAX )
      DOUBLE PRECISION G_L( CUTE_MMAX )
      DOUBLE PRECISION G_U( CUTE_MMAX )
      DOUBLE PRECISION LAM( CUTE_MMAX )

      INTEGER IERR, IPSOLVE

C32BIT      INTEGER IPROBLEM, IPCREATE
C64BIT      INTEGER*8 IPROBLEM, IPCREATE
C
      integer IDX_STYLE, NELE_JAC, NELE_HESS

      external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS
C
C     The following arrays are work space for the evaluation subroutines
C
      DOUBLE PRECISION DAT(CUTE_NMAX+CUTE_NZMAX)
      INTEGER IDAT(2*CUTE_NZMAX)

      REAL CALLS( 7 ), CPU( 2 )
      CHARACTER*10 PNAME
      CHARACTER*10 VNAMES( CUTE_NMAX ), GNAMES( CUTE_MMAX )
      DOUBLE PRECISION F
C
      logical equatn(CUTE_MMAX), linear(CUTE_MMAX)
      integer i, cnr_input
      logical efirst, lfirst, nvfrst, ex
      double precision init_val
C
C     Initialize the CUTEr interface and get the initial point
C
      cnr_input = 60
      efirst = .false.
      lfirst = .false.
      nvfrst = .false.

      open(cnr_input,file='OUTSDIF.d',status='old')

      call CSETUP(cnr_input, IOUT, N, M, X, X_L, X_U, CUTE_NMAX,
     1     equatn, linear, LAM, G_L, G_U, CUTE_MMAX,
     2     efirst, lfirst, nvfrst)
      close(cnr_input)
C
C     See if we want to set a different initial point
C
      inquire(file='INITPOINT.VAL', exist=ex)
      if (ex) then
         open(70, file='INITPOINT.VAL', status='old')
         read(70,'(d25.16)') init_val
         do i = 1, N
            X(i) = init_val
         enddo
         close(70)
      endif
C
C     Obtain the number of nonzeros in Jacobian and Hessian
C
      CALL CDIMSJ(NELE_JAC)
      NELE_JAC = NELE_JAC - N
      if (NELE_JAC.gt.CUTE_NZMAX) then
         write(*,*) 'NELE_JAC = ',NELE_JAC,' larger than CUTE_NZMAX = ',
     1        CUTE_NZMAX
         write(*,*) 'Increase CUTE_NZMAX'
         stop
      endif
      CALL CDIMSH(NELE_HESS)
      if (NELE_HESS.gt.CUTE_NZMAX) then
         write(*,*) 'NELE_HESS = ',NELE_HESS,
     1        ' larger than CUTE_NZMAX = ', CUTE_NZMAX
         write(*,*) 'Increase CUTE_NZMAX'
         stop
      endif
C
C     Get problem name.
C
      CALL CNAMES(N, M, PNAME, VNAMES, GNAMES)
C
C     Call IPOPT
C
      IDX_STYLE = 1
      IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,
     1     IDX_STYLE, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)
      if (IPROBLEM.eq.0) then
         write(*,*) 'Error creating Ipopt Problem.'
         stop
      endif
C
      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)
C
      call IPFREE(IPROBLEM)
C
C     Display CUTEr statistics
C
      CALL CREPRT( CALLS, CPU )
      WRITE ( IOUT, 2000 ) PNAME, N, M, CALLS(1), CALLS(2),
     .     CALLS(3), CALLS(4), CALLS(5), CALLS(6), CALLS(7),
     .     IERR, F, CPU(1), CPU(2)
c
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *     ,/,' Code used               :  IPOPT',    /
     *     ,' Problem                 :  ', A10,    /
     *     ,' # variables             =      ', I10 /
     *     ,' # constraints           =      ', I10 /
     *     ,' # objective functions   =      ', E15.7 /
     *     ,' # objective gradients   =      ', E15.7 /
     *     ,' # objective Hessians    =      ', E15.7 /
     *     ,' # Hessian-vector prdct  =      ', E15.7 /
     *     ,' # constraints functions =      ', E15.7 /
     *     ,' # constraints gradients =      ', E15.7 /
     *     ,' # constraints Hessians  =      ', E15.7 /
     *     ,' Exit code               =      ', I10 /
     *     ,' Final f                 = ', E15.7 /
     *     ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ,' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     ,/,66('*') / )

 9999 CONTINUE
      END



C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Eclipse Public License.
C*******************************************************************************
C
      subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
C
C*******************************************************************************
C
C    $Id$
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute objective function value to CUTEr problem
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
CA    Andreas Waechter      02/25/99
CA    Andreas Waechter      10/29/04 adapted for C++ version
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
CP   N         I    INT    number of variables in problem statement
CP   X         I    DP     point where F is to be evaluated
CP   NEW_X     I    INT    if 1, X has not been changed since last call
CP   F         O    DP     objective function value
CP   IDAT      P    INT    privat INT data for evaluation routines
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IERR      O    INT    set to nonzero value if error occurred
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
CS    COFG
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
      integer N
      double precision X(N)
      integer NEW_X
      double precision F
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision dummy
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
C
C     Call COFG to obtain value of objective function
C
      call COFG( N, X, F, dummy, .false.)

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Eclipse Public License.
C*******************************************************************************
C
      subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
C
C*******************************************************************************
C
C    $Id$
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute gradient of objective function to CUTEr problem
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
CA    Andreas Waechter      02/25/99
CA    Andreas Waechter      10/29/04 adapted for C++ version
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
CP   N         I    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where G is to be evaluated
CP   NEW_X     I    INT    if 1, X has not been changed since last call
CP   GRAD      O    DP     gradient of objective function
CP   IDAT      P    INT    privat INT data for evaluation routines
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IERR      O    INT    set to nonzero value if error occurred
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
CS    COFG
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
      integer N
      double precision X(N)
      integer NEW_X
      double precision GRAD(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision f
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
C
C     Call COFG to obtain gradient of objective function
C
      call COFG( N, X, f, GRAD, .true.)

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Eclipse Public License.
C*******************************************************************************
C
      subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
C
C*******************************************************************************
C
C    $Id$
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute values of constraints to CUTEr problem
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
CA    Andreas Waechter      02/25/99
CA    Andreas Waechter      07/01/99 BUG: problems if ineq not first
CA    Andreas Waechter      10/29/04 adapted for C++ version
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
CP   N         I    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where G is to be evaluated
CP   NEW_X     I    INT    if 1, X has not been changed since last call
CP   M         I    INT    number of constraints
CP   G         O    DP     values of constraints
CP   IDAT      P    INT    privat INT data for evaluation routines
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IERR      O    INT    set to nonzero value if error occurred
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
CS    CCFG
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
      integer N
      double precision X(N)
      integer NEW_X
      integer M
      double precision G(M)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision dummy
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
C
C     Call CCFG to obtain constraint values, but without slacks
C
      call CCFG(N, M, X, M, G, .FALSE., 1, 1, dummy, .FALSE.)

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Eclipse Public License.
C*******************************************************************************
C
      subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A,
     1     IDAT, DAT, IERR)
C
C*******************************************************************************
C
C    $Id$
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute Jacobian of constraints to CUTEr problem
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
CA    Andreas Waechter      02/25/99
CA    Andreas Waechter      10/29/04 adapted for C++ version
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
CP   TASK      I    INT     =0: Fill ACON and AVAR, don't use A
CP                         <>0: Fill A, don't use ACON, AVAR
CP   N         I    INT    number of variables in problem statement
CP   X         I    DP     point where A is to be evaluated
CP   NEW_X     I    INT    if 1, X has not been changed since last call
CP   M         I    INT    number of constraints
CP   NZ        I    INT    number of nonzero elements
CP                                     (size of A, AVAR, ACON)
CP   ACON      O    INT    (only TASK=0) row indices
CP   AVAR      O    INT    (only TASK=0) column indices
CP   A         O    DP     (only TASK<>0) values in Jacobian
CP   IDAT      P    INT    privat INT data for evaluation routines
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IERR      O    INT    set to nonzero value if error occurred
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
CS    CDIMSJ
CS    CCFSG
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
      double precision X(N)
      integer NEW_X
      integer M
      integer NZ
      double precision A(NZ)
      integer ACON(NZ)
      integer AVAR(NZ)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i, nele_jac
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
      if( TASK.eq.0 ) then
C
C     Get the nonzero structure
C
         do i = 1, N
            DAT(i) = 0.d0
         enddo
         call CCFSG(N, M, DAT(1), M, DAT(1), nele_jac,
     1        NZ, DAT(N+1), AVAR, ACON, .TRUE.)
      else
C
C     Get the values of nonzeros
C
         call CCFSG(N, M, X, M, DAT(1), nele_jac,
     1        NZ, A, IDAT(1), IDAT(1+NZ), .TRUE.)
      endif

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Eclipse Public License.
C*******************************************************************************
C

      subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM,
     1     NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
C
C*******************************************************************************
C
C    $Id$
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute Hessian of Lagrangian for CUTEr problem
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
CA    Andreas Waechter      03/23/00
CA    Andreas Waechter      10/29/04 adapted for C++ version
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
CP   TASK      I    INT     =0: Fill IRNH and ICNH, don't use HESS
CP                         <>0: Fill HESS, don't use IRNH, ICNH
CP   N         I    INT    number of variables in problem statement
CP   X         I    DP     point where A is to be evaluated
CP   NEW_X     I    INT    if 1, X has not been changed since last call
CP   OBJFACT   I    DP     weighting factor for objective function Hessian
CP   M         I    INT    number of constriants
CP   LAM       I    DP     weighting factors for the constraints
CP   NEW_LAM   I    INT    if 1, LAM has not been changed since last call
CP   NNZH      I    INT    number of nonzero elements
CP                                     (size of HESS, IRNH, ICNH)
CP   IRNH      O    INT    (only TASK=0) row indices
CP   ICNH      O    INT    (only TASK=0) column indices
CP   HESS      O    DP     (only TASK<>0) values in Hessian
CP   IDAT      P    INT    privat INT data for evaluation routines
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IERR      O    INT    set to nonzero value if error occurred
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
CS    CSH
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
      double precision X(N)
      integer NEW_X
      double precision OBJFACT
      integer M
      double precision LAM(M)
      integer NEW_LAM
      integer NNZH
      integer IRNH(NNZH)
      integer ICNH(NNZH)
      double precision HESS(NNZH)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i, nnzh2
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
      if( TASK.eq.0 ) then
C
C     Get the nonzero structure
C
         do i = 1, N
            DAT(i) = 0.d0
         enddo
         call CSH(N, M, DAT(1), M, DAT(1), nnzh2, NNZH, DAT(N+1),
     1        IRNH, ICNH)
      else
C
C     Call CSH to get the values
C
         if( OBJFACT.ne.0.d0 ) then

            if( OBJFACT.ne.1.d0 ) then
               do i = 1, M
                  DAT(i) = LAM(i)/OBJFACT
               enddo
               call CSH(N, M, X, M, DAT(1), nnzh2, NNZH, HESS,
     1              IDAT(1), IDAT(1+NNZH))
               do i = 1, NNZH
                  HESS(i) = HESS(i)*OBJFACT
               enddo
            else
               call CSH(N, M, X, M, LAM, nnzh2, NNZH, HESS,
     1              IDAT(1), IDAT(1+NNZH))
            endif

         else
C     now we have to call CSH twice, since we can't otherwise get rid of
C     the objective function entries
            do i = 1, M
               DAT(i) = 0.d0
            enddo
            call CSH(N, M, X, M, DAT(1), nnzh2, NNZH, DAT(1+M),
     1           IDAT(1), IDAT(1+NNZH))
            call CSH(N, M, X, M, LAM, nnzh2, NNZH, HESS,
     1           IDAT(1), IDAT(1+NNZH))
            do i = 1, NNZH
               HESS(i) = HESS(i) - DAT(M+i)
            enddo
         endif
      endif

 9999 continue
      return
      end
