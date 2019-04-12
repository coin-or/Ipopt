C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     $Id: basics.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Some basic routines, e.g. obtain collocation points, compute values
C     of basis polynomials.
C
C ============================================================================
C
      SUBROUTINE INITD2N(NCOL)
C !DEC$ ATTRIBUTES DLLEXPORT :: INITD2N
C
C ============================================================================
C
C     Initialization of the COMMON block /COLCOM/ in DAE2NLP.INC.
C     This routine has to be called once at the beginning.
C
C     Author:  Andreas Waechter    03-25-02
C
C     Meaning of the variables in COLCOM:
C
C     RHO       values of the collocations (Radau points in [0,1])
C     COEF      coefficients of the polynomial basis Omega;
C                 it is COEF(NCOL-i+1,q) = Omega^{(i)}_q(0)
C                 where Omega_q(0) = 0 and Omega^{(i)}_q(rho_r) = delta_{q,r}
C     OMEGAQ    value of the polynomial Omega at the collocation points RHO
C     OMEGA1    valud of the polynomial Omega at 1.d0
C
      implicit none
      integer NCOL              ! number of collocation points
C
      include 'DAE2NLP.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DAENLP/
C
      character*80 line
      double precision dummy
      integer i
C
C     Set collocation points
C
      if( NCOL.eq.1 ) then
         RHO(1) = 1.d0
      elseif( NCOL.eq.2 ) then
         RHO(1) = 0.3333333333333333d+00
         RHO(2) = 0.1000000000000000d+01
      elseif( NCOL.eq.3 ) then
         RHO(1) = 0.1550510257216822d+00
         RHO(2) = 0.6449489742783178d+00
         RHO(3) = 0.1000000000000000d+01
      elseif( NCOL.eq.4 ) then
         RHO(1) = 0.8858795951270395d-01
         RHO(2) = 0.4094668644407347d+00
         RHO(3) = 0.7876594617608471d+00
         RHO(4) = 0.1000000000000000d+01
      elseif( NCOL.eq.5 ) then
         RHO(1) = 0.5710419611451768d-01
         RHO(2) = 0.2768430136381238d+00
         RHO(3) = 0.5835904323689168d+00
         RHO(4) = 0.8602401356562195d+00
         RHO(5) = 0.1000000000000000d+01
      elseif( NCOL.eq.6 ) then
         RHO(1) = 0.3980985705146874d-01
         RHO(2) = 0.1980134178736082d+00
         RHO(3) = 0.4379748102473861d+00
         RHO(4) = 0.6954642733536361d+00
         RHO(5) = 0.9014649142011736d+00
         RHO(6) = 0.1000000000000000d+01
      else
         write(line,100) NCOL
 100     format(
     1   'initd2n: Error:  Invalid number of collocation points:',i4)
         call C_OUT(2,0,1,line)
         stop
      endif
C
C     Compute coefficients for the basis polynomial Omega
C
      do i = 1, NCOL
         call DCOPY(NCOL, 0.d0, 0, COEF(1,i), 1)
         COEF(NCOL-i+1,i) = 1.d0
         call SOLVEVM(NCOL, RHO, COEF(1,i))
      enddo
C
C     Compute values of basis polynomials at collocations points and 1.d0
C
      call GETPOLY(0, 1.d0, NCOL, OMEGA1, dummy)
      do i = 1, NCOL
         call GETPOLY(0, RHO(i), NCOL, OMEGAQ(1,I), dummy)
      enddo
      return
      end

C ============================================================================
C
      subroutine GETPOLY(TASK, T, NCOL, OMEGA, PHI)
C
C ============================================================================
C
C     Get value of basis polynomials at given point
C
C     Author:  Andreas Waechter    03-25-02
C
      implicit none
C
      include 'DAE2NLP.INC'
C
      integer TASK              ! 0: only OMEGA
                                ! 1: both OMEGA and PHI
      double precision T        ! point within [0,1] at which polynamial is
                                ! to be computed
      integer NCOL              ! number of collocation points
      double precision OMEGA(NCOL)
              ! omega(q) = sum( COEF(ncol-j+1,q)/(j!) s^{j-1}, j=1,ncol )
              !          = Omega_q(s) / s
      double precision PHI(NCOL)
              ! theta(q) = sum( COEF(ncol-j+1,q)/((j-1)!) s^{j-1}, j=1,ncol )
              !          = Omega_q^{(1)}(s)
              !          = Theta_q (s)

      double precision tmp(NCOLD2NMAX), sum
      integer j, q

      do j = 1, NCOL
         tmp(j) = T/dble(NCOL+1-j)
      enddo
C
C     Compute OMEGAs
C
      do q = 1, NCOL
         sum = COEF(1,q)
         do j = 2, NCOL
            sum = sum*tmp(j-1) + COEF(j,q)
         enddo
         OMEGA(q)=sum
      enddo
C
C     Compute PHIs
C
      if( TASK.eq.1 ) then
         do q = 1, NCOL
            sum = COEF(1,q)
            do j = 2, NCOL
               sum = sum*tmp(j) + COEF(j,q)
            enddo
            PHI(q)=sum
         enddo
      endif

      return
      end

C ============================================================================
C
      subroutine SOLVEVM( NCOL, RHO, F )
C
C ============================================================================
C
C     Solving Vandermonde System for polynomial interpolation
C     (based on Algorithm 4.6.1 in Golub and van Loan)
C
C     Author:  Andreas Waechter    03-25-02
C
      implicit none
      integer NCOL               ! number of collocation points
      double precision RHO(NCOL) ! collocation points
      double precision F(NCOL)   ! In:  RHS  Out: Solution
C
C     Here we solve
C
C     sum(i=0..NCOL-1, F_out(NCOL+1-i) * RHO(j)**(i-1)/(i-1)!)
C      = F_in(NCOL+1-i)
C
      integer k, i
      do k = 1, NCOL-1
         do i = NCOL, k+1, -1
            F(NCOL+1-i) = (F(NCOL+1-i)-F(NCOL+2-i))/(RHO(i)-RHO(i-k))
         enddo
      enddo
      do k = NCOL-1, 1, -1
         do i = k, NCOL-1
            F(NCOL+1-i) = F(NCOL+1-i) - F(NCOL-i)*RHO(k)
         enddo
      enddo
C     correct the factorials
      i = 1
      do k = 2, NCOL-1
         i = i * k
         F(NCOL-k) = F(NCOL-k) * dble(i)
      enddo

      return
      end

C ============================================================================
C
      SUBROUTINE APPROX (TASK, NZ, NY, NU, NE, NCOL, COEF, W, U, TI,
     1     IELE, T, OMEGA, PHI, ZVAL, YVAL, UVAL)
C
C ============================================================================
C
C     A.Waechter, 9-15-00
C
C       evaluate Z (and maybe Y and U) on time T based interpolation in
C          element IELE; results in ZVAL, YVAL, UVAL
C
C   variables
C
C     TASK     1: compute Z at collocation point with basis-function values
C                 OMEGA
C              2: compute Z at T (OMEGA will be overwritten)
C              3: compute Z,Y,U at T (OMEGA, PHI will be overwritten)
C              4: compute dZ at T (OMEGA, PHI will be overwritten)
C                    (result in ZVAL)
C     COEF     from INITCOLCOM
C     TI       element boundaries
C     W        [Z at beinning of element, Zdot at first col point, Y at
C               first col point, ..., Zdot at last col point, Y at last col
C               point]
C     U        [U at first col point, .., U at last col point]
C     IELE     number of element in which T is
C
      implicit none
      integer TASK, NZ, NY, NU, NE, NCOL
      double precision COEF(*)
      double precision W(NZ+NCOL*(NZ+NY))
      double precision U(NCOL*NU)
      double precision TI(NE+1)
      integer IELE
      double precision T, OMEGA(NCOL), PHI(NCOL)
      double precision ZVAL(NZ)
      double precision YVAL(NY)
      double precision UVAL(NU)

      double precision s, del_t, sum
      integer task2, lw, lu, i, j
C
C     If necessary compute RK basis
C
      if( TASK.ne.1 ) then
         s = (T - TI(IELE))/(TI(IELE+1) - TI(IELE))
         if( TASK.eq.2 ) then
            task2 = 0
         else
            task2 = 1
         endif
         call GETPOLY(task2, s, NCOL, OMEGA, PHI)
      endif
C
C     Compute values of Z
C
      if( TASK.ne.4 ) then
         del_t = T - TI(IELE)
         do i = 1, NZ
            lw = NZ + i
            sum = 0.d0
            do j = 1, NCOL
               sum = sum + OMEGA(j)*W(lw)
               lw = lw + NZ+NY
            enddo
            ZVAL(i) = sum*del_t + W(i)
         enddo
      endif
C
C     Compute values of Y and U
C
      if( TASK.eq.3 ) then
         do i = 1, NY
            lw = NZ + NZ + i
            sum = 0.d0
            do j = 1, NCOL
               sum = sum + PHI(j)*W(lw)
               lw = lw + NZ+NY
            enddo
            YVAL(i) = sum
         enddo
         do i = 1, NU
            lu = i
            sum = 0.d0
            do j = 1, NCOL
               sum = sum + PHI(j)*U(lu)
               lu = lu + NU
            enddo
            UVAL(i) = sum
         enddo
      endif
C
C     Compute values of dZ
C
      if( TASK.eq.4 ) then
         do i = 1, NZ
            lw = NZ + i
            sum = 0.d0
            do j = 1, NCOL
               sum = sum + PHI(j)*W(lw)
               lw = lw + NZ+NY
            enddo
            ZVAL(i) = sum
         enddo
      endif

      return
      end
