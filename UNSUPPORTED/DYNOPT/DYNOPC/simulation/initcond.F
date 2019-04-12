C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: initcond.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine INITCOND(TOL, Z0, DZ, Y0, U0, P0, CNR_NLEQ1S, IERR)
C
C     Determine consistent inital conditions using NLEQ1S
C
C     Authors:  Andreas Waechter    09-28-01
C
      implicit none
C
      double precision TOL      ! Tolerance for NLEQ1S
      double precision Z0(*)    ! initial conditions for Z (not changed)
      double precision DZ(*)    ! derivates for Z's, changed to obtain
                                ! consistent initial conditions
      double precision Y0(*)    ! changed to obtain consistent initial
                                ! conditions
      double precision U0(*)    ! values of U's (not changed)
      double precision P0(*)    ! values of P's (not changed)
      integer CNR_NLEQ1S        ! Fortran Channel number for NLEQ1S output
      integer IERR

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/

C     Variables for NLEQ1S
      integer n
      integer NFMAX
C      parameter( NFMAX = (NZMAX+NYMAX)**2 ) ! That's worst (=dense) case
      parameter( NFMAX = 100000 ) ! That's worst (=dense) case
      external fcn_nleq1s, jac_nleq1s
      double precision x(NZMAX+NYMAX), xscal(NZMAX+NYMAX)
      double precision rtol
      integer iopt(50)
      integer liwk, li2wk, lrwk
      parameter( liwk  = 8*(NZMAX+NYMAX) + 57 )
C      parameter( li2wk = 8*NFMAX + 5*(NZMAX+NYMAX) )
      parameter( li2wk = 6*NFMAX + 5*(NZMAX+NYMAX) ) ! (we use iopt(37)=1)
      parameter( lrwk  = 3*NFMAX + 11*(NZMAX+NYMAX) + 62 )
      integer iwk(liwk), i2wk(li2wk)
      double precision rwk(lrwk)

      integer i
C
C     First we need to set all control variables that do not have a specified
C     profile as fixed for the DAEMODEL interface - same for parameters
C     (see DYNAUX.INC)
C
      if( NU_PROF.ne.NU ) then
         do i = 1, NU
            if( IU_PROF(i).eq.0 ) then
               IU_PROF(i) = 1
               call DCOPY(NELE+1, U0(i), 0, U_PROF(i,1), NUMAX)
            endif
         enddo
         NU_PROF = NU
      endif

      if( NP_FIX.ne.NP ) then
         do i = 1, NP
            if( IP_FIX(i).eq.0 ) then
               IP_FIX(i) = 1
               P_FIX(i) = P0(i)
            endif
         enddo
         NP_FIX = NP
      endif
C
C     Copy user provided initial conditions for Z into DYNAUX common block
C
      call DCOPY(NZ, Z0, 1, ZINIT, 1)
C
C     Prepare parameters for NLEQ1S
C
      n = NZ + NY

C     Use default for most IOPT parameters
      do i = 1, 50
         iopt(i) = 0
         iwk (i) = 0
      enddo
      call DCOPY(50, 0.d0, 0, rwk, 1)
      if( CNR_NLEQ1S.gt.0 ) then
C         iopt(11) = 1           ! error and some warning messages
C         iopt(12) = CNR_NLEQ1S
C         iopt(13) = 1           ! standard iteration output
C         iopt(14) = CNR_NLEQ1S
      endif
      iopt(37) = 1              ! sparsity pattern does not change between
                                ! Newton iterations
      iopt(39) = 1              ! force strict convergence criterion

      rtol = TOL

C ============================================================================
C
C     Obtain consistent initial conditions from NLEQ1S
C
C ============================================================================

C     Prepare initial point for NLEQ1S
      call DCOPY(NZ, DZ  , 1, x(   1), 1)
      call DCOPY(NY, Y0  , 1, x(NZ+1), 1)

C     Don't know what should be put here...?????
      call DCOPY(n, 0.d0, 0, xscal, 1) 

C     Reset Flag for JAC_NLEQ1S
      JAC_INIT = .FALSE.

      write(*,90)
 90   format(/
     1     '>>>>> Starting computation of consistent initial conditions',
     2     /6x,'with NLEQ1S .'/)

      call NLEQ1S(n, NFMAX, FCN_NLEQ1S, JAC_NLEQ1S, x, xscal, rtol,
     1     iopt, IERR, liwk, iwk, li2wk, i2wk, lrwk, rwk)

      if( IERR.eq.0 .or. IERR.eq.5 .or. IERR.eq.4 ) then
         IERR = 0
         write(*,100)
 100     format(//
     1        '<<<<< NLEQ1S successfully determined ',
     2        'consistent initial conditions'/)
      else
         write(*,101) IERR, iwk(23)
 101     format(/
     1     6x,'***********************************************'/
     2     6x,'*     Problem with during initialization!     *'/
     3     6x,'*        NLEQ1S returned IERR = ',i4,'          *'/
     3     6x,'*            with IWK(23) = ',i4,'              *'/
     4     6x,'* Still need to decide what to do.  Now abort *'/
     4     6x,'***********************************************'/)
         stop
      endif

C     Copy initial conditions back into DZ and Y0
      call DCOPY(NZ, x(   1), 1, DZ, 1)
      call DCOPY(NY, x(NZ+1), 1, Y0, 1)


      return
      end

C ============================================================================

C $Id: initcond.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine FCN_NLEQ1S(N, X, F, IFAIL)
C
C     Compute residual for NLEQ1S using DAEMODEL_F.
C
C     Author:     Andreas Waechter     09-28-01
C
      implicit none
      include 'DYNAUX.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/
      integer N
      double precision X(N)
      double precision F(N)
      integer IFAIL
C
      double precision dummy(1)

      call DAEMODEL_F( NZ, NY, 0, 0, TI(1), ZINIT(1),
     1     X(1), X(NZ+1), dummy, dummy, F )
#ifndef WIN32
      call CPRINTF('%s','.')
#else
      write(*,100) '.'
  100 format(a1\)
#endif
      return
      end

C ============================================================================

C $Id: initcond.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine JAC_NLEQ1S(N, X, DFDX, IROW, ICOL, NFILL, IFAIL)
C
C     Compute Jacobian for NLEQ1S using DAEMODEL_DF.
C
C     Author:     Andreas Waechter     09-28-01
C
      implicit none
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/
      integer N
      double precision X(N)
      integer NFILL
      double precision DFDX(NFILL)
      integer IROW(NFILL), ICOL(NFILL)
      integer IFAIL

C     Store sparsity pattern
      integer NNZ, NNZMAX
C      parameter( NNZMAX = (NZMAX+NYMAX)**2 ) ! worst (=dense) case
      parameter( NNZMAX = 10*(NZMAX+NYMAX) ) ! maixmal 10 nonzeros per row
      integer IROWSTORE(NNZMAX)
      integer ICOLSTORE(NNZMAX)
      save NNZ, IROWSTORE, ICOLSTORE

C     work space for Model Jacobian
      integer LDDF
      parameter( LDDF = NZMAX + NYMAX )
      double precision DF(LDDF, 2*NZMAX + NYMAX)

      integer i, j
      double precision dummy

C     If this is first call, determine sparsity pattern and store it

      if( .not.JAC_INIT ) then
         JAC_INIT = .true.
         call DAEMODEL_DF_STRUC(NZ, NY, 0, 0, DF, LDDF)
C
         NNZ = 0
         do i = 1, NZ+NY
            do j = 1, NZ+NY
               if( DF(i,j+NZ).ne.0.d0 ) then
                  if( NNZ.eq.NNZMAX ) then
                     IFAIL = -2 ! NNZMAX too small
                     return
                  endif
                  NNZ = NNZ + 1
                  IROWSTORE(NNZ) = i
                  ICOLSTORE(NNZ) = j
               endif
            enddo
         enddo
         if( NNZ.gt.NFILL ) then
            IFAIL = -3          ! NFILL too small
         endif
      endif

      NFILL = NNZ

C     Get value of Jacobian
      call DAEMODEL_DF( NZ, NY, 0, 0, TI(1), ZINIT, X(1),
     1     X(NZ+1), dummy, dummy, DF, LDDF )

C     ... and copy nonzeros into what NLEQ1S wants
      do i = 1, NFILL
         IROW(i) = IROWSTORE(i)
         ICOL(i) = ICOLSTORE(i)
         DFDX(i) = DF(IROW(i),ICOL(i)+NZ)
      enddo
#ifndef WIN32
      call CPRINTF('%s','*')
#else
      write(*,100)"*"
  100 format(a1\)
#endif
      return
      end
