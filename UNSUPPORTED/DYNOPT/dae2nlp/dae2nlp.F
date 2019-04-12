C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine DAE2NLP(TASK, NZ, NY, NU, NP, NCOL, NAC, NELE, TI,
     1     ZINIT, X, IFLAG, IVAR, NNZA, A, VIN, VOUT, LAM,
     2     LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
C
C     $Id: dae2nlp.F 531 2004-03-11 01:31:07Z andreasw $
C
C     MAIN SUBROUTINE FOR DAE DECOMPOSITION
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
C     Version 1.1    ( 02-26-00 )  (implicit formulation of DAE,
C                                   additional constraints)
C     Version 2.0    ( 04-29-02 )  (released as part of IPOPT 2.0)
C
C     Description:
C     ============
C
C     This package provides an implementation of the discretization of a
C     system of Dynamic Algebraic Equations (DAE) using orthogonal
C     collocation.
C
C     Given an index one DAE system of the form
C
C     F(t,z(t),dz/dt(t),y(t),u(t),p) = 0,       z(t_0) = ZINIT     (DAE)
C
C     with differential variables z(t) of dimension NZ,
C          algebraic    variables y(t) of dimension NY,
C          control      variables u(t) of dimension NU,
C     time-independent parameters p    of dimension NP, and
C     initial conditions ZINIT,
C
C     these routines provide a discretization in NELE given finite elements
C     with element boundary TI(1), TI(2), ..., TI(NELE+1).  Within each
C     finite element, the differential variables are approximated by
C     polynomials of degree NCOL, and the algebraic and control variables are
C     approximated by polynomials of degree NCOL-1, where NCOL is the number
C     of collocation points.
C
C     Using this approximation, the system (DAE) is discretized as follows.
C     (For brevity, K = NCOL)
C
C     0 = Z(TI(1)) - ZINIT
C
C     0 = F(T(l,1),Z(T(l,1)),Zdot(T(l,1)),Y(T(l,1)),U(T(l,1)),P)   \
C                           .                                      |
C                           .                                      |  for l
C                           .                                      |    =
C     0 = F(T(l,1),Z(T(l,K)),Zdot(T(l,K)),Y(T(l,K)),U(T(l,K)),P)   | 1,..,NELE
C     0 = ADDCON(Z(T(l,1)),Zdot(T(l,1)),Y(T(l,1)),U(T(l,1)),...    |
C                Z(T(l,K)),Zdot(T(l,K)),Y(T(l,K)),U(T(l,K)))       |
C           continuity of Z across element boundaries              /
C
C     Here, ADDCON are additional linear constraints that can be posed
C     on the coefficients in each element l, for example in order to
C     reduce the order of the polynomial approximating the control
C     variables.  We assume, that the number of such additional
C     constraints is NAC per element.  With this distribution comes a
C     file addcon.f, which as an example implements constraints that
C     reduce the degree of the polynomial approximating the constraints.
C     If you want to use it, you first need initialize some COMMON block
C     by calling ADDCON_INIT.
C
C     More information can be found in "A reduced space interior point
C     strategy for optimization of differential algebraic systems" by
C     Arturo M. Cervantes, Andreas Waechter, Reha H. Tutuncu, Lorenz
C     T. Biegler in Computers and Chemical Engineering 24 (2000) 29-51.
C     However, since the publishers seemed not to be able to correct the
C     typoes that we marked you might be better of downloading this
C     paper from my web site (if it still exists) at
C     http://dynopt.cheme.cmu.edu/andreasw .
C     
C     It follows a description of the subroutines describing the model
C     and computing its derivatives that have to be provided by the
C     user.
C
C     - Computation of F:
C
C       DAEMODEL_F(NZ, NY, NU, NP, T, Z, ZDOT, Y, U, P, F)
C
C       Input:  T          time at which model is to be evaluated
C               Z(NZ)      values of differential variables Z
C               ZDOT(NZ)   values of Zdot
C               Y(NY)      values of algebraic variables Y
C               U(NU)      values of control variables U
C               P(NP)      values of parameters P
C       Output: F(NZ+NY)   violation of the DAE equations
C
C     - Sparsity structure of Jacobian DF
C
C       DAEMODEL_DF_STRUC(NZ, NY, NU, NP, DF, LDDF)
C
C       Output: LDDF       leading dimension of DF (at least NZ+NY)
C               DF(LDDF,2*NZ+NY+NU+NP)
C                          element nonzero, if according entry in Jacobian
C                          is structurally nonzero.
C                          (first dimension ordered like F;
C                           second dimension in order Z, ZDOT, Y, U, P)
C                          NOTE: DF is initialized to all 0.d0 on input
C
C     - Computation of Jacobian DF of F:
C
C       DAEMODEL_DF(NZ, NY, NU, NP, T, Z, ZDOT, Y, U, P, DF, LDDF)
C
C       Input:  see DAEMODEL_F
C       Output: LDDF       leading dimension of DF (at least NZ+NY)
C               DF(LDDF,2*NZ+NY+NU+NP)
C                          Jacobian of F
C                          (first dimension ordered like F;
C                           second dimension in order Z, ZDOT, Y, U, P)
C                          NOTE: DF is not initialized, but only non-zero
C                             elements (as identified before by
C                             DAEMODEL_DF_STRUC) have to be set in DAEMODEL_F
C
C     - Hessian-vector products
C
C       DAEMODEL_HV(NRHS, NZ, NY, NU, NP, T, Z, ZDOT, Y, U, P,
C                LAM, RHS, RES, LDRS)
C
C       Input:  T, Z, ZDOT, Y, P 
C                          see DAEMODEL_F
C               NRHS       number of vectors to be multiplied
C               LAM(NZ+NY) weights for individual constraint Hessians
C               LDRS       leading dimension of RHS and RES
C                            (at least 2*NZ+NY+NU+NP)
C               RHS(LDRS,NRHS)
C                          NRHS right hand sides (vectors which has to be
C                          postmuliplied with the Hessian)
C       Output: RES(LDRS,NRHS)
C                          RHS Hessian-vector products corresponding to RHS
C               RES(i,.) = sum( j=1..NZ+NY, LAM(j)*Hessian_F(j) )*RHS(i,.)
C
C
C     Below, the vector X of coefficient for the polynomial
C     representation is organized as follows:
C
C     X = [ Z   (TI(1)) ,
C
C           Zdot(T(l,1)), Y(T(l,1)), ... , Zdot(T(l,K)), Y(T(l,K)), Z(TI(l+1)),
C                                   (for each l = 1, NELE)
C
C
C           U(T(l,1)), U(T(l,2)), ... , U(T(l,K)),
C                                   (for each l = 1, NELE)
C
C           P ]
C
C     and thus is of dimension
C
C     NX = NZ + NELE*(NCOL*(NZ+NY) + NZ) + NCOL*NU + NP
C
C     Usage:
C     ======
C
C     1. Call with TASK = 0 to obtain storage and work space requirements.
C           This also initializes internal data structures.
C           (storage space means, that the values in those vectors must not
C            be changed between individual calls, whereas work space vectors
C            can be used somewhere else)
C
C     2. Call with TASK = 1 to initialize storage space and compute permutation
C           vector IVAR.  Even if no permutation is desired, this call is
C           necessary in order to initialize the structure of the Jacobian.
C
C     Further operations (for TASK > 1) can then be done in no particular order
C     except for the following restrictions:
C
C     -  If operations with the Jacobian (TASK > 3) are to be done, the
C        nonzeros of the Jacobian have to be computed beforehand with TASK = 3.
C        Note, that it is possible to store the nonzeros corresponding to
C        different points outside of this routine and perform TASK = 4,..,7
C        for those nonzeros in any order.
C
C     -  Before TASK = 9,10 can be executed, this routines must have been
C        called with TASK = 8 in order to factorize the elemental basis
C        matrices.
C        Note, that it is necesssary for TASK = 9,10 to provide the same
C        values of the nonzeros A as they were used for the last call with
C        TASK = 8.
C
C     -  It is possible to change the basis partition at a later point by
C        calling with TASk = 1 again.  However, since then internal data
C        structures change, nonzero elements and factorizations have to
C        be computed again.
C
C
C     Parameters:
C     ===========
C
C     TASK      (Integer, scalar). 
C               IN always includes NZ, NY, NU, NP, NCOL, NAC, NELE, IW, RW
C                  (IS, RS for TASK > 0)
C               OUT always includes IW, RW, (IS, RS for TASK > 0)
C
C               0: Initialize data structures; compute sizes of the vectors
C                  A, IS, RS, IW, RW. 
C                  Work space requirements are in that particular case
C                  LRW = (NZ+NY)*(2*NZ+NY+NU+NP) + NNZAC + NCOL*(NZ+NY+NU)
C                  LIW = NZ + NCOL*(NZ+NY+NU) + NP + 2*NNZAW ,
C                  where NNZAC is the number of nonzeros for the additional
C                  constraints (can be determined by calling ADDCON_NNZ).
C                  The user is asked to provide an estimate of the expected
C                  fill-in during the factorization of the elemental basis
C                  matrices in VIN(1) (e.g. VIN(1) = 5d0).
C                  The user also decides here whether a fast and less memory
C                  intensive factorization using the same pivot sequence in
C                  each element is desired.  If IFLAG = 0, fast factorizaton
C                  will be chosen, otherwise not.
C                     NOTE: THE FAST OPTION SEEMS NOT TO BE WORKING!!!
C                  IN : VIN(1), IFLAG
C                  OUT: NNZA, LRS, LIS, LRW, LIW
C
C               1: Determine basis partition vector IVAR and initialize data
C                  structures for Jacobian.
C                  If 1 <= IFLAG <= NELE, the elemental Jacobian for finite
C                  element IFLAG is evaluated and MA28AD is used to obtain
C                  a basis partition of this Jacobian.
C                  If IFLAG <= 0, the user provided elemental partition in
C                  the first NCOL*(NZ+NY+NU) entries of IVAR is taken.  Here,
C                  the columns corresponding to 
C                  [Zdot(T(l,1)), Y(T(l,1)), ... , Zdot(T(l,K)), Y(T(l,K),
C                   U(T(l,1)), ... , U(T(l,K))]
C                  are considered, i.e. for i = 1,..,NCOL*(NZ+NY)+NAC the column
C                  IVAR(i) in the above ordering is basic, and the remaining
C                  NCOL*NU-NAC columns nonbasic.
C                  For IFLAG > 0:
C                  IN : TI, X, IFLAG
C                  OUT: IVAR
C                  For IFLAG <= 0:
C                  IN : IVAR
C                  OUT: IVAR
C
C               2: Evaluate right hand side of collocation equations at X and
C                  provide result in VOUT.
C                  IN : ZINIT, TI, X
C                  OUT: VOUT
C
C               3: Evaluate nonzeros of Jacobian in store those in A.
C                  IN : TI, X
C                  OUT: A
C
C               4: Compute VOUT = C * VIN     (C: basis matrix)
C                  using nonzeros provided in A.
C                  IN : A, TI, VIN
C                  OUT: VOUT
C
C               5: Compute VOUT = C^T * VIN   (C: basis matrix)
C                  using nonzeros provided in A.
C                  IN : A, TI, VIN
C                  OUT: VOUT
C
C               6: Compute VOUT = N * VIN     (N: nonbasis matrix)
C                  using nonzeros provided in A.
C                  IN : A, TI, VIN
C                  OUT: VOUT
C
C               7: Compute VOUT = N^T * VIN   (N: nonbasis matrix)
C                  using nonzeros provided in A.
C                  IN : A, TI, VIN
C                  OUT: VOUT
C
C               8: Factorize basis matrix C with nonzeros provided in A.
C                  (Pivot sequence based on element IFLAG)
C                  IN : A, IFLAG
C                  OUT: -
C
C               9: Compute VOUT = C^{-1} * VIN
C                  using last factorzation (note: A needs to be the same as
C                  when last factorization was done.)
C                  IN : A, VIN
C                  OUT: VOUT
C
C              10: Compute VOUT = C^{-T} * VIN
C                  using last factorzation (note: A needs to be the same as
C                  when last factorization was done.)
C                  IN : A, VIN
C                  OUT: VOUT
C
C              11: Compute VOUT = Hessian LAG(X,LAM) * VIN
C                  IN : X, TI, VIN, LAM
C                  OUT: VOUT
C
C     NZ       (Integer, scalar) Number of differential variables
C
C     NY       (Integer, scalar) Number of algebraic variables
C
C     NU       (Integer, scalar) Number of control variables
C
C     NU       (Integer, scalar) Number of time-independent parameters
C
C     NCOL     (Integer, scalar) Number of collocation point per element
C
C     NAC      (Integer, scalar) Number of additional constraints
C
C     NELE     (Integer, scalar) Number of finite elements
C
C     TI       (Double precison, dimension NELE+1)
C                  Boundaries of finite elements (in increasing order)
C
C     ZINIT    (Double precison, dimension NZ)
C                  Initial conditions of differential variables
C
C     X        (Double precison, dimension NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
C                  Vector of coefficients for polynomical approximations
C                  (see above) ordered as described in "Description" above.
C
C     IFLAG    (Integer, scalar)
C                  For TASK = 0: Indicates whether fast factorization reusing
C                      same pivot sequence within each element is to be used
C                      (see TASK).
C                      DON'T USE IFLAG = 0!!!!
C                  For TASK = 1,8: Indicates at which finite element the
C                      partition or factorization is to be based on.
C
C     IVAR     (Integer, dimension NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
C                  Information about variable partition:
C                  X(i) is the value of XORIG(IVAR(i)) where XORIG is the
C                  unpartitioned vector of coefficients as described in
C                  "Description".  (See TASK = 1 for more information).
C
C     NNZA     (Integer, scalar)
C                  Dimension of A (Output only for TASK = 0)
C
C     A        (Double precision, dimension NNZA)
C                  Nonzero elements of Jacobian
C
C     VIN      (Double precision)
C                  For TASK = 0: (dimension 1), user-provided estimate of
C                      fill-in factor
C                  For TASK = 4,5,7,9,10:
C                      (dimension NZ+NELE*(NCOL*(NZ+NY)+NZ+NAC))
C                  For TASK = 6: (dimension NELE*NCOL*NU - NAC + NP)
C
C     VOUT     (Double precision)
C                  For TASK = 2,4,5,6,9,10:
C                      (dimension NZ+NELE*(NCOL*(NZ+NY)+NZ+NAC))
C                  For TASK = 7: (dimension NELE*NCOL*)NU - NAC + NP)
C
C     LAM      (Double precision, dimension NZ+NELE*(NCOL*(NZ+NY)+NZ+NAC) )
C                  For TASK = 11: vector of Lagrangian multipliers for
C                      Hessian of Lagrangian (ordered like right hand side)
C
C     LRS      (Integer, scalar)
C                  Length of RS (Output only for TASK = 0)
C
C     RS       (Double precision, dimension LRS)
C                  Storage space, must not be changed between calls
C
C     LIS      (Integer, scalar)
C                  Length of IS (Output only for TASK = 0)
C
C     IS       (Integer, dimension LIS)
C                  Storage space, must not be changed between calls
C
C     LRW      (Integer, scalar)
C                  Length of RW (Output only for TASK = 0)
C                  (For TASK = 0: at least
C                      (NZ+NY)*(2*NZ+NY+NU+NP) + NNZAC + NCOL*(NZ+NY+NU))
C
C     RW       (Double precision, dimension LRW)
C                  Work space
C
C     LIW      (Integer, scalar)
C                  Length of IW (Output only for TASK = 0)
C                  (For TASK = 0: at least NZ + NCOL*(NZ+NY+NU) + NP + 2*NNZAW)
C
C     IW       (Integer, dimension LIW)
C                  Work space
C
C     IERR     (Integer, scalar) always output
C                  Error flag:
C                      0: everything went fine (hopefully...)
C                      6: unable to find partition (elemental Jacobian seems
C                         rank deficient)
C                      7: IFLAG has invalid value
C                     20: factorization aborted, since one of the elemental
C                         basis matrices is singular
C                     21: The estimate of the fill-in factor (see TASK = 0)
C                         turns out to be too small
C                     22: Iterative refinement for MA48CD did not converge
C                     28: Trouble in MA28 (not good...)
C                     48: Trouble in MA48 (not good...)
C                     98: LRW too small
C                     99: LIW too small
C
      implicit none
C
C     Declaration of parameters
C
      integer TASK, NZ, NY, NU, NP, NCOL, NAC, NELE
      double precision TI(NELE+1)
      double precision ZINIT(NZ)
      double precision X(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      integer IFLAG
      integer IVAR(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      integer NNZA
      double precision A(*)    !(NNZA)
      double precision VIN(*)
      double precision VOUT(*)
      double precision LAM(NZ + NELE*(NCOL*(NZ+NY)+NZ+NAC))
      integer LRS
      double precision RS(LRS)
      integer LIS
      integer IS(LIS)
      integer LRW
c      double precision RW(LRW)
      double precision RW(*)
      integer LIW
c      integer IW(LIW)
      integer IW(*)
      integer IERR
C
C     Things to keep in mind:
C
      integer NNZAELE, NNZAC, LA48, P_IVARC, P_IR, P_JC, P_IFACT
      save    NNZAELE, NNZAC, LA48, P_IVARC, P_IR, P_JC, P_IFACT
      integer P_RFACT, P_IRC, P_IRN, P_JCC, P_JCN, P_C, P_N
      save    P_RFACT, P_IRC, P_IRN, P_JCC, P_JCN, P_C, P_N
      logical HAVE_PIV, FAST
      save    HAVE_PIV, FAST
      integer NNZFZZ, NNZFW, NNZFWZ(7), NNZFWY(7), NNZFWU(7), NNZCW
      save    NNZFZZ, NNZFW, NNZFWZ,    NNZFWY,    NNZFWU,    NNZCW
      integer NNZCWW(7), NNZFU, NNZFUZ(7), NNZFUY(7), NNZFUU(7)
      save    NNZCWW,    NNZFU, NNZFUZ,    NNZFUY,    NNZFUU
      integer NNZFPP, NNZCU, NNZCUU(7)
      save    NNZFPP, NNZCU, NNZCUU
      integer NNZFWD(7), NNZAW, NNZFUD(7), NNZAU
      save    NNZFWD   , NNZAW, NNZFUD   , NNZAU
C
C     local variables
C
      integer llrw, lliw, p_rwend, p_iwend, p_rw, p_iw, lw, lu, k, i
      integer l, p_ivarc1, lddf, p_df, p_zval, ttask, p_rhs
      integer nnzc, nkeep, p_res, p_ptmp, p_aac, p_irac, p_jcac
      integer p_wdummy, p_udummy, lw2, lu2
      character*1 trans

COLDADOLC
      integer p_zeros, p_hvtmp

      IERR = 0
      p_rwend = 0
      p_iwend = 0

      goto (10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100)
     1     TASK + 1

      IERR = 1
      goto 9999
C
C =============================================================================
C     BEGIN: Compute estimate for storage and work space demand; initialize
C =============================================================================
C
 10   continue
C
C     Estimate number of nonzeros in one element
C
C     First for the additional constraints
C
      if( NAC.gt.0 ) then
         call ADDCON_NNZ(NZ, NY, NU, NNZAC)
      else
         NNZAC = 0
      endif

      lddf     = NZ+NY
      p_df     = p_rwend
      p_aac    = p_df + lddf*(2*NZ+NY+NU+NP)
      p_wdummy = p_aac + NNZAC
      p_udummy = p_wdummy + NCOL*(NZ+NY)
      p_rwend  = p_udummy + NCOL*NU
      if( p_rwend.gt.LRW ) goto 9098

      p_iw    = p_iwend
      p_irac  = p_iw + NZ + NCOL*(NZ+NY+NU) + NP
      p_jcac  = p_irac + NNZAC
      p_iwend = p_jcac + NNZAC
      if( p_iwend.gt.LIW ) goto 9099

      call ESTIMNNZ(NZ, NY, NU, NP, NCOL, NAC, lddf, RW(p_df+1),
     1     NNZAELE, nnzc, NNZAC, IW(p_irac+1), IW(p_jcac+1),
     2     RW(p_aac+1), RW(p_wdummy+1), RW(p_udummy+1), IW(p_iw+1))
      NNZA = NELE * NNZAELE
C
C     Initialize our COMMON block :)
C
      call INITD2N(NCOL)
C
C     Initialize internal stuff and estimate storage space
C
      if( IFLAG.eq.0 ) then
         FAST = .true.
      else
         FAST = .false.
      endif

      HAVE_PIV = .false.
      if( vin(1).lt.2d0 ) then
         IERR = 7
         goto 9999
      endif
      LA48     = int(vin(1)*dble(nnzc))

      P_IVARC  = 0
COLD      P_IR     = P_IVARC + NZ + NELE*(NCOL*(NZ+NY+NU)+NZ)
      P_IR     = P_IVARC + NCOL*(NZ+NY+NU)
      P_JC     = P_IR    + NNZAELE + NCOL*NZ
      P_IFACT  = P_JC    + NNZAELE + NCOL*NZ
      nkeep    = 7 + 10*NCOL*(NZ+NY+NAC)
      if( FAST ) then
         LIS = P_IFACT + NNZAELE + LA48 + nkeep
      else
         LIS = P_IFACT + NNZAELE + NELE*(LA48 + nkeep)
      endif

      P_RFACT  = 0
      LRS      = P_RFACT + NELE*LA48
C
C     Estimate work spaces (task = 1 is probably the worst, so take that)
C
      LIW = max(40*NNZAELE + 13*NCOL*(NZ+NY+NU), 
     1        40*NNZAELE + 10*NCOL + 2*NNZAC)
      LRW = max(20*NNZAELE + NCOL*(NZ+NY+NU),
     1        NNZAELE + (NZ+NY)*(2*NZ+NY+NU+NP) + NZ + NNZAC +
     2        NCOL*(NZ+NY+NU))
C
      goto 9999
C
C =============================================================================
C     END: Compute estimate for storage and work space demand; initialize
C =============================================================================
C

C
C =============================================================================
C     BEGIN: Determine partition of variables; initialize structural data
C =============================================================================
C
 100  continue
C
C     if IFLAG <= 0, take user provided partition, otherwise call MA28AD
C
      if( IFLAG.le.0 ) then
         do i = 1, NCOL*(NZ+NY+NU)
            IS(P_IVARC+i) = IVAR(i)
         enddo
      elseif( IFLAG.le.NELE ) then
C
C     call MA28 to obtain partition for one individual element (IVARC)
C
         llrw    = max(20*NNZAELE + NCOL*(NZ+NY+NU),
     1        NNZAELE + (NZ+NY)*(2*NZ+NY+NU+NP) + NZ + NNZAC +
     2        NCOL*(NZ+NY+NU))
         p_rw    = p_rwend
         p_rwend = p_rw + llrw
         if( p_rwend.gt.LRW ) goto 9098

         lliw    = max(40*NNZAELE + 13*NCOL*(NZ+NY+NU), 
     1        40*NNZAELE + 10*NCOL + 2*NNZAC)
         p_iw    = p_iwend
         p_iwend = p_iw + lliw
         if( p_iwend.gt.LIW ) goto 9099

         call PARTELE(NZ, NY, NU, NP, NCOL, NAC, NELE, NNZAELE, NNZAC,
     1        TI, X, IFLAG, IS(P_IVARC+1), llrw, RW(p_rw+1), lliw,
     2        IW(p_iw+1), IERR)
         if( IERR.ne.0 ) goto 9999
         p_rwend = p_rw
         p_iwend = p_iw
      else
         IERR = 7
         goto 9999
      endif
C
C     Now that we have IVARC, compute overall IVAR
C
      do i = 1, NZ
         IVAR(i) = i
      enddo
      lw  = NZ
      lu  = NZ + NELE*(NCOL*(NZ+NY)+NAC+NZ)
      lw2 = NZ
      lu2 = NZ + NELE*(NCOL*(NZ+NY)+NZ)
      do k = 1, NELE
         do i = 1, NCOL*(NZ+NY)+NAC
            l = IS(P_IVARC+i)
            if( l.le.NCOL*(NZ+NY) ) then
               IVAR(lw+i) = lw2 + l
            else
               IVAR(lw+i) = lu2 + l - NCOL*(NZ+NY)
            endif
         enddo
         do i = 1, NCOL*NU - NAC
            l = IS(P_IVARC+NCOL*(NZ+NY)+NAC+i)
            if( l.le.NCOL*(NZ+NY) ) then
               IVAR(lu+i) = lw2 + l
            else
               IVAR(lu+i) = lu2 + l - NCOL*(NZ+NY)
            endif
         enddo
         lw  = lw  + NCOL*(NZ+NY) + NAC
         lu  = lu  + NCOL*NU - NAC
         lw2 = lw2 + NCOL*(NZ+NY)
         lu2 = lu2 + NCOL*NU

         do i = 1, NZ
            IVAR(lw+i) = lw2+i
         enddo
         lw  =  lw + NZ
         lw2 = lw2 + NZ
      enddo
      do i = 1, NP
         IVAR(lu+i) = lu2+i
      enddo
C
C     Now initialize structural data
C
      p_ivarc1 = p_iwend
      p_irac   = p_ivarc1 + NCOL*(NZ+NY+NU)
      p_jcac   = p_irac + NNZAC
      p_iwend  = p_jcac + NNZAC
      if( p_iwend.gt.LIW ) goto 9099

      lddf     = NZ+NY
      p_df     = p_rwend
      p_aac    = p_df + lddf*(2*NZ+NY+NU+NP)
      p_wdummy = p_aac + NNZAC
      p_udummy = NCOL*(NZ+NY)
      p_rwend  = p_udummy + NCOL*NU
      if( p_rwend.gt.LRW ) goto 9098

      do i = 1, NCOL*(NZ+NY+NU)
         IW(p_ivarc1+IS(P_IVARC+i)) = i ! inverse of IVARC
      enddo

      call STRUCTURE(NZ, NY, NU, NP, NCOL, NAC, IW(p_ivarc1+1), NNZAELE,
     1     IS(P_IR+1), IS(P_JC+1), lddf, RW(p_df+1), NNZAC,
     2     IW(p_irac+1), IW(p_jcac+1), RW(p_aac+1),
     3     NNZFZZ, NNZFW, NNZFWZ, NNZFWD, NNZFWY, NNZFWU, NNZAW, NNZCW,
     4     NNZCWW, NNZFU, NNZFUZ, NNZFUD, NNZFUY, NNZFUU, NNZFPP, NNZAU,
     5     NNZCU, NNZCUU, RW(p_wdummy+1), RW(p_udummy+1))
      P_IRC = P_IR
      P_IRN = P_IRC + NCOL*NNZFZZ + NNZFW + NNZAW + NNZCW
      P_JCC = P_JC
      P_JCN = P_JCC + NCOL*NNZFZZ + NNZFW + NNZAW + NNZCW
      P_C   = 0
      P_N   = P_C   + NELE*(NCOL*NNZFZZ + NNZFW + NNZAW)
C
C     Reset factorizaton flag
C
      HAVE_PIV = .false.
C
      goto 9999
C
C =============================================================================
C     END: Determine partition of variables; initialize structural data
C =============================================================================
C

C
C =============================================================================
C     BEGIN: Compute right hand side of system
C =============================================================================
C
 200  continue

      p_zval  = p_rwend
      p_rwend = p_zval  + NZ
      if( p_rwend.gt.LRW ) goto 9098

      call COLRHS(NZ, NY, NU, NP, NCOL, NAC, NELE, TI, ZINIT, X,
     1     VOUT, RW(p_zval+1))
C
      goto 9999
C
C =============================================================================
C     END: Compute right hand side of system
C =============================================================================
C

C
C =============================================================================
C     BEGIN: Evaluate Jacobian at point X
C =============================================================================
C
 300  continue
C
      lddf    = NZ+NY
      p_zval  = p_rwend
      p_df    = p_zval  + NZ
      p_aac   = p_df    + lddf*(2*NZ+NY+NU+NP)
      p_rwend = p_aac   + NNZAC
      if( p_rwend.gt.LRW ) goto 9098

      if( NAC.gt.0 ) then
         p_ivarc1 = p_iwend
         p_iwend  = p_ivarc1 + NCOL*(NZ+NY+NU)
      else
         p_ivarc1 = 0
      endif
      p_irac  = p_iwend
      p_jcac  = p_irac + NNZAC
      p_iwend = p_jcac + NNZAC
      if( p_iwend.gt.LIW ) goto 9099
      if( NAC.gt.0 ) then
         do i = 1, NCOL*(NZ+NY+NU)
            IW(p_ivarc1+IS(P_IVARC+i)) = i ! inverse of IVARC
         enddo
      endif

      call FILLJAC(NZ, NY, NU, NP, NCOL, NAC, NELE, TI, X,
     1     IS(P_IVARC+1), IW(p_ivarc1+1), NNZFZZ, NNZFW, NNZFWZ, NNZFWD,
     2     NNZFWY, NNZFWU, NNZFU, NNZFUZ, NNZFUD, NNZFUY, NNZFUU,
     3     NNZFPP, NNZAW, NNZAU, IS(P_IRC+1), IS(P_JCC+1),
     3     IS(P_IRN+1), IS(P_JCN+1), A(P_C+1), A(P_N+1),
     3     RW(p_zval+1), lddf, RW(p_df+1),
     4     NNZAC, IW(p_irac+1), IW(p_jcac+1), RW(p_aac+1))
C
      goto 9999
C
C =============================================================================
C     END: Evaluate Jacobian at point X
C =============================================================================
C

C
C =============================================================================
C     BEGIN: Compute products C * vin or C^T * vin
C =============================================================================
C
 400  continue
      trans = 'N'
      goto 550
C
 500  trans = 'T'
C
 550  continue
      
      call MULT_C(trans, NZ, NY, NU, NP, NCOL, NAC, NELE, TI,
     1     NNZFZZ, NNZFW, NNZAW, NNZCW, NNZCWW, IS(P_IRC+1),
     2     IS(P_JCC+1), A(P_C+1), VIN, VOUT)
C
      goto 9999
C
C =============================================================================
C     END: Compute products C * vin or C^T * vin
C =============================================================================
C

C
C =============================================================================
C     BEGIN: Compute products N * vin or N^T * vin
C =============================================================================
C
 600  continue
      trans = 'N'
      goto 750
C
 700  trans = 'T'
C
 750  continue
C
      call MULT_N(trans, NZ, NY, NU, NP, NCOL, NAC, NELE, TI,
     1     NNZFU, NNZAU, NNZFPP, NNZCU, NNZCUU, IS(P_IRN+1),
     2     IS(P_JCN+1), A(P_N+1), VIN, VOUT)
C
      goto 9999
C
C =============================================================================
C     END: Compute products N * vin or N^T * vin
C =============================================================================
C

C
C =============================================================================
C     BEGIN: Factorize basis matrix
C =============================================================================
C
 800  continue
C
      if( IFLAG.lt.1 .or. IFLAG.gt.NELE ) then
         IERR = 7
         goto 9999
      endif
C
C     Decide whether this is first factorization of that matrix
C
      if( HAVE_PIV ) then
         ttask = 2
      else
         ttask = 1
         HAVE_PIV = .true.
      endif
C
C     Do the factorization
C
      llrw    = LA48            ! worst estimate in case we initialize MA48
      p_rw    = p_rwend
      p_rwend = p_rw + llrw
      if( p_rwend.gt.LRW ) goto 9098

      lliw    = 9*(NCOL*(NZ+NY)+NAC) + LA48
                                ! worst estimate in case we initialize MA48
      p_iw    = p_iwend
      p_iwend = p_iw + lliw
      if( p_iwend.gt.LIW ) goto 9099

      call FACT_C(ttask, NZ, NY, NCOL, NAC, NELE, NNZFZZ, NNZFW, NNZAW,
     1     A(P_C+1), IS(P_IRC+1), IS(P_JCC+1), IFLAG, LA48, FAST,
     2     IS(P_IFACT+1), RS(P_RFACT+1), llrw, RW(p_rw+1),
     2     lliw, IW(p_iw+1), IERR)
C
      goto 9999
C
C =============================================================================
C     END: Factorize basis matrix
C =============================================================================
C

C
C =============================================================================
C     BEGIN: Compute products C^{-1} * vin or C^{-T} * vin
C =============================================================================
C
 900  continue
      trans = 'N'
      goto 1050
C
 1000 trans = 'T'
C
 1050 continue
C
      llrw    = 5*(NCOL*(NZ+NY)+NAC)
      p_rhs   = p_rwend
      p_rw    = p_rhs + NCOL*(NZ+NY)+NAC
      p_rwend = p_rw  + llrw
      if( p_rwend.gt.LRW ) goto 9098

      lliw    = NCOL*(NZ+NY)+NAC
      p_iw    = p_iwend
      p_iwend = p_iw  + lliw
      if( p_iwend.gt.LIW ) goto 9099

      call SOLVE_C(trans, NZ, NY, NU, NP, NCOL, NAC, NELE, TI,
     1     NNZFZZ, NNZFW, NNZAW, NNZCW, NNZCWW, IS(P_IRC+1),
     2     IS(P_JCC+1), A(P_C+1), LA48, FAST, IS(P_IFACT+1),
     3     RS(P_RFACT+1), VIN, VOUT, RW(p_rhs+1), llrw, RW(p_rw+1),
     4     lliw, IW(p_iw+1), IERR)
C
      goto 9999
C
C =============================================================================
C     END: Compute products C^{-1} * vin or C^{-T} * vin
C =============================================================================
C

C
C =============================================================================
C     BEGIN: Evaluate Hessian of Lagrangian times vector product
C =============================================================================
C
 1100 continue

#ifdef OLDADOLC
C
C     Leave in for possible debugging.
C
      p_zval  = p_rwend
      p_zeros = p_zval  + NZ
      p_hvtmp = p_zeros + max(NZ,NY,NU,NP)
      p_ptmp  = p_hvtmp + 2*NZ+NY+NU+NP
      p_rwend = p_ptmp  + NZ
      if( p_rwend.gt.LRW ) goto 9098

      call LAGVEC_OLD(NZ, NY, NU, NP, NCOL, NAC, NELE, TI, X, LAM, VIN,
     1     VOUT, RW(p_zval+1), RW(p_zeros+1), RW(p_hvtmp+1),
     2     RW(p_ptmp+1))
      goto 9999
#else

      p_zval  = p_rwend
      p_rhs   = p_zval  + NZ
      p_res   = p_rhs + (NCOL+1)*(2*NZ+NY+NU+NP)
      p_ptmp  = p_res + (NCOL+1)*(2*NZ+NY+NU+NP)
      p_rwend = p_ptmp  + NZ
      call LAGVEC(NZ, NY, NU, NP, NCOL, NAC, NELE, TI, X, LAM, VIN,
     1     VOUT, RW(p_zval+1), RW(p_rhs+1), RW(p_res+1),
     2     RW(p_ptmp+1))
      goto 9999
#endif
C
C =============================================================================
C     END: Evaluate Hessian of Lagrangian times vector product
C =============================================================================
C

C
C     Error Handling
C
 9098 IERR = 98
      goto 9999
 9099 IERR = 99
      goto 9999
C
 9999 continue
      return
      end
