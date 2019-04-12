C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine UPDATE_B(N, NIND, M, X, ITER, ERR, PZ, IVAR, NFIX,
     1     IFIX, NORIG, XORIG, CSCALE, RG, RGOLD, PYNRM, WBAR, ALPHA,
     2     G, GOLD, LAM, SKIP_UPDATE, B, W, C_SKIP,
     3     NLB, ILB, NUB, IUB, S_L, S_U, SIGMA_L, SIGMA_U,
     4     KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR, EV_F,
     5     EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: update_b.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Update reduced Hessian estimate
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
CP   N         I    INT    number of (free) variables; first NIND variables
CP                         are independent; remaining dependent
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of equality constraints / dependent variables
CP   X         I    DP     actual primal iterate
CP   ITER      I    INT    iteration counter (if ITER = 0 B is initialized)
CP   ERR       I    DP     actual KKT error
CP   PZ        I    DP     null space step from last iteration
CP                            (only for ITER > 0 )
CP                            (entries belong to independent
CP                             variables, these correspond to 's' in BFGS
CP                             formula, if multiplied by ALPHA)
CP   IVAR      I    INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP   NORIG     I    INT    total number of all variables (incl. fixed vars)
CP   XORIG     I    DP     actual iterate
CP                            XORIG is ordered in ORIGINAL order (i.e. not
CP                            partitioned into independent and dependent
CP                            variables)
CP   CSCALE    I    DP     scaling factors for constraints
CP   RG        I    DP     reduced gradient of objective function
CP   RGOLD     I    DP     reduced gradient of objective function from
CP                            last iteration
CP   PYNRM     I    DP     norm of old PY (for skipping rule)
CP   WBAR      I    DP     correction term (not for now!)
CP   ALPHA     I    DP     line search alpha from last iteration
CP   G         I    DP     actual gradient of objective function
CP                            (only for QQUASI < 0)
CP   GOLD      I    DP     gradient of objective function from last iteration
CP                            (only for QQUASI < 0)
CP   LAM       I    DP     Lagrangian multipliers
CP                            (only for QQUASI < 0)
CP   SKIP_UPDATE I/O DP    I: if true, do not do Quasi-Newton update
CP                         O: set to false
CP   B        I/O   DP     estimate for reduced Hessian
CP                            (packed upper triganular format!)
CP                            (for now only BFGS for orig Lagranian)
CP   W        I/O   DP     I: reduced primal-dual Hessian for barrier term (BB)
CP                         O: sum of B and BB
CP   C_SKIP    O    C*1    ='*' if update is skipped, otherwise = ' '
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. S_L(i) is slack for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. S_U(i) is slack for X(IUB(i)) )
CP   S_L       I    DP     slack variables for lower bounds
CP   S_U       I    DP     slack variables for upper bounds
CP   SIGMA_L   I    DP     primal-dual Hessian of lower bound barrier term
CP                            (NLB diagonal elements only)
CP   SIGMA_U   I    DP     primal-dual Hessian of upper bound barrier term
CP                            (NUB diagonal elements only)
CP   KCONSTR   I    INT    KCONSTR(1): LRS for CONSTR
CP                         KCONSTR(2): P_LRS for CONSTR
CP                         KCONSTR(3): LIS for CONSTR
CP                         KCONSTR(4): P_LIS for CONSTR
CP                         KCONSTR(5): LRW for CONSTR
CP                         KCONSTR(6): LIW for CONSTR
CP   LRS       I    INT    total length of RS
CP   RS       I/O   DP     DP storage space (all!)
CP   LIS       I    INT    total length of IS
CP   IS       I/O   INT    INT storage space (all!)
CP   LRW      I/O   INT    length of RW
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW      I/O   INT    length of IW
CP   IW       I/O   INT    can be used as INT work space but content will be
CP                            changed between calls
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
CP   EV_F      I    EXT    Subroutine for objective function
CP   EV_C      I    EXT    Subroutine for constraints
CP   EV_G      I    EXT    Subroutine for gradient of objective fuction
CP   EV_A      I    EXT    Subroutine for Jacobian
CP   EV_H      I    EXT    Subroutine for Lagrangian Hessian
CP   EV_HLV    I    EXT    Subroutine for Lagrangian Hessian-vector products
CP   EV_HOV    I    EXT    Subroutine for objective Hessian-vector products
CP   EV_HCV    I    EXT    Subroutine for constraint Hessian-vector products
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
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
CS    DDOT
CS    DCOPY
CS    DSCAL
CS    DAXPY
CS    DNRM2
CS    DSPMV
CS    DSPR
CS    GET_G
CS    COMPUTE_Y
CS    CONSTR
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
C     |QQUASI| = 1,3,5: BFGS update for original red. Hessian
C                2,4  : SR1 update for original red. Hessian
C      QQUASI > 0: Use multiplier free update
C             < 0: Use multipliers
C
      include 'IPOPT.INC'
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      integer NIND
      integer M
      double precision X(N)
      integer ITER
      double precision ERR
      double precision PZ(NIND)
      integer IVAR(N)
      integer NFIX
      integer IFIX(NFIX)
      integer NORIG
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      double precision RG(NIND)
      double precision RGOLD(NIND)
      double precision PYNRM
      double precision WBAR(NIND)
      double precision ALPHA
      double precision G(N)
      double precision GOLD(N)
      double precision LAM(M)
      logical SKIP_UPDATE
      double precision B(NIND*(NIND+1)/2)
      double precision W(NIND*(NIND+1)/2)
      character*1 C_SKIP
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      integer KCONSTR(6)
      integer LRS
      double precision RS(LRS)
      integer LIS
      integer IS(LIS)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR
      external EV_F
      external EV_C
      external EV_G
      external EV_A
      external EV_H
      external EV_HLV
      external EV_HOV
      external EV_HCV
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision DDOT, DNRM2
      integer i, k, p_iwend, p_rwend, p_y, p_bs, p_bdiag, p_gdiff
      integer p_pert, idummy, p_xorig, p_vout
      double precision dot, snorm, ynorm, pznrm, theta, bmin, tmp, tmp2
      double precision sBs, dummy, fact
      logical skip_rule, ex

      character*80 line(2)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
C
C     If no degrees of freedom, there is nothing to be done
C
      if( NIND.eq.0 ) then
         C_SKIP = 's'
         goto 9999
      endif

      p_iwend = 0
      p_rwend = 0

      C_SKIP = ' '

      if( ITER.eq.0 ) then
C
C     in first iteration initialize B
C
         if( QCG.ne.0 ) then
            goto 9020           ! unless this is done in GET_PZ_CG
         endif
C
C     If QBWARMSTART is set and file 'BWARM.DAT' exists, read initial
C     quasi-Newton estimate from file
C
         if( QBWARMSTART.eq.1 ) then
            inquire(file='BWARM.DAT', exist=ex)
            if( ex ) then
               open(80,file='BWARM.DAT',status='old',err=9008)
               read(80,1001,end=9009,err=9009) i
 1001          format(i16)
               if( i.ne.((NIND+1)*NIND)/2 ) then
                  write(line,100)
 100              format('update_b: Warning: Number of elements in',
     1                 ' BWARM.DAT for warm start',/,
     2      '          does not match.  Ignoring warm start option.')
                  call C_OUT(2,0,2,line)
                  close(80)
               else
                  do i = 1, ((NIND+1)*NIND)/2
                     read(80,1002,end=9009,err=9009) B(I)
 1002                format(d23.16)
                  enddo
                  close(80)
                  goto 9000
               endif
            else
               write(line,101)
 101           format('update_b:  File BWARM.DAT does not exist;',
     1              ' performing usual initialization.')
               call C_OUT(2,0,1,line)
            endif
         endif
C
C     This is the new version based on perturbation of diagonal...
C
         p_bdiag = p_rwend
         p_rwend = p_bdiag + NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif

         if( QINITB.ge.2 .and. QINITB.le.4 ) then

            p_gdiff = p_bdiag + NIND
            p_pert  = p_gdiff + N
            p_rwend = p_pert + N
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
C
C     Compute perturbation vector  (Z_0 * e )
C
            call DCOPY(NIND, 1d0, 0, RW(p_pert+M+1), 1)
            if( M.gt.0 ) then
               call CONSTR(5, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1              NORIG, XORIG, CSCALE, RW(p_pert+M+1), RW(p_pert+1),
     2              idummy, idummy,
     3              KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4              IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5              LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C,
     6              EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
               if( IERR.lt.0 ) then
                  write(line,*)
     1                 'update_b: Warning in CONSTR-5, IERR = ', IERR
                  call C_OUT(2,0,1,line)
               elseif( IERR.ne.0 ) then
                  write(line,*)
     1                 'update_b: Error in CONSTR-5, IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  goto 9999
               endif
               call DSCAL(M, -1d0, RW(p_pert+1), 1)
            endif
C
CTODO Cut length (need better theta!)
C
            theta = 1d-4
            call DSCAL(N, theta, RW(p_pert+1), 1)
            call DAXPY(N, 1d0, X, 1, RW(p_pert+1), 1)
C
C     Get gradient at perturbed point
C
            p_xorig = p_rwend
            p_rwend = p_xorig + NORIG
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
CTODO Check if xorigtmp is reallty necessary
            call DCOPY(NORIG, XORIG, 1, RW(p_xorig+1), 1)
            call GET_G(N, RW(p_pert+1), IVAR, NORIG, RW(p_xorig+1), M,
     1           CSCALE, 0, idummy, dummy, 0, idummy, dummy, 0.d0,
     1           RW(p_gdiff+1), LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     1           IW(p_iwend+1), IERR, EV_G, DAT, IDAT)
            if( IERR.ne.0 ) then
               write(line,*) 'update_b: get_g returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            p_rwend = p_pert
C
C     Compute diagonal elements of initial B
C
            call DAXPY(N, -1d0, G, 1, RW(p_gdiff+1), 1)
            if( M.gt.0 ) then
               call CONSTR(4, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1              NORIG, XORIG, CSCALE, RW(p_gdiff+1), RW(p_bdiag+1),
     2              idummy, idummy,
     3              KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4              IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5              LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C,
     6              EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
               if( IERR.lt.0 ) then
                  write(line,*)
     1                 'update_b: Warning in CONSTR-4, IERR = ', IERR
                  call C_OUT(2,0,1,line)
               elseif( IERR.ne.0 ) then
                  write(line,*)
     1                 'update_b: Error in CONSTR-4, IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  goto 9999
               endif
               call DSCAL(NIND, -1d0, RW(p_bdiag+1), 1)
               call DAXPY(NIND, 1d0, RW(p_gdiff+M+1), 1,
     1              RW(p_bdiag+1), 1)
            else
               call DCOPY(NIND, RW(p_gdiff+1), 1, RW(p_bdiag+1), 1)
            endif
            call DSCAL(NIND, 1d0/theta, RW(p_bdiag+1), 1)
            p_rwend = p_gdiff

         endif
C
C     Correct diagonal entries
C

CTODO WARM START
C
C     Implemented choices for QINITB
C
C     0: B_0 = identity
C     1: B_0 = 1d-6 * identity
C     2: B_0 = min of all abs(bdiag) * identity
C     3: B_0 = abs(bdiag)
C     4: B_0 = bdiag (! not for BFGS)
C     5: B_0 = identity, but re-scale before update (see p.200/201 in Nocedal)
C

CTODO Find good choice...(?)
         if( QINITB.eq.0 .or. QINITB.eq.5) then
            call DCOPY(NIND, 1d0, 0, RW(p_bdiag+1), 1)
         elseif( QINITB.eq.1 ) then
            call DCOPY(NIND, 1d-6, 0, RW(p_bdiag+1), 1)
         elseif( QINITB.eq.2 ) then
            bmin = dabs(RW(p_bdiag+1))
            do i = 2, NIND
               bmin = dmin1(bmin,dabs(RW(p_bdiag+i)))
            enddo
            bmin = dmax1(1d-6,bmin)
            call DCOPY(NIND, bmin, 0, RW(p_bdiag+1), 1)
         elseif( QINITB.eq.3 ) then
            do i = 1, NIND
               RW(p_bdiag+i) = dmax1(1d-6,dabs(RW(p_bdiag+i)))
            enddo
         elseif( QINITB.eq.4 ) then
            if( abs(QQUASI).ne.2 .and. abs(QQUASI).ne.4 ) then
               write(line,*)
     1              'update_b: Invalid choice of QINITB = ',QINITB
               call C_OUT(2,0,1,line)
               IERR = 4
               goto 9999
            endif
         endif
C
C     Initialize B
C
         call DCOPY( (NIND*(NIND+1))/2, 0.d0, 0, B, 1)
         k = 0
         do i = 1,NIND
            k = k + i
            B(k) = RW(p_bdiag+i)
         enddo
         p_rwend = p_bdiag

         goto 9000
      endif
C
C     Check if we want to switch from SR1 to BFGS
C
      if( abs(QQUASI).eq.3 .and. ERR.lt.QSR1TOL ) then
         call C_OUT(2,0,1,'Switching from BFGS to SR1!')
         if( QQUASI.eq.3 ) then
            QQUASI = 2
         else
            QQUASI = -2
         endif
      endif

      if( .not.SKIP_UPDATE ) then
C
C     BFGS
C
         if( abs(QQUASI).eq.1 .or. abs(QQUASI).eq.3 .or.
     1        abs(QQUASI).eq.5 ) then
C
C           compute y for BFGS formula
C
            p_y     = p_rwend
            p_rwend = p_y + NIND
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            call COMPUTE_Y(N, NIND, M, ITER, IVAR, NFIX, IFIX, NORIG,
     1           XORIG, CSCALE, RG, RGOLD, WBAR, ALPHA, G, GOLD,
     2           LAM, RW(p_y+1), KCONSTR,
     3           LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.gt.0 ) then
               write(line,*)
     1              'update_b: Error: compute_y ends with IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            elseif( IERR.ne.0 ) then
               write(line,*)
     1              'update_b: Warning: compute_y ends with IERR = ',
     2              IERR
               call C_OUT(2,0,1,line)
            endif
C
C     Add alpha* Z^T Sigma Z PZ in case overall reduced Hessian is approxiamted
C
            if( abs(QCG).eq.1 .and. (NLB.gt.0 .or. NUB.gt.0) ) then
               p_vout  = p_rwend
               p_rwend = p_vout + NIND
               if( p_rwend.gt.LRW ) then
                  IERR = 98
                  goto 9999
               endif
               call GET_ZWZV(4, N, NIND, M, ITER, IVAR, NFIX, IFIX,
     1              NORIG, XORIG, X, CSCALE, NLB, ILB, NUB, IUB,
     2              S_L, S_U, SIGMA_L, SIGMA_U, dummy, PZ, RW(p_vout+1),
     1              dummy, KCONSTR, LRS, RS, LIS, IS, LRW-p_rwend,
     1              RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR,
     2              EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     1              EV_HCV, DAT, IDAT)
               if( IERR.gt.0 ) then
                  write(line,*)
     1                 'update_b: Error: get_zwzv ends with IERR = ',
     1                 IERR
                  call C_OUT(2,0,1,line)
                  goto 9999
               elseif( IERR.ne.0 ) then
                  write(line,*)
     1                 'update_b: Warning: get_zwzv ends with IERR = ',
     2                 IERR
                  call C_OUT(2,0,1,line)
               endif
               call DAXPY(NIND, ALPHA, RW(p_vout+1), 1, RW(p_y+1), 1)
               p_rwend = p_vout
            endif
C
C           do the BFGS update
C
C           dot = ALPHA*PZ(1:NIND)'*y = s'*y
C
            dot = ALPHA*DDOT(NIND, PZ, 1, RW(p_y+1), 1)
            pznrm = DNRM2(NIND, PZ, 1)
            snorm = ALPHA*pznrm
            ynorm = DNRM2(NIND, RW(p_y+1), 1)
C
C     Check skipping criterion
C
CTODO Check skipping criterion!
            if( QSKIP.gt.0 ) then
               skip_rule = (PYNRM.gt.QSKIP*((0.98)**(dble(ITER)))*pznrm)
            else
               skip_rule = .false.
            endif

            if( (dot.gt.1.d-8*snorm*ynorm .or. abs(QQUASI).eq.5) .and.
     1           .not.skip_rule ) then

               if( abs(QCG).eq.1 .and. abs(QQUASI).eq.5) then
                  call C_OUT(2,0,1,
     1                 'update_b: Powell damping not implemented.')
                  IERR = 4
                  goto 9999
               endif
C
C     Do late correction for B_0
C
CTODO See if this is the best place to put this:
               if( ITER.eq.1 .and. QINITB.eq.5 .and. QCG.eq.0 ) then
                  if( abs(QQUASI).eq.5 ) then
                     call C_OUT(2,0,1,
     1  'update_b: Don''t know how to initialize for Powell damping...')
                     IERR = 4
                     goto 9999
                  endif
                  call DSCAL( ((NIND+1)*NIND)/2, ynorm**2/dot, B, 1)
               endif

               if( abs(QCG).ne.1 ) then
                                ! Regular BFGS update
                  p_bs = p_rwend
                  p_rwend = p_bs + NIND
                  if( p_rwend.gt.LRW ) then
                     IERR = 98
                     goto 9999
                  endif
C              bs = B*s (using symmetry of B)
                  call DSPMV('U', NIND, ALPHA, B, PZ, 1, 0.d0,
     1                 RW(p_bs+1), 1)

C              sBs = s'*bs
                  sBs = ALPHA*DDOT(NIND, PZ, 1, RW(p_bs+1), 1)
C
C     Do Powell damping if desired (s. p.540 Nocedal/Wright's book)
C
                  if (abs(QQUASI).eq.5 .and. dot.lt.0.2d0*sBs) then
                     theta = 0.8d0*sBs/(sBs - dot)
                     call DSCAL(NIND, theta, RW(p_y+1), 1)
                     call DAXPY(NIND, (1d0-theta), RW(p_bs+1), 1,
     1                    RW(p_y+1), 1)
                     dot = theta*dot + (1d0-theta)*sBs
                     if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
                        write(line,*)
     1                       'update_b: Powell damping theta = ',
     1                       theta
                        call C_OUT(2,2,1,line)
                        write(line,*) 'update_b: damped s^T y = ', dot
                        call C_OUT(2,2,1,line)
                     endif
                  endif

C              B_tmp = B + y*y'/dot
                  call DSPR('U', NIND, 1.d0/dot, RW(p_y+1), 1, B)

C              B_new = B_tmp -bs*bs'/dot
                  call DSPR('U', NIND, -1.d0/sBs, RW(p_bs+1), 1, B)

               else
                                ! use inverse BFGS formula
                  p_bs    = p_rwend
                  p_rwend = p_bs + NIND
                  if( p_rwend.gt.LRW ) then
                     IERR = 98
                     goto 9999
                  endif
C     bs = B * y
                  call DSPMV('U', NIND, 1.d0, B, RW(p_y+1), 1, 0.d0,
     1                 RW(p_bs+1), 1)
C     sBs = y' * B * y
                  sBs = DDOT(NIND, RW(p_y+1), 1, RW(p_bs+1), 1)
                  fact = ALPHA**2 * (sBs/dot + 1.d0)/dot
C     B <- B + fact s s^T
                  call DSPR('U', NIND, fact, PZ, 1, B)
C     B <- B - rho s (bs)^T - rho bs s^T
                  fact = -ALPHA/dot
                  call DSPR2('U', NIND, fact, RW(p_bs+1), 1, PZ, 1, B)
                  p_rwend = p_bs
               endif

            elseif( QSKIP.gt.0 .and. skip_rule ) then

               write(line,'(a,d12.5,a,d12.5)')
     1               'Skip BFGS, since PYNRM = ',PYNRM,
     1               ', lhs = ',QSKIP*((0.98)**(dble(ITER)))*pznrm
               call C_OUT(1,1,1,line)
               C_SKIP = '+'

            else
               write(line,700) ITER, dot
 700           format(' In Iter ',i5,
     1                ' no BFGS update since dot = ',
     1                d11.3)
               call C_OUT(1,2,0,line)
               call C_OUT(1,2,1,line)
               C_SKIP = '*'
            endif

            p_rwend = p_y
C
C     SR1
C
         elseif( abs(QQUASI).eq.2 .or. abs(QQUASI).eq.4 ) then

            p_y     = p_rwend
            p_bs    = p_y + NIND
            p_rwend = p_bs + NIND
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
C
C           compute y for SR1 formula (note: s = ALPHA*PZ)
C
            call COMPUTE_Y(N, NIND, M, ITER, IVAR, NFIX, IFIX, NORIG,
     1           XORIG, CSCALE, RG, RGOLD, WBAR, ALPHA, G, GOLD,
     2           LAM, RW(p_y+1), KCONSTR,
     3           LRS, RS, LIS, IS, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.gt.0 ) then
               write(line,*)
     1              'update_b: Error: compute_y ends with IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            elseif( IERR.ne.0 ) then
               write(line,*)
     1              'update_b: Warning: compute_y ends with IERR = ',
     2              IERR
               call C_OUT(2,0,1,line)
            endif
C           bs = B*s (using symmetry of B)
            call DSPMV('U', NIND, ALPHA, B, PZ, 1, 0.d0, RW(p_bs+1), 1)

C           bs = (B*s - y)
            call DAXPY(NIND, -1.d0, RW(p_y+1), 1, RW(p_bs+1), 1)

C           dot = s'*bs
            dot = ALPHA*DDOT(NIND,PZ,1,RW(p_bs+1),1)

C           snorm = |s|, ynorm = |y-Bs|
            snorm = ALPHA*DNRM2(NIND, PZ, 1)
            ynorm = DNRM2(NIND, RW(p_bs+1), 1)

            if( dabs(dot).gt.dmax1(1.d-300,1.d-8*snorm*ynorm) ) then
C            if( dabs(dot).gt.dmax1(1.d-300,1.d-4*snorm*ynorm) ) then

C
C     Do late correction for B_0
C
CTODO See if this is the best place to put this:
               if( ITER.eq.1 .and. QINITB.eq.5 ) then
                  tmp = DDOT(NIND, RW(p_y+1), 1, RW(p_y+1), 1)
                  tmp2 = ALPHA*DDOT(NIND, PZ, 1, RW(p_y+1), 1)
                  call DSCAL( ((NIND+1)*NIND)/2, tmp/tmp2, B, 1)
               endif

C              B_new = B_tmp -bs*bs'/dot
               call DSPR('U', NIND, -1.d0/dot, RW(p_bs+1), 1, B)

            else
               write(line,710) ITER, dot, snorm, ynorm
 710           format(' In Iter ',i5,
     1              ' no SR1 update since dot = ',
     1              3d11.3)
               call C_OUT(1,2,0,line)
c               call C_OUT(1,2,1,line)
               C_SKIP = '*'
            endif

            p_rwend = p_y

         elseif( QQUASI.ne.0 ) then
            write(line,*) 'Invalid Choice of QQUASI = ',QQUASI
            call C_OUT(2,0,1,line)
            stop
         endif

      else
         C_SKIP = '*'
      endif

 9000 continue
C
C     Add W = B + BB (oldW)
C
      if( QFULL.ne.1 ) then
         if( QQUASI.ne.0 .and. QCG.eq.0 ) then
            call DAXPY(INT(NIND*(NIND+1)/2), 1.d0, B, 1, W, 1)
         elseif( QCG.eq.0 ) then
            C_SKIP = '-'
         endif
      endif
 9020 continue

      SKIP_UPDATE = .false.

 9999 continue
      return

 9008 IERR = 8
      call C_OUT(2,0,1,'Error while trying to open BWARM.DAT')
      goto 9999
 9009 IERR = 8
      call C_OUT(2,0,1,'Error while reading BWARM.DAT')
      goto 9999
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine UPDATE_B_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1, lrw2, liw2

      LIW = 0
      LRW = 0

      if( QCG.eq.0 ) then
         if( QINITB.ge.2 .and. QINITB.le.4 ) then
            call GET_G_WS(N, M, NLB, NUB, NZA, lrw1, liw1)
            lrw1 = lrw1 + N
            if( M.gt.0 ) then
               call CONSTR_WS(N, M, NLB, NUB, NZA, lrw2, liw2,
     1              DAT, IDAT)
               lrw1 = max(lrw1, lrw2)
               liw1 = max(liw1, liw2)
            endif
            LIW = liw1
            LRW = 2*N+N-M + lrw1
         else
            LRW = N-M
         endif
      endif

      if( abs(QQUASI).eq.1 .or. abs(QQUASI).eq.3 .or.
     1     abs(QQUASI).eq.5 ) then

         call COMPUTE_Y_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
         if( abs(QCG).eq.1 .and. (NLB.gt.0 .or. NUB.gt.0) ) then
            call GET_ZWZV_WS(N, M, NLB, NUB, NZA, lrw2, liw2, DAT, IDAT)
            lrw1 = max(lrw1, lrw2+(N-M))
            liw1 = max(liw1, liw2)
         endif
         if( abs(QCG).ne.1 ) then
            lrw1 = max(lrw1, N-M)
         endif
         LRW = max(LRW, N-M+lrw1)
         LIW = max(LIW, liw1)
      else
         if( abs(QQUASI).eq.2 .or. abs(QQUASI).eq.4 ) then
            call COMPUTE_Y_WS(N, M, NLB, NUB, NZA, lrw1, liw1,
     1           DAT, IDAT)
            LRW = max(LRW, 2*(N-M)+lrw1)
            LIW = max(LIW, liw1)
         endif
      endif

      return
      end
