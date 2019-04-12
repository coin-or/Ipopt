C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

#ifdef OLDADOLC
C
C     Leave this in for possible debugging
C
      subroutine LAGVEC_OLD(NZ, NY, NU, NP, NCOL, NAC, NELE, TI, X, LAM,
     1     V, HV, zval, zeros, hvtmp, ptmp)
C
C     $Id: lagvec.F 531 2004-03-11 01:31:07Z andreasw $
C
C     Compute Hessian-Vector product
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      implicit none
      integer NZ, NY, NU, NP, NCOL, NAC, NELE
      double precision TI(NELE+1)
C      double precision X(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      double precision X(*)
      double precision LAM(NZ + NELE*(NCOL*(NZ+NY)+NZ+NAC))
C      double precision V(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      double precision V(*)     ! ordered like X
C      double precision HV(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      double precision HV(*)    ! ordered like X
      double precision zval(NZ) ! work space
      double precision zeros(*) ! work space of dimension max(NZ,NY,NU,NP)
C      double precision hvtmp(2*NZ+NY+NU+NP) ! work space
      double precision hvtmp(*) ! work space of dimension NZ+NY+NU+NP
      double precision ptmp(NZ) ! work space

      include 'DAE2NLP.INC'

      integer lw, lu, lp, llam, iele, i, r, k
      double precision h, hrho, trho, dummy, olrk

c      WRITE(*,*) 'Entering old version of LAGVEC.'

      call DCOPY(max(NZ,NY,NU,NP), 0d0, 0, zeros, 1)

      lw   = 1
      lu   = NZ + NELE*(NCOL*(NZ+NY)+NZ) + 1
      lp   = lu + NELE*NCOL*NU
      llam = NZ + 1
C
C     Initialize product to zero
C
      call DCOPY(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP, 0d0, 0, HV, 1)
C
C     Loop over all elements
C
      do iele = 1, NELE
C
C     Now loop over the collocation equations
C
         h = TI(iele+1) - TI(iele)
         do k = 1, NCOL
C
C     Compute values of Z at that collocation point
C
            hrho = h * RHO(k)
            trho = TI(iele) + hrho
            call APPROX(1, NZ, NY, NU, NELE, NCOL, COEF, X(lw), X(lu),
     1           TI, iele, trho, OMEGAQ(1,k), dummy, zval, dummy, dummy)
C
C     Compute p_0 's
C
            call DAEMODEL_HV_OLD(NZ, NY, NU, NP, trho, zval,
     1           X(lw+NZ+(k-1)*(NZ+NY)),
     1           X(lw+NZ+(k-1)*(NZ+NY)+NZ), X(lu+(k-1)*NU), X(lp),
     2           LAM(llam+(k-1)*(NZ+NY)), V(lw),
     3           V(lw+NZ+(k-1)*(NZ+NY)),
     3           V(lw+NZ+(k-1)*(NZ+NY)+NZ), V(lu+(k-1)*NU), V(lp),
     4           hvtmp(1), hvtmp(NZ+1), hvtmp(2*NZ+1), hvtmp(2*NZ+NY+1),
     5           hvtmp(2*NZ+NY+NU+1))
C
C     and add them to product
C
            call DCOPY(NZ, hvtmp, 1, ptmp, 1)
            call DAXPY(NZ, 1d0, hvtmp(NZ+1), 1,
     1           HV(lw+NZ+(k-1)*(NZ+NY)), 1)
            call DAXPY(NY, 1d0, hvtmp(2*NZ+1), 1,
     1           HV(lw+NZ+(k-1)*(NZ+NY)+NZ), 1)
            call DAXPY(NU, 1d0, hvtmp(2*NZ+NY+1), 1, HV(lu+(k-1)*NU), 1)
            call DAXPY(NP, 1d0, hvtmp(2*NZ+NY+NU+1), 1, HV(lp), 1)
C
C     Now compute the p_r 's
C
            do r = 1, NCOL
C
               call DAEMODEL_HV_OLD(NZ, NY, NU, NP, trho, zval,
     1              X(lw+NZ+(k-1)*(NZ+NY)), 
     1              X(lw+NZ+(k-1)*(NZ+NY)+NZ), X(lu+(k-1)*NU), X(lp),
     2              LAM(llam+(k-1)*(NZ+NY)), V(lw+NZ+(r-1)*(NZ+NY)),
     3              zeros, zeros, zeros, zeros,
     4              hvtmp(1), hvtmp(NZ+1), hvtmp(2*NZ+1),
     5              hvtmp(2*NZ+NY+1), hvtmp(2*NZ+NY+NU+1))
C
C     and add them to product
C
               olrk = hrho*OMEGAQ(r,k)
               call DAXPY(NZ, olrk, hvtmp, 1, ptmp, 1)
               call DAXPY(NZ, olrk, hvtmp(NZ+1), 1,
     1              HV(lw+NZ+(k-1)*(NZ+NY)), 1)
               call DAXPY(NY, olrk, hvtmp(2*NZ+1), 1,
     1              HV(lw+NZ+(k-1)*(NZ+NY)+NZ), 1)
               call DAXPY(NU, olrk, hvtmp(2*NZ+NY+1), 1,
     1              HV(lu+(k-1)*NU), 1)
               call DAXPY(NP, olrk, hvtmp(2*NZ+NY+NU+1), 1, HV(lp), 1)
C
            enddo
C
C     There are still the zdot's to be dealt with:
C
            call DAXPY(NZ, 1d0, ptmp, 1, HV(lw), 1)
            do r = 1, NCOL
C
               olrk = hrho*OMEGAQ(r,k)
               call DAXPY(NZ, olrk, ptmp, 1,
     1              HV(lw+NZ+(r-1)*(NZ+NY)), 1)
C
            enddo

         enddo
C
         llam = llam + NZ + NCOL*(NZ+NY) + NAC
C
         lw = lw + NCOL*(NZ+NY) + NZ
         lu = lu + NCOL*NU

      enddo

      return
      end
#else
C
C     New VERSION
C
      subroutine LAGVEC(NZ, NY, NU, NP, NCOL, NAC, NELE, TI, X, LAM, V,
     1     HV, zval, rhs, res, ptmp)
C
C     $Id: lagvec.F 531 2004-03-11 01:31:07Z andreasw $
C
C     Compute Hessian-Vector product
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      implicit none
      integer NZ, NY, NU, NP, NCOL, NAC, NELE
      double precision TI(NELE+1)
C      double precision X(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      double precision X(*)
      double precision LAM(NZ + NELE*(NCOL*(NZ+NY)+NZ+NAC))
C      double precision V(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      double precision V(*)     ! ordered like X
C      double precision HV(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      double precision HV(*)    ! ordered like X
      double precision zval(NZ) ! work space
      double precision rhs((NCOL+1)*(2*NZ+NY+NU+NP))
                                ! work space
C      double precision res((NCOL+1)*(2*NZ+NY+NU+NP))
      double precision res(*)
                                ! work space
      double precision ptmp(NZ) ! work space

      include 'DAE2NLP.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DAENLP/

      integer lw, lu, lp, llam, iele, i, r, k
      double precision h, hrho, trho, dummy, olrk

C
C     initialize RHS to zero (need to do this only once)
C
      call DCOPY((NCOL+1)*(2*NZ+NY+NU+NP), 0d0, 0, RHS, 1)

      lw   = 1
      lu   = NZ + NELE*(NCOL*(NZ+NY)+NZ) + 1
      lp   = lu + NELE*NCOL*NU
      llam = NZ + 1
C
C     Initialize product to zero
C
      call DCOPY(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP, 0d0, 0, HV, 1)
C
C     Loop over all elements
C
      do iele = 1, NELE
C
C     Now loop over the collocation equations
C
         h = TI(iele+1) - TI(iele)
         do k = 1, NCOL
C
C     Compute values of Z at that collocation point
C
            hrho = h * RHO(k)
            trho = TI(iele) + hrho
            call APPROX(1, NZ, NY, NU, NELE, NCOL, COEF, X(lw), X(lu),
     1           TI, iele, trho, OMEGAQ(1,k), dummy, zval, dummy, dummy)
C
C     Construct right hand side for computation of p_0 and p_r's
C
            call DCOPY(NZ, V(lw), 1, RHS(1), 1)
            call DCOPY(NZ, V(lw+NZ+(k-1)*(NZ+NY)), 1, RHS(1+NZ), 1)
            call DCOPY(NY, V(lw+NZ+(k-1)*(NZ+NY)+NZ), 1, RHS(1+2*NZ), 1)
            call DCOPY(NU, V(lu+(k-1)*NU), 1, RHS(1+2*NZ+NY), 1)
            call DCOPY(NP, V(lp), 1, RHS(1+2*NZ+NY+NU), 1)
            do r = 1, NCOL
               call DCOPY(NZ, V(lw+NZ+(r-1)*(NZ+NY)), 1,
     1              RHS(1+r*(2*NZ+NY+NU+NP)), 1)
            enddo
C
C     Compute the product
C
            call DAEMODEL_HV(NCOL+1, NZ, NY, NU, NP, trho, zval,
     1           X(lw+NZ+(k-1)*(NZ+NY)),
     1           X(lw+NZ+(k-1)*(NZ+NY)+NZ), X(lu+(k-1)*NU), X(lp),
     2           LAM(llam+(k-1)*(NZ+NY)), RHS, RES, 2*NZ+NY+NU+NP)
C
C     Add p_0 to product
C
            call DCOPY(NZ, RES(1), 1, ptmp, 1)
            call DAXPY(NZ, 1d0, RES(NZ+1), 1,
     1           HV(lw+NZ+(k-1)*(NZ+NY)), 1)
            call DAXPY(NY, 1d0, RES(2*NZ+1), 1,
     1           HV(lw+NZ+(k-1)*(NZ+NY)+NZ), 1)
            call DAXPY(NU, 1d0, RES(2*NZ+NY+1), 1, HV(lu+(k-1)*NU), 1)
            call DAXPY(NP, 1d0, RES(2*NZ+NY+NU+1), 1, HV(lp), 1)
C
C     Now add the p_r's to product
C
            do r = 1, NCOL
C
               olrk = hrho*OMEGAQ(r,k)
               call DAXPY(NZ, olrk,
     1              RES(1+r*(2*NZ+NY+NU+NP)), 1, ptmp, 1)
               call DAXPY(NZ, olrk,
     1              RES(NZ+1+r*(2*NZ+NY+NU+NP)), 1,
     1              HV(lw+NZ+(k-1)*(NZ+NY)), 1)
               call DAXPY(NY, olrk,
     1              RES(2*NZ+1+r*(2*NZ+NY+NU+NP)), 1,
     1              HV(lw+NZ+(k-1)*(NZ+NY)+NZ), 1)
               call DAXPY(NU, olrk,
     1              RES(2*NZ+NY+1+r*(2*NZ+NY+NU+NP)), 1,
     1              HV(lu+(k-1)*NU), 1)
               call DAXPY(NP, olrk,
     1              RES(2*NZ+NY+NU+1+r*(2*NZ+NY+NU+NP)), 1, HV(lp), 1)
C
            enddo
C
C     There are still the zdot's to be dealt with:
C
            call DAXPY(NZ, 1d0, ptmp, 1, HV(lw), 1)
            do r = 1, NCOL
C
               olrk = hrho*OMEGAQ(r,k)
               call DAXPY(NZ, olrk, ptmp, 1,
     1              HV(lw+NZ+(r-1)*(NZ+NY)), 1)
C
            enddo

         enddo
C
         llam = llam + NZ + NCOL*(NZ+NY) + NAC
C
         lw = lw + NCOL*(NZ+NY) + NZ
         lu = lu + NCOL*NU

      enddo

      return
      end
#endif
