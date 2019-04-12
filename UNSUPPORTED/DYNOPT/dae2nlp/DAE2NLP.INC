C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C
C     $Id: DAE2NLP.INC 531 2004-03-11 01:31:07Z andreasw $
C
C     Author:  Andreas Waechter
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA 15232
C
C     Internal values for the collocation
      integer NCOLD2NMAX
      parameter ( NCOLD2NMAX = 6 )
      double precision OMEGA1(NCOLD2NMAX)
      double precision OMEGAQ(NCOLD2NMAX,NCOLD2NMAX)
      double precision RHO(NCOLD2NMAX), COEF(NCOLD2NMAX,NCOLD2NMAX)
      common /DAENLP/ OMEGA1, OMEGAQ, RHO, COEF
      save /DAENLP/
