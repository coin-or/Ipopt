C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/
