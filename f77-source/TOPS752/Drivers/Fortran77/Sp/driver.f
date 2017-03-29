C
C
C        SRFTEST:  Portable Test Driver for SRFPACK
C                        07/09/98
C
C
C   This driver tests software package SRFPACK for construc-
C ting an interpolatory surface to data values at scattered
C points in the plane.  The software modules are tested on
C data for which the methods are exact, and maximum errors
C are printed.  These should be close to the machine pre-
C cision.  The final test using nonzero tension factors is
C included primarily to demonstrate how the modules are
C called.  Given the small number of nodes, large interpo-
C lation errors are expected in this case.
C
C   By default, tests are performed on a simple data set
C consisting of 12 nodes whose convex hull is the unit
C square with a .3 by .3 square constraint region in the
C center.  However, by enabling the READ statements below
C (C# in the first two columns), testing may be performed on
C an arbitrary set of up to NMAX nodes with up to NCMAX
C constraint curves.  (Refer to the PARAMETER statements be-
C low.)  Note that the choice of interpolation/extrapolation
C points assumes the nodes to approximately cover the unit
C square.  A data set consists of the following sequence of
C records:
C
C    N = Number of nodes (format I4):  3 to NMAX.
C    NCC = Number of constraint curves (format I4):  0 to
C          NCMAX.
C    (LCC(I), I = 1,NCC) = Indexes of the first node in each
C             constraint curve (format I4).  4 .LE. LCC(1)
C             and, for I .GT. 1, LCC(I-1) + 3 .LE. LCC(I)
C             .LE. N-2. (Each constraint curve has at least
C             three nodes.)
C    (X(I),Y(I), I = 1,N) = Nodal coordinates with non-
C                constraint nodes followed by the NCC
C                sequences of constraint nodes (format
C                2F13.8).
C
C   Various options may be changed by altering the DATA
C statements below.
C
C This driver must be linked to TRIPACK and SRFPACK.
C
      INTEGER I, IER, IFLGS, IMAX, IST, J, K, L, LIN, LNEW,
     .        LOUT, LPLT, LW, LWK, MAXITG, MAXITZ, N, N6,
     .        NA, NB, NCC, NCMAX, NCON, NI, NITG, NITZ, NJ,
     .        NLS, NMAX, NT
      LOGICAL DFLAG, SFLAG
      REAL    A1, A2, A3, DGMAX, DGMX, DMAX, DSM, DZMAX,
     .        DZMX, PLTSIZ, PX, PY, PZ, SVAL, TOL, XF, YF,
     .        ZTRUE, ZX, ZY
      REAL    AREAP, FT, VOLUME
C
      PARAMETER (NMAX=5000, NCMAX=5, N6=6*NMAX, LWK=2*NMAX,
     .           NI=13, NJ=13)
C
C Array storage:
C
      INTEGER LCC(NCMAX), LEND(NMAX), LIST(N6), LPTR(N6),
     .        NODES(LWK)
      REAL    GRAD(2,NMAX), P(NI), SIGMA(N6),
     .        X(NMAX), Y(NMAX), Z(NMAX), ZZ(NI,NJ)
C
C Plot size and number of contours for the contour plot.
C
      DATA    PLTSIZ/7.5/,  NCON/8/
C
C Default node set:
C
      DATA   N/12/, NCC/1/, LCC(1)/9/
      DATA    X(1),  X(2),  X(3),  X(4),  X(5),  X(6),  X(7)
     .       / 0.,    1.,    .5,   .15,   .85,    .5,   0./,
     .        X(8),  X(9), X(10), X(11), X(12)
     .       / 1.,   .35,   .65,   .65,   .35/,
     .        Y(1),  Y(2),  Y(3),  Y(4),  Y(5),  Y(6),  Y(7)
     .       / 0.,    0.,   .15,    .5,    .5,   .85,   1./,
     .        Y(8),  Y(9), Y(10), Y(11), Y(12)
     .       / 1.,   .35,   .35,   .65,   .65/
C
C Evaluation points:
C
      DATA    P(1),  P(2),  P(3),  P(4),  P(5),  P(6),  P(7)
     .       /-0.1,   0.0,   0.1,   0.2,   0.3,   0.4,  .5/,
     .        P(8),  P(9), P(10), P(11), P(12), P(13)
     .       / 0.6,   0.7,   0.8,   0.9,   1.0,   1.1/
C
C Test function for ZGRADG, GETSIG, and CRPLOT:
C
      FT(XF,YF) = (TANH(9.*(YF-XF)) + 1.)/9.
C
C Logical unit numbers for I/O:
C
      DATA    LIN/1/,  LOUT/2/,  LPLT/3/
      OPEN (LOUT,FILE='RES')
      OPEN (LPLT,FILE='RES.eps')
C
C *** Read triangulation parameters -- N, NCC, LCC, X, Y.
C
C#    OPEN (LIN,FILE='srftest.dat',STATUS='OLD')
C#    READ (LIN,100,ERR=30) N, NCC
      IF (N .LT. 3  .OR.  N .GT. NMAX  .OR.  NCC .LT. 0
     .    .OR.  NCC .GT. NCMAX) GO TO 31
C#    IF (NCC .GT. 0) READ (LIN,100,ERR=30)
C#   .                     (LCC(I), I = 1,NCC)
C#    READ (LIN,110,ERR=30) (X(I),Y(I), I = 1,N)
C#100 FORMAT (I4)
C#110 FORMAT (2F13.8)
C
C Set NLS to the number of nonconstraint nodes.
C
      IF (NCC .EQ. 0) THEN
        NLS = N
      ELSE
        NLS = LCC(1) - 1
      ENDIF
C
C Print a heading, and construct the triangulation.
C   NODES and GRAD are used as work space.
C
      WRITE (LOUT,400) N
      CALL TRMESH (N,X,Y, LIST,LPTR,LEND,LNEW,NODES,
     .             NODES(N+1),GRAD,IER)
      IF (IER .EQ. -2) THEN
        WRITE (LOUT,210)
        STOP
      ELSEIF (IER .EQ. -4) THEN
        WRITE (LOUT,212)
        STOP
      ELSEIF (IER .GT. 0) THEN
        WRITE (LOUT,214)
        STOP
      ENDIF
C
C *** Add the constraint curves.  Note that edges and
C     triangles are not removed from constraint regions.
C     ADDCST forces the inclusion of triangulation edges
C     connecting the sequences of constraint nodes.  If it
C     is necessary to alter the triangulation, the empty
C     circumcircle property is no longer satisfied.  NODES
C     is used as a work space array.
C
      LW = LWK
      CALL ADDCST (NCC,LCC,N,X,Y, LW,NODES,LIST,LPTR,
     .             LEND, IER)
      IF (IER .NE. 0) THEN
        WRITE (LOUT,220) IER
        STOP
      ENDIF
C
C *** Test VOLUME by comparing the values returned by AREAP
C                 applied to the triangulation boundary and
C                 VOLUME applied to the constant function
C                 Z = 1.
C
      CALL BNODES (N,LIST,LPTR,LEND, NODES,NB,NA,NT)
      A1 = AREAP(X,Y,NB,NODES)
      DO 1 K = 1,N
        Z(K) = 1.
        GRAD(1,K) = 0.
        GRAD(2,K) = 0.
    1   CONTINUE
      A2 = VOLUME(0,LCC,N,X,Y,Z,LIST,LPTR,LEND)
      A3 = VOLUME(NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND)
      WRITE (LOUT,300) NB, NA, NT, A1, A2, A3
C
C *** Test GRADG on the constant function with no tension.
C
      IFLGS = 0
      SIGMA(1) = 0.
      NITG = 1
      DGMX = 0.
      CALL GRADG (NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,
     .            IFLGS,SIGMA, NITG,DGMX,GRAD, IER)
      IF (IER .LT. 0) THEN
        WRITE (LOUT,230)
        STOP
      ENDIF
      DMAX = 0.
      DO 2 K = 1,N
        DMAX = MAX(DMAX,ABS(GRAD(1,K)),ABS(GRAD(2,K)))
    2   CONTINUE
      WRITE (LOUT,310) DMAX
C
C *** Test INTRC0 on the constant function (for which both
C                 interpolation and extrapolation are
C                 exact).
C
      DMAX = 0.
      IST = 1
      DO 4 J = 1,NJ
        PY = P(J)
        DO 3 I = 1,NI
          PX = P(I)
          CALL INTRC0 (PX,PY,NCC,LCC,N,X,Y,Z,LIST,
     .                 LPTR,LEND, IST, PZ,IER)
          DMAX = MAX(DMAX,ABS(PZ-1.))
    3     CONTINUE
    4   CONTINUE
      WRITE (LOUT,320) DMAX
C
C *** Test GRADL on a quadratic function.  The true partial
C                derivatives are stored in ZX and ZY.
C
      DO 5 K = 1,N
        Z(K) = ((X(K)+Y(K))**2)/4.
    5   CONTINUE
      DMAX = 0.
      DO 6 K = 1,N
        ZX = (X(K)+Y(K))/2.
        ZY = ZX
        CALL GRADL (K,NCC,LCC,N,X,Y,Z,LIST,LPTR,
     .              LEND, GRAD(1,K),GRAD(2,K),IER)
        DMAX = MAX(DMAX,ABS(ZX-GRAD(1,K)),
     .                  ABS(ZY-GRAD(2,K)))
    6   CONTINUE
      WRITE (LOUT,330) DMAX
C
C *** Test UNIF (and INTRC1) on the linear function (for
C               which both interpolation and extrapolation
C               are exact).
C
      DO 7 K = 1,N
        Z(K) = (X(K)+Y(K))/2.
        GRAD(1,K) = .5
        GRAD(2,K) = .5
    7   CONTINUE
      IFLGS = 0
      SIGMA(1) = 0.
      SFLAG = .FALSE.
      CALL UNIF (NCC,LCC,N,X,Y,Z,GRAD,LIST,LPTR,LEND,
     .           IFLGS,SIGMA,NI,NI,NJ,P,P,SFLAG,SVAL, ZZ,
     .           IER)
      DMAX = 0.
      DO 9 J = 1,NJ
        DO 8 I = 1,NI
          ZTRUE = (P(I)+P(J))/2.
          DMAX = MAX(DMAX,ABS(ZTRUE-ZZ(I,J)))
    8     CONTINUE
    9   CONTINUE
      WRITE (LOUT,340) DMAX, IER, NI*NJ
C
C *** Test GETSIG, ZINIT, and ZGRADG by computing tension
C     factors SIGMA, global gradient estimates GRAD, and
C     constraint node values (Z(I), I = LCC(1) to N) using
C     data values (Z(I), I = 1 to NLS) from test function
C     FT.  Gradient estimates are initialized to zeros.
C
      DO 10 K = 1,NLS
        Z(K) = FT(X(K),Y(K))
   10   CONTINUE
      DO 11 K = 1,N
        GRAD(1,K) = 0.
        GRAD(2,K) = 0.
   11   CONTINUE
C
C ZINIT is used to get good initial estimates of constraint-
C   node values for ZGRADG.
C
      CALL ZINIT (NCC,LCC,N,X,Y,LIST,LPTR,LEND, Z, IER)
C
C Initialize tension factors to zeros.  L is the length of
C   LIST, LPTR, and SIGMA.
C
      L = LNEW - 1
      DO 12 K = 1,L
        SIGMA(K) = 0.
   12   CONTINUE
C
C The global gradient estimates depend on the tension
C   factors, and the optimal tension factors depend on the
C   gradient estimates.  This requires an iteration around
C   the calls to ZGRADG/GRADG and GETSIG.  An appropriate
C   convergence test would involve an upper bound on the
C   relative change in the gradients between iterations.
C   However, a few iterations is usually sufficient to
C   produce good results.
C
C Parameters:
C
C   MAXITZ = Maximum number of Gauss-Seidel iterations to be
C            used by ZGRADG.
C   DZMAX = Tolerance defining convergence in ZGRADG.
C   MAXITG = Maximum number of Gauss-Seidel iterations to be
C            used by GRADG.
C   DGMAX = Tolerance defining convergence in GRADG.
C   TOL = Tolerance for GETSIG.
C   IMAX = Number of ZGRADG/GETSIG iterations.
C   IFLGS = Tension factor option:  0 means uniform tension.
C
      MAXITZ = 10
      DZMAX = 1.E-3
      MAXITG = 10
      DGMAX = 1.E-3
      TOL = 1.E-3
      IMAX = 3
      DO 13 I = 1,IMAX
        IFLGS = 0
        IF (I .GT. 1) IFLGS = 1
C
C Compute global gradient estimates and constraint node
C   values (ZGRADG).
C
        IF (NCC .GT. 0  .AND.  MAXITZ .GT. 0) THEN
          NITZ = MAXITZ
          DZMX = DZMAX
          CALL ZGRADG (NCC,LCC,N,X,Y,LIST,LPTR,LEND,
     .                 IFLGS,SIGMA, NITZ,DZMX,Z,
     .                 GRAD, IER)
          IF (IER .LT. 0) THEN
            WRITE (LOUT,240) IER, NITZ
            STOP
          ENDIF
        ELSE
          NITZ = 0
        ENDIF
C
C Improve the gradient estimates (GRADG).
C
        IF (MAXITG .GT. 0) THEN
          NITG = MAXITG
          DGMX = DGMAX
          CALL GRADG (NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,
     .                IFLGS,SIGMA, NITG,DGMX,GRAD, IER)
        ELSE
          NITG = 0
        ENDIF
C
C Compute tension factors SIGMA.
C
        CALL GETSIG (N,X,Y,Z,LIST,LPTR,LEND,GRAD,
     .               TOL, SIGMA, DSM,IER)
   13   CONTINUE
C
C Compute interpolated values and the maximum error on the
C   uniform grid.
C
      IFLGS = 1
      DFLAG = .FALSE.
      IST = 1
      DMAX = 0.
      DO 15 J = 1,NJ
        PY = P(J)
        DO 14 I = 1,NI
          PX = P(I)
          CALL INTRC1 (PX,PY,NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,
     .                 IFLGS,SIGMA,GRAD,DFLAG, IST, ZZ(I,J),
     .                 ZX,ZY,IER)
          IF (IER .LT. 0) THEN
            WRITE (LOUT,250) IER
            STOP
          ENDIF
C
C   Exclude extrapolation errors.
C
          IF (IER .GT. 0) GO TO 14
          ZTRUE = FT(PX,PY)
          DMAX = MAX(DMAX,ABS(ZTRUE-ZZ(I,J)))
   14     CONTINUE
   15   CONTINUE
      WRITE (LOUT,350) DMAX
C
C Print SIGMA.
C
      CALL SGPRNT (N,LOUT,LIST,LPTR,LEND,SIGMA)
C
C Create a contour plot.  NODES, X, and Y are used as work
C   space arrays.  It is assumed that NMAX >= 2.5*NI*NJ.
C
      CALL CRPLOT (LPLT,PLTSIZ,NI,NJ,P,P,ZZ,NCON,NODES,
     .             X,Y, IER)
      IF (IER .NE. 0) THEN
        WRITE (LOUT,260) IER
        STOP
      ENDIF
      WRITE (LOUT,360)
C
C Successful test.
C
      STOP
C
C Error reading the data set.
C
C# 30 WRITE (*,200)
      STOP
C
C Invalid value of N or NCC.
C
   31 WRITE (*,205) N, NCC
      STOP
C
C Heading:
C
  400 FORMAT (///1X,21X,'SRFPACK Test:  N =',I4///)
C
C Error message formats:
C
C#200 FORMAT (//5X,'*** Input data set invalid ***'/)
  205 FORMAT (//5X,'*** N or NCC is outside its valid ',
     .             'range:  N =',I5,', NCC = ',I4,' ***'/)
  210 FORMAT (//5X,'*** Error in TRMESH:  the first three ',
     .        'nodes are collinear ***'/)
  212 FORMAT (//5X,'*** Error in TRMESH:  invalid ',
     .        'triangulation ***'/)
  214 FORMAT (//5X,'*** Error in TRMESH:  duplicate nodes ',
     .        'encountered ***'/)
  220 FORMAT (//5X,'*** Error in ADDCST:  IER = ',I1,
     .        ' ***'/)
  230 FORMAT (//5X,'*** Error in GRADG (invalid flag ',
     .        'returned):  IER = ',I2,' ***'/)
  240 FORMAT (//5X,'*** Error in ZGRADG:  IER = ',
     .        I2,', NIT = ',I3,' ***'/)
  250 FORMAT (//5X,'*** Error in INTRC1 (invalid flag ',
     .        'returned):  IER = ',I2,' ***'/)
  260 FORMAT (//5X,'*** Error in CRPLOT:  IER = ',I1,
     .        ' ***'/)
C
C Test message formats:
C
  300 FORMAT (//5X,'Output from BNODES, AREA, and VOLUME'//
     .        5X,'BNODES:  ',I4,' Boundary Nodes,  ',I4,
     .        ' Arcs,  ',I4,' Triangles'/5X,
     .        'AREAP:   Area of convex hull = ',E15.8/5X,
     .        'VOLUME:  Area of convex hull = ',E15.8/5X,
     .        'VOLUME:  Area excluding constraint regions',
     .        ' = ',E15.8//)
  310 FORMAT (5X,'GRADG Test:  Maximum error = ',
     .        E15.8//)
  320 FORMAT (5X,'INTRC0 Test:  Maximum error = ',
     .        E15.8//)
  330 FORMAT (5X,'GRADL Test:  Maximum error = ',
     .        E15.8//)
  340 FORMAT (5X,'UNIF (INTRC1) Test:  Maximum ',
     .        'error = ',E15.8/5X,I2,' of ',I3,
     .        ' points required extrapolation'//)
  350 FORMAT (5X,'ZGRADG/GETSIG Test:  Maximum ',
     .        'interpolation error = ',E15.8//)
  360 FORMAT (/5X,'A contour plot file was ',
     .            'successfully created.')
      END
