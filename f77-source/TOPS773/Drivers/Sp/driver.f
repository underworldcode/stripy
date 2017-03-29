C
C                        SSRFTEST
C                        08/15/96
C
C
C   This program uses SSRFPACK to fit a smoothing or inter-
C polating surface under tension to a set of data values
C at arbitrarily distributed nodes on the surface of a
C sphere.  In addition to functioning as a test program,
C SSRFTEST serves as a demonstration driver and template
C for use of the software package.  Options are selected by
C altering flag values in the DATA statements below.
C
C
C INPUT Data:
C
C   An input data set (with file name data) must be
C provided.  This consists of the number of nodes N
C (Format I5) followed by N ordered triples (RLAT,RLON,W)
C (Format 3F15.9), where RLAT and RLON are latitude and lon-
C gitude in degrees, and W is the associated data value.
C Duplicate nodes (RLAT,RLON) are not allowed, and N must
C be in the range [3,NMAX] (or [7,NMAX] if MFLAG = 1), where
C NMAX is defined in the PARAMETER statement below.
C
C
C OUTPUT Data:
C
C   The fitting function F (selected by the method flag
C MFLAG) is evaluated on an NI by NJ rectangular grid
C defined by lines of latitude PLAT and longitude PLON
C uniformly spaced over the range of input data (the convex
C hull of the nodes).  The output data set (file name
C ssrftest.out) consists of the following records:  NI
C (Format I4), NJ (Format I4), a strictly increasing
C sequence of NI latitude values PLAT in radians (Format
C E15.8), a strictly increasing sequence of NJ longitude
C values PLON in radians (Format E15.8), and the array of
C NI*NJ function values FF (Format E15.8), where FF(K) =
C F(PLAT(I),PLON(J)) for K = (J-1)*NI + I.  The array of
C function values FF reflects the (optional) scaling of the
C data values.  FF is a discrete representation of F
C suitable for input to a surface plotting program.  NI and
C NJ are defined in the PARAMETER statement below.  Both
C must have value at least 2.
C
C   Two additional output files are created:  informative
C messages and error messages are written to ssrftest.msg,
C and tension factors are written to ssrftest.sig if
C VT = TRUE.
C
C
C OPTION Flags:
C
C   The type of surface fit is defined by the method flag
C MFLAG (0 to 3) and the tension factor option VT.
C
C   MFLAG = 0 for piecewise linear interpolation (INTRC0).
C             No tension is used regardless of the value of
C             VT in this case.
C   MFLAG = 1 for smooth interpolation (INTRC1) with local
C             gradient estimates (GRADL).
C   MFLAG = 2 for smooth interpolation (INTRC1) with global
C             gradient estimates (GRADG).
C   MFLAG = 3 for smoothing (SMSURF).
C
C   VT = TRUE if and only if variable tension (a set of non-
C             zero tension factors computed by GETSIG) is to
C             be used.
C
C   SCALE = TRUE iff the data values are to be scaled by
C           1/DW, where DW is the difference between the
C           largest and smallest values.
C
C This program must be linked to STRIPACK and SSRFPACK.
C
      INTEGER NMAX, NMX6, NI, NJ
      PARAMETER (NMAX=1500, NMX6=6*NMAX, NI=32, NJ=32)
C
C Array storage:  F is used iff MFLAG = 3, WT is used if
C                 MFLAG = 3 and (in any case) by TRMESH
C                 as work space, and only one element of
C                 SIGMA is used if VT = FALSE.
C
      INTEGER LIST(NMX6), LPTR(NMX6), LEND(NMAX),
     .        IWK(NMAX,2)
      REAL    RLAT(NMAX), RLON(NMAX), W(NMAX), X(NMAX),
     .        Y(NMAX), Z(NMAX), F(NMAX), WT(NMAX),
     .        GRAD(3,NMAX), SIGMA(NMX6),
     .        PLAT(NI), PLON(NJ), FF(NI,NJ)
C
      INTEGER I, IER, IFLGS, IST, ITER, ITGS, J, K, L, LIN,
     .        LMSG, LNEW, LOUT, LSIG, MAXIT, MFLAG, N, NEV,
     .        NITG, NXP
      LOGICAL SCALE, VT
      REAL    DGMAX, DGMX, DLAT, DLON, DSM, E, GSTOL, RLNMN,
     .        RLNMX, RLTMN, RLTMX, SF, SM, SMTOL, SUM, TOL,
     .        WMN, WMX, WTK
C
C Flag values:
C
      DATA    MFLAG/1/,  VT/.TRUE./,  SCALE/.TRUE./
C
C GRADG tolerance and maximum number of iterations.
C
      DATA    DGMAX/.01/,  MAXIT/10/
C
C GETSIG tolerance:
C
      DATA    TOL/.01/
C
C SMSURF parameters:
C
C   E = Positive expected squared error in a typical
C       (scaled) data value.
C   SM = Smoothing parameter (upper bound on the weighted
C        sum of squares of deviations from the data values).
C        If SM .LE. 0, the value is taken to be N.
C
      DATA    E/.01/,  SM/0./
C
C Number of GRADG/GETSIG or SMSURF/GETSIG iterations used
C   when MFLAG > 1 and VT = TRUE (at least 1):
C
      DATA    ITGS/3/
C
C Logical unit numbers for input, data output, message
C   output, and tension factors:
C
      DATA LIN/1/,  LOUT/2/,  LMSG/3/,  LSIG/4/
      OPEN (LIN,FILE='data',STATUS='OLD')
      OPEN (LOUT,FILE='RES',STATUS='NEW')
      OPEN (LMSG,FILE='MSG')
      IF (VT) OPEN (LSIG,FILE='SIG')
C
C Input Formats:
C
  300 FORMAT (I5)
  310 FORMAT (3F15.9)
C
C Output Formats:
C
  350 FORMAT (I4)
  360 FORMAT (E15.8)
C
C *** Read the input data and determine the range of
C     latitude, longitude, and data values.
C
      READ (LIN,300,ERR=160) N
      IF (N .LT. 3  .OR.  (MFLAG .EQ. 1  .AND.  N .LT. 7)
     .    .OR.  N .GT. NMAX) GO TO 150
      READ (LIN,310,ERR=160) RLAT(1), RLON(1), W(1)
      RLTMN = RLAT(1)
      RLTMX = RLAT(1)
      RLNMN = RLON(1)
      RLNMX = RLON(1)
      WMN = W(1)
      WMX = W(1)
      DO 10 K = 2,N
        READ (LIN,310,ERR=160) RLAT(K), RLON(K), W(K)
        IF (RLAT(K) .LT. RLTMN) RLTMN = RLAT(K)
        IF (RLAT(K) .GT. RLTMX) RLTMX = RLAT(K)
        IF (RLON(K) .LT. RLNMN) RLNMN = RLON(K)
        IF (RLON(K) .GT. RLNMX) RLNMX = RLON(K)
        IF (W(K) .LT. WMN) WMN = W(K)
        IF (W(K) .GT. WMX) WMX = W(K)
   10   CONTINUE
C
C   Print a description of the input data set.
C
      WRITE (LMSG,500) N, RLTMN, RLTMX, RLNMN, RLNMX,
     .                 WMN, WMX
      WRITE (*,500) N, RLTMN, RLTMX, RLNMN, RLNMX,
     .              WMN, WMX
  500 FORMAT (///26X,'SSRFTEST Output'//
     .        5X,'Number of nodes =',I5/
     .        5X,'Range of latitude = ',F6.2,' to ',F6.2/
     .        5X,'Range of longitude = ',F7.2,' to ',F7.2/
     .        5X,'Range of data values = ',F9.4,' to ',F9.4)
C
C *** Store PLAT and PLON (in radians) defining a uniform
C     rectangular grid of NEV evaluation points.
C
      NEV = NI*NJ
      SF = ATAN(1.)/45.
      RLTMN = SF*RLTMN
      RLTMX = SF*RLTMX
      DLAT = (RLTMX-RLTMN)/REAL(NI+1)
      DO 20 I = 1,NI
        PLAT(I) = RLTMN + REAL(I)*DLAT
   20   CONTINUE
      RLNMN = SF*RLNMN
      RLNMX = SF*RLNMX
      DLON = (RLNMX-RLNMN)/REAL(NJ+1)
      DO 25 J = 1,NJ
        PLON(J) = RLNMN + REAL(J)*DLON
   25   CONTINUE
C
C *** Convert the nodes (RLAT,RLON) to radians and then
C     to Cartesian coordinates (X,Y,Z) and create the
C     triangulation.
C
      DO 30 K = 1,N
        RLAT(K) = SF*RLAT(K)
        RLON(K) = SF*RLON(K)
   30   CONTINUE
      CALL TRANS (N,RLAT,RLON, X,Y,Z)
      CALL TRMESH (N,X,Y,Z, LIST,LPTR,LEND,LNEW,IWK,
     .             IWK(1,2),WT,IER)
      IF (IER .NE. 0) THEN
        WRITE (LMSG,420) IER
        WRITE (*,420) IER
        STOP
      ENDIF
C
C *** Scale the data if SCALE = TRUE.
C
      IF (SCALE  .AND.  WMX .GT. WMN) THEN
        SF = 1.0/(WMX-WMN)
        DO 40 K = 1,N
          W(K) = SF*W(K)
   40     CONTINUE
      ENDIF
C
C *** If MFLAG > 0, then store uniform tension factor
C     SIGMA(1) (if VT = FALSE), or initialize variable
C     tension factors to zeros (if VT = TRUE).
C
      IF (MFLAG .GT. 0) SIGMA(1) = 0.
      IF (MFLAG .GT. 0  .AND.  VT) THEN
        L = LNEW-1
        DO 50 K = 1,L
          SIGMA(K) = 0.
   50     CONTINUE
      ENDIF
      IF (MFLAG .EQ. 0) THEN
C
C *** C-0 interpolation (INTRC0).
C
        NXP = 0
        IST = 1
        DO 70 J = 1,NJ
          DO 60 I = 1,NI
            CALL INTRC0 (N,PLAT(I),PLON(J),X,Y,Z,W,LIST,
     .                   LPTR,LEND, IST, FF(I,J),IER)
            IF (IER .GT. 0) NXP = NXP + 1
            IF (IER .LT. 0) THEN
              WRITE (LMSG,430) I, J, IER
              WRITE (*,430) I, J, IER
              STOP
            ENDIF
   60       CONTINUE
   70     CONTINUE
        WRITE (LMSG,510) NEV, NXP
        WRITE (*,510) NEV, NXP
  510   FORMAT (//5X,'INTRC0:  Number of evaluation ',
     .          'points = ',I4/
     .          14X,'Number of extrapolation points =',I4)
      ELSEIF (MFLAG .EQ. 1) THEN
C
C *** C-1 interpolation (INTRC1) with local gradients GRADL.
C
C   Accumulate the sum of the numbers of nodes used in the
C     least squares fits in SUM.
C
        SUM = 0.
        DO 80 K = 1,N
          CALL GRADL (N,K,X,Y,Z,W,LIST,LPTR,LEND, GRAD(1,K),
     .                IER)
          IF (IER .LT. 0) THEN
            WRITE (LMSG,440) N, K, IER
            WRITE (*,440) N, K, IER
            STOP
          ENDIF
          SUM = SUM + REAL(IER)
   80     CONTINUE
        SUM = SUM/REAL(N)
        WRITE (LMSG,520) SUM
        WRITE (*,520) SUM
  520   FORMAT (//5X,'GRADL:  Average number of nodes used',
     .          ' in the least squares fits = ',F4.1)
        IF (VT) THEN
C
C   Compute tension factors SIGMA (GETSIG).
C
          CALL GETSIG (N,X,Y,Z,W,LIST,LPTR,LEND,GRAD,
     .                 TOL, SIGMA, DSM,IER)
          IF (IER .LT. 0) THEN
            WRITE (LMSG,460) IER
            WRITE (*,460) IER
            STOP
          ENDIF
          WRITE (LMSG,530) IER, DSM
          WRITE (*,530) IER, DSM
  530     FORMAT (//5X,'GETSIG:  ',I4,' tension factors ',
     .            'altered;  Max change = ',E9.3)
C
C   Write the tension factors to disk (SGPRNT).
C
          CALL SGPRNT (N,LSIG,LIST,LPTR,LEND,SIGMA)
        ENDIF
C
C   Compute interpolated values on the uniform grid (UNIF).
C
        IFLGS = 0
        IF (VT) IFLGS = 1
        CALL UNIF (N,X,Y,Z,W,LIST,LPTR,LEND,IFLGS,SIGMA,NI,
     .             NI,NJ,PLAT,PLON,1, GRAD, FF,IER)
        IF (IER .LT. 0) THEN
          WRITE (LMSG,470) IER
          WRITE (*,470) IER
          STOP
        ELSE
          WRITE (LMSG,540) NEV, IER
          WRITE (*,540) NEV, IER
        ENDIF
  540   FORMAT (//5X,'UNIF:  Number of evaluation ',
     .          'points = ',I4/
     .          12X,'Number of extrapolation points =',I4)
      ELSEIF (MFLAG .EQ. 2) THEN
C
C *** C-1 interpolation (INTRC1) with global gradients
C     GRADG.
C
C   Initialize gradients GRAD to zeros.
C
        DO 90 K = 1,N
          GRAD(1,K) = 0.
          GRAD(2,K) = 0.
          GRAD(3,K) = 0.
   90     CONTINUE
        IF (.NOT. VT) ITGS = 1
C
C   Loop on GRADG/GETSIG iterations.
C
        IFLGS = 0
        DO 100 ITER = 1,ITGS
          NITG = MAXIT
          DGMX = DGMAX
          CALL GRADG (N,X,Y,Z,W,LIST,LPTR,LEND,IFLGS,
     .                SIGMA, NITG,DGMX,GRAD, IER)
          IF (IER .LT. 0) THEN
            WRITE (LMSG,450) IER
            WRITE (*,450) IER
            STOP
          ENDIF
          WRITE (LMSG,550) DGMAX, DGMX, MAXIT, NITG, IER
          WRITE (*,550) DGMAX, DGMX, MAXIT, NITG, IER
  550     FORMAT (//5X,'GRADG:  Tolerance = ',E9.3,', Max ',
     .            'change = ',E9.3/13X,'MAXIT = ',I2,
     .            ', No. iterations = ',I2,', IER = ',I2)
          IF (VT) THEN
C
C   Compute tension factors SIGMA (GETSIG).  IFLGS > 0 iff
C     VT = TRUE.
C
            IFLGS = 1
            CALL GETSIG (N,X,Y,Z,W,LIST,LPTR,LEND,GRAD,
     .                   TOL, SIGMA, DSM,IER)
            IF (IER .LT. 0) THEN
              WRITE (LMSG,460) IER
              WRITE (*,460) IER
              STOP
            ENDIF
            WRITE (LMSG,530) IER, DSM
            WRITE (*,530) IER, DSM
C
C   Write the tension factors to disk (SGPRNT).
C
            CALL SGPRNT (N,LSIG,LIST,LPTR,LEND,SIGMA)
          ENDIF
  100     CONTINUE
C
C   Compute interpolated values on the uniform grid (UNIF).
C
        CALL UNIF (N,X,Y,Z,W,LIST,LPTR,LEND,IFLGS,SIGMA,NI,
     .             NI,NJ,PLAT,PLON,1, GRAD, FF,IER)
        IF (IER .LT. 0) THEN
          WRITE (LMSG,470) IER
          WRITE (*,470) IER
          STOP
        ELSE
          WRITE (LMSG,540) NEV, IER
          WRITE (*,540) NEV, IER
        ENDIF
      ELSEIF (MFLAG .EQ. 3) THEN
C
C *** C-1 smoothing method SMSURF.
C
        WTK = 1./E
C
C   Store the weights WT.
C
        DO 110 K = 1,N
          WT(K) = WTK
  110     CONTINUE
C
C   Compute and print SMSURF parameters.
C
        IF (SM .LE. 0.) SM = REAL(N)
        SMTOL = SQRT(2./SM)
        GSTOL = .05*E
        WRITE (LMSG,560) E, SM, GSTOL, SMTOL, WTK
        WRITE (*,560) E, SM, GSTOL, SMTOL, WTK
  560   FORMAT (///5X,'SMSURF parameters:  Expected ',
     .          'squared error = ',E15.8/
     .          25X,'Smoothing parameter SM = ',E15.8/
     .          25X,'Gauss-Seidel tolerance = ',E15.8/
     .          28X,'Smoothing tolerance = ',E15.8/
     .          40X,'Weights = ',E15.8)
        IF (.NOT. VT) ITGS = 1
C
C   Loop on SMSURF/GETSIG iterations.
C
        IFLGS = 0
        DO 120 ITER = 1,ITGS
          CALL SMSURF (N,X,Y,Z,W,LIST,LPTR,LEND,IFLGS,SIGMA,
     .                 WT,SM,SMTOL,GSTOL,-1, F,GRAD,IER)
          IF (IER .LT. 0) THEN
            WRITE (LMSG,480) IER
            WRITE (*,480) IER
            STOP
          ENDIF
          IF (IER .EQ. 1) THEN
            WRITE (LMSG,570)
            WRITE (*,570)
  570       FORMAT (//5X,'Inactive constraint in SMSURF.  ',
     .              'F is a constant function.')
          ENDIF
          IF (VT) THEN
C
C   Compute tension factors SIGMA (GETSIG).  IFLGS > 0 iff
C     VT = TRUE.
C
            IFLGS = 1
            CALL GETSIG (N,X,Y,Z,F,LIST,LPTR,LEND,GRAD,
     .                   TOL, SIGMA, DSM,IER)
            IF (IER .LT. 0) THEN
              WRITE (LMSG,460) IER
              WRITE (*,460) IER
              STOP
            ENDIF
            WRITE (LMSG,530) IER, DSM
            WRITE (*,530) IER, DSM
C
C   Write the tension factors to disk (SGPRNT).
C
            CALL SGPRNT (N,LSIG,LIST,LPTR,LEND,SIGMA)
          ENDIF
  120     CONTINUE
C
C   Compute interpolated values on the uniform grid (UNIF).
C
        CALL UNIF (N,X,Y,Z,F,LIST,LPTR,LEND,IFLGS,SIGMA,NI,
     .             NI,NJ,PLAT,PLON,1, GRAD, FF,IER)
        IF (IER .LT. 0) THEN
          WRITE (LMSG,470) IER
          WRITE (*,470) IER
          STOP
        ELSE
          WRITE (LMSG,540) NEV, IER
          WRITE (*,540) NEV, IER
        ENDIF
      ENDIF
C
C *** Create the output data set.
C
      WRITE (LOUT,350) NI
      WRITE (LOUT,350) NJ
      WRITE (LOUT,360) (PLAT(I), I = 1,NI)
      WRITE (LOUT,360) (PLON(J), J = 1,NJ)
      WRITE (LOUT,360) ((FF(I,J), I = 1,NI), J = 1,NJ)
C
C *** Normal termination.
C
      WRITE (LMSG,580)
      WRITE (*,580)
  580 FORMAT (//5X,'No error encountered.'/)
      STOP
C
C Invalid input data.
C
  150 WRITE (LMSG,400) N, NMAX
      WRITE (*,400) N, NMAX
  400 FORMAT (//5X,'N = ',I5,' is outside its valid range.',
     .        '  NMAX = ',I4/)
      STOP
  160 WRITE (LMSG,410)
      WRITE (*,410)
  410 FORMAT (//5X,'Error encountered in reading data.'/)
      STOP
C
C Error message Formats:
C
  420 FORMAT (//5X,'Error in TRMESH:  IER =',I5/)
  430 FORMAT (//5X,'Error in INTRC0:  I = ',I3,', J = ',I3,
     .        ', IER = ',I2/)
  440 FORMAT (//5X,'Error in GRADL:  K =',I5,', IER = ',I2/)
  450 FORMAT (//5X,'Error in GRADG:  IER = ',I2/)
  460 FORMAT (//5X,'Error in GETSIG:  IER = ',I2/)
  470 FORMAT (//5X,'Error in UNIF:  IER = ',I2/)
  480 FORMAT (//5X,'Error in SMSURF:  IER = ',I2/)
      END
