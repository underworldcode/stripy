C--**--CH465--772--P:RW--21:7:1999
C--**--CH464--772--C:F--21:7:1999
C
C
C             STRITEST:  STRIPACK Test Program
C                         07/30/98
C
C
C   This driver tests software package STRIPACK for con-
C structing a Delaunay triangulation and Voronoi diagram
C of a set of nodes on the surface of the unit sphere.
C All STRIPACK subprograms are tested.
C
C   By default, a triangulation is created from a set of N
C nodes consisting of the north pole and N-1 points uniform-
C ly distributed around the 60-degree parallel (with
C constant longitudinal separation).  However, by enabling
C the READ statements below (C# in the first two columns),
C testing may be performed on an arbitrary set of up to NMAX
C nodes (refer to the PARAMETER statement below).  A data
C set consists of the following sequence of records:
C
C    N = Number of nodes (3 to NMAX):  format I5.
C    (RLAT(I),RLON(I), I = 1,N) = Nodal coordinates in
C                 degrees latitude (-90 to 90) and degrees
C                 longitude (-180 to 180):  format 2F15.9.
C
C This program must be linked to STRIPACK.
C
C
      CHARACTER*80 TITLET, TITLEV
      INTEGER IER, IFLAG, K, KSUM, KT, LIN, LOUT, LP, LPL,
     .        LPLT, LPLV, LW, LWK, LNEW, N, N0, N1, N2, N3,
     .        NA, NB, NCOL, NMAX, NN, NROW, NT, NT6, NTMX,
     .        NV
      INTEGER NEARND
      LOGICAL INSIDE, NUMBR
      REAL    A, AL, AREA, DEL, ELAT, ELON, P(3), PLTSIZ,
     .        SC, V1(3), V2(3), V3(3), VLAT, VLON, VNRM
      REAL    AREAS
C
      PARAMETER (NMAX=100, NTMX=2*NMAX, NT6=6*NMAX,
     .           LWK=2*NMAX, NCOL=NMAX, NROW=9)
C
C Array storage for the triangulation, work space, and nodal
C   coordinates.
C
      INTEGER LIST(NT6), LPTR(NT6), LEND(NMAX), IWK(LWK)
      REAL    DS(NMAX), RLAT(NMAX), RLON(NMAX),
     .        X(NMAX), Y(NMAX), Z(NMAX)
C
C Array storage for the Voronoi diagram:  adjacency array,
C   boundary triangle list, triangle circumcenters, and
C   circumradii.
C
      INTEGER LISTC(NT6), LBTRI(6,NCOL)
      REAL    XC(NTMX), YC(NTMX), ZC(NTMX), RC(NTMX)
C
C Array storage for the triangle list.
C
      INTEGER LTRI(NROW,NTMX)
      INTEGER I1MACH
C
C Plot size.
C
      DATA    PLTSIZ/7.5/
C
C Logical unit numbers for I/O:
C
      DATA    LPLT/3/,  LPLV/4/
      LIN = I1MACH(1)
      LOUT = I1MACH(2)
      OPEN (LOUT,FILE='res')
      OPEN (LPLT,FILE='res_3.eps')
      OPEN (LPLV,FILE='res_4.eps')
C
C Store plot titles  They must be enclosed in parentheses.
C
      TITLET = '(Triangulation created by STRITEST)'
      TITLEV = '(Voronoi diagram created by STRITEST)'
C
C Generate the default set of nodes as latitudinal and lon-
C   gitudinal coordinates.  DEL is the separation in degrees
C   between the nodes on the 60-degree line of latitude.
C
      N = 9
      RLAT(1) = 90.
      RLON(1) = 0.
      DEL = 360./REAL(N-1)
      DO 1 K = 2,N
        RLAT(K) = 60.
        RLON(K) = REAL(K-2)*DEL
    1   CONTINUE
C
C *** Read a data set:  N, RLAT, and RLON.
C
C#    OPEN (LIN,FILE='stritest.dat',STATUS='OLD')
C#    READ (LIN,100,ERR=30) N
      IF (N .LT. 3  .OR.  N .GT. NMAX) GO TO 31
C#    READ (LIN,110,ERR=30) (RLAT(K),RLON(K), K = 1,N)
C#100 FORMAT (I5)
C#110 FORMAT (2F15.9)
C
C Print a heading.
C
      WRITE (LOUT,300)
  300 FORMAT (///1X,25X,'STRIPACK test'///)
C
C Set X and Y to the values of RLON and RLAT, respectively,
C   in radians.  (RLON and RLAT are saved for printing by
C   Subroutine TRPRNT).
C
      SC = ATAN(1.)/45.
      DO 2 K = 1,N
        X(K) = SC*RLON(K)
        Y(K) = SC*RLAT(K)
    2   CONTINUE
C
C *** Transform spherical coordinates X and Y to Cartesian
C       coordinates (X,Y,Z) on the unit sphere (X**2 +
C       Y**2 + Z**2 = 1).
C
      CALL TRANS (N,Y,X, X,Y,Z)
C
C *** Create the triangulation and test the error flag.
C
      CALL TRMESH (N,X,Y,Z, LIST,LPTR,LEND,LNEW,IWK,
     .             IWK(N+1),DS,IER)
      IF (IER .EQ. -2) THEN
        WRITE (LOUT,510)
        STOP
      ELSEIF (IER .GT. 0) THEN
        WRITE (LOUT,515)
        STOP
      ENDIF
C
C *** Print the spherical coordinates and adjacency informa-
C       tion on LOUT.  IFLAG > 0 indicates that RLON and
C       RLAT only are to be printed.
C
      IFLAG = 1
      CALL TRPRNT (N,RLON,RLAT,Z,IFLAG,LIST,LPTR,LEND,LOUT)
C
C *** Test TRLIST and TRLPRT by creating and printing a
C                 triangle list.
C
      CALL TRLIST (N,LIST,LPTR,LEND,NROW, NT,LTRI,IER)
      CALL TRLPRT (N,RLON,RLAT,Z,IFLAG,NROW,NT,LTRI,LOUT)
C
C *** Test TRPLOT by plotting the portion of the triangula-
C                 tion contained in the hemisphere centered
C                 at E = (ELAT,ELON), where ELAT and ELON
C                 are taken to be the center of the range of
C                 the nodal latitudes and longitudes.
C
      ELAT = RLAT(1)
      VLAT = ELAT
      ELON = RLON(1)
      VLON = ELON
      DO 3 K = 2,N
        IF (RLAT(K) .LT. ELAT) ELAT = RLAT(K)
        IF (RLAT(K) .GT. VLAT) VLAT = RLAT(K)
        IF (RLON(K) .LT. ELON) ELON = RLON(K)
        IF (RLON(K) .GT. VLON) VLON = RLON(K)
    3   CONTINUE
      ELAT = (ELAT+VLAT)/2.0
      ELON = (ELON+VLON)/2.0
      A = 90.0
      NUMBR = N .LE. 200
      CALL TRPLOT (LPLT,PLTSIZ,ELAT,ELON,A,N,X,Y,Z,LIST,
     .             LPTR,LEND,TITLET,NUMBR, IER)
      IF (IER .EQ. 0) THEN
        WRITE (LOUT,305)
  305   FORMAT (/5X,'A triangulation plot file was ',
     .              'successfully created.'/)
      ELSE
        WRITE (LOUT,590) IER
      ENDIF
C
C *** Test AREAS by computing and printing the area of the
C                convex hull of the nodes (sum of triangle
C                areas) relative to the total surface area
C                (4*Pi).
C
      AREA = 0.
      DO 4 KT = 1,NT
        N1 = LTRI(1,KT)
        N2 = LTRI(2,KT)
        N3 = LTRI(3,KT)
        V1(1) = X(N1)
        V1(2) = Y(N1)
        V1(3) = Z(N1)
        V2(1) = X(N2)
        V2(2) = Y(N2)
        V2(3) = Z(N2)
        V3(1) = X(N3)
        V3(2) = Y(N3)
        V3(3) = Z(N3)
        AREA = AREA + AREAS(V1,V2,V3)
    4   CONTINUE
      AREA = AREA/(16.0*ATAN(1.0))
      WRITE (LOUT,310) AREA
  310 FORMAT (//5X,'Area of convex hull relative to total ',
     .             'surface area = ',F8.2)
C
C *** Test BNODES.  The ordered sequence of boundary nodes
C                   is stored in IWK.
C
      CALL BNODES (N,LIST,LPTR,LEND, IWK,NB,NA,NT)
      WRITE (LOUT,320) NB, NA, NT
  320 FORMAT (//5X,'Output from BNODES:  ',I4,
     .        ' boundary nodes,  ',I4,
     .        ' arcs,  ',I4,' triangles.'//)
C
C *** Test GETNP by ordering the nodes on distance from N0
C                and verifying the ordering.  The sequence
C                of nodal indexes is stored in IWK, and
C                the values of an increasing function (the
C                negative cosine) of angular distance is
C                stored in DS.
C
      N0 = N/2
      IWK(1) = N0
      DS(1) = -1.0
      KSUM = N0
      DO 5 K = 2,N
        CALL GETNP (X,Y,Z,LIST,LPTR,LEND,K, IWK, DS(K),IER)
        IF (IER .NE. 0  .OR.  DS(K) .LT. DS(K-1)) THEN
          WRITE (LOUT,520)
          STOP
        ENDIF
        KSUM = KSUM + IWK(K)
    5   CONTINUE
C
C   Test for all nodal indexes included in IWK.
C
      IF (KSUM .NE. (N*(N+1))/2) THEN
        WRITE (LOUT,520)
        STOP
      ENDIF
C
C *** Test NEARND by verifying that the nearest node to K is
C                 node K for K = 1 to N.
C
      DO 6 K = 1,N
        P(1) = X(K)
        P(2) = Y(K)
        P(3) = Z(K)
        N0 = NEARND (P,1,N,X,Y,Z,LIST,LPTR,LEND, AL)
        IF (N0 .NE. K  .OR.  AL .GT. 1.E-3) THEN
          WRITE (LOUT,530)
          STOP
        ENDIF
    6   CONTINUE
C
C *** Test DELARC by removing a boundary arc if possible.
C                 The last two nodes define a boundary arc
C                 in the default data set.
C
      N1 = N-1
      N2 = N
      CALL DELARC (N,N1,N2, LIST,LPTR,LEND,LNEW, IER)
      IF (IER .EQ. 1  .OR.  IER .EQ. 4) THEN
        WRITE (LOUT,540) IER
        STOP
      ENDIF
      IF (IER .NE. 0) THEN
        WRITE (LOUT,330) N1, N2
  330   FORMAT (5X,'Subroutine DELARC not tested:'/
     .          5X,'Nodes ',I4,' and ',I4,' do not form a ',
     .             'removable boundary arc.'//)
      ELSE
        CALL TRMESH (N,X,Y,Z, LIST,LPTR,LEND,LNEW,IWK,
     .               IWK(N+1),DS,IER)
      ENDIF
C
C *** Test CRLIST, VRPLOT, and SCOORD by constructing and
C                 plotting the Voronoi diagram, and printing
C                 the Voronoi region boundary (ordered
C                 sequence of Voronoi vertices) associated
C                 with N0.
C
C     Note that the triangulation data structure
C       is altered if NB > 0.
C
      CALL CRLIST (N,NCOL,X,Y,Z,LIST,LEND, LPTR,LNEW,
     .             LBTRI, LISTC,NB,XC,YC,ZC,RC,IER)
      IF (IER .NE. 0) THEN
        WRITE (LOUT,550) IER
        STOP
      ENDIF
C
C   Use the same parameter values that were used for the
C     triangulation plot (except the output unit and title).
C
      NT = 2*N-4
      CALL VRPLOT (LPLV,PLTSIZ,ELAT,ELON,A,N,X,Y,Z,NT,LISTC,
     .             LPTR,LEND,XC,YC,ZC,TITLEV,NUMBR, IER)
      IF (IER .EQ. 0) THEN
        WRITE (LOUT,335)
  335   FORMAT (/5X,'A Voronoi diagram plot file was ',
     .              'successfully created.'/)
      ELSE
        WRITE (LOUT,600) IER
      ENDIF
      N0 = 1
      WRITE (LOUT,340) N0
  340 FORMAT (17X,'Voronoi region for node ',I4//
     .        5X,'Triangle',5X,'Latitude',5X,'Longitude',
     .        5X,'Circumradius'/)
C
C   Initialize for loop on Voronoi vertices (triangle cir-
C     cumcenters).  The number of vertices is accumulated
C     in NV, and the vertex indexes are stored in IWK.  The
C     vertices are converted to latitude and longitude in
C     degrees for printing.
C
      NV = 0
      LPL = LEND(N0)
      LP = LPL
    7 LP = LPTR(LP)
        KT = LISTC(LP)
        NV = NV + 1
        IWK(NV) = KT
        CALL SCOORD (XC(KT),YC(KT),ZC(KT), VLAT,VLON,VNRM)
        VLAT = VLAT/SC
        VLON = VLON/SC
        WRITE (LOUT,345) KT, VLAT, VLON, RC(KT)
        IF (LP .NE. LPL) GO TO 7
  345 FORMAT (9X,I4,1X,F12.6,2X,F12.6,5X,F12.6)
C
C *** Test INSIDE by checking for node N0 inside its
C                 Voronoi region.
C
      P(1) = X(N0)
      P(2) = Y(N0)
      P(3) = Z(N0)
      IF (.NOT. INSIDE(P,NTMX,XC,YC,ZC,NV,IWK, IER))
     .  WRITE (LOUT,560) N0
      IF (IER .NE. 0) WRITE (LOUT,565) IER
C
C *** Recreate the triangulation and test the error flag.
C
      CALL TRMESH (N,X,Y,Z, LIST,LPTR,LEND,LNEW,IWK,
     .             IWK(N+1),DS,IER)
      IF (IER .EQ. -2) THEN
        WRITE (LOUT,510)
        STOP
      ELSEIF (IER .GT. 0) THEN
        WRITE (LOUT,515)
        STOP
      ENDIF
C
C *** Test EDGE by forcing an edge between nodes N1=1 and
C               N2=N.  LW is the number of columns reserved
C               for a 2 by LW work space array (IWK).
C
      N1 = 1
      N2 = N
      LW = LWK/2
      CALL EDGE (N1,N2,X,Y,Z, LW,IWK,LIST,LPTR,LEND, IER)
      IF (IER .NE. 0  .AND.  IER .NE. 5) THEN
        WRITE (LOUT,570) IER
        STOP
      ENDIF
C
C *** Test DELNOD by removing nodes 4 to N (in reverse
C                 order).  LW is the number of columns
C                 reserved for a 2 by LW work space array
C                 (IWK).
C
      IF (N .LE. 3) THEN
        WRITE (LOUT,350)
  350   FORMAT (//5X,'Subroutine DELNOD not tested:  ',
     .          'N cannot be reduced below 3.'//)
      ELSE
        NN = N
    8   LW = LWK/2
          CALL DELNOD (NN, NN,X,Y,Z,LIST,LPTR,LEND,LNEW,LW,
     .                 IWK, IER)
          IF (IER .NE. 0  .AND.  IER .NE. 5) THEN
            WRITE (LOUT,580) IER
            STOP
          ENDIF
          IF (NN .GT. 3) GO TO 8
      ENDIF
C
C Successful test.
C
      WRITE (LOUT,360)
  360 FORMAT (//5X,'No errors encountered.'/)
      STOP
C
C Error reading the data set.
C
C#   30 WRITE (*,500)
C#      STOP
C
C Invalid value of N.
C
   31 WRITE (*,505) N
      STOP
C
C Error message formats:
C
C#500 FORMAT (//5X,'*** Input data set invalid ***'/)
  505 FORMAT (//5X,'*** N is outside its valid ',
     .             'range:  N =',I5,' ***'/)
  510 FORMAT (//5X,'*** Error in TRMESH:  the first three ',
     .        'nodes are collinear ***'/)
  515 FORMAT (//5X,'*** Error in TRMESH:  duplicate nodes ',
     .        'encountered ***'/)
  520 FORMAT (//5X,'*** Error in GETNP ***'/)
  530 FORMAT (//5X,'*** Error in NEARND ***'/)
  540 FORMAT (//5X,'*** Error in DELARC:  IER = ',I1,
     .        ' ***'/)
  550 FORMAT (//5X,'*** Error in CRLIST:  IER = ',I1,
     .        ' ***'/)
  560 FORMAT (//5X,'*** Error in INSIDE:  node ',I4,' is ',
     .        'not contained in its Voronoi region ***'/)
  565 FORMAT (//5X,'*** Error in INSIDE:  IER = ',I1,
     .        ' ***'/)
  570 FORMAT (//5X,'*** Error in EDGE:  IER = ',I1,' ***'/)
  580 FORMAT (//5X,'*** Error in DELNOD:  IER = ',I1,
     .        ' ***'/)
  590 FORMAT (//5X,'*** Error in TRPLOT:  IER = ',I1,
     .        ' ***'/)
  600 FORMAT (//5X,'*** Error in VRPLOT:  IER = ',I1,
     .        ' ***'/)
      END
