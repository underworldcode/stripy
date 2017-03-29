      SUBROUTINE ARCINT (B,X1,X2,Y1,Y2,H1,H2,HX1,HX2,HY1,
     .                   HY2,SIGMA,DFLAG, HP,HXP,HYP,IER)
      INTEGER IER
      LOGICAL DFLAG
      REAL    B, X1, X2, Y1, Y2, H1, H2, HX1, HX2, HY1,
     .        HY2, SIGMA, HP, HXP, HYP
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/18/96
C
C   Gigen a line segment P1-P2 containing a point P =
C (XP,YP), along with function values and partial deriva-
C tives at the endpoints, this subroutine computes an
C interpolated value and, optionally, a gradient at P.  The
C value and tangential gradient component at P are taken to
C be the value and derivative of the Hermite interpolatory
C tension spline H defined by the endpoint values and tan-
C gential gradient components.  The normal gradient compo-
C nent at P is obtained by linear interpolation applied to
C the normal components at the endpoints.
C
C On input:
C
C       B = Local coordinate of P with respect to P1-P2:
C           P = B*P1 + (1-B)*P2.  Note that B may be comput-
C           ed from the coordinates of P as <P2-P1,P2-P>/
C           <P2-P1,P2-P1>.
C
C       X1,X2,Y1,Y2 = Coordinates of a pair of distinct
C                     points P1 and P2.
C
C       H1,H2 = Values of the interpolant H at P1 and P2,
C               respectively.
C
C       HX1,HX2,HY1,HY2 = x and y partial derivatives of H
C                         at P1 and P2.
C
C       SIGMA = Tension factor associated with P1-P2.
C
C       DFLAG = Logical flag which specifies whether first
C               partial derivatives at P are to be computed:
C               DFLAG = .TRUE. if and only if partials are
C               to be returned.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       HP = Interpolated value at P unless IER < 0, in
C            which case HP is not defined.
C
C       HXP,HYP = x and y partial derivatives at P unless
C                 DFLAG = FALSE or IER < 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if B < 0 or B > 1 and thus HP is an
C                     extrapolated value.
C             IER = -1 if P1 and P2 coincide.
C
C SRFPACK module required by ARCINT:  SNHCSH
C
C Intrinsic functions called by ARCINT:  ABS, EXP
C
C***********************************************************
C
      REAL B1, B2, CM, CM2, CMM, D1, D2, DS, DUMMY, DX, DY,
     .     E, E1, E2, EMS, GN, GT, S, S1, S2, SB1, SB2,
     .     SBIG, SIG, SINH2, SM, SM2, TM, TM1, TM2, TP1,
     .     TP2, TS
      DATA SBIG/85./
C
      DX = X2 - X1
      DY = Y2 - Y1
      DS = DX*DX + DY*DY
      IF (DS .EQ. 0.) GO TO 1
      IER = 0
C
C Compute local coordinates B1 and B2, tangential deriva-
C   tives S1 and S2, slope S, and second differences D1 and
C   D2.  S1, S2, S, D1, and D2 are scaled by the separation
C   D between P1 and P2.
C
      B1 = B
      B2 = 1. - B1
      IF (B1 .LT. 0.  .OR.  B2 .LT. 0.) IER = 1
      S1 = HX1*DX + HY1*DY
      S2 = HX2*DX + HY2*DY
      S = H2 - H1
      D1 = S - S1
      D2 = S2 - S
C
C Compute HP and, if required, the scaled tangential grad-
C   ient component GT.
C
      SIG = ABS(SIGMA)
      IF (SIG .LT. 1.E-9) THEN
C
C SIG = 0:  use Hermite cubic interpolation.
C
        HP = H1 + B2*(S1 + B2*(D1 + B1*(D1 - D2)))
        IF (.NOT. DFLAG) RETURN
        GT = S1 + B2*(D1 + D2 + 3.*B1*(D1 - D2))
      ELSEIF (SIG .LE. .5) THEN
C
C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
C   cancellation error in the hyperbolic functions.
C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        E = SIG*SM - CMM - CMM
        HP = H1 + B2*S1 + ((CM*SM2-SM*CM2)*(D1+D2) + SIG*
     .                     (CM*CM2-(SM+SIG)*SM2)*D1)/(SIG*E)
        IF (.NOT. DFLAG) RETURN
        SINH2 = SM2 + SB2
        GT = S1 + ((CM*CM2-SM*SINH2)*(D1+D2) + SIG*
     .             (CM*SINH2-(SM+SIG)*CM2)*D1)/E
      ELSE
C
C SIG > .5:  use negative exponentials in order to avoid
C   overflow.  Note that EMS = EXP(-SIG).  In the case of
C   extrapolation (negative B1 or B2), H is approximated
C   by a linear function if -SIG*B1 or -SIG*B2 is large.
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          HP = H1 + B2*S
          IF (.NOT. DFLAG) RETURN
          GT = S
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1. - EMS
          TS = TM*TM
          TM1 = 1. - E1
          TM2 = 1. - E2
          E = TM*(SIG*(1.+EMS) - TM - TM)
          HP = H1 + B2*S + (TM*TM1*TM2*(D1+D2) + SIG*
     .                      ((E2*TM1*TM1-B1*TS)*D1 +
     .                       (E1*TM2*TM2-B2*TS)*D2))/(SIG*E)
          IF (.NOT. DFLAG) RETURN
          TP1 = 1. + E1
          TP2 = 1. + E2
          GT = S + (TM1*(TM*TP2-SIG*E2*TP1)*D1 -
     .              TM2*(TM*TP1-SIG*E1*TP2)*D2)/E
        ENDIF
      ENDIF
C
C Compute the gradient at P, (HXP,HYP) = (GT/D)T + (GN/D)N,
C   where T = (DX,DY)/D (unit tangent vector), N = (-DY,DX)/
C   D (unit normal), and the scaled normal component is GN =
C   B1<(HX1,HY1),N> + B2<(HX2,HY2),N>.
C
      GN = B1*(HY1*DX-HX1*DY) + B2*(HY2*DX-HX2*DY)
      HXP = (GT*DX - GN*DY)/DS
      HYP = (GT*DY + GN*DX)/DS
      RETURN
C
C P1 and P2 coincide.
C
    1 IER = -1
      RETURN
      END
      SUBROUTINE CNTOUR (NX,NY,X,Y,Z,CVAL,LC,NCMAX,IWK, XC,
     .                   YC,ILC,NC,IER)
      INTEGER NX, NY, LC, NCMAX, IWK(NX,*), ILC(NCMAX), NC,
     .        IER
      REAL    X(NX), Y(NY), Z(NX,NY), CVAL, XC(LC), YC(LC)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   04/28/90
C
C   Given a set of function values Z = F(X,Y) at the verti-
C ces of an NX by NY rectangular grid, this subroutine de-
C termines a set of contour lines associated with F = CVAL.
C A contour line is specified by an ordered sequence of
C points (XC,YC), each lying on a grid edge and computed
C from the linear interpolant of the function values at the
C endpoints of the edge.  The accuracy of the contour lines
C is thus directly related to the number of grid points.  If
C a contour line forms a closed curve, the first point coin-
C cides with the last point.  Otherwise, the first and last
C points lie on the grid boundary.
C
C   Note that the problem is ill-conditioned in the vicinity
C of a double zero of F-CVAL.  Thus, if a grid cell is
C crossed by two contour lines (all four sides intersected),
C three different configurations are possible, corresponding
C to a local minimum, a local maximum, or a saddle point.
C It is arbitrarily assumed in this case that the contour
C lines intersect, representing a saddle point.  Also, in
C order to treat the case of F = CVAL at a vertex in a con-
C sistent manner, this case is always treated as F > CVAL.
C Hence, if F takes on the same value at both ends of an
C edge, it is assumed that no contour line intersects that
C edge.  In particular, a constant function, including
C F = CVAL, results in no contour lines.
C
C On input:
C
C       NX = Number of grid points in the x direction.
C            NX .GE. 2.
C
C       NY = Number of grid points in the y direction.
C            NY .GE. 2.
C
C       X = Array of length NX containing a strictly in-
C           creasing sequence of values.
C
C       Y = Array of length NY containing a strictly in-
C           creasing sequence of values.
C
C       Z = Array of function values at the vertices of the
C           rectangular grid.  Z(I,J) = F(X(I),Y(J)) for
C           I = 1,...,NX and J = 1,...,NY.
C
C       CVAL = Constant function value defining a contour
C              line as the set of points (X,Y) such that
C              F(X,Y) = CVAL.
C
C       LC = Length of arrays XC and YC, and maximum allow-
C            able number of points defining contour lines.
C            LC = 2(NX-1)(NY-1) + (NX*NY+1)/2 is (probably
C            more than) sufficient.  LC .GE. 2.
C
C       NCMAX = Length of array ILC, and maximum allowable
C               number of contour lines.  NCMAX = (NX*NY+1)/
C               2 is sufficient.  NCMAX .GE. 1.
C
C The above parameters are not altered by this routine.
C
C       IWK = Integer array of length .GE. NX*(NY-1) to be
C             used as work space.
C
C       XC,YC = Arrays of length LC.
C
C       ILC = Integer array of length NCMAX.
C
C On output:
C
C       XC,YC = Arrays containing the coordinates of NC con-
C               tour lines.  For K = 1,...,NC, contour line
C               K is defined by the sequence of points with
C               indexes ILC(K-1)+1,...,ILC(K) where ILC(0) =
C               0.
C
C       ILC = Array containing the indexes (to XC and YC)
C             associated with the terminal point of contour
C             line K in position K for K = 1,...,NC (if NC
C             .GT. 0).
C
C       NC = Number of contour lines whose points are stored
C            in XC and YC.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and all
C                     contour lines were found.
C             IER = 1 if NX, NY, LC, or NCMAX is outside its
C                     valid range.  NC = 0 and XC, YC, and
C                     ILC are not altered in this case.
C             IER = 2 if X or Y is not strictly increasing.
C                     NC = 0 and XC, YC, and ILC are not
C                     altered in this case.
C             IER = K for K > LC, where K is the required
C                     length of XC and YC, if more storage
C                     space is required to complete the
C                     specification of contour line NC and/
C                     or additional contour lines up to a
C                     total of NCMAX.  NC .GE. 1 and ILC(NC)
C                     = LC in this case.
C            IER = -1 if more than NCMAX contour lines are
C                     present (more space is required in
C                     ILC).  NC = NCMAX, and LC may or may
C                     not be sufficient for the additional
C                     contour lines in this case.  (This is
C                     not determined.)
C
C   In the unlikely event of an internal failure, a message
C is printed on logical unit LUN (specified in the DATA
C statement below).  IER may be 0 in this case.
C
C Modules required by CNTOUR:  None
C
C***********************************************************
C
      INTEGER I, I1, I2, IB, IN, IND, ISID, ISIDB, ISIDN,
     .        J, J1, J2, JB, JN, K, LCON, LMX, LUN, NCMX,
     .        NCON, NI, NIM1, NJ, NJM1
      LOGICAL BDRY
      REAL    CV, W, XF, XN, XP, YF, YN, YP, Z1, Z2
      DATA    LUN/0/
C
C Store parameters in local variables.
C
      NI = NX
      NJ = NY
      NIM1 = NI - 1
      NJM1 = NJ - 1
      CV = CVAL
      LMX = LC
      NCMX = NCMAX
      NC = 0
C
C Test for invalid input parameters.
C
      IER = 1
      IF (NI .LT. 2  .OR.  NJ .LT. 2  .OR.  LMX .LT. 2  .OR.
     .    NCMX .LT. 1) RETURN
C
C Test for nonincreasing values of X or Y.
C
      IER = 2
      DO 1 I = 2,NI
        IF (X(I) .LE. X(I-1)) RETURN
    1   CONTINUE
      DO 2 J = 2,NJ
        IF (Y(J) .LE. Y(J-1)) RETURN
    2   CONTINUE
C
C Loop on grid cells, initializing edge indicators (stored
C   in IWK) to zeros.  For each cell, the indicator IND is a
C   4-bit integer with each bit corresponding to an edge of
C   the cell, and having value 1 iff the edge has been pro-
C   cessed.  Note that two IND values must be adjusted when
C   an interior edge is processed.  The cell sides (edges)
C   are numbered (1,2,4,8) in counterclockwise order start-
C   ing from the bottom.  This corresponds to an ordering of
C   the weighted IND bits from low order to high order.
C   Grid cells are identified with their lower left corners.
C
      DO 4 J = 1,NJM1
        DO 3 I = 1,NIM1
          IWK(I,J) = 0
    3     CONTINUE
    4   CONTINUE
C
C First determine open contours by looping on boundary edges
C   in counterclockwise order starting from the lower left.
C   For each unprocessed boundary edge intersected by a con-
C   tour line, the contour line is determined and IWK is up-
C   dated to reflect the edges intersected.  The boundary
C   cell (lower left corner) is indexed by (IB,JB) and the
C   boundary edge is specified by ISIDB.  NCON and LCON are
C   local variables containing the number of contour lines
C   encountered and the current length of XC and YC.
C
      NCON = 0
      LCON = 0
      ISIDB = 1
      IB = 1
      JB = 1
C
C Top of loop on boundary edges.  The edge has been
C   processed iff IND/ISIDB is odd.
C
    5 IND = IWK(IB,JB)
      IF (IND/ISIDB .NE. 2*((IND/ISIDB)/2)) GO TO 9
C
C Update the edge indicator and store the vertex indexes of
C   the endpoints of the edge.
C
      IWK(IB,JB) = IND + ISIDB
      IF (ISIDB .EQ. 1) THEN
        I1 = IB
        J1 = JB
        I2 = IB + 1
        J2 = JB
      ELSEIF (ISIDB .EQ. 2) THEN
        I1 = IB + 1
        J1 = JB
        I2 = IB + 1
        J2 = JB + 1
      ELSEIF (ISIDB .EQ. 4) THEN
        I1 = IB + 1
        J1 = JB + 1
        I2 = IB
        J2 = JB + 1
      ELSE
        I1 = IB
        J1 = JB + 1
        I2 = IB
        J2 = JB
      ENDIF
C
C Proceed to the next edge if there is no intersection.
C
      Z1 = Z(I1,J1)
      Z2 = Z(I2,J2)
      IF ((Z1 .LT. CV  .AND.  Z2 .LT. CV)  .OR.
     .    (Z1 .GE. CV  .AND.  Z2 .GE. CV)) GO TO 9
C
C Store the zero of the linear interpolant of Z1-CV and
C   Z2-CV as the first point of an open contour unless
C   NCMAX contour lines have been found or there is in-
C   sufficient space reserved for XC and YC.
C
      IF (NCON .EQ. NCMX) THEN
        IER = -1
        GO TO 16
      ENDIF
      NCON = NCON + 1
      LCON = LCON + 1
      W = (CV-Z1)/(Z2-Z1)
      XP = X(I1) + W*(X(I2)-X(I1))
      YP = Y(J1) + W*(Y(J2)-Y(J1))
      IF (LCON .LE. LMX) THEN
        XC(LCON) = XP
        YC(LCON) = YP
      ENDIF
C
C Initialize for loop on cells intersected by the open
C   contour line.
C
      I = IB
      J = JB
      ISID = ISIDB
C
C Traverse the contour line.  Cell (I,J) was entered on side
C   ISID = (I1,J1)->(I2,J2).  Find an exit edge E (unproces-
C   sed edge intersected by the contour) by looping on the
C   remaining three sides, starting with the side opposite
C   ISID.
C
    6 IND = IWK(I,J)
      DO 7 K = 1,3
        ISID = 2*ISID
        IF (K .NE. 2) ISID = 2*ISID
        IF (ISID .GT. 15) ISID = ISID/16
        IF (ISID .EQ. 1) THEN
          I1 = I
          J1 = J
          I2 = I + 1
          J2 = J
        ELSEIF (ISID .EQ. 2) THEN
          I1 = I + 1
          J1 = J
          I2 = I + 1
          J2 = J + 1
        ELSEIF (ISID .EQ. 4) THEN
          I1 = I + 1
          J1 = J + 1
          I2 = I
          J2 = J + 1
        ELSE
          I1 = I
          J1 = J + 1
          I2 = I
          J2 = J
        ENDIF
C
C Test for a 1 in bit position ISID of cell (I,J) and bypass
C   the edge if it has been previously encountered.
C
        IF (IND/ISID .NE. 2*((IND/ISID)/2)) GO TO 7
C
C Update IWK for edge E = (I1,J1)->(I2,J2).  (IN,JN) indexes
C   the cell which shares E with cell (I,J), and ISIDN is
C   the side number of E in (IN,JN).  BDRY is true iff E is
C   a boundary edge (with no neighboring cell).
C
        IWK(I,J) = IWK(I,J) + ISID
        IF (ISID .LE. 2) THEN
          IN = I1
          JN = J2 - 1
          ISIDN = 4*ISID
        ELSE
          IN = I1 - 1
          JN = J2
          ISIDN = ISID/4
        ENDIF
        BDRY = IN .EQ. 0  .OR.  IN .EQ. NI  .OR.
     .         JN .EQ. 0  .OR.  JN .EQ. NJ
        IF (.NOT. BDRY) IWK(IN,JN) = IWK(IN,JN) + ISIDN
C
C Exit the loop on sides if E is intersected by the contour.
C
        Z1 = Z(I1,J1)
        Z2 = Z(I2,J2)
        IF ((Z1 .LT. CV  .AND.  Z2 .GE. CV)  .OR.
     .      (Z1 .GE. CV  .AND.  Z2 .LT. CV)) GO TO 8
    7   CONTINUE
C*
C Error -- No exit point found.  Print a message and exit
C          the contour traversal loop.
C
      WRITE (LUN,100) NCON
  100 FORMAT (///5X,'Error in CNTOUR:  Contour line L ',
     .        'begins on the boundary'/5X,'and terminates ',
     .        'in the interior for L =',I4/)
      ILC(NCON) = LCON
      GO TO 9
C*
C Add the intersection point (XN,YN) to the list unless it
C   coincides with the previous point (XP,YP) or there is
C   not enough space in XC and YC.
C
    8 W = (CV-Z1)/(Z2-Z1)
      XN = X(I1) + W*(X(I2)-X(I1))
      YN = Y(J1) + W*(Y(J2)-Y(J1))
      IF (XN .NE. XP  .OR.  YN .NE. YP) THEN
        LCON = LCON + 1
        XP = XN
        YP = YN
        IF (LCON .LE. LMX) THEN
          XC(LCON) = XN
          YC(LCON) = YN
        ENDIF
      ENDIF
C
C Bottom of contour traversal loop.  If E is not a boundary
C   edge, reverse the edge direction (endpoint indexes) and
C   update the cell index and side number.
C
      IF (.NOT. BDRY) THEN
        I = I1
        J = J1
        I1 = I2
        J1 = J2
        I2 = I
        J2 = J
        I = IN
        J = JN
        ISID = ISIDN
        GO TO 6
      ENDIF
C
C Update ILC with a pointer to the end of the contour line.
C
      ILC(NCON) = LCON
C
C Bottom of loop on boundary edges.  Update the boundary
C   cell index and side number, and test for termination.
C
    9 IF (ISIDB .EQ. 1) THEN
        IF (IB .LT. NIM1) THEN
          IB = IB + 1
        ELSE
          ISIDB = 2
        ENDIF
      ELSEIF (ISIDB .EQ. 2) THEN
        IF (JB .LT. NJM1) THEN
          JB = JB + 1
        ELSE
          ISIDB = 4
        ENDIF
      ELSEIF (ISIDB .EQ. 4) THEN
        IF (IB .GT. 1) THEN
          IB = IB - 1
        ELSE
          ISIDB = 8
        ENDIF
      ELSE
        IF (JB .GT. 1) THEN
          JB = JB - 1
        ELSE
          ISIDB = 16
        ENDIF
      ENDIF
      IF (ISIDB .LT. 16) GO TO 5
C
C Determine closed contours by looping on interior edges --
C   the first two sides (bottom and right) of each cell,
C   excluding boundary edges.  The beginning cell is indexed
C   by (IB,JB), and the beginning side number is ISIDB.
C
      DO 15 JB = 1,NJM1
      DO 14 IB = 1,NIM1
      DO 13 ISIDB = 1,2
        IF (JB .EQ. 1  .AND.  ISIDB .EQ. 1) GO TO 13
        IF (IB .EQ. NIM1  .AND.  ISIDB .EQ. 2) GO TO 13
C
C Bypass the edge if it was previously encountered
C   (IND/ISIDB odd).
C
        IND = IWK(IB,JB)
        IF (IND/ISIDB .NE. 2*((IND/ISIDB)/2)) GO TO 13
C
C Determine the endpoint indexes of the beginning edge E =
C   (I1,J1)->(I2,J2), find the index (I,J) and side number
C   ISID of the cell which shares E with (IB,JB), and up-
C   date IWK.
C
        IF (ISIDB .EQ. 1) THEN
          I1 = IB
          J1 = JB
          I2 = IB + 1
          J2 = JB
          I = IB
          J = JB - 1
          ISID = 4
        ELSE
          I1 = IB + 1
          J1 = JB
          I2 = IB + 1
          J2 = JB + 1
          I = I1
          J = J1
          ISID = 8
        ENDIF
        IWK(IB,JB) = IND + ISIDB
        IWK(I,J) = IWK(I,J) + ISID
C
C Proceed to the next interior edge if there is no
C   intersection.
C
        Z1 = Z(I1,J1)
        Z2 = Z(I2,J2)
        IF ((Z1 .LT. CV  .AND.  Z2 .LT. CV)  .OR.
     .      (Z1 .GE. CV  .AND.  Z2 .GE. CV)) GO TO 13
C
C Store the intersection point as the first point of a
C   closed contour unless NCMAX contour lines have been
C   found or there is insufficient space in XC and YC.
C
        IF (NCON .EQ. NCMX) THEN
          IER = -1
          GO TO 16
        ENDIF
        NCON = NCON + 1
        LCON = LCON + 1
        W = (CV-Z1)/(Z2-Z1)
        XP = X(I1) + W*(X(I2)-X(I1))
        YP = Y(J1) + W*(Y(J2)-Y(J1))
        IF (LCON .LE. LMX) THEN
          XC(LCON) = XP
          YC(LCON) = YP
        ENDIF
        XF = XP
        YF = YP
C
C Traverse the contour line.  Cell (I,J) was entered on side
C   ISID = edge (I2,J2)->(I1,J1).  Reverse the edge direc-
C   tion.
C
   10   IN = I1
        JN = J1
        I1 = I2
        J1 = J2
        I2 = IN
        J2 = JN
        IND = IWK(I,J)
C
C Find an exit edge E by looping on the remaining three
C   sides, starting with the side opposite ISID.
C
        DO 11 K = 1,3
          ISID = 2*ISID
          IF (K .NE. 2) ISID = 2*ISID
          IF (ISID .GT. 15) ISID = ISID/16
          IF (ISID .EQ. 1) THEN
            I1 = I
            J1 = J
            I2 = I + 1
            J2 = J
          ELSEIF (ISID .EQ. 2) THEN
            I1 = I + 1
            J1 = J
            I2 = I + 1
            J2 = J + 1
          ELSEIF (ISID .EQ. 4) THEN
            I1 = I + 1
            J1 = J + 1
            I2 = I
            J2 = J + 1
          ELSE
            I1 = I
            J1 = J + 1
            I2 = I
            J2 = J
          ENDIF
C
C Bypass the edge if it has been previously encountered.
C
          IF (IND/ISID .NE. 2*((IND/ISID)/2)) GO TO 11
C
C Determine the index (IN,JN) and side number ISIDN of the
C   cell which shares edge E = (I1,J1)->(I2,J2) with cell
C   (I,J), and update IWK.
C
          IF (ISID .LE. 2) THEN
            IN = I1
            JN = J2 - 1
            ISIDN = 4*ISID
          ELSE
            IN = I1 - 1
            JN = J2
            ISIDN = ISID/4
          ENDIF
          IWK(I,J) = IWK(I,J) + ISID
          IWK(IN,JN) = IWK(IN,JN) + ISIDN
C
C Exit the loop on sides if E is intersected.
C
          Z1 = Z(I1,J1)
          Z2 = Z(I2,J2)
          IF ((Z1 .LT. CV  .AND.  Z2 .GE. CV)  .OR.
     .        (Z1 .GE. CV  .AND.  Z2 .LT. CV)) GO TO 12
   11     CONTINUE
C*
C Error -- No exit point found.  Print a message and exit
C          the contour traversal loop.
C
        WRITE (LUN,110) NCON
  110   FORMAT (///5X,'Error in CNTOUR:  Contour line L ',
     .          'is open but'/5X,'does not intersect the ',
     .          'boundary for L =',I4/)
        ILC(NCON) = LCON
        GO TO 13
C*
C Add the intersection point to the list unless it coincides
C   with the previous point or there is not enough space in
C   XC and YC.
C
   12   W = (CV-Z1)/(Z2-Z1)
        XN = X(I1) + W*(X(I2)-X(I1))
        YN = Y(J1) + W*(Y(J2)-Y(J1))
        IF (XN .NE. XP  .OR.  YN .NE. YP) THEN
          LCON = LCON + 1
          XP = XN
          YP = YN
          IF (LCON .LE. LMX) THEN
            XC(LCON) = XN
            YC(LCON) = YN
          ENDIF
        ENDIF
C
C Bottom of contour traversal loop.  If the next cell is not
C   the beginning cell, update the cell index and side num-
C   ber.
C
        IF (IN .NE. IB  .OR.  JN .NE. JB) THEN
          I = IN
          J = JN
          ISID = ISIDN
          GO TO 10
        ENDIF
C
C Add the first point as the last point (unless the first
C   and last points already coincide), and update ILC.
C
        IF (XP .NE. XF  .OR.  YP .NE. YF) THEN
          LCON = LCON + 1
          IF (LCON .LE. LMX) THEN
            XC(LCON) = XF
            YC(LCON) = YF
          ENDIF
        ENDIF
        ILC(NCON) = LCON
C
C Bottom of loop on interior edges.
C
   13   CONTINUE
   14   CONTINUE
   15   CONTINUE
      IER = 0
C
C Test for insufficient storage reserved for XC and YC.
C
   16 IF (LCON .GT. LMX) IER = LCON
      NC = NCON
      RETURN
      END
      SUBROUTINE COORDS (XP,YP,X1,X2,X3,Y1,Y2,Y3, B1,B2,
     .                   B3, IER)
      INTEGER IER
      REAL    XP, YP, X1, X2, X3, Y1, Y2, Y3, B1, B2, B3
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This subroutine computes the barycentric (areal) coordi-
C nates B1, B2, and B3 of a point P with respect to the tri-
C angle with vertices P1, P2, and P3:  the solution to the
C linear system defined by B1 + B2 + B3 = 1 and B1*P1 +
C B2*P2 + B3*P3 = P.  Note that B1 is a linear function
C of P which satisfies B1 = 1 at P = P1 and B1 = 0 on
C the triangle side P2-P3.  Also, B1 < 0 if and only if
C P is to the right of P2->P3 (and thus exterior to the
C triangle).  B2 and B3 satisfy similar properties.
C
C On input:
C
C       XP,YP = Cartesian coordinates of P.
C
C       X1,X2,X3,Y1,Y2,Y3 = Coordinates of the vertices of
C                           the triangle P1, P2, and P3.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       B1,B2,B3 = Barycentric coordinates unless IER = 1,
C                  in which case the coordinates are not
C                  defined.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if the vertices of the triangle are
C                     collinear.
C
C Modules required by COORDS:  None
C
C***********************************************************
C
      REAL A, PX, PY, XP1, XP2, XP3, YP1, YP2, YP3
C
      PX = XP
      PY = YP
C
C Compute components of the vectors P->P1, P->P2, and P->P3.
C
      XP1 = X1 - PX
      YP1 = Y1 - PY
      XP2 = X2 - PX
      YP2 = Y2 - PY
      XP3 = X3 - PX
      YP3 = Y3 - PY
C
C Compute subtriangle areas B1 = P->P2 X P->P3, B2 = P->P3 X
C   P->P1, and B3 = P->P1 X P->P2.
C
      B1 = XP2*YP3 - XP3*YP2
      B2 = XP3*YP1 - XP1*YP3
      B3 = XP1*YP2 - XP2*YP1
C
C Compute twice the signed area of the triangle.
C
      A = B1 + B2 + B3
      IF (A .EQ. 0.) GO TO 1
C
C Normalize the coordinates.
C
      B1 = B1/A
      B2 = B2/A
      B3 = B3/A
      IER = 0
      RETURN
C
C The vertices are collinear.
C
    1 IER = -1
      RETURN
      END
      SUBROUTINE CRPLOT (LUN,PLTSIZ,NX,NY,PX,PY,PZ,NCON,IWK,
     .                   XC,YC, IER)
      INTEGER LUN, NX, NY, NCON, IWK(*), IER
      REAL    PLTSIZ, PX(NX), PY(NY), PZ(NX,NY),
     .        XC(*), YC(*)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   04/12/97
C
C   Given a set of function values PZ = F(X,Y) at the ver-
C tices of an NX by NY rectangular grid, this subroutine
C creates a level-2 Encapsulated PostScript (EPS) file
C containing a contour plot of the piecewise bilinear inter-
C polant of the function values.
C
C   The accuracy of the contour lines increases with the
C number of grid points.  Refer to Subroutine CNTOUR for
C further details.
C
C
C On input:
C
C       LUN = Logical unit number in the range 0 to 99.
C             The unit should be opened with an appropriate
C             file name before the call to this routine.
C
C       PLTSIZ = Plot size in inches.  A window containing
C                the plot is mapped, with aspect ratio
C                preserved, to a rectangular viewport with
C                maximum side-length PLTSIZ.  The viewport
C                is centered on the 8.5 by 11 inch page, and
C                its boundary is drawn.  1.0 .LE. PLTSIZ
C                .LE. 7.5.
C
C       NX = Number of grid points in the x direction.
C            NX .GE. 2.
C
C       NY = Number of grid points in the y direction.
C            NY .GE. 2.
C
C       PX = Array of length NX containing a strictly in-
C            creasing sequence of values.
C
C       PY = Array of length NY containing a strictly in-
C            creasing sequence of values.
C
C       PZ = Array of function values at the vertices of the
C            rectangular grid.  PZ(I,J) = F(PX(I),PY(J)) for
C            I = 1,...,NX and J = 1,...,NY.
C
C       NCON = Number of contour values.  The contour values
C              are uniformly distributed over the range of
C              PZ values.  NCON .GE. 1.
C
C The above parameters are not altered by this routine.
C
C       IWK = Integer array of length at least 1.5*NX*NY to
C             be used as work space.
C
C       XC,YC = Real arrays of length at least 2.5*NX*NY to
C               be used as work space.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LUN, PLTSIZ, NX, NY, or NCON is
C                     outside its valid range.
C             IER = 2 if PX or PY is not strictly
C                     increasing.
C             IER = 3 if the range of PZ values has zero
C                     width (F is constant).
C             IER = 4 if an error was encountered in writing
C                     to unit LUN.
C             IER = 5 if an unexpected error flag was re-
C                     turned by Subroutine CNTOUR.  This
C                     should not occur.
C
C   In the unlikely event of an internal failure, a message
C is printed on the standard output device.  IER may be 0
C in this case.
C
C Module required by CRPLOT:  CNTOUR
C
C Intrinsic functions called by CRPLOT:  CHAR, REAL
C
C***********************************************************
C
      INTEGER I, IC, IERR, IH, IPX1, IPX2, IPY1, IPY2, IW,
     .        J, K, KV, LC, NC, NCMAX
      REAL    CVAL, DX, DY, DZ, PZIJ, R, SFX, SFY, T,
     .        TX, TY, ZMAX, ZMIN
C
C Local parameters:
C
C CVAL =      Contour value between ZMIN and ZMAX
C DX =        Window width PX(NX)-PX(1)N
C DY =        Window height PY(NY)-PY(1)
C DZ =        Interval between contour values:
C               (ZMAX-ZMIN)/(NCON+1)
C I,J =       Row and column indexes for PZ
C IC =        Index (for IWK) of a contour line associated
C               with contour value CVAL:  1 to NC
C IERR =      Error flag for calls to CNTOUR
C IH =        Height of the bounding box (viewport) in
C               points
C IPX1,IPY1 = X and y coordinates (in points) of the lower
C               left corner of the bounding box
C IPX2,IPY2 = X and y coordinates (in points) of the upper
C               right corner of the bounding box
C IW =        Width of the bounding box in points
C K =         Index (for XC and YC) of a point on a contour
C               line
C KV =        DO-loop index for loop on contour values
C LC =        Length of arrays XC and YC
C NC =        Number of contour lines associated with
C               contour value CVAL
C NCMAX =     Maximum allowable value of NC
C PZIJ =      PZ(I,J)
C R =         Aspect ratio DX/DY
C SFX,SFY =   Scale factors for mapping window coordinates
C               to viewport coordinates
C T =         Temporary variable
C TX,TY =     Translation vector for mapping window coordi-
C               nates to viewport coordinates
C ZMIN,ZMAX = Minimum and maximum of the PZ values
C
C
C Test for error 1.
C
      IF (LUN .LT. 0  .OR.  LUN .GT. 99  .OR.
     .    PLTSIZ .LT. 1.0  .OR.  PLTSIZ .GT. 7.5  .OR.
     .    NX .LT. 2  .OR.  NY .LT. 2  .OR.  NCON .LT. 1)
     .  GO TO 11
C
C Compute the aspect ratio of the window.
C
      DX = PX(NX) - PX(1)
      DY = PY(NY) - PY(1)
      IF (DX .EQ. 0.0  .OR.  DY .EQ. 0.0) GO TO 12
      R = DX/DY
C
C Compute the range of PZ values and the interval between
C   contour values.
C
      ZMIN = PZ(1,1)
      ZMAX = ZMIN
      DO 2 J = 1,NY
        DO 1 I = 1,NX
          PZIJ = PZ(I,J)
          IF (PZIJ .LT. ZMIN) ZMIN = PZIJ
          IF (PZIJ .GT. ZMAX) ZMAX = PZIJ
    1     CONTINUE
    2   CONTINUE
      DZ = (ZMAX-ZMIN)/REAL(NCON+1)
      IF (DZ .LE. 0.0) GO TO 13
C
C Compute the lower left (IPX1,IPY1) and upper right
C   (IPX2,IPY2) corner coordinates of the bounding box
C   (the viewport).  The coordinates, specified in default
C   user space units (points, at 72 points/inch with origin
C   at the lower left corner of the page), are chosen to
C   preserve the aspect ratio R, and to center the plot on
C   the 8.5 by 11 inch page.  The center of the page is
C   (306,396), and T = PLTSIZ/2 in points.
C
      T = 36.0*PLTSIZ
      IF (R .GE. 1.0) THEN
        IPX1 = 306 - NINT(T)
        IPX2 = 306 + NINT(T)
        IPY1 = 396 - NINT(T/R)
        IPY2 = 396 + NINT(T/R)
      ELSE
        IPX1 = 306 - NINT(T*R)
        IPX2 = 306 + NINT(T*R)
        IPY1 = 396 - NINT(T)
        IPY2 = 396 + NINT(T)
      ENDIF
C
C Output header comments.
C
      WRITE (LUN,100,ERR=14) IPX1, IPY1, IPX2, IPY2
  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/
     .        '%%BoundingBox:',4I4/
     .        '%%Title:  Contour Plot'/
     .        '%%Creator:  SRFPACK'/
     .        '%%EndComments')
C
C Draw the bounding box.
C
      WRITE (LUN,110,ERR=14) IPX1, IPY1
      WRITE (LUN,120,ERR=14) IPX1, IPY2
      WRITE (LUN,120,ERR=14) IPX2, IPY2
      WRITE (LUN,120,ERR=14) IPX2, IPY1
      WRITE (LUN,130,ERR=14)
      WRITE (LUN,140,ERR=14)
  110 FORMAT (2I4,' moveto')
  120 FORMAT (2I4,' lineto')
  130 FORMAT ('closepath')
  140 FORMAT ('stroke')
C
C Set up a mapping from the window to the viewport.
C
      IW = IPX2 - IPX1
      IH = IPY2 - IPY1
      SFX = REAL(IW)/DX
      SFY = REAL(IH)/DY
      TX = IPX1 - SFX*PX(1)
      TY = IPY1 - SFY*PY(1)
      WRITE (LUN,150,ERR=14) TX, TY, SFX, SFY
  150 FORMAT (2F12.6,' translate'/
     .        2F12.6,' scale')
C
C Set the line thickness to 2 points.  (Since the scale
C   factors are applied to everything, the width must be
C   specified in world coordinates.)
C
      T = 4.0/(SFX+SFY)
      WRITE (LUN,160,ERR=14) T
  160 FORMAT (F12.6,' setlinewidth')
C
C Compute parameters for CNTOUR:
C
C   NCMAX = Maximum allowable number of contour lines
C           associated with each contour value.
C   LC = Length of arrays XC and YC and maximum allowable
C        number of points defining all the contour lines
C        associated with a contour value.
C
      NCMAX = (NX*NY+1)/2
      LC = 2*(NX-1)*(NY-1) + NCMAX
C
C Loop on contour values CVAL uniformly spaced in the open
C   interval (ZMIN,ZMAX).
C
      CVAL = ZMIN
      DO 5 KV = 1,NCON
        CVAL = CVAL + DZ
C
C Compute a sequence of NC contour lines associated with
C   F = CVAL.  For IC = 1 to NC, IWK(IC) is the index (for
C   XC and YC) of the last point of contour IC.
C
        CALL CNTOUR (NX,NY,PX,PY,PZ,CVAL,LC,NCMAX,
     .               IWK(NCMAX+1), XC,YC,IWK,NC,IERR)
        IF (IERR .EQ. 2) GO TO 12
        IF (IERR .NE. 0) GO TO 15
C
C Draw the NC contours.
C
        IC = 0
        K = 0
    3   IC = IC + 1
          K = K + 1
C
C   Create a path consisting of contour IC.
C
          WRITE (LUN,170,ERR=14) XC(K), YC(K)
  170     FORMAT (2F12.6,' moveto')
    4     K = K + 1
            WRITE (LUN,180,ERR=14) XC(K), YC(K)
  180       FORMAT (2F12.6,' lineto')
            IF (K .NE. IWK(IC)) GO TO 4
C
C   Paint the path.
C
          WRITE (LUN,140,ERR=14)
          IF (IC .NE. NC) GO TO 3
    5   CONTINUE
C
C Output the showpage command and end-of-file indicator.
C
      WRITE (LUN,200,ERR=14)
  200 FORMAT ('showpage'/
     .        '%%EOF')
C
C HP's interpreters require a one-byte End-of-PostScript-Job
C   indicator (to eliminate a timeout error message):
C   ASCII 4.
C
      WRITE (LUN,210,ERR=14) CHAR(4)
  210 FORMAT (A1)
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   11 IER = 1
      RETURN
C
C PX or PY is not strictly increasing.
C
   12 IER = 2
      RETURN
C
C DZ = 0.
C
   13 IER = 3
      RETURN
C
C Error writing to unit LUN.
C
   14 IER = 4
      RETURN
C
C Error flag returned by CNTOUR.
C
   15 IER = 5
      RETURN
      END
      SUBROUTINE FVAL (XP,YP,X1,X2,X3,Y1,Y2,Y3,F1,F2,F3,
     .                 FX1,FX2,FX3,FY1,FY2,FY3,SIG1,SIG2,
     .                 SIG3, FP,IER)
      INTEGER IER
      REAL    XP, YP, X1, X2, X3, Y1, Y2, Y3, F1, F2,
     .        F3, FX1, FX2, FX3, FY1, FY2, FY3, SIG1,
     .        SIG2, SIG3, FP
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   03/18/90
C
C   Given function values and gradients at the three ver-
C tices of a triangle containing a point P, this routine
C computes the value of F at P where F interpolates the ver-
C tex data.  Along the triangle arcs, the interpolatory
C function F is the Hermite interpolatory tension spline de-
C fined by the values and tangential gradient components at
C the endpoints, and the derivative in the direction normal
C to the arc varies linearly between the normal gradient
C components at the endpoints.  A first-order C-1 blending
C method is used to extend F to the interior of the trian-
C gle.  Thus, since values and gradients on an arc depend
C only on the vertex data, the method results in C-1 contin-
C uity when used to interpolate over a triangulation.
C
C   The blending method consists of taking F(P) to be the
C weighted sum of the values at P of the three univariate
C Hermite interpolatory tension splines defined on the line
C segments which join the vertices to the opposite sides and
C pass through P.  The tension factors for these splines are
C obtained by linear interpolation between the pair of ten-
C sion factors associated with the triangle sides which join
C at the appropriate vertex.
C
C   A tension factor SIGMA associated with a Hermite interp-
C olatory tension spline is a nonnegative parameter which
C determines the curviness of the spline.  SIGMA = 0 results
C in a cubic spline, and the spline approaches the linear
C interpolant as SIGMA increases.
C
C On input:
C
C       XP,YP = Coordinates of a point P at which an interp-
C               olated value is to be computed.
C
C       X1,X2,X3,Y1,Y2,Y3 = Coordinates of the vertices of a
C                           triangle (V1,V2,V3) containing
C                           P.  V3 is strictly to the left
C                           of V1->V2.
C
C       F1,F2,F3 = Values of the interpolatory function at
C                  the vertices.
C
C       FX1,FX2,FX3 = x components of the gradients of F at
C                     the vertices.
C
C       FY1,FY2,FY3 = y components of the gradients of F at
C                     the vertices.
C
C       SIG1,SIG2,SIG3 = Tension factors associated with the
C                        arcs opposite V1, V2, and V3, re-
C                        spectively.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       FP = Interpolated value at P unless IER < 0, in
C            which case FP is not defined.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if P is not contained in the triangle.
C                     This may result from roundoff error
C                     when P lies near an arc, and the int-
C                     erpolated value FP is valid in that
C                     case.
C             IER = -1 if the triangle vertices are
C                      collinear.
C
C SRFPACK modules required by FVAL:  ARCINT, COORDS, SNHCSH
C
C***********************************************************
C
      INTEGER IERR
      REAL    B, B1, B2, B3, C1, C2, C3, DUM, FQ, FXQ, FYQ,
     .        H1, H2, H3, PX, PY, SIG, SUM, XQ, YQ
C
      PX = XP
      PY = YP
C
C F(P) = C1*H1(P) + C2*H2(P) + C3*H3(P) where C1, C2, and C3
C   are weight functions which sum to 1, and H1, H2, and H3
C   are Hermite interpolatory tension splines on the line
C   segments which join vertices to opposite sides and con-
C   tain P.
C
C Compute barycentric coordinates of P with respect to the
C   triangle.
C
      CALL COORDS (PX,PY,X1,X2,X3,Y1,Y2,Y3, B1,B2,B3,IER)
      IF (IER .NE. 0) RETURN
      IF (B1 .LT. 0.  .OR.  B2 .LT. 0.  .OR.  B3 .LT. 0.)
     .   IER = 1
C
C Compute the coefficients of the partial interpolants.
C   C1 = 1 on the side opposite V1, and C1 = 0 on the other
C   arcs.  Similarly for C2 and C3.
C
      C1 = B2*B3
      C2 = B3*B1
      C3 = B1*B2
      SUM = C1 + C2 + C3
      IF (SUM .EQ. 0.) THEN
C
C P coincides with a vertex.
C
        FP = B1*F1 + B2*F2 + B3*F3
        RETURN
      ENDIF
C
C Normalize the coefficients.
C
      C1 = C1/SUM
      C2 = C2/SUM
      C3 = C3/SUM
C
C For each vertex Vi, compute the intersection Q of the side
C   opposite Vi with the line defined by Vi and P, the value
C   and gradient at Q, and the partial interpolant value Hi
C   at P.
C
C   Side opposite V1:
C
      B = B2/(B2+B3)
      XQ = B*X2 + (1.-B)*X3
      YQ = B*Y2 + (1.-B)*Y3
      SIG = B*SIG3 + (1.-B)*SIG2
      CALL ARCINT (B,X2,X3,Y2,Y3,F2,F3,FX2,FX3,FY2,FY3,SIG1,
     .             .TRUE., FQ,FXQ,FYQ,IERR)
      CALL ARCINT (B1,X1,XQ,Y1,YQ,F1,FQ,FX1,FXQ,FY1,FYQ,SIG,
     .             .FALSE., H1,DUM,DUM,IERR)
C
C   Side opposite V2:
C
      B = B3/(B3+B1)
      XQ = B*X3 + (1.-B)*X1
      YQ = B*Y3 + (1.-B)*Y1
      SIG = B*SIG1 + (1.-B)*SIG3
      CALL ARCINT (B,X3,X1,Y3,Y1,F3,F1,FX3,FX1,FY3,FY1,SIG2,
     .             .TRUE., FQ,FXQ,FYQ,IERR)
      CALL ARCINT (B2,X2,XQ,Y2,YQ,F2,FQ,FX2,FXQ,FY2,FYQ,SIG,
     .             .FALSE., H2,DUM,DUM,IERR)
C
C   Side opposite V3:
C
      B = B1/(B1+B2)
      XQ = B*X1 + (1.-B)*X2
      YQ = B*Y1 + (1.-B)*Y2
      SIG = B*SIG2 + (1.-B)*SIG1
      CALL ARCINT (B,X1,X2,Y1,Y2,F1,F2,FX1,FX2,FY1,FY2,SIG3,
     .             .TRUE., FQ,FXQ,FYQ,IERR)
      CALL ARCINT (B3,X3,XQ,Y3,YQ,F3,FQ,FX3,FXQ,FY3,FYQ,SIG,
     .             .FALSE., H3,DUM,DUM,IERR)
C
C Accumulate the partial interpolant values.
C
      FP = C1*H1 + C2*H2 + C3*H3
      RETURN
      END
      SUBROUTINE GETSIG (N,X,Y,H,LIST,LPTR,LEND,HXHY,
     .                   TOL, SIGMA, DSMAX,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IER
      REAL    X(N), Y(N), H(N), HXHY(2,N), TOL, SIGMA(*),
     .        DSMAX
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/03/98
C
C   Given a triangulation of a set of nodes in the plane,
C along with data values H and gradients (HX,HY) at the
C nodes, this subroutine determines, for each triangulation
C arc, the smallest (nonnegative) tension factor SIGMA such
C that the Hermite interpolatory tension spline H(T), de-
C fined by SIGMA and the endpoint values and directional
C derivatives, preserves local shape properties of the data.
C In order to define the shape properties on an arc, it is
C convenient to map the arc to an interval (T1,T2).  Then,
C denoting the endpoint data values by H1,H2 and the deriva-
C tives by HP1,HP2, and letting S = (H2-H1)/(T2-T1), the
C data properties are
C
C       Monotonicity:  S, HP1, and HP2 are nonnegative or
C                        nonpositive,
C   and
C
C       Convexity:     HP1 .LE. S .LE. HP2  or  HP1 .GE. S
C                        .GE. HP2.
C
C The corresponding properties of H are constant sign of the
C first and second derivatives, respectively.  Note that,
C unless HP1 = S = HP2, infinite tension is required (and H
C is linear on the interval) if S = 0 in the case of mono-
C tonicity, or if HP1 = S or HP2 = S in the case of
C convexity.
C
C   Note that if gradients are to be computed by Subroutine
C GRADG or function values and gradients are computed by
C SMSURF, it may be desirable to alternate those computa-
C tions (which require tension factors) with calls to this
C subroutine.  This iterative procedure should terminate
C with a call to GETSIG in order to ensure that the shape
C properties are preserved, and convergence can be achieved
C (at the cost of optimality) by allowing only increases in
C tension factors (refer to the parameter descriptions for
C SIGMA, DSMAX, and IER).
C
C   Refer to functions SIG0, SIG1, and SIG2 for means of
C selecting minimum tension factors to preserve more general
C properties.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I)).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C       HXHY = Array dimensioned 2 by N whose columns con-
C              tain partial derivatives at the nodes (X
C              partials in the first row).  Refer to Subrou-
C              tines GRADC, GRADG, GRADL, and SMSURF.
C
C       TOL = Tolerance whose magnitude determines how close
C             each tension factor is to its optimal value
C             when nonzero finite tension is necessary and
C             sufficient to satisfy the constraint --
C             abs(TOL) is an upper bound on the magnitude
C             of the smallest (nonnegative) or largest (non-
C             positive) value of the first or second deriva-
C             tive of H in the interval.  Thus, the con-
C             straint is satisfied, but possibly with more
C             tension than necessary.
C
C The above parameters are not altered by this routine.
C
C       SIGMA = Array of length 2*NA = 6*(N-1)-2*NB, where
C               NA and NB are the numbers of arcs and boun-
C               dary nodes, respectively, containing minimum
C               values of the tension factors.  The tension
C               factors are associated with arcs in one-to-
C               one correspondence with LIST entries.  Note
C               that each arc N1-N2 has two LIST entries and
C               thus, the tension factor is stored in both
C               SIGMA(I) and SIGMA(J) where LIST(I) = N2 (in
C               the adjacency list for N1) and LIST(J) = N1
C               (in the list associated with N2).  SIGMA
C               should be set to all zeros if minimal ten-
C               sion is desired, and should be unchanged
C               from a previous call in order to ensure con-
C               vergence of the iterative procedure describ-
C               ed in the header comments.
C
C On output:
C
C       SIGMA = Array containing tension factors for which
C               H(T) preserves the local data properties on
C               each triangulation arc, with the restriction
C               that SIGMA(I) .LE. 85 for all I (unless the
C               input value is larger).  The factors are as
C               small as possible (within the tolerance) but
C               not less than their input values.  If infin-
C               ite tension is required on an arc, the cor-
C               responding factor is SIGMA(I) = 85 (and H
C               is an approximation to the linear inter-
C               polant on the arc), and if neither property
C               is satisfied by the data, then SIGMA(I) = 0
C               (assuming its input value is 0), and thus H
C               is cubic on the arc.
C
C       DSMAX = Maximum increase in a component of SIGMA
C               from its input value.
C
C       IER = Error indicator and information flag:
C             IER = I if no errors were encountered and I
C                     components of SIGMA were altered from
C                     their input values for I .GE. 0.
C             IER = -1 if N < 3.  SIGMA is not altered in
C                      this case.
C             IER = -2 if duplicate nodes were encountered.
C
C TRIPACK modules required by GETSIG:  LSTPTR, STORE
C
C SRFPACK module required by GETSIG:  SNHCSH
C
C Intrinsic functions called by GETSIG:  ABS, EXP, MAX, MIN,
C                                          SIGN, SQRT
C
C***********************************************************
C
      INTEGER LSTPTR
      REAL    STORE
      INTEGER ICNT, LP1, LP2, LPL, LUN, N1, N2, NIT, NM1
      REAL    A, C1, C2, COSHM, COSHMM, D0, D1, D1D2, D1PD2,
     .        D2, DMAX, DSIG, DSM, DT, DX, DY, E, EMS, EMS2,
     .        F, F0, FMAX, FNEG, FP, FTOL, RTOL, S, S1, S2,
     .        SBIG, SCM, SGN, SIG, SIGIN, SINHM, SSINH, SSM,
     .        STOL, T, T0, T1, T2, TM, TP1
C
      DATA SBIG/85./,  LUN/-1/
      NM1 = N - 1
      IF (NM1 .LT. 2) GO TO 11
C
C Compute an absolute tolerance FTOL = abs(TOL) and a
C   relative tolerance RTOL = 100*Macheps.
C
      FTOL = ABS(TOL)
      RTOL = 1.
    1 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 1
      RTOL = RTOL*200.
C
C Print a heading.
C
      IF (LUN .GE. 0) WRITE (LUN,100) N, FTOL
  100 FORMAT (///13X,'GETSIG:  N =',I4,', TOL = ',E10.3//)
C
C Initialize change counter ICNT and maximum change DSM for
C   the loop on arcs.
C
      ICNT = 0
      DSM = 0.
C
C Loop on arcs N1-N2 for which N2 > N1.  LPL points to the
C   last neighbor of N1.
C
      DO 10 N1 = 1,NM1
        LPL = LEND(N1)
        LP1 = LPL
C
C   Top of loop on neighbors N2 of N1.
C
    2   LP1 = LPTR(LP1)
        N2 = ABS(LIST(LP1))
        IF (N2 .LE. N1) GO TO 9
C
C Print a message and compute parameters for the arc:  DT =
C   arc length and SIGIN = input SIGMA value.
C
        IF (LUN .GE. 0) WRITE (LUN,110) N1, N2
  110   FORMAT (/1X,'Arc',I4,' -',I4)
        DX = X(N2) - X(N1)
        DY = Y(N2) - Y(N1)
        DT = SQRT(DX*DX + DY*DY)
        IF (DT .EQ. 0.) GO TO 12
        SIGIN = SIGMA(LP1)
        IF (SIGIN .GE. SBIG) GO TO 9
C
C Compute scaled directional derivatives S1,S2 at the end-
C   points (for the direction N1->N2), first difference S,
C   and second differences D1,D2.
C
        S1 = HXHY(1,N1)*DX + HXHY(2,N1)*DY
        S2 = HXHY(1,N2)*DX + HXHY(2,N2)*DY
        S = H(N2) - H(N1)
        D1 = S - S1
        D2 = S2 - S
        D1D2 = D1*D2
C
C Test for infinite tension required to satisfy either
C   property.
C
        SIG = SBIG
        IF ((D1D2 .EQ. 0.  .AND.  S1 .NE. S2)  .OR.
     .      (S .EQ. 0.  .AND.  S1*S2 .GT. 0.)) GO TO 8
C
C Test for SIGMA = 0 sufficient.  The data satisfies convex-
C   ity iff D1D2 .GE. 0, and D1D2 = 0 implies S1 = S = S2.
C
        SIG = 0.
        IF (D1D2 .LT. 0.) GO TO 4
        IF (D1D2 .EQ. 0.) GO TO 8
        T = MAX(D1/D2,D2/D1)
        IF (T .LE. 2.) GO TO 8
        TP1 = T + 1.
C
C Convexity:  find a zero of F(SIG) = SIG*COSHM(SIG)/
C   SINHM(SIG) - TP1.
C
C   F(0) = 2-T < 0, F(TP1) .GE. 0, the derivative of F
C     vanishes at SIG = 0, and the second derivative of F is
C     .2 at SIG = 0.  A quadratic approximation is used to
C     obtain a starting point for the Newton method.
C
        SIG = SQRT(10.*T-20.)
        NIT = 0
C
C   Top of loop:
C
    3   IF (SIG .LE. .5) THEN
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          T1 = COSHM/SINHM
          FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.)
        ELSE
C
C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
C     overflow with large SIG.
C
          EMS = EXP(-SIG)
          SSM = 1. - EMS*(EMS+SIG+SIG)
          T1 = (1.-EMS)*(1.-EMS)/SSM
          FP = T1 + SIG*(2.*SIG*EMS/SSM - T1*T1 + 1.)
        ENDIF
C
        F = SIG*T1 - TP1
        IF (LUN .GE. 0) WRITE (LUN,120) SIG, F, FP
  120   FORMAT (1X,'Convexity:  SIG = ',E15.8,
     .          ', F(SIG) = ',E15.8/1X,35X,'FP(SIG) = ',
     .          E15.8)
        NIT = NIT + 1
C
C   Test for convergence.
C
        IF (FP .LE. 0.) GO TO 8
        DSIG = -F/FP
        IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.  .AND.
     .      F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 8
C
C   Update SIG.
C
        SIG = SIG + DSIG
        GO TO 3
C
C Convexity cannot be satisfied.  Monotonicity can be satis-
C   fied iff S1*S .GE. 0 and S2*S .GE. 0 since S .NE. 0.
C
    4   IF (S1*S .LT. 0.  .OR.  S2*S .LT. 0.) GO TO 8
        T0 = 3.*S - S1 - S2
        D0 = T0*T0 - S1*S2
C
C SIGMA = 0 is sufficient for monotonicity iff S*T0 .GE. 0
C   or D0 .LE. 0.
C
        IF (D0 .LE. 0.  .OR.  S*T0 .GE. 0.) GO TO 8
C
C Monotonicity:  find a zero of F(SIG) = sign(S)*HP(R),
C   where HPP(R) = 0 and HP, HPP denote derivatives of H.
C   F has a unique zero, F(0) < 0, and F approaches
C   abs(S) as SIG increases.
C
C   Initialize parameters for the secant method.  The method
C     uses three points:  (SG0,F0), (SIG,F), and
C     (SNEG,FNEG), where SG0 and SNEG are defined implicitly
C     by DSIG = SIG - SG0 and DMAX = SIG - SNEG.
C
        SGN = SIGN(1.,S)
        SIG = SBIG
        FMAX = SGN*(SIG*S-S1-S2)/(SIG-2.)
        IF (FMAX .LE. 0.) GO TO 8
        STOL = RTOL*SIG
        F = FMAX
        F0 = SGN*D0/(3.*(D1-D2))
        FNEG = F0
        DSIG = SIG
        DMAX = SIG
        D1PD2 = D1 + D2
        NIT = 0
C
C   Top of loop:  compute the change in SIG by linear
C     interpolation.
C
    5   DSIG = -F*DSIG/(F-F0)
        IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130   FORMAT (1X,'Monotonicity:  DSIG = ',E15.8)
        IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .       DSIG*DMAX .GT. 0. ) GO TO 7
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
        IF (ABS(DSIG) .LT. STOL/2.) DSIG = -SIGN(STOL/2.,
     .                              DMAX)
C
C   Update SIG, F0, and F.
C
        SIG = SIG + DSIG
        F0 = F
        IF (SIG .LE. .5) THEN
C
C   Use approximations to the hyperbolic functions designed
C     to avoid cancellation error with small SIG.
C
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          C1 = SIG*COSHM*D2 - SINHM*D1PD2
          C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
          A = C2 - C1
          E = SIG*SINHM - COSHMM - COSHMM
        ELSE
C
C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
C     overflow with large SIG.
C
          EMS = EXP(-SIG)
          EMS2 = EMS + EMS
          TM = 1. - EMS
          SSINH = TM*(1.+EMS)
          SSM = SSINH - SIG*EMS2
          SCM = TM*TM
          C1 = SIG*SCM*D2 - SSM*D1PD2
          C2 = SIG*SSINH*D2 - SCM*D1PD2
C
C   R is in (0,1) and well-defined iff HPP(T1)*HPP(T2) < 0.
C
          F = FMAX
          IF (C1*(SIG*SCM*D1 - SSM*D1PD2) .GE. 0.) GO TO 6
          A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
          IF (A*(C2+C1) .LT. 0.) GO TO 6
          E = SIG*SSINH - SCM - SCM
        ENDIF
C
        F = (SGN*(E*S2-C2) + SQRT(A*(C2+C1)))/E
C
C   Update the number of iterations NIT.
C
    6   NIT = NIT + 1
        IF (LUN .GE. 0) WRITE (LUN,140) NIT, SIG, F
  140   FORMAT (1X,11X,I2,' -- SIG = ',E15.8,', F = ',
     .          E15.8)
C
C   Test for convergence.
C
        STOL = RTOL*SIG
        IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.  .AND.
     .      F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 8
        DMAX = DMAX + DSIG
        IF (F0*F .GT. 0.  .AND.  ABS(F) .GE. ABS(F0))
     .     GO TO 7
        IF (F0*F .LE. 0.) THEN
C
C   F and F0 have opposite signs.  Update (SNEG,FNEG) to
C     (SG0,F0) so that F and FNEG always have opposite
C     signs.  If SIG is closer to SNEG than SG0 and abs(F)
C     < abs(FNEG), then swap (SNEG,FNEG) with (SG0,F0).
C
          T1 = DMAX
          T2 = FNEG
          DMAX = DSIG
          FNEG = F0
          IF ( ABS(DSIG) .GT. ABS(T1)  .AND.
     .         ABS(F) .LT. ABS(T2) ) THEN
C
            DSIG = T1
            F0 = T2
          ENDIF
        ENDIF
        GO TO 5
C
C   Bottom of loop:  F0*F > 0 and the new estimate would
C     be outside of the bracketing interval of length
C     abs(DMAX).  Reset (SG0,F0) to (SNEG,FNEG).
C
    7   DSIG = DMAX
        F0 = FNEG
        GO TO 5
C
C  Update SIGMA, ICNT, and DSM if necessary.
C
    8   SIG = MIN(SIG,SBIG)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(LP1) = SIG
          LP2 = LSTPTR(LEND(N2),N1,LIST,LPTR)
          SIGMA(LP2) = SIG
          ICNT = ICNT + 1
          DSM = MAX(DSM,SIG-SIGIN)
        ENDIF
C
C Bottom of loop on neighbors N2 of N1.
C
    9   IF (LP1 .NE. LPL) GO TO 2
   10   CONTINUE
C
C No errors encountered.
C
      DSMAX = DSM
      IER = ICNT
      RETURN
C
C N < 3
C
   11 DSMAX = 0.
      IER = -1
      RETURN
C
C Nodes N1 and N2 coincide.
C
   12 DSMAX = DSM
      IER = -2
      RETURN
      END
      SUBROUTINE GIVENS ( A,B, C,S)
      REAL A, B, C, S
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This subroutine constructs the Givens plane rotation,
C
C           ( C  S)
C       G = (     ) , where C*C + S*S = 1,
C           (-S  C)
C
C which zeros the second component of the vector (A,B)**T
C (transposed).  Subroutine ROTATE may be called to apply
C the transformation to a 2 by N matrix.
C
C   This routine is identical to Subroutine SROTG from the
C LINPACK BLAS (Basic Linear Algebra Subroutines).
C
C On input:
C
C       A,B = Components of the vector defining the rota-
C             tion.  These are overwritten by values R
C             and Z (described below) which define C and S.
C
C On output:
C
C       A = Signed Euclidean norm R of the input vector:
C           R = +/-SQRT(A*A + B*B)
C
C       B = Value Z such that:
C             C = SQRT(1-Z*Z) and S=Z if ABS(Z) .LE. 1, and
C             C = 1/Z and S = SQRT(1-C*C) if ABS(Z) > 1.
C
C       C = +/-(A/R) or 1 if R = 0.
C
C       S = +/-(B/R) or 0 if R = 0.
C
C Modules required by GIVENS:  None
C
C Intrinsic functions called by GIVENS:  ABS, SQRT
C
C***********************************************************
C
      REAL AA, BB, R, U, V
C
C Local parameters:
C
C AA,BB = Local copies of A and B
C R =     C*A + S*B = +/-SQRT(A*A+B*B)
C U,V =   Variables used to scale A and B for computing R
C
      AA = A
      BB = B
      IF (ABS(AA) .LE. ABS(BB)) GO TO 1
C
C ABS(A) > ABS(B).
C
      U = AA + AA
      V = BB/U
      R = SQRT(.25 + V*V) * U
      C = AA/R
      S = V * (C + C)
C
C Note that R has the sign of A, C > 0, and S has
C   SIGN(A)*SIGN(B).
C
      B = S
      A = R
      RETURN
C
C ABS(A) .LE. ABS(B).
C
    1 IF (BB .EQ. 0.) GO TO 2
      U = BB + BB
      V = AA/U
C
C Store R in A.
C
      A = SQRT(.25 + V*V) * U
      S = BB/A
      C = V * (S + S)
C
C Note that R has the sign of B, S > 0, and C has
C   SIGN(A)*SIGN(B).
C
      B = 1.
      IF (C .NE. 0.) B = 1./C
      RETURN
C
C A = B = 0.
C
    2 C = 1.
      S = 0.
      RETURN
      END
      SUBROUTINE GRADC (K,NCC,LCC,N,X,Y,Z,LIST,LPTR,
     .                  LEND, DX,DY,DXX,DXY,DYY,IER)
      INTEGER K, NCC, LCC(*), N, LIST(*), LPTR(*),
     .        LEND(N), IER
      REAL    X(N), Y(N), Z(N), DX, DY, DXX, DXY, DYY
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/22/97
C
C   Given a Delaunay triangulation of N points in the plane
C with associated data values Z, this subroutine estimates
C first and second partial derivatives at node K.  The der-
C ivatives are taken to be the partials at K of a cubic
C function which interpolates Z(K) and fits the data values
C at a set of nearby nodes in a weighted least squares
C sense.  A Marquardt stabilization factor is used if neces-
C sary to ensure a well-conditioned system.  Thus, a unique
C solution exists if there are at least 10 noncollinear
C nodes.
C
C   The triangulation may include constraints introduced by
C Subroutine ADDCST, in which case the derivative estimates
C are influenced by the nonconvex geometry of the domain.
C Refer to Subroutine GETNP.  If data values at the con-
C straint nodes are not known, Subroutine ZGRADL, which
C computes approximate data values at constraint nodes along
C with gradients, should be called in place of this routine.
C
C   Subroutine GRADL uses a quadratic polynomial instead of
C the cubic and may be more accurate if the nodal distribu-
C tion is sparse.  Another alternative routine, GRADG,
C employs a global method to compute the first partial
C derivatives at all of the nodes at once.  That method
C is usually more efficient (when all first partials are
C needed) and may be more accurate, depending on the data.
C
C On input:
C
C       K = Index of the node at which derivatives are to be
C           estimated.  1 .LE. K .LE. N.
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index of the
C             first node of constraint I in LCC(I).  For I =
C             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1.
C
C       N = Number of nodes in the triangulation.
C           N .GE. 10.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       Z = Array of length N containing data values associ-
C           ated with the nodes.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       DX,DY = Estimated first partial derivatives at node
C               K unless IER < 0.
C
C       DXX,DXY,DYY = Estimated second partial derivatives
C                     at node K unless IER < 0.
C
C       IER = Error indicator:
C             IER = L > 0 if no errors were encountered and
C                         L nodes (including node K) were
C                         employed in the least squares fit.
C             IER = -1 if K, NCC, an LCC entry, or N is
C                      outside its valid range on input.
C             IER = -2 if all nodes are collinear.
C
C TRIPACK modules required by GRADC:  GETNP, INTSEC
C
C SRFPACK modules required by GRADC:  GIVENS, ROTATE, SETRO3
C
C Intrinsic functions called by GRADC:  ABS, MIN, REAL, SQRT
C
C***********************************************************
C
      INTEGER   LMN, LMX
      PARAMETER (LMN=14,  LMX=30)
      INTEGER   I, IERR, J, JP1, KK, L, LMAX, LMIN, LM1,
     .          LNP, NP, NPTS(LMX)
      REAL      A(10,10), C, DIST(LMX), DMIN, DS, DTOL, RIN,
     .          RS, RTOL, S, SF, SFC, SFS, STF, SUM, W, XK,
     .          YK, ZK
      DATA      RTOL/1.E-5/, DTOL/.01/
C
C Local parameters:
C
C A =         Transpose of the augmented regression matrix
C C =         First component of the plane rotation deter-
C               mined by Subroutine GIVENS
C DIST =      Array containing the distances between K and
C               the elements of NPTS (refer to GETNP)
C DMIN =      Minimum of the magnitudes of the diagonal
C               elements of the regression matrix after
C               zeros are introduced below the diagonal
C DS =        Squared distance between nodes K and NPTS(LNP)
C DTOL =      Tolerance for detecting an ill-conditioned
C               system.  The system is accepted when DMIN/W
C               .GE. DTOL.
C I =         DO-loop index
C IERR =      Error flag for calls to GETNP
C J =         DO-loop index
C JP1 =       J+1
C KK =        Local copy of K
C L =         Number of columns of A**T to which a rotation
C               is applied
C LMAX,LMIN = Min(LMX,N), Min(LMN,N)
C LMN,LMX =   Minimum and maximum values of LNP for N
C               sufficiently large.  In most cases LMN-1
C               nodes are used in the fit.  4 .LE. LMN .LE.
C               LMX.
C LM1 =       LMIN-1 or LNP-1
C LNP =       Length of NPTS
C NP =        Element of NPTS to be added to the system
C NPTS =      Array containing the indexes of a sequence of
C               nodes ordered by distance from K.  NPTS(1)=K
C               and the first LNP-1 elements of NPTS are
C               used in the least squares fit.  Unless LNP
C               exceeds LMAX, NPTS(LNP) determines R.
C RIN =       Inverse of the distance R between node K and
C               NPTS(LNP) or some point further from K than
C               NPTS(LMAX) if NPTS(LMAX) is used in the fit.
C               R is a radius of influence which enters into
C               the weight W.
C RS =        R*R
C RTOL =      Tolerance for determining R.  If the relative
C               change in DS between two elements of NPTS is
C               not greater than RTOL, they are treated as
C               being the same distance from node K.
C S =         Second component of the plane rotation deter-
C               mined by Subroutine GIVENS
C SF =        Scale factor for the linear terms (columns 8
C               and 9) in the least squares fit -- inverse
C               of the root-mean-square distance between K
C               and the nodes (other than K) in the least
C               squares fit
C SFS =       Scale factor for the quadratic terms (columns
C               5, 6, and 7) in the least squares fit --
C               SF*SF
C SFC =       Scale factor for the cubic terms (first 4
C               columns) in the least squares fit -- SF**3
C STF =       Marquardt stabilization factor used to damp
C               out the first 4 solution components (third
C               partials of the cubic) when the system is
C               ill-conditioned.  As STF increases, the
C               fitting function approaches a quadratic
C               polynomial.
C SUM =       Sum of squared distances between node K and
C               the nodes used in the least squares fit
C W =         Weight associated with a row of the augmented
C               regression matrix -- 1/D - 1/R, where D < R
C               and D is the distance between K and a node
C               entering into the least squares fit
C XK,YK,ZK =  Coordinates and data value associated with K
C
      KK = K
C
C Test for errors and initialize LMIN and LMAX.
C
      IF (KK .LT. 1  .OR.  KK .GT. N  .OR.  NCC .LT. 0
     .    .OR.  N .LT. 10) GO TO 13
      LMIN = MIN(LMN,N)
      LMAX = MIN(LMX,N)
C
C Compute NPTS, DIST, LNP, SF, SFS, SFC, and RIN --
C
C   Set NPTS to the closest LMIN-1 nodes to K.
C
      SUM = 0.
      NPTS(1) = KK
      DIST(1) = 0.
      LM1 = LMIN - 1
      DO 1 LNP = 2,LM1
        CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,
     .              LNP, NPTS,DIST, IERR)
        IF (IERR .NE. 0) GO TO 13
        DS = DIST(LNP)**2
        SUM = SUM + DS
    1   CONTINUE
C
C Add additional nodes to NPTS until the relative increase
C   in DS is at least RTOL.
C
      DO 3 LNP = LMIN,LMAX
        CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,
     .              LNP, NPTS,DIST, IERR)
        RS = DIST(LNP)**2
        IF ((RS-DS)/DS .LE. RTOL) GO TO 2
        IF (LNP .GT. 10) GO TO 4
    2   SUM = SUM + RS
    3   CONTINUE
C
C Use all LMAX nodes in the least squares fit.  RS is
C   arbitrarily increased by 10 per cent.
C
      RS = 1.1*RS
      LNP = LMAX + 1
C
C There are LNP-2 equations corresponding to nodes NPTS(2),
C   ...,NPTS(LNP-1).
C
    4 SFS = REAL(LNP-2)/SUM
      SF = SQRT(SFS)
      SFC = SF*SFS
      RIN = 1./SQRT(RS)
      XK = X(KK)
      YK = Y(KK)
      ZK = Z(KK)
C
C A Q-R decomposition is used to solve the least squares
C   system.  The transpose of the augmented regression
C   matrix is stored in A with columns (rows of A) defined
C   as follows:  1-4 are the cubic terms, 5-7 are the quad-
C   ratic terms with coefficients DXX/2, DXY, and DYY/2,
C   8 and 9 are the linear terms with coefficients DX and
C   DY, and the last column is the right hand side.
C
C Set up the first 9 equations and zero out the lower tri-
C   angle with Givens rotations.
C
      DO 6 I = 1,9
        NP = NPTS(I+1)
        W = 1./DIST(I+1) - RIN
        CALL SETRO3 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,
     .               SFC,W, A(1,I))
        IF (I .EQ. 1) GO TO 6
        DO 5 J = 1,I-1
          JP1 = J + 1
          L = 10 - J
          CALL GIVENS (A(J,J),A(J,I),C,S)
          CALL ROTATE (L,C,S,A(JP1,J),A(JP1,I))
    5     CONTINUE
    6   CONTINUE
C
C Add the additional equations to the system using
C   the last column of A.  I .LE. LNP.
C
      I = 11
    7   IF (I .LT. LNP) THEN
          NP = NPTS(I)
          W = 1./DIST(I) - RIN
          CALL SETRO3 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,
     .                 SFC,W, A(1,10))
          DO 8 J = 1,9
            JP1 = J + 1
            L = 10 - J
            CALL GIVENS (A(J,J),A(J,10),C,S)
            CALL ROTATE (L,C,S,A(JP1,J),A(JP1,10))
    8       CONTINUE
          I = I + 1
          GO TO 7
        ENDIF
C
C Test the system for ill-conditioning.
C
      DMIN = MIN( ABS(A(1,1)),ABS(A(2,2)),ABS(A(3,3)),
     .            ABS(A(4,4)),ABS(A(5,5)),ABS(A(6,6)),
     .            ABS(A(7,7)),ABS(A(8,8)),ABS(A(9,9)) )
      IF (DMIN/W .GE. DTOL) GO TO 12
      IF (LNP .LE. LMAX) THEN
C
C   Add another node to the system and increase R.  Note
C     that I = LNP.
C
        LNP = LNP + 1
        IF (LNP .LE. LMAX) THEN
          CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,
     .                LNP, NPTS,DIST, IERR)
          RS = DIST(LNP)**2
        ENDIF
        RIN = 1./SQRT(1.1*RS)
        GO TO 7
      ENDIF
C
C Stabilize the system by damping third partials -- add
C   multiples of the first four unit vectors to the first
C   four equations.
C
      STF = W
      DO 11 I = 1,4
        A(I,10) = STF
        DO 9 J = I+1,10
          A(J,10) = 0.
    9     CONTINUE
        DO 10 J = I,9
          JP1 = J + 1
          L = 10 - J
          CALL GIVENS (A(J,J),A(J,10),C,S)
          CALL ROTATE (L,C,S,A(JP1,J),A(JP1,10))
   10     CONTINUE
   11   CONTINUE
C
C Test the damped system for ill-conditioning.
C
      DMIN = MIN( ABS(A(5,5)),ABS(A(6,6)),ABS(A(7,7)),
     .            ABS(A(8,8)),ABS(A(9,9)) )
      IF (DMIN/W .LT. DTOL) GO TO 14
C
C Solve the 9 by 9 triangular system for the last 5
C   components (first and second partial derivatives).
C
   12 DY = A(10,9)/A(9,9)
      DX = (A(10,8)-A(9,8)*DY)/A(8,8)
      DYY = (A(10,7)-A(8,7)*DX-A(9,7)*DY)/A(7,7)
      DXY = (A(10,6)-A(7,6)*DYY-A(8,6)*DX-A(9,6)*DY)/A(6,6)
      DXX = (A(10,5)-A(6,5)*DXY-A(7,5)*DYY-A(8,5)*DX-
     .       A(9,5)*DY)/A(5,5)
C
C Scale the solution components.
C
      DX = SF*DX
      DY = SF*DY
      DXX = 2.*SFS*DXX
      DXY = SFS*DXY
      DYY = 2.*SFS*DYY
      IER = LNP - 1
      RETURN
C
C Invalid input parameter.
C
   13 IER = -1
      RETURN
C
C No unique solution due to collinear nodes.
C
   14 IER = -2
      RETURN
      END
      SUBROUTINE GRADG (NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,
     .                  IFLGS,SIGMA, NIT,DGMAX,GRAD, IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        IFLGS, NIT, IER
      REAL    X(N), Y(N), Z(N), SIGMA(*), DGMAX, GRAD(2,N)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/12/94
C
C   Given a triangulation of N nodes in the plane, along
C with data values at the nodes and tension factors associ-
C ated with the arcs, this subroutine employs a global
C method to compute estimated gradients at the nodes.  The
C method consists of minimizing a quadratic functional Q(G)
C over vectors G of length 2N (gradient components), where
C Q is an approximation to the linearized curvature over the
C triangulation of a C-1 bivariate function F(X,Y) which
C interpolates the nodal values and gradients.
C
C   The restriction of F to an arc of the triangulation is
C taken to be the Hermite interpolatory tension spline
C defined by the data values and tangential gradient compo-
C nents at the endpoints of the arc, and Q is the sum over
C the triangulation arcs, excluding interior constraint
C arcs, of the linearized curvatures of F along the arcs --
C the integrals over the arcs of D2F(T)**2, where D2F(T) is
C the second derivative of F with respect to distance T
C along the arc.
C
C   Subroutines INTRC1 and UNIF may be called to evaluate F
C at arbitrary points.  The interpolant F is further de-
C scribed in Subroutines FVAL and TVAL, and Q is identical
C to the functional Q1 described in Subroutine SMSURF.
C
C   The minimization problem corresponds to an order 2N
C symmetric positive definite sparse linear system which is
C solved for the X and Y partial derivatives by the block
C Gauss-Seidel method with N blocks of order 2.
C
C   If constraints are present and data values at the con-
C straint nodes are not known, Subroutine ZGRADG, which
C computes approximate data values at constraint nodes
C along with the gradients, should be called in place of
C this routine.
C
C   An alternative method, Subroutine GRADC or GRADL, com-
C putes a local approximation to the partials at a single
C node and may be more accurate, depending on the data
C values and distribution of nodes (neither method emerged
C as superior in tests for accuracy).  If all gradients are
C required and a uniform tension factor SIGMA = 0 is used,
C GRADG is significantly faster than either GRADC or GRADL.
C
C On input:
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index of the
C             first node of constraint I in LCC(I).  For I =
C             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       Z = Array of length N containing data values at the
C           nodes.  Z(I) is associated with (X(I),Y(I)).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C       IFLGS = Tension factor option:
C               IFLGS .LE. 0 if a single uniform tension
C                            factor is to be used.
C               IFLGS .GE. 1 if variable tension is desired.
C
C       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               routines GETSIG, SIG0, SIG1, and SIG2.
C
C The above parameters are not altered by this routine.
C
C       NIT = Maximum number of Gauss-Seidel iterations to
C             be employed.  This maximum will likely be
C             achieved if DGMAX is smaller than the machine
C             precision.  Note that complete convergence is
C             not necessary to achieve maximum accuracy of
C             the interpolant.  For SIGMA = 0, optimal ef-
C             ficiency was achieved in testing with DGMAX =
C             0, and NIT = 3 or 4.  NIT > 0.
C
C       DGMAX = Nonnegative convergence criterion.  The
C               method is terminated when the maximum change
C               in a gradient between iterations is at most
C               DGMAX.  The change in a gradient is taken to
C               be the Euclidean norm of the difference rel-
C               ative to 1 plus the norm of the old value.
C               DGMAX = 1.E-3 is sufficient for effective
C               convergence.
C
C       GRAD = 2 by N array whose columns contain initial
C              estimates of the partial derivatives.  Zero
C              vectors are sufficient.
C
C On output:
C
C       NIT = Number of Gauss-Seidel iterations employed.
C
C       DGMAX = Maximum relative change in a gradient at the
C               last iteration.
C
C       GRAD = Estimated X and Y partial derivatives at the
C              nodes with X partials in the first row.  Grad
C              is not altered if IER = -1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     convergence criterion was achieved.
C             IER = 1 if no errors were encountered but con-
C                     vergence was not achieved within NIT
C                     iterations.
C             IER = -1 if NCC, an LCC entry, N, NIT, or
C                      DGMAX is outside its valid range on
C                      input.
C             IER = -2 if all nodes are collinear or the
C                      triangulation data structure is in-
C                      valid.
C             IER = -3 if duplicate nodes were encountered.
C
C SRFPACK modules required by GRADG:  GRCOEF, SNHCSH
C
C Intrinsic functions called by GRADG:  ABS, MAX, SQRT
C
C***********************************************************
C
      INTEGER I, IFL, IFRST, ILAST, ITER, J, K, KBAK, KFOR,
     .        LCC1, LP, LPJ, LPL, MAXIT, NB, NN
      REAL    A11, A12, A22, D, DCUB, DELX, DELXS, DELY,
     .        DELYS, DET, DF, DGMX, DSQ, DZX, DZY, R1, R2,
     .        SDF, SIG, T, TOL, XK, YK, ZK, ZXK, ZYK
C
      NN = N
      IFL = IFLGS
      MAXIT = NIT
      TOL = DGMAX
C
C Test for errors in input parameters.
C
      IF (NCC .LT. 0  .OR.  MAXIT .LT. 1  .OR.  TOL .LT. 0.)
     .  GO TO 9
      LCC1 = NN+1
      IF (NCC .EQ. 0) THEN
        IF (NN .LT. 3) GO TO 9
      ELSE
        DO 1 I = NCC,1,-1
          IF (LCC1-LCC(I) .LT. 3) GO TO 9
          LCC1 = LCC(I)
    1     CONTINUE
        IF (LCC1 .LT. 1) GO TO 9
      ENDIF
C
C Initialize iteration count and SIG (overwritten if
C   IFLGS > 0).
C
      ITER = 0
      SIG = SIGMA(1)
C
C Top of iteration loop:  If K is a constraint node, I
C   indexes the constraint containing node K, IFRST and
C   ILAST are the first and last nodes of constraint I,
C   and (KBAK,K,KFOR) is a subsequence of constraint I.
C
    2 IF (ITER .EQ. MAXIT) GO TO 8
      DGMX = 0.
      I = 0
      IFRST = 1
      ILAST = LCC1-1
      KBAK = 0
      KFOR = 0
C
C Loop on nodes.
C
      DO 7 K = 1,NN
        IF (K .GE. LCC1) THEN
          IF (K .GT. ILAST) THEN
            I = I + 1
            IFRST = K
            IF (I .LT. NCC) THEN
              ILAST = LCC(I+1) - 1
            ELSE
              ILAST = NN
            ENDIF
            KBAK = ILAST
            KFOR = K + 1
          ELSE
            KBAK = K - 1
            IF (K .LT. ILAST) THEN
              KFOR = K + 1
            ELSE
              KFOR = IFRST
            ENDIF
          ENDIF
        ENDIF
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        ZXK = GRAD(1,K)
        ZYK = GRAD(2,K)
C
C   Initialize components of the order 2 system for the
C     change (DZX,DZY) in the K-th solution components
C     (symmetric matrix in A and residual in R).
C
        A11 = 0.
        A12 = 0.
        A22 = 0.
        R1 = 0.
        R2 = 0.
C
C   Loop on neighbors J of node K.
C
        LPL = LEND(K)
        LPJ = LPL
    3   LPJ = LPTR(LPJ)
          J = ABS(LIST(LPJ))
C
C   Arc K-J lies in a constraint region and is bypassed iff
C     K and J are nodes in the same constraint and J follows
C     KFOR and precedes KBAK as a neighbor of K.
C
          IF (K .LT. LCC1  .OR.  J .LT. IFRST  .OR.
     .        J .GT. ILAST) GO TO 5
          IF (J .EQ. KBAK  .OR.  J .EQ. KFOR) GO TO 5
          LP = LPJ
C
    4     LP = LPTR(LP)
            NB = ABS(LIST(LP))
            IF (NB .EQ. KBAK) GO TO 6
            IF (NB .NE. KFOR) GO TO 4
C
C   Compute parameters associated with edge
C     K->J, and test for duplicate nodes.
C
    5     DELX = X(J) - XK
          DELY = Y(J) - YK
          DELXS = DELX*DELX
          DELYS = DELY*DELY
          DSQ = DELXS + DELYS
          D = SQRT(DSQ)
          DCUB = D*DSQ
          IF (D .EQ. 0.) GO TO 11
          IF (IFL .GE. 1) SIG = SIGMA(LPJ)
          CALL GRCOEF (SIG,DCUB, DF,SDF)
C
C   Update the system components for node J.  The contribu-
C     tion from edge K->J is weighted by 1/D, where D is
C     the arc length.
C
          A11 = A11 + DF*DELXS/D
          A12 = A12 + DF*DELX*DELY/D
          A22 = A22 + DF*DELYS/D
          T = ((DF+SDF)*(Z(J)-ZK) - DF*(ZXK*DELX + ZYK*DELY)
     .          - SDF*(GRAD(1,J)*DELX + GRAD(2,J)*DELY))/D
          R1 = R1 + T*DELX
          R2 = R2 + T*DELY
C
C   Bottom of loop on neighbors.
C
    6     IF (LPJ .NE. LPL) GO TO 3
C
C   Solve the system associated with the K-th block.
C
        DET = A11*A22 - A12*A12
        IF (DET .EQ. 0.  .OR.  A11 .EQ. 0.) GO TO 10
        DZY = (A11*R2 - A12*R1)/DET
        DZX = (R1 - A12*DZY)/A11
C
C   Update the partials at node K and the maximum relative
C     change DGMX.
C
        GRAD(1,K) = ZXK + DZX
        GRAD(2,K) = ZYK + DZY
        DGMX = MAX(DGMX,SQRT(DZX*DZX+DZY*DZY)/
     .             (1.+SQRT(ZXK*ZXK+ZYK*ZYK)))
    7   CONTINUE
C
C   Increment ITER and test for convergence.
C
      ITER = ITER + 1
      IF (DGMX .GT. TOL) GO TO 2
C
C Method converged.
C
      NIT = ITER
      DGMAX = DGMX
      IER = 0
      RETURN
C
C Method failed to converge within NIT iterations.
C
    8 DGMAX = DGMX
      IER = 1
      RETURN
C
C Invalid input parameter.
C
    9 NIT = 0
      DGMAX = 0.
      IER = -1
      RETURN
C
C Node K and its neighbors are collinear, resulting in a
C   singular system.
C
   10 NIT = 0
      DGMAX = DGMX
      IER = -2
      RETURN
C
C Nodes K and J coincide.
C
   11 NIT = 0
      DGMAX = DGMX
      IER = -3
      RETURN
      END
      SUBROUTINE GRADL (K,NCC,LCC,N,X,Y,Z,LIST,LPTR,
     .                  LEND, DX,DY,IER)
      INTEGER K, NCC, LCC(*), N, LIST(*), LPTR(*),
     .        LEND(N), IER
      REAL    X(N), Y(N), Z(N), DX, DY
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/22/97
C
C   Given a Delaunay triangulation of N points in the plane
C with associated data values Z, this subroutine estimates
C X and Y partial derivatives at node K.  The derivatives
C are taken to be the partials at K of a quadratic function
C which interpolates Z(K) and fits the data values at a set
C of nearby nodes in a weighted least squares sense. A Mar-
C quardt stabilization factor is used if necessary to ensure
C a well-conditioned system.  Thus, a unique solution exists
C if there are at least 6 noncollinear nodes.
C
C   The triangulation may include constraints introduced by
C Subroutine ADDCST, in which case the gradient estimates
C are influenced by the nonconvex geometry of the domain.
C Refer to Subroutine GETNP.  If data values at the con-
C straint nodes are not known, Subroutine ZGRADL, which
C computes approximate data values at constraint nodes along
C with gradients, should be called in place of this routine.
C
C   Subroutine GRADC uses a cubic polynomial instead of the
C quadratic and is generally more accurate than this routine
C if the nodal distribution is sufficiently dense.  Another
C alternative routine, GRADG, employs a global method to
C compute the partial derivatives at all of the nodes at
C once.  That method is usually more efficient (when all
C partials are needed) and may be more accurate, depending
C on the data.
C
C On input:
C
C       K = Index of the node at which derivatives are to be
C           estimated.  1 .LE. K .LE. N.
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index of the
C             first node of constraint I in LCC(I).  For I =
C             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1.
C
C       N = Number of nodes in the triangulation.  N .GE. 6.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       Z = Array of length N containing data values associ-
C           ated with the nodes.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       DX,DY = Estimated partial derivatives at node K
C               unless IER < 0.
C
C       IER = Error indicator:
C             IER = L > 0 if no errors were encountered and
C                         L nodes (including node K) were
C                         employed in the least squares fit.
C             IER = -1 if K, NCC, an LCC entry, or N is
C                      outside its valid range on input.
C             IER = -2 if all nodes are collinear.
C
C TRIPACK modules required by GRADL:  GETNP, INTSEC
C
C SRFPACK modules required by GRADL:  GIVENS, ROTATE, SETRO1
C
C Intrinsic functions called by GRADL:  ABS, MIN, REAL, SQRT
C
C***********************************************************
C
      INTEGER   LMN, LMX
      PARAMETER (LMN=10,  LMX=30)
      INTEGER   I, IERR, J, JP1, KK, L, LMAX, LMIN, LM1,
     .          LNP, NP, NPTS(LMX)
      REAL      A(6,6), C, DIST(LMX), DMIN, DS, DTOL, RIN,
     .          RS, RTOL, S, SF, SFS, STF, SUM, W, XK, YK,
     .          ZK
      DATA      RTOL/1.E-5/, DTOL/.01/
C
C Local parameters:
C
C A =         Transpose of the augmented regression matrix
C C =         First component of the plane rotation deter-
C               mined by Subroutine GIVENS
C DIST =      Array containing the distances between K and
C               the elements of NPTS (refer to GETNP)
C DMIN =      Minimum of the magnitudes of the diagonal
C               elements of the regression matrix after
C               zeros are introduced below the diagonal
C DS =        Squared distance between nodes K and NPTS(LNP)
C DTOL =      Tolerance for detecting an ill-conditioned
C               system.  The system is accepted when DMIN/W
C               .GE. DTOL
C I =         DO-loop index
C IERR =      Error flag for calls to GETNP
C J =         DO-loop index
C JP1 =       J+1
C KK =        Local copy of K
C L =         Number of columns of A**T to which a rotation
C               is applied
C LMAX,LMIN = Min(LMX,N), Min(LMN,N)
C LMN,LMX =   Minimum and maximum values of LNP for N
C               sufficiently large.  In most cases LMN-1
C               nodes are used in the fit.  4 .LE. LMN .LE.
C               LMX.
C LM1 =       LMIN-1 or LNP-1
C LNP =       Length of NPTS
C NP =        Element of NPTS to be added to the system
C NPTS =      Array containing the indexes of a sequence of
C               nodes ordered by distance from K.  NPTS(1)=K
C               and the first LNP-1 elements of NPTS are
C               used in the least squares fit.  Unless LNP
C               exceeds LMAX, NPTS(LNP) determines R.
C RIN =       Inverse of the distance R between node K and
C               NPTS(LNP) or some point further from K than
C               NPTS(LMAX) if NPTS(LMAX) is used in the fit.
C               R is a radius of influence which enters into
C               the weight W.
C RS =        R*R
C RTOL =      Tolerance for determining R.  If the relative
C               change in DS between two elements of NPTS is
C               not greater than RTOL, they are treated as
C               being the same distance from node K
C S =         Second component of the plane rotation deter-
C               mined by Subroutine GIVENS
C SF =        Scale factor for the linear terms (columns 4
C               and 5) in the least squares fit -- inverse
C               of the root-mean-square distance between K
C               and the nodes (other than K) in the least
C               squares fit.
C SFS =       Scale factor for the quadratic terms (first 3
C               columns) in the least squares fit -- SF*SF.
C STF =       Marquardt stabilization factor used to damp
C               out the first 3 solution components (second
C               partials of the quadratic) when the system
C               is ill-conditioned.  As STF increases, the
C               fitting function approaches a linear
C SUM =       Sum of squared distances between node K and
C               the nodes used in the least squares fit
C W =         Weight associated with a row of the augmented
C               regression matrix -- 1/R - 1/D, where D < R
C               and D is the distance between K and a node
C               entering into the least squares fit.
C XK,YK,ZK =  Coordinates and data value associated with K
C
      KK = K
C
C Test for errors and initialize LMIN and LMAX.
C
      IF (KK .LT. 1  .OR.  KK .GT. N  .OR.  NCC .LT. 0
     .    .OR.  N .LT. 6) GO TO 13
      LMIN = MIN(LMN,N)
      LMAX = MIN(LMX,N)
C
C Compute NPTS, DIST, LNP, SF, SFS, and RIN --
C
C   Set NPTS to the closest LMIN-1 nodes to K.
C
      SUM = 0.
      NPTS(1) = KK
      DIST(1) = 0.
      LM1 = LMIN - 1
      DO 1 LNP = 2,LM1
        CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,
     .              LNP, NPTS,DIST, IERR)
        IF (IERR .NE. 0) GO TO 13
        DS = DIST(LNP)**2
        SUM = SUM + DS
    1   CONTINUE
C
C Add additional nodes to NPTS until the relative increase
C   in DS is at least RTOL.
C
      DO 3 LNP = LMIN,LMAX
        CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,
     .              LNP, NPTS,DIST, IERR)
        RS = DIST(LNP)**2
        IF ((RS-DS)/DS .LE. RTOL) GO TO 2
        IF (LNP .GT. 6) GO TO 4
    2   SUM = SUM + RS
    3   CONTINUE
C
C Use all LMAX nodes in the least squares fit.  RS is
C   arbitrarily increased by 10 per cent.
C
      RS = 1.1*RS
      LNP = LMAX + 1
C
C There are LNP-2 equations corresponding to nodes NPTS(2),
C   ...,NPTS(LNP-1).
C
    4 SFS = REAL(LNP-2)/SUM
      SF = SQRT(SFS)
      RIN = 1./SQRT(RS)
      XK = X(KK)
      YK = Y(KK)
      ZK = Z(KK)
C
C A Q-R decomposition is used to solve the least squares
C   system.  The transpose of the augmented regression
C   matrix is stored in A with columns (rows of A) defined
C   as follows:  1-3 are the quadratic terms, 4 and 5 are
C   the linear terms with coefficients DX and DY, and the
C   last column is the right hand side.
C
C Set up the first 5 equations and zero out the lower tri-
C   angle with Givens rotations.
C
      DO 6 I = 1,5
        NP = NPTS(I+1)
        W = 1./DIST(I+1) - RIN
        CALL SETRO1 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,
     .               W, A(1,I))
        IF (I .EQ. 1) GO TO 6
        DO 5 J = 1,I-1
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS (A(J,J),A(J,I),C,S)
          CALL ROTATE (L,C,S,A(JP1,J),A(JP1,I))
    5     CONTINUE
    6   CONTINUE
C
C Add the additional equations to the system using
C   the last column of A.  I .LE. LNP.
C
      I = 7
    7   IF (I .LT. LNP) THEN
          NP = NPTS(I)
          W = 1./DIST(I) - RIN
          CALL SETRO1 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,
     .                 W, A(1,6))
          DO 8 J = 1,5
            JP1 = J + 1
            L = 6 - J
            CALL GIVENS (A(J,J),A(J,6),C,S)
            CALL ROTATE (L,C,S,A(JP1,J),A(JP1,6))
    8       CONTINUE
          I = I + 1
          GO TO 7
        ENDIF
C
C Test the system for ill-conditioning.
C
      DMIN = MIN( ABS(A(1,1)),ABS(A(2,2)),ABS(A(3,3)),
     .            ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN/W .GE. DTOL) GO TO 12
      IF (LNP .LE. LMAX) THEN
C
C   Add another node to the system and increase R.  Note
C     that I = LNP.
C
        LNP = LNP + 1
        IF (LNP .LE. LMAX) THEN
          CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,
     .                LNP, NPTS,DIST, IERR)
          RS = DIST(LNP)**2
        ENDIF
        RIN = 1./SQRT(1.1*RS)
        GO TO 7
      ENDIF
C
C Stabilize the system by damping second partials -- add
C   multiples of the first three unit vectors to the first
C   three equations.
C
      STF = W
      DO 11 I = 1,3
        A(I,6) = STF
        DO 9 J = I+1,6
          A(J,6) = 0.
    9     CONTINUE
        DO 10 J = I,5
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS (A(J,J),A(J,6),C,S)
          CALL ROTATE (L,C,S,A(JP1,J),A(JP1,6))
   10     CONTINUE
   11   CONTINUE
C
C Test the damped system for ill-conditioning.
C
      DMIN = MIN( ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN/W .LT. DTOL) GO TO 14
C
C Solve the 2 by 2 triangular system for the partial
C   derivatives.
C
   12 DY = A(6,5)/A(5,5)
      DX = SF*(A(6,4) - A(5,4)*DY)/A(4,4)
      DY = SF*DY
      IER = LNP - 1
      RETURN
C
C Invalid input parameter.
C
   13 IER = -1
      RETURN
C
C No unique solution due to collinear nodes.
C
   14 IER = -2
      RETURN
      END
      SUBROUTINE GRCOEF (SIGMA,DCUB, D,SD)
      REAL SIGMA, DCUB, D, SD
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/18/96
C
C   This subroutine computes factors involved in the linear
C system solved by Subroutines GRADG and SMSGS.
C
C On input:
C
C       SIGMA = Nonnegative tension factor associated with a
C               triangulation arc.
C
C       DCUB = Cube of the positive arc length.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       D = Diagonal factor.  D = SIG*(SIG*COSHM(SIG) -
C           SINHM(SIG))/(E*DCUB), where E = SIG*SINH(SIG) -
C           2*COSHM(SIG).  D > 0.
C
C       SD = Off-diagonal factor.  SD = SIG*SINHM(SIG)/
C            (E*DCUB).  SD > 0.
C
C SRFPACK module required by GRCOEF:  SNHCSH
C
C Intrinsic function called by GRCOEF:  EXP
C
C***********************************************************
C
      REAL COSHM, COSHMM, E, EMS, SCM, SIG, SINHM, SSINH,
     .     SSM
C
      SIG = SIGMA
      IF (SIG .LT. 1.E-9) THEN
C
C SIG = 0:  cubic interpolant.
C
        D = 4./DCUB
        SD = 2./DCUB
      ELSEIF (SIG .LE. .5) THEN
C
C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
C                      cancellation error in the hyperbolic
C                      functions when SIGMA is small.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        E = (SIG*SINHM - COSHMM - COSHMM)*DCUB
        D = SIG*(SIG*COSHM-SINHM)/E
        SD = SIG*SINHM/E
      ELSE
C
C SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
C            to avoid overflow when SIGMA is large.
C
        EMS = EXP(-SIG)
        SSINH = 1. - EMS*EMS
        SSM = SSINH - 2.*SIG*EMS
        SCM = (1.-EMS)*(1.-EMS)
        E = (SIG*SSINH - SCM - SCM)*DCUB
        D = SIG*(SIG*SCM-SSM)/E
        SD = SIG*SSM/E
      ENDIF
      RETURN
      END
      SUBROUTINE INTRC0 (PX,PY,NCC,LCC,N,X,Y,Z,LIST,LPTR,
     .                   LEND, IST, PZ,IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        IST, IER
      REAL    PX, PY, X(N), Y(N), Z(N), PZ
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/03/98
C
C   Given a triangulation of a set of nodes in the plane,
C along with data values at the nodes, this subroutine com-
C putes the value at P = (PX,PY) of the piecewise linear
C function which interpolates the data values.  The surface
C is extended in a continuous fashion beyond the boundary of
C the triangulation, allowing extrapolation.
C
C On input:
C
C       PX,PY = Coordinates of the point P at which the sur-
C               face is to be evaluated.
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index of the
C             first node of constraint I in LCC(I).  For I =
C             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       Z = Array of length N containing data values at the
C           nodes.  Refer to Subroutine ZGRADL.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C The above parameters are not altered by this routine.
C
C       IST = Index of the starting node in the search for a
C             triangle containing P.  1 .LE. IST .LE. N.
C             The output value of IST from a previous call
C             may be a good choice.
C
C On output:
C
C       IST = Index of one of the vertices of the triangle
C             containing P (or a boundary node which is vis-
C             ible from P) unless IER < 0.
C
C       PZ = Value of the interpolatory surface at P, or
C            zero if IER < 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and P is
C                     contained in a triangle but not in a
C                     constraint region.
C             IER = 1 if no errors were encountered and P
C                     lies in a constraint region triangle.
C                     PZ is effectively an extrapolated
C                     value in this case.
C             IER = 2 if no errors were encountered and P is
C                     exterior to the triangulation.  PZ is
C                     an extrapolated value in this case.
C             IER = -1 if NCC, N, or IST is outside its
C                      valid range on input.  LCC is not
C                      tested for validity.
C             IER = -2 if the nodes are collinear.
C
C TRIPACK modules required by INTRC0:  CRTRI, JRAND, LEFT,
C                                        LSTPTR, TRFIND
C
C SRFPACK module required by INTRC0:  COORDS
C
C***********************************************************
C
      LOGICAL CRTRI
      INTEGER I1, I2, I3, IERR, LPL, N1, N2
      REAL    B1, B2, B3, DP, X1, X2, XP, Y1, Y2, YP
C
      XP = PX
      YP = PY
      PZ = 0.
C
C Test for invalid input parameters.
C
      IF (NCC .LT. 0  .OR.  N .LT. 3  .OR.  IST .LT. 1
     .    .OR.  IST .GT. N) THEN
        IER = -1
        RETURN
      ENDIF
C
C Find a triangle (I1,I2,I3) containing P, or a pair of
C   visible boundary nodes I1 and I2.
C
      CALL TRFIND (IST,XP,YP,N,X,Y,LIST,LPTR,LEND, I1,I2,I3)
      IF (I1 .EQ. 0) THEN
        IER = -2
        RETURN
      ENDIF
      IST = I1
      IF (I3 .EQ. 0) GO TO 1
C
C P is in a triangle.  Compute its barycentric coordinates.
C
      CALL COORDS (XP,YP,X(I1),X(I2),X(I3),Y(I1),Y(I2),
     .             Y(I3), B1,B2,B3,IERR)
      IF (IERR .NE. 0) THEN
        IER = -2
        RETURN
      ENDIF
C
C Compute an interpolated value.
C
      PZ = B1*Z(I1) + B2*Z(I2) + B3*Z(I3)
      IER = 0
C
      IF (CRTRI(NCC,LCC,I1,I2,I3)) THEN
        IER = 1
      ELSE
        IER = 0
      ENDIF
      RETURN
C
C P is exterior to the triangulation.  Extrapolate to P by
C   extending the interpolatory surface as a constant
C   beyond the boundary:  PZ is the function value at Q
C   where Q is the closest boundary point to P.
C
C Determine Q by traversing the boundary starting from the
C   rightmost visible node I1.
C
    1 IER = 2
      N2 = I1
C
C Top of loop:
C
C   Set N1 to the last neighbor of N2, and compute the dot
C     product DP = (N2->N1,N2->P).  P FORWARD N2->N1 iff
C     DP > 0.
C
    2 LPL = LEND(N2)
      N1 = -LIST(LPL)
      X1 = X(N1)
      Y1 = Y(N1)
      X2 = X(N2)
      Y2 = Y(N2)
      DP = (X1-X2)*(XP-X2) + (Y1-Y2)*(YP-Y2)
      IF (DP .LE. 0.) THEN
C
C   N2 is the closest boundary point to P.
C
        PZ = Z(N2)
        RETURN
      ENDIF
C
C   P FORWARD N2->N1.  Test for P FORWARD N1->N2.
C
      IF ((XP-X1)*(X2-X1) + (YP-Y1)*(Y2-Y1) .GT. 0.) THEN
C
C   The closest boundary point to P lies on N2-N1.  Compute
C     its local coordinates with respect to N2-N1.
C
        B1 = DP/( (X2-X1)**2 + (Y2-Y1)**2 )
        B2 = 1. - B1
        PZ = B1*Z(N1) + B2*Z(N2)
        RETURN
      ENDIF
C
C   Bottom of boundary traversal loop.
C
      N2 = N1
      GO TO 2
      END
      SUBROUTINE INTRC1 (PX,PY,NCC,LCC,N,X,Y,Z,LIST,LPTR,
     .                   LEND,IFLGS,SIGMA,GRAD,
     .                   DFLAG, IST, PZ,PZX,PZY,IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        IFLGS, IST, IER
      LOGICAL DFLAG
      REAL    PX, PY, X(N), Y(N), Z(N), SIGMA(*), GRAD(2,N),
     .        PZ, PZX, PZY
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/03/98
C
C   Given a triangulation of a set of nodes in the plane,
C along with data values and estimated gradients at the
C nodes, this subroutine computes the value and, optionally,
C the first partial derivatives at P = (PX,PY) of a C-1
C (once-continuously differentiable) function F which int-
C erpolates the data values and gradients.  Extrapolation to
C a point exterior to the triangulation is accomplished by
C extending the surface in such a way that F is C-1 over the
C entire plane.
C
C   Subroutine FVAL is used to evaluate an interpolatory
C surface under tension, while Subroutine TVAL is called in
C the case of no tension (IFLGS .LE. 0 and SIGMA(1) = 0).
C However, the surface under tension is well-defined with
C SIGMA = 0, and, in this case, the two interpolants are
C identical on triangulation arcs and outside the triangula-
C tion.  The use of FVAL with no tension can be forced (at
C a cost in efficiency) by setting IFLGS = 1 and storing
C zeros in all components of SIGMA.  Note, however, that
C first partial derivatives are only available from TVAL
C (and at points outside the triangulation);  i.e., a proce-
C dure for differentiating the surface under tension has not
C been implemented.
C
C   A set of interpolated values at the vertices of a rec-
C tangular grid can be obtained by a single call to
C Subroutine UNIF.  Subroutine INTRC0 provides for evalua-
C tion of the piecewise linear interpolatory surface.
C
C On input:
C
C       PX,PY = Coordinates of the point P at which the sur-
C               face is to be evaluated.
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index of the
C             first node of constraint I in LCC(I).  For I =
C             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       Z = Array of length N containing data values at the
C           nodes.  Refer to Subroutines ZGRADG and ZGRADL.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C       IFLGS = Tension factor option:
C               IFLGS .LE. 0 if a single uniform tension
C                            factor is to be used.
C               IFLGS .GE. 1 if variable tension is desired.
C
C       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               routines FVAL, GETSIG, SIG0, SIG1, and SIG2.
C
C       GRAD = 2 by N array whose columns contain estimated
C              gradients at the nodes with X partial deriva-
C              tives in the first row and Y partials in the
C              second.  Refer to Subroutines GRADC, GRADG,
C              GRADL, SMSURF, ZGRADG, and ZGRADL.
C
C       DFLAG = Logical flag which specifies whether first
C               partial derivatives at P are to be computed:
C               DFLAG = TRUE if and only if partials are
C               to be computed by TVAL.  This option is only
C               valid for IFLGS .LE. 0 and SIGMA(1) = 0 (and
C               for points outside the triangulation).
C
C The above parameters are not altered by this routine.
C
C       IST = Index of the starting node in the search for a
C             triangle containing P.  1 .LE. IST .LE. N.
C             The output value of IST from a previous call
C             may be a good choice.
C
C On output:
C
C       IST = Index of one of the vertices of the triangle
C             containing P (or a boundary node which is vis-
C             ible from P) unless IER = -1 or IER = -2.
C
C       PZ = Value of the interpolatory surface at P, or
C            zero if IER < 0.
C
C       PZX,PZY = X and Y partials at P if DFLAG = .TRUE.
C                 and IER .GE. 0, unaltered otherwise.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and P is
C                     contained in a triangle but not in a
C                     constraint region.
C             IER = 1 if no errors were encountered and P
C                     lies in a constraint region triangle.
C                     PZ is effectively an extrapolated
C                     value in this case.
C             IER = 2 if no errors were encountered and P is
C                     exterior to the triangulation.  PZ is
C                     an extrapolated value in this case.
C             IER = -1 if NCC, N, or IST is outside its
C                      valid range on input.  LCC is not
C                      tested for validity.
C             IER = -2 if the nodes are collinear.
C             IER = -3 if P is contained in a triangle and
C                      DFLAG = TRUE, but IFLGS > 0 or
C                      SIGMA(1) .NE. 0.
C
C TRIPACK modules required by INTRC1:  CRTRI, JRAND, LEFT,
C                                        LSTPTR, TRFIND
C
C SRFPACK modules required by INTRC1:  ARCINT, COORDS, FVAL,
C                                        SNHCSH, TVAL
C
C Intrinsic function called by INTRC1:  SQRT
C
C***********************************************************
C
      INTEGER LSTPTR
      LOGICAL CRTRI
      INTEGER I1, I2, I3, IERR, LP, LPL, N1, N2, N3
      LOGICAL TENSN
      REAL    A1, A2, B1, B2, C1, C2, D, D1, D2, D3, DP,
     .        DP1, DP3, F1, F2, R1, R12, R2, SIG, SIG1,
     .        SIG2, SIG3, T, T1, T2, X1, X2, X3, X12, X23,
     .        X2P, XP, XQ, XQP, Y1, Y2, Y3, Y12, Y23, Y2P,
     .        YP, YQ, YQP, Z1, Z2, Z3, ZQ, ZX1, ZX2, ZX3,
     .        ZXQ, ZY1, ZY2, ZY3, ZYQ
C
      XP = PX
      YP = PY
      PZ = 0.
C
C Test for invalid input parameters.
C
      IF (NCC .LT. 0  .OR.  N .LT. 3  .OR.  IST .LT. 1
     .    .OR.  IST .GT. N) THEN
        IER = -1
        RETURN
      ENDIF
C
C Find a triangle (I1,I2,I3) containing P, or a pair of
C   visible boundary nodes I1 and I2.
C
      CALL TRFIND (IST,XP,YP,N,X,Y,LIST,LPTR,LEND, I1,I2,I3)
      IF (I1 .EQ. 0) THEN
        IER = -2
        RETURN
      ENDIF
      IST = I1
      TENSN = IFLGS .GE. 1  .OR.  SIGMA(1) .NE. 0
      IF (I3 .EQ. 0) GO TO 1
      IF (DFLAG  .AND.  TENSN) THEN
        IER = -3
        RETURN
      ENDIF
C
C P is in a triangle.  Store local parameters for the
C   call to FVAL or TVAL.
C
      X1 = X(I1)
      Y1 = Y(I1)
      X2 = X(I2)
      Y2 = Y(I2)
      X3 = X(I3)
      Y3 = Y(I3)
      Z1 = Z(I1)
      Z2 = Z(I2)
      Z3 = Z(I3)
      ZX1 = GRAD(1,I1)
      ZX2 = GRAD(1,I2)
      ZX3 = GRAD(1,I3)
      ZY1 = GRAD(2,I1)
      ZY2 = GRAD(2,I2)
      ZY3 = GRAD(2,I3)
      IF (TENSN) THEN
C
C Set SIG1, SIG2, and SIG3 to the tension factors associated
C   with the sides opposite I1, I2, and I3, respectively,
C   and compute a value from FVAL.
C
        IF (IFLGS .LE. 0) THEN
          SIG1 = SIGMA(1)
          SIG2 = SIG1
          SIG3 = SIG1
        ELSE
          LP = LSTPTR(LEND(I2),I3,LIST,LPTR)
          SIG1 = SIGMA(LP)
          LP = LSTPTR(LEND(I3),I1,LIST,LPTR)
          SIG2 = SIGMA(LP)
          LP = LSTPTR(LEND(I1),I2,LIST,LPTR)
          SIG3 = SIGMA(LP)
        ENDIF
        CALL FVAL (XP,YP,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,ZX1,ZX2,
     .             ZX3,ZY1,ZY2,ZY3,SIG1,SIG2,SIG3, PZ,IERR)
        IF (IERR .LT. 0) THEN
          IER = -2
          RETURN
        ENDIF
      ELSE
C
C Compute an interpolated value from TVAL for no tension.
C
        CALL TVAL (XP,YP,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,ZX1,ZX2,
     .             ZX3,ZY1,ZY2,ZY3,DFLAG, PZ,PZX,PZY,IERR)
        IF (IERR .NE. 0) THEN
          IER = -2
          RETURN
        ENDIF
      ENDIF
C
      IF (CRTRI(NCC,LCC,I1,I2,I3)) THEN
        IER = 1
      ELSE
        IER = 0
      ENDIF
      RETURN
C
C P is exterior to the triangulation.  Extrapolate to P by
C   passing a linear function of one variable through the
C   value and directional derivative (in the direction
C   Q->P) of the interpolatory surface F at Q, where Q is
C   the closest boundary point to P.
C
C Determine Q by traversing the boundary starting from the
C   rightmost visible node I1.
C
    1 IER = 2
      N2 = I1
C
C Top of loop:
C
C   Set N1 to the last neighbor of N2, and compute the dot
C     product DP = (N2->N1,N2->P).  P FORWARD N2->N1 iff
C     DP > 0.
C
    2 LPL = LEND(N2)
      N1 = -LIST(LPL)
      X1 = X(N1)
      Y1 = Y(N1)
      X2 = X(N2)
      Y2 = Y(N2)
      DP = (X1-X2)*(XP-X2) + (Y1-Y2)*(YP-Y2)
      IF (DP .LE. 0.) THEN
C
C   N2 is the closest boundary point to P:  P lies in a
C     wedge with sides orthogonal to N1-N2 and N2-N3, where
C     N3 is the first neighbor of N2.  The linear interpo-
C     lant must be modified by a correction term which
C     provides for continuity of the derivative across the
C     sides of the wedge.
C
        LP = LPTR(LPL)
        N3 = LIST(LP)
        ZX1 = GRAD(1,N1)
        ZX2 = GRAD(1,N2)
        ZX3 = GRAD(1,N3)
        ZY1 = GRAD(2,N1)
        ZY2 = GRAD(2,N2)
        ZY3 = GRAD(2,N3)
        X12 = X2-X1
        Y12 = Y2-Y1
        X23 = X(N3)-X2
        Y23 = Y(N3)-Y2
        X2P = XP-X2
        Y2P = YP-Y2
        DP1 = -DP
        DP3 = X23*X2P + Y23*Y2P
        D2 = X2P*X2P + Y2P*Y2P
        D1 = SQRT((X12*X12 + Y12*Y12)*D2)
        D3 = SQRT((X23*X23 + Y23*Y23)*D2)
        D = (X12*Y23 - Y12*X23)**2
        T1 = DP3*( X12*(ZY2-ZY1)-Y12*(ZX2-ZX1) )/D1
        T2 = DP1*( X23*(ZY3-ZY2)-Y23*(ZX3-ZX2) )/D3
        T = DP1*DP3*(T1+T2)
        PZ = Z(N2) + X2P*ZX2 + Y2P*ZY2
        IF (D .NE. 0.) PZ = PZ - T/D
        IF (DFLAG) THEN
          T = T/D2
          D1 = DP3*(T1+T2+T2)
          D2 = DP1*(T1+T1+T2)
          PZX = ZX2
          IF (D .NE. 0.) PZX = PZX + (T*X2P - D1*X12 -
     .                                D2*X23)/D
          PZY = ZY2
          IF (D .NE. 0.) PZY = PZY + (T*Y2P - D1*Y12 -
     .                                D2*Y23)/D
        ENDIF
        RETURN
      ENDIF
C
C   P FORWARD N2->N1.  Test for P FORWARD N1->N2.
C
      IF ((XP-X1)*(X2-X1) + (YP-Y1)*(Y2-Y1) .LE. 0.) THEN
C
C   Bottom of boundary traversal loop.
C
        N2 = N1
        GO TO 2
      ENDIF
C
C The closest boundary point Q lies on N2-N1.  Store par-
C   tials at N1 and N2, and compute Q and its barycentric
C   coordinates R1 and R2.
C
      ZX1 = GRAD(1,N1)
      ZY1 = GRAD(2,N1)
      ZX2 = GRAD(1,N2)
      ZY2 = GRAD(2,N2)
      X12 = X2-X1
      Y12 = Y2-Y1
      D2 = X12*X12 + Y12*Y12
      R1 = DP/D2
      R2 = 1. - R1
      XQ = R1*X1 + R2*X2
      YQ = R1*Y1 + R2*Y2
      IF (TENSN) THEN
C
C   Set SIG to the tension factor associated with N1-N2 and
C     compute an interpolated value ZQ at Q from FVAL.
C
        IF (IFLGS .LE. 0) THEN
          SIG = SIGMA(1)
        ELSE
          SIG = SIGMA(LPL)
        ENDIF
        CALL ARCINT (R1,X1,X2,Y1,Y2,Z(N1),Z(N2),ZX1,ZX2,
     .               ZY1,ZY2,SIG,.TRUE., ZQ,ZXQ,ZYQ,IERR)
C
C   Compute the extrapolated value at P.
C
        XQP = XP-XQ
        YQP = YP-YQ
        PZ = ZQ + ZXQ*XQP + ZYQ*YQP
        IF (DFLAG) THEN
          T = ((ZX2-ZX1)*XQP + (ZY2-ZY1)*YQP)/D2
          PZX = ZXQ + X12*T
          PZY = ZYQ + Y12*T
        ENDIF
      ELSE
C
C   Compute the cardinal function values and interpolated
C     value at Q associated with TVAL.
C
        R12 = R1*R2
        F1 = R1*R12
        F2 = R2*R12
        A1 = R1 + (F1-F2)
        A2 = R2 - (F1-F2)
        B1 = X12*F1
        B2 = -X12*F2
        C1 = Y12*F1
        C2 = -Y12*F2
        ZQ = A1*Z(N1) + A2*Z(N2) + B1*ZX1 + B2*ZX2 +
     .       C1*ZY1 + C2*ZY2
C
C   Compute the extrapolated value at P.
C
        XQP = XP-XQ
        YQP = YP-YQ
        T1 = R1*ZX1 + R2*ZX2
        T2 = R1*ZY1 + R2*ZY2
        PZ = ZQ + T1*XQP + T2*YQP
        IF (DFLAG) THEN
          T = (3.*R12*(2.*(Z(N2)-Z(N1)) - X12*(ZX1+ZX2) -
     .         Y12*(ZY1+ZY2)) + (ZX2-ZX1)*XQP +
     .         (ZY2-ZY1)*YQP)/D2
          PZX = T1 + X12*T
          PZY = T2 + Y12*T
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE ROTATE (N,C,S, X,Y )
      INTEGER N
      REAL    C, S, X(N), Y(N)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C                                                ( C  S)
C   This subroutine applies the Givens rotation  (     )  to
C                                                (-S  C)
C                    (X(1) ... X(N))
C the 2 by N matrix  (             ) .
C                    (Y(1) ... Y(N))
C
C   This routine is identical to Subroutine SROT from the
C LINPACK BLAS (Basic Linear Algebra Subroutines).
C
C On input:
C
C       N = Number of columns to be rotated.
C
C       C,S = Elements of the Givens rotation.  Refer to
C             Subroutine GIVENS.
C
C The above parameters are not altered by this routine.
C
C       X,Y = Arrays of length .GE. N containing the compo-
C             nents of the vectors to be rotated.
C
C On output:
C
C       X,Y = Arrays containing the rotated vectors (not
C             altered if N < 1).
C
C Modules required by ROTATE:  None
C
C***********************************************************
C
      INTEGER I
      REAL    XI, YI
C
      DO 1 I = 1,N
        XI = X(I)
        YI = Y(I)
        X(I) = C*XI + S*YI
        Y(I) = -S*XI + C*YI
    1   CONTINUE
      RETURN
      END
      SUBROUTINE SETRO1 (XK,YK,ZK,XI,YI,ZI,S1,S2,W, ROW)
      REAL XK, YK, ZK, XI, YI, ZI, S1, S2, W, ROW(6)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This subroutine sets up the I-th row of an augmented re-
C gression matrix for a weighted least squares fit of a
C quadratic function Q(X,Y) to a set of data values Z, where
C Q(XK,YK) = ZK.  The first three columns (quadratic terms)
C are scaled by S2, and the fourth and fifth columns (lin-
C ear terms) are scaled by S1.
C
C On input:
C
C       XK,YK = Coordinates of node K.
C
C       ZK = Data value at node K to be interpolated by Q.
C
C       XI,YI,ZI = Coordinates and data value at node I.
C
C       S1,S2 = Scale factors.
C
C       W = Weight associated with node I.
C
C The above parameters are not altered by this routine.
C
C       ROW = Array of length 6.
C
C On output:
C
C       ROW = Array containing a row of the augmented re-
C             gression matrix.
C
C Modules required by SETRO1:  None
C
C***********************************************************
C
      REAL DX, DY, W1, W2
C
      DX = XI - XK
      DY = YI - YK
      W1 = S1*W
      W2 = S2*W
      ROW(1) = DX*DX*W2
      ROW(2) = DX*DY*W2
      ROW(3) = DY*DY*W2
      ROW(4) = DX*W1
      ROW(5) = DY*W1
      ROW(6) = (ZI - ZK)*W
      RETURN
      END
      SUBROUTINE SETRO2 (XK,YK,ZK,XI,YI,ZI,S1,S2,W, ROW)
      REAL XK, YK, ZK, XI, YI, ZI, S1, S2, W, ROW(7)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This subroutine sets up the I-th row of an augmented re-
C gression matrix for a weighted least squares fit of a
C quadratic function Q(X,Y) to a set of data values Z.  The
C first 3 columns (quadratic terms) are scaled by S2, and
C the fourth and fifth columns (linear terms) are scaled by
C S1.
C
C On input:
C
C       XK,YK = Coordinates of node K.
C
C       ZK = Data value at node K to be interpolated by Q
C            (Q(XK,YK) = ZK) if the constant term, ROW(6),
C            is to be ignored, or zero if Q(XK,YK) is to be
C            a parameter (coefficient of ROW(6)).
C
C       XI,YI,ZI = Coordinates and data value at node I.
C
C       S1,S2 = Scale factors.
C
C       W = Weight associated with node I.
C
C The above parameters are not altered by this routine.
C
C       ROW = Array of length 7.
C
C On output:
C
C       ROW = Array containing a row of the augmented re-
C             gression matrix.
C
C Modules required by SETRO2:  None
C
C***********************************************************
C
      REAL DX, DY, W1, W2
C
      DX = XI - XK
      DY = YI - YK
      W1 = S1*W
      W2 = S2*W
      ROW(1) = DX*DX*W2
      ROW(2) = DX*DY*W2
      ROW(3) = DY*DY*W2
      ROW(4) = DX*W1
      ROW(5) = DY*W1
      ROW(6) = W
      ROW(7) = (ZI - ZK)*W
      RETURN
      END
      SUBROUTINE SETRO3 (XK,YK,ZK,XI,YI,ZI,S1,S2,S3,W, ROW)
      REAL XK, YK, ZK, XI, YI, ZI, S1, S2, S3, W, ROW(10)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   01/25/97
C
C   This subroutine sets up the I-th row of an augmented re-
C gression matrix for a weighted least squares fit of a
C cubic function f(x,y) to a set of data values z, where
C f(XK,YK) = ZK.  The first four columns (cubic terms) are
C scaled by S3, the next three columns (quadratic terms)
C are scaled by S2, and the eighth and ninth columns (lin-
C ear terms) are scaled by S1.
C
C On input:
C
C       XK,YK = Coordinates of node K.
C
C       ZK = Data value at node K to be interpolated by f.
C
C       XI,YI,ZI = Coordinates and data value at node I.
C
C       S1,S2,S3 = Scale factors.
C
C       W = Weight associated with node I.
C
C The above parameters are not altered by this routine.
C
C       ROW = Array of length 10.
C
C On output:
C
C       ROW = Array containing a row of the augmented re-
C             gression matrix.
C
C Modules required by SETRO3:  None
C
C***********************************************************
C
      REAL DX, DY, W1, W2, W3
C
      DX = XI - XK
      DY = YI - YK
      W1 = S1*W
      W2 = S2*W
      W3 = S3*W
      ROW(1) = DX*DX*DX*W3
      ROW(2) = DX*DX*DY*W3
      ROW(3) = DX*DY*DY*W3
      ROW(4) = DY*DY*DY*W3
      ROW(5) = DX*DX*W2
      ROW(6) = DX*DY*W2
      ROW(7) = DY*DY*W2
      ROW(8) = DX*W1
      ROW(9) = DY*W1
      ROW(10) = (ZI - ZK)*W
      RETURN
      END
      SUBROUTINE SGPRNT (N,LUNIT,LIST,LPTR,LEND,SIGMA)
      INTEGER N, LUNIT, LIST(*), LPTR(*), LEND(N)
      REAL    SIGMA(*)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/03/98
C
C   Given a triangulation of a set of nodes in the plane,
C along with an array of tension factors associated with the
C triangulation arcs, this subroutine prints the list of
C arcs (with tension factors) ordered by endpoint nodal in-
C dexes.  An arc is identified with its smaller endpoint
C index:  N1-N2, where N1 < N2.
C
C On input:
C
C       N = Number of nodes in the triangulation.  3 .LE. N
C           .LE. 9999.
C
C       LUNIT = Logical unit for output.  0 .LE. LUNIT .LE.
C               99.  Output is printed on unit 6 if LUNIT is
C               outside its valid range.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C       SIGMA = Array of length 2*NA = 6*(N-1)-2*NB, where
C               NA and NB are the numbers of arcs and boun-
C               dary nodes, respectively, containing tension
C               factors associated with arcs in one-to-one
C               correspondence with LIST entries.  Note that
C               each arc N1-N2 has two LIST entries and
C               thus, SIGMA(I) and SIGMA(J) should be iden-
C               tical, where LIST(I) = N2 (in the adjacency
C               list for N1) and LIST(J) = N1 (in the list
C               associated with N2).  Both SIGMA(I) and
C               SIGMA(J) are printed if they are not iden-
C               tical.
C
C None of the parameters are altered by this routine.
C
C TRIPACK module required by SGPRNT:  LSTPTR
C
C Intrinsic function called by SGPRNT:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER LP1, LP2, LPL, LUN, N1, N2, NA, NAT, NB, NE,
     .        NL, NLMAX, NM1, NMAX
      LOGICAL ERROR
      REAL    SIG
      DATA NMAX/9999/,  NLMAX/60/
C
      LUN = LUNIT
      IF (LUN .LT. 0  .OR.  LUN .GT. 99) LUN = 6
C
C Print a heading, test for invalid N, and initialize coun-
C   ters:
C
C NL = Number of lines printed on the current page
C NA = Number of arcs encountered
C NE = Number of errors in SIGMA encountered
C NB = Number of boundary nodes encountered
C
      WRITE (LUN,100) N
      IF (N .LT. 3  .OR.  N .GT. NMAX) GO TO 4
      NL = 6
      NA = 0
      NE = 0
      NB = 0
C
C Outer loop on nodes N1.  LPL points to the last neighbor
C   of N1.
C
      NM1 = N - 1
      DO 3 N1 = 1,NM1
        LPL = LEND(N1)
        IF (LIST(LPL) .LT. 0) NB = NB + 1
        LP1 = LPL
C
C Inner loop on neighbors N2 of N1 such that N1 < N2.
C
    1   LP1 = LPTR(LP1)
          N2 = ABS(LIST(LP1))
          IF (N2 .LT. N1) GO TO 2
          NA = NA + 1
          SIG = SIGMA(LP1)
C
C   Test for an invalid SIGMA entry.
C
          LP2 = LSTPTR (LEND(N2),N1,LIST,LPTR)
          ERROR = SIGMA(LP2) .NE. SIG
          IF (ERROR) NE = NE + 1
C
C   Print a line and update the counters.
C
          IF (.NOT. ERROR) WRITE (LUN,110) N1, N2, SIG
          IF (ERROR) WRITE (LUN,120) N1, N2, SIG, SIGMA(LP2)
          NL = NL + 1
          IF (NL .GE. NLMAX) THEN
            WRITE (LUN,130)
            NL = 1
          ENDIF
C
C Bottom of loop on neighbors N2 of N1.
C
    2     IF (LP1 .NE. LPL) GO TO 1
    3   CONTINUE
      LPL = LEND(N)
      IF (LIST(LPL) .LT. 0) NB = NB + 1
C
C Test for errors in SIGMA.
C
      IF (NE .GT. 0) WRITE (LUN,200) NE
C
C Print NA and test for an invalid triangulation.
C
      WRITE (LUN,140) NA
      NAT = 3*NM1 - NB
      IF (NAT .NE. NA) WRITE (LUN,210) NAT
      RETURN
C
C N is outside its valid range.
C
    4 WRITE (LUN,220) NMAX
      RETURN
C
C Print formats:
C
  100 FORMAT (///14X,'Tension Factors,  N =',I5,
     .        ' Nodes'//1X,18X,'N1',5X,'N2',8X,'Tension'//)
  110 FORMAT (1X,16X,I4,3X,I4,5X,F12.8)
  120 FORMAT (1X,16X,I4,3X,I4,5X,F12.8,3X,F12.8,' *')
  130 FORMAT (///)
  140 FORMAT (//1X,10X,'NA =',I5,' Arcs')
C
C Error messages:
C
  200 FORMAT (//1X,10X,'*',I5,' Errors in SIGMA')
  210 FORMAT (/1X,10X,'*** Error in triangulation:  ',
     .        '3N-NB-3 = ',I5,' ***')
  220 FORMAT (1X,10X,'*** N is outside its valid range:  ',
     .        'NMAX = ',I4,' ***')
      END
      REAL FUNCTION SIG0 (N1,N2,N,X,Y,H,LIST,LPTR,LEND,HXHY,
     .                    IFLGB,HBND,TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGB,
     .        IFLGS, IER
      REAL    X(N), Y(N), H(N), HXHY(2,N), HBND, TOL,
     .        SIGMA(*)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/03/98
C
C   Given a triangulation of a set of nodes in the plane,
C along with data values H and gradients (HX,HY) at the
C nodes, this function determines the smallest tension fac-
C tor SIG0 such that the Hermite interpolatory tension
C spline H(T), defined by SIG0 and the endpoint values and
C directional derivatives associated with an arc N1-N2, is
C bounded (either above or below) by HBND for all T in
C (T1,T2), where (T1,T2) denotes an interval corresponding
C to the arc.
C
C On input:
C
C       N1,N2 = Nodal indexes of the endpoints of an arc for
C               which the tension factor is to be computed.
C               The indexes must be distinct and lie in the
C               range 1 to N, and if IFLGS .GE. 1, they must
C               correspond to adjacent nodes in the triangu-
C               lation.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I)).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C       HXHY = Array dimensioned 2 by N whose columns con-
C              partial derivatives at the nodes (X partials
C              in the first row).  Refer to Subroutines
C              GRADC, GRADG, GRADL, and SMSURF.
C
C       IFLGB = Bound option indicator:
C               IFLGB = -1 if HBND is a lower bound on H.
C               IFLGB = 1 if HBND is an upper bound on H.
C
C       HBND = Bound on H.  HBND .LE. min(H1,H2) if IFLGB =
C              -1 and HBND .GE. max(H1,H2) if IFLGB = 1,
C              where H1 and H2 are the data values at the
C              endpoints of the arc N1-N2.
C
C       TOL = Tolerance whose magnitude determines how close
C             SIG0 is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy the constraint.  For a lower bound,
C             SIG0 is chosen so that HBND .LE. HMIN .LE.
C             HBND + abs(TOL), where HMIN is the minimum
C             value of H on the arc, and for an upper bound,
C             the maximum of H satisfies HBND - abs(TOL)
C             .LE. HMAX .LE. HBND.  Thus, the constraint is
C             satisfied but possibly with more tension than
C             necessary.
C
C       IFLGS = Tension array option indicator:
C               IFLGS .LE. 0 if SIGMA is not to be used.
C               IFLGS .GE. 1 if SIGMA is to be updated by
C                            storing SIG0 in the appropriate
C                            locations.
C
C The above parameters are not altered by this function.
C
C       SIGMA = Dummy array of length 1 (IFLGS .LE. 0) or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               routine GETSIG.
C
C On output:
C
C       SIGMA = Tension factor array updated with the new
C               value if and only if IFLGS .GE. 1 and IER
C               .GE. 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     constraint can be satisfied with fin-
C                     ite tension.
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     the constraint (e.g., IFLGB = -1, HBND
C                     = H(T1), and the derivative of H at T1
C                     is negative).
C             IER = -1 if N1, N2, N, or IFLGB is outside its
C                      valid range.
C             IER = -2 if nodes N1 and N2 coincide or IFLGS
C                      .GE. 1 and the nodes are not adja-
C                      cent.
C             IER = -3 if HBND is outside its valid range.
C
C       SIG0 = Minimum tension factor defined above unless
C              IER < 0, in which case SIG0 = -1.  If IER
C              = 1, SIG0 is set to 85, resulting in an
C              approximation to the linear interpolant of
C              the endpoint values.
C
C TRIPACK module required by SIG0:  STORE
C
C SRFPACK module required by SIG0:  SNHCSH
C
C Intrinsic functions called by SIG0:  ABS, EXP, LOG, MAX,
C                                        MIN, REAL, SIGN,
C                                        SQRT
C
C***********************************************************
C
      REAL    STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    A, AA, A0, B, B0, BND, C, C1, C2, COSHM,
     .        COSHMM, D, D0, D1PD2, D2, DMAX, DSIG, DX, DY,
     .        E, EMS, F, F0, FMAX, FNEG, FTOL, H1, H2, R,
     .        RF, RSIG, RTOL, S, S1, S2, SBIG, SCM, SIG,
     .        SINHM, SNEG, SSINH, SSM, STOL, T, T0, T1, T2,
     .        TM
      DATA SBIG/85./,  LUN/-1/
      RF = REAL(IFLGB)
      BND = HBND
C
C Print a heading.
C
      IF (LUN .GE. 0) THEN
        IF (RF .LT. 0.) WRITE (LUN,100) N1, N2, BND
        IF (RF .GT. 0.) WRITE (LUN,110) N1, N2, BND
      ENDIF
  100 FORMAT (//1X,'SIG0 -- N1 =',I4,', N2 =',I4,
     .        ', Lower bound = ',E15.8)
  110 FORMAT (//1X,'SIG0 -- N1 =',I4,', N2 =',I4,
     .        ', Upper bound = ',E15.8)
C
C Test for errors and store local parameters.
C
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.
     .    MAX(N1,N2,3) .GT. N  .OR.  ABS(RF) .NE. 1.)
     .   GO TO 11
      IER = -2
      IF (IFLGS .GT. 0) THEN
C
C   Set LP1 and LP2 to the pointers to N2 as a neighbor of
C     N1 and N1 as a neighbor of N2, respectively.
C
        LPL = LEND(N1)
        LP1 = LPTR(LPL)
    1   IF (LIST(LP1) .EQ. N2) GO TO 2
          LP1 = LPTR(LP1)
          IF (LP1 .NE. LPL) GO TO 1
        IF (ABS(LIST(LP1)) .NE. N2) GO TO 11
C
    2   LPL = LEND(N2)
        LP2 = LPTR(LPL)
    3   IF (LIST(LP2) .EQ. N1) GO TO 4
          LP2 = LPTR(LP2)
          IF (LP2 .NE. LPL) GO TO 3
        IF (ABS(LIST(LP2)) .NE. N1) GO TO 11
      ENDIF
C
C Test for arc length DT = SQRT(DX**2+DY**2) = 0.
C
    4 DX = X(N2) - X(N1)
      DY = Y(N2) - Y(N1)
      IF (DX .EQ. 0.  .AND.  DY .EQ. 0.) GO TO 11
C
C Store endpoint values and test for a valid constraint.
C
      H1 = H(N1)
      H2 = H(N2)
      IER = -3
      IF ((RF .LT. 0.  .AND.  MIN(H1,H2) .LT. BND)  .OR.
     .    (RF .GT. 0.  .AND.  BND .LT. MAX(H1,H2)))
     .   GO TO 11
C
C Compute scaled directional derivatives S1,S2 at the end-
C   points (for the direction N1->N2) and test for infinite
C   tension required.
C
      S1 = HXHY(1,N1)*DX + HXHY(2,N1)*DY
      S2 = HXHY(1,N2)*DX + HXHY(2,N2)*DY
      IER = 1
      SIG = SBIG
      IF ((H1 .EQ. BND  .AND.  RF*S1 .GT. 0.)  .OR.
     .    (H2 .EQ. BND  .AND.  RF*S2 .LT. 0.)) GO TO 10
C
C Test for SIG = 0 sufficient.
C
      IER = 0
      SIG = 0.
      IF (RF*S1 .LE. 0.  .AND.  RF*S2 .GE. 0.) GO TO 10
C
C   Compute first difference S and coefficients A0 and B0
C     of the Hermite cubic interpolant H0(T) = H2 - (S2*R +
C     B0*R**2 + A0*R**3), where R = (T2-T)/DT.
C
      S = H2 - H1
      T0 = 3.*S - S1 - S2
      A0 = 3.*(S-T0)
      B0 = T0 - S2
      D0 = T0*T0 - S1*S2
C
C   H0 has local extrema in (T1,T2) iff S1*S2 < 0 or
C     (T0*(S1+S2) < 0 and D0 .GE. 0).
C
      IF (S1*S2 .GE. 0.  .AND.  (T0*(S1+S2) .GE. 0.  .OR.
     .    D0 .LT. 0.)) GO TO 10
      IF (A0 .EQ. 0.) THEN
C
C   H0 is quadratic and has an extremum at R = -S2/(2*B0).
C     H0(R) = H2 + S2**2/(4*B0).  Note that A0 = 0 implies
C     2*B0 = S1-S2, and S1*S2 < 0 implies B0 .NE. 0.
C     Also, the extremum is a min iff HBND is a lower bound.
C
        F0 = (BND - H2 - S2*S2/(4.*B0))*RF
      ELSE
C
C   A0 .NE. 0 and H0 has extrema at R = (-B0 +/- SQRT(D0))/
C     A0 = S2/(-B0 -/+ SQRT(D0)), where the negative root
C     corresponds to a min.  The expression for R is chosen
C     to avoid cancellation error.  H0(R) = H2 + (S2*B0 +
C     2*D0*R)/(3*A0).
C
        T = -B0 - SIGN(SQRT(D0),B0)
        R = T/A0
        IF (RF*B0 .GT. 0.) R = S2/T
        F0 = (BND - H2 - (S2*B0+2.*D0*R)/(3.*A0))*RF
      ENDIF
C
C   F0 .GE. 0 iff SIG = 0 is sufficient to satisfy the
C     constraint.
C
      IF (F0 .GE. 0.) GO TO 10
C
C Find a zero of F(SIG) = (BND-H(R))*RF where the derivative
C   of H, HP, vanishes at R.  F is a nondecreasing function,
C   F(0) < 0, and F = FMAX for SIG sufficiently large.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
C   where SG0 and SNEG are defined implicitly by DSIG = SIG
C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
C   SNEG is initialized to a sufficiently large value that
C   FNEG > 0.  This value is used only if the initial value
C   of F is negative.
C
      FMAX = MAX(1.E-3,MIN(ABS(H1-BND),ABS(H2-BND)))
      T = MAX(ABS(H1-BND),ABS(H2-BND))
      SIG = MAX(ABS(S1),ABS(S2))/T
      DMAX = SIG*(1.-T/FMAX)
      SNEG = SIG - DMAX
      IF (LUN .GE. 0) WRITE (LUN,120) SIG, SNEG, F0, FMAX
  120 FORMAT (1X,8X,'SIG = ',E15.8,', SNEG = ',E15.8/
     .        1X,9X,'F0 = ',E15.8,', FMAX = ',E15.8/)
      DSIG = SIG
      FNEG = FMAX
      D2 = S2 - S
      D1PD2 = S2 - S1
      NIT = 0
C
C Compute an absolute tolerance FTOL = abs(TOL) and a
C   relative tolerance RTOL = 100*Macheps.
C
      FTOL = ABS(TOL)
      RTOL = 1.
    5 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 5
      RTOL = RTOL*200.
C
C Top of loop:  compute F.
C
    6 EMS = EXP(-SIG)
      IF (SIG .LE. .5) THEN
C
C   Use approximations designed to avoid cancellation error
C     (associated with small SIG) in the modified hyperbolic
C     functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        AA = A/EMS
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
C
C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
C     overflow.
C
        TM = 1. - EMS
        SSINH = TM*(1.+EMS)
        SSM = SSINH - 2.*SIG*EMS
        SCM = TM*TM
        C1 = SIG*SCM*D2 - SSM*D1PD2
        C2 = SIG*SSINH*D2 - SCM*D1PD2
        AA = 2.*(SIG*TM*D2 + (TM-SIG)*D1PD2)
        A = EMS*AA
        E = SIG*SSINH - SCM - SCM
      ENDIF
C
C   HP(R) = (S2 - (C1*SINH(SIG*R) - C2*COSHM(SIG*R))/E)/DT
C     = 0 for ESR = (-B +/- SQRT(D))/A = C/(-B -/+ SQRT(D))
C     where ESR = exp(SIG*R), A = C2-C1, D = B**2 - A*C, and
C     B and C are defined below.
C
      B = E*S2 - C2
      C = C2 + C1
      D = B*B - A*C
      F = 0.
      IF (AA*C .EQ. 0.  .AND.  B .EQ. 0.) GO TO 7
      F = FMAX
      IF (D .LT. 0.) GO TO 7
      T1 = SQRT(D)
      T = -B - SIGN(T1,B)
      RSIG = 0.
      IF (RF*B .LT. 0.  .AND.  AA .NE. 0.) THEN
        IF (T/AA .GT. 0.) RSIG = SIG + LOG(T/AA)
      ENDIF
      IF ((RF*B .GT. 0.  .OR.  AA .EQ. 0.)  .AND.
     .    C/T .GT. 0.) RSIG = LOG(C/T)
      IF ((RSIG .LE. 0.  .OR.  RSIG .GE. SIG)  .AND.
     .    B .NE. 0.) GO TO 7
C
C   H(R) = H2 - (B*SIG*R + C1 + RF*SQRT(D))/(SIG*E).
C
      F = (BND - H2 + (B*RSIG+C1+RF*T1)/(SIG*E))*RF
C
C   Update the number of iterations NIT.
C
    7 NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,130) NIT, SIG, F
  130 FORMAT (1X,3X,I2,' -- SIG = ',E15.8,', F = ',
     .        E15.8)
      IF (F0*F .LT. 0.) THEN
C
C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F and
C     FNEG always have opposite signs.  If SIG is closer to
C     SNEG than SG0, then swap (SNEG,FNEG) with (SG0,F0).
C
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF (ABS(DSIG) .GT. ABS(T1)) THEN
C
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
C
C   Test for convergence.
C
      STOL = RTOL*SIG
      IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 10
C
C   Test for F0 = F = FMAX or F < 0 on the first iteration.
C
      IF (F0 .NE. F  .AND.  (NIT .GT. 1  .OR.  F .GT. 0.))
     .   GO TO 9
C
C   F*F0 > 0 and either the new estimate would be outside
C     of the bracketing interval of length abs(DMAX) or
C     F < 0 on the first iteration.  Reset (SG0,F0) to
C     (SNEG,FNEG).
C
    8 DSIG = DMAX
      F0 = FNEG
C
C   Compute the change in SIG by linear interpolation
C     between (SG0,F0) and (SIG,F).
C
    9 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,140) DSIG
  140 FORMAT (1X,8X,'DSIG = ',E15.8)
      IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .     DSIG*DMAX .GT. 0. ) GO TO 8
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
      IF (ABS(DSIG) .LT. STOL/2.) DSIG = -SIGN(STOL/2.,DMAX)
C
C   Bottom of loop:  Update SIG, DMAX, and F0.
C
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 6
C
C No errors encountered.
C
   10 SIG0 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
C
C Error termination.
C
   11 SIG0 = -1.
      RETURN
      END
      REAL FUNCTION SIG1 (N1,N2,N,X,Y,H,LIST,LPTR,LEND,HXHY,
     .                    IFLGB,HPBND,TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGB,
     .        IFLGS, IER
      REAL    X(N), Y(N), H(N), HXHY(2,N), HPBND, TOL,
     .        SIGMA(*)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/03/98
C
C   Given a triangulation of a set of nodes in the plane,
C along with data values H and gradients (HX,HY) at the
C nodes, this function determines the smallest tension fac-
C tor SIG1 such that the first derivative HP(T) of the
C Hermite interpolatory tension spline H(T), defined by SIG1
C and the endpoint values and directional derivatives asso-
C ciated with an arc N1-N2, is bounded (either above or
C below) by HPBND for all T in (T1,T2), where (T1,T2) de-
C notes an interval corresponding to the arc.
C
C On input:
C
C       N1,N2 = Nodal indexes of the endpoints of an arc for
C               which the tension factor is to be computed.
C               The indexes must be distinct and lie in the
C               range 1 to N, and if IFLGS .GE. 1, they must
C               correspond to adjacent nodes in the triangu-
C               lation.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I)).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C       HXHY = Array dimensioned 2 by N whose columns con-
C              partial derivatives at the nodes (X partials
C              in the first row).  Refer to Subroutines
C              GRADC, GRADG, GRADL, and SMSURF.
C
C       IFLGB = Bound option indicator:
C               IFLGB = -1 if HPBND is a lower bound on HP.
C               IFLGB = 1 if HPBND is an upper bound on HP.
C
C       HPBND = Bound on HP.  HPBND .LE. min(HP1,HP2,S) if
C               IFLGB = -1 and HPBND .GE. max(HP1,HP2,S) if
C               IFLGB = 1, where HP1 and HP2 are the direc-
C               tional derivatives at the endpoints of the
C               arc N1-N2, and S is the slope of the linear
C               interpolant of the endpoint data values.
C
C       TOL = Tolerance whose magnitude determines how close
C             SIG1 is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy the constraint.  For a lower bound,
C             SIG1 is chosen so that HPBND .LE. HPMIN .LE.
C             HPBND + abs(TOL), where HPMIN is the minimum
C             value of HP on the arc.  For an upper bound,
C             the maximum of HP satisfies HPBND - abs(TOL)
C             .LE. HPMAX .LE. HPBND.  Thus, the constraint
C             is satisfied but possibly with more tension
C             than necessary.
C
C       IFLGS = Tension array option indicator:
C               IFLGS .LE. 0 if SIGMA is not to be used.
C               IFLGS .GE. 1 if SIGMA is to be updated by
C                            storing SIG1 in the appropriate
C                            locations.
C
C The above parameters are not altered by this function.
C
C       SIGMA = Dummy array of length 1 (IFLGS .LE. 0) or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               routine GETSIG.
C
C On output:
C
C       SIGMA = Tension factor array updated with the new
C               value if and only if IFLGS .GE. 1 and IER
C               .GE. 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     constraint can be satisfied with fin-
C                     ite tension.
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     the constraint (e.g., IFLGB = -1,
C                     HPBND = S, and HP1 > S).
C             IER = -1 if N1, N2, N, or IFLGB is outside its
C                      valid range.
C             IER = -2 if nodes N1 and N2 coincide or IFLGS
C                      .GE. 1 and the nodes are not adja-
C                      cent.
C             IER = -3 if HPBND is outside its valid range.
C
C       SIG1 = Minimum tension factor defined above unless
C              IER < 0, in which case SIG1 = -1.  If IER
C              = 1, SIG1 is set to 85, resulting in an
C              approximation to the linear interpolant of
C              the endpoint values.
C
C TRIPACK module required by SIG1:  STORE
C
C SRFPACK module required by SIG1:  SNHCSH
C
C Intrinsic functions called by SIG1:  ABS, EXP, MAX, MIN,
C                                        REAL, SIGN, SQRT
C
C***********************************************************
C
      REAL    STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    A, A0, B0, BND, C0, C1, C2, COSHM, COSHMM, D0,
     .        D1, D1PD2, D2, DMAX, DSIG, DT, DX, DY, E, EMS,
     .        EMS2, F, F0, FMAX, FNEG, FTOL, RF, RTOL, S,
     .        S1, S2, SBIG, SIG, SINH, SINHM, STOL, T0, T1,
     .        T2, TM
      DATA SBIG/85./,  LUN/-1/
      RF = REAL(IFLGB)
      BND = HPBND
C
C Print a heading.
C
      IF (LUN .GE. 0) THEN
        IF (RF .LT. 0.) WRITE (LUN,100) N1, N2, BND
        IF (RF .GT. 0.) WRITE (LUN,110) N1, N2, BND
      ENDIF
  100 FORMAT (//1X,'SIG1 -- N1 =',I4,', N2 =',I4,
     .        ', Lower bound = ',E15.8)
  110 FORMAT (//1X,'SIG1 -- N1 =',I4,', N2 =',I4,
     .        ', Upper bound = ',E15.8)
C
C Test for errors and store local parameters.
C
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.
     .    MAX(N1,N2,3) .GT. N  .OR.  ABS(RF) .NE. 1.)
     .   GO TO 10
      IER = -2
      IF (IFLGS .GT. 0) THEN
C
C   Set LP1 and LP2 to the pointers to N2 as a neighbor of
C     N1 and N1 as a neighbor of N2, respectively.
C
        LPL = LEND(N1)
        LP1 = LPTR(LPL)
    1   IF (LIST(LP1) .EQ. N2) GO TO 2
          LP1 = LPTR(LP1)
          IF (LP1 .NE. LPL) GO TO 1
        IF (ABS(LIST(LP1)) .NE. N2) GO TO 10
C
    2   LPL = LEND(N2)
        LP2 = LPTR(LPL)
    3   IF (LIST(LP2) .EQ. N1) GO TO 4
          LP2 = LPTR(LP2)
          IF (LP2 .NE. LPL) GO TO 3
        IF (ABS(LIST(LP2)) .NE. N1) GO TO 10
      ENDIF
C
C Test for arc length DT = SQRT(DX**2+DY**2) = 0.
C
    4 DX = X(N2) - X(N1)
      DY = Y(N2) - Y(N1)
      IF (DX .EQ. 0.  .AND.  DY .EQ. 0.) GO TO 10
C
C Compute first difference S and scaled directional deriva-
C   tives S1,S2 at the endpoints (for the direction N1->N2).
C
      S = H(N2) - H(N1)
      S1 = HXHY(1,N1)*DX + HXHY(2,N1)*DY
      S2 = HXHY(1,N2)*DX + HXHY(2,N2)*DY
C
C Test for a valid constraint.
C
      IER = -3
      IF ((RF .LT. 0.  .AND.  MIN(S1,S2,S) .LT. BND)  .OR.
     .    (RF .GT. 0.  .AND.  BND .LT. MAX(S1,S2,S)))
     .   GO TO 10
C
C Test for infinite tension required.
C
      IER = 1
      SIG = SBIG
      IF (S .EQ. BND  .AND.  (S1 .NE. S  .OR.  S2 .NE. S))
     .   GO TO 9
C
C Test for SIG = 0 sufficient.  The Hermite cubic interpo-
C   lant H0 has derivative HP0(T) = (S2 + 2*B0*R + A0*R**2)/
C   DT, where R = (T2-T)/DT.
C
      IER = 0
      SIG = 0.
      T0 = 3.*S - S1 - S2
      B0 = T0 - S2
      C0 = T0 - S1
      A0 = -B0 - C0
C
C   HP0(R) has an extremum (at R = -B0/A0) in (0,1) iff
C     B0*C0 > 0 and the third derivative of H0 has the
C     sign of A0.
C
      IF (B0*C0 .LE. 0.  .OR.  A0*RF .GT. 0.) GO TO 9
C
C   A0*RF < 0 and HP0(R) = -D0/(DT*A0) at R = -B0/A0.
C
      DT = SQRT(DX*DX + DY*DY)
      D0 = T0*T0 - S1*S2
      F0 = (BND + D0/(A0*DT))*RF
      IF (F0 .GE. 0.) GO TO 9
C
C Find a zero of F(SIG) = (BND-HP(R))*RF, where HP has an
C   extremum at R.  F has a unique zero, F(0) = F0 < 0, and
C   F = (BND-S)*RF > 0 for SIG sufficiently large.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
C   where SG0 and SNEG are defined implicitly by DSIG = SIG
C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
C   SIG is initialized to the zero of (BND - (SIG*S-S1-S2)/
C   (DT*(SIG-2.)))*RF -- a value for which F(SIG) .GE. 0 and
C   F(SIG) = 0 for SIG sufficiently large that 2*SIG is in-
C   significant relative to exp(SIG).
C
      FMAX = (BND-S/DT)*RF
      SIG = 2. - A0/(3.*(DT*BND-S))
      IF (LUN .GE. 0) WRITE (LUN,120) F0, FMAX, SIG
  120 FORMAT (1X,9X,'F0 = ',E15.8,', FMAX = ',E15.8/
     .        1X,8X,'SIG = ',E15.8/)
      IF (STORE(SIG*EXP(-SIG)+.5) .EQ. .5) GO TO 9
      DSIG = SIG
      DMAX = -2.*SIG
      FNEG = FMAX
      D1 = S - S1
      D2 = S2 - S
      D1PD2 = D1 + D2
      NIT = 0
C
C Compute an absolute tolerance FTOL = abs(TOL), and a
C   relative tolerance RTOL = 100*Macheps.
C
      FTOL = ABS(TOL)
      RTOL = 1.
    5 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 5
      RTOL = RTOL*200.
C
C Top of loop:  compute F.
C
    6 IF (SIG .LE. .5) THEN
C
C   Use approximations designed to avoid cancellation
C     error (associated with small SIG) in the modified
C     hyperbolic functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
C
C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
C     overflow.
C
        EMS = EXP(-SIG)
        EMS2 = EMS + EMS
        TM = 1. - EMS
        SINH = TM*(1.+EMS)
        SINHM = SINH - SIG*EMS2
        COSHM = TM*TM
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*SINH*D2 - COSHM*D1PD2
        A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
        E = SIG*SINH - COSHM - COSHM
      ENDIF
C
C   The second derivative HPP of H(R) has a zero at exp(SIG*
C     R) = SQRT((C2+C1)/A) and R is in (0,1) and well-
C     defined iff HPP(T1)*HPP(T2) < 0.
C
      F = FMAX
      T1 = A*(C2+C1)
      IF (T1 .GE. 0.) THEN
        IF (C1*(SIG*COSHM*D1 - SINHM*D1PD2) .LT. 0.) THEN
C
C   HP(R) = (B+SIGN(A)*SQRT(A*C))/(DT*E) at the critical
C     value of R, where A = C2-C1, B = E*S2-C2, and C = C2 +
C     C1.  Note that RF*A < 0.
C
          F = (BND - (E*S2-C2 - RF*SQRT(T1))/(DT*E))*RF
        ENDIF
      ENDIF
C
C   Update the number of iterations NIT.
C
      NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,130) NIT, SIG, F
  130 FORMAT (1X,3X,I2,' -- SIG = ',E15.8,', F = ',
     .        E15.8)
      IF (F0*F .LT. 0.) THEN
C
C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F
C     and FNEG always have opposite signs.  If SIG is closer
C     to SNEG than SG0 and abs(F) < abs(FNEG), then swap
C     (SNEG,FNEG) with (SG0,F0).
C
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF ( ABS(DSIG) .GT. ABS(T1)  .AND.
     .       ABS(F) .LT. ABS(T2) ) THEN
C
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
C
C   Test for convergence.
C
      STOL = RTOL*SIG
      IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 9
      IF (F0*F .LT. 0.  .OR.  ABS(F) .LT. ABS(F0)) GO TO 8
C
C   F*F0 > 0 and the new estimate would be outside of the
C     bracketing interval of length abs(DMAX).  Reset
C     (SG0,F0) to (SNEG,FNEG).
C
    7 DSIG = DMAX
      F0 = FNEG
C
C   Compute the change in SIG by linear interpolation
C     between (SG0,F0) and (SIG,F).
C
    8 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,140) DSIG
  140 FORMAT (1X,8X,'DSIG = ',E15.8)
      IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .     DSIG*DMAX .GT. 0. ) GO TO 7
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
      IF (ABS(DSIG) .LT. STOL/2.) DSIG = -SIGN(STOL/2.,DMAX)
C
C   Bottom of loop:  update SIG, DMAX, and F0.
C
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 6
C
C No errors encountered.
C
    9 SIG1 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
C
C Error termination.
C
   10 SIG1 = -1.
      RETURN
      END
      REAL FUNCTION SIG2 (N1,N2,N,X,Y,H,LIST,LPTR,LEND,HXHY,
     .                    TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGS,
     .        IER
      REAL    X(N), Y(N), H(N), HXHY(2,N), TOL, SIGMA(*)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/22/97
C
C   Given a triangulation of a set of nodes in the plane,
C along with data values H and gradients (HX,HY) at the
C nodes, this function determines the smallest tension fac-
C tor SIG2 such that the Hermite interpolatory tension
C spline H(T), defined by SIG2 and the endpoint values and
C directional derivatives associated with an arc N1-N2,
C preserves convexity (or concavity) of the data:
C
C   HP1 .LE. S .LE. HP2 implies HPP(T) .GE. 0, and
C   HP1 .GE. S .GE. HP2 implies HPP(T) .LE. 0
C
C for all T in the open interval (T1,T2) corresponding to
C the arc, where HP1 and HP2 are the derivative values of H
C at the endpoints, S is the slope of the linear interpolant
C of the endpoint data values, and HPP denotes the second
C derivative of H.  Note, however, that infinite tension is
C required if HP1 = S or HP2 = S (unless HP1 = HP2 = S).
C
C On input:
C
C       N1,N2 = Nodal indexes of the endpoints of an arc for
C               which the tension factor is to be computed.
C               The indexes must be distinct and lie in the
C               range 1 to N, and if IFLGS .GE. 1, they must
C               correspond to adjacent nodes in the triangu-
C               lation.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I)).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C       HXHY = Array dimensioned 2 by N whose columns con-
C              partial derivatives at the nodes (X partials
C              in the first row).  Refer to Subroutines
C              GRADC, GRADG, GRADL, and SMSURF.
C
C       TOL = Tolerance whose magnitude determines how close
C             SIG2 is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy convexity or concavity.  In the case
C             convexity, SIG2 is chosen so that 0 .LE.
C             HPPMIN .LE. abs(TOL), where HPPMIN is the
C             minimum value of HPP on the arc.  In the case
C             of concavity, the maximum value of HPP satis-
C             fies -abs(TOL) .LE. HPPMAX .LE. 0.  Thus, the
C             constraint is satisfied but possibly with more
C             tension than necessary.
C
C       IFLGS = Tension array option indicator:
C               IFLGS .LE. 0 if SIGMA is not to be used.
C               IFLGS .GE. 1 if SIGMA is to be updated by
C                            storing SIG2 in the appropriate
C                            locations.
C
C The above parameters are not altered by this function.
C
C       SIGMA = Dummy array of length 1 (IFLGS .LE. 0) or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               routine GETSIG.
C
C On output:
C
C       SIGMA = Tension factor array updated with the new
C               value if and only if IFLGS .GE. 1 and IER
C               .GE. 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and fin-
C                     ite tension is sufficient to satisfy
C                     convexity (or concavity).
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     convexity.
C             IER = 2 if the data does not satisfy convexity
C                     or concavity.
C             IER = -1 if N1, N2, or N is outside its valid
C                      range.
C             IER = -2 if nodes N1 and N2 coincide or IFLGS
C                      .GE. 1 and the nodes are not adja-
C                      cent.
C
C       SIG2 = Minimum tension factor defined above unless
C              IER < 0, in which case SIG2 = -1.  If IER
C              = 1, SIG2 is set to 85, resulting in an
C              approximation to the linear interpolant of
C              the endpoint values.  If IER = 2, SIG2 = 0,
C              resulting in the Hermite cubic interpolant.
C
C TRIPACK module required by SIG2:  STORE
C
C SRFPACK module required by SIG2:  SNHCSH
C
C Intrinsic functions called by SIG2:  ABS, EXP, MAX, MIN,
C                                        SQRT
C
C***********************************************************
C
      REAL    STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    COSHM, D1, D1D2, D2, DSIG, DUMMY, DX, DY, EMS,
     .        F, FP, FTOL, RTOL, S, SBIG, SIG, SINHM, SSM,
     .        T, T1, TP1
      DATA SBIG/85./,  LUN/-1/
C
C Print a heading.
C
      IF (LUN .GE. 0) WRITE (LUN,100) N1, N2
  100 FORMAT (//1X,'SIG2 -- N1 =',I4,', N2 =',I4)
C
C Test for errors and store local parameters.
C
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.
     .    MAX(N1,N2,3) .GT. N) GO TO 8
      IER = -2
      IF (IFLGS .GT. 0) THEN
C
C   Set LP1 and LP2 to the pointers to N2 as a neighbor of
C     N1 and N1 as a neighbor of N2, respectively.
C
        LPL = LEND(N1)
        LP1 = LPTR(LPL)
    1   IF (LIST(LP1) .EQ. N2) GO TO 2
          LP1 = LPTR(LP1)
          IF (LP1 .NE. LPL) GO TO 1
        IF (ABS(LIST(LP1)) .NE. N2) GO TO 8
C
    2   LPL = LEND(N2)
        LP2 = LPTR(LPL)
    3   IF (LIST(LP2) .EQ. N1) GO TO 4
          LP2 = LPTR(LP2)
          IF (LP2 .NE. LPL) GO TO 3
        IF (ABS(LIST(LP2)) .NE. N1) GO TO 8
      ENDIF
C
C Test for arc length DT = SQRT(DX**2+DY**2) = 0.
C
    4 DX = X(N2) - X(N1)
      DY = Y(N2) - Y(N1)
      IF (DX .EQ. 0.  .AND.  DY .EQ. 0.) GO TO 8
C
C Compute first and second differences and test for infinite
C   tension required.
C
      S = H(N2) - H(N1)
      D1 = S - HXHY(1,N1)*DX - HXHY(2,N1)*DY
      D2 = HXHY(1,N2)*DX + HXHY(2,N2)*DY - S
      D1D2 = D1*D2
      IER = 1
      SIG = SBIG
      IF (D1D2 .EQ. 0.  .AND.  D1 .NE. D2) GO TO 7
C
C Test for a valid constraint.
C
      IER = 2
      SIG = 0.
      IF (D1D2 .LT. 0.) GO TO 7
C
C Test for SIG = 0 sufficient.
C
      IER = 0
      IF (D1D2 .EQ. 0.) GO TO 7
      T = MAX(D1/D2,D2/D1)
      IF (T .LE. 2.) GO TO 7
C
C Find a zero of F(SIG) = SIG*COSHM(SIG)/SINHM(SIG) - (T+1).
C   Since the derivative of F vanishes at the origin, a
C   quadratic approximation is used to obtain an initial
C   estimate for the Newton method.
C
      TP1 = T + 1.
      SIG = SQRT(10.*T-20.)
      NIT = 0
C
C   Compute an absolute tolerance FTOL = abs(TOL) and a
C     relative tolerance RTOL = 100*Macheps.
C
      FTOL = ABS(TOL)
      RTOL = 1.
    5 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 5
      RTOL = RTOL*200.
C
C Top of loop:  evaluate F and its derivative FP.
C
    6 IF (SIG .LE. .5) THEN
C
C   Use approximations designed to avoid cancellation error
C     in the hyperbolic functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,DUMMY)
        T1 = COSHM/SINHM
        FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.)
      ELSE
C
C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
C     overflow.
C
        EMS = EXP(-SIG)
        SSM = 1. - EMS*(EMS+SIG+SIG)
        T1 = (1.-EMS)*(1.-EMS)/SSM
        FP = T1 + SIG*(2.*SIG*EMS/SSM - T1*T1 + 1.)
      ENDIF
C
      F = SIG*T1 - TP1
C
C   Update the number of iterations NIT.
C
      NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,110) NIT, SIG, F, FP
  110 FORMAT (1X,3X,I2,' -- SIG = ',E15.8,', F = ',
     .        E15.8/1X,31X,'FP = ',E15.8)
C
C   Test for convergence.
C
      IF (FP .LE. 0.) GO TO 7
      DSIG = -F/FP
      IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 7
C
C   Bottom of loop:  update SIG.
C
      SIG = SIG + DSIG
      GO TO 6
C
C No errors encountered.
C
    7 SIG2 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
C
C Error termination.
C
    8 SIG2 = -1.
      RETURN
      END
      SUBROUTINE SMSGS (NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,
     .                  IFLGS,SIGMA,W,P, NIT,DFMAX,F,
     .                  FXFY, IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        IFLGS, NIT, IER
      REAL    X(N), Y(N), Z(N), SIGMA(*), W(N), P, DFMAX,
     .        F(N), FXFY(2,N)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   08/26/91
C
C   This subroutine employs the block Gauss-Seidel method
C (3 by 3 blocks) to solve the order 3N symmetric positive
C definite linear system associated with minimizing the
C quadratic functional Q(F,FX,FY) described in Subroutine
C SMSURF.
C
C   Note that small relative changes in F can cause large
C relative changes in FX and FY, resulting in an ill-
C conditioned system.  However, good F values should be
C achieved with a small number of iterations, and the
C gradients (with fixed F) can then be improved by a call
C to Subroutine GRADG.
C
C On input:
C
C   NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,IFLGS,SIGMA,W =
C Parameters defined in Subroutine SMSURF.
C
C       P = Positive smoothing parameter defining Q.
C
C The above parameters are not altered by this routine.
C
C       NIT = Maximum number of Gauss-Seidel iterations to
C             be employed.  This maximum will likely be
C             achieved if DFMAX is smaller than the machine
C             precision.  NIT .GE. 0.
C
C       DFMAX = Nonnegative convergence criterion.  The
C               method is terminated when the maximum
C               change in a solution F-component between
C               iterations is at most DFMAX.  The change in
C               a component is taken to be the absolute
C               difference relative to 1 plus the old value.
C
C       F = Initial estimate of the first N solution compo-
C           nents.
C
C       FXFY = 2 by N array containing initial estimates of
C              the last 2N solution components.
C
C On output:
C
C       NIT = Number of Gauss-Seidel iterations employed.
C
C       DFMAX = Maximum relative change in a solution F-
C               component at the last iteration.
C
C       F = First N solution components -- function values
C           at the nodes.
C
C       FXFY = Last 2N solution components -- gradients at
C              the nodes with X partials in the first row.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     convergence criterion was achieved.
C             IER = 1 if no errors were encountered but con-
C                     vergence was not achieved within NIT
C                     iterations.
C             IER = -1 if NCC, N, P, NIT, or DFMAX is out-
C                      side its valid range on input.  F
C                      and FXFY are not altered in this
C                      case.  LCC is not tested for valid-
C                      ity.
C             IER = -2 if all nodes are collinear or the
C                      triangulation data structure is not
C                      valid.
C             IER = -3 if duplicate nodes were encountered.
C
C SRFPACK modules required by SMSGS:  GRCOEF, SNHCSH
C
C Intrinsic functions called by SMSGS:  ABS, MAX, SQRT
C
C***********************************************************
C
      INTEGER I, IFL, IFRST, ILAST, ITER, ITMAX, J, K, KBAK,
     .        KFOR, LCC1, LP, LPJ, LPL, LPLJ, NB, NN
      REAL    C11, C12, C13, C22, C23, C33, CC22, CC23,
     .        CC33, DCUB, DET, DF, DFMX, DFX, DFY, DSQ, DX,
     .        DXS, DXDY, DY, DYS, FK, FXJ, FXK, FYJ, FYK,
     .        PP, R1, R2, R3, RR2, RR3, SIG, T, T1, T2, T3,
     .        TOL, TRMX, TRMY, XK, YK
C
      NN = N
      IFL = IFLGS
      PP = P
      ITMAX = NIT
      TOL = DFMAX
      IF (NCC .EQ. 0) THEN
        LCC1 = NN+1
      ELSE
        LCC1 = LCC(1)
      ENDIF
C
C Test for errors in input and initialize the iteration
C   count ITER, tension factor SIG, and output value of
C   DFMAX.
C
      IF (NCC .LT. 0  .OR.  NN .LT. 3  .OR.  PP .LE. 0.
     .    .OR.  ITMAX .LT. 0  .OR.  TOL .LT. 0.) GO TO 8
      ITER = 0
      SIG = SIGMA(1)
      DFMX = 0.
C
C Top of iteration loop:  If K is a constraint node, I
C   indexes the constraint containing node K, IFRST and
C   ILAST are the first and last nodes of constraint I,
C   and (KBAK,K,KFOR) is a subsequence of constraint I.
C
    1 IF (ITER .EQ. ITMAX) GO TO 7
      DFMX = 0.
      I = 0
      ILAST = LCC1-1
      KBAK = 0
      KFOR = 0
C
C   Loop on nodes.
C
      DO 6 K = 1,NN
        IF (K .GE. LCC1) THEN
          IF (K .GT. ILAST) THEN
            I = I + 1
            IFRST = K
            IF (I .LT. NCC) THEN
              ILAST = LCC(I+1) - 1
            ELSE
              ILAST = NN
            ENDIF
            KBAK = ILAST
            KFOR = K + 1
          ELSE
            KBAK = K - 1
            IF (K .LT. ILAST) THEN
              KFOR = K + 1
            ELSE
              KFOR = IFRST
            ENDIF
          ENDIF
        ENDIF
        XK = X(K)
        YK = Y(K)
        FK = F(K)
        FXK = FXFY(1,K)
        FYK = FXFY(2,K)
C
C   Initialize components of the order 3 system for the
C     change (DF,DFX,DFY) in the K-th solution components
C     (symmetric matrix in C and residual in R).
C
        C11 = PP*W(K)
        C12 = 0.
        C13 = 0.
        C22 = 0.
        C23 = 0.
        C33 = 0.
        R1 = C11*(Z(K)-FK)
        R2 = 0.
        R3 = 0.
C
C   Loop on neighbors J of node K.
C
        LPL = LEND(K)
        LPJ = LPL
    2   LPJ = LPTR(LPJ)
          J = ABS(LIST(LPJ))
C
C   Arc K-J lies in a constraint region and is bypassed iff
C     K and J are nodes in the same constraint and J follows
C     KFOR and precedes KBAK as a neighbor of K.  Also, K-J
C     is bypassed if it is both a constraint arc and a
C     boundary arc of the triangulation.
C
          IF (K .LT. LCC1  .OR.  J .LT. IFRST  .OR.
     .        J .GT. ILAST) GO TO 4
          IF (J .EQ. KBAK  .OR.  J .EQ. KFOR) THEN
            LPLJ = LEND(J)
            IF (LIST(LPL) .EQ. -J  .OR.
     .          LIST(LPLJ) .EQ. -K) THEN
              GO TO 5
            ELSE
              GO TO 4
            ENDIF
          ENDIF
          LP = LPJ
C
    3     LP = LPTR(LP)
            NB = ABS(LIST(LP))
            IF (NB .EQ. KBAK) GO TO 5
            IF (NB .NE. KFOR) GO TO 3
C
C   Compute parameters associated with edge K->J, and test
C     for duplicate nodes.
C
    4     DX = X(J) - XK
          DY = Y(J) - YK
          DXS = DX*DX
          DXDY = DX*DY
          DYS = DY*DY
          DSQ = DXS + DYS
          DCUB = DSQ*SQRT(DSQ)
          IF (DCUB .EQ. 0.) GO TO 10
          IF (IFL .GE. 1) SIG = SIGMA(LPJ)
          CALL GRCOEF (SIG,DCUB, T3,T2)
          T1 = T2 + T3
C
C   T1 = SIG*SIG*COSHM/(DCUB*E), T2 = SIG*SINHM/(DCUB*E),
C     and T3 = SIG*(SIG*COSHM-SINHM)/(DCUB*E) for E =
C     SIG*SINH - 2*COSHM.
C
          T = T1*(FK-F(J))
          FXJ = FXFY(1,J)
          FYJ = FXFY(2,J)
C
C   Update the system components for node J.
C
          C11 = C11 + T1 + T1
          C12 = C12 + T1*DX
          C13 = C13 + T1*DY
          C22 = C22 + T3*DXS
          C23 = C23 + T3*DXDY
          C33 = C33 + T3*DYS
          R1 = R1 - T - T - T1*(DX*(FXK+FXJ) + DY*(FYK+FYJ))
          TRMX = T3*FXK + T2*FXJ
          TRMY = T3*FYK + T2*FYJ
          R2 = R2 - T*DX - TRMX*DXS - TRMY*DXDY
          R3 = R3 - T*DY - TRMX*DXDY - TRMY*DYS
C
C   Bottom of loop on neighbors.
C
    5     IF (LPJ .NE. LPL) GO TO 2
C
C   Solve the system associated with the K-th block.
C
        CC22 = C11*C22 - C12*C12
        CC23 = C11*C23 - C12*C13
        CC33 = C11*C33 - C13*C13
        RR2 = C11*R2 - C12*R1
        RR3 = C11*R3 - C13*R1
        DET = CC22*CC33 - CC23*CC23
        IF (DET .EQ. 0.  .OR.  CC22 .EQ. 0.  .OR.
     .      C11 .EQ. 0.) GO TO 9
        DFY = (CC22*RR3 - CC23*RR2)/DET
        DFX = (RR2 - CC23*DFY)/CC22
        DF = (R1 - C12*DFX - C13*DFY)/C11
C
C   Update the solution components for node K and the
C     maximum relative change in F.
C
        F(K) = FK + DF
        FXFY(1,K) = FXK + DFX
        FXFY(2,K) = FYK + DFY
        DFMX = MAX(DFMX,ABS(DF)/(1.+ABS(FK)))
    6   CONTINUE
C
C   Increment ITER and test for convergence.
C
      ITER = ITER + 1
      IF (DFMX .GT. TOL) GO TO 1
C
C Method converged.
C
      NIT = ITER
      DFMAX = DFMX
      IER = 0
      RETURN
C
C Method failed to converge within NIT iterations.
C
    7 DFMAX = DFMX
      IER = 1
      RETURN
C
C Invalid input parameter.
C
    8 NIT = 0
      DFMAX = 0.
      IER = -1
      RETURN
C
C Node K and its neighbors are collinear, resulting in a
C   singular system.
C
    9 NIT = 0
      DFMAX = DFMX
      IER = -2
      RETURN
C
C Nodes J and K coincide.
C
   10 NIT = 0
      DFMAX = DFMX
      IER = -3
      RETURN
      END
      SUBROUTINE SMSURF (NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,
     .                   IFLGS,SIGMA,W,SM,SMTOL,GSTOL, F,
     .                   FXFY,IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        IFLGS, IER
      REAL    X(N), Y(N), Z(N), SIGMA(*), W(N), SM, SMTOL,
     .        GSTOL, F(N), FXFY(2,N)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/03/98
C
C   Given a triangulation of N nodes in the plane, along
C with data values Z at the nodes and tension factors SIGMA
C associated with the arcs, this subroutine determines a
C set of nodal function values F and gradients (FX,FY) such
C that a quadratic functional Q1(F,FX,FY) is minimized sub-
C ject to the constraint Q2(F) .LE. SM for Q2(F) = (Z-F)**T*
C W*(Z-F), where W is a diagonal matrix of positive weights.
C The functional Q1 is an approximation to the linearized
C curvature over the triangulation of a C-1 bivariate func-
C tion F(X,Y) which interpolates the nodal values and
C gradients.  Subroutines INTRC1 and UNIF may be called to
C evaluate F at arbitrary points.
C
C   The smoothing procedure is an extension of the method
C for cubic spline smoothing due to C. Reinsch -- Numer.
C Math., 10 (1967) and 16 (1971).  Refer to Subroutines FVAL
C and TVAL for a further description of the interpolant F.
C Letting D1F(T) and D2F(T) denote first and second deriva-
C tives of F with respect to a parameter T varying along a
C triangulation arc, Q1 is the sum over the triangulation
C arcs, excluding interior constraint arcs, of the integrals
C of
C
C      D2F(T)**2 + [(SIGMA/L)*(D1F(T)-S)]**2 ,
C
C where L denotes arc length, SIGMA is the appropriate ten-
C sion factor, and S is the slope of the linear function of
C T which interpolates the values of F at the endpoints of
C the arc.  Introducing a smoothing parameter P, and assum-
C ing the constraint is active, the problem is equivalent to
C minimizing
C
C      Q(P,F,FX,FY) = Q1(F,FX,FY) + P*(Q2(F)-SM) .
C
C The secant method is used to find a zero of
C
C      G(P) = 1/SQRT(Q2) - 1/SQRT(SM) ,
C
C where F(P) satisfies the order 3N symmetric positive def-
C inite linear system obtained by setting the gradient of Q
C (treated as a function of F, FX, and FY) to zero.  The
C linear system is solved by the Gauss-Seidel method.
C
C   Note that the method can also be used to select grad-
C ients for the interpolation problem (F = Z, SM = 0, and P
C infinite).  This is achieved by a call to Subroutine
C GRADG.
C
C On input:
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index of the
C             first node of constraint I in LCC(I).  For I =
C             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       Z = Array of length N containing data values at the
C           nodes.  Z(I) is associated with (X(I),Y(I)).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C       IFLGS = Tension factor option:
C               IFLGS .LE. 0 if a single uniform tension
C                            factor is to be used.
C               IFLGS .GE. 1 if variable tension is desired.
C
C       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               routines GETSIG, SIG0, SIG1, and SIG2.
C
C       W = Array of length N containing positive weights
C           associated with the data values.  The recommend-
C           ed value of W(I) is 1/DZ**2, where DZ is the
C           standard deviation associated with Z(I).  DZ**2
C           is the expected value of the squared error in
C           the measurement of Z(I).  (The mean error is
C           assumed to be zero.)
C
C       SM = Positive parameter specifying an upper bound on
C            Q2(F).  Note that F(X,Y) is linear (and Q2(F)
C            is minimized) if SM is sufficiently large that
C            the constraint is not active.  It is recommend-
C            ed that SM satisfy N-SQRT(2N) .LE. SM .LE. N+
C            SQRT(2N).
C
C       SMTOL = Parameter in the open interval (0,1) speci-
C               fying the relative error allowed in satisfy-
C               ing the constraint -- the constraint is
C               assumed to be satisfied if SM*(1-SMTOL) .LE.
C               Q2 .LE. SM*(1+SMTOL).  A reasonable value
C               for SMTOL is SQRT(2/N).
C
C       GSTOL = Nonnegative tolerance defining the conver-
C               gence criterion for the Gauss-Seidel method.
C               Refer to parameter DFMAX in Subroutine
C               SMSGS.  A recommended value is .05*DU**2,
C               where DU is an average standard deviation
C               in the data values.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       F = Array of length N containing nodal function val-
C           ues unless IER < 0.
C
C       FXFY = 2 by N array whose columns contain partial
C              derivatives of F at the nodes unless IER < 0,
C              with FX in the first row, FY in the second.
C
C       IER = Error indicator and information flag:
C             IER = 0 if no errors were encountered and the
C                     constraint is active -- Q2(F) is ap-
C                     proximately equal to SM.
C             IER = 1 if no errors were encountered but the
C                     constraint is not active -- F, FX, and
C                     FY are the values and partials of a
C                     linear function which minimizes Q2(F),
C                     and Q1 = 0.
C             IER = -1 if NCC, an LCC entry, N, W, SM,
C                      SMTOL, or GSTOL is outside its
C                      valid range on input.
C             IER = -2 if all nodes are collinear or the
C                      triangulation data structure is not
C                      valid.
C             IER = -3 if duplicate nodes were encountered.
C
C SRFPACK modules required by SMSURF:  GRCOEF, SMSGS, SNHCSH
C
C Intrinsic functions called by SMSURF:  ABS, SQRT
C
C***********************************************************
C
      INTEGER I, IERR, ITER, LCCIP1, LUN, NIT, NITMAX, NN
      REAL    C11, C12, C13, C22, C23, C33, CC22, CC23,
     .        CC33, DET, DFMAX, DMAX, DP, F0, FX, FY, G, G0,
     .        GNEG, P, Q2, Q2MAX, Q2MIN, R1, R2, R3, RR2,
     .        RR3, S, TOL, WI, WIXI, WIYI, WIZI, XI, YI
C
      DATA NITMAX/40/,  LUN/-1/
C
C LUN = Logical unit on which diagnostic messages are print-
C       ed (unless LUN < 0).  For each secant iteration, the
C       following values are printed:  P, G(P), NIT, DFMAX,
C       and DP, where NIT denotes the number of Gauss-Seidel
C       iterations used in the computation of G, DFMAX de-
C       notes the maximum relative change in a solution
C       component in the last Gauss-Seidel iteration, and
C       DP is the change in P computed by linear interpola-
C       tion between the current point (P,G) and a previous
C       point.
C
      NN = N
      TOL = GSTOL
C
C Test for errors in input parameters.
C
      IER = -1
      IF (NCC .LT. 0  .OR.  SM .LE. 0.  .OR.  SMTOL .LE. 0.
     .    .OR.  SMTOL .GE. 1.  .OR.  TOL .LE. 0.) RETURN
      IF (NCC .EQ. 0) THEN
        IF (NN .LT. 3) RETURN
      ELSE
        LCCIP1 = NN+1
        DO 1 I = NCC,1,-1
          IF (LCCIP1-LCC(I) .LT. 3) RETURN
          LCCIP1 = LCC(I)
    1     CONTINUE
        IF (LCCIP1 .LT. 1) RETURN
      ENDIF
C
C Compute the components of the 3 by 3 system (normal
C   equations) for the weighted least squares linear fit.
C
      C11 = 0.
      C12 = 0.
      C13 = 0.
      C22 = 0.
      C23 = 0.
      C33 = 0.
      R1 = 0.
      R2 = 0.
      R3 = 0.
      DO 2 I = 1,NN
        WI = W(I)
        IF (WI .LE. 0.) RETURN
        XI = X(I)
        YI = Y(I)
        WIXI = WI*XI
        WIYI = WI*YI
        WIZI = WI*Z(I)
        C11 = C11 + WIXI*XI
        C12 = C12 + WIXI*YI
        C13 = C13 + WIXI
        C22 = C22 + WIYI*YI
        C23 = C23 + WIYI
        C33 = C33 + WI
        R1 = R1 + WIZI*XI
        R2 = R2 + WIZI*YI
        R3 = R3 + WIZI
    2   CONTINUE
C
C Solve the system for (FX,FY,F0) where (FX,FY) is the
C   gradient (constant) and F0 = F(0,0).
C
      CC22 = C11*C22 - C12*C12
      CC23 = C11*C23 - C12*C13
      CC33 = C11*C33 - C13*C13
      RR2 = C11*R2 - C12*R1
      RR3 = C11*R3 - C13*R1
      DET = CC22*CC33 - CC23*CC23
      IER = -2
      IF (DET .EQ. 0.  .OR.  CC22 .EQ. 0.  .OR.
     .    C11 .EQ. 0.) RETURN
      F0 = (CC22*RR3 - CC23*RR2)/DET
      FY = (RR2 - CC23*F0)/CC22
      FX = (R1 - C12*FY - C13*F0)/C11
C
C Compute nodal values and gradients, and accumulate Q2 =
C   (Z-F)**T*W*(Z-F).
C
      Q2 = 0.
      DO 3 I = 1,NN
        F(I) = FX*X(I) + FY*Y(I) + F0
        FXFY(1,I) = FX
        FXFY(2,I) = FY
        Q2 = Q2 + W(I)*(Z(I)-F(I))**2
    3   CONTINUE
C
C Compute bounds on Q2 defined by SMTOL, and test for the
C   constraint satisfied by the linear fit.
C
      Q2MIN = SM*(1.-SMTOL)
      Q2MAX = SM*(1.+SMTOL)
      IF (Q2 .LE. Q2MAX) THEN
C
C   The constraint is satisfied by a planar surface.
C
        IER = 1
        IF (LUN .GE. 0) WRITE (LUN,100)
  100   FORMAT (///1X,'SMSURF:  The constraint is not ',
     .          'active and the surface is linear.'/)
        RETURN
      ENDIF
C
C Compute G0 = G(0) and print a heading.
C
      IER = 0
      S = 1./SQRT(SM)
      G0 = 1./SQRT(Q2) - S
      IF (LUN .GE. 0) WRITE (LUN,110) SM, TOL, NITMAX, G0
  110 FORMAT (///1X,'SMSURF -- SM = ',E10.4,', GSTOL = ',
     .        E7.1,', NITMAX = ',I2,', G(0) = ',E15.8)
C
C G(P) is strictly increasing and concave, and G(0) < 0.
C   Initialize parameters for the secant method.  The method
C   uses three points:  (P0,G0), (P,G), and (PNEG,GNEG),
C   where P0 and PNEG are defined implicitly by DP = P - P0
C   and DMAX = P - PNEG.
C
      P = 10.*SM
      DP = P
      DMAX = 0.
      ITER = 0
C
C Top of loop -- compute G.
C
    4 NIT = NITMAX
      DFMAX = TOL
      CALL SMSGS (NCC,LCC,NN,X,Y,Z,LIST,LPTR,LEND,IFLGS,
     .            SIGMA,W,P, NIT,DFMAX,F,FXFY, IERR)
      IF (IERR .LT. 0) THEN
        IER = IERR
        RETURN
      ENDIF
      Q2 = 0.
      DO 5 I = 1,NN
        Q2 = Q2 + W(I)*(Z(I)-F(I))**2
    5   CONTINUE
      G = 1./SQRT(Q2) - S
      ITER = ITER + 1
      IF (LUN .GE. 0) WRITE (LUN,120) ITER, P, G, NIT, DFMAX
  120 FORMAT (/1X,I2,' -- P = ',E15.8,', G = ',E15.8,
     .        ', NIT = ',I2,', DFMAX = ',E12.6)
C
C   Test for convergence.
C
      IF (G .EQ. G0  .OR.  (Q2MIN .LE. Q2  .AND.
     .                      Q2 .LE. Q2MAX)) RETURN
      IF (DMAX .NE. 0.  .OR.  G .GT. 0.) GO TO 6
C
C   Increase P until G(P) > 0.
C
      P = 10.*P
      DP = P
      GO TO 4
C
C   A bracketing interval [P0,P] has been found.
C
    6 IF (G0*G .LE. 0.) THEN
C
C   G0*G < 0.  Update (PNEG,GNEG) to (P0,G0) so that G
C     and GNEG always have opposite signs.
C
        DMAX = DP
        GNEG = G0
      ENDIF
C
C   Compute the change in P by linear interpolation between
C     (P0,G0) and (P,G).
C
    7 DP = -G*DP/(G-G0)
      IF (LUN .GE. 0) WRITE (LUN,130) DP
  130 FORMAT (1X,5X,'DP = ',E15.8)
      IF (ABS(DP) .GT. ABS(DMAX)) THEN
C
C   G0*G > 0 and the new estimate would be outside of the
C     bracketing interval of length abs(DMAX).  Reset
C     (P0,G0) to (PNEG,GNEG).
C
        DP = DMAX
        G0 = GNEG
        GO TO 7
      ENDIF
C
C   Bottom of loop -- update P, DMAX, and G0.
C
      P = P + DP
      DMAX = DMAX + DP
      G0 = G
      GO TO 4
      END
      SUBROUTINE SNHCSH (X, SINHM,COSHM,COSHMM)
      REAL X, SINHM, COSHM, COSHMM
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   03/18/90
C
C   This subroutine computes approximations to the modified
C hyperbolic functions defined below with relative error
C bounded by 4.7E-12 for a floating point number system with
C sufficient precision.  For IEEE standard single precision,
C the relative error is less than 1.E-5 for all x.
C
C   Note that the 13-digit constants in the data statements
C below may not be acceptable to all compilers.
C
C On input:
C
C       X = Point at which the functions are to be
C           evaluated.
C
C X is not altered by this routine.
C
C On output:
C
C       SINHM = sinh(X) - X.
C
C       COSHM = cosh(X) - 1.
C
C       COSHMM = cosh(X) - 1 - X*X/2.
C
C Modules required by SNHCSH:  None
C
C Intrinsic functions called by SNHCSH:  ABS, EXP
C
C***********************************************************
C
      REAL AX, C1, C2, C3, C4, EXPX, F, XC, XS, XSD2, XSD4
C
      DATA C1/.1666666666659E0/,
     .     C2/.8333333431546E-2/,
     .     C3/.1984107350948E-3/,
     .     C4/.2768286868175E-5/
      AX = ABS(X)
      XS = AX*AX
      IF (AX .LE. .5) THEN
C
C Approximations for small X:
C
        XC = X*XS
        SINHM = XC*(((C4*XS+C3)*XS+C2)*XS+C1)
        XSD4 = .25*XS
        XSD2 = XSD4 + XSD4
        F = (((C4*XSD4+C3)*XSD4+C2)*XSD4+C1)*XSD4
        COSHMM = XSD2*F*(F+2.)
        COSHM = COSHMM + XSD2
      ELSE
C
C Approximations for large X:
C
        EXPX = EXP(AX)
        SINHM = -(((1./EXPX+AX)+AX)-EXPX)/2.
        IF (X .LT. 0.) SINHM = -SINHM
        COSHM = ((1./EXPX-2.)+EXPX)/2.
        COSHMM = COSHM - XS/2.
      ENDIF
      RETURN
      END
      REAL FUNCTION TRVOL (X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3)
      REAL X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This function computes the integral over a triangle of
C the linear (planar) surface which interpolates data
C values at the vertices.
C
C On input:
C
C       X1,X2,X3 = X coordinates of the vertices of the tri-
C                  angle in counterclockwise order.
C
C       Y1,Y2,Y3 = Y coordinates of the vertices of the tri-
C                  angle in one-to-one correspondence with
C                  X1, X2, and X3.
C
C       Z1,Z2,Z3 = Data values at the vertices (X1,Y1),
C                  (X2,Y2), (X3,Y3), respectively.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       TRVOL = Integral over the triangle of the linear
C               interpolant.  Note that TRVOL will have
C               the wrong sign if the vertices are speci-
C               fied in clockwise order.
C
C Modules required by TRVOL:  None
C
C***********************************************************
C
      REAL AREA
C
      AREA = (X2-X1)*(Y3-Y1) - (X3-X1)*(Y2-Y1)
C
C AREA is twice the (signed) area of the triangle.
C TRVOL is the mean of the data values times the area of the
C   triangle.
C
      TRVOL = (Z1 + Z2 + Z3)*AREA/6.
      RETURN
      END
      SUBROUTINE TVAL (X,Y,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,ZX1,
     .                 ZX2,ZX3,ZY1,ZY2,ZY3,DFLAG, F,FX,FY,
     .                 IER)
      INTEGER IER
      LOGICAL DFLAG
      REAL    X, Y, X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3,
     .        ZX1, ZX2, ZX3, ZY1, ZY2, ZY3, F, FX, FY
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/12/90
C
C   Given function values and first partial derivatives at
C the vertices of a triangle, along with a point P in the
C triangle, this subroutine computes an interpolated value
C F(P) and, optionally, the first partial derivatives of F
C at P.
C
C   The interpolant F of the vertex values and gradients is
C the Clough-Tocher finite element.  F is cubic in each of
C the three subtriangles of equal area obtained by joining
C the vertices to the barycenter, but has only quadratic
C precision (exact for values and partials from a quadratic
C polynomial).  Along each triangle side, F is the Hermite
C cubic interpolant of the endpoint values and tangential
C gradient components, and the normal gradient component of
C F varies linearly between the interpolated endpoint nor-
C mal components.  Thus, since values and first partials on
C a triangle side depend only on the endpoint data, the
C method results in a C-1 interpolant over a triangulation.
C Second derivatives are discontinuous across subtriangle
C boundaries.
C
C   The computational procedure, due to Charles Lawson, has
C the following operation counts:  62 adds, 54 multiplies,
C 8 divides, and 6 compares for an interpolated value, and
C 170 adds, 142 multiplies, 14 divides, and 6 compares for
C both a value and a pair of first partial derivatives.
C
C On input:
C
C       X,Y = Coordinates of the point P at which F is to
C             be evaluated.
C
C       X1,X2,X3 = X coordinates of the vertices of the tri-
C                  angle in counterclockwise order.
C
C       Y1,Y2,Y3 = Y coordinates of the vertices of the tri-
C                  angle in one-to-one correspondence with
C                  X1, X2, and X3.
C
C       Z1,Z2,Z3 = Data values at the vertices (X1,Y1),
C                  (X2,Y2), (X3,Y3), respectively.
C
C       ZX1,ZX2,ZX3 = X-derivative values at the vertices.
C
C       ZY1,ZY2,ZY3 = Y-derivative values at the vertices.
C
C       DFLAG = Logical flag which specifies whether first
C               partial derivatives at P are to be computed:
C               DFLAG = .TRUE. if and only if partials are
C               to be returned.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       F = Value of the interpolatory function at P if
C           IER = 0, or zero if IER = 1.  Note that, if
C           P is not contained in the triangle, F is an
C           extrapolated value.
C
C       FX,FY = Partial derivatives of F at P if DFLAG =
C               .TRUE. and IER = 0, unaltered otherwise.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if the vertices of the triangle are
C                     collinear.
C
C Modules required by TVAL:  None
C
C***********************************************************
C
      INTEGER I, IP1, IP2, IP3
      REAL    A(3), AREA, AX(3), AY(3), B(3), BX(3), BY(3),
     .        C(3), CX(3), CY(3), C1, C2, FF(3), G(3),
     .        GX(3), GY(3), P(3), PHI(3), PHIX(3), PHIY(3),
     .        PX(3), PY(3), Q(3), QX(3), QY(3), R(3), RMIN,
     .        RO(3), ROX(3), ROY(3), RX(3), RY(3), SL(3),
     .        U(3), V(3), XP, YP
C
C Local parameters:
C
C A(K) =            Cardinal function whose coefficient is
C                     Z(K)
C AREA =            Twice the area of the triangle
C AX(K),AY(K) =     X,Y partials of A(K) -- cardinal
C                     functions for FX and FY
C B(K) =            Twice the cardinal function whose
C                     coefficient is ZX(K)
C BX(K),BY(K) =     X,Y partials of B(K)
C C(K) =            Twice the cardinal function whose
C                     coefficient is ZY(K)
C CX(K),CY(K) =     X,Y partials of C(K)
C C1,C2 =           Factors for computing RO
C FF(K) =           Factors for computing G, GX, and GY --
C                     constant
C G(K) =            Factors for computing the cardinal
C                     functions -- cubic
C GX(K),GY(K) =     X,Y partials of G(K)
C I =               DO-loop index
C IP1,IP2,IP3 =     Permuted indexes for computing RO, ROX,
C                     and ROY
C P(K) =            G(K) + PHI(K)
C PHI(K)            R(K-1)*R(K+1) -- quadratic
C PHIX(K),PHIY(K) = X,Y partials of PHI(K)
C PX(K),PY(K) =     X,Y partials of P(K)
C Q(K) =            G(K) - PHI(K)
C QX(K),QY(K) =     X,Y partials of Q(K)
C R(K) =            K-th barycentric coordinate
C RMIN =            Min(R1,R2,R3)
C RO(K) =           Factors for computing G -- cubic
C                     correction terms
C ROX(K),ROY(K) =   X,Y partials of RO(K)
C RX(K),RY(K) =     X,Y partial derivatives of R(K)
C SL(K) =           Square of the length of the side
C                     opposite vertex K
C U(K) =            X-component of the vector representing
C                     the side opposite vertex K
C V(K) =            Y-component of the vector representing
C                     the side opposite vertex K
C XP,YP =           X-X1, Y-Y1
C
      U(1) = X3 - X2
      U(2) = X1 - X3
      U(3) = X2 - X1
C
      V(1) = Y3 - Y2
      V(2) = Y1 - Y3
      V(3) = Y2 - Y1
C
      DO 1 I = 1,3
        SL(I) = U(I)*U(I) + V(I)*V(I)
    1   CONTINUE
C
C AREA = 3->1 X 3->2.
C
      AREA = U(1)*V(2) - U(2)*V(1)
      IF (AREA .EQ. 0.) GO TO 9
      IER = 0
C
C R(1) = (2->3 X 2->P)/AREA, R(2) = (1->P X 1->3)/AREA,
C   R(3) = (1->2 X 1->P)/AREA.
C
      R(1) = (U(1)*(Y-Y2) - V(1)*(X-X2))/AREA
      XP = X - X1
      YP = Y - Y1
      R(2) = (U(2)*YP - V(2)*XP)/AREA
      R(3) = (U(3)*YP - V(3)*XP)/AREA
C
      PHI(1) = R(2)*R(3)
      PHI(2) = R(3)*R(1)
      PHI(3) = R(1)*R(2)
C
      IF (R(1) .GT. R(2)  .OR.  R(1) .GT. R(3)) GO TO 3
      RMIN = R(1)
      IP1 = 1
      IP2 = 2
      IP3 = 3
      GO TO 5
    3 IF (R(2) .GT. R(3)) GO TO 4
      RMIN = R(2)
      IP1 = 2
      IP2 = 3
      IP3 = 1
      GO TO 5
    4 RMIN = R(3)
      IP1 = 3
      IP2 = 1
      IP3 = 2
C
    5 C1 = RMIN*RMIN/2.
      C2 = RMIN/3.
      RO(IP1) = (PHI(IP1) + 5.*C1/3.)*R(IP1) - C1
      RO(IP2) = C1*(R(IP3) - C2)
      RO(IP3) = C1*(R(IP2) - C2)
C
      FF(1) = 3.*(SL(2)-SL(3))/SL(1)
      FF(2) = 3.*(SL(3)-SL(1))/SL(2)
      FF(3) = 3.*(SL(1)-SL(2))/SL(3)
C
      G(1) = (R(2)-R(3))*PHI(1) + FF(1)*RO(1) - RO(2)+RO(3)
      G(2) = (R(3)-R(1))*PHI(2) + FF(2)*RO(2) - RO(3)+RO(1)
      G(3) = (R(1)-R(2))*PHI(3) + FF(3)*RO(3) - RO(1)+RO(2)
C
      DO 6 I = 1,3
        P(I) = G(I) + PHI(I)
        Q(I) = G(I) - PHI(I)
    6   CONTINUE
C
      A(1) = R(1) + G(3) - G(2)
      A(2) = R(2) + G(1) - G(3)
      A(3) = R(3) + G(2) - G(1)
C
      B(1) = U(3)*P(3) + U(2)*Q(2)
      B(2) = U(1)*P(1) + U(3)*Q(3)
      B(3) = U(2)*P(2) + U(1)*Q(1)
C
      C(1) = V(3)*P(3) + V(2)*Q(2)
      C(2) = V(1)*P(1) + V(3)*Q(3)
      C(3) = V(2)*P(2) + V(1)*Q(1)
C
C F is a linear combination of the cardinal functions.
C
      F = A(1)*Z1 + A(2)*Z2 + A(3)*Z3 + (B(1)*ZX1 + B(2)*ZX2
     .    + B(3)*ZX3 + C(1)*ZY1 + C(2)*ZY2 + C(3)*ZY3)/2.
      IF (.NOT. DFLAG) RETURN
C
C Compute FX and FY.
C
      DO 7 I = 1,3
        RX(I) = -V(I)/AREA
        RY(I) = U(I)/AREA
    7   CONTINUE
C
      PHIX(1) = R(2)*RX(3) + RX(2)*R(3)
      PHIY(1) = R(2)*RY(3) + RY(2)*R(3)
      PHIX(2) = R(3)*RX(1) + RX(3)*R(1)
      PHIY(2) = R(3)*RY(1) + RY(3)*R(1)
      PHIX(3) = R(1)*RX(2) + RX(1)*R(2)
      PHIY(3) = R(1)*RY(2) + RY(1)*R(2)
C
      ROX(IP1) = RX(IP1)*(PHI(IP1) + 5.*C1) +
     .           R(IP1)*(PHIX(IP1) - RX(IP1))
      ROY(IP1) = RY(IP1)*(PHI(IP1) + 5.*C1) +
     .           R(IP1)*(PHIY(IP1) - RY(IP1))
      ROX(IP2) = RX(IP1)*(PHI(IP2) - C1) + C1*RX(IP3)
      ROY(IP2) = RY(IP1)*(PHI(IP2) - C1) + C1*RY(IP3)
      ROX(IP3) = RX(IP1)*(PHI(IP3) - C1) + C1*RX(IP2)
      ROY(IP3) = RY(IP1)*(PHI(IP3) - C1) + C1*RY(IP2)
C
      GX(1) = (RX(2) - RX(3))*PHI(1) + (R(2) - R(3))*PHIX(1)
     .        + FF(1)*ROX(1) - ROX(2) + ROX(3)
      GY(1) = (RY(2) - RY(3))*PHI(1) + (R(2) - R(3))*PHIY(1)
     .        + FF(1)*ROY(1) - ROY(2) + ROY(3)
      GX(2) = (RX(3) - RX(1))*PHI(2) + (R(3) - R(1))*PHIX(2)
     .        + FF(2)*ROX(2) - ROX(3) + ROX(1)
      GY(2) = (RY(3) - RY(1))*PHI(2) + (R(3) - R(1))*PHIY(2)
     .        + FF(2)*ROY(2) - ROY(3) + ROY(1)
      GX(3) = (RX(1) - RX(2))*PHI(3) + (R(1) - R(2))*PHIX(3)
     .        + FF(3)*ROX(3) - ROX(1) + ROX(2)
      GY(3) = (RY(1) - RY(2))*PHI(3) + (R(1) - R(2))*PHIY(3)
     .        + FF(3)*ROY(3) - ROY(1) + ROY(2)
C
      DO 8 I = 1,3
        PX(I) = GX(I) + PHIX(I)
        PY(I) = GY(I) + PHIY(I)
        QX(I) = GX(I) - PHIX(I)
        QY(I) = GY(I) - PHIY(I)
    8   CONTINUE
C
      AX(1) = RX(1) + GX(3) - GX(2)
      AY(1) = RY(1) + GY(3) - GY(2)
      AX(2) = RX(2) + GX(1) - GX(3)
      AY(2) = RY(2) + GY(1) - GY(3)
      AX(3) = RX(3) + GX(2) - GX(1)
      AY(3) = RY(3) + GY(2) - GY(1)
C
      BX(1) = U(3)*PX(3) + U(2)*QX(2)
      BY(1) = U(3)*PY(3) + U(2)*QY(2)
      BX(2) = U(1)*PX(1) + U(3)*QX(3)
      BY(2) = U(1)*PY(1) + U(3)*QY(3)
      BX(3) = U(2)*PX(2) + U(1)*QX(1)
      BY(3) = U(2)*PY(2) + U(1)*QY(1)
C
      CX(1) = V(3)*PX(3) + V(2)*QX(2)
      CY(1) = V(3)*PY(3) + V(2)*QY(2)
      CX(2) = V(1)*PX(1) + V(3)*QX(3)
      CY(2) = V(1)*PY(1) + V(3)*QY(3)
      CX(3) = V(2)*PX(2) + V(1)*QX(1)
      CY(3) = V(2)*PY(2) + V(1)*QY(1)
C
C FX and FY are linear combinations of the cardinal
C   functions.
C
      FX = AX(1)*Z1 + AX(2)*Z2 + AX(3)*Z3 + (BX(1)*ZX1 +
     .     BX(2)*ZX2 + BX(3)*ZX3 + CX(1)*ZY1 + CX(2)*ZY2 +
     .     CX(3)*ZY3)/2.
      FY = AY(1)*Z1 + AY(2)*Z2 + AY(3)*Z3 + (BY(1)*ZX1 +
     .     BY(2)*ZX2 + BY(3)*ZX3 + CY(1)*ZY1 + CY(2)*ZY2 +
     .     CY(3)*ZY3)/2.
      RETURN
C
C The vertices are collinear.
C
    9 IER = 1
      F = 0.
      RETURN
      END
      SUBROUTINE UNIF (NCC,LCC,N,X,Y,Z,GRAD,LIST,LPTR,LEND,
     .                 IFLGS,SIGMA,NROW,NX,NY,PX,PY,SFLAG,
     .                 SVAL, ZZ,IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        IFLGS, NROW, NX, NY, IER
      LOGICAL SFLAG
      REAL    X(N), Y(N), Z(N), GRAD(2,N), SIGMA(*), PX(NX),
     .        PY(NY), SVAL, ZZ(NROW,NY)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/22/97
C
C   Given a Delaunay triangulation of a set of points in the
C plane with associated data values and gradients, this sub-
C routine interpolates the data to a set of rectangular grid
C points for such applications as contouring.  Extrapolation
C is performed at grid points exterior to the triangulation,
C and the interpolant is once-continuously differentiable
C over the entire plane.  Refer to Subroutine INTRC1 for
C further details.
C
C On input:
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index of the
C             first node of constraint I in LCC(I).  For I =
C             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       Z = Array of length N containing data values at the
C           nodes.  Refer to Subroutines ZGRADG and ZGRADL.
C
C       GRAD = 2 by N array whose columns contain estimated
C              gradients at the nodes with X partial deriva-
C              tives in the first row and Y partials in the
C              second.  Refer to Subroutines GRADC, GRADG,
C              GRADL, SMSURF, ZGRADG, and ZGRADL.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C       IFLGS = Tension factor option:
C               IFLGS .LE. 0 if a single uniform tension
C                            factor is to be used.
C               IFLGS .GE. 1 if variable tension is desired.
C
C       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               routines FVAL, GETSIG, SIG0, SIG1, and SIG2.
C
C       NROW = Number of rows in the dimension statement of
C              ZZ.
C
C       NX,NY = Number of rows and columns, respectively, in
C               the rectangular grid.  1 .LE. NX .LE. NROW,
C               and 1 .LE. NY.
C
C       PX,PY = Arrays of length NX and NY, respectively,
C               containing the coordinates of the grid
C               lines.
C
C       SFLAG = Special value flag:
C               SFLAG = .FALSE. if special values are not to
C                               be used (ZZ contains only
C                               interpolated or extrapolated
C                               values.
C               SFLAG = .TRUE. if SVAL is to be stored in ZZ
C                              elements corresponding to
C                              grid points which lie in a
C                              constraint region.
C
C       SVAL = Special value for grid points lying in a con-
C              straint region, or dummy parameter if SFLAG =
C              .FALSE.
C
C The above parameters are not altered by this routine.
C
C       ZZ = NROW by NCOL array for some NCOL .GE. NY.
C
C On output:
C
C       ZZ = Interpolated values at the grid points (or
C            special values) if IER .GE. 0.  ZZ(I,J) =
C            F(PX(I),PY(J)) for I = 1,...,NX and J = 1,...,
C            NY, where F is the interpolatory surface.
C
C       IER = Error indicator:
C             IER .GE. 0 if no errors were encountered.
C                        IER contains the number of grid
C                        points exterior to the triangula-
C                        tion or contained in a constraint
C                        region triangle (extrapolated
C                        values).
C             IER = -1 if NCC, N, NROW, NX, or NY is
C                      outside its valid range on input.
C                      LCC is not tested for validity.
C             IER = -2 if the nodes are collinear or the
C                      triangulation is invalid.
C
C TRIPACK modules required by UNIF:  CRTRI, JRAND, LEFT,
C                                      LSTPTR, TRFIND
C
C SRFPACK modules required by UNIF:  ARCINT, COORDS, FVAL,
C                                      INTRC1, SNHCSH, TVAL
C
C***********************************************************
C
      INTEGER I, IERR, IST, J, NEX, NI, NJ, NST
      LOGICAL DFLAG, SFL
      REAL    DUM
      DATA    DFLAG/.FALSE./,  NST/1/
C
C Local parameters:
C
C DFLAG = Derivative flag for INTRC1
C DUM =   Dummy INTRC1 parameter
C I,J =   DO-loop indexes
C IERR =  Error flag for calls to INTRC1
C IST =   Parameter for INTRC1
C NEX =   Number of grid points exterior to the triangula-
C           tion boundary (number of extrapolated values)
C NI,NJ = Local copies of NX and NY
C NST =   Initial value for IST
C SFL =   Local copy of SFLAG
C
      NI = NX
      NJ = NY
      IF (NCC .LT. 0  .OR.  N .LT. 3  .OR.  NI .LT. 1  .OR.
     .    NI .GT. NROW  .OR.  NJ .LT. 1) GO TO 3
      SFL = SFLAG
      IST = NST
C
C Compute interpolated values.
C
      NEX = 0
      DO 2 J = 1,NJ
        DO 1 I = 1,NI
          CALL INTRC1 (PX(I),PY(J),NCC,LCC,N,X,Y,Z,LIST,
     .                 LPTR,LEND,IFLGS,SIGMA,GRAD,
     .                 DFLAG, IST, ZZ(I,J),DUM,DUM,IERR)
          IF (IERR .LT. 0) GO TO 4
          IF (IERR .GT. 0) NEX = NEX + 1
          IF (SFL  .AND. IERR .EQ. 1) ZZ(I,J) = SVAL
    1     CONTINUE
    2   CONTINUE
      IER = NEX
      RETURN
C
C Invalid input parameter.
C
    3 IER = -1
      RETURN
C
C Triangulation nodes are collinear.
C
    4 IER = -2
      RETURN
      END
      REAL FUNCTION VOLUME (NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N)
      REAL    X(N), Y(N), Z(N)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   08/26/91
C
C   Given a triangulation of a set of N nodes, along with
C data values at the nodes, this function computes the int-
C egral over a region R of the piecewise linear interpolant
C of the data values.  R is the convex hull of the nodes
C with constraint regions excluded.
C
C On input:
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index of the
C             first node of constraint I in LCC(I).  For I =
C             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       Z = Array of length N containing data values at the
C           nodes.  Refer to Subroutine ZGRADL.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       VOLUME = Sum of the volumes of the linear interpo-
C                lants on the non-constraint triangles, or
C                zero if a parameter is outside its valid
C                range on input.
C
C SRFPACK module required by VOLUME:  TRVOL
C
C Intrinsic function called by VOLUME:  ABS
C
C***********************************************************
C
      REAL    TRVOL
      INTEGER I, ILAST, LCC1, LP2, LP3, LPL, N1, N2, N3,
     .        NM2, NN
      REAL    SUM, XN1, YN1, ZN1
C
C Test for invalid input parameters.
C
      IF (NCC .LT. 0) GO TO 5
      NN = N
      LCC1 = NN+1
      IF (NCC .EQ. 0) THEN
        IF (NN .LT. 3) GO TO 5
      ELSE
        DO 1 I = NCC,1,-1
          IF (LCC1-LCC(I) .LT. 3) GO TO 5
          LCC1 = LCC(I)
    1     CONTINUE
        IF (LCC1 .LT. 1) GO TO 5
      ENDIF
C
C Initialize for loop on triangles (N1,N2,N3) such that N2
C   and N3 have larger indexes than N1.  SUM contains the
C   accumulated volume, I is the index of the constraint
C   containing N1 if N1 is a constraint node, and ILAST is
C   the last node of constraint I.
C
      I = 0
      ILAST = LCC1 - 1
      SUM = 0.
      NM2 = NN - 2
      DO 4 N1 = 1,NM2
        XN1 = X(N1)
        YN1 = Y(N1)
        ZN1 = Z(N1)
        IF (N1 .GT. ILAST) THEN
          I = I + 1
          IF (I .LT. NCC) THEN
            ILAST = LCC(I+1) - 1
          ELSE
            ILAST = NN
          ENDIF
        ENDIF
C
C Top of loop on neighbors of N1.
C
        LPL = LEND(N1)
        LP2 = LPL
    2   LP2 = LPTR(LP2)
          N2 = LIST(LP2)
          LP3 = LPTR(LP2)
          N3 = ABS(LIST(LP3))
          IF (N2 .LT. N1  .OR.  N3 .LT. N1) GO TO 3
C
C   (N1,N2,N3) lies in a constraint region iff the vertices
C     are nodes of the same constraint and N2 < N3.
C
          IF (N1 .LT. LCC1  .OR.  N2 .GT. N3  .OR.
     .        N3 .GT. ILAST) THEN
            SUM = SUM + TRVOL(XN1,X(N2),X(N3),YN1,Y(N2),
     .                        Y(N3),ZN1,Z(N2),Z(N3))
          ENDIF
C
C   Bottom of loop on neighbors.
C
    3     IF (LP2 .NE. LPL) GO TO 2
    4   CONTINUE
C
      VOLUME = SUM
      RETURN
C
C Invalid input parameter.
C
    5 VOLUME = 0.
      RETURN
      END
      SUBROUTINE ZGRADG (NCC,LCC,N,X,Y,LIST,LPTR,LEND,
     .                   IFLGS,SIGMA, NIT,DZMAX,Z,
     .                   GRAD, IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        IFLGS, NIT, IER
      REAL    X(N), Y(N), SIGMA(*), DZMAX, Z(N), GRAD(2,N)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/22/97
C
C   Given a triangulation of N nodes, along with data values
C at non-constraint nodes, this subroutine employs a global
C method to compute estimated gradients at the nodes and
C approximate data values at the constraint nodes.  For
C NCN = N-LCC(1)+1 constraint nodes, the method consists of
C minimizing a quadratic functional Q(U) over vectors U of
C length 2N+NCN containing gradients and values at con-
C straint nodes.  Q is taken to be the sum over the
C triangulation arcs, excluding interior constraint arcs,
C of the linearized curvature (integral of squared second
C derivative) of the Hermite interpolatory tension spline
C defined by the data values and tangential gradient comp-
C onents at the endpoints of the arc.
C
C   This minimization problem corresponds to a symmetric
C positive definite sparse linear system which is solved by
C a block Gauss-Seidel method with N blocks of order 2, or
C order 3 for constraint nodes.
C
C   An alternative method (Subroutine ZGRADL) computes a
C local approximation to the gradient and data value (if not
C specified) at a single node.  The relative speed and
C accuracy of the two methods depends on the distribution
C of the nodes.  Relative accuracy also depends on the data
C values.
C
C   Note that the call to ZGRADG can be followed by a call
C to GRADG in order to compute improved gradient estimates
C with fixed data values.  This is recommended for improved
C efficiency.
C
C On input:
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC > 0.
C
C       LCC = Array of length NCC containing the index of
C             the first node of constraint I in LCC(I).  For
C             I = 1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1, and LCC(1) .GE. 4.
C
C       N = Number of nodes in the triangulation.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C       IFLGS = Tension factor option:
C               IFLGS .LE. 0 if a single uniform tension
C                            factor is to be used.
C               IFLGS .GE. 1 if variable tension is desired.
C
C       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               routines GETSIG, SIG0, SIG1, and SIG2.
C
C The above parameters are not altered by this routine.
C
C       NIT = Maximum number of Gauss-Seidel iterations to
C             be employed.  This maximum will likely be
C             achieved if DZMAX is smaller than the machine
C             precision.  Note that complete convergence is
C             not necessary to achieve maximum accuracy of
C             the interpolant.  NIT > 0.
C
C       DZMAX - Nonnegative convergence criterion.  The
C               method is terminated when the maximum change
C               in a solution Z-component between iterations
C               is at most DZMAX.  The change in a solution
C               component is taken to be the magnitude of
C               the difference relative to 1 plus the magni-
C               tude of the previous value.
C
C       Z = Array of length N containing data values in the
C           first LCC(1)-1 locations and initial solution
C           estimates in the remaining locations.  Zeros are
C           sufficient, but Subroutine ZINIT may be called
C           to provide better initial estimates.
C
C       GRAD = 2 by N array whose columns contain initial
C              estimates of the gradients with X partial
C              derivatives in the first row, Y partials in
C              the second.  Zeros are sufficient.
C
C On output:
C
C       NIT = Number of Gauss-Seidel iterations employed.
C
C       DZMAX = Maximum relative change in a solution Z-
C               component at the last iteration.
C
C       Z = Array updated with approximate data values in
C           the last NCN = N-LCC(1)+1 locations if IER .GE.
C           0.  Z is not altered if IER = -1.
C
C       GRAD = Estimated gradients at the nodes if IER .GE.
C              0.  GRAD is not altered if IER = -1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     convergence criterion was achieved.
C             IER = 1 if no errors were encountered but con-
C                     vergence was not achieved within NIT
C                     iterations.
C             IER = -1 if NCC, an LCC entry, N, NIT, or
C                      DZMAX is outside its valid range
C                      on input.
C             IER = -2 if all nodes are collinear or the
C                      triangulation data structure is in-
C                      valid.
C             IER = -3 if duplicate nodes were encountered.
C
C SRFPACK modules required by ZGRADG:  GRCOEF, SNHCSH
C
C Intrinsic functions called by ZGRADG:  ABS, MAX, SQRT
C
C***********************************************************
C
      INTEGER I, IFL, IFRST, ILAST, ITER, J, JN, K, KBAK,
     .        KFOR, LCC1, LP, LPF, LPJ, LPL, LPN, MAXIT, NB,
     .        NN
      REAL    A11, A12, A13, A22, A23, A33, AREAJ, AREAN,
     .        AREAP, D, DCUB, DF, DSQ, DX, DXS, DY, DYS, DZ,
     .        DZJ, DZK, DZMX, DZX, DZY, R1, R2, R3, SDF,
     .        SIG, T, TOL, W, XK, YK, ZK, ZXK, ZYK
C
      NN = N
      IFL = IFLGS
      MAXIT = NIT
      TOL = DZMAX
C
C Test for errors in input parameters.
C
      IF (NCC .LE. 0  .OR.  MAXIT .LT. 1  .OR.  TOL .LT. 0.)
     .  GO TO 9
      LCC1 = NN+1
      DO 1 I = NCC,1,-1
        IF (LCC1-LCC(I) .LT. 3) GO TO 9
        LCC1 = LCC(I)
    1   CONTINUE
      IF (LCC1 .LT. 4) GO TO 9
C
C Initialize iteration count and SIG (overwritten if
C   IFLGS > 0).
C
      ITER = 0
      SIG = SIGMA(1)
C
C Top of iteration loop:  If K is a constraint node, I
C   indexes the constraint containing node K, IFRST and
C   ILAST are the first and last nodes of constraint I, and
C   (KBAK,K,KFOR) is a subsequence of constraint I.
C
    2 IF (ITER .EQ. MAXIT) GO TO 8
      DZMX = 0.
      I = 0
      ILAST = LCC1-1
      KBAK = 0
      KFOR = 0
C
C Loop on nodes.
C
      DO 7 K = 1,NN
        IF (K .GE. LCC1) THEN
          IF (K .GT. ILAST) THEN
            I = I + 1
            IFRST = K
            IF (I .LT. NCC) THEN
              ILAST = LCC(I+1) - 1
            ELSE
              ILAST = NN
            ENDIF
            KBAK = ILAST
            KFOR = K + 1
          ELSE
            KBAK = K-1
            IF (K .LT. ILAST) THEN
              KFOR = K+1
            ELSE
              KFOR = IFRST
            ENDIF
          ENDIF
        ENDIF
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        ZXK = GRAD(1,K)
        ZYK = GRAD(2,K)
C
C Initialize components of the 2 by 2 (or 3 by 3) block --
C   symmetric matrix in A and residual in R.  The unknowns
C   are ordered (DZX,DZY,DZ).
C
        A11 = 0.
        A12 = 0.
        A13 = 0.
        A22 = 0.
        A23 = 0.
        A33 = 0.
        R1 = 0.
        R2 = 0.
        R3 = 0.
C
C Loop on neighbors J of node K.  The equation associated
C   with K->J (and hence its contribution to the functional)
C   is weighted by AREAJ/D, where AREAJ is twice the sum of
C   the areas of the triangles containing K-J (excluding
C   those which lie in a constraint region) and D is the arc
C   length.  JN is the neighbor of K following J.  AREAP is
C   to the right of K->J and AREAN is to the left.
C
        LPL = LEND(K)
        J = LIST(LPL)
        LPF = LPTR(LPL)
        JN = LIST(LPF)
        AREAN = 0.
        IF (J .GT. 0) AREAN = (X(J)-XK)*(Y(JN)-YK) -
     .                        (Y(J)-YK)*(X(JN)-XK)
        LPN = LPF
C
C Top of loop:  LPF and LPL point to the first and last
C   neighbors of K, and LPN points to JN.
C
    3   LPJ = LPN
          LPN = LPTR(LPN)
          J = JN
          AREAP = AREAN
          JN = ABS(LIST(LPN))
C
C Arc K-J lies in a constraint region and is bypassed iff K
C   and J are nodes in the same constraint and J follows
C   KFOR and precedes KBAK as a neighbor of K.
C
          IF (K .LT. LCC1  .OR.  J .LT. IFRST  .OR.
     .        J .GT. ILAST) GO TO 5
          IF (J .EQ. KBAK) AREAP = 0.
          IF (J .EQ. KBAK  .OR.  J .EQ. KFOR) GO TO 5
C
          LP = LPN
    4     NB = ABS(LIST(LP))
            IF (NB .EQ. KFOR) GO TO 5
            IF (NB .EQ. KBAK) GO TO 6
            LP = LPTR(LP)
            GO TO 4
C
C   Compute parameters associated with the edge K->J, and
C     test for duplicate nodes.  Note that AREAJ = 0 and
C     K->J is bypassed if K-J is both a constraint arc and
C     a boundary arc of the triangulation.
C
    5     DX = X(J) - XK
          DY = Y(J) - YK
          AREAN = 0.
          IF (LIST(LPL) .NE. -J  .AND.  J .NE. KFOR) AREAN =
     .      DX*(Y(JN)-YK) - DY*(X(JN)-XK)
          AREAJ = AREAP + AREAN
          IF (AREAJ .EQ. 0.) GO TO 6
          DXS = DX*DX
          DYS = DY*DY
          DSQ = DXS + DYS
          D = SQRT(DSQ)
          DCUB = D*DSQ
          IF (D .EQ. 0.) GO TO 11
          IF (IFL .GE. 1) SIG = SIGMA(LPJ)
          CALL GRCOEF (SIG,DCUB, DF,SDF)
          W = AREAJ/D
C
C   Update the 2 by 2 system components for node J.
C
          A11 = A11 + DF*DXS*W
          A12 = A12 + DF*DX*DY*W
          A22 = A22 + DF*DYS*W
          DZ = Z(J) - ZK
          DZJ = GRAD(1,J)*DX + GRAD(2,J)*DY
          DZK = ZXK*DX + ZYK*DY
          T = ( (DF+SDF)*DZ - SDF*DZJ - DF*DZK )*W
          R1 = R1 + T*DX
          R2 = R2 + T*DY
          IF (K .GE. LCC1) THEN
C
C   K is a constraint node.  Update the remaining components.
C
            W = (DF+SDF)*W
            A13 = A13 + DX*W
            A23 = A23 + DY*W
            A33 = A33 + 2.0*W
            R3 = R3 + (2.0*DZ - DZJ - DZK)*W
          ENDIF
C
C   Bottom of loop on J.
C
    6     IF (LPN .NE. LPF) GO TO 3
C
C Solve the linear system associated with the K-th block.
C
        A22 = A11*A22 - A12*A12
        R2 = A11*R2 - A12*R1
        IF (A11 .EQ. 0.  .OR.  A22 .EQ. 0.) GO TO 10
        IF (K .GE. LCC1) THEN
          A23 = A11*A23 - A12*A13
          A33 = A22*(A11*A33 - A13*A13) - A23*A23
          R3 = A22*(A11*R3 - A13*R1) - A23*R2
          IF (A33 .EQ. 0.) GO TO 10
          DZ = R3/A33
        ENDIF
        DZY = (R2 - A23*DZ)/A22
        DZX = (R1 - A12*DZY - A13*DZ)/A11
C
C Update the solution components for node K and the maxi-
C   mum relative change DZMX.
C
        GRAD(1,K) = ZXK + DZX
        GRAD(2,K) = ZYK + DZY
        IF (K .GE. LCC1) THEN
          Z(K) = ZK + DZ
          DZMX = MAX(DZMX,ABS(DZ)/(1.+ABS(ZK)))
        ENDIF
    7   CONTINUE
C
C Increment ITER and test for convergence.
C
      ITER = ITER + 1
      IF (DZMX .GT. TOL) GO TO 2
C
C Method converged.
C
      NIT = ITER
      DZMAX = DZMX
      IER = 0
      RETURN
C
C Method failed to converge within NIT iterations.
C
    8 DZMAX = DZMX
      IER = 1
      RETURN
C
C Invalid input parameter.
C
    9 NIT = 0
      DZMAX = 0.
      IER = -1
      RETURN
C
C Node K and its neighbors are collinear, resulting in a
C   singular system.
C
   10 NIT = 0
      DZMAX = DZMX
      IER = -2
      RETURN
C
C Nodes J and K coincide.
C
   11 NIT = 0
      DZMAX = DZMX
      IER = -3
      RETURN
      END
      SUBROUTINE ZGRADL (K,NCC,LCC,N,X,Y,LIST,LPTR,
     .                   LEND, NDV,Z,NPTS,DS, DX,DY,IER)
      INTEGER K, NCC, LCC(*), N, LIST(*), LPTR(*),
     .        LEND(N), NDV, NPTS(*), IER
      REAL    X(N), Y(N), Z(N), DS(*), DX, DY
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/22/97
C
C   Given a Delaunay triangulation of N nodes, along with
C data values at non-constraint nodes, this subroutine com-
C putes an estimated gradient and, if K is a constraint
C node, an approximate data value, at node K.  The values
C are taken from a quadratic function which fits the data
C values at a set of non-constraint nodes close to K in a
C weighted least squares sense.  If K is not a constraint
C node, the fitting function interpolates the data value at
C node K.  If there are fewer than six data values (non-
C constraint nodes), a linear fitting function is used.
C Also, a Marquardt stabilization factor is used if neces-
C sary to ensure a well-conditioned system.  Thus, a unique
C solution exists unless the non-constraint nodes are col-
C linear.
C
C   An alternative routine, ZGRADG, employs a global method
C to compute values at constraint nodes and gradients at all
C of the nodes at once.  The relative speed and accuracy of
C the two methods depends on whether or not constraints are
C present and on the distribution of the nodes.  Relative
C accuracy also depends on the data values.
C
C   This subroutine may be used for the following purposes:
C
C 1)  to compute gradient estimates and constraint node
C       values for INTRC1 or UNIF,
C 2)  to fill in missing constraint node values for linear
C       interpolation (INTRC0 and VOLUME) or smoothing
C       (SMSURF), and
C 3)  to provide initial estimates for GRADG or ZGRADG
C       (probably a waste of computing time).
C
C If data values at the constraint nodes are known, Subrou-
C tine GRADC or GRADL should be used in place of this
C routine.  If there are no constraints this routine differs
C from GRADL only in providing more flexibility (refer to
C NDV below).
C
C On input:
C
C       K = Index of the node at which a gradient is to be
C           computed.  1 .LE. K .LE. N.
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index of the
C             first node of constraint I in LCC(I).  For I =
C             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1, and LCC(1) .GE. 4.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C The above parameters are not altered by this routine.
C
C       NDV = Number of data values (including Z(K) if K is
C             not a constraint node) to be used in the least
C             squares fit, unless more equations are re-
C             quired for stability.  3 .LE. NDV .LT. LCC(1).
C             Note that a linear fitting function will be
C             used if NDV < 6 on input.  A reasonable value
C             is NDV = 9.
C
C       Z = Array of length N containing data values in the
C           first LCC(1)-1 locations.
C
C       NPTS,DS = Arrays of length at least Min(L+1,N) where
C                 L is defined below.  (Length N is suffi-
C                 cient.)
C
C On output:
C
C       NDV = Number of data values (non-constraint nodes)
C             used in the least squares fit.
C
C       Z = Array updated with an approximate data value at
C           node K in Z(K) if K .GE. LCC(1) and IER = 0.
C
C       NPTS = Array containing the indexes of the ordered
C              sequence of L closest nodes to node K (with K
C              in the first position) unless IER .NE. 0.  L
C              is the smallest integer such that the se-
C              quence contains NDV (output value) non-
C              constraint nodes.  NPTS(L+1) = 0 if L < N.
C
C       DS = Array containing the distance between node K
C            and NPTS(I) in DS(I) for I = 1,...,L unless
C            IER .NE. 0.  Distance is measured within the
C            non-constraint region (refer to Subroutine
C            GETNP).
C
C       DX,DY = Estimated X and Y partial derivatives at
C               node K unless IER .NE. 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if K, NCC, an LCC entry, N, or NDV is
C                     outside its valid range on input.
C             IER = 2 if all non-constraint nodes are col-
C                     linear.
C
C TRIPACK modules required by ZGRADL:  GETNP, INTSEC
C
C SRFPACK modules required by ZGRADL:  GIVENS, ROTATE,
C                                        SETRO2
C
C Intrinsic functions called by ZGRADL:  ABS, MIN
C
C***********************************************************
C
      INTEGER I, IERR, IR, IROW1, J, JP1, JR, KK, L, LCC1,
     .        LNP, LR, ND, NDMIN, NP, NPAR, NPM1, NPP1
      REAL    A(7,7), C, DMIN, DTOL, RFAC, RIN, S, SF, SFS,
     .        STF, W, XK, YK, ZK
      LOGICAL INIT, STABL
      DATA    RFAC/1.05/,  DTOL/.01/
C
C Store parameters in local variables, test for errors, and
C   initialize switches.
C
      KK = K
      IF (NCC .GT. 0) THEN
        LCC1 = LCC(1)
      ELSE
        LCC1 = N+1
      ENDIF
      NDMIN = NDV
      IF (KK .LT. 1  .OR.  KK .GT. N  .OR.  NCC .LT. 0
     .    .OR.  LCC1 .LT. 4  .OR.  NDMIN .LT. 3  .OR.
     .    NDMIN .GE. LCC1) GO TO 13
      XK = X(KK)
      YK = Y(KK)
      ZK = 0.
      IF (KK .LT. LCC1) ZK = Z(KK)
      INIT = .FALSE.
      STABL = .FALSE.
C
C Set NPTS to the closest LNP nodes to K, where LNP is the
C   smallest integer such that NPTS contains NDMIN non-
C   constraint nodes.  ND is the number of non-constraint
C   nodes currently in NPTS.
C
      LNP = 1
      NPTS(1) = KK
      DS(1) = 0.
      ND = 0
      IF (KK .LT. LCC1) ND = 1
C
C   Get a new non-constraint node.
C
    1 LNP = LNP + 1
      CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,LNP, NPTS,
     .            DS, IERR)
      IF (IERR .NE. 0) GO TO 13
      IF (NPTS(LNP) .GE. LCC1) GO TO 1
      ND = ND + 1
      IF (ND .LT. NDMIN) GO TO 1
C
C Compute an inverse radius of influence to be used in the
C   weights, and test the initialization switch -- INIT =
C   .TRUE. iff A has been initialized with the first NPAR
C   equations.
C
      RIN = 1./(RFAC*DS(LNP))
      IF (INIT) GO TO 5
C
C A Q-R decomposition is used to solve the least squares
C   system.  For a quadratic fit there are NPAR = 5 or
C   NPAR = 6 parameters, depending on whether or not K is
C   a constraint node.  (At least 6 data values are needed
C   in either case.)  The transpose of the augmented regres-
C   sion matrix is stored in A with columns (rows of A) de-
C   fined as follows -- 1-3 are the quadratic terms, 4 and 5
C   are the linear terms with coefficients DX and DY, column
C   6 is the constant term with coefficient Z(K) (extraneous
C   if NPAR = 5), and the last column is the right hand
C   side.  In the case of a linear fit, the first 3 columns
C   are ignored and the first 3 rows are omitted.  The lin-
C   ear terms are scaled by SF = 1/DMAX, where DMAX is the
C   maximum distance between K and a non-constraint node in
C   NPTS, and the quadratic terms are scaled by SF**2.
C
      SF = 1./DS(LNP)
      SFS = SF*SF
      IROW1 = 1
      IF (ND .LT. 6) IROW1 = 4
      NPAR = 5
      IF (KK .GE. LCC1) NPAR = 6
      NPM1 = NPAR - 1
      NPP1 = NPAR + 1
C
C Set up the first NPAR equations and zero out the lower
C   triangle (upper triangle of A) with Givens rotations --
C
      L = 1
      DO 4 IR = IROW1,NPAR
    2   L = L + 1
          NP = NPTS(L)
          IF (NP .GE. LCC1) GO TO 2
        W = 1./DS(L) - RIN
        CALL SETRO2 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,
     .               W, A(1,IR))
        IF (IR .EQ. IROW1) GO TO 4
        DO 3 JR = IROW1,IR-1
          JP1 = JR + 1
          LR = 7 - JR
          CALL GIVENS (A(JR,JR),A(JR,IR),C,S)
          CALL ROTATE (LR,C,S,A(JP1,JR),A(JP1,IR))
    3     CONTINUE
    4   CONTINUE
      INIT = .TRUE.
C
C Incorporate additional equations into the system using the
C   last column of A (or next to last if NPAR = 5).
C
    5 IF (L .EQ. LNP) GO TO 7
        L = L + 1
        NP = NPTS(L)
        IF (NP .GE. LCC1) GO TO 5
        W = 1./DS(L) - RIN
        CALL SETRO2 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,
     .               W, A(1,NPP1))
        DO 6 JR = IROW1,NPAR
          JP1 = JR + 1
          LR = 7 - JR
          CALL GIVENS (A(JR,JR),A(JR,NPP1),C,S)
          CALL ROTATE (LR,C,S,A(JP1,JR),A(JP1,NPP1))
    6     CONTINUE
        GO TO 5
C
C Test the system for ill-conditioning.
C
    7 DMIN = ABS(A(NPAR,NPAR))
      DO 8 I = IROW1,NPM1
        DMIN = MIN(DMIN,ABS(A(I,I)))
    8   CONTINUE
      IF (DMIN/W .GE. DTOL) GO TO 12
      IF (ND .LT. LCC1) THEN
C
C   Add another equation to the system and increase the
C     radius R.
C
        NDMIN = NDMIN + 1
        GO TO 1
      ENDIF
C
C The system is ill-conditioned and all non-constraint nodes
C   have been used.  Stabilize the system by damping out the
C   second partials (coefficients of the quadratic terms)
C   unless the system has already been stabilized or a
C   linear fitting function is being used.  Add multiples
C   of the first 3 unit vectors to the first 3 equations.
C
      IF (STABL  .OR.  IROW1 .EQ. 4) GO TO 14
      STF = W
      DO 11 I = 1,3
        A(I,NPP1) = STF
        DO 9 J = I+1,7
          A(J,NPP1) = 0.
    9     CONTINUE
        DO 10 JR = I,NPAR
          JP1 = JR + 1
          LR = 7 - JR
          CALL GIVENS (A(JR,JR),A(JR,NPP1),C,S)
          CALL ROTATE (LR,C,S,A(JP1,JR),A(JP1,NPP1))
   10     CONTINUE
   11   CONTINUE
      STABL = .TRUE.
      GO TO 7
C
C Solve the 2 by 2 (or 3 by 3 if K is a constraint node)
C   lower triangular system.
C
   12 ZK = 0.
      IF (KK .GE. LCC1) ZK = A(7,6)/A(6,6)
      DY = (A(7,5) - A(6,5)*ZK)/A(5,5)
      DX = SF*(A(7,4) - A(6,4)*ZK - A(5,4)*DY)/A(4,4)
      DY = SF*DY
      IF (KK .GE. LCC1) Z(K) = ZK
      NDV = ND
      IF (L .LT. N) NPTS(L+1) = 0
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   13 NDV = 0
      IER = 1
      RETURN
C
C No unique solution due to collinear non-constraint nodes.
C
   14 NDV = 0
      IER = 2
      RETURN
      END
      SUBROUTINE ZINIT (NCC,LCC,N,X,Y,LIST,LPTR,
     .                  LEND, Z, IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        IER
      REAL    X(N), Y(N), Z(N)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   08/27/91
C
C   Given a triangulation of N nodes, along with data values
C at non-constraint nodes, this subroutine computes approxi-
C mate data values at the constraint nodes.  The approximate
C values are intended only to serve as initial estimates for
C Subroutine ZGRADG which computes refined estimates.
C
C   For each subsequence (KM2,KM1,K) of a constraint, the
C approximate value at node KM1 is taken to be the closest-
C point value (data value at the closest non-constraint
C node) at KM1 averaged with the value at KM1 of the linear
C interpolant (along the constraint boundary) of the approx-
C imate value at KM2 and the closest-point value at K.
C
C On input:
C
C       NCC = Number of constraint curves (refer to TRIPACK
C             Subroutine ADDCST).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index of the
C             first node of constraint I in LCC(I).  For I =
C             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
C             LCC(NCC+1) = N+1, and LCC(1) .GE. 4.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRIPACK
C                        Subroutine TRMESH.
C
C The above parameters are not altered by this routine.
C
C       Z = Array of length N containing data values in the
C           first LCC(1)-1 locations.
C
C On output:
C
C       Z = Array updated with approximate data values in
C           the last N-LCC(1)+1 locations if IER = 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if NCC, N, or an LCC entry is outside
C                     its valid range on input.
C
C TRIPACK modules required by ZINIT:  GETNP, INTSEC
C
C Intrinsic functions called by ZINIT:  ABS, SQRT
C
C***********************************************************
C
      INTEGER   LMAX
      PARAMETER (LMAX=12)
      INTEGER   I, IERR, IFRST, ILAST, ILSTM1, K, KM1, KM2,
     .          KN, LCC1, LNP, LP, LPL, NPTS(LMAX)
      REAL      D, DMIN, DS(LMAX), H1, H2, XK, YK, ZN
C
C Test for errors in input parameters.  (LCC is tested by
C   Subroutine GETNP.)
C
      IER = 1
      IF (NCC .GT. 0) THEN
        LCC1 = LCC(1)
      ELSE
        LCC1 = N+1
      ENDIF
      IF (NCC .LT. 0  .OR.  LCC1 .LT. 4) RETURN
C
C Outer loop on constraint I with first and last nodes IFRST
C   and ILAST.
C
      DO 6 I = 1,NCC
        IFRST = LCC(I)
        IF (I .LT. NCC) THEN
          ILAST = LCC(I+1) - 1
        ELSE
          ILAST = N
        ENDIF
C
C Initialize Z(ILAST) with the data value at the closest
C   non-constraint node to ILAST.  Unless the LMAX closest
C   nodes to ILAST (including ILAST) are all constraint
C   nodes, NPTS is set to the closest LNP nodes (with
C   distance measured in the non-constraint region), where
C   LNP is the smallest integer such that NPTS contains a
C   non-constraint node.  The value at LCC(1)-1 is used if
C   LMAX is too small.
C
        LNP = 1
        NPTS(1) = ILAST
        DS(1) = 0.
    1   LNP = LNP + 1
          CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,
     .                LNP, NPTS,DS, IERR)
          IF (IERR .NE. 0) RETURN
          KN = NPTS(LNP)
          IF (KN .GE. LCC1  .AND.  LNP .LT. LMAX) GO TO 1
        IF (KN .GE. LCC1) KN = LCC1-1
        Z(ILAST) = Z(KN)
C
C Loop on constraint nodes K.  LPL points to the last
C   neighbor of K.  At each step, Z(K) is set to the
C   closest-point value at K, and Z(KM1) is set to the
C   (final) approximate data value at KM1 (except when
C   K = IFRST).
C
        KM1 = ILAST
        ILSTM1 = ILAST - 1
        DO 5 K = IFRST,ILSTM1
          XK = X(K)
          YK = Y(K)
          LPL = LEND(K)
C
C   Set LP to point to KM1 as a neighbor of K.
C
          LP = LPL
    2     LP = LPTR(LP)
            IF (ABS(LIST(LP)) .NE. KM1) GO TO 2
C
C   Initialize for loop on non-constraint node neighbors of
C     K.  If K has no such neighbors, the closest non-
C     constraint node to K is (implicitly) taken to be the
C     closest non-constraint node to KM1.
C
          DMIN = -1.
          ZN = Z(KM1)
    3     LP = LPTR(LP)
            KN = ABS(LIST(LP))
            IF (KN .EQ. K+1) GO TO 4
            IF (KN .GE. LCC1) GO TO 3
            D = (X(KN)-XK)**2 + (Y(KN)-YK)**2
            IF (DMIN .GE. 0.  .AND.  DMIN .LT. D) GO TO 3
            DMIN = D
            ZN = Z(KN)
            GO TO 3
C
C   ZN is the closest-point value at K.  Set H2 to the arc
C     length of KM1-K, and compute Z(KM1) if K > IFRST.
C     (H1 is the arc length of KM2-KM1).
C
    4     H2 = SQRT( (XK-X(KM1))**2 + (YK-Y(KM1))**2 )
          IF (K .NE. IFRST) Z(KM1) = .5*( Z(KM1) +
     .      (H1*ZN+H2*Z(KM2))/(H1+H2) )
          Z(K) = ZN
C
C   Bottom of loop on K.
C
          H1 = H2
          KM2 = KM1
          KM1 = K
    5     CONTINUE
C
C For K = ILAST, the closest-point value has already been
C   computed.
C
        H2 = SQRT ( (X(ILAST)-X(ILSTM1))**2 +
     .              (Y(ILAST)-Y(ILSTM1))**2 )
        Z(ILSTM1) = .5*( Z(ILSTM1) + (H1*Z(ILAST)+H2*Z(KM2))
     .                              /(H1+H2) )
C
C Compute the final value at ILAST.
C
        H1 = H2
        H2 = SQRT ( (X(IFRST)-X(ILAST))**2 +
     .              (Y(IFRST)-Y(ILAST))**2 )
        Z(ILAST) = .5*(Z(ILAST) + (H1*Z(IFRST)+H2*Z(ILSTM1))
     .                           /(H1+H2))
    6   CONTINUE
C
C No errors encountered.
C
      IER = 0
      RETURN
      END
