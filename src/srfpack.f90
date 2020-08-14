      include 'tripack.f90'
      SUBROUTINE ARCINT (B,X1,X2,Y1,Y2,H1,H2,HX1,HX2,HY1,&
     &                   HY2,SIGMA,DFLAG, HP,HXP,HYP,IER)
      INTEGER IER
      LOGICAL DFLAG
      REAL    B, X1, X2, Y1, Y2, H1, H2, HX1, HX2, HY1,&
     &        HY2, SIGMA, HP, HXP, HYP
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/18/96
!
!   Gigen a line segment P1-P2 containing a point P =
! (XP,YP), along with function values and partial deriva-
! tives at the endpoints, this subroutine computes an
! interpolated value and, optionally, a gradient at P.  The
! value and tangential gradient component at P are taken to
! be the value and derivative of the Hermite interpolatory
! tension spline H defined by the endpoint values and tan-
! gential gradient components.  The normal gradient compo-
! nent at P is obtained by linear interpolation applied to
! the normal components at the endpoints.
!
! On input:
!
!       B = Local coordinate of P with respect to P1-P2:
!           P = B*P1 + (1-B)*P2.  Note that B may be comput-
!           ed from the coordinates of P as <P2-P1,P2-P>/
!           <P2-P1,P2-P1>.
!
!       X1,X2,Y1,Y2 = Coordinates of a pair of distinct
!                     points P1 and P2.
!
!       H1,H2 = Values of the interpolant H at P1 and P2,
!               respectively.
!
!       HX1,HX2,HY1,HY2 = x and y partial derivatives of H
!                         at P1 and P2.
!
!       SIGMA = Tension factor associated with P1-P2.
!
!       DFLAG = Logical flag which specifies whether first
!               partial derivatives at P are to be computed:
!               DFLAG = .TRUE. if and only if partials are
!               to be returned.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       HP = Interpolated value at P unless IER < 0, in
!            which case HP is not defined.
!
!       HXP,HYP = x and y partial derivatives at P unless
!                 DFLAG = FALSE or IER < 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if B < 0 or B > 1 and thus HP is an
!                     extrapolated value.
!             IER = -1 if P1 and P2 coincide.
!
! SRFPACK module required by ARCINT:  SNHCSH
!
! Intrinsic functions called by ARCINT:  ABS, EXP
!
!***********************************************************
!
      REAL B1, B2, CM, CM2, CMM, D1, D2, DS, DUMMY, DX, DY,&
     &     E, E1, E2, EMS, GN, GT, S, S1, S2, SB1, SB2,&
     &     SBIG, SIG, SINH2, SM, SM2, TM, TM1, TM2, TP1,&
     &     TP2, TS
      DATA SBIG/85./
!
      DX = X2 - X1
      DY = Y2 - Y1
      DS = DX*DX + DY*DY
      IF (DS .EQ. 0.) GO TO 1
      IER = 0
!
! Compute local coordinates B1 and B2, tangential deriva-
!   tives S1 and S2, slope S, and second differences D1 and
!   D2.  S1, S2, S, D1, and D2 are scaled by the separation
!   D between P1 and P2.
!
      B1 = B
      B2 = 1. - B1
      IF (B1 .LT. 0.  .OR.  B2 .LT. 0.) IER = 1
      S1 = HX1*DX + HY1*DY
      S2 = HX2*DX + HY2*DY
      S = H2 - H1
      D1 = S - S1
      D2 = S2 - S
!
! Compute HP and, if required, the scaled tangential grad-
!   ient component GT.
!
      SIG = ABS(SIGMA)
      IF (SIG .LT. 1.E-9) THEN
!
! SIG = 0:  use Hermite cubic interpolation.
!
        HP = H1 + B2*(S1 + B2*(D1 + B1*(D1 - D2)))
        IF (.NOT. DFLAG) RETURN
        GT = S1 + B2*(D1 + D2 + 3.*B1*(D1 - D2))
      ELSEIF (SIG .LE. .5) THEN
!
! 0 .LT. SIG .LE. .5:  use approximations designed to avoid
!   cancellation error in the hyperbolic functions.
!
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        E = SIG*SM - CMM - CMM
        HP = H1 + B2*S1 + ((CM*SM2-SM*CM2)*(D1+D2) + SIG*&
     &                     (CM*CM2-(SM+SIG)*SM2)*D1)/(SIG*E)
        IF (.NOT. DFLAG) RETURN
        SINH2 = SM2 + SB2
        GT = S1 + ((CM*CM2-SM*SINH2)*(D1+D2) + SIG*&
     &             (CM*SINH2-(SM+SIG)*CM2)*D1)/E
      ELSE
!
! SIG > .5:  use negative exponentials in order to avoid
!   overflow.  Note that EMS = EXP(-SIG).  In the case of
!   extrapolation (negative B1 or B2), H is approximated
!   by a linear function if -SIG*B1 or -SIG*B2 is large.
!
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
          HP = H1 + B2*S + (TM*TM1*TM2*(D1+D2) + SIG*&
     &                      ((E2*TM1*TM1-B1*TS)*D1 +&
     &                       (E1*TM2*TM2-B2*TS)*D2))/(SIG*E)
          IF (.NOT. DFLAG) RETURN
          TP1 = 1. + E1
          TP2 = 1. + E2
          GT = S + (TM1*(TM*TP2-SIG*E2*TP1)*D1 -&
     &              TM2*(TM*TP1-SIG*E1*TP2)*D2)/E
        ENDIF
      ENDIF
!
! Compute the gradient at P, (HXP,HYP) = (GT/D)T + (GN/D)N,
!   where T = (DX,DY)/D (unit tangent vector), N = (-DY,DX)/
!   D (unit normal), and the scaled normal component is GN =
!   B1<(HX1,HY1),N> + B2<(HX2,HY2),N>.
!
      GN = B1*(HY1*DX-HX1*DY) + B2*(HY2*DX-HX2*DY)
      HXP = (GT*DX - GN*DY)/DS
      HYP = (GT*DY + GN*DX)/DS
      RETURN
!
! P1 and P2 coincide.
!
    1 IER = -1
      RETURN
      END
      SUBROUTINE CNTOUR (NX,NY,X,Y,Z,CVAL,LC,NCMAX,IWK, XC,&
     &                   YC,ILC,NC,IER)
      INTEGER NX, NY, LC, NCMAX, IWK(NX,*), ILC(NCMAX), NC,&
     &        IER
      REAL    X(NX), Y(NY), Z(NX,NY), CVAL, XC(LC), YC(LC)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   04/28/90
!
!   Given a set of function values Z = F(X,Y) at the verti-
! ces of an NX by NY rectangular grid, this subroutine de-
! termines a set of contour lines associated with F = CVAL.
! A contour line is specified by an ordered sequence of
! points (XC,YC), each lying on a grid edge and computed
! from the linear interpolant of the function values at the
! endpoints of the edge.  The accuracy of the contour lines
! is thus directly related to the number of grid points.  If
! a contour line forms a closed curve, the first point coin-
! cides with the last point.  Otherwise, the first and last
! points lie on the grid boundary.
!
!   Note that the problem is ill-conditioned in the vicinity
! of a double zero of F-CVAL.  Thus, if a grid cell is
! crossed by two contour lines (all four sides intersected),
! three different configurations are possible, corresponding
! to a local minimum, a local maximum, or a saddle point.
! It is arbitrarily assumed in this case that the contour
! lines intersect, representing a saddle point.  Also, in
! order to treat the case of F = CVAL at a vertex in a con-
! sistent manner, this case is always treated as F > CVAL.
! Hence, if F takes on the same value at both ends of an
! edge, it is assumed that no contour line intersects that
! edge.  In particular, a constant function, including
! F = CVAL, results in no contour lines.
!
! On input:
!
!       NX = Number of grid points in the x direction.
!            NX .GE. 2.
!
!       NY = Number of grid points in the y direction.
!            NY .GE. 2.
!
!       X = Array of length NX containing a strictly in-
!           creasing sequence of values.
!
!       Y = Array of length NY containing a strictly in-
!           creasing sequence of values.
!
!       Z = Array of function values at the vertices of the
!           rectangular grid.  Z(I,J) = F(X(I),Y(J)) for
!           I = 1,...,NX and J = 1,...,NY.
!
!       CVAL = Constant function value defining a contour
!              line as the set of points (X,Y) such that
!              F(X,Y) = CVAL.
!
!       LC = Length of arrays XC and YC, and maximum allow-
!            able number of points defining contour lines.
!            LC = 2(NX-1)(NY-1) + (NX*NY+1)/2 is (probably
!            more than) sufficient.  LC .GE. 2.
!
!       NCMAX = Length of array ILC, and maximum allowable
!               number of contour lines.  NCMAX = (NX*NY+1)/
!               2 is sufficient.  NCMAX .GE. 1.
!
! The above parameters are not altered by this routine.
!
!       IWK = Integer array of length .GE. NX*(NY-1) to be
!             used as work space.
!
!       XC,YC = Arrays of length LC.
!
!       ILC = Integer array of length NCMAX.
!
! On output:
!
!       XC,YC = Arrays containing the coordinates of NC con-
!               tour lines.  For K = 1,...,NC, contour line
!               K is defined by the sequence of points with
!               indexes ILC(K-1)+1,...,ILC(K) where ILC(0) =
!               0.
!
!       ILC = Array containing the indexes (to XC and YC)
!             associated with the terminal point of contour
!             line K in position K for K = 1,...,NC (if NC
!             .GT. 0).
!
!       NC = Number of contour lines whose points are stored
!            in XC and YC.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and all
!                     contour lines were found.
!             IER = 1 if NX, NY, LC, or NCMAX is outside its
!                     valid range.  NC = 0 and XC, YC, and
!                     ILC are not altered in this case.
!             IER = 2 if X or Y is not strictly increasing.
!                     NC = 0 and XC, YC, and ILC are not
!                     altered in this case.
!             IER = K for K > LC, where K is the required
!                     length of XC and YC, if more storage
!                     space is required to complete the
!                     specification of contour line NC and/
!                     or additional contour lines up to a
!                     total of NCMAX.  NC .GE. 1 and ILC(NC)
!                     = LC in this case.
!            IER = -1 if more than NCMAX contour lines are
!                     present (more space is required in
!                     ILC).  NC = NCMAX, and LC may or may
!                     not be sufficient for the additional
!                     contour lines in this case.  (This is
!                     not determined.)
!
!   In the unlikely event of an internal failure, a message
! is printed on logical unit LUN (specified in the DATA
! statement below).  IER may be 0 in this case.
!
! Modules required by CNTOUR:  None
!
!***********************************************************
!
      INTEGER I, I1, I2, IB, IN, IND, ISID, ISIDB, ISIDN,&
     &        J, J1, J2, JB, JN, K, LCON, LMX, LUN, NCMX,&
     &        NCON, NI, NIM1, NJ, NJM1
      LOGICAL BDRY
      REAL    CV, W, XF, XN, XP, YF, YN, YP, Z1, Z2
      DATA    LUN/0/
!
! Store parameters in local variables.
!
      NI = NX
      NJ = NY
      NIM1 = NI - 1
      NJM1 = NJ - 1
      CV = CVAL
      LMX = LC
      NCMX = NCMAX
      NC = 0
!
! Test for invalid input parameters.
!
      IER = 1
      IF (NI .LT. 2  .OR.  NJ .LT. 2  .OR.  LMX .LT. 2  .OR.&
     &    NCMX .LT. 1) RETURN
!
! Test for nonincreasing values of X or Y.
!
      IER = 2
      DO 1 I = 2,NI
        IF (X(I) .LE. X(I-1)) RETURN
    1   CONTINUE
      DO 2 J = 2,NJ
        IF (Y(J) .LE. Y(J-1)) RETURN
    2   CONTINUE
!
! Loop on grid cells, initializing edge indicators (stored
!   in IWK) to zeros.  For each cell, the indicator IND is a
!   4-bit integer with each bit corresponding to an edge of
!   the cell, and having value 1 iff the edge has been pro-
!   cessed.  Note that two IND values must be adjusted when
!   an interior edge is processed.  The cell sides (edges)
!   are numbered (1,2,4,8) in counterclockwise order start-
!   ing from the bottom.  This corresponds to an ordering of
!   the weighted IND bits from low order to high order.
!   Grid cells are identified with their lower left corners.
!
      DO 4 J = 1,NJM1
        DO 3 I = 1,NIM1
          IWK(I,J) = 0
    3     CONTINUE
    4   CONTINUE
!
! First determine open contours by looping on boundary edges
!   in counterclockwise order starting from the lower left.
!   For each unprocessed boundary edge intersected by a con-
!   tour line, the contour line is determined and IWK is up-
!   dated to reflect the edges intersected.  The boundary
!   cell (lower left corner) is indexed by (IB,JB) and the
!   boundary edge is specified by ISIDB.  NCON and LCON are
!   local variables containing the number of contour lines
!   encountered and the current length of XC and YC.
!
      NCON = 0
      LCON = 0
      ISIDB = 1
      IB = 1
      JB = 1
!
! Top of loop on boundary edges.  The edge has been
!   processed iff IND/ISIDB is odd.
!
    5 IND = IWK(IB,JB)
      IF (IND/ISIDB .NE. 2*((IND/ISIDB)/2)) GO TO 9
!
! Update the edge indicator and store the vertex indexes of
!   the endpoints of the edge.
!
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
!
! Proceed to the next edge if there is no intersection.
!
      Z1 = Z(I1,J1)
      Z2 = Z(I2,J2)
      IF ((Z1 .LT. CV  .AND.  Z2 .LT. CV)  .OR.&
     &    (Z1 .GE. CV  .AND.  Z2 .GE. CV)) GO TO 9
!
! Store the zero of the linear interpolant of Z1-CV and
!   Z2-CV as the first point of an open contour unless
!   NCMAX contour lines have been found or there is in-
!   sufficient space reserved for XC and YC.
!
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
!
! Initialize for loop on cells intersected by the open
!   contour line.
!
      I = IB
      J = JB
      ISID = ISIDB
!
! Traverse the contour line.  Cell (I,J) was entered on side
!   ISID = (I1,J1)->(I2,J2).  Find an exit edge E (unproces-
!   sed edge intersected by the contour) by looping on the
!   remaining three sides, starting with the side opposite
!   ISID.
!
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
!
! Test for a 1 in bit position ISID of cell (I,J) and bypass
!   the edge if it has been previously encountered.
!
        IF (IND/ISID .NE. 2*((IND/ISID)/2)) GO TO 7
!
! Update IWK for edge E = (I1,J1)->(I2,J2).  (IN,JN) indexes
!   the cell which shares E with cell (I,J), and ISIDN is
!   the side number of E in (IN,JN).  BDRY is true iff E is
!   a boundary edge (with no neighboring cell).
!
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
        BDRY = IN .EQ. 0  .OR.  IN .EQ. NI  .OR.&
     &         JN .EQ. 0  .OR.  JN .EQ. NJ
        IF (.NOT. BDRY) IWK(IN,JN) = IWK(IN,JN) + ISIDN
!
! Exit the loop on sides if E is intersected by the contour.
!
        Z1 = Z(I1,J1)
        Z2 = Z(I2,J2)
        IF ((Z1 .LT. CV  .AND.  Z2 .GE. CV)  .OR.&
     &      (Z1 .GE. CV  .AND.  Z2 .LT. CV)) GO TO 8
    7   CONTINUE
!*
! Error -- No exit point found.  Print a message and exit
!          the contour traversal loop.
!
      WRITE (LUN,100) NCON
  100 FORMAT (///5X,'Error in CNTOUR:  Contour line L ',&
     &        'begins on the boundary'/5X,'and terminates ',&
     &        'in the interior for L =',I4/)
      ILC(NCON) = LCON
      GO TO 9
!*
! Add the intersection point (XN,YN) to the list unless it
!   coincides with the previous point (XP,YP) or there is
!   not enough space in XC and YC.
!
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
!
! Bottom of contour traversal loop.  If E is not a boundary
!   edge, reverse the edge direction (endpoint indexes) and
!   update the cell index and side number.
!
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
!
! Update ILC with a pointer to the end of the contour line.
!
      ILC(NCON) = LCON
!
! Bottom of loop on boundary edges.  Update the boundary
!   cell index and side number, and test for termination.
!
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
!
! Determine closed contours by looping on interior edges --
!   the first two sides (bottom and right) of each cell,
!   excluding boundary edges.  The beginning cell is indexed
!   by (IB,JB), and the beginning side number is ISIDB.
!
      DO 15 JB = 1,NJM1
      DO 14 IB = 1,NIM1
      DO 13 ISIDB = 1,2
        IF (JB .EQ. 1  .AND.  ISIDB .EQ. 1) GO TO 13
        IF (IB .EQ. NIM1  .AND.  ISIDB .EQ. 2) GO TO 13
!
! Bypass the edge if it was previously encountered
!   (IND/ISIDB odd).
!
        IND = IWK(IB,JB)
        IF (IND/ISIDB .NE. 2*((IND/ISIDB)/2)) GO TO 13
!
! Determine the endpoint indexes of the beginning edge E =
!   (I1,J1)->(I2,J2), find the index (I,J) and side number
!   ISID of the cell which shares E with (IB,JB), and up-
!   date IWK.
!
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
!
! Proceed to the next interior edge if there is no
!   intersection.
!
        Z1 = Z(I1,J1)
        Z2 = Z(I2,J2)
        IF ((Z1 .LT. CV  .AND.  Z2 .LT. CV)  .OR.&
     &      (Z1 .GE. CV  .AND.  Z2 .GE. CV)) GO TO 13
!
! Store the intersection point as the first point of a
!   closed contour unless NCMAX contour lines have been
!   found or there is insufficient space in XC and YC.
!
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
!
! Traverse the contour line.  Cell (I,J) was entered on side
!   ISID = edge (I2,J2)->(I1,J1).  Reverse the edge direc-
!   tion.
!
   10   IN = I1
        JN = J1
        I1 = I2
        J1 = J2
        I2 = IN
        J2 = JN
        IND = IWK(I,J)
!
! Find an exit edge E by looping on the remaining three
!   sides, starting with the side opposite ISID.
!
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
!
! Bypass the edge if it has been previously encountered.
!
          IF (IND/ISID .NE. 2*((IND/ISID)/2)) GO TO 11
!
! Determine the index (IN,JN) and side number ISIDN of the
!   cell which shares edge E = (I1,J1)->(I2,J2) with cell
!   (I,J), and update IWK.
!
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
!
! Exit the loop on sides if E is intersected.
!
          Z1 = Z(I1,J1)
          Z2 = Z(I2,J2)
          IF ((Z1 .LT. CV  .AND.  Z2 .GE. CV)  .OR.&
     &        (Z1 .GE. CV  .AND.  Z2 .LT. CV)) GO TO 12
   11     CONTINUE
!*
! Error -- No exit point found.  Print a message and exit
!          the contour traversal loop.
!
        WRITE (LUN,110) NCON
  110   FORMAT (///5X,'Error in CNTOUR:  Contour line L ',&
     &          'is open but'/5X,'does not intersect the ',&
     &          'boundary for L =',I4/)
        ILC(NCON) = LCON
        GO TO 13
!*
! Add the intersection point to the list unless it coincides
!   with the previous point or there is not enough space in
!   XC and YC.
!
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
!
! Bottom of contour traversal loop.  If the next cell is not
!   the beginning cell, update the cell index and side num-
!   ber.
!
        IF (IN .NE. IB  .OR.  JN .NE. JB) THEN
          I = IN
          J = JN
          ISID = ISIDN
          GO TO 10
        ENDIF
!
! Add the first point as the last point (unless the first
!   and last points already coincide), and update ILC.
!
        IF (XP .NE. XF  .OR.  YP .NE. YF) THEN
          LCON = LCON + 1
          IF (LCON .LE. LMX) THEN
            XC(LCON) = XF
            YC(LCON) = YF
          ENDIF
        ENDIF
        ILC(NCON) = LCON
!
! Bottom of loop on interior edges.
!
   13   CONTINUE
   14   CONTINUE
   15   CONTINUE
      IER = 0
!
! Test for insufficient storage reserved for XC and YC.
!
   16 IF (LCON .GT. LMX) IER = LCON
      NC = NCON
      RETURN
      END
      SUBROUTINE COORDS (XP,YP,X1,X2,X3,Y1,Y2,Y3, B1,B2,&
     &                   B3, IER)
      INTEGER IER
      REAL    XP, YP, X1, X2, X3, Y1, Y2, Y3, B1, B2, B3
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88
!
!   This subroutine computes the barycentric (areal) coordi-
! nates B1, B2, and B3 of a point P with respect to the tri-
! angle with vertices P1, P2, and P3:  the solution to the
! linear system defined by B1 + B2 + B3 = 1 and B1*P1 +
! B2*P2 + B3*P3 = P.  Note that B1 is a linear function
! of P which satisfies B1 = 1 at P = P1 and B1 = 0 on
! the triangle side P2-P3.  Also, B1 < 0 if and only if
! P is to the right of P2->P3 (and thus exterior to the
! triangle).  B2 and B3 satisfy similar properties.
!
! On input:
!
!       XP,YP = Cartesian coordinates of P.
!
!       X1,X2,X3,Y1,Y2,Y3 = Coordinates of the vertices of
!                           the triangle P1, P2, and P3.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       B1,B2,B3 = Barycentric coordinates unless IER = 1,
!                  in which case the coordinates are not
!                  defined.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if the vertices of the triangle are
!                     collinear.
!
! Modules required by COORDS:  None
!
!***********************************************************
!
      REAL A, PX, PY, XP1, XP2, XP3, YP1, YP2, YP3
!
      PX = XP
      PY = YP
!
! Compute components of the vectors P->P1, P->P2, and P->P3.
!
      XP1 = X1 - PX
      YP1 = Y1 - PY
      XP2 = X2 - PX
      YP2 = Y2 - PY
      XP3 = X3 - PX
      YP3 = Y3 - PY
!
! Compute subtriangle areas B1 = P->P2 X P->P3, B2 = P->P3 X
!   P->P1, and B3 = P->P1 X P->P2.
!
      B1 = XP2*YP3 - XP3*YP2
      B2 = XP3*YP1 - XP1*YP3
      B3 = XP1*YP2 - XP2*YP1
!
! Compute twice the signed area of the triangle.
!
      A = B1 + B2 + B3
      IF (A .EQ. 0.) GO TO 1
!
! Normalize the coordinates.
!
      B1 = B1/A
      B2 = B2/A
      B3 = B3/A
      IER = 0
      RETURN
!
! The vertices are collinear.
!
    1 IER = -1
      RETURN
      END
      SUBROUTINE CRPLOT (LUN,PLTSIZ,NX,NY,PX,PY,PZ,NCON,IWK,&
     &                   XC,YC, IER)
      INTEGER LUN, NX, NY, NCON, IWK(*), IER
      REAL    PLTSIZ, PX(NX), PY(NY), PZ(NX,NY),&
     &        XC(*), YC(*)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   04/12/97
!
!   Given a set of function values PZ = F(X,Y) at the ver-
! tices of an NX by NY rectangular grid, this subroutine
! creates a level-2 Encapsulated PostScript (EPS) file
! containing a contour plot of the piecewise bilinear inter-
! polant of the function values.
!
!   The accuracy of the contour lines increases with the
! number of grid points.  Refer to Subroutine CNTOUR for
! further details.
!
!
! On input:
!
!       LUN = Logical unit number in the range 0 to 99.
!             The unit should be opened with an appropriate
!             file name before the call to this routine.
!
!       PLTSIZ = Plot size in inches.  A window containing
!                the plot is mapped, with aspect ratio
!                preserved, to a rectangular viewport with
!                maximum side-length PLTSIZ.  The viewport
!                is centered on the 8.5 by 11 inch page, and
!                its boundary is drawn.  1.0 .LE. PLTSIZ
!                .LE. 7.5.
!
!       NX = Number of grid points in the x direction.
!            NX .GE. 2.
!
!       NY = Number of grid points in the y direction.
!            NY .GE. 2.
!
!       PX = Array of length NX containing a strictly in-
!            creasing sequence of values.
!
!       PY = Array of length NY containing a strictly in-
!            creasing sequence of values.
!
!       PZ = Array of function values at the vertices of the
!            rectangular grid.  PZ(I,J) = F(PX(I),PY(J)) for
!            I = 1,...,NX and J = 1,...,NY.
!
!       NCON = Number of contour values.  The contour values
!              are uniformly distributed over the range of
!              PZ values.  NCON .GE. 1.
!
! The above parameters are not altered by this routine.
!
!       IWK = Integer array of length at least 1.5*NX*NY to
!             be used as work space.
!
!       XC,YC = Real arrays of length at least 2.5*NX*NY to
!               be used as work space.
!
! On output:
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if LUN, PLTSIZ, NX, NY, or NCON is
!                     outside its valid range.
!             IER = 2 if PX or PY is not strictly
!                     increasing.
!             IER = 3 if the range of PZ values has zero
!                     width (F is constant).
!             IER = 4 if an error was encountered in writing
!                     to unit LUN.
!             IER = 5 if an unexpected error flag was re-
!                     turned by Subroutine CNTOUR.  This
!                     should not occur.
!
!   In the unlikely event of an internal failure, a message
! is printed on the standard output device.  IER may be 0
! in this case.
!
! Module required by CRPLOT:  CNTOUR
!
! Intrinsic functions called by CRPLOT:  CHAR, REAL
!
!***********************************************************
!
      INTEGER I, IC, IERR, IH, IPX1, IPX2, IPY1, IPY2, IW,&
     &        J, K, KV, LC, NC, NCMAX
      REAL    CVAL, DX, DY, DZ, PZIJ, R, SFX, SFY, T,&
     &        TX, TY, ZMAX, ZMIN
!
! Local parameters:
!
! CVAL =      Contour value between ZMIN and ZMAX
! DX =        Window width PX(NX)-PX(1)N
! DY =        Window height PY(NY)-PY(1)
! DZ =        Interval between contour values:
!               (ZMAX-ZMIN)/(NCON+1)
! I,J =       Row and column indexes for PZ
! IC =        Index (for IWK) of a contour line associated
!               with contour value CVAL:  1 to NC
! IERR =      Error flag for calls to CNTOUR
! IH =        Height of the bounding box (viewport) in
!               points
! IPX1,IPY1 = X and y coordinates (in points) of the lower
!               left corner of the bounding box
! IPX2,IPY2 = X and y coordinates (in points) of the upper
!               right corner of the bounding box
! IW =        Width of the bounding box in points
! K =         Index (for XC and YC) of a point on a contour
!               line
! KV =        DO-loop index for loop on contour values
! LC =        Length of arrays XC and YC
! NC =        Number of contour lines associated with
!               contour value CVAL
! NCMAX =     Maximum allowable value of NC
! PZIJ =      PZ(I,J)
! R =         Aspect ratio DX/DY
! SFX,SFY =   Scale factors for mapping window coordinates
!               to viewport coordinates
! T =         Temporary variable
! TX,TY =     Translation vector for mapping window coordi-
!               nates to viewport coordinates
! ZMIN,ZMAX = Minimum and maximum of the PZ values
!
!
! Test for error 1.
!
      IF (LUN .LT. 0  .OR.  LUN .GT. 99  .OR.&
     &    PLTSIZ .LT. 1.0  .OR.  PLTSIZ .GT. 7.5  .OR.&
     &    NX .LT. 2  .OR.  NY .LT. 2  .OR.  NCON .LT. 1)&
     &  GO TO 11
!
! Compute the aspect ratio of the window.
!
      DX = PX(NX) - PX(1)
      DY = PY(NY) - PY(1)
      IF (DX .EQ. 0.0  .OR.  DY .EQ. 0.0) GO TO 12
      R = DX/DY
!
! Compute the range of PZ values and the interval between
!   contour values.
!
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
!
! Compute the lower left (IPX1,IPY1) and upper right
!   (IPX2,IPY2) corner coordinates of the bounding box
!   (the viewport).  The coordinates, specified in default
!   user space units (points, at 72 points/inch with origin
!   at the lower left corner of the page), are chosen to
!   preserve the aspect ratio R, and to center the plot on
!   the 8.5 by 11 inch page.  The center of the page is
!   (306,396), and T = PLTSIZ/2 in points.
!
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
!
! Output header comments.
!
      WRITE (LUN,100,ERR=14) IPX1, IPY1, IPX2, IPY2
  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/&
     &        '%%BoundingBox:',4I4/&
     &        '%%Title:  Contour Plot'/&
     &        '%%Creator:  SRFPACK'/&
     &        '%%EndComments')
!
! Draw the bounding box.
!
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
!
! Set up a mapping from the window to the viewport.
!
      IW = IPX2 - IPX1
      IH = IPY2 - IPY1
      SFX = REAL(IW)/DX
      SFY = REAL(IH)/DY
      TX = IPX1 - SFX*PX(1)
      TY = IPY1 - SFY*PY(1)
      WRITE (LUN,150,ERR=14) TX, TY, SFX, SFY
  150 FORMAT (2F12.6,' translate'/&
     &        2F12.6,' scale')
!
! Set the line thickness to 2 points.  (Since the scale
!   factors are applied to everything, the width must be
!   specified in world coordinates.)
!
      T = 4.0/(SFX+SFY)
      WRITE (LUN,160,ERR=14) T
  160 FORMAT (F12.6,' setlinewidth')
!
! Compute parameters for CNTOUR:
!
!   NCMAX = Maximum allowable number of contour lines
!           associated with each contour value.
!   LC = Length of arrays XC and YC and maximum allowable
!        number of points defining all the contour lines
!        associated with a contour value.
!
      NCMAX = (NX*NY+1)/2
      LC = 2*(NX-1)*(NY-1) + NCMAX
!
! Loop on contour values CVAL uniformly spaced in the open
!   interval (ZMIN,ZMAX).
!
      CVAL = ZMIN
      DO 5 KV = 1,NCON
        CVAL = CVAL + DZ
!
! Compute a sequence of NC contour lines associated with
!   F = CVAL.  For IC = 1 to NC, IWK(IC) is the index (for
!   XC and YC) of the last point of contour IC.
!
        CALL CNTOUR (NX,NY,PX,PY,PZ,CVAL,LC,NCMAX,&
     &               IWK(NCMAX+1), XC,YC,IWK,NC,IERR)
        IF (IERR .EQ. 2) GO TO 12
        IF (IERR .NE. 0) GO TO 15
!
! Draw the NC contours.
!
        IC = 0
        K = 0
    3   IC = IC + 1
          K = K + 1
!
!   Create a path consisting of contour IC.
!
          WRITE (LUN,170,ERR=14) XC(K), YC(K)
  170     FORMAT (2F12.6,' moveto')
    4     K = K + 1
            WRITE (LUN,180,ERR=14) XC(K), YC(K)
  180       FORMAT (2F12.6,' lineto')
            IF (K .NE. IWK(IC)) GO TO 4
!
!   Paint the path.
!
          WRITE (LUN,140,ERR=14)
          IF (IC .NE. NC) GO TO 3
    5   CONTINUE
!
! Output the showpage command and end-of-file indicator.
!
      WRITE (LUN,200,ERR=14)
  200 FORMAT ('showpage'/&
     &        '%%EOF')
!
! HP's interpreters require a one-byte End-of-PostScript-Job
!   indicator (to eliminate a timeout error message):
!   ASCII 4.
!
      WRITE (LUN,210,ERR=14) CHAR(4)
  210 FORMAT (A1)
!
! No error encountered.
!
      IER = 0
      RETURN
!
! Invalid input parameter.
!
   11 IER = 1
      RETURN
!
! PX or PY is not strictly increasing.
!
   12 IER = 2
      RETURN
!
! DZ = 0.
!
   13 IER = 3
      RETURN
!
! Error writing to unit LUN.
!
   14 IER = 4
      RETURN
!
! Error flag returned by CNTOUR.
!
   15 IER = 5
      RETURN
      END
      SUBROUTINE FVAL (XP,YP,X1,X2,X3,Y1,Y2,Y3,F1,F2,F3,&
     &                 FX1,FX2,FX3,FY1,FY2,FY3,SIG1,SIG2,&
     &                 SIG3, FP,IER)
      INTEGER IER
      REAL    XP, YP, X1, X2, X3, Y1, Y2, Y3, F1, F2,&
     &        F3, FX1, FX2, FX3, FY1, FY2, FY3, SIG1,&
     &        SIG2, SIG3, FP
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   03/18/90
!
!   Given function values and gradients at the three ver-
! tices of a triangle containing a point P, this routine
! computes the value of F at P where F interpolates the ver-
! tex data.  Along the triangle arcs, the interpolatory
! function F is the Hermite interpolatory tension spline de-
! fined by the values and tangential gradient components at
! the endpoints, and the derivative in the direction normal
! to the arc varies linearly between the normal gradient
! components at the endpoints.  A first-order C-1 blending
! method is used to extend F to the interior of the trian-
! gle.  Thus, since values and gradients on an arc depend
! only on the vertex data, the method results in C-1 contin-
! uity when used to interpolate over a triangulation.
!
!   The blending method consists of taking F(P) to be the
! weighted sum of the values at P of the three univariate
! Hermite interpolatory tension splines defined on the line
! segments which join the vertices to the opposite sides and
! pass through P.  The tension factors for these splines are
! obtained by linear interpolation between the pair of ten-
! sion factors associated with the triangle sides which join
! at the appropriate vertex.
!
!   A tension factor SIGMA associated with a Hermite interp-
! olatory tension spline is a nonnegative parameter which
! determines the curviness of the spline.  SIGMA = 0 results
! in a cubic spline, and the spline approaches the linear
! interpolant as SIGMA increases.
!
! On input:
!
!       XP,YP = Coordinates of a point P at which an interp-
!               olated value is to be computed.
!
!       X1,X2,X3,Y1,Y2,Y3 = Coordinates of the vertices of a
!                           triangle (V1,V2,V3) containing
!                           P.  V3 is strictly to the left
!                           of V1->V2.
!
!       F1,F2,F3 = Values of the interpolatory function at
!                  the vertices.
!
!       FX1,FX2,FX3 = x components of the gradients of F at
!                     the vertices.
!
!       FY1,FY2,FY3 = y components of the gradients of F at
!                     the vertices.
!
!       SIG1,SIG2,SIG3 = Tension factors associated with the
!                        arcs opposite V1, V2, and V3, re-
!                        spectively.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       FP = Interpolated value at P unless IER < 0, in
!            which case FP is not defined.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if P is not contained in the triangle.
!                     This may result from roundoff error
!                     when P lies near an arc, and the int-
!                     erpolated value FP is valid in that
!                     case.
!             IER = -1 if the triangle vertices are
!                      collinear.
!
! SRFPACK modules required by FVAL:  ARCINT, COORDS, SNHCSH
!
!***********************************************************
!
      INTEGER IERR
      REAL    B, B1, B2, B3, C1, C2, C3, DUM, FQ, FXQ, FYQ,&
     &        H1, H2, H3, PX, PY, SIG, SUM, XQ, YQ
!
      PX = XP
      PY = YP
!
! F(P) = C1*H1(P) + C2*H2(P) + C3*H3(P) where C1, C2, and C3
!   are weight functions which sum to 1, and H1, H2, and H3
!   are Hermite interpolatory tension splines on the line
!   segments which join vertices to opposite sides and con-
!   tain P.
!
! Compute barycentric coordinates of P with respect to the
!   triangle.
!
      CALL COORDS (PX,PY,X1,X2,X3,Y1,Y2,Y3, B1,B2,B3,IER)
      IF (IER .NE. 0) RETURN
      IF (B1 .LT. 0.  .OR.  B2 .LT. 0.  .OR.  B3 .LT. 0.)&
     &   IER = 1
!
! Compute the coefficients of the partial interpolants.
!   C1 = 1 on the side opposite V1, and C1 = 0 on the other
!   arcs.  Similarly for C2 and C3.
!
      C1 = B2*B3
      C2 = B3*B1
      C3 = B1*B2
      SUM = C1 + C2 + C3
      IF (SUM .EQ. 0.) THEN
!
! P coincides with a vertex.
!
        FP = B1*F1 + B2*F2 + B3*F3
        RETURN
      ENDIF
!
! Normalize the coefficients.
!
      C1 = C1/SUM
      C2 = C2/SUM
      C3 = C3/SUM
!
! For each vertex Vi, compute the intersection Q of the side
!   opposite Vi with the line defined by Vi and P, the value
!   and gradient at Q, and the partial interpolant value Hi
!   at P.
!
!   Side opposite V1:
!
      B = B2/(B2+B3)
      XQ = B*X2 + (1.-B)*X3
      YQ = B*Y2 + (1.-B)*Y3
      SIG = B*SIG3 + (1.-B)*SIG2
      CALL ARCINT (B,X2,X3,Y2,Y3,F2,F3,FX2,FX3,FY2,FY3,SIG1,&
     &             .TRUE., FQ,FXQ,FYQ,IERR)
      CALL ARCINT (B1,X1,XQ,Y1,YQ,F1,FQ,FX1,FXQ,FY1,FYQ,SIG,&
     &             .FALSE., H1,DUM,DUM,IERR)
!
!   Side opposite V2:
!
      B = B3/(B3+B1)
      XQ = B*X3 + (1.-B)*X1
      YQ = B*Y3 + (1.-B)*Y1
      SIG = B*SIG1 + (1.-B)*SIG3
      CALL ARCINT (B,X3,X1,Y3,Y1,F3,F1,FX3,FX1,FY3,FY1,SIG2,&
     &             .TRUE., FQ,FXQ,FYQ,IERR)
      CALL ARCINT (B2,X2,XQ,Y2,YQ,F2,FQ,FX2,FXQ,FY2,FYQ,SIG,&
     &             .FALSE., H2,DUM,DUM,IERR)
!
!   Side opposite V3:
!
      B = B1/(B1+B2)
      XQ = B*X1 + (1.-B)*X2
      YQ = B*Y1 + (1.-B)*Y2
      SIG = B*SIG2 + (1.-B)*SIG1
      CALL ARCINT (B,X1,X2,Y1,Y2,F1,F2,FX1,FX2,FY1,FY2,SIG3,&
     &             .TRUE., FQ,FXQ,FYQ,IERR)
      CALL ARCINT (B3,X3,XQ,Y3,YQ,F3,FQ,FX3,FXQ,FY3,FYQ,SIG,&
     &             .FALSE., H3,DUM,DUM,IERR)
!
! Accumulate the partial interpolant values.
!
      FP = C1*H1 + C2*H2 + C3*H3
      RETURN
      END
      SUBROUTINE GETSIG (N,X,Y,H,LIST,LPTR,LEND,HXHY,&
     &                   TOL, SIGMA, DSMAX,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IER
      REAL    X(N), Y(N), H(N), HXHY(2,N), TOL, SIGMA(*),&
     &        DSMAX
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/03/98
!
!   Given a triangulation of a set of nodes in the plane,
! along with data values H and gradients (HX,HY) at the
! nodes, this subroutine determines, for each triangulation
! arc, the smallest (nonnegative) tension factor SIGMA such
! that the Hermite interpolatory tension spline H(T), de-
! fined by SIGMA and the endpoint values and directional
! derivatives, preserves local shape properties of the data.
! In order to define the shape properties on an arc, it is
! convenient to map the arc to an interval (T1,T2).  Then,
! denoting the endpoint data values by H1,H2 and the deriva-
! tives by HP1,HP2, and letting S = (H2-H1)/(T2-T1), the
! data properties are
!
!       Monotonicity:  S, HP1, and HP2 are nonnegative or
!                        nonpositive,
!   and
!
!       Convexity:     HP1 .LE. S .LE. HP2  or  HP1 .GE. S
!                        .GE. HP2.
!
! The corresponding properties of H are constant sign of the
! first and second derivatives, respectively.  Note that,
! unless HP1 = S = HP2, infinite tension is required (and H
! is linear on the interval) if S = 0 in the case of mono-
! tonicity, or if HP1 = S or HP2 = S in the case of
! convexity.
!
!   Note that if gradients are to be computed by Subroutine
! GRADG or function values and gradients are computed by
! SMSURF, it may be desirable to alternate those computa-
! tions (which require tension factors) with calls to this
! subroutine.  This iterative procedure should terminate
! with a call to GETSIG in order to ensure that the shape
! properties are preserved, and convergence can be achieved
! (at the cost of optimality) by allowing only increases in
! tension factors (refer to the parameter descriptions for
! SIGMA, DSMAX, and IER).
!
!   Refer to functions SIG0, SIG1, and SIG2 for means of
! selecting minimum tension factors to preserve more general
! properties.
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.
!
!       H = Array of length N containing data values at the
!           nodes.  H(I) is associated with (X(I),Y(I)).
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
!       HXHY = Array dimensioned 2 by N whose columns con-
!              tain partial derivatives at the nodes (X
!              partials in the first row).  Refer to Subrou-
!              tines GRADC, GRADG, GRADL, and SMSURF.
!
!       TOL = Tolerance whose magnitude determines how close
!             each tension factor is to its optimal value
!             when nonzero finite tension is necessary and
!             sufficient to satisfy the constraint --
!             abs(TOL) is an upper bound on the magnitude
!             of the smallest (nonnegative) or largest (non-
!             positive) value of the first or second deriva-
!             tive of H in the interval.  Thus, the con-
!             straint is satisfied, but possibly with more
!             tension than necessary.
!
! The above parameters are not altered by this routine.
!
!       SIGMA = Array of length 2*NA = 6*(N-1)-2*NB, where
!               NA and NB are the numbers of arcs and boun-
!               dary nodes, respectively, containing minimum
!               values of the tension factors.  The tension
!               factors are associated with arcs in one-to-
!               one correspondence with LIST entries.  Note
!               that each arc N1-N2 has two LIST entries and
!               thus, the tension factor is stored in both
!               SIGMA(I) and SIGMA(J) where LIST(I) = N2 (in
!               the adjacency list for N1) and LIST(J) = N1
!               (in the list associated with N2).  SIGMA
!               should be set to all zeros if minimal ten-
!               sion is desired, and should be unchanged
!               from a previous call in order to ensure con-
!               vergence of the iterative procedure describ-
!               ed in the header comments.
!
! On output:
!
!       SIGMA = Array containing tension factors for which
!               H(T) preserves the local data properties on
!               each triangulation arc, with the restriction
!               that SIGMA(I) .LE. 85 for all I (unless the
!               input value is larger).  The factors are as
!               small as possible (within the tolerance) but
!               not less than their input values.  If infin-
!               ite tension is required on an arc, the cor-
!               responding factor is SIGMA(I) = 85 (and H
!               is an approximation to the linear inter-
!               polant on the arc), and if neither property
!               is satisfied by the data, then SIGMA(I) = 0
!               (assuming its input value is 0), and thus H
!               is cubic on the arc.
!
!       DSMAX = Maximum increase in a component of SIGMA
!               from its input value.
!
!       IER = Error indicator and information flag:
!             IER = I if no errors were encountered and I
!                     components of SIGMA were altered from
!                     their input values for I .GE. 0.
!             IER = -1 if N < 3.  SIGMA is not altered in
!                      this case.
!             IER = -2 if duplicate nodes were encountered.
!
! TRIPACK modules required by GETSIG:  LSTPTR, STORE
!
! SRFPACK module required by GETSIG:  SNHCSH
!
! Intrinsic functions called by GETSIG:  ABS, EXP, MAX, MIN,
!                                          SIGN, SQRT
!
!***********************************************************
!
      INTEGER LSTPTR
      REAL    STORE
      INTEGER ICNT, LP1, LP2, LPL, LUN, N1, N2, NIT, NM1
      REAL    A, C1, C2, COSHM, COSHMM, D0, D1, D1D2, D1PD2,&
     &        D2, DMAX, DSIG, DSM, DT, DX, DY, E, EMS, EMS2,&
     &        F, F0, FMAX, FNEG, FP, FTOL, RTOL, S, S1, S2,&
     &        SBIG, SCM, SGN, SIG, SIGIN, SINHM, SSINH, SSM,&
     &        STOL, T, T0, T1, T2, TM, TP1
!
      DATA SBIG/85./,  LUN/-1/
      NM1 = N - 1
      IF (NM1 .LT. 2) GO TO 11
!
! Compute an absolute tolerance FTOL = abs(TOL) and a
!   relative tolerance RTOL = 100*Macheps.
!
      FTOL = ABS(TOL)
      RTOL = 1.
    1 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 1
      RTOL = RTOL*200.
!
! Print a heading.
!
      IF (LUN .GE. 0) WRITE (LUN,100) N, FTOL
  100 FORMAT (///13X,'GETSIG:  N =',I4,', TOL = ',E10.3//)
!
! Initialize change counter ICNT and maximum change DSM for
!   the loop on arcs.
!
      ICNT = 0
      DSM = 0.
!
! Loop on arcs N1-N2 for which N2 > N1.  LPL points to the
!   last neighbor of N1.
!
      DO 10 N1 = 1,NM1
        LPL = LEND(N1)
        LP1 = LPL
!
!   Top of loop on neighbors N2 of N1.
!
    2   LP1 = LPTR(LP1)
        N2 = ABS(LIST(LP1))
        IF (N2 .LE. N1) GO TO 9
!
! Print a message and compute parameters for the arc:  DT =
!   arc length and SIGIN = input SIGMA value.
!
        IF (LUN .GE. 0) WRITE (LUN,110) N1, N2
  110   FORMAT (/1X,'Arc',I4,' -',I4)
        DX = X(N2) - X(N1)
        DY = Y(N2) - Y(N1)
        DT = SQRT(DX*DX + DY*DY)
        IF (DT .EQ. 0.) GO TO 12
        SIGIN = SIGMA(LP1)
        IF (SIGIN .GE. SBIG) GO TO 9
!
! Compute scaled directional derivatives S1,S2 at the end-
!   points (for the direction N1->N2), first difference S,
!   and second differences D1,D2.
!
        S1 = HXHY(1,N1)*DX + HXHY(2,N1)*DY
        S2 = HXHY(1,N2)*DX + HXHY(2,N2)*DY
        S = H(N2) - H(N1)
        D1 = S - S1
        D2 = S2 - S
        D1D2 = D1*D2
!
! Test for infinite tension required to satisfy either
!   property.
!
        SIG = SBIG
        IF ((D1D2 .EQ. 0.  .AND.  S1 .NE. S2)  .OR.&
     &      (S .EQ. 0.  .AND.  S1*S2 .GT. 0.)) GO TO 8
!
! Test for SIGMA = 0 sufficient.  The data satisfies convex-
!   ity iff D1D2 .GE. 0, and D1D2 = 0 implies S1 = S = S2.
!
        SIG = 0.
        IF (D1D2 .LT. 0.) GO TO 4
        IF (D1D2 .EQ. 0.) GO TO 8
        T = MAX(D1/D2,D2/D1)
        IF (T .LE. 2.) GO TO 8
        TP1 = T + 1.
!
! Convexity:  find a zero of F(SIG) = SIG*COSHM(SIG)/
!   SINHM(SIG) - TP1.
!
!   F(0) = 2-T < 0, F(TP1) .GE. 0, the derivative of F
!     vanishes at SIG = 0, and the second derivative of F is
!     .2 at SIG = 0.  A quadratic approximation is used to
!     obtain a starting point for the Newton method.
!
        SIG = SQRT(10.*T-20.)
        NIT = 0
!
!   Top of loop:
!
    3   IF (SIG .LE. .5) THEN
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          T1 = COSHM/SINHM
          FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.)
        ELSE
!
!   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
!     overflow with large SIG.
!
          EMS = EXP(-SIG)
          SSM = 1. - EMS*(EMS+SIG+SIG)
          T1 = (1.-EMS)*(1.-EMS)/SSM
          FP = T1 + SIG*(2.*SIG*EMS/SSM - T1*T1 + 1.)
        ENDIF
!
        F = SIG*T1 - TP1
        IF (LUN .GE. 0) WRITE (LUN,120) SIG, F, FP
  120   FORMAT (1X,'Convexity:  SIG = ',E15.8,&
     &          ', F(SIG) = ',E15.8/1X,35X,'FP(SIG) = ',&
     &          E15.8)
        NIT = NIT + 1
!
!   Test for convergence.
!
        IF (FP .LE. 0.) GO TO 8
        DSIG = -F/FP
        IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.  .AND.&
     &      F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 8
!
!   Update SIG.
!
        SIG = SIG + DSIG
        GO TO 3
!
! Convexity cannot be satisfied.  Monotonicity can be satis-
!   fied iff S1*S .GE. 0 and S2*S .GE. 0 since S .NE. 0.
!
    4   IF (S1*S .LT. 0.  .OR.  S2*S .LT. 0.) GO TO 8
        T0 = 3.*S - S1 - S2
        D0 = T0*T0 - S1*S2
!
! SIGMA = 0 is sufficient for monotonicity iff S*T0 .GE. 0
!   or D0 .LE. 0.
!
        IF (D0 .LE. 0.  .OR.  S*T0 .GE. 0.) GO TO 8
!
! Monotonicity:  find a zero of F(SIG) = sign(S)*HP(R),
!   where HPP(R) = 0 and HP, HPP denote derivatives of H.
!   F has a unique zero, F(0) < 0, and F approaches
!   abs(S) as SIG increases.
!
!   Initialize parameters for the secant method.  The method
!     uses three points:  (SG0,F0), (SIG,F), and
!     (SNEG,FNEG), where SG0 and SNEG are defined implicitly
!     by DSIG = SIG - SG0 and DMAX = SIG - SNEG.
!
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
!
!   Top of loop:  compute the change in SIG by linear
!     interpolation.
!
    5   DSIG = -F*DSIG/(F-F0)
        IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130   FORMAT (1X,'Monotonicity:  DSIG = ',E15.8)
        IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.&
     &       DSIG*DMAX .GT. 0. ) GO TO 7
!
!   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
!     Note that DSIG and DMAX have opposite signs.
!
        IF (ABS(DSIG) .LT. STOL/2.) DSIG = -SIGN(STOL/2.,&
     &                              DMAX)
!
!   Update SIG, F0, and F.
!
        SIG = SIG + DSIG
        F0 = F
        IF (SIG .LE. .5) THEN
!
!   Use approximations to the hyperbolic functions designed
!     to avoid cancellation error with small SIG.
!
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          C1 = SIG*COSHM*D2 - SINHM*D1PD2
          C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
          A = C2 - C1
          E = SIG*SINHM - COSHMM - COSHMM
        ELSE
!
!   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
!     overflow with large SIG.
!
          EMS = EXP(-SIG)
          EMS2 = EMS + EMS
          TM = 1. - EMS
          SSINH = TM*(1.+EMS)
          SSM = SSINH - SIG*EMS2
          SCM = TM*TM
          C1 = SIG*SCM*D2 - SSM*D1PD2
          C2 = SIG*SSINH*D2 - SCM*D1PD2
!
!   R is in (0,1) and well-defined iff HPP(T1)*HPP(T2) < 0.
!
          F = FMAX
          IF (C1*(SIG*SCM*D1 - SSM*D1PD2) .GE. 0.) GO TO 6
          A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
          IF (A*(C2+C1) .LT. 0.) GO TO 6
          E = SIG*SSINH - SCM - SCM
        ENDIF
!
        F = (SGN*(E*S2-C2) + SQRT(A*(C2+C1)))/E
!
!   Update the number of iterations NIT.
!
    6   NIT = NIT + 1
        IF (LUN .GE. 0) WRITE (LUN,140) NIT, SIG, F
  140   FORMAT (1X,11X,I2,' -- SIG = ',E15.8,', F = ',&
     &          E15.8)
!
!   Test for convergence.
!
        STOL = RTOL*SIG
        IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.  .AND.&
     &      F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 8
        DMAX = DMAX + DSIG
        IF (F0*F .GT. 0.  .AND.  ABS(F) .GE. ABS(F0))&
     &     GO TO 7
        IF (F0*F .LE. 0.) THEN
!
!   F and F0 have opposite signs.  Update (SNEG,FNEG) to
!     (SG0,F0) so that F and FNEG always have opposite
!     signs.  If SIG is closer to SNEG than SG0 and abs(F)
!     < abs(FNEG), then swap (SNEG,FNEG) with (SG0,F0).
!
          T1 = DMAX
          T2 = FNEG
          DMAX = DSIG
          FNEG = F0
          IF ( ABS(DSIG) .GT. ABS(T1)  .AND.&
     &         ABS(F) .LT. ABS(T2) ) THEN
!
            DSIG = T1
            F0 = T2
          ENDIF
        ENDIF
        GO TO 5
!
!   Bottom of loop:  F0*F > 0 and the new estimate would
!     be outside of the bracketing interval of length
!     abs(DMAX).  Reset (SG0,F0) to (SNEG,FNEG).
!
    7   DSIG = DMAX
        F0 = FNEG
        GO TO 5
!
!  Update SIGMA, ICNT, and DSM if necessary.
!
    8   SIG = MIN(SIG,SBIG)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(LP1) = SIG
          LP2 = LSTPTR(LEND(N2),N1,LIST,LPTR)
          SIGMA(LP2) = SIG
          ICNT = ICNT + 1
          DSM = MAX(DSM,SIG-SIGIN)
        ENDIF
!
! Bottom of loop on neighbors N2 of N1.
!
    9   IF (LP1 .NE. LPL) GO TO 2
   10   CONTINUE
!
! No errors encountered.
!
      DSMAX = DSM
      IER = ICNT
      RETURN
!
! N < 3
!
   11 DSMAX = 0.
      IER = -1
      RETURN
!
! Nodes N1 and N2 coincide.
!
   12 DSMAX = DSM
      IER = -2
      RETURN
      END
      SUBROUTINE GIVENS ( A,B, C,S)
      REAL A, B, C, S
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88
!
!   This subroutine constructs the Givens plane rotation,
!
!           ( C  S)
!       G = (     ) , where C*C + S*S = 1,
!           (-S  C)
!
! which zeros the second component of the vector (A,B)**T
! (transposed).  Subroutine ROTATE may be called to apply
! the transformation to a 2 by N matrix.
!
!   This routine is identical to Subroutine SROTG from the
! LINPACK BLAS (Basic Linear Algebra Subroutines).
!
! On input:
!
!       A,B = Components of the vector defining the rota-
!             tion.  These are overwritten by values R
!             and Z (described below) which define C and S.
!
! On output:
!
!       A = Signed Euclidean norm R of the input vector:
!           R = +/-SQRT(A*A + B*B)
!
!       B = Value Z such that:
!             C = SQRT(1-Z*Z) and S=Z if ABS(Z) .LE. 1, and
!             C = 1/Z and S = SQRT(1-C*C) if ABS(Z) > 1.
!
!       C = +/-(A/R) or 1 if R = 0.
!
!       S = +/-(B/R) or 0 if R = 0.
!
! Modules required by GIVENS:  None
!
! Intrinsic functions called by GIVENS:  ABS, SQRT
!
!***********************************************************
!
      REAL AA, BB, R, U, V
!
! Local parameters:
!
! AA,BB = Local copies of A and B
! R =     C*A + S*B = +/-SQRT(A*A+B*B)
! U,V =   Variables used to scale A and B for computing R
!
      AA = A
      BB = B
      IF (ABS(AA) .LE. ABS(BB)) GO TO 1
!
! ABS(A) > ABS(B).
!
      U = AA + AA
      V = BB/U
      R = SQRT(.25 + V*V) * U
      C = AA/R
      S = V * (C + C)
!
! Note that R has the sign of A, C > 0, and S has
!   SIGN(A)*SIGN(B).
!
      B = S
      A = R
      RETURN
!
! ABS(A) .LE. ABS(B).
!
    1 IF (BB .EQ. 0.) GO TO 2
      U = BB + BB
      V = AA/U
!
! Store R in A.
!
      A = SQRT(.25 + V*V) * U
      S = BB/A
      C = V * (S + S)
!
! Note that R has the sign of B, S > 0, and C has
!   SIGN(A)*SIGN(B).
!
      B = 1.
      IF (C .NE. 0.) B = 1./C
      RETURN
!
! A = B = 0.
!
    2 C = 1.
      S = 0.
      RETURN
      END


      
      SUBROUTINE GRADC (K,NCC,LCC,N,X,Y,Z,LIST,LPTR,&
     &                  LEND, DX,DY,DXX,DXY,DYY,IER)
      INTEGER K, NCC, LCC(*), N, LIST(*), LPTR(*),&
     &        LEND(N), IER
      REAL    X(N), Y(N), Z(N), DX, DY, DXX, DXY, DYY
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/22/97
!
!   Given a Delaunay triangulation of N points in the plane
! with associated data values Z, this subroutine estimates
! first and second partial derivatives at node K.  The der-
! ivatives are taken to be the partials at K of a cubic
! function which interpolates Z(K) and fits the data values
! at a set of nearby nodes in a weighted least squares
! sense.  A Marquardt stabilization factor is used if neces-
! sary to ensure a well-conditioned system.  Thus, a unique
! solution exists if there are at least 10 noncollinear
! nodes.
!
!   The triangulation may include constraints introduced by
! Subroutine ADDCST, in which case the derivative estimates
! are influenced by the nonconvex geometry of the domain.
! Refer to Subroutine GETNP.  If data values at the con-
! straint nodes are not known, Subroutine ZGRADL, which
! computes approximate data values at constraint nodes along
! with gradients, should be called in place of this routine.
!
!   Subroutine GRADL uses a quadratic polynomial instead of
! the cubic and may be more accurate if the nodal distribu-
! tion is sparse.  Another alternative routine, GRADG,
! employs a global method to compute the first partial
! derivatives at all of the nodes at once.  That method
! is usually more efficient (when all first partials are
! needed) and may be more accurate, depending on the data.
!
! On input:
!
!       K = Index of the node at which derivatives are to be
!           estimated.  1 .LE. K .LE. N.
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC .GE. 0.
!
!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index of the
!             first node of constraint I in LCC(I).  For I =
!             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1.
!
!       N = Number of nodes in the triangulation.
!           N .GE. 10.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       Z = Array of length N containing data values associ-
!           ated with the nodes.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       DX,DY = Estimated first partial derivatives at node
!               K unless IER < 0.
!
!       DXX,DXY,DYY = Estimated second partial derivatives
!                     at node K unless IER < 0.
!
!       IER = Error indicator:
!             IER = L > 0 if no errors were encountered and
!                         L nodes (including node K) were
!                         employed in the least squares fit.
!             IER = -1 if K, NCC, an LCC entry, or N is
!                      outside its valid range on input.
!             IER = -2 if all nodes are collinear.
!
! TRIPACK modules required by GRADC:  GETNP, INTSEC
!
! SRFPACK modules required by GRADC:  GIVENS, ROTATE, SETRO3
!
! Intrinsic functions called by GRADC:  ABS, MIN, REAL, SQRT
!
!***********************************************************
!
      INTEGER   LMN, LMX
      PARAMETER (LMN=14,  LMX=30)
      INTEGER   I, IERR, J, JP1, KK, L, LMAX, LMIN, LM1,&
     &          LNP, NP, NPTS(LMX)
      REAL      A(10,10), C, DIST(LMX), DMIN, DS, DTOL, RIN,&
     &          RS, RTOL, S, SF, SFC, SFS, STF, SUM, W, XK,&
     &          YK, ZK
      DATA      RTOL/1.E-5/, DTOL/.01/
!
! Local parameters:
!
! A =         Transpose of the augmented regression matrix
! C =         First component of the plane rotation deter-
!               mined by Subroutine GIVENS
! DIST =      Array containing the distances between K and
!               the elements of NPTS (refer to GETNP)
! DMIN =      Minimum of the magnitudes of the diagonal
!               elements of the regression matrix after
!               zeros are introduced below the diagonal
! DS =        Squared distance between nodes K and NPTS(LNP)
! DTOL =      Tolerance for detecting an ill-conditioned
!               system.  The system is accepted when DMIN/W
!               .GE. DTOL.
! I =         DO-loop index
! IERR =      Error flag for calls to GETNP
! J =         DO-loop index
! JP1 =       J+1
! KK =        Local copy of K
! L =         Number of columns of A**T to which a rotation
!               is applied
! LMAX,LMIN = Min(LMX,N), Min(LMN,N)
! LMN,LMX =   Minimum and maximum values of LNP for N
!               sufficiently large.  In most cases LMN-1
!               nodes are used in the fit.  4 .LE. LMN .LE.
!               LMX.
! LM1 =       LMIN-1 or LNP-1
! LNP =       Length of NPTS
! NP =        Element of NPTS to be added to the system
! NPTS =      Array containing the indexes of a sequence of
!               nodes ordered by distance from K.  NPTS(1)=K
!               and the first LNP-1 elements of NPTS are
!               used in the least squares fit.  Unless LNP
!               exceeds LMAX, NPTS(LNP) determines R.
! RIN =       Inverse of the distance R between node K and
!               NPTS(LNP) or some point further from K than
!               NPTS(LMAX) if NPTS(LMAX) is used in the fit.
!               R is a radius of influence which enters into
!               the weight W.
! RS =        R*R
! RTOL =      Tolerance for determining R.  If the relative
!               change in DS between two elements of NPTS is
!               not greater than RTOL, they are treated as
!               being the same distance from node K.
! S =         Second component of the plane rotation deter-
!               mined by Subroutine GIVENS
! SF =        Scale factor for the linear terms (columns 8
!               and 9) in the least squares fit -- inverse
!               of the root-mean-square distance between K
!               and the nodes (other than K) in the least
!               squares fit
! SFS =       Scale factor for the quadratic terms (columns
!               5, 6, and 7) in the least squares fit --
!               SF*SF
! SFC =       Scale factor for the cubic terms (first 4
!               columns) in the least squares fit -- SF**3
! STF =       Marquardt stabilization factor used to damp
!               out the first 4 solution components (third
!               partials of the cubic) when the system is
!               ill-conditioned.  As STF increases, the
!               fitting function approaches a quadratic
!               polynomial.
! SUM =       Sum of squared distances between node K and
!               the nodes used in the least squares fit
! W =         Weight associated with a row of the augmented
!               regression matrix -- 1/D - 1/R, where D < R
!               and D is the distance between K and a node
!               entering into the least squares fit
! XK,YK,ZK =  Coordinates and data value associated with K
!
      KK = K
!
! Test for errors and initialize LMIN and LMAX.
!
      IF (KK .LT. 1  .OR.  KK .GT. N  .OR.  NCC .LT. 0&
         &.OR.  N .LT. 10) GO TO 13
      LMIN = MIN(LMN,N)
      LMAX = MIN(LMX,N)
!
! Compute NPTS, DIST, LNP, SF, SFS, SFC, and RIN --
!
!   Set NPTS to the closest LMIN-1 nodes to K.
!
      SUM = 0.
      NPTS(1) = KK
      DIST(1) = 0.
      LM1 = LMIN - 1
      DO 1 LNP = 2,LM1
        CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,&
     &              LNP, NPTS,DIST, IERR)
        IF (IERR .NE. 0) GO TO 13
        DS = DIST(LNP)**2
        SUM = SUM + DS
    1   CONTINUE
!
! Add additional nodes to NPTS until the relative increase
!   in DS is at least RTOL.
!
      DO 3 LNP = LMIN,LMAX
        CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,&
     &              LNP, NPTS,DIST, IERR)
        RS = DIST(LNP)**2
        IF ((RS-DS)/DS .LE. RTOL) GO TO 2
        IF (LNP .GT. 10) GO TO 4
    2   SUM = SUM + RS
    3   CONTINUE
!
! Use all LMAX nodes in the least squares fit.  RS is
!   arbitrarily increased by 10 per cent.
!
      RS = 1.1*RS
      LNP = LMAX + 1
!
! There are LNP-2 equations corresponding to nodes NPTS(2),
!   ...,NPTS(LNP-1).
!
    4 SFS = REAL(LNP-2)/SUM
      SF = SQRT(SFS)
      SFC = SF*SFS
      RIN = 1./SQRT(RS)
      XK = X(KK)
      YK = Y(KK)
      ZK = Z(KK)
!
! A Q-R decomposition is used to solve the least squares
!   system.  The transpose of the augmented regression
!   matrix is stored in A with columns (rows of A) defined
!   as follows:  1-4 are the cubic terms, 5-7 are the quad-
!   ratic terms with coefficients DXX/2, DXY, and DYY/2,
!   8 and 9 are the linear terms with coefficients DX and
!   DY, and the last column is the right hand side.
!
! Set up the first 9 equations and zero out the lower tri-
!   angle with Givens rotations.
!
      DO 6 I = 1,9
        NP = NPTS(I+1)
        W = 1./DIST(I+1) - RIN
        CALL SETRO3 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,&
     &               SFC,W, A(1,I))
        IF (I .EQ. 1) GO TO 6
        DO 5 J = 1,I-1
          JP1 = J + 1
          L = 10 - J
          CALL GIVENS (A(J,J),A(J,I),C,S)
          CALL ROTATE (L,C,S,A(JP1,J),A(JP1,I))
    5     CONTINUE
    6   CONTINUE
!
! Add the additional equations to the system using
!   the last column of A.  I .LE. LNP.
!
      I = 11
    7   IF (I .LT. LNP) THEN
          NP = NPTS(I)
          W = 1./DIST(I) - RIN
          CALL SETRO3 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,&
     &                 SFC,W, A(1,10))
          DO 8 J = 1,9
            JP1 = J + 1
            L = 10 - J
            CALL GIVENS (A(J,J),A(J,10),C,S)
            CALL ROTATE (L,C,S,A(JP1,J),A(JP1,10))
    8       CONTINUE
          I = I + 1
          GO TO 7
        ENDIF
!
! Test the system for ill-conditioning.
!
      DMIN = MIN( ABS(A(1,1)),ABS(A(2,2)),ABS(A(3,3)),&
     &            ABS(A(4,4)),ABS(A(5,5)),ABS(A(6,6)),&
     &            ABS(A(7,7)),ABS(A(8,8)),ABS(A(9,9)) )
      IF (DMIN/W .GE. DTOL) GO TO 12
      IF (LNP .LE. LMAX) THEN
!
!   Add another node to the system and increase R.  Note
!     that I = LNP.
!
        LNP = LNP + 1
        IF (LNP .LE. LMAX) THEN
          CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,&
     &                LNP, NPTS,DIST, IERR)
          RS = DIST(LNP)**2
        ENDIF
        RIN = 1./SQRT(1.1*RS)
        GO TO 7
      ENDIF
!
! Stabilize the system by damping third partials -- add
!   multiples of the first four unit vectors to the first
!   four equations.
!
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
!
! Test the damped system for ill-conditioning.
!
      DMIN = MIN( ABS(A(5,5)),ABS(A(6,6)),ABS(A(7,7)),&
     &            ABS(A(8,8)),ABS(A(9,9)) )
      IF (DMIN/W .LT. DTOL) GO TO 14
!
! Solve the 9 by 9 triangular system for the last 5
!   components (first and second partial derivatives).
!
   12 DY = A(10,9)/A(9,9)
      DX = (A(10,8)-A(9,8)*DY)/A(8,8)
      DYY = (A(10,7)-A(8,7)*DX-A(9,7)*DY)/A(7,7)
      DXY = (A(10,6)-A(7,6)*DYY-A(8,6)*DX-A(9,6)*DY)/A(6,6)
      DXX = (A(10,5)-A(6,5)*DXY-A(7,5)*DYY-A(8,5)*DX-&
     &       A(9,5)*DY)/A(5,5)
!
! Scale the solution components.
!
      DX = SF*DX
      DY = SF*DY
      DXX = 2.*SFS*DXX
      DXY = SFS*DXY
      DYY = 2.*SFS*DYY
      IER = LNP - 1
      RETURN
!
! Invalid input parameter.
!
   13 IER = -1
      RETURN
!
! No unique solution due to collinear nodes.
!
   14 IER = -2
      RETURN
      END
      SUBROUTINE GRADG (NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,&
     &                  IFLGS,SIGMA, NIT,DGMAX,GRAD, IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),&
     &        IFLGS, NIT, IER
      REAL    X(N), Y(N), Z(N), SIGMA(*), DGMAX, GRAD(2,N)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/12/94
!
!   Given a triangulation of N nodes in the plane, along
! with data values at the nodes and tension factors associ-
! ated with the arcs, this subroutine employs a global
! method to compute estimated gradients at the nodes.  The
! method consists of minimizing a quadratic functional Q(G)
! over vectors G of length 2N (gradient components), where
! Q is an approximation to the linearized curvature over the
! triangulation of a C-1 bivariate function F(X,Y) which
! interpolates the nodal values and gradients.
!
!   The restriction of F to an arc of the triangulation is
! taken to be the Hermite interpolatory tension spline
! defined by the data values and tangential gradient compo-
! nents at the endpoints of the arc, and Q is the sum over
! the triangulation arcs, excluding interior constraint
! arcs, of the linearized curvatures of F along the arcs --
! the integrals over the arcs of D2F(T)**2, where D2F(T) is
! the second derivative of F with respect to distance T
! along the arc.
!
!   Subroutines INTRC1 and UNIF may be called to evaluate F
! at arbitrary points.  The interpolant F is further de-
! scribed in Subroutines FVAL and TVAL, and Q is identical
! to the functional Q1 described in Subroutine SMSURF.
!
!   The minimization problem corresponds to an order 2N
! symmetric positive definite sparse linear system which is
! solved for the X and Y partial derivatives by the block
! Gauss-Seidel method with N blocks of order 2.
!
!   If constraints are present and data values at the con-
! straint nodes are not known, Subroutine ZGRADG, which
! computes approximate data values at constraint nodes
! along with the gradients, should be called in place of
! this routine.
!
!   An alternative method, Subroutine GRADC or GRADL, com-
! putes a local approximation to the partials at a single
! node and may be more accurate, depending on the data
! values and distribution of nodes (neither method emerged
! as superior in tests for accuracy).  If all gradients are
! required and a uniform tension factor SIGMA = 0 is used,
! GRADG is significantly faster than either GRADC or GRADL.
!
! On input:
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC .GE. 0.
!
!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index of the
!             first node of constraint I in LCC(I).  For I =
!             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       Z = Array of length N containing data values at the
!           nodes.  Z(I) is associated with (X(I),Y(I)).
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
!       IFLGS = Tension factor option:
!               IFLGS .LE. 0 if a single uniform tension
!                            factor is to be used.
!               IFLGS .GE. 1 if variable tension is desired.
!
!       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
!               array containing tension factors associated
!               with arcs in one-to-one correspondence with
!               LIST entries (IFLGS .GE. 1).  Refer to Sub-
!               routines GETSIG, SIG0, SIG1, and SIG2.
!
! The above parameters are not altered by this routine.
!
!       NIT = Maximum number of Gauss-Seidel iterations to
!             be employed.  This maximum will likely be
!             achieved if DGMAX is smaller than the machine
!             precision.  Note that complete convergence is
!             not necessary to achieve maximum accuracy of
!             the interpolant.  For SIGMA = 0, optimal ef-
!             ficiency was achieved in testing with DGMAX =
!             0, and NIT = 3 or 4.  NIT > 0.
!
!       DGMAX = Nonnegative convergence criterion.  The
!               method is terminated when the maximum change
!               in a gradient between iterations is at most
!               DGMAX.  The change in a gradient is taken to
!               be the Euclidean norm of the difference rel-
!               ative to 1 plus the norm of the old value.
!               DGMAX = 1.E-3 is sufficient for effective
!               convergence.
!
!       GRAD = 2 by N array whose columns contain initial
!              estimates of the partial derivatives.  Zero
!              vectors are sufficient.
!
! On output:
!
!       NIT = Number of Gauss-Seidel iterations employed.
!
!       DGMAX = Maximum relative change in a gradient at the
!               last iteration.
!
!       GRAD = Estimated X and Y partial derivatives at the
!              nodes with X partials in the first row.  Grad
!              is not altered if IER = -1.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and the
!                     convergence criterion was achieved.
!             IER = 1 if no errors were encountered but con-
!                     vergence was not achieved within NIT
!                     iterations.
!             IER = -1 if NCC, an LCC entry, N, NIT, or
!                      DGMAX is outside its valid range on
!                      input.
!             IER = -2 if all nodes are collinear or the
!                      triangulation data structure is in-
!                      valid.
!             IER = -3 if duplicate nodes were encountered.
!
! SRFPACK modules required by GRADG:  GRCOEF, SNHCSH
!
! Intrinsic functions called by GRADG:  ABS, MAX, SQRT
!
!***********************************************************
!
      INTEGER I, IFL, IFRST, ILAST, ITER, J, K, KBAK, KFOR,&
     &        LCC1, LP, LPJ, LPL, MAXIT, NB, NN
      REAL    A11, A12, A22, D, DCUB, DELX, DELXS, DELY,&
     &        DELYS, DET, DF, DGMX, DSQ, DZX, DZY, R1, R2,&
     &        SDF, SIG, T, TOL, XK, YK, ZK, ZXK, ZYK
!
      NN = N
      IFL = IFLGS
      MAXIT = NIT
      TOL = DGMAX
!
! Test for errors in input parameters.
!
      IF (NCC .LT. 0  .OR.  MAXIT .LT. 1  .OR.  TOL .LT. 0.)&
     &  GO TO 9
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
!
! Initialize iteration count and SIG (overwritten if
!   IFLGS > 0).
!
      ITER = 0
      SIG = SIGMA(1)
!
! Top of iteration loop:  If K is a constraint node, I
!   indexes the constraint containing node K, IFRST and
!   ILAST are the first and last nodes of constraint I,
!   and (KBAK,K,KFOR) is a subsequence of constraint I.
!
    2 IF (ITER .EQ. MAXIT) GO TO 8
      DGMX = 0.
      I = 0
      IFRST = 1
      ILAST = LCC1-1
      KBAK = 0
      KFOR = 0
!
! Loop on nodes.
!
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
!
!   Initialize components of the order 2 system for the
!     change (DZX,DZY) in the K-th solution components
!     (symmetric matrix in A and residual in R).
!
        A11 = 0.
        A12 = 0.
        A22 = 0.
        R1 = 0.
        R2 = 0.
!
!   Loop on neighbors J of node K.
!
        LPL = LEND(K)
        LPJ = LPL
    3   LPJ = LPTR(LPJ)
          J = ABS(LIST(LPJ))
!
!   Arc K-J lies in a constraint region and is bypassed iff
!     K and J are nodes in the same constraint and J follows
!     KFOR and precedes KBAK as a neighbor of K.
!
          IF (K .LT. LCC1  .OR.  J .LT. IFRST  .OR.&
     &        J .GT. ILAST) GO TO 5
          IF (J .EQ. KBAK  .OR.  J .EQ. KFOR) GO TO 5
          LP = LPJ
!
    4     LP = LPTR(LP)
            NB = ABS(LIST(LP))
            IF (NB .EQ. KBAK) GO TO 6
            IF (NB .NE. KFOR) GO TO 4
!
!   Compute parameters associated with edge
!     K->J, and test for duplicate nodes.
!
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
!
!   Update the system components for node J.  The contribu-
!     tion from edge K->J is weighted by 1/D, where D is
!     the arc length.
!
          A11 = A11 + DF*DELXS/D
          A12 = A12 + DF*DELX*DELY/D
          A22 = A22 + DF*DELYS/D
          T = ((DF+SDF)*(Z(J)-ZK) - DF*(ZXK*DELX + ZYK*DELY)&
     &          - SDF*(GRAD(1,J)*DELX + GRAD(2,J)*DELY))/D
          R1 = R1 + T*DELX
          R2 = R2 + T*DELY
!
!   Bottom of loop on neighbors.
!
    6     IF (LPJ .NE. LPL) GO TO 3
!
!   Solve the system associated with the K-th block.
!
        DET = A11*A22 - A12*A12
        IF (DET .EQ. 0.  .OR.  A11 .EQ. 0.) GO TO 10
        DZY = (A11*R2 - A12*R1)/DET
        DZX = (R1 - A12*DZY)/A11
!
!   Update the partials at node K and the maximum relative
!     change DGMX.
!
        GRAD(1,K) = ZXK + DZX
        GRAD(2,K) = ZYK + DZY
        DGMX = MAX(DGMX,SQRT(DZX*DZX+DZY*DZY)/&
     &             (1.+SQRT(ZXK*ZXK+ZYK*ZYK)))
    7   CONTINUE
!
!   Increment ITER and test for convergence.
!
      ITER = ITER + 1
      IF (DGMX .GT. TOL) GO TO 2
!
! Method converged.
!
      NIT = ITER
      DGMAX = DGMX
      IER = 0
      RETURN
!
! Method failed to converge within NIT iterations.
!
    8 DGMAX = DGMX
      IER = 1
      RETURN
!
! Invalid input parameter.
!
    9 NIT = 0
      DGMAX = 0.
      IER = -1
      RETURN
!
! Node K and its neighbors are collinear, resulting in a
!   singular system.
!
   10 NIT = 0
      DGMAX = DGMX
      IER = -2
      RETURN
!
! Nodes K and J coincide.
!
   11 NIT = 0
      DGMAX = DGMX
      IER = -3
      RETURN
      END
      SUBROUTINE GRADL (K,NCC,LCC,N,X,Y,Z,LIST,LPTR,&
     &                  LEND, DX,DY,IER)
      INTEGER K, NCC, LCC(*), N, LIST(*), LPTR(*),&
     &        LEND(N), IER
      REAL    X(N), Y(N), Z(N), DX, DY
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/22/97
!
!   Given a Delaunay triangulation of N points in the plane
! with associated data values Z, this subroutine estimates
! X and Y partial derivatives at node K.  The derivatives
! are taken to be the partials at K of a quadratic function
! which interpolates Z(K) and fits the data values at a set
! of nearby nodes in a weighted least squares sense. A Mar-
! quardt stabilization factor is used if necessary to ensure
! a well-conditioned system.  Thus, a unique solution exists
! if there are at least 6 noncollinear nodes.
!
!   The triangulation may include constraints introduced by
! Subroutine ADDCST, in which case the gradient estimates
! are influenced by the nonconvex geometry of the domain.
! Refer to Subroutine GETNP.  If data values at the con-
! straint nodes are not known, Subroutine ZGRADL, which
! computes approximate data values at constraint nodes along
! with gradients, should be called in place of this routine.
!
!   Subroutine GRADC uses a cubic polynomial instead of the
! quadratic and is generally more accurate than this routine
! if the nodal distribution is sufficiently dense.  Another
! alternative routine, GRADG, employs a global method to
! compute the partial derivatives at all of the nodes at
! once.  That method is usually more efficient (when all
! partials are needed) and may be more accurate, depending
! on the data.
!
! On input:
!
!       K = Index of the node at which derivatives are to be
!           estimated.  1 .LE. K .LE. N.
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC .GE. 0.
!
!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index of the
!             first node of constraint I in LCC(I).  For I =
!             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1.
!
!       N = Number of nodes in the triangulation.  N .GE. 6.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       Z = Array of length N containing data values associ-
!           ated with the nodes.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       DX,DY = Estimated partial derivatives at node K
!               unless IER < 0.
!
!       IER = Error indicator:
!             IER = L > 0 if no errors were encountered and
!                         L nodes (including node K) were
!                         employed in the least squares fit.
!             IER = -1 if K, NCC, an LCC entry, or N is
!                      outside its valid range on input.
!             IER = -2 if all nodes are collinear.
!
! TRIPACK modules required by GRADL:  GETNP, INTSEC
!
! SRFPACK modules required by GRADL:  GIVENS, ROTATE, SETRO1
!
! Intrinsic functions called by GRADL:  ABS, MIN, REAL, SQRT
!
!***********************************************************
!
      INTEGER   LMN, LMX
      PARAMETER (LMN=10,  LMX=30)
      INTEGER   I, IERR, J, JP1, KK, L, LMAX, LMIN, LM1,&
     &          LNP, NP, NPTS(LMX)
      REAL      A(6,6), C, DIST(LMX), DMIN, DS, DTOL, RIN,&
     &          RS, RTOL, S, SF, SFS, STF, SUM, W, XK, YK,&
     &          ZK
      DATA      RTOL/1.E-5/, DTOL/.01/
!
! Local parameters:
!
! A =         Transpose of the augmented regression matrix
! C =         First component of the plane rotation deter-
!               mined by Subroutine GIVENS
! DIST =      Array containing the distances between K and
!               the elements of NPTS (refer to GETNP)
! DMIN =      Minimum of the magnitudes of the diagonal
!               elements of the regression matrix after
!               zeros are introduced below the diagonal
! DS =        Squared distance between nodes K and NPTS(LNP)
! DTOL =      Tolerance for detecting an ill-conditioned
!               system.  The system is accepted when DMIN/W
!               .GE. DTOL
! I =         DO-loop index
! IERR =      Error flag for calls to GETNP
! J =         DO-loop index
! JP1 =       J+1
! KK =        Local copy of K
! L =         Number of columns of A**T to which a rotation
!               is applied
! LMAX,LMIN = Min(LMX,N), Min(LMN,N)
! LMN,LMX =   Minimum and maximum values of LNP for N
!               sufficiently large.  In most cases LMN-1
!               nodes are used in the fit.  4 .LE. LMN .LE.
!               LMX.
! LM1 =       LMIN-1 or LNP-1
! LNP =       Length of NPTS
! NP =        Element of NPTS to be added to the system
! NPTS =      Array containing the indexes of a sequence of
!               nodes ordered by distance from K.  NPTS(1)=K
!               and the first LNP-1 elements of NPTS are
!               used in the least squares fit.  Unless LNP
!               exceeds LMAX, NPTS(LNP) determines R.
! RIN =       Inverse of the distance R between node K and
!               NPTS(LNP) or some point further from K than
!               NPTS(LMAX) if NPTS(LMAX) is used in the fit.
!               R is a radius of influence which enters into
!               the weight W.
! RS =        R*R
! RTOL =      Tolerance for determining R.  If the relative
!               change in DS between two elements of NPTS is
!               not greater than RTOL, they are treated as
!               being the same distance from node K
! S =         Second component of the plane rotation deter-
!               mined by Subroutine GIVENS
! SF =        Scale factor for the linear terms (columns 4
!               and 5) in the least squares fit -- inverse
!               of the root-mean-square distance between K
!               and the nodes (other than K) in the least
!               squares fit.
! SFS =       Scale factor for the quadratic terms (first 3
!               columns) in the least squares fit -- SF*SF.
! STF =       Marquardt stabilization factor used to damp
!               out the first 3 solution components (second
!               partials of the quadratic) when the system
!               is ill-conditioned.  As STF increases, the
!               fitting function approaches a linear
! SUM =       Sum of squared distances between node K and
!               the nodes used in the least squares fit
! W =         Weight associated with a row of the augmented
!               regression matrix -- 1/R - 1/D, where D < R
!               and D is the distance between K and a node
!               entering into the least squares fit.
! XK,YK,ZK =  Coordinates and data value associated with K
!
      KK = K
!
! Test for errors and initialize LMIN and LMAX.
!
      IF (KK .LT. 1  .OR.  KK .GT. N  .OR.  NCC .LT. 0&
         &.OR.  N .LT. 6) GO TO 13
      LMIN = MIN(LMN,N)
      LMAX = MIN(LMX,N)
!
! Compute NPTS, DIST, LNP, SF, SFS, and RIN --
!
!   Set NPTS to the closest LMIN-1 nodes to K.
!
      SUM = 0.
      NPTS(1) = KK
      DIST(1) = 0.
      LM1 = LMIN - 1
      DO 1 LNP = 2,LM1
        CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,&
     &              LNP, NPTS,DIST, IERR)
        IF (IERR .NE. 0) GO TO 13
        DS = DIST(LNP)**2
        SUM = SUM + DS
    1   CONTINUE
!
! Add additional nodes to NPTS until the relative increase
!   in DS is at least RTOL.
!
      DO 3 LNP = LMIN,LMAX
        CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,&
     &              LNP, NPTS,DIST, IERR)
        RS = DIST(LNP)**2
        IF ((RS-DS)/DS .LE. RTOL) GO TO 2
        IF (LNP .GT. 6) GO TO 4
    2   SUM = SUM + RS
    3   CONTINUE
!
! Use all LMAX nodes in the least squares fit.  RS is
!   arbitrarily increased by 10 per cent.
!
      RS = 1.1*RS
      LNP = LMAX + 1
!
! There are LNP-2 equations corresponding to nodes NPTS(2),
!   ...,NPTS(LNP-1).
!
    4 SFS = REAL(LNP-2)/SUM
      SF = SQRT(SFS)
      RIN = 1./SQRT(RS)
      XK = X(KK)
      YK = Y(KK)
      ZK = Z(KK)
!
! A Q-R decomposition is used to solve the least squares
!   system.  The transpose of the augmented regression
!   matrix is stored in A with columns (rows of A) defined
!   as follows:  1-3 are the quadratic terms, 4 and 5 are
!   the linear terms with coefficients DX and DY, and the
!   last column is the right hand side.
!
! Set up the first 5 equations and zero out the lower tri-
!   angle with Givens rotations.
!
      DO 6 I = 1,5
        NP = NPTS(I+1)
        W = 1./DIST(I+1) - RIN
        CALL SETRO1 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,&
     &               W, A(1,I))
        IF (I .EQ. 1) GO TO 6
        DO 5 J = 1,I-1
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS (A(J,J),A(J,I),C,S)
          CALL ROTATE (L,C,S,A(JP1,J),A(JP1,I))
    5     CONTINUE
    6   CONTINUE
!
! Add the additional equations to the system using
!   the last column of A.  I .LE. LNP.
!
      I = 7
    7   IF (I .LT. LNP) THEN
          NP = NPTS(I)
          W = 1./DIST(I) - RIN
          CALL SETRO1 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,&
     &                 W, A(1,6))
          DO 8 J = 1,5
            JP1 = J + 1
            L = 6 - J
            CALL GIVENS (A(J,J),A(J,6),C,S)
            CALL ROTATE (L,C,S,A(JP1,J),A(JP1,6))
    8       CONTINUE
          I = I + 1
          GO TO 7
        ENDIF
!
! Test the system for ill-conditioning.
!
      DMIN = MIN( ABS(A(1,1)),ABS(A(2,2)),ABS(A(3,3)),&
     &            ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN/W .GE. DTOL) GO TO 12
      IF (LNP .LE. LMAX) THEN
!
!   Add another node to the system and increase R.  Note
!     that I = LNP.
!
        LNP = LNP + 1
        IF (LNP .LE. LMAX) THEN
          CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,&
     &                LNP, NPTS,DIST, IERR)
          RS = DIST(LNP)**2
        ENDIF
        RIN = 1./SQRT(1.1*RS)
        GO TO 7
      ENDIF
!
! Stabilize the system by damping second partials -- add
!   multiples of the first three unit vectors to the first
!   three equations.
!
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
!
! Test the damped system for ill-conditioning.
!
      DMIN = MIN( ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN/W .LT. DTOL) GO TO 14
!
! Solve the 2 by 2 triangular system for the partial
!   derivatives.
!
   12 DY = A(6,5)/A(5,5)
      DX = SF*(A(6,4) - A(5,4)*DY)/A(4,4)
      DY = SF*DY
      IER = LNP - 1
      RETURN
!
! Invalid input parameter.
!
   13 IER = -1
      RETURN
!
! No unique solution due to collinear nodes.
!
   14 IER = -2
      RETURN
      END
      SUBROUTINE GRCOEF (SIGMA,DCUB, D,SD)
      REAL SIGMA, DCUB, D, SD
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/18/96
!
!   This subroutine computes factors involved in the linear
! system solved by Subroutines GRADG and SMSGS.
!
! On input:
!
!       SIGMA = Nonnegative tension factor associated with a
!               triangulation arc.
!
!       DCUB = Cube of the positive arc length.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       D = Diagonal factor.  D = SIG*(SIG*COSHM(SIG) -
!           SINHM(SIG))/(E*DCUB), where E = SIG*SINH(SIG) -
!           2*COSHM(SIG).  D > 0.
!
!       SD = Off-diagonal factor.  SD = SIG*SINHM(SIG)/
!            (E*DCUB).  SD > 0.
!
! SRFPACK module required by GRCOEF:  SNHCSH
!
! Intrinsic function called by GRCOEF:  EXP
!
!***********************************************************
!
      REAL COSHM, COSHMM, E, EMS, SCM, SIG, SINHM, SSINH,&
     &     SSM
!
      SIG = SIGMA
      IF (SIG .LT. 1.E-9) THEN
!
! SIG = 0:  cubic interpolant.
!
        D = 4./DCUB
        SD = 2./DCUB
      ELSEIF (SIG .LE. .5) THEN
!
! 0 .LT. SIG .LE. .5:  use approximations designed to avoid
!                      cancellation error in the hyperbolic
!                      functions when SIGMA is small.
!
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        E = (SIG*SINHM - COSHMM - COSHMM)*DCUB
        D = SIG*(SIG*COSHM-SINHM)/E
        SD = SIG*SINHM/E
      ELSE
!
! SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
!            to avoid overflow when SIGMA is large.
!
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
      SUBROUTINE INTRC0 (PX,PY,NCC,LCC,N,X,Y,Z,LIST,LPTR,&
     &                   LEND, IST, PZ,IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),&
     &        IST, IER
      REAL    PX, PY, X(N), Y(N), Z(N), PZ
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/03/98
!
!   Given a triangulation of a set of nodes in the plane,
! along with data values at the nodes, this subroutine com-
! putes the value at P = (PX,PY) of the piecewise linear
! function which interpolates the data values.  The surface
! is extended in a continuous fashion beyond the boundary of
! the triangulation, allowing extrapolation.
!
! On input:
!
!       PX,PY = Coordinates of the point P at which the sur-
!               face is to be evaluated.
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC .GE. 0.
!
!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index of the
!             first node of constraint I in LCC(I).  For I =
!             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       Z = Array of length N containing data values at the
!           nodes.  Refer to Subroutine ZGRADL.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
! The above parameters are not altered by this routine.
!
!       IST = Index of the starting node in the search for a
!             triangle containing P.  1 .LE. IST .LE. N.
!             The output value of IST from a previous call
!             may be a good choice.
!
! On output:
!
!       IST = Index of one of the vertices of the triangle
!             containing P (or a boundary node which is vis-
!             ible from P) unless IER < 0.
!
!       PZ = Value of the interpolatory surface at P, or
!            zero if IER < 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and P is
!                     contained in a triangle but not in a
!                     constraint region.
!             IER = 1 if no errors were encountered and P
!                     lies in a constraint region triangle.
!                     PZ is effectively an extrapolated
!                     value in this case.
!             IER = 2 if no errors were encountered and P is
!                     exterior to the triangulation.  PZ is
!                     an extrapolated value in this case.
!             IER = -1 if NCC, N, or IST is outside its
!                      valid range on input.  LCC is not
!                      tested for validity.
!             IER = -2 if the nodes are collinear.
!
! TRIPACK modules required by INTRC0:  CRTRI, JRAND, LEFT,
!                                        LSTPTR, TRFIND
!
! SRFPACK module required by INTRC0:  COORDS
!
!***********************************************************
!
      LOGICAL CRTRI
      INTEGER I1, I2, I3, IERR, LPL, N1, N2
      REAL    B1, B2, B3, DP, X1, X2, XP, Y1, Y2, YP
!
      XP = PX
      YP = PY
      PZ = 0.
!
! Test for invalid input parameters.
!
      IF (NCC .LT. 0  .OR.  N .LT. 3  .OR.  IST .LT. 1&
         &.OR.  IST .GT. N) THEN
        IER = -1
        RETURN
      ENDIF
!
! Find a triangle (I1,I2,I3) containing P, or a pair of
!   visible boundary nodes I1 and I2.
!
      CALL TRFIND (IST,XP,YP,N,X,Y,LIST,LPTR,LEND, I1,I2,I3)
      IF (I1 .EQ. 0) THEN
        IER = -2
        RETURN
      ENDIF
      IST = I1
      IF (I3 .EQ. 0) GO TO 1
!
! P is in a triangle.  Compute its barycentric coordinates.
!
      CALL COORDS (XP,YP,X(I1),X(I2),X(I3),Y(I1),Y(I2),&
     &             Y(I3), B1,B2,B3,IERR)
      IF (IERR .NE. 0) THEN
        IER = -2
        RETURN
      ENDIF
!
! Compute an interpolated value.
!
      PZ = B1*Z(I1) + B2*Z(I2) + B3*Z(I3)
      IER = 0
!
      IF (CRTRI(NCC,LCC,I1,I2,I3)) THEN
        IER = 1
      ELSE
        IER = 0
      ENDIF
      RETURN
!
! P is exterior to the triangulation.  Extrapolate to P by
!   extending the interpolatory surface as a constant
!   beyond the boundary:  PZ is the function value at Q
!   where Q is the closest boundary point to P.
!
! Determine Q by traversing the boundary starting from the
!   rightmost visible node I1.
!
    1 IER = 2
      N2 = I1
!
! Top of loop:
!
!   Set N1 to the last neighbor of N2, and compute the dot
!     product DP = (N2->N1,N2->P).  P FORWARD N2->N1 iff
!     DP > 0.
!
    2 LPL = LEND(N2)
      N1 = -LIST(LPL)
      X1 = X(N1)
      Y1 = Y(N1)
      X2 = X(N2)
      Y2 = Y(N2)
      DP = (X1-X2)*(XP-X2) + (Y1-Y2)*(YP-Y2)
      IF (DP .LE. 0.) THEN
!
!   N2 is the closest boundary point to P.
!
        PZ = Z(N2)
        RETURN
      ENDIF
!
!   P FORWARD N2->N1.  Test for P FORWARD N1->N2.
!
      IF ((XP-X1)*(X2-X1) + (YP-Y1)*(Y2-Y1) .GT. 0.) THEN
!
!   The closest boundary point to P lies on N2-N1.  Compute
!     its local coordinates with respect to N2-N1.
!
        B1 = DP/( (X2-X1)**2 + (Y2-Y1)**2 )
        B2 = 1. - B1
        PZ = B1*Z(N1) + B2*Z(N2)
        RETURN
      ENDIF
!
!   Bottom of boundary traversal loop.
!
      N2 = N1
      GO TO 2
      END
      SUBROUTINE INTRC1 (PX,PY,NCC,LCC,N,X,Y,Z,LIST,LPTR,&
     &                   LEND,IFLGS,SIGMA,GRAD,&
     &                   DFLAG, IST, PZ,PZX,PZY,IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),&
     &        IFLGS, IST, IER
      LOGICAL DFLAG
      REAL    PX, PY, X(N), Y(N), Z(N), SIGMA(*), GRAD(2,N),&
     &        PZ, PZX, PZY
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/03/98
!
!   Given a triangulation of a set of nodes in the plane,
! along with data values and estimated gradients at the
! nodes, this subroutine computes the value and, optionally,
! the first partial derivatives at P = (PX,PY) of a C-1
! (once-continuously differentiable) function F which int-
! erpolates the data values and gradients.  Extrapolation to
! a point exterior to the triangulation is accomplished by
! extending the surface in such a way that F is C-1 over the
! entire plane.
!
!   Subroutine FVAL is used to evaluate an interpolatory
! surface under tension, while Subroutine TVAL is called in
! the case of no tension (IFLGS .LE. 0 and SIGMA(1) = 0).
! However, the surface under tension is well-defined with
! SIGMA = 0, and, in this case, the two interpolants are
! identical on triangulation arcs and outside the triangula-
! tion.  The use of FVAL with no tension can be forced (at
! a cost in efficiency) by setting IFLGS = 1 and storing
! zeros in all components of SIGMA.  Note, however, that
! first partial derivatives are only available from TVAL
! (and at points outside the triangulation);  i.e., a proce-
! dure for differentiating the surface under tension has not
! been implemented.
!
!   A set of interpolated values at the vertices of a rec-
! tangular grid can be obtained by a single call to
! Subroutine UNIF.  Subroutine INTRC0 provides for evalua-
! tion of the piecewise linear interpolatory surface.
!
! On input:
!
!       PX,PY = Coordinates of the point P at which the sur-
!               face is to be evaluated.
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC .GE. 0.
!
!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index of the
!             first node of constraint I in LCC(I).  For I =
!             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       Z = Array of length N containing data values at the
!           nodes.  Refer to Subroutines ZGRADG and ZGRADL.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
!       IFLGS = Tension factor option:
!               IFLGS .LE. 0 if a single uniform tension
!                            factor is to be used.
!               IFLGS .GE. 1 if variable tension is desired.
!
!       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
!               array containing tension factors associated
!               with arcs in one-to-one correspondence with
!               LIST entries (IFLGS .GE. 1).  Refer to Sub-
!               routines FVAL, GETSIG, SIG0, SIG1, and SIG2.
!
!       GRAD = 2 by N array whose columns contain estimated
!              gradients at the nodes with X partial deriva-
!              tives in the first row and Y partials in the
!              second.  Refer to Subroutines GRADC, GRADG,
!              GRADL, SMSURF, ZGRADG, and ZGRADL.
!
!       DFLAG = Logical flag which specifies whether first
!               partial derivatives at P are to be computed:
!               DFLAG = TRUE if and only if partials are
!               to be computed by TVAL.  This option is only
!               valid for IFLGS .LE. 0 and SIGMA(1) = 0 (and
!               for points outside the triangulation).
!
! The above parameters are not altered by this routine.
!
!       IST = Index of the starting node in the search for a
!             triangle containing P.  1 .LE. IST .LE. N.
!             The output value of IST from a previous call
!             may be a good choice.
!
! On output:
!
!       IST = Index of one of the vertices of the triangle
!             containing P (or a boundary node which is vis-
!             ible from P) unless IER = -1 or IER = -2.
!
!       PZ = Value of the interpolatory surface at P, or
!            zero if IER < 0.
!
!       PZX,PZY = X and Y partials at P if DFLAG = .TRUE.
!                 and IER .GE. 0, unaltered otherwise.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and P is
!                     contained in a triangle but not in a
!                     constraint region.
!             IER = 1 if no errors were encountered and P
!                     lies in a constraint region triangle.
!                     PZ is effectively an extrapolated
!                     value in this case.
!             IER = 2 if no errors were encountered and P is
!                     exterior to the triangulation.  PZ is
!                     an extrapolated value in this case.
!             IER = -1 if NCC, N, or IST is outside its
!                      valid range on input.  LCC is not
!                      tested for validity.
!             IER = -2 if the nodes are collinear.
!             IER = -3 if P is contained in a triangle and
!                      DFLAG = TRUE, but IFLGS > 0 or
!                      SIGMA(1) .NE. 0.
!
! TRIPACK modules required by INTRC1:  CRTRI, JRAND, LEFT,
!                                        LSTPTR, TRFIND
!
! SRFPACK modules required by INTRC1:  ARCINT, COORDS, FVAL,
!                                        SNHCSH, TVAL
!
! Intrinsic function called by INTRC1:  SQRT
!
!***********************************************************
!
      INTEGER LSTPTR
      LOGICAL CRTRI
      INTEGER I1, I2, I3, IERR, LP, LPL, N1, N2, N3
      LOGICAL TENSN
      REAL    A1, A2, B1, B2, C1, C2, D, D1, D2, D3, DP,&
     &        DP1, DP3, F1, F2, R1, R12, R2, SIG, SIG1,&
     &        SIG2, SIG3, T, T1, T2, X1, X2, X3, X12, X23,&
     &        X2P, XP, XQ, XQP, Y1, Y2, Y3, Y12, Y23, Y2P,&
     &        YP, YQ, YQP, Z1, Z2, Z3, ZQ, ZX1, ZX2, ZX3,&
     &        ZXQ, ZY1, ZY2, ZY3, ZYQ
!
      XP = PX
      YP = PY
      PZ = 0.
!
! Test for invalid input parameters.
!
      IF (NCC .LT. 0  .OR.  N .LT. 3  .OR.  IST .LT. 1&
         &.OR.  IST .GT. N) THEN
        IER = -1
        RETURN
      ENDIF
!
! Find a triangle (I1,I2,I3) containing P, or a pair of
!   visible boundary nodes I1 and I2.
!
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
!
! P is in a triangle.  Store local parameters for the
!   call to FVAL or TVAL.
!
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
!
! Set SIG1, SIG2, and SIG3 to the tension factors associated
!   with the sides opposite I1, I2, and I3, respectively,
!   and compute a value from FVAL.
!
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
        CALL FVAL (XP,YP,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,ZX1,ZX2,&
     &             ZX3,ZY1,ZY2,ZY3,SIG1,SIG2,SIG3, PZ,IERR)
        IF (IERR .LT. 0) THEN
          IER = -2
          RETURN
        ENDIF
      ELSE
!
! Compute an interpolated value from TVAL for no tension.
!
        CALL TVAL (XP,YP,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,ZX1,ZX2,&
     &             ZX3,ZY1,ZY2,ZY3,DFLAG, PZ,PZX,PZY,IERR)
        IF (IERR .NE. 0) THEN
          IER = -2
          RETURN
        ENDIF
      ENDIF
!
      IF (CRTRI(NCC,LCC,I1,I2,I3)) THEN
        IER = 1
      ELSE
        IER = 0
      ENDIF
      RETURN
!
! P is exterior to the triangulation.  Extrapolate to P by
!   passing a linear function of one variable through the
!   value and directional derivative (in the direction
!   Q->P) of the interpolatory surface F at Q, where Q is
!   the closest boundary point to P.
!
! Determine Q by traversing the boundary starting from the
!   rightmost visible node I1.
!
    1 IER = 2
      N2 = I1
!
! Top of loop:
!
!   Set N1 to the last neighbor of N2, and compute the dot
!     product DP = (N2->N1,N2->P).  P FORWARD N2->N1 iff
!     DP > 0.
!
    2 LPL = LEND(N2)
      N1 = -LIST(LPL)
      X1 = X(N1)
      Y1 = Y(N1)
      X2 = X(N2)
      Y2 = Y(N2)
      DP = (X1-X2)*(XP-X2) + (Y1-Y2)*(YP-Y2)
      IF (DP .LE. 0.) THEN
!
!   N2 is the closest boundary point to P:  P lies in a
!     wedge with sides orthogonal to N1-N2 and N2-N3, where
!     N3 is the first neighbor of N2.  The linear interpo-
!     lant must be modified by a correction term which
!     provides for continuity of the derivative across the
!     sides of the wedge.
!
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
          IF (D .NE. 0.) PZX = PZX + (T*X2P - D1*X12 -&
     &                                D2*X23)/D
          PZY = ZY2
          IF (D .NE. 0.) PZY = PZY + (T*Y2P - D1*Y12 -&
     &                                D2*Y23)/D
        ENDIF
        RETURN
      ENDIF
!
!   P FORWARD N2->N1.  Test for P FORWARD N1->N2.
!
      IF ((XP-X1)*(X2-X1) + (YP-Y1)*(Y2-Y1) .LE. 0.) THEN
!
!   Bottom of boundary traversal loop.
!
        N2 = N1
        GO TO 2
      ENDIF
!
! The closest boundary point Q lies on N2-N1.  Store par-
!   tials at N1 and N2, and compute Q and its barycentric
!   coordinates R1 and R2.
!
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
!
!   Set SIG to the tension factor associated with N1-N2 and
!     compute an interpolated value ZQ at Q from FVAL.
!
        IF (IFLGS .LE. 0) THEN
          SIG = SIGMA(1)
        ELSE
          SIG = SIGMA(LPL)
        ENDIF
        CALL ARCINT (R1,X1,X2,Y1,Y2,Z(N1),Z(N2),ZX1,ZX2,&
     &               ZY1,ZY2,SIG,.TRUE., ZQ,ZXQ,ZYQ,IERR)
!
!   Compute the extrapolated value at P.
!
        XQP = XP-XQ
        YQP = YP-YQ
        PZ = ZQ + ZXQ*XQP + ZYQ*YQP
        IF (DFLAG) THEN
          T = ((ZX2-ZX1)*XQP + (ZY2-ZY1)*YQP)/D2
          PZX = ZXQ + X12*T
          PZY = ZYQ + Y12*T
        ENDIF
      ELSE
!
!   Compute the cardinal function values and interpolated
!     value at Q associated with TVAL.
!
        R12 = R1*R2
        F1 = R1*R12
        F2 = R2*R12
        A1 = R1 + (F1-F2)
        A2 = R2 - (F1-F2)
        B1 = X12*F1
        B2 = -X12*F2
        C1 = Y12*F1
        C2 = -Y12*F2
        ZQ = A1*Z(N1) + A2*Z(N2) + B1*ZX1 + B2*ZX2 +&
     &       C1*ZY1 + C2*ZY2
!
!   Compute the extrapolated value at P.
!
        XQP = XP-XQ
        YQP = YP-YQ
        T1 = R1*ZX1 + R2*ZX2
        T2 = R1*ZY1 + R2*ZY2
        PZ = ZQ + T1*XQP + T2*YQP
        IF (DFLAG) THEN
          T = (3.*R12*(2.*(Z(N2)-Z(N1)) - X12*(ZX1+ZX2) -&
     &         Y12*(ZY1+ZY2)) + (ZX2-ZX1)*XQP +&
     &         (ZY2-ZY1)*YQP)/D2
          PZX = T1 + X12*T
          PZY = T2 + Y12*T
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE ROTATE (N,C,S, X,Y )
      INTEGER N
      REAL    C, S, X(N), Y(N)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88
!
!                                                ( C  S)
!   This subroutine applies the Givens rotation  (     )  to
!                                                (-S  C)
!                    (X(1) ... X(N))
! the 2 by N matrix  (             ) .
!                    (Y(1) ... Y(N))
!
!   This routine is identical to Subroutine SROT from the
! LINPACK BLAS (Basic Linear Algebra Subroutines).
!
! On input:
!
!       N = Number of columns to be rotated.
!
!       C,S = Elements of the Givens rotation.  Refer to
!             Subroutine GIVENS.
!
! The above parameters are not altered by this routine.
!
!       X,Y = Arrays of length .GE. N containing the compo-
!             nents of the vectors to be rotated.
!
! On output:
!
!       X,Y = Arrays containing the rotated vectors (not
!             altered if N < 1).
!
! Modules required by ROTATE:  None
!
!***********************************************************
!
      INTEGER I
      REAL    XI, YI
!
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
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88
!
!   This subroutine sets up the I-th row of an augmented re-
! gression matrix for a weighted least squares fit of a
! quadratic function Q(X,Y) to a set of data values Z, where
! Q(XK,YK) = ZK.  The first three columns (quadratic terms)
! are scaled by S2, and the fourth and fifth columns (lin-
! ear terms) are scaled by S1.
!
! On input:
!
!       XK,YK = Coordinates of node K.
!
!       ZK = Data value at node K to be interpolated by Q.
!
!       XI,YI,ZI = Coordinates and data value at node I.
!
!       S1,S2 = Scale factors.
!
!       W = Weight associated with node I.
!
! The above parameters are not altered by this routine.
!
!       ROW = Array of length 6.
!
! On output:
!
!       ROW = Array containing a row of the augmented re-
!             gression matrix.
!
! Modules required by SETRO1:  None
!
!***********************************************************
!
      REAL DX, DY, W1, W2
!
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
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88
!
!   This subroutine sets up the I-th row of an augmented re-
! gression matrix for a weighted least squares fit of a
! quadratic function Q(X,Y) to a set of data values Z.  The
! first 3 columns (quadratic terms) are scaled by S2, and
! the fourth and fifth columns (linear terms) are scaled by
! S1.
!
! On input:
!
!       XK,YK = Coordinates of node K.
!
!       ZK = Data value at node K to be interpolated by Q
!            (Q(XK,YK) = ZK) if the constant term, ROW(6),
!            is to be ignored, or zero if Q(XK,YK) is to be
!            a parameter (coefficient of ROW(6)).
!
!       XI,YI,ZI = Coordinates and data value at node I.
!
!       S1,S2 = Scale factors.
!
!       W = Weight associated with node I.
!
! The above parameters are not altered by this routine.
!
!       ROW = Array of length 7.
!
! On output:
!
!       ROW = Array containing a row of the augmented re-
!             gression matrix.
!
! Modules required by SETRO2:  None
!
!***********************************************************
!
      REAL DX, DY, W1, W2
!
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
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   01/25/97
!
!   This subroutine sets up the I-th row of an augmented re-
! gression matrix for a weighted least squares fit of a
! cubic function f(x,y) to a set of data values z, where
! f(XK,YK) = ZK.  The first four columns (cubic terms) are
! scaled by S3, the next three columns (quadratic terms)
! are scaled by S2, and the eighth and ninth columns (lin-
! ear terms) are scaled by S1.
!
! On input:
!
!       XK,YK = Coordinates of node K.
!
!       ZK = Data value at node K to be interpolated by f.
!
!       XI,YI,ZI = Coordinates and data value at node I.
!
!       S1,S2,S3 = Scale factors.
!
!       W = Weight associated with node I.
!
! The above parameters are not altered by this routine.
!
!       ROW = Array of length 10.
!
! On output:
!
!       ROW = Array containing a row of the augmented re-
!             gression matrix.
!
! Modules required by SETRO3:  None
!
!***********************************************************
!
      REAL DX, DY, W1, W2, W3
!
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
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/03/98
!
!   Given a triangulation of a set of nodes in the plane,
! along with an array of tension factors associated with the
! triangulation arcs, this subroutine prints the list of
! arcs (with tension factors) ordered by endpoint nodal in-
! dexes.  An arc is identified with its smaller endpoint
! index:  N1-N2, where N1 < N2.
!
! On input:
!
!       N = Number of nodes in the triangulation.  3 .LE. N
!           .LE. 9999.
!
!       LUNIT = Logical unit for output.  0 .LE. LUNIT .LE.
!               99.  Output is printed on unit 6 if LUNIT is
!               outside its valid range.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
!       SIGMA = Array of length 2*NA = 6*(N-1)-2*NB, where
!               NA and NB are the numbers of arcs and boun-
!               dary nodes, respectively, containing tension
!               factors associated with arcs in one-to-one
!               correspondence with LIST entries.  Note that
!               each arc N1-N2 has two LIST entries and
!               thus, SIGMA(I) and SIGMA(J) should be iden-
!               tical, where LIST(I) = N2 (in the adjacency
!               list for N1) and LIST(J) = N1 (in the list
!               associated with N2).  Both SIGMA(I) and
!               SIGMA(J) are printed if they are not iden-
!               tical.
!
! None of the parameters are altered by this routine.
!
! TRIPACK module required by SGPRNT:  LSTPTR
!
! Intrinsic function called by SGPRNT:  ABS
!
!***********************************************************
!
      INTEGER LSTPTR
      INTEGER LP1, LP2, LPL, LUN, N1, N2, NA, NAT, NB, NE,&
     &        NL, NLMAX, NM1, NMAX
      LOGICAL ERROR
      REAL    SIG
      DATA NMAX/9999/,  NLMAX/60/
!
      LUN = LUNIT
      IF (LUN .LT. 0  .OR.  LUN .GT. 99) LUN = 6
!
! Print a heading, test for invalid N, and initialize coun-
!   ters:
!
! NL = Number of lines printed on the current page
! NA = Number of arcs encountered
! NE = Number of errors in SIGMA encountered
! NB = Number of boundary nodes encountered
!
      WRITE (LUN,100) N
      IF (N .LT. 3  .OR.  N .GT. NMAX) GO TO 4
      NL = 6
      NA = 0
      NE = 0
      NB = 0
!
! Outer loop on nodes N1.  LPL points to the last neighbor
!   of N1.
!
      NM1 = N - 1
      DO 3 N1 = 1,NM1
        LPL = LEND(N1)
        IF (LIST(LPL) .LT. 0) NB = NB + 1
        LP1 = LPL
!
! Inner loop on neighbors N2 of N1 such that N1 < N2.
!
    1   LP1 = LPTR(LP1)
          N2 = ABS(LIST(LP1))
          IF (N2 .LT. N1) GO TO 2
          NA = NA + 1
          SIG = SIGMA(LP1)
!
!   Test for an invalid SIGMA entry.
!
          LP2 = LSTPTR (LEND(N2),N1,LIST,LPTR)
          ERROR = SIGMA(LP2) .NE. SIG
          IF (ERROR) NE = NE + 1
!
!   Print a line and update the counters.
!
          IF (.NOT. ERROR) WRITE (LUN,110) N1, N2, SIG
          IF (ERROR) WRITE (LUN,120) N1, N2, SIG, SIGMA(LP2)
          NL = NL + 1
          IF (NL .GE. NLMAX) THEN
            WRITE (LUN,130)
            NL = 1
          ENDIF
!
! Bottom of loop on neighbors N2 of N1.
!
    2     IF (LP1 .NE. LPL) GO TO 1
    3   CONTINUE
      LPL = LEND(N)
      IF (LIST(LPL) .LT. 0) NB = NB + 1
!
! Test for errors in SIGMA.
!
      IF (NE .GT. 0) WRITE (LUN,200) NE
!
! Print NA and test for an invalid triangulation.
!
      WRITE (LUN,140) NA
      NAT = 3*NM1 - NB
      IF (NAT .NE. NA) WRITE (LUN,210) NAT
      RETURN
!
! N is outside its valid range.
!
    4 WRITE (LUN,220) NMAX
      RETURN
!
! Print formats:
!
  100 FORMAT (///14X,'Tension Factors,  N =',I5,&
     &        ' Nodes'//1X,18X,'N1',5X,'N2',8X,'Tension'//)
  110 FORMAT (1X,16X,I4,3X,I4,5X,F12.8)
  120 FORMAT (1X,16X,I4,3X,I4,5X,F12.8,3X,F12.8,' *')
  130 FORMAT (///)
  140 FORMAT (//1X,10X,'NA =',I5,' Arcs')
!
! Error messages:
!
  200 FORMAT (//1X,10X,'*',I5,' Errors in SIGMA')
  210 FORMAT (/1X,10X,'*** Error in triangulation:  ',&
     &        '3N-NB-3 = ',I5,' ***')
  220 FORMAT (1X,10X,'*** N is outside its valid range:  ',&
     &        'NMAX = ',I4,' ***')
      END
      REAL FUNCTION SIG0 (N1,N2,N,X,Y,H,LIST,LPTR,LEND,HXHY,&
     &                    IFLGB,HBND,TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGB,&
     &        IFLGS, IER
      REAL    X(N), Y(N), H(N), HXHY(2,N), HBND, TOL,&
     &        SIGMA(*)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/03/98
!
!   Given a triangulation of a set of nodes in the plane,
! along with data values H and gradients (HX,HY) at the
! nodes, this function determines the smallest tension fac-
! tor SIG0 such that the Hermite interpolatory tension
! spline H(T), defined by SIG0 and the endpoint values and
! directional derivatives associated with an arc N1-N2, is
! bounded (either above or below) by HBND for all T in
! (T1,T2), where (T1,T2) denotes an interval corresponding
! to the arc.
!
! On input:
!
!       N1,N2 = Nodal indexes of the endpoints of an arc for
!               which the tension factor is to be computed.
!               The indexes must be distinct and lie in the
!               range 1 to N, and if IFLGS .GE. 1, they must
!               correspond to adjacent nodes in the triangu-
!               lation.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.
!
!       H = Array of length N containing data values at the
!           nodes.  H(I) is associated with (X(I),Y(I)).
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
!       HXHY = Array dimensioned 2 by N whose columns con-
!              partial derivatives at the nodes (X partials
!              in the first row).  Refer to Subroutines
!              GRADC, GRADG, GRADL, and SMSURF.
!
!       IFLGB = Bound option indicator:
!               IFLGB = -1 if HBND is a lower bound on H.
!               IFLGB = 1 if HBND is an upper bound on H.
!
!       HBND = Bound on H.  HBND .LE. min(H1,H2) if IFLGB =
!              -1 and HBND .GE. max(H1,H2) if IFLGB = 1,
!              where H1 and H2 are the data values at the
!              endpoints of the arc N1-N2.
!
!       TOL = Tolerance whose magnitude determines how close
!             SIG0 is to its optimal value when nonzero
!             finite tension is necessary and sufficient to
!             satisfy the constraint.  For a lower bound,
!             SIG0 is chosen so that HBND .LE. HMIN .LE.
!             HBND + abs(TOL), where HMIN is the minimum
!             value of H on the arc, and for an upper bound,
!             the maximum of H satisfies HBND - abs(TOL)
!             .LE. HMAX .LE. HBND.  Thus, the constraint is
!             satisfied but possibly with more tension than
!             necessary.
!
!       IFLGS = Tension array option indicator:
!               IFLGS .LE. 0 if SIGMA is not to be used.
!               IFLGS .GE. 1 if SIGMA is to be updated by
!                            storing SIG0 in the appropriate
!                            locations.
!
! The above parameters are not altered by this function.
!
!       SIGMA = Dummy array of length 1 (IFLGS .LE. 0) or
!               array containing tension factors associated
!               with arcs in one-to-one correspondence with
!               LIST entries (IFLGS .GE. 1).  Refer to Sub-
!               routine GETSIG.
!
! On output:
!
!       SIGMA = Tension factor array updated with the new
!               value if and only if IFLGS .GE. 1 and IER
!               .GE. 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and the
!                     constraint can be satisfied with fin-
!                     ite tension.
!             IER = 1 if no errors were encountered but in-
!                     finite tension is required to satisfy
!                     the constraint (e.g., IFLGB = -1, HBND
!                     = H(T1), and the derivative of H at T1
!                     is negative).
!             IER = -1 if N1, N2, N, or IFLGB is outside its
!                      valid range.
!             IER = -2 if nodes N1 and N2 coincide or IFLGS
!                      .GE. 1 and the nodes are not adja-
!                      cent.
!             IER = -3 if HBND is outside its valid range.
!
!       SIG0 = Minimum tension factor defined above unless
!              IER < 0, in which case SIG0 = -1.  If IER
!              = 1, SIG0 is set to 85, resulting in an
!              approximation to the linear interpolant of
!              the endpoint values.
!
! TRIPACK module required by SIG0:  STORE
!
! SRFPACK module required by SIG0:  SNHCSH
!
! Intrinsic functions called by SIG0:  ABS, EXP, LOG, MAX,
!                                        MIN, REAL, SIGN,
!                                        SQRT
!
!***********************************************************
!
      REAL    STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    A, AA, A0, B, B0, BND, C, C1, C2, COSHM,&
     &        COSHMM, D, D0, D1PD2, D2, DMAX, DSIG, DX, DY,&
     &        E, EMS, F, F0, FMAX, FNEG, FTOL, H1, H2, R,&
     &        RF, RSIG, RTOL, S, S1, S2, SBIG, SCM, SIG,&
     &        SINHM, SNEG, SSINH, SSM, STOL, T, T0, T1, T2,&
     &        TM
      DATA SBIG/85./,  LUN/-1/
      RF = REAL(IFLGB)
      BND = HBND
!
! Print a heading.
!
      IF (LUN .GE. 0) THEN
        IF (RF .LT. 0.) WRITE (LUN,100) N1, N2, BND
        IF (RF .GT. 0.) WRITE (LUN,110) N1, N2, BND
      ENDIF
  100 FORMAT (//1X,'SIG0 -- N1 =',I4,', N2 =',I4,&
     &        ', Lower bound = ',E15.8)
  110 FORMAT (//1X,'SIG0 -- N1 =',I4,', N2 =',I4,&
     &        ', Upper bound = ',E15.8)
!
! Test for errors and store local parameters.
!
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.&
     &    MAX(N1,N2,3) .GT. N  .OR.  ABS(RF) .NE. 1.)&
     &   GO TO 11
      IER = -2
      IF (IFLGS .GT. 0) THEN
!
!   Set LP1 and LP2 to the pointers to N2 as a neighbor of
!     N1 and N1 as a neighbor of N2, respectively.
!
        LPL = LEND(N1)
        LP1 = LPTR(LPL)
    1   IF (LIST(LP1) .EQ. N2) GO TO 2
          LP1 = LPTR(LP1)
          IF (LP1 .NE. LPL) GO TO 1
        IF (ABS(LIST(LP1)) .NE. N2) GO TO 11
!
    2   LPL = LEND(N2)
        LP2 = LPTR(LPL)
    3   IF (LIST(LP2) .EQ. N1) GO TO 4
          LP2 = LPTR(LP2)
          IF (LP2 .NE. LPL) GO TO 3
        IF (ABS(LIST(LP2)) .NE. N1) GO TO 11
      ENDIF
!
! Test for arc length DT = SQRT(DX**2+DY**2) = 0.
!
    4 DX = X(N2) - X(N1)
      DY = Y(N2) - Y(N1)
      IF (DX .EQ. 0.  .AND.  DY .EQ. 0.) GO TO 11
!
! Store endpoint values and test for a valid constraint.
!
      H1 = H(N1)
      H2 = H(N2)
      IER = -3
      IF ((RF .LT. 0.  .AND.  MIN(H1,H2) .LT. BND)  .OR.&
     &    (RF .GT. 0.  .AND.  BND .LT. MAX(H1,H2)))&
     &   GO TO 11
!
! Compute scaled directional derivatives S1,S2 at the end-
!   points (for the direction N1->N2) and test for infinite
!   tension required.
!
      S1 = HXHY(1,N1)*DX + HXHY(2,N1)*DY
      S2 = HXHY(1,N2)*DX + HXHY(2,N2)*DY
      IER = 1
      SIG = SBIG
      IF ((H1 .EQ. BND  .AND.  RF*S1 .GT. 0.)  .OR.&
     &    (H2 .EQ. BND  .AND.  RF*S2 .LT. 0.)) GO TO 10
!
! Test for SIG = 0 sufficient.
!
      IER = 0
      SIG = 0.
      IF (RF*S1 .LE. 0.  .AND.  RF*S2 .GE. 0.) GO TO 10
!
!   Compute first difference S and coefficients A0 and B0
!     of the Hermite cubic interpolant H0(T) = H2 - (S2*R +
!     B0*R**2 + A0*R**3), where R = (T2-T)/DT.
!
      S = H2 - H1
      T0 = 3.*S - S1 - S2
      A0 = 3.*(S-T0)
      B0 = T0 - S2
      D0 = T0*T0 - S1*S2
!
!   H0 has local extrema in (T1,T2) iff S1*S2 < 0 or
!     (T0*(S1+S2) < 0 and D0 .GE. 0).
!
      IF (S1*S2 .GE. 0.  .AND.  (T0*(S1+S2) .GE. 0.  .OR.&
     &    D0 .LT. 0.)) GO TO 10
      IF (A0 .EQ. 0.) THEN
!
!   H0 is quadratic and has an extremum at R = -S2/(2*B0).
!     H0(R) = H2 + S2**2/(4*B0).  Note that A0 = 0 implies
!     2*B0 = S1-S2, and S1*S2 < 0 implies B0 .NE. 0.
!     Also, the extremum is a min iff HBND is a lower bound.
!
        F0 = (BND - H2 - S2*S2/(4.*B0))*RF
      ELSE
!
!   A0 .NE. 0 and H0 has extrema at R = (-B0 +/- SQRT(D0))/
!     A0 = S2/(-B0 -/+ SQRT(D0)), where the negative root
!     corresponds to a min.  The expression for R is chosen
!     to avoid cancellation error.  H0(R) = H2 + (S2*B0 +
!     2*D0*R)/(3*A0).
!
        T = -B0 - SIGN(SQRT(D0),B0)
        R = T/A0
        IF (RF*B0 .GT. 0.) R = S2/T
        F0 = (BND - H2 - (S2*B0+2.*D0*R)/(3.*A0))*RF
      ENDIF
!
!   F0 .GE. 0 iff SIG = 0 is sufficient to satisfy the
!     constraint.
!
      IF (F0 .GE. 0.) GO TO 10
!
! Find a zero of F(SIG) = (BND-H(R))*RF where the derivative
!   of H, HP, vanishes at R.  F is a nondecreasing function,
!   F(0) < 0, and F = FMAX for SIG sufficiently large.
!
! Initialize parameters for the secant method.  The method
!   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
!   where SG0 and SNEG are defined implicitly by DSIG = SIG
!   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
!   SNEG is initialized to a sufficiently large value that
!   FNEG > 0.  This value is used only if the initial value
!   of F is negative.
!
      FMAX = MAX(1.E-3,MIN(ABS(H1-BND),ABS(H2-BND)))
      T = MAX(ABS(H1-BND),ABS(H2-BND))
      SIG = MAX(ABS(S1),ABS(S2))/T
      DMAX = SIG*(1.-T/FMAX)
      SNEG = SIG - DMAX
      IF (LUN .GE. 0) WRITE (LUN,120) SIG, SNEG, F0, FMAX
  120 FORMAT (1X,8X,'SIG = ',E15.8,', SNEG = ',E15.8/&
     &        1X,9X,'F0 = ',E15.8,', FMAX = ',E15.8/)
      DSIG = SIG
      FNEG = FMAX
      D2 = S2 - S
      D1PD2 = S2 - S1
      NIT = 0
!
! Compute an absolute tolerance FTOL = abs(TOL) and a
!   relative tolerance RTOL = 100*Macheps.
!
      FTOL = ABS(TOL)
      RTOL = 1.
    5 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 5
      RTOL = RTOL*200.
!
! Top of loop:  compute F.
!
    6 EMS = EXP(-SIG)
      IF (SIG .LE. .5) THEN
!
!   Use approximations designed to avoid cancellation error
!     (associated with small SIG) in the modified hyperbolic
!     functions.
!
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        AA = A/EMS
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
!
!   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
!     overflow.
!
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
!
!   HP(R) = (S2 - (C1*SINH(SIG*R) - C2*COSHM(SIG*R))/E)/DT
!     = 0 for ESR = (-B +/- SQRT(D))/A = C/(-B -/+ SQRT(D))
!     where ESR = exp(SIG*R), A = C2-C1, D = B**2 - A*C, and
!     B and C are defined below.
!
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
      IF ((RF*B .GT. 0.  .OR.  AA .EQ. 0.)  .AND.&
     &    C/T .GT. 0.) RSIG = LOG(C/T)
      IF ((RSIG .LE. 0.  .OR.  RSIG .GE. SIG)  .AND.&
     &    B .NE. 0.) GO TO 7
!
!   H(R) = H2 - (B*SIG*R + C1 + RF*SQRT(D))/(SIG*E).
!
      F = (BND - H2 + (B*RSIG+C1+RF*T1)/(SIG*E))*RF
!
!   Update the number of iterations NIT.
!
    7 NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,130) NIT, SIG, F
  130 FORMAT (1X,3X,I2,' -- SIG = ',E15.8,', F = ',&
     &        E15.8)
      IF (F0*F .LT. 0.) THEN
!
!   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F and
!     FNEG always have opposite signs.  If SIG is closer to
!     SNEG than SG0, then swap (SNEG,FNEG) with (SG0,F0).
!
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF (ABS(DSIG) .GT. ABS(T1)) THEN
!
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
!
!   Test for convergence.
!
      STOL = RTOL*SIG
      IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.  .AND.&
     &    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 10
!
!   Test for F0 = F = FMAX or F < 0 on the first iteration.
!
      IF (F0 .NE. F  .AND.  (NIT .GT. 1  .OR.  F .GT. 0.))&
     &   GO TO 9
!
!   F*F0 > 0 and either the new estimate would be outside
!     of the bracketing interval of length abs(DMAX) or
!     F < 0 on the first iteration.  Reset (SG0,F0) to
!     (SNEG,FNEG).
!
    8 DSIG = DMAX
      F0 = FNEG
!
!   Compute the change in SIG by linear interpolation
!     between (SG0,F0) and (SIG,F).
!
    9 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,140) DSIG
  140 FORMAT (1X,8X,'DSIG = ',E15.8)
      IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.&
     &     DSIG*DMAX .GT. 0. ) GO TO 8
!
!   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
!     Note that DSIG and DMAX have opposite signs.
!
      IF (ABS(DSIG) .LT. STOL/2.) DSIG = -SIGN(STOL/2.,DMAX)
!
!   Bottom of loop:  Update SIG, DMAX, and F0.
!
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 6
!
! No errors encountered.
!
   10 SIG0 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
!
! Error termination.
!
   11 SIG0 = -1.
      RETURN
      END
      REAL FUNCTION SIG1 (N1,N2,N,X,Y,H,LIST,LPTR,LEND,HXHY,&
     &                    IFLGB,HPBND,TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGB,&
     &        IFLGS, IER
      REAL    X(N), Y(N), H(N), HXHY(2,N), HPBND, TOL,&
     &        SIGMA(*)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/03/98
!
!   Given a triangulation of a set of nodes in the plane,
! along with data values H and gradients (HX,HY) at the
! nodes, this function determines the smallest tension fac-
! tor SIG1 such that the first derivative HP(T) of the
! Hermite interpolatory tension spline H(T), defined by SIG1
! and the endpoint values and directional derivatives asso-
! ciated with an arc N1-N2, is bounded (either above or
! below) by HPBND for all T in (T1,T2), where (T1,T2) de-
! notes an interval corresponding to the arc.
!
! On input:
!
!       N1,N2 = Nodal indexes of the endpoints of an arc for
!               which the tension factor is to be computed.
!               The indexes must be distinct and lie in the
!               range 1 to N, and if IFLGS .GE. 1, they must
!               correspond to adjacent nodes in the triangu-
!               lation.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.
!
!       H = Array of length N containing data values at the
!           nodes.  H(I) is associated with (X(I),Y(I)).
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
!       HXHY = Array dimensioned 2 by N whose columns con-
!              partial derivatives at the nodes (X partials
!              in the first row).  Refer to Subroutines
!              GRADC, GRADG, GRADL, and SMSURF.
!
!       IFLGB = Bound option indicator:
!               IFLGB = -1 if HPBND is a lower bound on HP.
!               IFLGB = 1 if HPBND is an upper bound on HP.
!
!       HPBND = Bound on HP.  HPBND .LE. min(HP1,HP2,S) if
!               IFLGB = -1 and HPBND .GE. max(HP1,HP2,S) if
!               IFLGB = 1, where HP1 and HP2 are the direc-
!               tional derivatives at the endpoints of the
!               arc N1-N2, and S is the slope of the linear
!               interpolant of the endpoint data values.
!
!       TOL = Tolerance whose magnitude determines how close
!             SIG1 is to its optimal value when nonzero
!             finite tension is necessary and sufficient to
!             satisfy the constraint.  For a lower bound,
!             SIG1 is chosen so that HPBND .LE. HPMIN .LE.
!             HPBND + abs(TOL), where HPMIN is the minimum
!             value of HP on the arc.  For an upper bound,
!             the maximum of HP satisfies HPBND - abs(TOL)
!             .LE. HPMAX .LE. HPBND.  Thus, the constraint
!             is satisfied but possibly with more tension
!             than necessary.
!
!       IFLGS = Tension array option indicator:
!               IFLGS .LE. 0 if SIGMA is not to be used.
!               IFLGS .GE. 1 if SIGMA is to be updated by
!                            storing SIG1 in the appropriate
!                            locations.
!
! The above parameters are not altered by this function.
!
!       SIGMA = Dummy array of length 1 (IFLGS .LE. 0) or
!               array containing tension factors associated
!               with arcs in one-to-one correspondence with
!               LIST entries (IFLGS .GE. 1).  Refer to Sub-
!               routine GETSIG.
!
! On output:
!
!       SIGMA = Tension factor array updated with the new
!               value if and only if IFLGS .GE. 1 and IER
!               .GE. 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and the
!                     constraint can be satisfied with fin-
!                     ite tension.
!             IER = 1 if no errors were encountered but in-
!                     finite tension is required to satisfy
!                     the constraint (e.g., IFLGB = -1,
!                     HPBND = S, and HP1 > S).
!             IER = -1 if N1, N2, N, or IFLGB is outside its
!                      valid range.
!             IER = -2 if nodes N1 and N2 coincide or IFLGS
!                      .GE. 1 and the nodes are not adja-
!                      cent.
!             IER = -3 if HPBND is outside its valid range.
!
!       SIG1 = Minimum tension factor defined above unless
!              IER < 0, in which case SIG1 = -1.  If IER
!              = 1, SIG1 is set to 85, resulting in an
!              approximation to the linear interpolant of
!              the endpoint values.
!
! TRIPACK module required by SIG1:  STORE
!
! SRFPACK module required by SIG1:  SNHCSH
!
! Intrinsic functions called by SIG1:  ABS, EXP, MAX, MIN,
!                                        REAL, SIGN, SQRT
!
!***********************************************************
!
      REAL    STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    A, A0, B0, BND, C0, C1, C2, COSHM, COSHMM, D0,&
     &        D1, D1PD2, D2, DMAX, DSIG, DT, DX, DY, E, EMS,&
     &        EMS2, F, F0, FMAX, FNEG, FTOL, RF, RTOL, S,&
     &        S1, S2, SBIG, SIG, SINH, SINHM, STOL, T0, T1,&
     &        T2, TM
      DATA SBIG/85./,  LUN/-1/
      RF = REAL(IFLGB)
      BND = HPBND
!
! Print a heading.
!
      IF (LUN .GE. 0) THEN
        IF (RF .LT. 0.) WRITE (LUN,100) N1, N2, BND
        IF (RF .GT. 0.) WRITE (LUN,110) N1, N2, BND
      ENDIF
  100 FORMAT (//1X,'SIG1 -- N1 =',I4,', N2 =',I4,&
     &        ', Lower bound = ',E15.8)
  110 FORMAT (//1X,'SIG1 -- N1 =',I4,', N2 =',I4,&
     &        ', Upper bound = ',E15.8)
!
! Test for errors and store local parameters.
!
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.&
     &    MAX(N1,N2,3) .GT. N  .OR.  ABS(RF) .NE. 1.)&
     &   GO TO 10
      IER = -2
      IF (IFLGS .GT. 0) THEN
!
!   Set LP1 and LP2 to the pointers to N2 as a neighbor of
!     N1 and N1 as a neighbor of N2, respectively.
!
        LPL = LEND(N1)
        LP1 = LPTR(LPL)
    1   IF (LIST(LP1) .EQ. N2) GO TO 2
          LP1 = LPTR(LP1)
          IF (LP1 .NE. LPL) GO TO 1
        IF (ABS(LIST(LP1)) .NE. N2) GO TO 10
!
    2   LPL = LEND(N2)
        LP2 = LPTR(LPL)
    3   IF (LIST(LP2) .EQ. N1) GO TO 4
          LP2 = LPTR(LP2)
          IF (LP2 .NE. LPL) GO TO 3
        IF (ABS(LIST(LP2)) .NE. N1) GO TO 10
      ENDIF
!
! Test for arc length DT = SQRT(DX**2+DY**2) = 0.
!
    4 DX = X(N2) - X(N1)
      DY = Y(N2) - Y(N1)
      IF (DX .EQ. 0.  .AND.  DY .EQ. 0.) GO TO 10
!
! Compute first difference S and scaled directional deriva-
!   tives S1,S2 at the endpoints (for the direction N1->N2).
!
      S = H(N2) - H(N1)
      S1 = HXHY(1,N1)*DX + HXHY(2,N1)*DY
      S2 = HXHY(1,N2)*DX + HXHY(2,N2)*DY
!
! Test for a valid constraint.
!
      IER = -3
      IF ((RF .LT. 0.  .AND.  MIN(S1,S2,S) .LT. BND)  .OR.&
     &    (RF .GT. 0.  .AND.  BND .LT. MAX(S1,S2,S)))&
     &   GO TO 10
!
! Test for infinite tension required.
!
      IER = 1
      SIG = SBIG
      IF (S .EQ. BND  .AND.  (S1 .NE. S  .OR.  S2 .NE. S))&
     &   GO TO 9
!
! Test for SIG = 0 sufficient.  The Hermite cubic interpo-
!   lant H0 has derivative HP0(T) = (S2 + 2*B0*R + A0*R**2)/
!   DT, where R = (T2-T)/DT.
!
      IER = 0
      SIG = 0.
      T0 = 3.*S - S1 - S2
      B0 = T0 - S2
      C0 = T0 - S1
      A0 = -B0 - C0
!
!   HP0(R) has an extremum (at R = -B0/A0) in (0,1) iff
!     B0*C0 > 0 and the third derivative of H0 has the
!     sign of A0.
!
      IF (B0*C0 .LE. 0.  .OR.  A0*RF .GT. 0.) GO TO 9
!
!   A0*RF < 0 and HP0(R) = -D0/(DT*A0) at R = -B0/A0.
!
      DT = SQRT(DX*DX + DY*DY)
      D0 = T0*T0 - S1*S2
      F0 = (BND + D0/(A0*DT))*RF
      IF (F0 .GE. 0.) GO TO 9
!
! Find a zero of F(SIG) = (BND-HP(R))*RF, where HP has an
!   extremum at R.  F has a unique zero, F(0) = F0 < 0, and
!   F = (BND-S)*RF > 0 for SIG sufficiently large.
!
! Initialize parameters for the secant method.  The method
!   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
!   where SG0 and SNEG are defined implicitly by DSIG = SIG
!   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
!   SIG is initialized to the zero of (BND - (SIG*S-S1-S2)/
!   (DT*(SIG-2.)))*RF -- a value for which F(SIG) .GE. 0 and
!   F(SIG) = 0 for SIG sufficiently large that 2*SIG is in-
!   significant relative to exp(SIG).
!
      FMAX = (BND-S/DT)*RF
      SIG = 2. - A0/(3.*(DT*BND-S))
      IF (LUN .GE. 0) WRITE (LUN,120) F0, FMAX, SIG
  120 FORMAT (1X,9X,'F0 = ',E15.8,', FMAX = ',E15.8/&
     &        1X,8X,'SIG = ',E15.8/)
      IF (STORE(SIG*EXP(-SIG)+.5) .EQ. .5) GO TO 9
      DSIG = SIG
      DMAX = -2.*SIG
      FNEG = FMAX
      D1 = S - S1
      D2 = S2 - S
      D1PD2 = D1 + D2
      NIT = 0
!
! Compute an absolute tolerance FTOL = abs(TOL), and a
!   relative tolerance RTOL = 100*Macheps.
!
      FTOL = ABS(TOL)
      RTOL = 1.
    5 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 5
      RTOL = RTOL*200.
!
! Top of loop:  compute F.
!
    6 IF (SIG .LE. .5) THEN
!
!   Use approximations designed to avoid cancellation
!     error (associated with small SIG) in the modified
!     hyperbolic functions.
!
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
!
!   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
!     overflow.
!
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
!
!   The second derivative HPP of H(R) has a zero at exp(SIG*
!     R) = SQRT((C2+C1)/A) and R is in (0,1) and well-
!     defined iff HPP(T1)*HPP(T2) < 0.
!
      F = FMAX
      T1 = A*(C2+C1)
      IF (T1 .GE. 0.) THEN
        IF (C1*(SIG*COSHM*D1 - SINHM*D1PD2) .LT. 0.) THEN
!
!   HP(R) = (B+SIGN(A)*SQRT(A*C))/(DT*E) at the critical
!     value of R, where A = C2-C1, B = E*S2-C2, and C = C2 +
!     C1.  Note that RF*A < 0.
!
          F = (BND - (E*S2-C2 - RF*SQRT(T1))/(DT*E))*RF
        ENDIF
      ENDIF
!
!   Update the number of iterations NIT.
!
      NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,130) NIT, SIG, F
  130 FORMAT (1X,3X,I2,' -- SIG = ',E15.8,', F = ',&
     &        E15.8)
      IF (F0*F .LT. 0.) THEN
!
!   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F
!     and FNEG always have opposite signs.  If SIG is closer
!     to SNEG than SG0 and abs(F) < abs(FNEG), then swap
!     (SNEG,FNEG) with (SG0,F0).
!
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF ( ABS(DSIG) .GT. ABS(T1)  .AND.&
     &       ABS(F) .LT. ABS(T2) ) THEN
!
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
!
!   Test for convergence.
!
      STOL = RTOL*SIG
      IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.  .AND.&
     &    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 9
      IF (F0*F .LT. 0.  .OR.  ABS(F) .LT. ABS(F0)) GO TO 8
!
!   F*F0 > 0 and the new estimate would be outside of the
!     bracketing interval of length abs(DMAX).  Reset
!     (SG0,F0) to (SNEG,FNEG).
!
    7 DSIG = DMAX
      F0 = FNEG
!
!   Compute the change in SIG by linear interpolation
!     between (SG0,F0) and (SIG,F).
!
    8 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,140) DSIG
  140 FORMAT (1X,8X,'DSIG = ',E15.8)
      IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.&
     &     DSIG*DMAX .GT. 0. ) GO TO 7
!
!   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
!     Note that DSIG and DMAX have opposite signs.
!
      IF (ABS(DSIG) .LT. STOL/2.) DSIG = -SIGN(STOL/2.,DMAX)
!
!   Bottom of loop:  update SIG, DMAX, and F0.
!
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 6
!
! No errors encountered.
!
    9 SIG1 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
!
! Error termination.
!
   10 SIG1 = -1.
      RETURN
      END
      REAL FUNCTION SIG2 (N1,N2,N,X,Y,H,LIST,LPTR,LEND,HXHY,&
     &                    TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGS,&
     &        IER
      REAL    X(N), Y(N), H(N), HXHY(2,N), TOL, SIGMA(*)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/22/97
!
!   Given a triangulation of a set of nodes in the plane,
! along with data values H and gradients (HX,HY) at the
! nodes, this function determines the smallest tension fac-
! tor SIG2 such that the Hermite interpolatory tension
! spline H(T), defined by SIG2 and the endpoint values and
! directional derivatives associated with an arc N1-N2,
! preserves convexity (or concavity) of the data:
!
!   HP1 .LE. S .LE. HP2 implies HPP(T) .GE. 0, and
!   HP1 .GE. S .GE. HP2 implies HPP(T) .LE. 0
!
! for all T in the open interval (T1,T2) corresponding to
! the arc, where HP1 and HP2 are the derivative values of H
! at the endpoints, S is the slope of the linear interpolant
! of the endpoint data values, and HPP denotes the second
! derivative of H.  Note, however, that infinite tension is
! required if HP1 = S or HP2 = S (unless HP1 = HP2 = S).
!
! On input:
!
!       N1,N2 = Nodal indexes of the endpoints of an arc for
!               which the tension factor is to be computed.
!               The indexes must be distinct and lie in the
!               range 1 to N, and if IFLGS .GE. 1, they must
!               correspond to adjacent nodes in the triangu-
!               lation.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.
!
!       H = Array of length N containing data values at the
!           nodes.  H(I) is associated with (X(I),Y(I)).
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
!       HXHY = Array dimensioned 2 by N whose columns con-
!              partial derivatives at the nodes (X partials
!              in the first row).  Refer to Subroutines
!              GRADC, GRADG, GRADL, and SMSURF.
!
!       TOL = Tolerance whose magnitude determines how close
!             SIG2 is to its optimal value when nonzero
!             finite tension is necessary and sufficient to
!             satisfy convexity or concavity.  In the case
!             convexity, SIG2 is chosen so that 0 .LE.
!             HPPMIN .LE. abs(TOL), where HPPMIN is the
!             minimum value of HPP on the arc.  In the case
!             of concavity, the maximum value of HPP satis-
!             fies -abs(TOL) .LE. HPPMAX .LE. 0.  Thus, the
!             constraint is satisfied but possibly with more
!             tension than necessary.
!
!       IFLGS = Tension array option indicator:
!               IFLGS .LE. 0 if SIGMA is not to be used.
!               IFLGS .GE. 1 if SIGMA is to be updated by
!                            storing SIG2 in the appropriate
!                            locations.
!
! The above parameters are not altered by this function.
!
!       SIGMA = Dummy array of length 1 (IFLGS .LE. 0) or
!               array containing tension factors associated
!               with arcs in one-to-one correspondence with
!               LIST entries (IFLGS .GE. 1).  Refer to Sub-
!               routine GETSIG.
!
! On output:
!
!       SIGMA = Tension factor array updated with the new
!               value if and only if IFLGS .GE. 1 and IER
!               .GE. 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and fin-
!                     ite tension is sufficient to satisfy
!                     convexity (or concavity).
!             IER = 1 if no errors were encountered but in-
!                     finite tension is required to satisfy
!                     convexity.
!             IER = 2 if the data does not satisfy convexity
!                     or concavity.
!             IER = -1 if N1, N2, or N is outside its valid
!                      range.
!             IER = -2 if nodes N1 and N2 coincide or IFLGS
!                      .GE. 1 and the nodes are not adja-
!                      cent.
!
!       SIG2 = Minimum tension factor defined above unless
!              IER < 0, in which case SIG2 = -1.  If IER
!              = 1, SIG2 is set to 85, resulting in an
!              approximation to the linear interpolant of
!              the endpoint values.  If IER = 2, SIG2 = 0,
!              resulting in the Hermite cubic interpolant.
!
! TRIPACK module required by SIG2:  STORE
!
! SRFPACK module required by SIG2:  SNHCSH
!
! Intrinsic functions called by SIG2:  ABS, EXP, MAX, MIN,
!                                        SQRT
!
!***********************************************************
!
      REAL    STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    COSHM, D1, D1D2, D2, DSIG, DUMMY, DX, DY, EMS,&
     &        F, FP, FTOL, RTOL, S, SBIG, SIG, SINHM, SSM,&
     &        T, T1, TP1
      DATA SBIG/85./,  LUN/-1/
!
! Print a heading.
!
      IF (LUN .GE. 0) WRITE (LUN,100) N1, N2
  100 FORMAT (//1X,'SIG2 -- N1 =',I4,', N2 =',I4)
!
! Test for errors and store local parameters.
!
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.&
     &    MAX(N1,N2,3) .GT. N) GO TO 8
      IER = -2
      IF (IFLGS .GT. 0) THEN
!
!   Set LP1 and LP2 to the pointers to N2 as a neighbor of
!     N1 and N1 as a neighbor of N2, respectively.
!
        LPL = LEND(N1)
        LP1 = LPTR(LPL)
    1   IF (LIST(LP1) .EQ. N2) GO TO 2
          LP1 = LPTR(LP1)
          IF (LP1 .NE. LPL) GO TO 1
        IF (ABS(LIST(LP1)) .NE. N2) GO TO 8
!
    2   LPL = LEND(N2)
        LP2 = LPTR(LPL)
    3   IF (LIST(LP2) .EQ. N1) GO TO 4
          LP2 = LPTR(LP2)
          IF (LP2 .NE. LPL) GO TO 3
        IF (ABS(LIST(LP2)) .NE. N1) GO TO 8
      ENDIF
!
! Test for arc length DT = SQRT(DX**2+DY**2) = 0.
!
    4 DX = X(N2) - X(N1)
      DY = Y(N2) - Y(N1)
      IF (DX .EQ. 0.  .AND.  DY .EQ. 0.) GO TO 8
!
! Compute first and second differences and test for infinite
!   tension required.
!
      S = H(N2) - H(N1)
      D1 = S - HXHY(1,N1)*DX - HXHY(2,N1)*DY
      D2 = HXHY(1,N2)*DX + HXHY(2,N2)*DY - S
      D1D2 = D1*D2
      IER = 1
      SIG = SBIG
      IF (D1D2 .EQ. 0.  .AND.  D1 .NE. D2) GO TO 7
!
! Test for a valid constraint.
!
      IER = 2
      SIG = 0.
      IF (D1D2 .LT. 0.) GO TO 7
!
! Test for SIG = 0 sufficient.
!
      IER = 0
      IF (D1D2 .EQ. 0.) GO TO 7
      T = MAX(D1/D2,D2/D1)
      IF (T .LE. 2.) GO TO 7
!
! Find a zero of F(SIG) = SIG*COSHM(SIG)/SINHM(SIG) - (T+1).
!   Since the derivative of F vanishes at the origin, a
!   quadratic approximation is used to obtain an initial
!   estimate for the Newton method.
!
      TP1 = T + 1.
      SIG = SQRT(10.*T-20.)
      NIT = 0
!
!   Compute an absolute tolerance FTOL = abs(TOL) and a
!     relative tolerance RTOL = 100*Macheps.
!
      FTOL = ABS(TOL)
      RTOL = 1.
    5 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 5
      RTOL = RTOL*200.
!
! Top of loop:  evaluate F and its derivative FP.
!
    6 IF (SIG .LE. .5) THEN
!
!   Use approximations designed to avoid cancellation error
!     in the hyperbolic functions.
!
        CALL SNHCSH (SIG, SINHM,COSHM,DUMMY)
        T1 = COSHM/SINHM
        FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.)
      ELSE
!
!   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
!     overflow.
!
        EMS = EXP(-SIG)
        SSM = 1. - EMS*(EMS+SIG+SIG)
        T1 = (1.-EMS)*(1.-EMS)/SSM
        FP = T1 + SIG*(2.*SIG*EMS/SSM - T1*T1 + 1.)
      ENDIF
!
      F = SIG*T1 - TP1
!
!   Update the number of iterations NIT.
!
      NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,110) NIT, SIG, F, FP
  110 FORMAT (1X,3X,I2,' -- SIG = ',E15.8,', F = ',&
     &        E15.8/1X,31X,'FP = ',E15.8)
!
!   Test for convergence.
!
      IF (FP .LE. 0.) GO TO 7
      DSIG = -F/FP
      IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.  .AND.&
     &    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 7
!
!   Bottom of loop:  update SIG.
!
      SIG = SIG + DSIG
      GO TO 6
!
! No errors encountered.
!
    7 SIG2 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
!
! Error termination.
!
    8 SIG2 = -1.
      RETURN
      END
      SUBROUTINE SMSGS (NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,&
     &                  IFLGS,SIGMA,W,P, NIT,DFMAX,F,&
     &                  FXFY, IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),&
     &        IFLGS, NIT, IER
      REAL    X(N), Y(N), Z(N), SIGMA(*), W(N), P, DFMAX,&
     &        F(N), FXFY(2,N)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   08/26/91
!
!   This subroutine employs the block Gauss-Seidel method
! (3 by 3 blocks) to solve the order 3N symmetric positive
! definite linear system associated with minimizing the
! quadratic functional Q(F,FX,FY) described in Subroutine
! SMSURF.
!
!   Note that small relative changes in F can cause large
! relative changes in FX and FY, resulting in an ill-
! conditioned system.  However, good F values should be
! achieved with a small number of iterations, and the
! gradients (with fixed F) can then be improved by a call
! to Subroutine GRADG.
!
! On input:
!
!   NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,IFLGS,SIGMA,W =
! Parameters defined in Subroutine SMSURF.
!
!       P = Positive smoothing parameter defining Q.
!
! The above parameters are not altered by this routine.
!
!       NIT = Maximum number of Gauss-Seidel iterations to
!             be employed.  This maximum will likely be
!             achieved if DFMAX is smaller than the machine
!             precision.  NIT .GE. 0.
!
!       DFMAX = Nonnegative convergence criterion.  The
!               method is terminated when the maximum
!               change in a solution F-component between
!               iterations is at most DFMAX.  The change in
!               a component is taken to be the absolute
!               difference relative to 1 plus the old value.
!
!       F = Initial estimate of the first N solution compo-
!           nents.
!
!       FXFY = 2 by N array containing initial estimates of
!              the last 2N solution components.
!
! On output:
!
!       NIT = Number of Gauss-Seidel iterations employed.
!
!       DFMAX = Maximum relative change in a solution F-
!               component at the last iteration.
!
!       F = First N solution components -- function values
!           at the nodes.
!
!       FXFY = Last 2N solution components -- gradients at
!              the nodes with X partials in the first row.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and the
!                     convergence criterion was achieved.
!             IER = 1 if no errors were encountered but con-
!                     vergence was not achieved within NIT
!                     iterations.
!             IER = -1 if NCC, N, P, NIT, or DFMAX is out-
!                      side its valid range on input.  F
!                      and FXFY are not altered in this
!                      case.  LCC is not tested for valid-
!                      ity.
!             IER = -2 if all nodes are collinear or the
!                      triangulation data structure is not
!                      valid.
!             IER = -3 if duplicate nodes were encountered.
!
! SRFPACK modules required by SMSGS:  GRCOEF, SNHCSH
!
! Intrinsic functions called by SMSGS:  ABS, MAX, SQRT
!
!***********************************************************
!
      INTEGER I, IFL, IFRST, ILAST, ITER, ITMAX, J, K, KBAK,&
     &        KFOR, LCC1, LP, LPJ, LPL, LPLJ, NB, NN
      REAL    C11, C12, C13, C22, C23, C33, CC22, CC23,&
     &        CC33, DCUB, DET, DF, DFMX, DFX, DFY, DSQ, DX,&
     &        DXS, DXDY, DY, DYS, FK, FXJ, FXK, FYJ, FYK,&
     &        PP, R1, R2, R3, RR2, RR3, SIG, T, T1, T2, T3,&
     &        TOL, TRMX, TRMY, XK, YK
!
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
!
! Test for errors in input and initialize the iteration
!   count ITER, tension factor SIG, and output value of
!   DFMAX.
!
      IF (NCC .LT. 0  .OR.  NN .LT. 3  .OR.  PP .LE. 0.&
         &.OR.  ITMAX .LT. 0  .OR.  TOL .LT. 0.) GO TO 8
      ITER = 0
      SIG = SIGMA(1)
      DFMX = 0.
!
! Top of iteration loop:  If K is a constraint node, I
!   indexes the constraint containing node K, IFRST and
!   ILAST are the first and last nodes of constraint I,
!   and (KBAK,K,KFOR) is a subsequence of constraint I.
!
    1 IF (ITER .EQ. ITMAX) GO TO 7
      DFMX = 0.
      I = 0
      ILAST = LCC1-1
      KBAK = 0
      KFOR = 0
!
!   Loop on nodes.
!
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
!
!   Initialize components of the order 3 system for the
!     change (DF,DFX,DFY) in the K-th solution components
!     (symmetric matrix in C and residual in R).
!
        C11 = PP*W(K)
        C12 = 0.
        C13 = 0.
        C22 = 0.
        C23 = 0.
        C33 = 0.
        R1 = C11*(Z(K)-FK)
        R2 = 0.
        R3 = 0.
!
!   Loop on neighbors J of node K.
!
        LPL = LEND(K)
        LPJ = LPL
    2   LPJ = LPTR(LPJ)
          J = ABS(LIST(LPJ))
!
!   Arc K-J lies in a constraint region and is bypassed iff
!     K and J are nodes in the same constraint and J follows
!     KFOR and precedes KBAK as a neighbor of K.  Also, K-J
!     is bypassed if it is both a constraint arc and a
!     boundary arc of the triangulation.
!
          IF (K .LT. LCC1  .OR.  J .LT. IFRST  .OR.&
     &        J .GT. ILAST) GO TO 4
          IF (J .EQ. KBAK  .OR.  J .EQ. KFOR) THEN
            LPLJ = LEND(J)
            IF (LIST(LPL) .EQ. -J  .OR.&
     &          LIST(LPLJ) .EQ. -K) THEN
              GO TO 5
            ELSE
              GO TO 4
            ENDIF
          ENDIF
          LP = LPJ
!
    3     LP = LPTR(LP)
            NB = ABS(LIST(LP))
            IF (NB .EQ. KBAK) GO TO 5
            IF (NB .NE. KFOR) GO TO 3
!
!   Compute parameters associated with edge K->J, and test
!     for duplicate nodes.
!
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
!
!   T1 = SIG*SIG*COSHM/(DCUB*E), T2 = SIG*SINHM/(DCUB*E),
!     and T3 = SIG*(SIG*COSHM-SINHM)/(DCUB*E) for E =
!     SIG*SINH - 2*COSHM.
!
          T = T1*(FK-F(J))
          FXJ = FXFY(1,J)
          FYJ = FXFY(2,J)
!
!   Update the system components for node J.
!
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
!
!   Bottom of loop on neighbors.
!
    5     IF (LPJ .NE. LPL) GO TO 2
!
!   Solve the system associated with the K-th block.
!
        CC22 = C11*C22 - C12*C12
        CC23 = C11*C23 - C12*C13
        CC33 = C11*C33 - C13*C13
        RR2 = C11*R2 - C12*R1
        RR3 = C11*R3 - C13*R1
        DET = CC22*CC33 - CC23*CC23
        IF (DET .EQ. 0.  .OR.  CC22 .EQ. 0.  .OR.&
     &      C11 .EQ. 0.) GO TO 9
        DFY = (CC22*RR3 - CC23*RR2)/DET
        DFX = (RR2 - CC23*DFY)/CC22
        DF = (R1 - C12*DFX - C13*DFY)/C11
!
!   Update the solution components for node K and the
!     maximum relative change in F.
!
        F(K) = FK + DF
        FXFY(1,K) = FXK + DFX
        FXFY(2,K) = FYK + DFY
        DFMX = MAX(DFMX,ABS(DF)/(1.+ABS(FK)))
    6   CONTINUE
!
!   Increment ITER and test for convergence.
!
      ITER = ITER + 1
      IF (DFMX .GT. TOL) GO TO 1
!
! Method converged.
!
      NIT = ITER
      DFMAX = DFMX
      IER = 0
      RETURN
!
! Method failed to converge within NIT iterations.
!
    7 DFMAX = DFMX
      IER = 1
      RETURN
!
! Invalid input parameter.
!
    8 NIT = 0
      DFMAX = 0.
      IER = -1
      RETURN
!
! Node K and its neighbors are collinear, resulting in a
!   singular system.
!
    9 NIT = 0
      DFMAX = DFMX
      IER = -2
      RETURN
!
! Nodes J and K coincide.
!
   10 NIT = 0
      DFMAX = DFMX
      IER = -3
      RETURN
      END
      SUBROUTINE SMSURF (NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND,&
     &                   IFLGS,SIGMA,W,SM,SMTOL,GSTOL, F,&
     &                   FXFY,IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),&
     &        IFLGS, IER
      REAL    X(N), Y(N), Z(N), SIGMA(*), W(N), SM, SMTOL,&
     &        GSTOL, F(N), FXFY(2,N)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/03/98
!
!   Given a triangulation of N nodes in the plane, along
! with data values Z at the nodes and tension factors SIGMA
! associated with the arcs, this subroutine determines a
! set of nodal function values F and gradients (FX,FY) such
! that a quadratic functional Q1(F,FX,FY) is minimized sub-
! ject to the constraint Q2(F) .LE. SM for Q2(F) = (Z-F)**T*
! W*(Z-F), where W is a diagonal matrix of positive weights.
! The functional Q1 is an approximation to the linearized
! curvature over the triangulation of a C-1 bivariate func-
! tion F(X,Y) which interpolates the nodal values and
! gradients.  Subroutines INTRC1 and UNIF may be called to
! evaluate F at arbitrary points.
!
!   The smoothing procedure is an extension of the method
! for cubic spline smoothing due to C. Reinsch -- Numer.
! Math., 10 (1967) and 16 (1971).  Refer to Subroutines FVAL
! and TVAL for a further description of the interpolant F.
! Letting D1F(T) and D2F(T) denote first and second deriva-
! tives of F with respect to a parameter T varying along a
! triangulation arc, Q1 is the sum over the triangulation
! arcs, excluding interior constraint arcs, of the integrals
! of
!
!      D2F(T)**2 + [(SIGMA/L)*(D1F(T)-S)]**2 ,
!
! where L denotes arc length, SIGMA is the appropriate ten-
! sion factor, and S is the slope of the linear function of
! T which interpolates the values of F at the endpoints of
! the arc.  Introducing a smoothing parameter P, and assum-
! ing the constraint is active, the problem is equivalent to
! minimizing
!
!      Q(P,F,FX,FY) = Q1(F,FX,FY) + P*(Q2(F)-SM) .
!
! The secant method is used to find a zero of
!
!      G(P) = 1/SQRT(Q2) - 1/SQRT(SM) ,
!
! where F(P) satisfies the order 3N symmetric positive def-
! inite linear system obtained by setting the gradient of Q
! (treated as a function of F, FX, and FY) to zero.  The
! linear system is solved by the Gauss-Seidel method.
!
!   Note that the method can also be used to select grad-
! ients for the interpolation problem (F = Z, SM = 0, and P
! infinite).  This is achieved by a call to Subroutine
! GRADG.
!
! On input:
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC .GE. 0.
!
!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index of the
!             first node of constraint I in LCC(I).  For I =
!             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       Z = Array of length N containing data values at the
!           nodes.  Z(I) is associated with (X(I),Y(I)).
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
!       IFLGS = Tension factor option:
!               IFLGS .LE. 0 if a single uniform tension
!                            factor is to be used.
!               IFLGS .GE. 1 if variable tension is desired.
!
!       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
!               array containing tension factors associated
!               with arcs in one-to-one correspondence with
!               LIST entries (IFLGS .GE. 1).  Refer to Sub-
!               routines GETSIG, SIG0, SIG1, and SIG2.
!
!       W = Array of length N containing positive weights
!           associated with the data values.  The recommend-
!           ed value of W(I) is 1/DZ**2, where DZ is the
!           standard deviation associated with Z(I).  DZ**2
!           is the expected value of the squared error in
!           the measurement of Z(I).  (The mean error is
!           assumed to be zero.)
!
!       SM = Positive parameter specifying an upper bound on
!            Q2(F).  Note that F(X,Y) is linear (and Q2(F)
!            is minimized) if SM is sufficiently large that
!            the constraint is not active.  It is recommend-
!            ed that SM satisfy N-SQRT(2N) .LE. SM .LE. N+
!            SQRT(2N).
!
!       SMTOL = Parameter in the open interval (0,1) speci-
!               fying the relative error allowed in satisfy-
!               ing the constraint -- the constraint is
!               assumed to be satisfied if SM*(1-SMTOL) .LE.
!               Q2 .LE. SM*(1+SMTOL).  A reasonable value
!               for SMTOL is SQRT(2/N).
!
!       GSTOL = Nonnegative tolerance defining the conver-
!               gence criterion for the Gauss-Seidel method.
!               Refer to parameter DFMAX in Subroutine
!               SMSGS.  A recommended value is .05*DU**2,
!               where DU is an average standard deviation
!               in the data values.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       F = Array of length N containing nodal function val-
!           ues unless IER < 0.
!
!       FXFY = 2 by N array whose columns contain partial
!              derivatives of F at the nodes unless IER < 0,
!              with FX in the first row, FY in the second.
!
!       IER = Error indicator and information flag:
!             IER = 0 if no errors were encountered and the
!                     constraint is active -- Q2(F) is ap-
!                     proximately equal to SM.
!             IER = 1 if no errors were encountered but the
!                     constraint is not active -- F, FX, and
!                     FY are the values and partials of a
!                     linear function which minimizes Q2(F),
!                     and Q1 = 0.
!             IER = -1 if NCC, an LCC entry, N, W, SM,
!                      SMTOL, or GSTOL is outside its
!                      valid range on input.
!             IER = -2 if all nodes are collinear or the
!                      triangulation data structure is not
!                      valid.
!             IER = -3 if duplicate nodes were encountered.
!
! SRFPACK modules required by SMSURF:  GRCOEF, SMSGS, SNHCSH
!
! Intrinsic functions called by SMSURF:  ABS, SQRT
!
!***********************************************************
!
      INTEGER I, IERR, ITER, LCCIP1, LUN, NIT, NITMAX, NN
      REAL    C11, C12, C13, C22, C23, C33, CC22, CC23,&
     &        CC33, DET, DFMAX, DMAX, DP, F0, FX, FY, G, G0,&
     &        GNEG, P, Q2, Q2MAX, Q2MIN, R1, R2, R3, RR2,&
     &        RR3, S, TOL, WI, WIXI, WIYI, WIZI, XI, YI
!
      DATA NITMAX/40/,  LUN/-1/
!
! LUN = Logical unit on which diagnostic messages are print-
!       ed (unless LUN < 0).  For each secant iteration, the
!       following values are printed:  P, G(P), NIT, DFMAX,
!       and DP, where NIT denotes the number of Gauss-Seidel
!       iterations used in the computation of G, DFMAX de-
!       notes the maximum relative change in a solution
!       component in the last Gauss-Seidel iteration, and
!       DP is the change in P computed by linear interpola-
!       tion between the current point (P,G) and a previous
!       point.
!
      NN = N
      TOL = GSTOL
!
! Test for errors in input parameters.
!
      IER = -1
      IF (NCC .LT. 0  .OR.  SM .LE. 0.  .OR.  SMTOL .LE. 0.&
         &.OR.  SMTOL .GE. 1.  .OR.  TOL .LE. 0.) RETURN
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
!
! Compute the components of the 3 by 3 system (normal
!   equations) for the weighted least squares linear fit.
!
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
!
! Solve the system for (FX,FY,F0) where (FX,FY) is the
!   gradient (constant) and F0 = F(0,0).
!
      CC22 = C11*C22 - C12*C12
      CC23 = C11*C23 - C12*C13
      CC33 = C11*C33 - C13*C13
      RR2 = C11*R2 - C12*R1
      RR3 = C11*R3 - C13*R1
      DET = CC22*CC33 - CC23*CC23
      IER = -2
      IF (DET .EQ. 0.  .OR.  CC22 .EQ. 0.  .OR.&
     &    C11 .EQ. 0.) RETURN
      F0 = (CC22*RR3 - CC23*RR2)/DET
      FY = (RR2 - CC23*F0)/CC22
      FX = (R1 - C12*FY - C13*F0)/C11
!
! Compute nodal values and gradients, and accumulate Q2 =
!   (Z-F)**T*W*(Z-F).
!
      Q2 = 0.
      DO 3 I = 1,NN
        F(I) = FX*X(I) + FY*Y(I) + F0
        FXFY(1,I) = FX
        FXFY(2,I) = FY
        Q2 = Q2 + W(I)*(Z(I)-F(I))**2
    3   CONTINUE
!
! Compute bounds on Q2 defined by SMTOL, and test for the
!   constraint satisfied by the linear fit.
!
      Q2MIN = SM*(1.-SMTOL)
      Q2MAX = SM*(1.+SMTOL)
      IF (Q2 .LE. Q2MAX) THEN
!
!   The constraint is satisfied by a planar surface.
!
        IER = 1
        IF (LUN .GE. 0) WRITE (LUN,100)
  100   FORMAT (///1X,'SMSURF:  The constraint is not ',&
     &          'active and the surface is linear.'/)
        RETURN
      ENDIF
!
! Compute G0 = G(0) and print a heading.
!
      IER = 0
      S = 1./SQRT(SM)
      G0 = 1./SQRT(Q2) - S
      IF (LUN .GE. 0) WRITE (LUN,110) SM, TOL, NITMAX, G0
  110 FORMAT (///1X,'SMSURF -- SM = ',E10.4,', GSTOL = ',&
     &        E7.1,', NITMAX = ',I2,', G(0) = ',E15.8)
!
! G(P) is strictly increasing and concave, and G(0) < 0.
!   Initialize parameters for the secant method.  The method
!   uses three points:  (P0,G0), (P,G), and (PNEG,GNEG),
!   where P0 and PNEG are defined implicitly by DP = P - P0
!   and DMAX = P - PNEG.
!
      P = 10.*SM
      DP = P
      DMAX = 0.
      ITER = 0
!
! Top of loop -- compute G.
!
    4 NIT = NITMAX
      DFMAX = TOL
      CALL SMSGS (NCC,LCC,NN,X,Y,Z,LIST,LPTR,LEND,IFLGS,&
     &            SIGMA,W,P, NIT,DFMAX,F,FXFY, IERR)
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
  120 FORMAT (/1X,I2,' -- P = ',E15.8,', G = ',E15.8,&
     &        ', NIT = ',I2,', DFMAX = ',E12.6)
!
!   Test for convergence.
!
      IF (G .EQ. G0  .OR.  (Q2MIN .LE. Q2  .AND.&
     &                      Q2 .LE. Q2MAX)) RETURN
      IF (DMAX .NE. 0.  .OR.  G .GT. 0.) GO TO 6
!
!   Increase P until G(P) > 0.
!
      P = 10.*P
      DP = P
      GO TO 4
!
!   A bracketing interval [P0,P] has been found.
!
    6 IF (G0*G .LE. 0.) THEN
!
!   G0*G < 0.  Update (PNEG,GNEG) to (P0,G0) so that G
!     and GNEG always have opposite signs.
!
        DMAX = DP
        GNEG = G0
      ENDIF
!
!   Compute the change in P by linear interpolation between
!     (P0,G0) and (P,G).
!
    7 DP = -G*DP/(G-G0)
      IF (LUN .GE. 0) WRITE (LUN,130) DP
  130 FORMAT (1X,5X,'DP = ',E15.8)
      IF (ABS(DP) .GT. ABS(DMAX)) THEN
!
!   G0*G > 0 and the new estimate would be outside of the
!     bracketing interval of length abs(DMAX).  Reset
!     (P0,G0) to (PNEG,GNEG).
!
        DP = DMAX
        G0 = GNEG
        GO TO 7
      ENDIF
!
!   Bottom of loop -- update P, DMAX, and G0.
!
      P = P + DP
      DMAX = DMAX + DP
      G0 = G
      GO TO 4
      END
      SUBROUTINE SNHCSH (X, SINHM,COSHM,COSHMM)
      REAL X, SINHM, COSHM, COSHMM
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   03/18/90
!
!   This subroutine computes approximations to the modified
! hyperbolic functions defined below with relative error
! bounded by 4.7E-12 for a floating point number system with
! sufficient precision.  For IEEE standard single precision,
! the relative error is less than 1.E-5 for all x.
!
!   Note that the 13-digit constants in the data statements
! below may not be acceptable to all compilers.
!
! On input:
!
!       X = Point at which the functions are to be
!           evaluated.
!
! X is not altered by this routine.
!
! On output:
!
!       SINHM = sinh(X) - X.
!
!       COSHM = cosh(X) - 1.
!
!       COSHMM = cosh(X) - 1 - X*X/2.
!
! Modules required by SNHCSH:  None
!
! Intrinsic functions called by SNHCSH:  ABS, EXP
!
!***********************************************************
!
      REAL AX, C1, C2, C3, C4, EXPX, F, XC, XS, XSD2, XSD4
!
      DATA C1/.1666666666659E0/,&
     &     C2/.8333333431546E-2/,&
     &     C3/.1984107350948E-3/,&
     &     C4/.2768286868175E-5/
      AX = ABS(X)
      XS = AX*AX
      IF (AX .LE. .5) THEN
!
! Approximations for small X:
!
        XC = X*XS
        SINHM = XC*(((C4*XS+C3)*XS+C2)*XS+C1)
        XSD4 = .25*XS
        XSD2 = XSD4 + XSD4
        F = (((C4*XSD4+C3)*XSD4+C2)*XSD4+C1)*XSD4
        COSHMM = XSD2*F*(F+2.)
        COSHM = COSHMM + XSD2
      ELSE
!
! Approximations for large X:
!
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
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88
!
!   This function computes the integral over a triangle of
! the linear (planar) surface which interpolates data
! values at the vertices.
!
! On input:
!
!       X1,X2,X3 = X coordinates of the vertices of the tri-
!                  angle in counterclockwise order.
!
!       Y1,Y2,Y3 = Y coordinates of the vertices of the tri-
!                  angle in one-to-one correspondence with
!                  X1, X2, and X3.
!
!       Z1,Z2,Z3 = Data values at the vertices (X1,Y1),
!                  (X2,Y2), (X3,Y3), respectively.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       TRVOL = Integral over the triangle of the linear
!               interpolant.  Note that TRVOL will have
!               the wrong sign if the vertices are speci-
!               fied in clockwise order.
!
! Modules required by TRVOL:  None
!
!***********************************************************
!
      REAL AREA
!
      AREA = (X2-X1)*(Y3-Y1) - (X3-X1)*(Y2-Y1)
!
! AREA is twice the (signed) area of the triangle.
! TRVOL is the mean of the data values times the area of the
!   triangle.
!
      TRVOL = (Z1 + Z2 + Z3)*AREA/6.
      RETURN
      END
      SUBROUTINE TVAL (X,Y,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,ZX1,&
     &                 ZX2,ZX3,ZY1,ZY2,ZY3,DFLAG, F,FX,FY,&
     &                 IER)
      INTEGER IER
      LOGICAL DFLAG
      REAL    X, Y, X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3,&
     &        ZX1, ZX2, ZX3, ZY1, ZY2, ZY3, F, FX, FY
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   06/12/90
!
!   Given function values and first partial derivatives at
! the vertices of a triangle, along with a point P in the
! triangle, this subroutine computes an interpolated value
! F(P) and, optionally, the first partial derivatives of F
! at P.
!
!   The interpolant F of the vertex values and gradients is
! the Clough-Tocher finite element.  F is cubic in each of
! the three subtriangles of equal area obtained by joining
! the vertices to the barycenter, but has only quadratic
! precision (exact for values and partials from a quadratic
! polynomial).  Along each triangle side, F is the Hermite
! cubic interpolant of the endpoint values and tangential
! gradient components, and the normal gradient component of
! F varies linearly between the interpolated endpoint nor-
! mal components.  Thus, since values and first partials on
! a triangle side depend only on the endpoint data, the
! method results in a C-1 interpolant over a triangulation.
! Second derivatives are discontinuous across subtriangle
! boundaries.
!
!   The computational procedure, due to Charles Lawson, has
! the following operation counts:  62 adds, 54 multiplies,
! 8 divides, and 6 compares for an interpolated value, and
! 170 adds, 142 multiplies, 14 divides, and 6 compares for
! both a value and a pair of first partial derivatives.
!
! On input:
!
!       X,Y = Coordinates of the point P at which F is to
!             be evaluated.
!
!       X1,X2,X3 = X coordinates of the vertices of the tri-
!                  angle in counterclockwise order.
!
!       Y1,Y2,Y3 = Y coordinates of the vertices of the tri-
!                  angle in one-to-one correspondence with
!                  X1, X2, and X3.
!
!       Z1,Z2,Z3 = Data values at the vertices (X1,Y1),
!                  (X2,Y2), (X3,Y3), respectively.
!
!       ZX1,ZX2,ZX3 = X-derivative values at the vertices.
!
!       ZY1,ZY2,ZY3 = Y-derivative values at the vertices.
!
!       DFLAG = Logical flag which specifies whether first
!               partial derivatives at P are to be computed:
!               DFLAG = .TRUE. if and only if partials are
!               to be returned.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       F = Value of the interpolatory function at P if
!           IER = 0, or zero if IER = 1.  Note that, if
!           P is not contained in the triangle, F is an
!           extrapolated value.
!
!       FX,FY = Partial derivatives of F at P if DFLAG =
!               .TRUE. and IER = 0, unaltered otherwise.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if the vertices of the triangle are
!                     collinear.
!
! Modules required by TVAL:  None
!
!***********************************************************
!
      INTEGER I, IP1, IP2, IP3
      REAL    A(3), AREA, AX(3), AY(3), B(3), BX(3), BY(3),&
     &        C(3), CX(3), CY(3), C1, C2, FF(3), G(3),&
     &        GX(3), GY(3), P(3), PHI(3), PHIX(3), PHIY(3),&
     &        PX(3), PY(3), Q(3), QX(3), QY(3), R(3), RMIN,&
     &        RO(3), ROX(3), ROY(3), RX(3), RY(3), SL(3),&
     &        U(3), V(3), XP, YP
!
! Local parameters:
!
! A(K) =            Cardinal function whose coefficient is
!                     Z(K)
! AREA =            Twice the area of the triangle
! AX(K),AY(K) =     X,Y partials of A(K) -- cardinal
!                     functions for FX and FY
! B(K) =            Twice the cardinal function whose
!                     coefficient is ZX(K)
! BX(K),BY(K) =     X,Y partials of B(K)
! C(K) =            Twice the cardinal function whose
!                     coefficient is ZY(K)
! CX(K),CY(K) =     X,Y partials of C(K)
! C1,C2 =           Factors for computing RO
! FF(K) =           Factors for computing G, GX, and GY --
!                     constant
! G(K) =            Factors for computing the cardinal
!                     functions -- cubic
! GX(K),GY(K) =     X,Y partials of G(K)
! I =               DO-loop index
! IP1,IP2,IP3 =     Permuted indexes for computing RO, ROX,
!                     and ROY
! P(K) =            G(K) + PHI(K)
! PHI(K)            R(K-1)*R(K+1) -- quadratic
! PHIX(K),PHIY(K) = X,Y partials of PHI(K)
! PX(K),PY(K) =     X,Y partials of P(K)
! Q(K) =            G(K) - PHI(K)
! QX(K),QY(K) =     X,Y partials of Q(K)
! R(K) =            K-th barycentric coordinate
! RMIN =            Min(R1,R2,R3)
! RO(K) =           Factors for computing G -- cubic
!                     correction terms
! ROX(K),ROY(K) =   X,Y partials of RO(K)
! RX(K),RY(K) =     X,Y partial derivatives of R(K)
! SL(K) =           Square of the length of the side
!                     opposite vertex K
! U(K) =            X-component of the vector representing
!                     the side opposite vertex K
! V(K) =            Y-component of the vector representing
!                     the side opposite vertex K
! XP,YP =           X-X1, Y-Y1
!
      U(1) = X3 - X2
      U(2) = X1 - X3
      U(3) = X2 - X1
!
      V(1) = Y3 - Y2
      V(2) = Y1 - Y3
      V(3) = Y2 - Y1
!
      DO 1 I = 1,3
        SL(I) = U(I)*U(I) + V(I)*V(I)
    1   CONTINUE
!
! AREA = 3->1 X 3->2.
!
      AREA = U(1)*V(2) - U(2)*V(1)
      IF (AREA .EQ. 0.) GO TO 9
      IER = 0
!
! R(1) = (2->3 X 2->P)/AREA, R(2) = (1->P X 1->3)/AREA,
!   R(3) = (1->2 X 1->P)/AREA.
!
      R(1) = (U(1)*(Y-Y2) - V(1)*(X-X2))/AREA
      XP = X - X1
      YP = Y - Y1
      R(2) = (U(2)*YP - V(2)*XP)/AREA
      R(3) = (U(3)*YP - V(3)*XP)/AREA
!
      PHI(1) = R(2)*R(3)
      PHI(2) = R(3)*R(1)
      PHI(3) = R(1)*R(2)
!
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
!
    5 C1 = RMIN*RMIN/2.
      C2 = RMIN/3.
      RO(IP1) = (PHI(IP1) + 5.*C1/3.)*R(IP1) - C1
      RO(IP2) = C1*(R(IP3) - C2)
      RO(IP3) = C1*(R(IP2) - C2)
!
      FF(1) = 3.*(SL(2)-SL(3))/SL(1)
      FF(2) = 3.*(SL(3)-SL(1))/SL(2)
      FF(3) = 3.*(SL(1)-SL(2))/SL(3)
!
      G(1) = (R(2)-R(3))*PHI(1) + FF(1)*RO(1) - RO(2)+RO(3)
      G(2) = (R(3)-R(1))*PHI(2) + FF(2)*RO(2) - RO(3)+RO(1)
      G(3) = (R(1)-R(2))*PHI(3) + FF(3)*RO(3) - RO(1)+RO(2)
!
      DO 6 I = 1,3
        P(I) = G(I) + PHI(I)
        Q(I) = G(I) - PHI(I)
    6   CONTINUE
!
      A(1) = R(1) + G(3) - G(2)
      A(2) = R(2) + G(1) - G(3)
      A(3) = R(3) + G(2) - G(1)
!
      B(1) = U(3)*P(3) + U(2)*Q(2)
      B(2) = U(1)*P(1) + U(3)*Q(3)
      B(3) = U(2)*P(2) + U(1)*Q(1)
!
      C(1) = V(3)*P(3) + V(2)*Q(2)
      C(2) = V(1)*P(1) + V(3)*Q(3)
      C(3) = V(2)*P(2) + V(1)*Q(1)
!
! F is a linear combination of the cardinal functions.
!
      F = A(1)*Z1 + A(2)*Z2 + A(3)*Z3 + (B(1)*ZX1 + B(2)*ZX2&
     &    + B(3)*ZX3 + C(1)*ZY1 + C(2)*ZY2 + C(3)*ZY3)/2.
      IF (.NOT. DFLAG) RETURN
!
! Compute FX and FY.
!
      DO 7 I = 1,3
        RX(I) = -V(I)/AREA
        RY(I) = U(I)/AREA
    7   CONTINUE
!
      PHIX(1) = R(2)*RX(3) + RX(2)*R(3)
      PHIY(1) = R(2)*RY(3) + RY(2)*R(3)
      PHIX(2) = R(3)*RX(1) + RX(3)*R(1)
      PHIY(2) = R(3)*RY(1) + RY(3)*R(1)
      PHIX(3) = R(1)*RX(2) + RX(1)*R(2)
      PHIY(3) = R(1)*RY(2) + RY(1)*R(2)
!
      ROX(IP1) = RX(IP1)*(PHI(IP1) + 5.*C1) +&
     &           R(IP1)*(PHIX(IP1) - RX(IP1))
      ROY(IP1) = RY(IP1)*(PHI(IP1) + 5.*C1) +&
     &           R(IP1)*(PHIY(IP1) - RY(IP1))
      ROX(IP2) = RX(IP1)*(PHI(IP2) - C1) + C1*RX(IP3)
      ROY(IP2) = RY(IP1)*(PHI(IP2) - C1) + C1*RY(IP3)
      ROX(IP3) = RX(IP1)*(PHI(IP3) - C1) + C1*RX(IP2)
      ROY(IP3) = RY(IP1)*(PHI(IP3) - C1) + C1*RY(IP2)
!
      GX(1) = (RX(2) - RX(3))*PHI(1) + (R(2) - R(3))*PHIX(1)&
     &        + FF(1)*ROX(1) - ROX(2) + ROX(3)
      GY(1) = (RY(2) - RY(3))*PHI(1) + (R(2) - R(3))*PHIY(1)&
     &        + FF(1)*ROY(1) - ROY(2) + ROY(3)
      GX(2) = (RX(3) - RX(1))*PHI(2) + (R(3) - R(1))*PHIX(2)&
     &        + FF(2)*ROX(2) - ROX(3) + ROX(1)
      GY(2) = (RY(3) - RY(1))*PHI(2) + (R(3) - R(1))*PHIY(2)&
     &        + FF(2)*ROY(2) - ROY(3) + ROY(1)
      GX(3) = (RX(1) - RX(2))*PHI(3) + (R(1) - R(2))*PHIX(3)&
     &        + FF(3)*ROX(3) - ROX(1) + ROX(2)
      GY(3) = (RY(1) - RY(2))*PHI(3) + (R(1) - R(2))*PHIY(3)&
     &        + FF(3)*ROY(3) - ROY(1) + ROY(2)
!
      DO 8 I = 1,3
        PX(I) = GX(I) + PHIX(I)
        PY(I) = GY(I) + PHIY(I)
        QX(I) = GX(I) - PHIX(I)
        QY(I) = GY(I) - PHIY(I)
    8   CONTINUE
!
      AX(1) = RX(1) + GX(3) - GX(2)
      AY(1) = RY(1) + GY(3) - GY(2)
      AX(2) = RX(2) + GX(1) - GX(3)
      AY(2) = RY(2) + GY(1) - GY(3)
      AX(3) = RX(3) + GX(2) - GX(1)
      AY(3) = RY(3) + GY(2) - GY(1)
!
      BX(1) = U(3)*PX(3) + U(2)*QX(2)
      BY(1) = U(3)*PY(3) + U(2)*QY(2)
      BX(2) = U(1)*PX(1) + U(3)*QX(3)
      BY(2) = U(1)*PY(1) + U(3)*QY(3)
      BX(3) = U(2)*PX(2) + U(1)*QX(1)
      BY(3) = U(2)*PY(2) + U(1)*QY(1)
!
      CX(1) = V(3)*PX(3) + V(2)*QX(2)
      CY(1) = V(3)*PY(3) + V(2)*QY(2)
      CX(2) = V(1)*PX(1) + V(3)*QX(3)
      CY(2) = V(1)*PY(1) + V(3)*QY(3)
      CX(3) = V(2)*PX(2) + V(1)*QX(1)
      CY(3) = V(2)*PY(2) + V(1)*QY(1)
!
! FX and FY are linear combinations of the cardinal
!   functions.
!
      FX = AX(1)*Z1 + AX(2)*Z2 + AX(3)*Z3 + (BX(1)*ZX1 +&
     &     BX(2)*ZX2 + BX(3)*ZX3 + CX(1)*ZY1 + CX(2)*ZY2 +&
     &     CX(3)*ZY3)/2.
      FY = AY(1)*Z1 + AY(2)*Z2 + AY(3)*Z3 + (BY(1)*ZX1 +&
     &     BY(2)*ZX2 + BY(3)*ZX3 + CY(1)*ZY1 + CY(2)*ZY2 +&
     &     CY(3)*ZY3)/2.
      RETURN
!
! The vertices are collinear.
!
    9 IER = 1
      F = 0.
      RETURN
      END
      SUBROUTINE UNIF (NCC,LCC,N,X,Y,Z,GRAD,LIST,LPTR,LEND,&
     &                 IFLGS,SIGMA,NROW,NX,NY,PX,PY,SFLAG,&
     &                 SVAL, ZZ,IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),&
     &        IFLGS, NROW, NX, NY, IER
      LOGICAL SFLAG
      REAL    X(N), Y(N), Z(N), GRAD(2,N), SIGMA(*), PX(NX),&
     &        PY(NY), SVAL, ZZ(NROW,NY)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/22/97
!
!   Given a Delaunay triangulation of a set of points in the
! plane with associated data values and gradients, this sub-
! routine interpolates the data to a set of rectangular grid
! points for such applications as contouring.  Extrapolation
! is performed at grid points exterior to the triangulation,
! and the interpolant is once-continuously differentiable
! over the entire plane.  Refer to Subroutine INTRC1 for
! further details.
!
! On input:
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC .GE. 0.
!
!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index of the
!             first node of constraint I in LCC(I).  For I =
!             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       Z = Array of length N containing data values at the
!           nodes.  Refer to Subroutines ZGRADG and ZGRADL.
!
!       GRAD = 2 by N array whose columns contain estimated
!              gradients at the nodes with X partial deriva-
!              tives in the first row and Y partials in the
!              second.  Refer to Subroutines GRADC, GRADG,
!              GRADL, SMSURF, ZGRADG, and ZGRADL.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
!       IFLGS = Tension factor option:
!               IFLGS .LE. 0 if a single uniform tension
!                            factor is to be used.
!               IFLGS .GE. 1 if variable tension is desired.
!
!       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
!               array containing tension factors associated
!               with arcs in one-to-one correspondence with
!               LIST entries (IFLGS .GE. 1).  Refer to Sub-
!               routines FVAL, GETSIG, SIG0, SIG1, and SIG2.
!
!       NROW = Number of rows in the dimension statement of
!              ZZ.
!
!       NX,NY = Number of rows and columns, respectively, in
!               the rectangular grid.  1 .LE. NX .LE. NROW,
!               and 1 .LE. NY.
!
!       PX,PY = Arrays of length NX and NY, respectively,
!               containing the coordinates of the grid
!               lines.
!
!       SFLAG = Special value flag:
!               SFLAG = .FALSE. if special values are not to
!                               be used (ZZ contains only
!                               interpolated or extrapolated
!                               values.
!               SFLAG = .TRUE. if SVAL is to be stored in ZZ
!                              elements corresponding to
!                              grid points which lie in a
!                              constraint region.
!
!       SVAL = Special value for grid points lying in a con-
!              straint region, or dummy parameter if SFLAG =
!              .FALSE.
!
! The above parameters are not altered by this routine.
!
!       ZZ = NROW by NCOL array for some NCOL .GE. NY.
!
! On output:
!
!       ZZ = Interpolated values at the grid points (or
!            special values) if IER .GE. 0.  ZZ(I,J) =
!            F(PX(I),PY(J)) for I = 1,...,NX and J = 1,...,
!            NY, where F is the interpolatory surface.
!
!       IER = Error indicator:
!             IER .GE. 0 if no errors were encountered.
!                        IER contains the number of grid
!                        points exterior to the triangula-
!                        tion or contained in a constraint
!                        region triangle (extrapolated
!                        values).
!             IER = -1 if NCC, N, NROW, NX, or NY is
!                      outside its valid range on input.
!                      LCC is not tested for validity.
!             IER = -2 if the nodes are collinear or the
!                      triangulation is invalid.
!
! TRIPACK modules required by UNIF:  CRTRI, JRAND, LEFT,
!                                      LSTPTR, TRFIND
!
! SRFPACK modules required by UNIF:  ARCINT, COORDS, FVAL,
!                                      INTRC1, SNHCSH, TVAL
!
!***********************************************************
!
      INTEGER I, IERR, IST, J, NEX, NI, NJ, NST
      LOGICAL DFLAG, SFL
      REAL    DUM
      DATA    DFLAG/.FALSE./,  NST/1/
!
! Local parameters:
!
! DFLAG = Derivative flag for INTRC1
! DUM =   Dummy INTRC1 parameter
! I,J =   DO-loop indexes
! IERR =  Error flag for calls to INTRC1
! IST =   Parameter for INTRC1
! NEX =   Number of grid points exterior to the triangula-
!           tion boundary (number of extrapolated values)
! NI,NJ = Local copies of NX and NY
! NST =   Initial value for IST
! SFL =   Local copy of SFLAG
!
      NI = NX
      NJ = NY
      IF (NCC .LT. 0  .OR.  N .LT. 3  .OR.  NI .LT. 1  .OR.&
     &    NI .GT. NROW  .OR.  NJ .LT. 1) GO TO 3
      SFL = SFLAG
      IST = NST
!
! Compute interpolated values.
!
      NEX = 0
      DO 2 J = 1,NJ
        DO 1 I = 1,NI
          CALL INTRC1 (PX(I),PY(J),NCC,LCC,N,X,Y,Z,LIST,&
     &                 LPTR,LEND,IFLGS,SIGMA,GRAD,&
     &                 DFLAG, IST, ZZ(I,J),DUM,DUM,IERR)
          IF (IERR .LT. 0) GO TO 4
          IF (IERR .GT. 0) NEX = NEX + 1
          IF (SFL  .AND. IERR .EQ. 1) ZZ(I,J) = SVAL
    1     CONTINUE
    2   CONTINUE
      IER = NEX
      RETURN
!
! Invalid input parameter.
!
    3 IER = -1
      RETURN
!
! Triangulation nodes are collinear.
!
    4 IER = -2
      RETURN
      END
      REAL FUNCTION VOLUME (NCC,LCC,N,X,Y,Z,LIST,LPTR,LEND)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N)
      REAL    X(N), Y(N), Z(N)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   08/26/91
!
!   Given a triangulation of a set of N nodes, along with
! data values at the nodes, this function computes the int-
! egral over a region R of the piecewise linear interpolant
! of the data values.  R is the convex hull of the nodes
! with constraint regions excluded.
!
! On input:
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC .GE. 0.
!
!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index of the
!             first node of constraint I in LCC(I).  For I =
!             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       Z = Array of length N containing data values at the
!           nodes.  Refer to Subroutine ZGRADL.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       VOLUME = Sum of the volumes of the linear interpo-
!                lants on the non-constraint triangles, or
!                zero if a parameter is outside its valid
!                range on input.
!
! SRFPACK module required by VOLUME:  TRVOL
!
! Intrinsic function called by VOLUME:  ABS
!
!***********************************************************
!
      REAL    TRVOL
      INTEGER I, ILAST, LCC1, LP2, LP3, LPL, N1, N2, N3,&
     &        NM2, NN
      REAL    SUM, XN1, YN1, ZN1
!
! Test for invalid input parameters.
!
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
!
! Initialize for loop on triangles (N1,N2,N3) such that N2
!   and N3 have larger indexes than N1.  SUM contains the
!   accumulated volume, I is the index of the constraint
!   containing N1 if N1 is a constraint node, and ILAST is
!   the last node of constraint I.
!
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
!
! Top of loop on neighbors of N1.
!
        LPL = LEND(N1)
        LP2 = LPL
    2   LP2 = LPTR(LP2)
          N2 = LIST(LP2)
          LP3 = LPTR(LP2)
          N3 = ABS(LIST(LP3))
          IF (N2 .LT. N1  .OR.  N3 .LT. N1) GO TO 3
!
!   (N1,N2,N3) lies in a constraint region iff the vertices
!     are nodes of the same constraint and N2 < N3.
!
          IF (N1 .LT. LCC1  .OR.  N2 .GT. N3  .OR.&
     &        N3 .GT. ILAST) THEN
            SUM = SUM + TRVOL(XN1,X(N2),X(N3),YN1,Y(N2),&
     &                        Y(N3),ZN1,Z(N2),Z(N3))
          ENDIF
!
!   Bottom of loop on neighbors.
!
    3     IF (LP2 .NE. LPL) GO TO 2
    4   CONTINUE
!
      VOLUME = SUM
      RETURN
!
! Invalid input parameter.
!
    5 VOLUME = 0.
      RETURN
      END
      SUBROUTINE ZGRADG (NCC,LCC,N,X,Y,LIST,LPTR,LEND,&
     &                   IFLGS,SIGMA, NIT,DZMAX,Z,&
     &                   GRAD, IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),&
     &        IFLGS, NIT, IER
      REAL    X(N), Y(N), SIGMA(*), DZMAX, Z(N), GRAD(2,N)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/22/97
!
!   Given a triangulation of N nodes, along with data values
! at non-constraint nodes, this subroutine employs a global
! method to compute estimated gradients at the nodes and
! approximate data values at the constraint nodes.  For
! NCN = N-LCC(1)+1 constraint nodes, the method consists of
! minimizing a quadratic functional Q(U) over vectors U of
! length 2N+NCN containing gradients and values at con-
! straint nodes.  Q is taken to be the sum over the
! triangulation arcs, excluding interior constraint arcs,
! of the linearized curvature (integral of squared second
! derivative) of the Hermite interpolatory tension spline
! defined by the data values and tangential gradient comp-
! onents at the endpoints of the arc.
!
!   This minimization problem corresponds to a symmetric
! positive definite sparse linear system which is solved by
! a block Gauss-Seidel method with N blocks of order 2, or
! order 3 for constraint nodes.
!
!   An alternative method (Subroutine ZGRADL) computes a
! local approximation to the gradient and data value (if not
! specified) at a single node.  The relative speed and
! accuracy of the two methods depends on the distribution
! of the nodes.  Relative accuracy also depends on the data
! values.
!
!   Note that the call to ZGRADG can be followed by a call
! to GRADG in order to compute improved gradient estimates
! with fixed data values.  This is recommended for improved
! efficiency.
!
! On input:
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC > 0.
!
!       LCC = Array of length NCC containing the index of
!             the first node of constraint I in LCC(I).  For
!             I = 1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1, and LCC(1) .GE. 4.
!
!       N = Number of nodes in the triangulation.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
!       IFLGS = Tension factor option:
!               IFLGS .LE. 0 if a single uniform tension
!                            factor is to be used.
!               IFLGS .GE. 1 if variable tension is desired.
!
!       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
!               array containing tension factors associated
!               with arcs in one-to-one correspondence with
!               LIST entries (IFLGS .GE. 1).  Refer to Sub-
!               routines GETSIG, SIG0, SIG1, and SIG2.
!
! The above parameters are not altered by this routine.
!
!       NIT = Maximum number of Gauss-Seidel iterations to
!             be employed.  This maximum will likely be
!             achieved if DZMAX is smaller than the machine
!             precision.  Note that complete convergence is
!             not necessary to achieve maximum accuracy of
!             the interpolant.  NIT > 0.
!
!       DZMAX - Nonnegative convergence criterion.  The
!               method is terminated when the maximum change
!               in a solution Z-component between iterations
!               is at most DZMAX.  The change in a solution
!               component is taken to be the magnitude of
!               the difference relative to 1 plus the magni-
!               tude of the previous value.
!
!       Z = Array of length N containing data values in the
!           first LCC(1)-1 locations and initial solution
!           estimates in the remaining locations.  Zeros are
!           sufficient, but Subroutine ZINIT may be called
!           to provide better initial estimates.
!
!       GRAD = 2 by N array whose columns contain initial
!              estimates of the gradients with X partial
!              derivatives in the first row, Y partials in
!              the second.  Zeros are sufficient.
!
! On output:
!
!       NIT = Number of Gauss-Seidel iterations employed.
!
!       DZMAX = Maximum relative change in a solution Z-
!               component at the last iteration.
!
!       Z = Array updated with approximate data values in
!           the last NCN = N-LCC(1)+1 locations if IER .GE.
!           0.  Z is not altered if IER = -1.
!
!       GRAD = Estimated gradients at the nodes if IER .GE.
!              0.  GRAD is not altered if IER = -1.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and the
!                     convergence criterion was achieved.
!             IER = 1 if no errors were encountered but con-
!                     vergence was not achieved within NIT
!                     iterations.
!             IER = -1 if NCC, an LCC entry, N, NIT, or
!                      DZMAX is outside its valid range
!                      on input.
!             IER = -2 if all nodes are collinear or the
!                      triangulation data structure is in-
!                      valid.
!             IER = -3 if duplicate nodes were encountered.
!
! SRFPACK modules required by ZGRADG:  GRCOEF, SNHCSH
!
! Intrinsic functions called by ZGRADG:  ABS, MAX, SQRT
!
!***********************************************************
!
      INTEGER I, IFL, IFRST, ILAST, ITER, J, JN, K, KBAK,&
     &        KFOR, LCC1, LP, LPF, LPJ, LPL, LPN, MAXIT, NB,&
     &        NN
      REAL    A11, A12, A13, A22, A23, A33, AREAJ, AREAN,&
     &        AREAP, D, DCUB, DF, DSQ, DX, DXS, DY, DYS, DZ,&
     &        DZJ, DZK, DZMX, DZX, DZY, R1, R2, R3, SDF,&
     &        SIG, T, TOL, W, XK, YK, ZK, ZXK, ZYK
!
      NN = N
      IFL = IFLGS
      MAXIT = NIT
      TOL = DZMAX
!
! Test for errors in input parameters.
!
      IF (NCC .LE. 0  .OR.  MAXIT .LT. 1  .OR.  TOL .LT. 0.)&
     &  GO TO 9
      LCC1 = NN+1
      DO 1 I = NCC,1,-1
        IF (LCC1-LCC(I) .LT. 3) GO TO 9
        LCC1 = LCC(I)
    1   CONTINUE
      IF (LCC1 .LT. 4) GO TO 9
!
! Initialize iteration count and SIG (overwritten if
!   IFLGS > 0).
!
      ITER = 0
      SIG = SIGMA(1)
!
! Top of iteration loop:  If K is a constraint node, I
!   indexes the constraint containing node K, IFRST and
!   ILAST are the first and last nodes of constraint I, and
!   (KBAK,K,KFOR) is a subsequence of constraint I.
!
    2 IF (ITER .EQ. MAXIT) GO TO 8
      DZMX = 0.
      I = 0
      ILAST = LCC1-1
      KBAK = 0
      KFOR = 0
!
! Loop on nodes.
!
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
!
! Initialize components of the 2 by 2 (or 3 by 3) block --
!   symmetric matrix in A and residual in R.  The unknowns
!   are ordered (DZX,DZY,DZ).
!
        A11 = 0.
        A12 = 0.
        A13 = 0.
        A22 = 0.
        A23 = 0.
        A33 = 0.
        R1 = 0.
        R2 = 0.
        R3 = 0.
!
! Loop on neighbors J of node K.  The equation associated
!   with K->J (and hence its contribution to the functional)
!   is weighted by AREAJ/D, where AREAJ is twice the sum of
!   the areas of the triangles containing K-J (excluding
!   those which lie in a constraint region) and D is the arc
!   length.  JN is the neighbor of K following J.  AREAP is
!   to the right of K->J and AREAN is to the left.
!
        LPL = LEND(K)
        J = LIST(LPL)
        LPF = LPTR(LPL)
        JN = LIST(LPF)
        AREAN = 0.
        IF (J .GT. 0) AREAN = (X(J)-XK)*(Y(JN)-YK) -&
     &                        (Y(J)-YK)*(X(JN)-XK)
        LPN = LPF
!
! Top of loop:  LPF and LPL point to the first and last
!   neighbors of K, and LPN points to JN.
!
    3   LPJ = LPN
          LPN = LPTR(LPN)
          J = JN
          AREAP = AREAN
          JN = ABS(LIST(LPN))
!
! Arc K-J lies in a constraint region and is bypassed iff K
!   and J are nodes in the same constraint and J follows
!   KFOR and precedes KBAK as a neighbor of K.
!
          IF (K .LT. LCC1  .OR.  J .LT. IFRST  .OR.&
     &        J .GT. ILAST) GO TO 5
          IF (J .EQ. KBAK) AREAP = 0.
          IF (J .EQ. KBAK  .OR.  J .EQ. KFOR) GO TO 5
!
          LP = LPN
    4     NB = ABS(LIST(LP))
            IF (NB .EQ. KFOR) GO TO 5
            IF (NB .EQ. KBAK) GO TO 6
            LP = LPTR(LP)
            GO TO 4
!
!   Compute parameters associated with the edge K->J, and
!     test for duplicate nodes.  Note that AREAJ = 0 and
!     K->J is bypassed if K-J is both a constraint arc and
!     a boundary arc of the triangulation.
!
    5     DX = X(J) - XK
          DY = Y(J) - YK
          AREAN = 0.
          IF (LIST(LPL) .NE. -J  .AND.  J .NE. KFOR) AREAN =&
     &      DX*(Y(JN)-YK) - DY*(X(JN)-XK)
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
!
!   Update the 2 by 2 system components for node J.
!
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
!
!   K is a constraint node.  Update the remaining components.
!
            W = (DF+SDF)*W
            A13 = A13 + DX*W
            A23 = A23 + DY*W
            A33 = A33 + 2.0*W
            R3 = R3 + (2.0*DZ - DZJ - DZK)*W
          ENDIF
!
!   Bottom of loop on J.
!
    6     IF (LPN .NE. LPF) GO TO 3
!
! Solve the linear system associated with the K-th block.
!
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
!
! Update the solution components for node K and the maxi-
!   mum relative change DZMX.
!
        GRAD(1,K) = ZXK + DZX
        GRAD(2,K) = ZYK + DZY
        IF (K .GE. LCC1) THEN
          Z(K) = ZK + DZ
          DZMX = MAX(DZMX,ABS(DZ)/(1.+ABS(ZK)))
        ENDIF
    7   CONTINUE
!
! Increment ITER and test for convergence.
!
      ITER = ITER + 1
      IF (DZMX .GT. TOL) GO TO 2
!
! Method converged.
!
      NIT = ITER
      DZMAX = DZMX
      IER = 0
      RETURN
!
! Method failed to converge within NIT iterations.
!
    8 DZMAX = DZMX
      IER = 1
      RETURN
!
! Invalid input parameter.
!
    9 NIT = 0
      DZMAX = 0.
      IER = -1
      RETURN
!
! Node K and its neighbors are collinear, resulting in a
!   singular system.
!
   10 NIT = 0
      DZMAX = DZMX
      IER = -2
      RETURN
!
! Nodes J and K coincide.
!
   11 NIT = 0
      DZMAX = DZMX
      IER = -3
      RETURN
      END
      SUBROUTINE ZGRADL (K,NCC,LCC,N,X,Y,LIST,LPTR,&
     &                   LEND, NDV,Z,NPTS,DS, DX,DY,IER)
      INTEGER K, NCC, LCC(*), N, LIST(*), LPTR(*),&
     &        LEND(N), NDV, NPTS(*), IER
      REAL    X(N), Y(N), Z(N), DS(*), DX, DY
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/22/97
!
!   Given a Delaunay triangulation of N nodes, along with
! data values at non-constraint nodes, this subroutine com-
! putes an estimated gradient and, if K is a constraint
! node, an approximate data value, at node K.  The values
! are taken from a quadratic function which fits the data
! values at a set of non-constraint nodes close to K in a
! weighted least squares sense.  If K is not a constraint
! node, the fitting function interpolates the data value at
! node K.  If there are fewer than six data values (non-
! constraint nodes), a linear fitting function is used.
! Also, a Marquardt stabilization factor is used if neces-
! sary to ensure a well-conditioned system.  Thus, a unique
! solution exists unless the non-constraint nodes are col-
! linear.
!
!   An alternative routine, ZGRADG, employs a global method
! to compute values at constraint nodes and gradients at all
! of the nodes at once.  The relative speed and accuracy of
! the two methods depends on whether or not constraints are
! present and on the distribution of the nodes.  Relative
! accuracy also depends on the data values.
!
!   This subroutine may be used for the following purposes:
!
! 1)  to compute gradient estimates and constraint node
!       values for INTRC1 or UNIF,
! 2)  to fill in missing constraint node values for linear
!       interpolation (INTRC0 and VOLUME) or smoothing
!       (SMSURF), and
! 3)  to provide initial estimates for GRADG or ZGRADG
!       (probably a waste of computing time).
!
! If data values at the constraint nodes are known, Subrou-
! tine GRADC or GRADL should be used in place of this
! routine.  If there are no constraints this routine differs
! from GRADL only in providing more flexibility (refer to
! NDV below).
!
! On input:
!
!       K = Index of the node at which a gradient is to be
!           computed.  1 .LE. K .LE. N.
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC .GE. 0.
!
!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index of the
!             first node of constraint I in LCC(I).  For I =
!             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1, and LCC(1) .GE. 4.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
! The above parameters are not altered by this routine.
!
!       NDV = Number of data values (including Z(K) if K is
!             not a constraint node) to be used in the least
!             squares fit, unless more equations are re-
!             quired for stability.  3 .LE. NDV .LT. LCC(1).
!             Note that a linear fitting function will be
!             used if NDV < 6 on input.  A reasonable value
!             is NDV = 9.
!
!       Z = Array of length N containing data values in the
!           first LCC(1)-1 locations.
!
!       NPTS,DS = Arrays of length at least Min(L+1,N) where
!                 L is defined below.  (Length N is suffi-
!                 cient.)
!
! On output:
!
!       NDV = Number of data values (non-constraint nodes)
!             used in the least squares fit.
!
!       Z = Array updated with an approximate data value at
!           node K in Z(K) if K .GE. LCC(1) and IER = 0.
!
!       NPTS = Array containing the indexes of the ordered
!              sequence of L closest nodes to node K (with K
!              in the first position) unless IER .NE. 0.  L
!              is the smallest integer such that the se-
!              quence contains NDV (output value) non-
!              constraint nodes.  NPTS(L+1) = 0 if L < N.
!
!       DS = Array containing the distance between node K
!            and NPTS(I) in DS(I) for I = 1,...,L unless
!            IER .NE. 0.  Distance is measured within the
!            non-constraint region (refer to Subroutine
!            GETNP).
!
!       DX,DY = Estimated X and Y partial derivatives at
!               node K unless IER .NE. 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if K, NCC, an LCC entry, N, or NDV is
!                     outside its valid range on input.
!             IER = 2 if all non-constraint nodes are col-
!                     linear.
!
! TRIPACK modules required by ZGRADL:  GETNP, INTSEC
!
! SRFPACK modules required by ZGRADL:  GIVENS, ROTATE,
!                                        SETRO2
!
! Intrinsic functions called by ZGRADL:  ABS, MIN
!
!***********************************************************
!
      INTEGER I, IERR, IR, IROW1, J, JP1, JR, KK, L, LCC1,&
     &        LNP, LR, ND, NDMIN, NP, NPAR, NPM1, NPP1
      REAL    A(7,7), C, DMIN, DTOL, RFAC, RIN, S, SF, SFS,&
     &        STF, W, XK, YK, ZK
      LOGICAL INIT, STABL
      DATA    RFAC/1.05/,  DTOL/.01/
!
! Store parameters in local variables, test for errors, and
!   initialize switches.
!
      KK = K
      IF (NCC .GT. 0) THEN
        LCC1 = LCC(1)
      ELSE
        LCC1 = N+1
      ENDIF
      NDMIN = NDV
      IF (KK .LT. 1  .OR.  KK .GT. N  .OR.  NCC .LT. 0&
         &.OR.  LCC1 .LT. 4  .OR.  NDMIN .LT. 3  .OR.&
     &    NDMIN .GE. LCC1) GO TO 13
      XK = X(KK)
      YK = Y(KK)
      ZK = 0.
      IF (KK .LT. LCC1) ZK = Z(KK)
      INIT = .FALSE.
      STABL = .FALSE.
!
! Set NPTS to the closest LNP nodes to K, where LNP is the
!   smallest integer such that NPTS contains NDMIN non-
!   constraint nodes.  ND is the number of non-constraint
!   nodes currently in NPTS.
!
      LNP = 1
      NPTS(1) = KK
      DS(1) = 0.
      ND = 0
      IF (KK .LT. LCC1) ND = 1
!
!   Get a new non-constraint node.
!
    1 LNP = LNP + 1
      CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,LNP, NPTS,&
     &            DS, IERR)
      IF (IERR .NE. 0) GO TO 13
      IF (NPTS(LNP) .GE. LCC1) GO TO 1
      ND = ND + 1
      IF (ND .LT. NDMIN) GO TO 1
!
! Compute an inverse radius of influence to be used in the
!   weights, and test the initialization switch -- INIT =
!   .TRUE. iff A has been initialized with the first NPAR
!   equations.
!
      RIN = 1./(RFAC*DS(LNP))
      IF (INIT) GO TO 5
!
! A Q-R decomposition is used to solve the least squares
!   system.  For a quadratic fit there are NPAR = 5 or
!   NPAR = 6 parameters, depending on whether or not K is
!   a constraint node.  (At least 6 data values are needed
!   in either case.)  The transpose of the augmented regres-
!   sion matrix is stored in A with columns (rows of A) de-
!   fined as follows -- 1-3 are the quadratic terms, 4 and 5
!   are the linear terms with coefficients DX and DY, column
!   6 is the constant term with coefficient Z(K) (extraneous
!   if NPAR = 5), and the last column is the right hand
!   side.  In the case of a linear fit, the first 3 columns
!   are ignored and the first 3 rows are omitted.  The lin-
!   ear terms are scaled by SF = 1/DMAX, where DMAX is the
!   maximum distance between K and a non-constraint node in
!   NPTS, and the quadratic terms are scaled by SF**2.
!
      SF = 1./DS(LNP)
      SFS = SF*SF
      IROW1 = 1
      IF (ND .LT. 6) IROW1 = 4
      NPAR = 5
      IF (KK .GE. LCC1) NPAR = 6
      NPM1 = NPAR - 1
      NPP1 = NPAR + 1
!
! Set up the first NPAR equations and zero out the lower
!   triangle (upper triangle of A) with Givens rotations --
!
      L = 1
      DO 4 IR = IROW1,NPAR
    2   L = L + 1
          NP = NPTS(L)
          IF (NP .GE. LCC1) GO TO 2
        W = 1./DS(L) - RIN
        CALL SETRO2 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,&
     &               W, A(1,IR))
        IF (IR .EQ. IROW1) GO TO 4
        DO 3 JR = IROW1,IR-1
          JP1 = JR + 1
          LR = 7 - JR
          CALL GIVENS (A(JR,JR),A(JR,IR),C,S)
          CALL ROTATE (LR,C,S,A(JP1,JR),A(JP1,IR))
    3     CONTINUE
    4   CONTINUE
      INIT = .TRUE.
!
! Incorporate additional equations into the system using the
!   last column of A (or next to last if NPAR = 5).
!
    5 IF (L .EQ. LNP) GO TO 7
        L = L + 1
        NP = NPTS(L)
        IF (NP .GE. LCC1) GO TO 5
        W = 1./DS(L) - RIN
        CALL SETRO2 (XK,YK,ZK,X(NP),Y(NP),Z(NP),SF,SFS,&
     &               W, A(1,NPP1))
        DO 6 JR = IROW1,NPAR
          JP1 = JR + 1
          LR = 7 - JR
          CALL GIVENS (A(JR,JR),A(JR,NPP1),C,S)
          CALL ROTATE (LR,C,S,A(JP1,JR),A(JP1,NPP1))
    6     CONTINUE
        GO TO 5
!
! Test the system for ill-conditioning.
!
    7 DMIN = ABS(A(NPAR,NPAR))
      DO 8 I = IROW1,NPM1
        DMIN = MIN(DMIN,ABS(A(I,I)))
    8   CONTINUE
      IF (DMIN/W .GE. DTOL) GO TO 12
      IF (ND .LT. LCC1) THEN
!
!   Add another equation to the system and increase the
!     radius R.
!
        NDMIN = NDMIN + 1
        GO TO 1
      ENDIF
!
! The system is ill-conditioned and all non-constraint nodes
!   have been used.  Stabilize the system by damping out the
!   second partials (coefficients of the quadratic terms)
!   unless the system has already been stabilized or a
!   linear fitting function is being used.  Add multiples
!   of the first 3 unit vectors to the first 3 equations.
!
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
!
! Solve the 2 by 2 (or 3 by 3 if K is a constraint node)
!   lower triangular system.
!
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
!
! Invalid input parameter.
!
   13 NDV = 0
      IER = 1
      RETURN
!
! No unique solution due to collinear non-constraint nodes.
!
   14 NDV = 0
      IER = 2
      RETURN
      END
      SUBROUTINE ZINIT (NCC,LCC,N,X,Y,LIST,LPTR,&
     &                  LEND, Z, IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),&
     &        IER
      REAL    X(N), Y(N), Z(N)
!
!***********************************************************
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   08/27/91
!
!   Given a triangulation of N nodes, along with data values
! at non-constraint nodes, this subroutine computes approxi-
! mate data values at the constraint nodes.  The approximate
! values are intended only to serve as initial estimates for
! Subroutine ZGRADG which computes refined estimates.
!
!   For each subsequence (KM2,KM1,K) of a constraint, the
! approximate value at node KM1 is taken to be the closest-
! point value (data value at the closest non-constraint
! node) at KM1 averaged with the value at KM1 of the linear
! interpolant (along the constraint boundary) of the approx-
! imate value at KM2 and the closest-point value at K.
!
! On input:
!
!       NCC = Number of constraint curves (refer to TRIPACK
!             Subroutine ADDCST).  NCC .GE. 0.
!
!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index of the
!             first node of constraint I in LCC(I).  For I =
!             1 to NCC, LCC(I+1)-LCC(I) .GE. 3, where
!             LCC(NCC+1) = N+1, and LCC(1) .GE. 4.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations, followed by NCC se-
!             quences of constraint nodes.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRIPACK
!                        Subroutine TRMESH.
!
! The above parameters are not altered by this routine.
!
!       Z = Array of length N containing data values in the
!           first LCC(1)-1 locations.
!
! On output:
!
!       Z = Array updated with approximate data values in
!           the last N-LCC(1)+1 locations if IER = 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if NCC, N, or an LCC entry is outside
!                     its valid range on input.
!
! TRIPACK modules required by ZINIT:  GETNP, INTSEC
!
! Intrinsic functions called by ZINIT:  ABS, SQRT
!
!***********************************************************
!
      INTEGER   LMAX
      PARAMETER (LMAX=12)
      INTEGER   I, IERR, IFRST, ILAST, ILSTM1, K, KM1, KM2,&
     &          KN, LCC1, LNP, LP, LPL, NPTS(LMAX)
      REAL      D, DMIN, DS(LMAX), H1, H2, XK, YK, ZN
!
! Test for errors in input parameters.  (LCC is tested by
!   Subroutine GETNP.)
!
      IER = 1
      IF (NCC .GT. 0) THEN
        LCC1 = LCC(1)
      ELSE
        LCC1 = N+1
      ENDIF
      IF (NCC .LT. 0  .OR.  LCC1 .LT. 4) RETURN
!
! Outer loop on constraint I with first and last nodes IFRST
!   and ILAST.
!
      DO 6 I = 1,NCC
        IFRST = LCC(I)
        IF (I .LT. NCC) THEN
          ILAST = LCC(I+1) - 1
        ELSE
          ILAST = N
        ENDIF
!
! Initialize Z(ILAST) with the data value at the closest
!   non-constraint node to ILAST.  Unless the LMAX closest
!   nodes to ILAST (including ILAST) are all constraint
!   nodes, NPTS is set to the closest LNP nodes (with
!   distance measured in the non-constraint region), where
!   LNP is the smallest integer such that NPTS contains a
!   non-constraint node.  The value at LCC(1)-1 is used if
!   LMAX is too small.
!
        LNP = 1
        NPTS(1) = ILAST
        DS(1) = 0.
    1   LNP = LNP + 1
          CALL GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,&
     &                LNP, NPTS,DS, IERR)
          IF (IERR .NE. 0) RETURN
          KN = NPTS(LNP)
          IF (KN .GE. LCC1  .AND.  LNP .LT. LMAX) GO TO 1
        IF (KN .GE. LCC1) KN = LCC1-1
        Z(ILAST) = Z(KN)
!
! Loop on constraint nodes K.  LPL points to the last
!   neighbor of K.  At each step, Z(K) is set to the
!   closest-point value at K, and Z(KM1) is set to the
!   (final) approximate data value at KM1 (except when
!   K = IFRST).
!
        KM1 = ILAST
        ILSTM1 = ILAST - 1
        DO 5 K = IFRST,ILSTM1
          XK = X(K)
          YK = Y(K)
          LPL = LEND(K)
!
!   Set LP to point to KM1 as a neighbor of K.
!
          LP = LPL
    2     LP = LPTR(LP)
            IF (ABS(LIST(LP)) .NE. KM1) GO TO 2
!
!   Initialize for loop on non-constraint node neighbors of
!     K.  If K has no such neighbors, the closest non-
!     constraint node to K is (implicitly) taken to be the
!     closest non-constraint node to KM1.
!
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
!
!   ZN is the closest-point value at K.  Set H2 to the arc
!     length of KM1-K, and compute Z(KM1) if K > IFRST.
!     (H1 is the arc length of KM2-KM1).
!
    4     H2 = SQRT( (XK-X(KM1))**2 + (YK-Y(KM1))**2 )
          IF (K .NE. IFRST) Z(KM1) = .5*( Z(KM1) +&
     &      (H1*ZN+H2*Z(KM2))/(H1+H2) )
          Z(K) = ZN
!
!   Bottom of loop on K.
!
          H1 = H2
          KM2 = KM1
          KM1 = K
    5     CONTINUE
!
! For K = ILAST, the closest-point value has already been
!   computed.
!
        H2 = SQRT ( (X(ILAST)-X(ILSTM1))**2 +&
     &              (Y(ILAST)-Y(ILSTM1))**2 )
        Z(ILSTM1) = .5*( Z(ILSTM1) + (H1*Z(ILAST)+H2*Z(KM2))&
     &                              /(H1+H2) )
!
! Compute the final value at ILAST.
!
        H1 = H2
        H2 = SQRT ( (X(IFRST)-X(ILAST))**2 +&
     &              (Y(IFRST)-Y(ILAST))**2 )
        Z(ILAST) = .5*(Z(ILAST) + (H1*Z(IFRST)+H2*Z(ILSTM1))&
     &                           /(H1+H2))
    6   CONTINUE
!
! No errors encountered.
!
      IER = 0
      RETURN
      END
