      include 'stripack.f'
      SUBROUTINE APLYR (X,Y,Z,CX,SX,CY,SY, XP,YP,ZP)
      REAL X, Y, Z, CX, SX, CY, SY, XP, YP, ZP
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This subroutine applies the rotation R defined by Sub-
C routine CONSTR to the unit vector (X Y Z)**T, i,e. (X,Y,Z)
C is rotated to (XP,YP,ZP).  If (XP,YP,ZP) lies in the
C southern hemisphere (ZP < 0), (XP,YP) are set to the
C coordinates of the nearest point of the equator, ZP re-
C maining unchanged.
C
C On input:
C
C       X,Y,Z = Coordinates of a point on the unit sphere.
C
C       CX,SX,CY,SY = Elements of the rotation defined by
C                     Subroutine CONSTR.
C
C Input parameters are not altered except as noted below.
C
C On output:
C
C       XP,YP,ZP = Coordinates of the rotated point on the
C                  sphere unless ZP < 0, in which case
C                  (XP,YP,0) is the closest point of the
C                  equator to the rotated point.  Storage
C                  for XP, YP, and ZP may coincide with
C                  storage for X, Y, and Z, respectively,
C                  if the latter need not be saved.
C
C Modules required by APLYR:  None
C
C Intrinsic function called by APLYR:  SQRT
C
C***********************************************************
C
      REAL T
C
C Local parameter:
C
C T = Temporary variable
C
      T = SX*Y + CX*Z
      YP = CX*Y - SX*Z
      ZP = SY*X + CY*T
      XP = CY*X - SY*T
      IF (ZP .GE. 0.) RETURN
C
C Move (XP,YP,ZP) to the equator.
C
      T = SQRT(XP*XP + YP*YP)
      IF (T .EQ. 0.) GO TO 1
      XP = XP/T
      YP = YP/T
      RETURN
C
C Move the south pole to an arbitrary point of the equator.
C
    1 XP = 1.
      YP = 0.
      RETURN
      END
      SUBROUTINE APLYRT (G1P,G2P,CX,SX,CY,SY, G)
      REAL G1P, G2P, CX, SX, CY, SY, G(3)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This subroutine applies the inverse (transpose) of the
C rotation defined by Subroutine CONSTR to the vector
C (G1P G2P 0)**T, i.e., the gradient (G1P,G2P,0) in the rot-
C ated coordinate system is mapped to (G1,G2,G3) in the
C original coordinate system.
C
C On input:
C
C       G1P,G2P = X and Y components, respectively, of the
C                 gradient in the rotated coordinate system.
C
C       CX,SX,CY,SY = Elements of the rotation R constructed
C                     by Subroutine CONSTR.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       G = X, Y, and Z components (in that order) of the
C           inverse rotation applied to (G1P,G2P,0) --
C           gradient in the original coordinate system.
C
C Modules required by APLYRT:  None
C
C***********************************************************
C
      REAL T
C
C Local parameters:
C
C T = Temporary variable
C
      T = SY*G1P
      G(1) = CY*G1P
      G(2) = CX*G2P - SX*T
      G(3) = -SX*G2P - CX*T
      RETURN
      END
      SUBROUTINE ARCINT (P,P1,P2,F1,F2,G1,G2,SIGMA, F,G,GN)
      REAL    P(3), P1(3), P2(3), F1, F2, G1(3), G2(3),
     .        SIGMA, F, G(3), GN
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   Given 3 points P, P1, and P2 lying on a common geodesic
C of the unit sphere with P between P1 and P2, along with
C data values and gradients at P1 and P2, this subroutine
C computes an interpolated value F and a gradient vector G
C AT P.  F and the tangential component of G are taken to be
C the value and derivative (with respect to arc-length) of
C a Hermite interpolatory tension spline defined by the end-
C point values and tangential gradient components.  The nor-
C mal component of G is obtained by linear interpolation of
C the normal components of the gradients at P1 and P2.
C
C On input:
C
C       P = Cartesian coordinates of a point lying on the
C           arc defined by P1 and P2.  P(1)**2 + P(2)**2 +
C           P(3)**2 = 1.
C
C       P1,P2 = Coordinates of distinct points on the unit
C               sphere defining an arc with length less than
C               180 degrees.
C
C       F1,F2 = Data values associated with P1 and P2,
C               respectively.
C
C       G1,G2 = Gradient vectors associated with P1 and P2.
C               G1 and G2 are orthogonal to P1 and P2,
C               respectively.
C
C       SIGMA = Tension factor associated with P1-P2.
C
C The above parameters are not altered by this routine.
C
C       G = Array of length 3.
C
C On output:
C
C       F = Interpolated value at P.
C
C       G = Interpolated gradient at P.
C
C       GN = Normal component of G with the direction
C            P1 X P2 taken to be positive.  The extrapola-
C            tion procedure requires this component.
C
C   For each vector V, V(1), V(2), and V(3) contain X, Y,
C and Z components, respectively.
C
C SSRFPACK modules required by ARCINT:  ARCLEN, SNHCSH
C
C Intrinsic functions called by ARCINT:  ABS, EXP, SQRT
C
C***********************************************************
C
      REAL    ARCLEN
      INTEGER I, LUN
      REAL    A, AL, B1, B2, CM, CMM, CM2, DUMMY, D1, D2, E,
     .        EMS, E1, E2, GT, S, SB1, SB2, SIG, SINH,
     .        SINH2, SM, SM2, TAU1, TAU2, TM, TM1, TM2, TP1,
     .        TP2, TS, UN(3), UNORM
      DATA    LUN/6/
C
C Local parameters:
C
C A =         Angle in radians (arc-length) between P1 and
C               P2
C AL =        Arc-length between P1 and P
C B1,B2 =     Local coordinates of P with respect to P1-P2
C CM,CMM =    Coshm(SIG) and Coshmm(SIG) -- refer to SNHCSH
C CM2 =       Coshm(SB2)
C DUMMY =     Dummy parameter for SNHCSH
C D1,D2 =     Scaled second differences
C E =         CM**2 - SM*Sinh = SIG*SM - 2*CMM (scaled by
C               2*EMS if SIG > .5)
C EMS =       Exp(-SIG)
C E1,E2 =     Exp(-SB1), Exp(-SB2)
C GT =        Tangential component of G -- component in the
C               direction UN X P
C I =         DO-loop index
C LUN =       Logical unit for error messages
C S =         Slope:  (F2-F1)/A
C SB1,SB2 =   SIG*B1, SIG*B2
C SIG =       Abs(SIGMA)
C SINH =      Sinh(SIGMA)
C SINH2 =     Sinh(SB2)
C SM,SM2 =    Sinhm(SIG), Sinhm(SB2)
C TAU1,TAU2 = Tangential derivatives (components of G1,G2)
C               at P1 and P2
C TM =        1-EMS
C TM1,TM2 =   1-E1, 1-E2
C TP1,TP2 =   1+E1, 1+E2
C TS =        TM**2
C UN =        Unit normal to the plane of P, P1, and P2
C UNORM =     Euclidean norm of P1 X P2 -- used to normalize
C               UN
C
C
C Compute unit normal UN.
C
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.) GO TO 2
C
C Normalize UN.
C
      DO 1 I = 1,3
        UN(I) = UN(I)/UNORM
    1   CONTINUE
C
C Compute tangential derivatives at the endpoints:
C   TAU1 = (G1,UN X P1) = (G1,P2)/UNORM and
C   TAU2 = (G2,UN X P2) = -(G2,P1)/UNORM.
C
      TAU1 = (G1(1)*P2(1) + G1(2)*P2(2) + G1(3)*P2(3))/UNORM
      TAU2 =-(G2(1)*P1(1) + G2(2)*P1(2) + G2(3)*P1(3))/UNORM
C
C Compute arc-lengths A, AL.
C
      A = ARCLEN(P1,P2)
      IF (A .EQ. 0.) GO TO 2
      AL = ARCLEN(P1,P)
C
C Compute local coordinates, slope, and second differences.
C
      B2 = AL/A
      B1 = 1. - B2
      S = (F2-F1)/A
      D1 = S - TAU1
      D2 = TAU2 - S
C
C Test the range of SIGMA.
C
      SIG = ABS(SIGMA)
      IF (SIG .LT. 1.E-9) THEN
C
C Hermite cubic interpolation.
C
        F = F1 + AL*(TAU1 + B2*(D1 + B1*(D1 - D2)))
        GT = TAU1 + B2*(D1 + D2 + 3.*B1*(D1 - D2))
      ELSEIF (SIG .LE. .5) THEN
C
C 0 < SIG .LE. .5.  Use approximations designed to avoid
C   cancellation error in the hyperbolic functions.
C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        SINH = SM + SIG
        SINH2 = SM2 + SB2
        E = SIG*SM - CMM - CMM
        F = F1 + AL*TAU1 + A*((CM*SM2-SM*CM2)*(D1+D2) + SIG*
     .                        (CM*CM2-SINH*SM2)*D1)/(SIG*E)
        GT = TAU1 + ((CM*CM2-SM*SINH2)*(D1+D2) + SIG*
     .               (CM*SINH2-SINH*CM2)*D1)/E
      ELSE
C
C SIG > .5.  Use negative exponentials in order to avoid
C   overflow.  Note that EMS = EXP(-SIG).
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        E1 = EXP(-SB1)
        E2 = EXP(-SB2)
        EMS = E1*E2
        TM = 1. - EMS
        TS = TM*TM
        TM1 = 1. - E1
        TM2 = 1. - E2
        E = TM*(SIG*(1.+EMS) - TM - TM)
        F = F1 + AL*S + A*(TM*TM1*TM2*(D1+D2) + SIG*
     .                     ((E2*TM1*TM1-B1*TS)*D1 +
     .                      (E1*TM2*TM2-B2*TS)*D2))/(SIG*E)
        TP1 = 1. + E1
        TP2 = 1. + E2
        GT = S + (TM1*(TM*TP2-SIG*E2*TP1)*D1 -
     .            TM2*(TM*TP1-SIG*E1*TP2)*D2)/E
      ENDIF
C
C Compute GN.
C
      GN = B1*(UN(1)*G1(1) + UN(2)*G1(2) + UN(3)*G1(3)) +
     .     B2*(UN(1)*G2(1) + UN(2)*G2(2) + UN(3)*G2(3))
C
C Compute G = GT*(UN X P) + GN*UN.
C
      G(1) = GT*(UN(2)*P(3) - UN(3)*P(2)) + GN*UN(1)
      G(2) = GT*(UN(3)*P(1) - UN(1)*P(3)) + GN*UN(2)
      G(3) = GT*(UN(1)*P(2) - UN(2)*P(1)) + GN*UN(3)
      RETURN
C
C P1 X P2 = 0.  Print an error message and terminate
C   processing.
C
    2 WRITE (LUN,100) (P1(I),I=1,3), (P2(I),I=1,3)
  100 FORMAT ('1','ERROR IN ARCINT -- P1 = ',2(F9.6,',  '),
     .        F9.6/1X,19X,'P2 = ',2(F9.6,',  '),F9.6)
      STOP
      END
      REAL FUNCTION ARCLEN (P,Q)
      REAL P(3), Q(3)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This function computes the arc-length (angle in radians)
C between a pair of points on the unit sphere.
C
C On input:
C
C       P,Q = Arrays of length 3 containing the X, Y, and Z
C             coordinates (in that order) of points on the
C             unit sphere.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       ARCLEN = Angle in radians between the unit vectors
C                P and Q.  0 .LE. ARCLEN .LE. PI.
C
C Modules required by ARCLEN:  None
C
C Intrinsic functions called by ARCLEN:  ATAN, SQRT
C
C***********************************************************
C
      INTEGER I
      REAL    D
C
C Local parameters:
C
C D = Euclidean norm squared of P+Q
C I = DO-loop index
C
      D = 0.
      DO 1 I = 1,3
        D = D + (P(I) + Q(I))**2
    1   CONTINUE
      IF (D .EQ. 0.) THEN
C
C P and Q are separated by 180 degrees.
C
        ARCLEN = 4.*ATAN(1.)
      ELSEIF (D .GE. 4.) THEN
C
C P and Q coincide.
C
        ARCLEN = 0.
      ELSE
        ARCLEN = 2.*ATAN(SQRT((4.-D)/D))
      ENDIF
      RETURN
      END
      SUBROUTINE CONSTR (XK,YK,ZK, CX,SX,CY,SY)
      REAL XK, YK, ZK, CX, SX, CY, SY
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This subroutine constructs the elements of a 3 by 3
C orthogonal matrix R which rotates a point (XK,YK,ZK) on
C the unit sphere to the north pole, i.e.,
C
C      (XK)     (CY  0 -SY)   (1   0   0)   (XK)     (0)
C  R * (YK)  =  ( 0  1   0) * (0  CX -SX) * (YK)  =  (0)
C      (ZK)     (SY  0  CY)   (0  SX  CX)   (ZK)     (1)
C
C On input:
C
C       XK,YK,ZK = Components of a unit vector to be
C                  rotated to (0,0,1).
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       CX,SX,CY,SY = Elements of R:  CX,SX define a rota-
C                     tion about the X-axis and CY,SY define
C                     a rotation about the Y-axis.
C
C Modules required by CONSTR:  None
C
C Intrinsic function called by CONSTR:  SQRT
C
C***********************************************************
C
      CY = SQRT(YK*YK + ZK*ZK)
      SY = XK
      IF (CY .NE. 0.) THEN
        CX = ZK/CY
        SX = YK/CY
      ELSE
C
C (XK,YK,ZK) lies on the X-axis.
C
        CX = 1.
        SX = 0.
      ENDIF
      RETURN
      END
      REAL FUNCTION FVAL (B1,B2,B3,V1,V2,V3,F1,F2,F3,G1,G2,
     .                    G3,SIG1,SIG2,SIG3)
      REAL B1, B2, B3, V1(3), V2(3), V3(3), F1, F2, F3,
     .     G1(3), G2(3), G3(3), SIG1, SIG2, SIG3
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   Given data values and gradients at the three vertices of
C a spherical triangle containing a point P, this routine
C computes the value of F at P where F interpolates the ver-
C tex data.  Along the triangle sides, the interpolatory
C function F is the Hermite interpolatory tension spline
C defined by the values and tangential gradient components
C at the endpoints, and the gradient component normal to the
C triangle side varies linearly with respect to arc-length
C between the normal gradient components at the endpoints.
C A first-order C-1 blending method is used on the underly-
C ing planar triangle.  Since values and gradients on an arc
C depend only on the vertex data, the method results in C-1
C continuity when used to interpolate over a triangulation.
C
C   The blending method consists of taking F(P) to be a
C weighted sum of the values at PP of three univariate Her-
C mite interpolatory tension splines defined on the line
C segments which join the vertices to the opposite sides and
C pass through PP:  the central projection of P onto the
C underlying planar triangle.  The tension factors for these
C splines are obtained by linear interpolation between the
C pair of tension factors associated with the triangle sides
C which join at the appropriate vertex.
C
C   A tension factor SIGMA associated with a Hermite interp-
C olatory tension spline is a nonnegative parameter which
C determines the curviness of the spline.  SIGMA = 0 results
C in a cubic spline, and the spline approaches the linear
C interpolant as SIGMA increases.
C
C On input:
C
C       B1,B2,B3 = Barycentric coordinates of PP with re-
C                  spect to the (planar) underlying triangle
C                  (V1,V2,V3), where PP is the central
C                  projection of P onto this triangle.
C
C       V1,V2,V3 = Cartesian coordinates of the vertices of
C                  a spherical triangle containing P.  V3
C                  Left V1->V2.
C
C       F1,F2,F3 = Data values associated with the vertices.
C
C       G1,G2,G3 = Gradients associated with the vertices.
C                  Gi is orthogonal to Vi for i = 1,2,3.
C
C       SIG1,SIG2,SIG3 = Tension factors associated with the
C                        triangle sides opposite V1, V2, and
C                        V3, respectively.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       FVAL = Interpolated value at P.
C
C Each vector V above contains X, Y, and Z components in
C   V(1), V(2), and V(3), respectively.
C
C SSRFPACK modules required by FVAL:  ARCINT, ARCLEN, HVAL
C
C Intrinsic function called by FVAL:  SQRT
C
C***********************************************************
C
      REAL    HVAL
      INTEGER I
      REAL    C1, C2, C3, DS, DUM, DV, F, G(3),
     .        Q1(3), Q2(3), Q3(3), SIG, SUM, S1, S2, S3,
     .        U1(3), U2(3), U3(3), U1N, U2N, U3N, VAL
C
C Local parameters:
C
C C1,C2,C3 =    Coefficients (weight functions) of partial
C                 interpolants.  C1 = 1 on the edge opposite
C                 V1 and C1 = 0 on the other edges.  Simi-
C                 larly for C2 and C3.  C1+C2+C3 = 1.
C DS =          Directional derivative (scaled by distnace)
C                 at U1, U2, or U3:  DS = (G,U1-V1)/U1N =
C                 -(G,V1)/U1N on side opposite V1, where G/
C                 U1N (plus an orthogonal component) is the
C                 projection of G onto the planar triangle
C DUM =         Dummy variable for calls to ARCINT
C DV =          Directional derivatives (scaled by distance)
C                 at a vertex:  D1 = (G1,U1-V1) = (G1,U1)
C F,G =         Value and gradient at Q1 Q2, or Q3 obtained
C                 by interpolation along one of the arcs of
C                 the spherical triangle
C I =           DO-loop index
C Q1,Q2,Q3 =    Central projections of U1, U2, and U3 onto
C                 the sphere and thus lying on an arc of the
C                 spherical triangle
C SIG =         Tension factor for a side-vertex (partial)
C                 interpolant:  obtained by linear interpo-
C                 lation applied to triangle side tensions
C SUM =         Quantity used to normalize C1, C2, and C3
C S1,S2,S3 =    Sums of pairs of barycentric coordinates:
C                 used to compute U1, U2, U3, and SIG
C U1,U2,U3 =    Points on the boundary of the planar trian-
C                 gle and lying on the lines containing PP
C                 and the vertices.  U1 is opposite V1, etc.
C U1N,U2N,U3N = Quantities used to compute Q1, Q2, and Q3
C                 (magnitudes of U1, U2, and U3)
C VAL =         Local variable used to accumulate the con-
C                 tributions to FVAL
C
C
C Compute weight functions C1, C2, and C3.
C
      C1 = B2*B3
      C2 = B3*B1
      C3 = B1*B2
      SUM = C1 + C2 + C3
      IF (SUM .LE. 0.) THEN
C
C P coincides with a vertex.
C
        FVAL = B1*F1 + B2*F2 + B3*F3
        RETURN
      ENDIF
C
C Normalize C1, C2, and C3.
C
      C1 = C1/SUM
      C2 = C2/SUM
      C3 = C3/SUM
C
C Compute (S1,S2,S3), (U1,U2,U3) and (U1N,U2N,U3N).
C
      S1 = B2 + B3
      S2 = B3 + B1
      S3 = B1 + B2
      U1N = 0.
      U2N = 0.
      U3N = 0.
      DO 1 I = 1,3
        U1(I) = (B2*V2(I) + B3*V3(I))/S1
        U2(I) = (B3*V3(I) + B1*V1(I))/S2
        U3(I) = (B1*V1(I) + B2*V2(I))/S3
        U1N = U1N + U1(I)*U1(I)
        U2N = U2N + U2(I)*U2(I)
        U3N = U3N + U3(I)*U3(I)
    1   CONTINUE
C
C Compute Q1, Q2, and Q3.
C
      U1N = SQRT(U1N)
      U2N = SQRT(U2N)
      U3N = SQRT(U3N)
      DO 2 I = 1,3
        Q1(I) = U1(I)/U1N
        Q2(I) = U2(I)/U2N
        Q3(I) = U3(I)/U3N
    2   CONTINUE
C
C Compute interpolated value (VAL) at P by looping on
C   triangle sides.
C
      VAL = 0.
C
C Contribution from side opposite V1:
C
C   Compute value and gradient at Q1 by interpolating
C     between V2 and V3.
C
      CALL ARCINT (Q1,V2,V3,F2,F3,G2,G3,SIG1, F,G,DUM)
C
C   Add in the contribution.
C
      DV = G1(1)*U1(1) + G1(2)*U1(2) + G1(3)*U1(3)
      DS = -(G(1)*V1(1) + G(2)*V1(2) + G(3)*V1(3))/U1N
      SIG = (B2*SIG3 + B3*SIG2)/S1
      VAL = VAL + C1*HVAL(B1,F1,F,DV,DS,SIG)
C
C Contribution from side opposite V2:
C
C   Compute value and gradient at Q2 by interpolating
C     between V3 and V1.
C
      CALL ARCINT (Q2,V3,V1,F3,F1,G3,G1,SIG2, F,G,DUM)
C
C   Add in the contribution.
C
      DV = G2(1)*U2(1) + G2(2)*U2(2) + G2(3)*U2(3)
      DS = -(G(1)*V2(1) + G(2)*V2(2) + G(3)*V2(3))/U2N
      SIG = (B3*SIG1 + B1*SIG3)/S2
      VAL = VAL + C2*HVAL(B2,F2,F,DV,DS,SIG)
C
C Contribution from side opposite V3:
C
C   Compute interpolated value and gradient at Q3
C     by interpolating between V1 and V2.
C
      CALL ARCINT (Q3,V1,V2,F1,F2,G1,G2,SIG3, F,G,DUM)
C
C   Add in the final contribution.
C
      DV = G3(1)*U3(1) + G3(2)*U3(2) + G3(3)*U3(3)
      DS = -(G(1)*V3(1) + G(2)*V3(2) + G(3)*V3(3))/U3N
      SIG = (B1*SIG2 + B2*SIG1)/S3
      FVAL = VAL + C3*HVAL(B3,F3,F,DV,DS,SIG)
      RETURN
      END
      SUBROUTINE GETSIG (N,X,Y,Z,H,LIST,LPTR,LEND,GRAD,
     .                   TOL, SIGMA, DSMAX,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IER
      REAL    X(N), Y(N), Z(N), H(N), GRAD(3,N), TOL,
     .        SIGMA(*), DSMAX
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values H and gradients GRAD at the
C nodes, this subroutine determines, for each triangulation
C arc, the smallest (nonnegative) tension factor SIGMA such
C that the Hermite interpolatory tension spline H(A), de-
C fined by SIGMA and the endpoint values and directional
C derivatives, preserves local shape properties of the data.
C In order to define the shape properties on an arc, it is
C convenient to map the arc to an interval (A1,A2).  Then,
C denoting the endpoint data values by H1,H2 and the deriva-
C tives (tangential gradient components) by HP1,HP2, and
C letting S = (H2-H1)/(A2-A1), the data properties are
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
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
C           for I = 1,...,N.
C
C       LIST,LPTR,LEND = Data structure defining the tri-
C                        angulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       GRAD = Array dimensioned 3 by N whose columns con-
C              tain gradients at the nodes.  GRAD( ,J) must
C              be orthogonal to node J:  GRAD(1,J)*X(J) +
C              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0..  Refer
C              to Subroutines GRADG, GRADL, and SMSURF.
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
C               H(A) preserves the local data properties on
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
C STRIPACK modules required by GETSIG:  LSTPTR, STORE
C
C SSRFPACK modules required by GETSIG:  ARCLEN, SNHCSH
C
C Intrinsic functions called by GETSIG:  ABS, EXP, MAX, MIN,
C                                          SIGN, SQRT
C
C***********************************************************
C
      INTEGER LSTPTR
      REAL    ARCLEN, STORE
      INTEGER ICNT, LP1, LP2, LPL, LUN, N1, N2, NIT, NM1
      REAL    A, AL, C1, C2, COSHM, COSHMM, D0, D1, D1D2,
     .        D1PD2, D2, DMAX, DSIG, DSM, E, EMS, EMS2, F,
     .        F0, FMAX, FNEG, FP, FTOL, P1(3), P2(3), RTOL,
     .        S, S1, S2, SBIG, SCM, SGN, SIG, SIGIN, SINHM,
     .        SSINH, SSM, STOL, T, T0, T1, T2, TM, TP1,
     .        UN(3), UNORM
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
  100 FORMAT ('1',13X,'GETSIG -- N =',I4,', TOL = ',E10.3//)
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
C Print a message and compute parameters for the arc:
C   nodal coordinates P1 and P2, arc-length AL,
C   UNORM = magnitude of P1 X P2, and
C   SIGIN = input SIGMA value.
C
        IF (LUN .GE. 0) WRITE (LUN,110) N1, N2
  110   FORMAT (/1X,'ARC',I4,' -',I4)
        P1(1) = X(N1)
        P1(2) = Y(N1)
        P1(3) = Z(N1)
        P2(1) = X(N2)
        P2(2) = Y(N2)
        P2(3) = Z(N2)
        AL = ARCLEN(P1,P2)
        UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
        UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
        UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
        UNORM = SQRT(UN(1)*UN(1)+UN(2)*UN(2)+UN(3)*UN(3))
        IF (UNORM .EQ. 0.  .OR.  AL .EQ. 0.) GO TO 12
        SIGIN = SIGMA(LP1)
        IF (SIGIN .GE. SBIG) GO TO 9
C
C Compute scaled directional derivatives S1,S2 at the end-
C   points (for the direction N1->N2), first difference S,
C   and second differences D1,D2.
C
        S1 = AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +
     .               GRAD(3,N1)*P2(3))/UNORM
        S2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +
     .            GRAD(3,N2)*P1(3))/UNORM
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
C Convexity:  find a zero of F(SIG) = SIG*coshm(SIG)/
C   sinhm(SIG) - TP1.
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
C   Scale sinhm and coshm by 2*exp(-SIG) in order to avoid
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
  120   FORMAT (1X,'CONVEXITY -- SIG = ',E15.8,
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
  130   FORMAT (1X,'MONOTONICITY -- DSIG = ',E15.8)
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
C   Scale sinhm and coshm by 2*exp(-SIG) in order to avoid
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
C N < 3.
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
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
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
      IF (ABS(AA) .GT. ABS(BB)) THEN
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
      ELSEIF (BB .NE. 0.) THEN
C
C ABS(A) .LE. ABS(B).
C
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
      ELSE
C
C A = B = 0.
C
        C = 1.
        S = 0.
      ENDIF
      RETURN
      END
      SUBROUTINE GRADG (N,X,Y,Z,F,LIST,LPTR,LEND,IFLGS,
     .                  SIGMA, NIT,DGMAX,GRAD, IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, NIT, IER
      REAL    X(N), Y(N), Z(N), F(N), SIGMA(*), DGMAX,
     .        GRAD(3,N)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/24/96
C
C   Given a triangulation of N nodes on the unit sphere with
C data values F at the nodes and tension factors SIGMA asso-
C ciated with the arcs, this routine uses a global method
C to compute estimated gradients at the nodes.  The method
C consists of minimizing a quadratic functional Q(G) over
C the N-vector G of gradients, where Q approximates the
C linearized curvature of the restriction to arcs of the
C interpolatory function F defined by Function FVAL.  The
C restriction of F to an arc of the triangulation is the
C Hermite interpolatory tension spline defined by the data
C values and tangential gradient components at the endpoints
C of the arc.  Letting D1F(A) and D2F(A) denote first and
C second derivatives of F with respect to a parameter A var-
C ying along a triangulation arc, Q is the sum of integrals
C over the arcs of D2F(A)**2 + ((SIGMA/L)*(D1F(A)-S))**2,
C where L denotes arc-length, SIGMA is the appropriate ten-
C sion factor, and S is the slope of the linear function of
C A which interpolates the values of F at the endpoints of
C the arc.
C
C   Since the gradient at node K lies in the plane tangent
C to the sphere surface at K, it is effectively defined by
C only two components -- its X and Y components in the coor-
C dinate system obtained by rotating K to the north pole.
C Thus, the minimization problem corresponds to an order-2N
C symmetric positive-definite sparse linear system which is
C solved by a block Gauss-Seidel method with 2 by 2 blocks.
C
C   An alternative method, Subroutine GRADL, computes a
C local approximation to the gradient at a single node and,
C although less efficient when all gradients are needed, was
C found to be generally more accurate (in the case of uni-
C form zero tension) when the nodal distribution is very
C dense, varies greatly, or does not cover the sphere.
C GRADG, on the other hand, was found to be slightly more
C accurate on a uniform distribution of 514 nodes.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing Cartesian
C               coordinates of the nodes.  X(I)**2 + Y(I)**2
C               + Z(I)**2 = 1 for I = 1,...,N.
C
C       F = Array of length N containing data values at the
C           nodes.  F(I) is associated with (X(I),Y(I),Z(I))
C           for I = 1,...,N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
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
C               programs GETSIG, SIG0, SIG1, and SIG2.
C
C The above parameters are not altered by this routine.
C
C       NIT = Maximum number of Gauss-Seidel iterations to
C             be applied.  This maximum will likely be a-
C             chieved if DGMAX is smaller than the machine
C             precision.  NIT .GE. 0.
C
C       DGMAX = Nonnegative convergence criterion.  The
C               method is terminated when the maximum change
C               in a gradient between iterations is at most
C               DGMAX.  The change in a gradient is taken to
C               be the Euclidean norm of the difference (in
C               the rotated coordinate system) relative to 1
C               plus the norm of the old gradient value.
C
C       GRAD = 3 by N array whose columns contain initial
C              solution estimates (zero vectors are suffici-
C              ent).  GRAD(I,J) contains component I of the
C              gradient at node J for I = 1,2,3 (X,Y,Z) and
C              J = 1,...,N.  GRAD( ,J) must be orthogonal to
C              node J -- GRAD(1,J)*X(J) + GRAD(2,J)*Y(J) +
C              GRAD(3,J)*Z(J) = 0.
C
C On output:
C
C       NIT = Number of Gauss-Seidel iterations employed.
C
C       DGMAX = Maximum change in a gradient at the last
C               iteration.
C
C       GRAD = Estimated gradients.  See the description
C              under input parameters.  GRAD is not changed
C              if IER = -1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     convergence criterion was achieved.
C             IER = 1 if no errors were encountered but con-
C                     vergence was not achieved within NIT
C                     iterations.
C             IER = -1 if N or DGMAX is outside its valid
C                      range or NIT .LT. 0 on input.
C             IER = -2 if all nodes are collinear or the
C                      triangulation is invalid.
C             IER = -3 if duplicate nodes were encountered.
C
C SSRFPACK modules required by GRADG:  APLYRT, CONSTR,
C                                        GRCOEF, SNHCSH
C
C Intrinsic functions called by GRADG:  ATAN, MAX, SQRT
C
C***********************************************************
C
      INTEGER IFL, ITER, J, K, LPJ, LPL, MAXIT, NN
      REAL    ALFA, A11, A12, A22, CX, CY, D, DEN, DET,
     .        DGK(3), DGMX, DG1, DG2, FK, G1, G2, G3, R1,
     .        R2, SD, SIG, SINAL, SX, SY, T, TOL, XK, YK,
     .        ZK, XJ, YJ, ZJ, XS, YS
C
C Local parameters:
C
C ALFA =        Arc-length between nodes K and J
C A11,A12,A22 = Matrix components of the 2 by 2 block A*DG
C                 = R where A is symmetric, (DG1,DG2,0) is
C                 the change in the gradient at K, and R is
C                 the residual
C CX,CY =       Components of a rotation mapping K to the
C                 north pole (0,0,1)
C D =           Function of SIG computed by GRCOEF -- factor
C                 in the order-2 system
C DEN =         ALFA*SINAL**2 -- factor in the 2 by 2 system
C DET =         Determinant of the order-2 matrix
C DGK =         Change in GRAD( ,K) from the previous esti-
C                 mate in the original coordinate system
C DGMX =        Maximum change in a gradient between itera-
C                 tions
C DG1,DG2 =     Solution of the 2 by 2 system -- first 2
C                 components of DGK in the rotated coordi-
C                 nate system
C FK =          Data value F(K)
C G1,G2,G3 =    Components of GRAD( ,K)
C IFL =         Local copy of IFLGS
C ITER =        Number of iterations used
C J =           Neighbor of K
C K =           DO-loop and node index
C LPJ =         LIST pointer of node J as a neighbor of K
C LPL =         Pointer to the last neighbor of K
C MAXIT =       Input value of NIT
C NN =          Local copy of N
C R1,R2 =       Components of the residual -- derivatives of
C                 Q with respect to the components of the
C                 gradient at node K
C SD =          Function of SIG computed by GRCOEF -- factor
C                 in the order-2 system
C SIG =         Tension factor associated with ARC K-J
C SINAL =       SIN(ALFA) -- magnitude of the vector cross
C                 product between nodes K and J
C SX,SY =       Components of a rotation mapping K to the
C                 north pole (0,0,1)
C T =           Temporary storage for factors in the system
C                 components
C TOL =         Local copy of DGMAX
C XK,YK,ZK =    Coordinates of node K -- X(K), Y(K), Z(K)
C XJ,YJ,ZJ =    Coordinates of node J in the rotated coor-
C                 dinate system
C XS,YS =       XJ**2, YJ**2
C
      NN = N
      IFL = IFLGS
      MAXIT = NIT
      TOL = DGMAX
C
C Test for errors in input, and initialize iteration count,
C   tension factor, and output value of DGMAX.
C
      IF (NN .LT. 3  .OR.  MAXIT .LT. 0  .OR.  TOL .LT. 0.)
     .   GO TO 11
      ITER = 0
      SIG = SIGMA(1)
      DGMX = 0.
C
C Top of iteration loop.
C
    1 IF (ITER .EQ. MAXIT) GO TO 4
      DGMX = 0.
C
C Loop on nodes.
C
      DO 3 K = 1,NN
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        FK = F(K)
        G1 = GRAD(1,K)
        G2 = GRAD(2,K)
        G3 = GRAD(3,K)
C
C   Construct the rotation mapping node K to the north pole.
C
        CALL CONSTR (XK,YK,ZK, CX,SX,CY,SY)
C
C   Initialize components of the 2 by 2 system for the
C     change (DG1,DG2,0) in the K-th solution components
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
    2   LPJ = LPTR(LPJ)
          J = ABS(LIST(LPJ))
C
C   Compute the coordinates of J in the rotated system.
C
          T = SX*Y(J) + CX*Z(J)
          YJ = CX*Y(J) - SX*Z(J)
          ZJ = SY*X(J) + CY*T
          XJ = CY*X(J) - SY*T
C
C   Compute arc-length ALFA between J and K, SINAL =
C     SIN(ALFA), and DEN = ALFA*SIN(ALFA)**2.
C
          ALFA = 2.*ATAN(SQRT((1.-ZJ)/(1.+ZJ)))
          XS = XJ*XJ
          YS = YJ*YJ
          SINAL = SQRT(XS+YS)
          DEN = ALFA*(XS+YS)
C
C   Test for coincident nodes and compute functions of SIG:
C     D = SIG*(SIG*COSHM-SINHM)/E and SD = SIG*SINHM/E for
C     E = SIG*SINH-2*COSHM.
C
          IF (DEN .EQ. 0.) GO TO 13
          IF (IFL .GE. 1) SIG = SIGMA(LPJ)
          CALL GRCOEF (SIG, D,SD)
C
C   Update the system components for node J.
C
          T = D/DEN
          A11 = A11 + T*XS
          A12 = A12 + T*XJ*YJ
          A22 = A22 + T*YS
          T = (D+SD)*(FK-F(J))/(ALFA*ALFA*SINAL) +
     .        ( D*(G1*X(J) + G2*Y(J) + G3*Z(J)) -
     .          SD*(GRAD(1,J)*XK + GRAD(2,J)*YK +
     .                 GRAD(3,J)*ZK) )/DEN
          R1 = R1 - T*XJ
          R2 = R2 - T*YJ
C
C   Bottom of loop on neighbors.
C
          IF (LPJ .NE. LPL) GO TO 2
C
C   Solve the 2 by 2 system and update DGMAX.
C
        DET = A11*A22 - A12*A12
        IF (DET .EQ. 0.  .OR.  A11 .EQ. 0.) GO TO 12
        DG2 = (A11*R2 - A12*R1)/DET
        DG1 = (R1 - A12*DG2)/A11
        DGMX = MAX(DGMX,SQRT(DG1*DG1+DG2*DG2)/
     .             (1.+SQRT(G1*G1+G2*G2+G3*G3)))
C
C   Rotate (DG1,DG2,0) back to the original coordinate
C     system and update GRAD( ,K).
C
        CALL APLYRT (DG1,DG2,CX,SX,CY,SY, DGK)
        GRAD(1,K) = G1 + DGK(1)
        GRAD(2,K) = G2 + DGK(2)
        GRAD(3,K) = G3 + DGK(3)
    3   CONTINUE
C
C   Increment ITER and test for convergence.
C
      ITER = ITER + 1
      IF (DGMX .GT. TOL) GO TO 1
C
C The method converged.
C
      NIT = ITER
      DGMAX = DGMX
      IER = 0
      RETURN
C
C The method failed to converge within NIT iterations.
C
    4 DGMAX = DGMX
      IER = 1
      RETURN
C
C Invalid input parameter.
C
   11 NIT = 0
      DGMAX = 0.
      IER = -1
      RETURN
C
C Node K and its neighbors are collinear.
C
   12 NIT = 0
      DGMAX = DGMX
      IER = -2
      RETURN
C
C Nodes K and J coincide.
C
   13 NIT = 0
      DGMAX = DGMX
      IER = -3
      RETURN
      END
      SUBROUTINE GRADL (N,K,X,Y,Z,W,LIST,LPTR,LEND, G,IER)
      INTEGER N, K, LIST(*), LPTR(*), LEND(N), IER
      REAL    X(N), Y(N), Z(N), W(N), G(3)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/24/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere with their associated data values W, this routine
C estimates a gradient vector at node K as follows:  the
C coordinate system is rotated so that K becomes the north
C pole, node K and a set of nearby nodes are projected
C orthogonally onto the X-Y plane (in the new coordinate
C system), a quadratic is fitted in a weighted least squares
C sense to the data values at the projected nodes such that
C the value (associated with K) at (0,0) is interpolated, X
C and Y-partial derivative estimates DX and DY are computed
C by differentiating the quadratic at (0,0), and the esti-
C mated gradient G is obtained by rotating (DX,DY,0) back to
C the original coordinate system.  Note that G lies in the
C plane tangent to the sphere at node K, i.e., G is orthogo-
C nal to the unit vector represented by node K.  A Marquardt
C stabilization factor is used if necessary to ensure a
C well-conditioned least squares system, and a unique solu-
C tion exists unless the nodes are collinear.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 7.
C
C       K = Node at which the gradient is sought.  1 .LE. K
C           .LE. N.
C
C       X,Y,Z = Arrays containing the Cartesian coordinates
C               of the nodes.
C
C       W = Array containing the data values at the nodes.
C           W(I) is associated with (X(I),Y(I),Z(I)) for
C           I = 1,...,N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       G = X, Y, and Z components (in that order) of the
C           estimated gradient at node K unless IER < 0.
C
C       IER = Error indicator:
C             IER .GE. 6 if no errors were encountered.
C                        IER contains the number of nodes
C                        (including K) used in the least
C                        squares fit.
C             IER = -1 if N or K is outside its valid range.
C             IER = -2 if the least squares system has no
C                      unique solution due to duplicate or
C                      collinear nodes.
C
C STRIPACK module required by GRADL:  GETNP
C
C SSRFPACK modules required by GRADL:  APLYR, APLYRT,
C                                        CONSTR, GIVENS,
C                                        ROTATE, SETUP
C
C Intrinsic functions called by GRADL:  ABS, MIN, REAL, SQRT
C
C***********************************************************
C
      INTEGER   LMN, LMX
      PARAMETER (LMN=10,  LMX=30)
      INTEGER I, IERR, IM1, IP1, J, JP1, KK, L, LM1, LMAX,
     .        LMIN, LNP, NN, NP, NPTS(LMX)
      REAL    A(6,6), AV, AVSQ, C, CX, CY, DF, DMIN, DTOL,
     .        DX, DY, RF, RIN, RTOL, S, SF, SUM, SX, SY,
     .        WK, WT, XP, YP, ZP
C
      DATA    RTOL/1.E-6/, DTOL/.01/, SF/1./
C
C Local parameters:
C
C A =         Transpose of the (upper triangle of the) aug-
C               mented regression matrix
C AV =        Root-mean-square distance (in the rotated
C               coordinate system) between the origin and
C               the nodes (other than K) in the least
C               squares fit.  The first 3 columns of A**T
C               are scaled by 1/AVSQ, the next 2 by 1/AV.
C AVSQ =      AV*AV:  accumulated in SUM
C C,S =       Components of the plane rotation used to
C               triangularize the regression matrix
C CX,SX =     Components of a plane rotation about the X-
C               axis which, together with CY and SY, define
C               a mapping from node K to the north pole
C               (0,0,1)
C CY,SY =     Components of a plane rotation about the Y-
C               axis
C DF =        Negative Z component (in the rotated coordi-
C               nate system) of an element NP of NPTS --
C               increasing function of the angular distance
C               between K and NP.  DF lies in the interval
C               (-1,1).
C DMIN =      Minimum of the magnitudes of the diagonal
C               elements of the triangularized regression
C               matrix
C DTOL =      Tolerance for detecting an ill-conditioned
C               system (DMIN is required to be at least
C               DTOL)
C DX,DY =     X and Y components of the estimated gradient
C               in the rotated coordinate system
C I,J =       Loop indexes
C IERR =      Error flag for calls to GETNP (not checked)
C IM1,IP1 =   I-1, I+1
C JP1 =       J+1
C KK =        Local copy of K
C L =         Number of columns of A**T to which a rotation
C               is applied
C LM1 =       LMIN-1
C LMIN,LMAX = Min(LMN,N), Min(LMX,N)
C LMN,LMX =   Minimum and maximum values of LNP for N
C               sufficiently large.  In most cases LMN-1
C               nodes are used in the fit.  7 .LE. LMN .LE.
C               LMX.
C LNP =       Length of NPTS or LMAX+1
C NN =        Local copy of N
C NP =        Element of NPTS to be added to the system
C NPTS =      Array containing the indexes of a sequence of
C               nodes ordered by angular distance from K.
C               NPTS(1)=K and the first LNP-1 elements of
C               NPTS are used in the least squares fit.
C               unless LNP = LMAX+1, NPTS(LNP) determines R
C               (see RIN).
C RF =        Value of DF associated with NPTS(LNP) unless
C               LNP = LMAX+1 (see RIN)
C RIN =       Inverse of a radius of influence R which
C               enters into WT:  R = 1+RF unless all ele-
C               ments of NPTS are used in the fit (LNP =
C               LMAX+1), in which case R is the distance
C               function associated with some point more
C               distant from K than NPTS(LMAX)
C RTOL =      Tolerance for determining LNP (and hence R):
C               if the increase in DF between two successive
C               elements of NPTS is less than RTOL, they are
C               treated as being the same distance from node
C               K and an additional node is added
C SF =        Marquardt stabilization factor used to damp
C               out the first 3 solution components (second
C               partials of the quadratic) when the system
C               is ill-conditioned.  Increasing SF results
C               in more damping (a more nearly linear fit).
C SUM =       Sum of squared Euclidean distances (in the
C               rotated coordinate system) between the
C               origin and the nodes used in the least
C               squares fit
C WK =        W(K) -- data value at node K
C WT =        Weight for the equation coreesponding to NP:
C               WT = (R-D)/(R*D) = 1/D - RIN, where D = 1-ZP
C               is associated with NP
C XP,YP,ZP =  Coordinates of NP in the rotated coordinate
C               system unless ZP < 0, in which case
C               (XP,YP,0) lies on the equator
C
      NN = N
      KK = K
      WK = W(KK)
C
C Check for errors and initialize LMIN, LMAX.
C
      IF (NN .LT. 7  .OR.  KK .LT. 1  .OR.  KK .GT. NN)
     .   GO TO 13
      LMIN = MIN(LMN,NN)
      LMAX = MIN(LMX,NN)
C
C Compute NPTS, LNP, AVSQ, AV, and R.
C   Set NPTS to the closest LMIN-1 nodes to K.  DF contains
C   the negative Z component (in the rotated coordinate
C   system) of the new node on return from GETNP.
C
      SUM = 0.
      NPTS(1) = KK
      LM1 = LMIN - 1
      DO 1 LNP = 2,LM1
        CALL GETNP (X,Y,Z,LIST,LPTR,LEND,LNP, NPTS, DF,IERR)
        SUM = SUM + 1. - DF*DF
    1   CONTINUE
C
C   Add additional nodes to NPTS until the increase in
C     R = 1+RF is at least RTOL.
C
      DO 2 LNP = LMIN,LMAX
        CALL GETNP (X,Y,Z,LIST,LPTR,LEND,LNP, NPTS, RF,IERR)
        IF (RF-DF .GE. RTOL) GO TO 3
        SUM = SUM + 1. - RF*RF
    2   CONTINUE
C
C   Use all LMAX nodes in the least squares fit.  R is
C     arbitrarily increased by 5 percent.
C
      RF = 1.05*RF + .05
      LNP = LMAX + 1
C
C   There are LNP-2 equations corresponding to nodes
C     NPTS(2),...,NPTS(LNP-1).
C
    3 AVSQ = SUM/REAL(LNP-2)
      AV = SQRT(AVSQ)
      RIN = 1./(1.+RF)
C
C Construct the rotation.
C
      CALL CONSTR (X(KK),Y(KK),Z(KK), CX,SX,CY,SY)
C
C Set up the first 5 equations of the augmented regression
C   matrix (transposed) as the columns of A, and zero out
C   the lower triangle (upper triangle of A) with Givens
C   rotations.
C
      DO 5 I = 1,5
        NP = NPTS(I+1)
        CALL APLYR (X(NP),Y(NP),Z(NP),CX,SX,CY,SY, XP,YP,ZP)
        WT = 1./(1.-ZP) - RIN
        CALL SETUP (XP,YP,W(NP),WK,AV,AVSQ,WT, A(1,I))
        IF (I .EQ. 1) GO TO 5
        IM1 = I - 1
        DO 4 J = 1,IM1
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS ( A(J,J),A(J,I), C,S)
          CALL ROTATE (L,C,S, A(JP1,J),A(JP1,I) )
    4     CONTINUE
    5   CONTINUE
C
C Add the additional equations to the system using
C   the last column of A.  I .LE. LNP.
C
      I = 7
    6   IF (I .EQ. LNP) GO TO 8
        NP = NPTS(I)
        CALL APLYR (X(NP),Y(NP),Z(NP),CX,SX,CY,SY, XP,YP,ZP)
        WT = 1./(1.-ZP) - RIN
        CALL SETUP (XP,YP,W(NP),WK,AV,AVSQ,WT, A(1,6))
        DO 7 J = 1,5
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS ( A(J,J),A(J,6), C,S)
          CALL ROTATE (L,C,S, A(JP1,J),A(JP1,6) )
    7     CONTINUE
        I = I + 1
        GO TO 6
C
C Test the system for ill-conditioning.
C
    8 DMIN = MIN( ABS(A(1,1)),ABS(A(2,2)),ABS(A(3,3)),
     .            ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN .GE. DTOL) GO TO 12
      IF (LNP .LE. LMAX) THEN
C
C Add another node to the system and increase R.
C   I = LNP.
C
        LNP = LNP + 1
        IF (LNP .LE. LMAX) CALL GETNP (X,Y,Z,LIST,LPTR,LEND,
     .                                 LNP,NPTS, RF,IERR)
        RIN = 1./(1.05*(1.+RF))
        GO TO 6
      ENDIF
C
C Stabilize the system by damping second partials.  Add
C   multiples of the first three unit vectors to the first
C   three equations.
C
      DO 11 I = 1,3
        A(I,6) = SF
        IP1 = I + 1
        DO 9 J = IP1,6
          A(J,6) = 0.
    9     CONTINUE
        DO 10 J = I,5
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS ( A(J,J),A(J,6), C,S)
          CALL ROTATE (L,C,S, A(JP1,J),A(JP1,6) )
   10     CONTINUE
   11   CONTINUE
C
C Test the linear portion of the stabilized system for
C   ill-conditioning.
C
      DMIN = MIN( ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN .LT. DTOL) GO TO 14
C
C Solve the 2 by 2 triangular system for the estimated
C   partial derivatives.
C
   12 DY = A(6,5)/A(5,5)
      DX = (A(6,4) - A(5,4)*DY)/A(4,4)/AV
      DY = DY/AV
C
C Rotate the gradient (DX,DY,0) back into the original
C   coordinate system.
C
      CALL APLYRT (DX,DY,CX,SX,CY,SY, G)
      IER = LNP - 1
      RETURN
C
C N or K is outside its valid range.
C
   13 IER = -1
      RETURN
C
C No unique solution due to collinear nodes.
C
   14 IER = -2
      RETURN
      END
      SUBROUTINE GRCOEF (SIGMA, D,SD)
      REAL SIGMA, D, SD
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   This subroutine computes factors involved in the linear
C systems solved by Subroutines GRADG and SMSGS.
C
C On input:
C
C       SIGMA = Nonnegative tension factor associated with a
C               triangulation arc.
C
C SIGMA is not altered by this routine.
C
C On output:
C
C       D = Diagonal factor.  D = SIG*(SIG*Coshm(SIG) -
C           Sinhm(SIG))/E where E = SIG*Sinh(SIG) - 2*
C           Coshm(SIG).  D > 0, and D = 4 at SIG = 0.
C
C       SD = Off-diagonal factor.  SD = SIG*Sinhm(SIG)/E.
C            SD > 0, and SD = 2 at SIG = 0.
C
C SSRFPACK module required by GRCOEF:  SNHCSH
C
C Intrinsic function called by GRCOEF:  EXP
C
C***********************************************************
C
      REAL COSHM, COSHMM, E, EMS, SCM, SIG, SINHM, SSINH,
     .     SSM
      SIG = SIGMA
      IF (SIG .LT. 1.E-9) THEN
C
C Cubic function:
C
        D = 4.
        SD = 2.
      ELSEIF (SIG .LE. .5) THEN
C
C 0 < SIG .LE. .5.
C
C Use approximations designed to avoid cancellation error
C   in the hyperbolic functions when SIGMA is small.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        E = SIG*SINHM - COSHMM - COSHMM
        D = SIG*(SIG*COSHM-SINHM)/E
        SD = SIG*SINHM/E
      ELSE
C
C SIG > .5.
C
C Scale SINHM, COSHM, and E by 2*EXP(-SIG) in order to
C   avoid overflow when SIGMA is large.
C
        EMS = EXP(-SIG)
        SSINH = 1. - EMS*EMS
        SSM = SSINH - 2.*SIG*EMS
        SCM = (1.-EMS)*(1.-EMS)
        E = SIG*SSINH - SCM - SCM
        D = SIG*(SIG*SCM-SSM)/E
        SD = SIG*SSM/E
      ENDIF
      RETURN
      END
      REAL FUNCTION HVAL (B,H1,H2,HP1,HP2,SIGMA)
      REAL B, H1, H2, HP1, HP2, SIGMA
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   Given a line segment P1-P2 containing a point P, along
C with values and derivatives at the endpoints, this func-
C tion returns the value H(P), where H is the Hermite inter-
C polatory tension spline defined by the endpoint data.
C
C On input:
C
C       B = Local coordinate of P with respect to P1-P2:
C           P = B*P1 + (1-B)*P2, and thus B = d(P,P2)/
C           d(P1,P2), where d(P1,P2) is the distance between
C           P1 and P2.  B < 0 or B > 1 results in extrapola-
C           tion.
C
C       H1,H2 = Values interpolated at P1 and P2, respec-
C               tively.
C
C       HP1,HP2 = Products of d(P1,P2) with first order der-
C                 ivatives at P1 and P2, respectively.  HP1
C                 may, for example, be the scalar product of
C                 P2-P1 with a gradient at P1.
C
C       SIGMA = Nonnegative tension factor associated with
C               the spline.  SIGMA = 0 corresponds to a
C               cubic spline, and H approaches the linear
C               interpolant of H1 and H2 as SIGMA increases.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       HVAL = Interpolated value H(P).
C
C SSRFPACK module required by HVAL:  SNHCSH
C
C Intrinsic functions called by HVAL:  ABS, EXP
C
C***********************************************************
C
      REAL B1, B2, CM, CM2, CMM, D1, D2, DUMMY, E, E1, E2,
     .     EMS, S, SB1, SB2, SIG, SM, SM2, TM, TM1, TM2, TS
      B1 = B
      B2 = 1. - B1
C
C Compute slope S and second differences D1 and D2 scaled
C   by the separation between P1 and P2.
C
      S = H2 - H1
      D1 = S - HP1
      D2 = HP2 - S
C
C Test the range of SIGMA.
C
      SIG = ABS(SIGMA)
      IF (SIG .LT. 1.E-9) THEN
C
C Hermite cubic interpolation:
C
        HVAL = H1 + B2*(HP1 + B2*(D1 + B1*(D1 - D2)))
      ELSEIF (SIG .LE. .5) THEN
C
C 0 < SIG .LE. .5.  Use approximations designed to avoid
C   cancellation error in the hyperbolic functions.
C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        E = SIG*SM - CMM - CMM
        HVAL = H1 + B2*HP1 + ((CM*SM2-SM*CM2)*(D1+D2) + SIG*
     .                     (CM*CM2-(SM+SIG)*SM2)*D1)/(SIG*E)
      ELSE
C
C SIG > .5.  Use negative exponentials in order to avoid
C   overflow.  Note that EMS = EXP(-SIG).
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        E1 = EXP(-SB1)
        E2 = EXP(-SB2)
        EMS = E1*E2
        TM = 1. - EMS
        TS = TM*TM
        TM1 = 1. - E1
        TM2 = 1. - E2
        E = TM*(SIG*(1.+EMS) - TM - TM)
        HVAL = H1 + B2*S + (TM*TM1*TM2*(D1+D2) + SIG*
     .                      ((E2*TM1*TM1-B1*TS)*D1 +
     .                       (E1*TM2*TM2-B2*TS)*D2))/(SIG*E)
      ENDIF
      RETURN
      END
      SUBROUTINE INTRC0 (N,PLAT,PLON,X,Y,Z,W,LIST,LPTR,
     .                   LEND, IST, PW,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IST, IER
      REAL    PLAT, PLON, X(N), Y(N), Z(N), W(N), PW
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/24/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values at the nodes, this sub-
C routine computes the value at a point P of a continuous
C function which interpolates the data values.  The interp-
C olatory function is linear on each underlying triangle
C (planar triangle with the same vertices as a spherical
C triangle).  If P is not contained in a triangle, an ex-
C trapolated value is taken to be the interpolated value at
C the nearest point of the triangulation boundary.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       PLAT,PLON = Latitude and longitude of P in radians.
C
C       X,Y,Z = Arrays containing Cartesian coordinates of
C               the nodes.
C
C       W = Array containing data values at the nodes.  W(I)
C           is associated with (X(I),Y(I),Z(I)) for I =
C           1,...,N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       IST = Index of the starting node in the search for a
C             triangle containing P.  1 .LE. IST .LE. N.
C             The output value of IST from a previous call
C             may be a good choice.
C
C Input parameters other than IST are not altered by this
C   routine.
C
C On output:
C
C       IST = Index of one of the vertices of the triangle
C             containing P (or nearest P) unless IER = -1
C             or IER = -2.
C
C       PW = Value of the interpolatory function at P if
C            IER .GE. 0.
C
C       IER = Error indicator:
C             IER = 0 if interpolation was performed
C                     successfully.
C             IER = 1 if extrapolation was performed
C                     successfully.
C             IER = -1 if N < 3 or IST is outside its valid
C                      range.
C             IER = -2 if the nodes are collinear.
C             IER = -3 if P is not in a triangle and the
C                      angle between P and the nearest boun-
C                      dary point is at least 90 degrees.
C
C STRIPACK modules required by INTRC0:  JRAND, LSTPTR,
C                                         STORE, TRFIND
C
C Intrinsic functions called by INTRC0:  COS, SIN
C
C***********************************************************
C
      INTEGER I1, I2, I3, LP, N1, N2
      REAL    B1, B2, B3, P(3), PTN1, PTN2, S12, SUM
C
C Local parameters:
C
C B1,B2,B3 = Barycentric coordinates of the central projec-
C              tion of P onto the underlying planar trian-
C              gle, or (B1 and B2) projection of Q onto the
C              underlying line segment N1-N2 when P is
C              exterior.  Unnormalized coordinates are
C              computed by TRFIND when P is in a triangle.
C I1,I2,I3 = Vertex indexes returned by TRFIND
C LP =       LIST pointer to N1 as a neighbor of N2 or N2
C              as a neighbor of N1
C N1,N2 =    Endpoints of a boundary arc which is visible
C              from P when P is not contained in a triangle
C P =        Cartesian coordinates of P
C PTN1 =     Scalar product (P,N1)
C PTN2 =     Scalar product (P,N2)
C S12 =      Scalar product (N1,N2)
C SUM =      Quantity used to normalize the barycentric
C              coordinates
C
      IF (N .LT. 3  .OR.  IST .LT. 1  .OR.  IST .GT. N)
     .    GO TO 11
C
C Transform (PLAT,PLON) to Cartesian coordinates.
C
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
C
C Find the vertex indexes of a triangle containing P.
C
      CALL TRFIND(IST,P,N,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,
     .            I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 12
      IST = I1
      IF (I3 .NE. 0) THEN
C
C P is contained in the triangle (I1,I2,I3).  Normalize the
C   barycentric coordinates.
C
        SUM = B1 + B2 + B3
        B1 = B1/SUM
        B2 = B2/SUM
        B3 = B3/SUM
        PW = B1*W(I1) + B2*W(I2) + B3*W(I3)
        IER = 0
        RETURN
      ENDIF
C
C P is exterior to the triangulation, and I1 and I2 are
C   boundary nodes which are visible from P.  Set PW to the
C   interpolated value at Q, where Q is the closest boundary
C   point to P.
C
C Traverse the boundary starting from the rightmost visible
C   node I1.
C
      N1 = I1
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
      IF (I1 .NE. I2) GO TO 2
C
C All boundary nodes are visible from P.  Find a boundary
C   arc N1->N2 such that P Left (N2 X N1)->N1.
C
C Counterclockwise boundary traversal:
C   Set N2 to the first neighbor of N1.
C
    1 LP = LEND(N1)
        LP = LPTR(LP)
        N2 = LIST(LP)
C
C Compute inner products (N1,N2) and (P,N2), and compute
C   B2 = DET(P,N1,N2 X N1).
C
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN2 = P(1)*X(N2) + P(2)*Y(N2) + P(3)*Z(N2)
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 2
C
C P Right (N2 X N1)->N1 -- Iterate.
C
        N1 = N2
        I1 = N1
        PTN1 = PTN2
        GO TO 1
C
C P Left (N2 X N1)->N1, where N2 is the first neighbor of P1.
C   Clockwise boundary traversal:
C
    2 N2 = N1
        PTN2 = PTN1
C
C Set N1 to the last neighbor of N2 and test for
C   termination.
C
        LP = LEND(N2)
        N1 = -LIST(LP)
        IF (N1 .EQ. I1) GO TO 13
C
C Compute inner products (N1,N2) and (P,N1).
C
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
C
C Compute B2 = DET(P,N1,N2 X N1) = DET(Q,N1,N2 X N1)*(P,Q).
C
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 2
C
C Compute B1 = DET(P,N2 X N1,N2) = DET(Q,N2 X N1,N2)*(P,Q).
C
      B1 = PTN1 - S12*PTN2
      IF (B1 .LE. 0.) THEN
C
C Q = N2.
C
        PW = W(N2)
      ELSE
C
C P Strictly Left (N2 X N1)->N2 and P Strictly Left
C   N1->(N2 X N1).  Thus Q lies on the interior of N1->N2.
C   Normalize the coordinates and compute PW.
C
        SUM = B1 + B2
        PW = (B1*W(N1) + B2*W(N2))/SUM
      ENDIF
      IER = 1
      RETURN
C
C N or IST is outside its valid range.
C
   11 IER = -1
      RETURN
C
C Collinear nodes.
C
   12 IER = -2
      RETURN
C
C The angular distance between P and the closest boundary
C   point to P is at least 90 degrees.
C
   13 IER = -3
      RETURN
      END
      SUBROUTINE INTRC1 (N,PLAT,PLON,X,Y,Z,F,LIST,LPTR,LEND,
     .                   IFLGS,SIGMA,IFLGG,GRAD, IST, FP,
     .                   IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, IFLGG,
     .        IST, IER
      REAL    PLAT, PLON, X(N), Y(N), Z(N), F(N), SIGMA(*),
     .        GRAD(3,N), FP
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/25/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values and gradients at the nodes,
C this routine computes a value F(P), where F interpolates
C the nodal data and is once-continuously differentiable
C over the convex hull of the nodes.  Refer to Function FVAL
C for further details.  If P is not contained in a triangle,
C an extrapolated value is computed by extending F beyond
C the boundary in a continuous fashion.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3
C           and N .GE. 7 if IFLGG .LE. 0.
C
C       PLAT,PLON = Latitude and longitude in radians of the
C                   point P at which F is to be evaluated.
C
C       X,Y,Z = Arrays of length N containing Cartesian
C               coordinates of the nodes.
C
C       F = Array of length N containing values of F at the
C           nodes:  F(I) = F(X(I),Y(I),Z(I)).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
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
C               programs FVAL, GETSIG, SIG0, SIG1, and SIG2.
C
C       IFLGG = Gradient option:
C               IFLGG .LE. 0 if INTRC1 is to provide grad-
C                            ient estimates as needed (from
C                            GRADL).
C               IFLGG .GE. 1 if gradients are user-provided
C                            in GRAD.  This is more effici-
C                            ent if INTRC1 is to be called
C                            several times.
C
C       GRAD = 3 by N array whose I-th column contains
C              an estimated gradient at node I if IFLGG .GE.
C              1, or unused dummy parameter if IFLGG .LE. 0.
C              Refer to Subroutines GRADL and GRADG.
C
C       IST = Index of the starting node in the search for a
C             triangle containing P.  The output value of
C             IST from a previous call may be a good choice.
C             1 .LE. IST .LE. N.
C
C Input parameters other than IST are not altered by this
C   routine.
C
C On output:
C
C       IST = Index of one of the vertices of the triangle
C             containing P (or a boundary node if P is not
C             contained in a triangle) unless IER = -1 or
C             IER = -2.
C
C       FP = Value of F at P unless IER < 0, in which case
C            FP is not defined.
C
C       IER = Error indicator and information flag:
C             IER = 0 if no errors were encountered and P is
C                     contained in a triangle.
C             IER = 1 if no errors were encountered and
C                     extrapolation was required.
C             IER = -1 if N or IST is outside its valid
C                      range.
C             IER = -2 if the nodes are collinear.
C             IER = -3 if the angular distance between P and
C                      the nearest point of the triangula-
C                      tion is at least 90 degrees.
C
C STRIPACK modules required by INTRC1: JRAND, LSTPTR, STORE,
C                                        TRFIND
C                    (and optionally)  GETNP if IFLGG .LE. 0
C
C SSRFPACK modules required by INTRC1:  ARCINT, ARCLEN,
C                                         FVAL, HVAL, SNHCSH
C              (and if IFLGG .LE. 0)  APLYR, APLYRT, CONSTR,
C                                       GIVENS, GRADL,
C                                       ROTATE, SETUP
C
C Intrinsic functions called by INTRC1:  COS, SIN, SQRT
C
C***********************************************************
C
      INTEGER LSTPTR
      REAL    ARCLEN, FVAL
      INTEGER I, IERR, I1, I2, I3, LP, N1, N2, NN
      REAL    A, B1, B2, B3, DUM(3), FQ, GQ(3), GQN, G1(3),
     .        G2(3), G3(3), P(3), P1(3), P2(3), P3(3), PTGQ,
     .        PTN1, PTN2, Q(3), QNORM, S1, S2, S3, S12, SUM
C
C Local parameters:
C
C A =        Angular separation between P and Q
C B1,B2,B3 = Barycentric coordinates of the central projec-
C              tion of P onto the underlying planar triangle,
C              or (B1 and B2) projection of Q onto the
C              underlying line segment N1-N2 when P is
C              exterior.  Unnormalized coordinates are
C              computed by TRFIND when P is in a triangle.
C DUM =      Dummy parameter for ARCINT
C FQ,GQ =    Interpolated value and gradient at Q
C GQN =      Negative of the component of GQ in the direction
C              Q->P
C G1,G2,G3 = Gradients at I1, I2, and I3, or (G1 and G2) at
C              N1 and N2
C I =        DO-loop index
C IERR =     Error flag for calls to GRADL
C I1,I2,I3 = Vertex indexes returned by TRFIND
C LP =       LIST pointer
C N1,N2 =    Indexes of the endpoints of a boundary arc when
C              P is exterior (not contained in a triangle)
C NN =       Local copy of N
C P =        Cartesian coordinates of P
C P1,P2,P3 = Cartesian coordinates of the vertices I1, I2,
C              and I3, or (P1 and P2) coordinates of N1 and
C              N2 if P is exterior
C PTGQ =     Scalar product (P,GQ) -- factor of the component
C              of GQ in the direction Q->P
C PTN1 =     Scalar product (P,N1) -- factor of B1 and B2
C PTN2 =     Scalar product (P,N2) -- factor of B1 and B2
C Q =        Closest boundary point to P when P is exterior
C QNORM =    Factor used to normalize Q
C S1,S2,S3 = Tension factors associated with the triangle
C              sides opposite I1, I2, and I3, or (S1) the
C              boundary arc N1-N2
C S12 =      Scalar product (N1,N2) -- factor of B1 and B2
C SUM =      Quantity used to normalize the barycentric
C              coordinates
C
      NN = N
      IF (NN .LT. 3  .OR.  (IFLGG .LE. 0  .AND.  NN .LT. 7)
     .    .OR.  IST .LT. 1  .OR.  IST .GT. NN) GO TO 11
C
C Transform (PLAT,PLON) to Cartesian coordinates.
C
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
C
C Locate P with respect to the triangulation.
C
      CALL TRFIND (IST,P,NN,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,
     .             I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 12
      IST = I1
      IF (I3 .NE. 0) THEN
C
C P is contained in the triangle (I1,I2,I3).  Store the
C   vertex coordinates, gradients, and tension factors in
C   local variables.
C
        P1(1) = X(I1)
        P1(2) = Y(I1)
        P1(3) = Z(I1)
        P2(1) = X(I2)
        P2(2) = Y(I2)
        P2(3) = Z(I2)
        P3(1) = X(I3)
        P3(2) = Y(I3)
        P3(3) = Z(I3)
        IF (IFLGG .GT. 0) THEN
C
C   Gradients are user-provided.
C
          DO 1 I = 1,3
            G1(I) = GRAD(I,I1)
            G2(I) = GRAD(I,I2)
            G3(I) = GRAD(I,I3)
    1       CONTINUE
        ELSE
C
C   Compute gradient estimates at the vertices.
C
          CALL GRADL (NN,I1,X,Y,Z,F,LIST,LPTR,LEND, G1,IERR)
          IF (IERR .LT. 0) GO TO 12
          CALL GRADL (NN,I2,X,Y,Z,F,LIST,LPTR,LEND, G2,IERR)
          IF (IERR .LT. 0) GO TO 12
          CALL GRADL (NN,I3,X,Y,Z,F,LIST,LPTR,LEND, G3,IERR)
          IF (IERR .LT. 0) GO TO 12
        ENDIF
C
        IF (IFLGS .GT. 0) THEN
C
C   Variable tension:
C
          LP = LSTPTR(LEND(I2),I3,LIST,LPTR)
          S1 = SIGMA(LP)
          LP = LSTPTR(LEND(I3),I1,LIST,LPTR)
          S2 = SIGMA(LP)
          LP = LSTPTR(LEND(I1),I2,LIST,LPTR)
          S3 = SIGMA(LP)
        ELSE
C
C   Uniform tension:
C
          S1 = SIGMA(1)
          S2 = S1
          S3 = S1
        ENDIF
C
C Normalize the coordinates.
C
        SUM = B1 + B2 + B3
        B1 = B1/SUM
        B2 = B2/SUM
        B3 = B3/SUM
        FP = FVAL(B1,B2,B3,P1,P2,P3,F(I1),F(I2),F(I3),G1,
     .            G2,G3,S1,S2,S3)
        IER = 0
        RETURN
      ENDIF
C
C P is exterior to the triangulation, and I1 and I2 are
C   boundary nodes which are visible from P.  Extrapolate to
C   P by linear (with respect to arc-length) interpolation
C   of the value and directional derivative (gradient comp-
C   onent in the direction Q->P) of the interpolatory
C   surface at Q where Q is the closest boundary point to P.
C
C Determine Q by traversing the boundary starting from I1.
C
      N1 = I1
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
      IF (I1 .NE. I2) GO TO 3
C
C All boundary nodes are visible from P.  Find a boundary
C   arc N1->N2 such that P Left (N2 X N1)->N1.
C
C Counterclockwise boundary traversal:
C   Set N2 to the first neighbor of N1.
C
    2 LP = LEND(N1)
        LP = LPTR(LP)
        N2 = LIST(LP)
C
C Compute inner products (N1,N2) and (P,N2), and compute
C   B2 = Det(P,N1,N2 X N1).
C
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN2 = P(1)*X(N2) + P(2)*Y(N2) + P(3)*Z(N2)
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 3
C
C P Right (N2 X N1)->N1:  iterate.
C
        N1 = N2
        I1 = N1
        PTN1 = PTN2
        GO TO 2
C
C P Left (N2 X N1)->N1 where N2 is the first neighbor of N1.
C   Clockwise boundary traversal:
C
    3 N2 = N1
        PTN2 = PTN1
C
C Set N1 to the last neighbor of N2 and test for
C   termination.
C
        LP = LEND(N2)
        N1 = -LIST(LP)
        IF (N1 .EQ. I1) GO TO 13
C
C Compute inner products (N1,N2) and (P,N1).
C
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
C
C Compute B2 = Det(P,N1,N2 X N1) = Det(Q,N1,N2 X N1)*(P,Q).
C
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 3
C
C Compute B1 = Det(P,N2 X N1,N2) = Det(Q,N2 X N1,N2)*(P,Q).
C
      B1 = PTN1 - S12*PTN2
      IF (B1 .LE. 0.) THEN
C
C Q = N2.  Store value, coordinates, and gradient at Q.
C
        FQ = F(N2)
        Q(1) = X(N2)
        Q(2) = Y(N2)
        Q(3) = Z(N2)
        IF (IFLGG .GT. 0) THEN
          DO 4 I = 1,3
            GQ(I) = GRAD(I,N2)
    4       CONTINUE
        ELSE
          CALL GRADL (NN,N2,X,Y,Z,F,LIST,LPTR,LEND, GQ,IERR)
          IF (IERR .LT. 0) GO TO 12
        ENDIF
C
C Extrapolate to P:  FP = FQ + A*(GQ,Q X (PXQ)/SIN(A)),
C   where A is the angular separation between Q and P,
C   and Sin(A) is the magnitude of P X Q.
C
        A = ARCLEN(Q,P)
        PTGQ = P(1)*GQ(1) + P(2)*GQ(2) + P(3)*GQ(3)
        FP = FQ
        IF (A .NE. 0.) FP = FP + PTGQ*A/SIN(A)
        IER = 1
        RETURN
      ENDIF
C
C P Strictly Left (N2 X N1)->N2 and P Strictly Left
C   N1->(N2 X N1).  Thus Q lies on the interior of N1->N2.
C   Store coordinates of N1 and N2 in local variables.
C
      P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
C
C Compute the central projection of Q onto P2-P1 and
C   normalize to obtain Q.
C
      QNORM = 0.
      DO 5 I = 1,3
        Q(I) = B1*P1(I) + B2*P2(I)
        QNORM = QNORM + Q(I)*Q(I)
    5   CONTINUE
      QNORM = SQRT(QNORM)
      DO 6 I = 1,3
        Q(I) = Q(I)/QNORM
    6   CONTINUE
C
C Store gradients at N1 and N2 and tension factor S1.
C
      IF (IFLGG .GT. 0) THEN
        DO 7 I = 1,3
          G1(I) = GRAD(I,N1)
          G2(I) = GRAD(I,N2)
    7     CONTINUE
      ELSE
        CALL GRADL (NN,N1,X,Y,Z,F,LIST,LPTR,LEND, G1,IERR)
        IF (IERR .LT. 0) GO TO 12
        CALL GRADL (NN,N2,X,Y,Z,F,LIST,LPTR,LEND, G2,IERR)
        IF (IERR .LT. 0) GO TO 12
      ENDIF
C
      IF (IFLGS .LE. 0) S1 = SIGMA(1)
      IF (IFLGS .GE. 1) S1 = SIGMA(LP)
C
C Compute an interpolated value and normal gradient
C   component at Q.
C
      CALL ARCINT (Q,P1,P2,F(N1),F(N2),G1,G2,S1, FQ,DUM,GQN)
C
C Extrapolate to P:  the normal gradient component GQN is
C   the negative of the component in the direction Q->P.
C
      FP = FQ - GQN*ARCLEN(Q,P)
      IER = 1
      RETURN
C
C N or IST is outside its valid range.
C
   11 IER = -1
      RETURN
C
C Collinear nodes encountered.
C
   12 IER = -2
      RETURN
C
C The distance between P and the closest boundary point
C   is at least 90 degrees.
C
   13 IER = -3
      RETURN
      END
      SUBROUTINE ROTATE (N,C,S, X,Y )
      INTEGER N
      REAL    C, S, X(N), Y(N)
C
C***********************************************************
C
C                                              From SSRFPACK
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
      SUBROUTINE SETUP (XI,YI,WI,WK,S1,S2,WT, ROW)
      REAL XI, YI, WI, WK, S1, S2, WT, ROW(6)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This subroutine sets up the I-th row of an augmented
C regression matrix for a weighted least squares fit of a
C quadratic function Q(X,Y) to a set of data values Wi,
C where Q(0,0) = Wk.  The first 3 columns (quadratic terms)
C are scaled by 1/S2 and the fourth and fifth columns (lin-
C ear terms) are scaled by 1/S1.
C
C On input:
C
C       XI,YI = Coordinates of node I.
C
C       WI = Data value at node I.
C
C       WK = Data value interpolated by Q at the origin.
C
C       S1,S2 = Inverse scale factors.
C
C       WT = Weight factor corresponding to the I-th
C            equation.
C
C       ROW = Array of length 6.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       ROW = Array containing a row of the augmented re-
C             gression matrix.
C
C Modules required by SETUP:  None
C
C***********************************************************
C
      REAL W1, W2
C
C Local parameters:
C
C W1 = Weighted scale factor for the linear terms
C W2 = Weighted scale factor for the quadratic terms
C
      W1 = WT/S1
      W2 = WT/S2
      ROW(1) = XI*XI*W2
      ROW(2) = XI*YI*W2
      ROW(3) = YI*YI*W2
      ROW(4) = XI*W1
      ROW(5) = YI*W1
      ROW(6) = (WI-WK)*WT
      RETURN
      END
      SUBROUTINE SGPRNT (N,LUNIT,LIST,LPTR,LEND,SIGMA)
      INTEGER N, LUNIT, LIST(*), LPTR(*), LEND(N)
      REAL    SIGMA(*)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/21/98
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with an array of tension factors associated
C with the triangulation arcs, this subroutine prints the
C list of arcs (with tension factors) ordered by endpoint
C nodal indexes.  An arc is identified with its smaller
C endpoint index:  N1-N2, where N1 < N2.
C
C   This routine is identical to the similarly named routine
C in SRFPACK.
C
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
C                        gulation.  Refer to STRIPACK
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
C STRIPACK module required by SGPRNT:  LSTPTR
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
      DATA NMAX/9999/,  NLMAX/58/
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
      IF (NB .NE. 0) THEN
        NAT = 3*NM1 - NB
      ELSE
        NAT = 3*N - 6
      ENDIF
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
  100 FORMAT ('1',14X,'TENSION FACTORS,  N =',I5,
     .        ' NODES'//1X,18X,'N1',5X,'N2',8X,'TENSION'//)
  110 FORMAT (1X,16X,I4,3X,I4,5X,F12.8)
  120 FORMAT (1X,16X,I4,3X,I4,5X,F12.8,3X,F12.8,' *')
  130 FORMAT ('1')
  140 FORMAT (//1X,10X,'NA =',I5,' ARCS')
C
C Error messages:
C
  200 FORMAT (//1X,10X,'*',I5,' ERRORS IN SIGMA')
  210 FORMAT (/1X,10X,'*** ERROR IN TRIANGULATION -- ',
     .        '3N-NB-3 = ',I5,' ***')
  220 FORMAT (1X,10X,'*** N IS OUT OF RANGE -- NMAX = ',
     .        I4,' ***')
      END
      REAL FUNCTION SIG0 (N1,N2,N,X,Y,Z,H,LIST,LPTR,LEND,
     .                    GRAD,IFLGB,HBND,TOL,
     .                    IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGB,
     .        IFLGS, IER
      REAL    X(N), Y(N), Z(N), H(N), GRAD(3,N), HBND, TOL,
     .        SIGMA(*)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values H and gradients GRAD at the
C nodes, this function determines the smallest tension fac-
C tor SIG0 such that the Hermite interpolatory tension
C spline H(A), defined by SIG0 and the endpoint values and
C directional derivatives associated with an arc N1-N2, is
C bounded (either above or below) by HBND for all A in
C (A1,A2), where (A1,A2) denotes an interval corresponding
C to the arc and A is the arc-length.
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
C       X,Y,Z = Arrays of length N containing coordinates of
C               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
C           for I = 1 to N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       GRAD = Array dimensioned 3 by N whose columns con-
C              tain gradients at the nodes.  GRAD( ,J) must
C              be orthogonal to node J:  GRAD(1,J)*X(J) +
C              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0.  Refer
C              to Subroutines GRADG, GRADL, and SMSURF.
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
C       SIGMA = Dummy parameter (IFLGS .LE. 0) or array con-
C               taining tension factors associated with arcs
C               in one-to-one correspondence with LIST
C               entries (IFLGS .GE. 1).  Refer to Subroutine
C               GETSIG.
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
C                     = H(A1), and the directional deriva-
C                     tive of H at A1 is negative).
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
C STRIPACK module required by SIG0:  STORE
C
C SSRFPACK modules required by SIG0:  ARCLEN, SNHCSH
C
C Intrinsic functions called by SIG0:  ABS, EXP, LOG, MAX,
C                                        MIN, REAL, SIGN,
C                                        SQRT
C
C***********************************************************
C
      REAL    ARCLEN, STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    A, A0, AA, AL, B, B0, BND, C, C1, C2, COSHM,
     .        COSHMM, D, D0, D1PD2, D2, DMAX, DSIG, E, EMS,
     .        F, F0, FMAX, FNEG, FTOL, H1, H2, P1(3), P2(3),
     .        R, RF, RSIG, RTOL, S, S1, S2, SBIG, SCM, SIG,
     .        SINHM, SNEG, SSINH, SSM, STOL, T, T0, T1, T2,
     .        TM, UN(3), UNORM
C
      DATA SBIG/85./,  LUN/-1/
      RF = REAL(IFLGB)
      BND = HBND
C
C Print a heading.
C
      IF (LUN .GE. 0  .AND.  RF .LT. 0.) WRITE (LUN,100) N1,
     .                                   N2, BND
      IF (LUN .GE. 0  .AND.  RF .GT. 0.) WRITE (LUN,110) N1,
     .                                   N2, BND
  100 FORMAT (//1X,'SIG0 -- N1 =',I4,', N2 =',I4,
     .        ', LOWER BOUND = ',E15.8)
  110 FORMAT (//1X,'SIG0 -- N1 =',I4,', N2 =',I4,
     .        ', UPPER BOUND = ',E15.8)
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
C Store nodal coordinates P1 and P2, compute arc-length AL
C   and unit normal UN = (P1 X P2)/UNORM, and test for
C   coincident nodes.
C
    4 P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
      AL = ARCLEN(P1,P2)
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.  .OR.  AL .EQ. 0.) GO TO 11
C
C Store endpoint data values and test for valid constraint.
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
      S1 = AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +
     .         GRAD(3,N1)*P2(3))/UNORM
      S2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +
     .          GRAD(3,N2)*P1(3))/UNORM
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
C     of the Hermite cubic interpolant H0(A) = H2 - (S2*R +
C     B0*R**2 + (A0/3)*R**3), where R(A) = (A2-A)/AL.
C
      S = H2 - H1
      T0 = 3.*S - S1 - S2
      A0 = 3.*(S-T0)
      B0 = T0 - S2
      D0 = T0*T0 - S1*S2
C
C   H0 has local extrema in (A1,A2) iff S1*S2 < 0 or
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
C   HP(R) = (S2 - (C1*sinh(SIG*R) - C2*coshm(SIG*R))/E)/DT
C     = 0 for ESR = (-B +/- sqrt(D))/A = C/(-B -/+ sqrt(D))
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
C   H(R) = H2 - (B*SIG*R + C1 + RF*sqrt(D))/(SIG*E).
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
      REAL FUNCTION SIG1 (N1,N2,N,X,Y,Z,H,LIST,LPTR,LEND,GRAD,
     .                    IFLGB,HPBND,TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGB,
     .        IFLGS, IER
      REAL    X(N), Y(N), Z(N), H(N), GRAD(3,N), HPBND, TOL,
     .        SIGMA(*)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values H and gradients GRAD at the
C nodes, this function determines the smallest tension fac-
C tor SIG1 such that the first derivative HP(A) of the
C Hermite interpolatory tension spline H(A), defined by SIG1
C and the endpoint values and directional derivatives asso-
C ciated with an arc N1-N2, is bounded (either above or
C below) by HPBND for all A in (A1,A2), where (A1,A2) de-
C notes an interval corresponding to the arc and A denotes
C arc-length.
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
C       X,Y,Z = Arrays of length N containing coordinates of
C               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
C           for I = 1 to N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       GRAD = Array dimensioned 3 by N whose columns con-
C              gradients at the nodes.  GRAD( ,J) must be
C              orthogonal to node J:  GRAD(1,J)*X(J) +
C              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0.  Refer
C              to Subroutines GRADG, GRADL, and SMSURF.
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
C       SIGMA = Dummy parameter (IFLGS .LE. 0) or array
C               containing tension factors associated with
C               arcs in one-to-one correspondence with LIST
C               entries (IFLGS .GE. 1).  Refer to Subroutine
C               GETSIG.
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
C STRIPACK module required by SIG1:  STORE
C
C SSRFPACK modules required by SIG1:  ARCLEN, SNHCSH
C
C Intrinsic functions called by SIG1:   ABS, EXP, MAX, MIN,
C                                         REAL, SIGN, SQRT
C
C***********************************************************
C
      REAL    ARCLEN, STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    A, A0, AL, B0, BND, C0, C1, C2, COSHM, COSHMM,
     .        D0, D1, D1PD2, D2, DMAX, DSIG, E, EMS, EMS2,
     .        F, F0, FMAX, FNEG, FTOL, P1(3), P2(3), RF,
     .        RTOL, S, S1, S2, SBIG, SIG, SINH, SINHM, STOL,
     .        T0, T1, T2, TM, UN(3), UNORM
C
      DATA SBIG/85./,  LUN/-1/
      RF = REAL(IFLGB)
      BND = HPBND
C
C Print a heading.
C
      IF (LUN .GE. 0  .AND.  RF .LT. 0.) WRITE (LUN,100) N1,
     .                                   N2, BND
      IF (LUN .GE. 0  .AND.  RF .GT. 0.) WRITE (LUN,110) N1,
     .                                   N2, BND
  100 FORMAT (//1X,'SIG1 -- N1 =',I4,', N2 =',I4,
     .        ', LOWER BOUND = ',E15.8)
  110 FORMAT (//1X,'SIG1 -- N1 =',I4,', N2 =',I4,
     .        ', UPPER BOUND = ',E15.8)
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
C Store nodal coordinates P1 and P2, compute arc-length AL
C   and unit normal UN = (P1 X P2)/UNORM, and test for
C   coincident nodes.
C
    4 P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
      AL = ARCLEN(P1,P2)
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.  .OR.  AL .EQ. 0.) GO TO 11
C
C Compute first difference S and scaled directional deriva-
C   tives S1,S2 at the endpoints (for the direction N1->N2).
C
      S = H(N2) - H(N1)
      S1 = AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +
     .         GRAD(3,N1)*P2(3))/UNORM
      S2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +
     .          GRAD(3,N2)*P1(3))/UNORM
C
C Test for a valid constraint.
C
      IER = -3
      IF ((RF .LT. 0.  .AND.  MIN(S1,S2,S) .LT. BND)  .OR.
     .    (RF .GT. 0.  .AND.  BND .LT. MAX(S1,S2,S)))
     .   GO TO 11
C
C Test for infinite tension required.
C
      IER = 1
      SIG = SBIG
      IF (S .EQ. BND  .AND.  (S1 .NE. S  .OR.  S2 .NE. S))
     .   GO TO 10
C
C Test for SIG = 0 sufficient.  The Hermite cubic interpo-
C   lant H0 has derivative HP0(T) = (S2 + 2*B0*R + A0*R**2)/
C   AL, where R = (T2-T)/AL.
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
      IF (B0*C0 .LE. 0.  .OR.  A0*RF .GT. 0.) GO TO 10
C
C   A0*RF < 0 and HP0(R) = -D0/(DT*A0) at R = -B0/A0.
C
      D0 = T0*T0 - S1*S2
      F0 = (BND + D0/(A0*AL))*RF
      IF (F0 .GE. 0.) GO TO 10
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
C   (AL*(SIG-2.)))*RF -- a value for which F(SIG) .GE. 0 and
C   F(SIG) = 0 for SIG sufficiently large that 2*SIG is in-
C   significant relative to exp(SIG).
C
      FMAX = (BND-S/AL)*RF
      SIG = 2. - A0/(3.*(AL*BND-S))
      IF (LUN .GE. 0) WRITE (LUN,120) F0, FMAX, SIG
  120 FORMAT (1X,9X,'F0 = ',E15.8,', FMAX = ',E15.8/
     .        1X,8X,'SIG = ',E15.8/)
      IF (STORE(SIG*EXP(-SIG)+.5) .EQ. .5) GO TO 10
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
C   HP(R) = (B+SIGN(A)*SQRT(A*C))/(AL*E) at the critical
C     value of R, where A = C2-C1, B = E*S2-C2, and C = C2 +
C     C1.  Note that RF*A < 0.
C
          F = (BND - (E*S2-C2 - RF*SQRT(T1))/(AL*E))*RF
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
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 10
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
   10 SIG1 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
C
C Error termination.
C
   11 SIG1 = -1.
      RETURN
      END
      REAL FUNCTION SIG2 (N1,N2,N,X,Y,Z,H,LIST,LPTR,LEND,
     .                    GRAD,TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGS,
     .        IER
      REAL    X(N), Y(N), Z(N), H(N), GRAD(3,N), TOL,
     .        SIGMA(*)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/25/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values H and gradients GRAD at the
C nodes, this function determines the smallest tension fac-
C tor SIG2 such that the Hermite interpolatory tension
C spline H(A), defined by SIG2 and the endpoint values and
C directional derivatives associated with an arc N1-N2,
C preserves convexity (or concavity) of the data:
C
C   HP1 .LE. S .LE. HP2 implies HPP(A) .GE. 0, and
C   HP1 .GE. S .GE. HP2 implies HPP(A) .LE. 0
C
C for all A in the open interval (A1,A2) corresponding to
C the arc, where HP1 and HP2 are the derivative values of H
C at the endpoints, S is the slope of the linear interpolant
C of the endpoint data values, HPP denotes the second deriv-
C ative of H, and A is arc-length.  Note, however, that
C infinite tension is required if HP1 = S or HP2 = S (unless
C HP1 = HP2 = S).
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
C       X,Y,Z = Arrays of length N containing coordinates of
C               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
C           for I = 1 to N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       GRAD = Array dimensioned 3 by N whose columns con-
C              gradients at the nodes.  GRAD( ,J) must be
C              orthogonal to node J:  GRAD(1,J)*X(J) +
C              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0.  Refer
C              to Subroutines GRADG, GRADL, and SMSURF.
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
C       SIGMA = Dummy parameter (IFLGS .LE. 0) or array
C               containing tension factors associated with
C               arcs in one-to-one correspondence with LIST
C               entries (IFLGS .GE. 1).  Refer to Subroutine
C               GETSIG.
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
C STRIPACK module required by SIG2:  STORE
C
C SSRFPACK modules required by SIG2:  ARCLEN, SNHCSH
C
C Intrinsic functions called by SIG2:  ABS, EXP, MAX, MIN,
C                                        SQRT
C
C***********************************************************
C
      REAL    ARCLEN, STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    AL, COSHM, D1, D1D2, D2, DSIG, DUMMY, EMS, F,
     .        FP, FTOL, P1(3), P2(3), RTOL, S, SBIG, SIG,
     .        SINHM, SSM, T, T1, TP1, UN(3), UNORM
C
      DATA SBIG/85./,  LUN/-1/
C
C Print a heading.
C
      IF (LUN .GE. 0) WRITE (LUN,100) N1, N2
  100 FORMAT (//1X,'SIG2 -- N1 =',I4,', N2 =',I4)
C
C Test for errors and set local parameters.
C
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.
     .    MAX(N1,N2,3) .GT. N) GO TO 11
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
C Store nodal coordinates P1 and P2, compute arc-length AL
C   and unit normal UN = (P1 X P2)/UNORM, and test for
C   coincident nodes.
C
    4 P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
      AL = ARCLEN(P1,P2)
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.  .OR.  AL .EQ. 0.) GO TO 11
C
C Compute first and second differences and test for infinite
C   tension required.
C
      S = H(N2) - H(N1)
      D1 = S - AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +
     .             GRAD(3,N1)*P2(3))/UNORM
      D2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +
     .          GRAD(3,N2)*P1(3))/UNORM - S
      D1D2 = D1*D2
      IER = 1
      SIG = SBIG
      IF (D1D2 .EQ. 0.  .AND.  D1 .NE. D2) GO TO 10
C
C Test for a valid constraint.
C
      IER = 2
      SIG = 0.
      IF (D1D2 .LT. 0.) GO TO 10
C
C Test for SIG = 0 sufficient.
C
      IER = 0
      IF (D1D2 .EQ. 0.) GO TO 10
      T = MAX(D1/D2,D2/D1)
      IF (T .LE. 2.) GO TO 10
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
      IF (FP .LE. 0.) GO TO 10
      DSIG = -F/FP
      IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 10
C
C   Bottom of loop:  update SIG.
C
      SIG = SIG + DSIG
      GO TO 6
C
C No errors encountered.
C
   10 SIG2 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
C
C Error termination.
C
   11 SIG2 = -1.
      RETURN
      END





      SUBROUTINE SMSGS (N,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,
     .                  SIGMA,W,P, NIT,DFMAX,F,GRAD, IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, NIT, IER
      REAL    X(N), Y(N), Z(N), U(N), SIGMA(*), W(N), P,
     .        DFMAX, F(N), GRAD(3,N)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/25/96
C
C   This subroutine solves the symmetric positive definite
C linear system associated with minimizing the quadratic
C functional Q(F,FX,FY,FZ) described in Subroutine SMSURF.
C Since the gradient at node K lies in the plane tangent to
C the sphere surface at K, it is effectively defined by only
C two components -- its X and Y components in the coordinate
C system obtained by rotating K to the north pole.  Thus,
C the minimization problem corresponds to an order-3N system
C which is solved by the block Gauss-Seidel method with 3 by
C 3 blocks.
C
C On input:
C
C       N,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,SIGMA,W = Parameters
C           as described in Subroutine SMSURF.
C
C       P = Positive smoothing parameter defining Q.
C
C The above parameters are not altered by this routine.
C
C       NIT = Maximum number of iterations to be used.  This
C             maximum will likely be achieved if DFMAX is
C             smaller than the machine precision.  NIT .GE.
C             0.
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
C       GRAD = 3 by N array containing initial estimates of
C              the last 3N solution components (the gradi-
C              ent with FX, FY, and FZ in rows 1, 2, and 3,
C              respectively).
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
C       GRAD = Last 3N solution components -- gradients at
C              the nodes.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     convergence criterion was achieved.
C             IER = 1 if no errors were encountered but con-
C                     vergence was not achieved within NIT
C                     iterations.
C             IER = -1 if N, P, NIT, or DFMAX is outside its
C                      valid range on input.  F and GRAD are
C                      not altered in this case.
C             IER = -2 if all nodes are collinear or the
C                      triangulation is invalid.
C             IER = -3 if duplicate nodes were encountered.
C
C SSRFPACK modules required by SMSGS:  APLYRT, CONSTR,
C                                        GRCOEF, SNHCSH
C
C Intrinsic functions called by SMSGS:  ABS, ATAN, MAX, SQRT
C
C***********************************************************
C
      INTEGER IFL, ITER, ITMAX, J, K, LPJ, LPL, NN
      REAL    ALFA, ALFSQ, C11, C12, C13, C22, C23, C33,
     .        CC22, CC23, CC33, CX, CY, DEN1, DEN2, DET, DF,
     .        DFMX, DGK(3), DGX, DGY, FK, G1, G2, G3, GJK,
     .        GKJ, PP, R1, R2, R3, RR2, RR3, SIG, SINAL, SX,
     .        SY, T, T1, T2, T3, T4, T5, T6, TOL, XJ, XK,
     .        XS, YJ, YK, YS, ZJ, ZK
C
      NN = N
      IFL = IFLGS
      PP = P
      ITMAX = NIT
      TOL = DFMAX
C
C Test for errors in input and initialize iteration count,
C   tension factor, and output value of DFMAX.
C
      IF (NN .LT. 3  .OR.  PP .LE. 0.  .OR.  ITMAX .LT. 0
     .    .OR.  TOL .LT. 0.) GO TO 5
      ITER = 0
      SIG = SIGMA(1)
      DFMX = 0.
C
C Top of iteration loop.
C
    1 IF (ITER .EQ. ITMAX) GO TO 4
      DFMX = 0.
C
C   Loop on nodes.
C
      DO 3 K = 1,NN
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        FK = F(K)
        G1 = GRAD(1,K)
        G2 = GRAD(2,K)
        G3 = GRAD(3,K)
C
C   Construct the rotation mapping node K to the north pole.
C
        CALL CONSTR (XK,YK,ZK, CX,SX,CY,SY)
C
C   Initialize components of the order-3 system for the
C     change (DF,DGX,DGY) in the K-th solution components.
C
        C11 = PP*W(K)
        C12 = 0.
        C13 = 0.
        C22 = 0.
        C23 = 0.
        C33 = 0.
        R1 = C11*(U(K)-FK)
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
C   Compute the coordinates of J in the rotated system.
C
          T = SX*Y(J) + CX*Z(J)
          YJ = CX*Y(J) - SX*Z(J)
          ZJ = SY*X(J) + CY*T
          XJ = CY*X(J) - SY*T
C
C   Compute arc-length ALFA between K and J, ALFSQ = ALFA*
C     ALFA, SINAL = SIN(ALFA), DEN1 = ALFA*SIN(ALFA)**2, and
C     DEN2 = ALFSQ*SINAL.
C
          ALFA = 2.*ATAN(SQRT((1.-ZJ)/(1.+ZJ)))
          ALFSQ = ALFA*ALFA
          XS = XJ*XJ
          YS = YJ*YJ
          SINAL = SQRT(XS+YS)
          DEN1 = ALFA*(XS+YS)
          DEN2 = ALFSQ*SINAL
C
C   Test for coincident nodes and compute functions of SIG:
C     T1 = SIG*SIG*COSHM/E, T2 = SIG*SINHM/E, and T3 = SIG*
C     (SIG*COSHM-SINHM)/E for E = SIG*SINH - 2*COSHM.
C
          IF (DEN1 .EQ. 0.) GO TO 7
          IF (IFL .GE. 1) SIG = SIGMA(LPJ)
          CALL GRCOEF (SIG, T3,T2)
          T1 = T2 + T3
C
C   Update system components for node J.
C
          T4 = 2.*T1/(ALFA*ALFSQ)
          T5 = T1/DEN2
          T6 = T3/DEN1
          C11 = C11 + T4
          C12 = C12 + T5*XJ
          C13 = C13 + T5*YJ
          C22 = C22 + T6*XS
          C23 = C23 + T6*XJ*YJ
          C33 = C33 + T6*YS
          GKJ = G1*X(J) + G2*Y(J) + G3*Z(J)
          GJK = GRAD(1,J)*XK + GRAD(2,J)*YK + GRAD(3,J)*ZK
          R1 = R1 + T4*(F(J)-FK) + T5*(GJK-GKJ)
          T = T5*(F(J)-FK) - T6*GKJ + T2*GJK/DEN1
          R2 = R2 + T*XJ
          R3 = R3 + T*YJ
C
C   Bottom of loop on neighbors.
C
          IF (LPJ .NE. LPL) GO TO 2
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
     .      C11 .EQ. 0.) GO TO 6
        DGY = (CC22*RR3 - CC23*RR2)/DET
        DGX = (RR2 - CC23*DGY)/CC22
        DF = (R1 - C12*DGX - C13*DGY)/C11
C
C   Rotate (DGX,DGY,0) back to the original coordinate
C     system, and update GRAD( ,K), F(K), and DFMX.
C
        CALL APLYRT (DGX,DGY,CX,SX,CY,SY, DGK)
        GRAD(1,K) = G1 + DGK(1)
        GRAD(2,K) = G2 + DGK(2)
        GRAD(3,K) = G3 + DGK(3)
        F(K) = FK + DF
        DFMX = MAX(DFMX,ABS(DF)/(1.+ABS(FK)))
    3   CONTINUE
C
C   Increment ITER and test for convergence.
C
      ITER = ITER + 1
      IF (DFMX .GT. TOL) GO TO 1
C
C The method converged.
C
      NIT = ITER
      DFMAX = DFMX
      IER = 0
      RETURN
C
C The method failed to converge within NIT iterations.
C
    4 DFMAX = DFMX
      IER = 1
      RETURN
C
C Invalid input parameter.
C
    5 NIT = 0
      DFMAX = 0.
      IER = -1
      RETURN
C
C Node K and its neighbors are collinear.
C
    6 NIT = 0
      DFMAX = DFMX
      IER = -2
      RETURN
C
C Nodes J and K coincide.
C
    7 NIT = 0
      DFMAX = DFMX
      IER = -3
      RETURN
      END
      SUBROUTINE SMSURF (N,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,
     .                   SIGMA,W,SM,SMTOL,GSTOL,LPRNT, F,
     .                   GRAD,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, LPRNT,
     .        IER
      REAL    X(N), Y(N), Z(N), U(N), SIGMA(*), W(N), SM,
     .        SMTOL, GSTOL, F(N), GRAD(3,N)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/21/98
C
C   Given a triangulation of N nodes on the unit sphere with
C data values U at the nodes and tension factors SIGMA
C associated with the arcs, this routine determines a set of
C nodal function values F and gradients GRAD = (FX,FY,FZ)
C such that a quadratic functional Q1(F,GRAD) is minimized
C subject to the constraint Q2(F) .LE. SM for Q2(F) =
C (U-F)**T*W*(U-F), where W is a diagonal matrix of positive
C weights.  The functional Q1 is an approximation to the
C linearized curvature over the triangulation of a C-1 fun-
C ction F(V), V a unit vector, which interpolates the nodal
C values and gradients.  Subroutines INTRC1 and UNIF may be
C called to evaluate F at arbitrary points.
C
C   The smoothing procedure is an extension of the method
C for cubic spline smoothing due to C. Reinsch -- Numer.
C Math., 10 (1967) and 16 (1971).  Refer to Function FVAL
C for a further description of the interpolant F.  Letting
C D1F(T) and D2F(T) denote first and second derivatives of F
C with respect to a parameter T varying along a triangula-
C tion arc, Q1 is the sum of integrals over the arcs of
C D2F(T)**2 + ((SIGMA/L)*(D1F(T)-S))**2 where L denotes arc-
C length, SIGMA is the appropriate tension factor, and S is
C the slope of the linear function of T which interpolates
C the values of F at the endpoints of the arc.  Introducing
C a smoothing parameter P, and assuming the constraint is
C active, the problem is equivalent to minimizing Q(P,F,
C GRAD) = Q1(F,GRAD) + P*(Q2(F)-SM).  The secant method is
C used to find a zero of G(P) = 1/SQRT(Q2) - 1/SQRT(SM)
C where F(P) satisfies the order-3N symmetric positive def-
C inite linear system obtained by setting the gradient of Q
C (treated as a function of F and GRAD with GRAD tangent to
C the sphere surface) to zero.  The linear system is solved
C by the block Gauss-Seidel method (refer to SMSGS).
C
C   Note that the method can also be used to select grad-
C ients for the interpolation problem (F = U, SM = 0, and P
C infinite).  This is achieved by a call to Subroutine
C GRADG.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing Cartesian
C               coordinates of the nodes.
C
C       U = Array of length N containing data values at the
C           nodes.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
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
C               programs GETSIG, SIG0, SIG1, and SIG2.
C
C       W = Array of length N containing positive weights
C           associated with the data values.  The recommend-
C           ed value of W(I) is 1/DU**2 where DU is the
C           standard deviation associated with U(I).  DU**2
C           is the expected value of the squared error in
C           the measurement of U(I).  (The mean error is
C           assumed to be zero.)
C
C       SM = Positive parameter specifying an upper bound on
C            Q2(F).  Note that F is constant (and Q2(F)
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
C       LPRNT = Logical unit on which diagnostic messages
C               are printed, or negative integer specifying
C               no diagnostics.  For each secant iteration,
C               the following values are printed:  P, G(P),
C               NIT, DFMAX, and DP, where NIT denotes the
C               number of Gauss-Seidel iterations used in
C               the computation of G, DFMAX denotes the max-
C               imum relative change in a solution component
C               in the last Gauss-Seidel iteration, and DP
C               is the change in P computed by linear inter-
C               polation between the current point (P,G) and
C               a previous point.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       F = Array of length N containing nodal function val-
C           ues unless IER < 0.
C
C       GRAD = 3 by N array whose columns contain gradients
C              of F at the nodes unless IER < 0.
C
C       IER = Error indicator and information flag:
C             IER = 0 if no errors were encountered and the
C                     constraint is active -- Q2(F) is ap-
C                     proximately equal to SM.
C             IER = 1 if no errors were encountered but the
C                     constraint is not active -- F and GRAD
C                     are the values and gradients of a con-
C                     stant function which minimizes Q2(F),
C                     and Q1 = 0.
C             IER = 2 if the constraint could not be satis-
C                     fied to within SMTOL due to
C                     ill-conditioned linear systems.
C             IER = -1 if N, W, SM, SMTOL, or GSTOL is out-
C                      side its valid range on input.
C             IER = -2 if all nodes are collinear or the
C                      triangulation is invalid.
C             IER = -3 if duplicate nodes were encountered.
C
C SSRFPACK modules required by SMSURF:  APLYRT, CONSTR,
C                                         GRCOEF, SMSGS,
C                                         SNHCSH
C
C Intrinsic functions called by SMSURF:  ABS, SQRT
C
C***********************************************************
C
      INTEGER I, IERR, ITER, ITMAX, LUN, NIT, NITMAX, NN
      REAL    C, DFMAX, DMAX, DP, G, G0, GNEG, P, Q2, Q2MAX,
     .        Q2MIN, S, SUMW, TOL, WI
C
C Local parameters:
C
C ITMAX = Maximum number of secant iterations.
C LUN = Local copy of LPRNT.
C NITMAX = Maximum number of Gauss-Seidel iterations for
C          each secant iteration.
C NN = Local copy of N.
C TOL = Local copy of GSTOL.
C
      DATA ITMAX/50/,  NITMAX/40/
C
      NN = N
      TOL = GSTOL
      LUN = LPRNT
      IF (LUN .GT. 99) LUN = -1
C
C Test for errors and initialize F to the weighted least
C   squares fit of a constant function to the data.
C
      IER = -1
      IF (NN .LT. 3  .OR.  SM .LE. 0.  .OR.  SMTOL .LE. 0.
     .    .OR.  SMTOL .GE. 1.  .OR.  TOL .LE. 0.) RETURN
      C = 0.
      SUMW = 0.
      DO 1 I = 1,NN
        WI = W(I)
        IF (WI .LE. 0.) RETURN
        C = C + WI*U(I)
        SUMW = SUMW + WI
    1   CONTINUE
      C = C/SUMW
C
C Compute nodal values and gradients, and accumulate Q2 =
C   (U-F)**T*W*(U-F).
C
      Q2 = 0.
      DO 2 I = 1,NN
        F(I) = C
        GRAD(1,I) = 0.
        GRAD(2,I) = 0.
        GRAD(3,I) = 0.
        Q2 = Q2 + W(I)*(U(I)-F(I))**2
    2   CONTINUE
C
C Compute bounds on Q2 defined by SMTOL, and test for the
C   constraint satisfied by the constant fit.
C
      Q2MIN = SM*(1.-SMTOL)
      Q2MAX = SM*(1.+SMTOL)
      IF (Q2 .LE. Q2MAX) THEN
C
C The constraint is satisfied by a constant function.
C
        IER = 1
        IF (LUN .GE. 0) WRITE (LUN,100)
  100   FORMAT (///1X,'SMSURF -- THE CONSTRAINT IS NOT ',
     .          'ACTIVE AND THE FITTING FCN IS CONSTANT.')
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
C G(P) is strictly increasing and concave, and G(0) .LT. 0.
C   Initialize parameters for the secant method.  The method
C   uses three points -- (P0,G0), (P,G), and (PNEG,GNEG)
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
    3 NIT = NITMAX
      DFMAX = TOL
      CALL SMSGS (NN,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,SIGMA,W,
     .            P, NIT,DFMAX,F,GRAD, IERR)
      IF (IERR .LT. 0) IER = IERR
C
C   IERR = -1 in SMSGS could be caused by P = 0 as a result
C     of inaccurate solutions to ill-conditioned systems.
C
      IF (IERR .EQ. -1) IER = 2
      IF (IERR .LT. 0) RETURN
      Q2 = 0.
      DO 4 I = 1,NN
        Q2 = Q2 + W(I)*(U(I)-F(I))**2
    4   CONTINUE
      G = 1./SQRT(Q2) - S
      ITER = ITER + 1
      IF (LUN .GE. 0) WRITE (LUN,120) ITER, P, G, NIT, DFMAX
  120 FORMAT (/1X,I2,' -- P = ',E15.8,', G = ',E15.8,
     .        ', NIT = ',I2,', DFMAX = ',E12.6)
C
C   Test for convergence.
C
      IF (Q2MIN .LE. Q2  .AND.  Q2 .LE. Q2MAX) RETURN
      IF (ITER .GE. ITMAX) THEN
        IER = 2
        RETURN
      ENDIF
      IF (DMAX .EQ. 0.  .AND.  G .LE. 0.) THEN
C
C   Increase P until G(P) > 0.
C
        P = 10.*P
        DP = P
        GO TO 3
      ENDIF
C
C   A bracketing interval [P0,P] has been found.
C
      IF (G0*G .LE. 0.) THEN
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
    5 DP = -G*DP/(G-G0)
      IF (LUN .GE. 0) WRITE (LUN,130) DP
  130 FORMAT (1X,5X,'DP = ',E15.8)
      IF (ABS(DP) .GT. ABS(DMAX)) THEN
C
C   G0*G .GT. 0 and the new estimate would be outside of the
C     bracketing interval of length ABS(DMAX).  Reset
C     (P0,G0) to (PNEG,GNEG).
C
        DP = DMAX
        G0 = GNEG
        GO TO 5
      ENDIF
C
C   Bottom of loop -- update P, DMAX, and G0.
C
      P = P + DP
      DMAX = DMAX + DP
      G0 = G
      GO TO 3
      END
      SUBROUTINE SNHCSH (X, SINHM,COSHM,COSHMM)
      REAL X, SINHM, COSHM, COSHMM
C
C***********************************************************
C
C                                              From SSRFPACK
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
      SUBROUTINE UNIF (N,X,Y,Z,F,LIST,LPTR,LEND,IFLGS,SIGMA,
     .                 NROW,NI,NJ,PLAT,PLON,IFLGG, GRAD, FF,
     .                 IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, NROW, NI,
     .        NJ, IFLGG, IER
      REAL    X(N), Y(N), Z(N), F(N), SIGMA(*), PLAT(NI),
     .        PLON(NJ), GRAD(3,N), FF(NROW,NJ)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/25/96
C
C   Given a Delaunay triangulation of a set of nodes on the
C unit sphere, along with data values and tension factors
C associated with the triangulation arcs, this routine
C interpolates the data values to a uniform grid for such
C applications as contouring.  The interpolant is once con-
C tinuously differentiable.  Extrapolation is performed at
C grid points exterior to the triangulation when the nodes
C do not cover the entire sphere.
C
C On input:
C
C       N = Number of nodes.  N .GE. 3 and N .GE. 7 if
C           IFLAG .NE. 1.
C
C       X,Y,Z = Arrays containing Cartesian coordinates of
C               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1
C               for I = 1 to N.
C
C       F = Array containing data values.  F(I) is associ-
C           ated with (X(I),Y(I),Z(I)).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
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
C               programs GETSIG, SIG0, SIG1, and SIG2.
C
C       NROW = Number of rows in the dimension statement of
C              FF.
C
C       NI,NJ = Number of rows and columns in the uniform
C               grid.  1 .LE. NI .LE. NROW and 1 .LE. NJ.
C
C       PLAT,PLON = Arrays of length NI and NJ, respective-
C                   ly, containing the latitudes and
C                   longitudes of the grid lines.
C
C       IFLGG = Option indicator:
C               IFLGG = 0 if gradient estimates at the ver-
C                         tices of a triangle are to be
C                         recomputed for each grid point in
C                         the triangle and not saved.
C               IFLGG = 1 if gradient estimates are input in
C                         GRAD.
C               IFLGG = 2 if gradient estimates are to be
C                         computed once for each node (by
C                         GRADL) and saved in GRAD.
C
C The above parameters are not altered by this routine.
C
C       GRAD = 3 by N array whose columns contain the X, Y,
C              and Z components (in that order) of the grad-
C              ients at the nodes if IFLGG = 1, array of
C              sufficient size if IFLGG = 2, or dummy para-
C              meter if IFLGG = 0.
C
C Gradient estimates may be computed by Subroutines GRADL or
C   GRADG if IFLGG = 1.
C
C       FF = NROW by NCOL array with NROW .GE. NI and NCOL
C            .GE. NJ.
C
C On output:
C
C       GRAD = Array containing estimated gradients as de-
C              fined above if IFLGG = 2 and IER .GE. 0.
C              GRAD is not altered if IFLGG .NE. 2.
C
C       FF = Interpolated values at the grid points if IER
C            .GE. 0.  FF(I,J) = F(PLAT(I),PLON(J)) for I =
C            1,...,NI and J = 1,...,NJ.
C
C       IER = Error indicator:
C             IER = K if no errors were encountered and K
C                     grid points required extrapolation for
C                     K .GE. 0.
C             IER = -1 if N, NI, NJ, or IFLGG is outside its
C                      valid range.
C             IER = -2 if the nodes are collinear.
C             IER = -3 if extrapolation failed due to the
C                      uniform grid extending too far beyond
C                      the triangulation boundary.
C
C STRIPACK modules required by UNIF:  GETNP, JRAND, LSTPTR,
C                                       STORE, TRFIND
C
C SSRFPACK modules required by UNIF:  APLYR, APLYRT, ARCINT,
C                                       ARCLEN, CONSTR,
C                                       FVAL, GIVENS, GRADL,
C                                       HVAL, INTRC1,
C                                       ROTATE, SETUP,
C                                       SNHCSH
C
C***********************************************************
C
      INTEGER I, J, IERR, IFL, IST, NEX, NN, NST, NX, NY
      DATA    NST/1/
C
C Local parameters:
C
C I,J =   DO-loop indexes
C IERR =  Error flag for calls to GRADL and INTRC1
C IFL =   Local copy of IFLGG
C IST =   Parameter for INTRC1
C NEX =   Number of grid points exterior to the triangula-
C           tion boundary (number of extrapolated values)
C NN =    Local copy of N
C NST =   Initial value for IST
C NX,NY = Local copies of NI and NJ
C
      NN = N
      NX = NI
      NY = NJ
      IFL = IFLGG
      IF (NX .LT. 1  .OR.  NX .GT. NROW  .OR.  NY .LT. 1
     .   .OR.  IFL .LT. 0  .OR.  IFL .GT. 2) GO TO 4
      IST = NST
      IF (IFL .EQ. 2) THEN
C
C Compute gradient estimates at the nodes.
C
        DO 1 I = 1,NN
          CALL GRADL (NN,I,X,Y,Z,F,LIST,LPTR,
     .                LEND, GRAD(1,I),IERR)
          IF (IERR .LT. 0) GO TO 5
    1     CONTINUE
        IFL = 1
      ENDIF
C
C Compute uniform grid points and interpolated values.
C
      NEX = 0
      DO 3 J = 1,NY
        DO 2 I = 1,NX
          CALL INTRC1 (NN,PLAT(I),PLON(J),X,Y,Z,F,LIST,LPTR,
     .                 LEND,IFLGS,SIGMA,IFL,
     .                 GRAD, IST, FF(I,J),IERR)
          IF (IERR .LT. 0) GO TO 5
          NEX = NEX + IERR
    2     CONTINUE
    3   CONTINUE
      IER = NEX
      RETURN
C
C NI, NJ, or IFLGG is outside its valid range.
C
    4 IER = -1
      RETURN
C
C Error in GRADL or INTRC1.
C
    5 IER = IERR
      RETURN
      END
      SUBROUTINE interp_linear(n,ns,olats,olons,x,y,z,datain,lst,
     .               lptr,lend,odata,edata,ierr)

      INTEGER n, ns, ierr
      INTEGER lst(6*(n-2)), lptr(6*(n-2)), lend(n)
      REAL olats(ns), olons(ns), odata(ns)
      REAL datain(n), x(n), y(n), z(n)
      INTEGER i, ierr1, ist, edata(ns)

      ist = 1
      ierr = 0

      DO i=1,ns
         CALL intrc0(n,olats(i),olons(i),x,y,z,datain,lst,lptr,lend,
     .               ist,odata(i),ierr1)

         edata(i) = ierr1

         IF (ierr1 .lt. 0) THEN
           !print *,n,'warning: ierr = ',ierr1,' in interp_n'
           !print *,olats(n), olons(n), npts
           !stop
           ierr = ierr + ierr1
         ENDIF
      ENDDO
      END
      SUBROUTINE interp_cubic(n,ns,olats,olons,x,y,z,datain,lst,
     .               lptr,lend,iflgs,sigma,iflgg,grad,odata,edata,ierr)

      INTEGER n, ns, ierr
      INTEGER lst(6*(n-2)), lptr(6*(n-2)), lend(n)
      REAL olats(ns), olons(ns), odata(ns)
      REAL datain(n), x(n), y(n), z(n), sigma(6*(n-2)), grad(3,n)
      INTEGER i, ierr1, ist, iflgs, iflgg, edata(ns)

      ist = 1
      ierr = 0

      DO i=1,ns
         CALL intrc1(n,olats(i),olons(i),x,y,z,datain,lst,lptr,lend,
     .                iflgs,sigma,iflgg,grad,ist,odata(i),ierr1)

         edata(i) = ierr1

         IF (ierr1 .lt. 0) THEN
           !print *,n,'warning: ierr = ',ierr1,' in interp_n'
           !print *,olats(n), olons(n), npts
           !stop
           ierr = ierr + ierr1
         ENDIF
      ENDDO
      END