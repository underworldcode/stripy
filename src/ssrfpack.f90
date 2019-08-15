      include 'stripack.f90'
      SUBROUTINE APLYR (X,Y,Z,CX,SX,CY,SY, XP,YP,ZP)
      REAL X, Y, Z, CX, SX, CY, SY, XP, YP, ZP
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   05/09/92
!
!   This subroutine applies the rotation R defined by Sub-
! routine CONSTR to the unit vector (X Y Z)**T, i,e. (X,Y,Z)
! is rotated to (XP,YP,ZP).  If (XP,YP,ZP) lies in the
! southern hemisphere (ZP < 0), (XP,YP) are set to the
! coordinates of the nearest point of the equator, ZP re-
! maining unchanged.
!
! On input:
!
!       X,Y,Z = Coordinates of a point on the unit sphere.
!
!       CX,SX,CY,SY = Elements of the rotation defined by
!                     Subroutine CONSTR.
!
! Input parameters are not altered except as noted below.
!
! On output:
!
!       XP,YP,ZP = Coordinates of the rotated point on the
!                  sphere unless ZP < 0, in which case
!                  (XP,YP,0) is the closest point of the
!                  equator to the rotated point.  Storage
!                  for XP, YP, and ZP may coincide with
!                  storage for X, Y, and Z, respectively,
!                  if the latter need not be saved.
!
! Modules required by APLYR:  None
!
! Intrinsic function called by APLYR:  SQRT
!
!***********************************************************
!
      REAL T
!
! Local parameter:
!
! T = Temporary variable
!
      T = SX*Y + CX*Z
      YP = CX*Y - SX*Z
      ZP = SY*X + CY*T
      XP = CY*X - SY*T
      IF (ZP .GE. 0.) RETURN
!
! Move (XP,YP,ZP) to the equator.
!
      T = SQRT(XP*XP + YP*YP)
      IF (T .EQ. 0.) GO TO 1
      XP = XP/T
      YP = YP/T
      RETURN
!
! Move the south pole to an arbitrary point of the equator.
!
    1 XP = 1.
      YP = 0.
      RETURN
      END
      SUBROUTINE APLYRT (G1P,G2P,CX,SX,CY,SY, G)
      REAL G1P, G2P, CX, SX, CY, SY, G(3)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   05/09/92
!
!   This subroutine applies the inverse (transpose) of the
! rotation defined by Subroutine CONSTR to the vector
! (G1P G2P 0)**T, i.e., the gradient (G1P,G2P,0) in the rot-
! ated coordinate system is mapped to (G1,G2,G3) in the
! original coordinate system.
!
! On input:
!
!       G1P,G2P = X and Y components, respectively, of the
!                 gradient in the rotated coordinate system.
!
!       CX,SX,CY,SY = Elements of the rotation R constructed
!                     by Subroutine CONSTR.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       G = X, Y, and Z components (in that order) of the
!           inverse rotation applied to (G1P,G2P,0) --
!           gradient in the original coordinate system.
!
! Modules required by APLYRT:  None
!
!***********************************************************
!
      REAL T
!
! Local parameters:
!
! T = Temporary variable
!
      T = SY*G1P
      G(1) = CY*G1P
      G(2) = CX*G2P - SX*T
      G(3) = -SX*G2P - CX*T
      RETURN
      END
      SUBROUTINE ARCINT (P,P1,P2,F1,F2,G1,G2,SIGMA, F,G,GN)
      REAL    P(3), P1(3), P2(3), F1, F2, G1(3), G2(3),&
     &        SIGMA, F, G(3), GN
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/21/96
!
!   Given 3 points P, P1, and P2 lying on a common geodesic
! of the unit sphere with P between P1 and P2, along with
! data values and gradients at P1 and P2, this subroutine
! computes an interpolated value F and a gradient vector G
! AT P.  F and the tangential component of G are taken to be
! the value and derivative (with respect to arc-length) of
! a Hermite interpolatory tension spline defined by the end-
! point values and tangential gradient components.  The nor-
! mal component of G is obtained by linear interpolation of
! the normal components of the gradients at P1 and P2.
!
! On input:
!
!       P = Cartesian coordinates of a point lying on the
!           arc defined by P1 and P2.  P(1)**2 + P(2)**2 +
!           P(3)**2 = 1.
!
!       P1,P2 = Coordinates of distinct points on the unit
!               sphere defining an arc with length less than
!               180 degrees.
!
!       F1,F2 = Data values associated with P1 and P2,
!               respectively.
!
!       G1,G2 = Gradient vectors associated with P1 and P2.
!               G1 and G2 are orthogonal to P1 and P2,
!               respectively.
!
!       SIGMA = Tension factor associated with P1-P2.
!
! The above parameters are not altered by this routine.
!
!       G = Array of length 3.
!
! On output:
!
!       F = Interpolated value at P.
!
!       G = Interpolated gradient at P.
!
!       GN = Normal component of G with the direction
!            P1 X P2 taken to be positive.  The extrapola-
!            tion procedure requires this component.
!
!   For each vector V, V(1), V(2), and V(3) contain X, Y,
! and Z components, respectively.
!
! SSRFPACK modules required by ARCINT:  ARCLEN, SNHCSH
!
! Intrinsic functions called by ARCINT:  ABS, EXP, SQRT
!
!***********************************************************
!
      REAL    ARCLEN
      INTEGER I, LUN
      REAL    A, AL, B1, B2, CM, CMM, CM2, DUMMY, D1, D2, E,&
     &        EMS, E1, E2, GT, S, SB1, SB2, SIG, SINH,&
     &        SINH2, SM, SM2, TAU1, TAU2, TM, TM1, TM2, TP1,&
     &        TP2, TS, UN(3), UNORM
      DATA    LUN/6/
!
! Local parameters:
!
! A =         Angle in radians (arc-length) between P1 and
!               P2
! AL =        Arc-length between P1 and P
! B1,B2 =     Local coordinates of P with respect to P1-P2
! CM,CMM =    Coshm(SIG) and Coshmm(SIG) -- refer to SNHCSH
! CM2 =       Coshm(SB2)
! DUMMY =     Dummy parameter for SNHCSH
! D1,D2 =     Scaled second differences
! E =         CM**2 - SM*Sinh = SIG*SM - 2*CMM (scaled by
!               2*EMS if SIG > .5)
! EMS =       Exp(-SIG)
! E1,E2 =     Exp(-SB1), Exp(-SB2)
! GT =        Tangential component of G -- component in the
!               direction UN X P
! I =         DO-loop index
! LUN =       Logical unit for error messages
! S =         Slope:  (F2-F1)/A
! SB1,SB2 =   SIG*B1, SIG*B2
! SIG =       Abs(SIGMA)
! SINH =      Sinh(SIGMA)
! SINH2 =     Sinh(SB2)
! SM,SM2 =    Sinhm(SIG), Sinhm(SB2)
! TAU1,TAU2 = Tangential derivatives (components of G1,G2)
!               at P1 and P2
! TM =        1-EMS
! TM1,TM2 =   1-E1, 1-E2
! TP1,TP2 =   1+E1, 1+E2
! TS =        TM**2
! UN =        Unit normal to the plane of P, P1, and P2
! UNORM =     Euclidean norm of P1 X P2 -- used to normalize
!               UN
!
!
! Compute unit normal UN.
!
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.) GO TO 2
!
! Normalize UN.
!
      DO 1 I = 1,3
        UN(I) = UN(I)/UNORM
    1   CONTINUE
!
! Compute tangential derivatives at the endpoints:
!   TAU1 = (G1,UN X P1) = (G1,P2)/UNORM and
!   TAU2 = (G2,UN X P2) = -(G2,P1)/UNORM.
!
      TAU1 = (G1(1)*P2(1) + G1(2)*P2(2) + G1(3)*P2(3))/UNORM
      TAU2 =-(G2(1)*P1(1) + G2(2)*P1(2) + G2(3)*P1(3))/UNORM
!
! Compute arc-lengths A, AL.
!
      A = ARCLEN(P1,P2)
      IF (A .EQ. 0.) GO TO 2
      AL = ARCLEN(P1,P)
!
! Compute local coordinates, slope, and second differences.
!
      B2 = AL/A
      B1 = 1. - B2
      S = (F2-F1)/A
      D1 = S - TAU1
      D2 = TAU2 - S
!
! Test the range of SIGMA.
!
      SIG = ABS(SIGMA)
      IF (SIG .LT. 1.E-9) THEN
!
! Hermite cubic interpolation.
!
        F = F1 + AL*(TAU1 + B2*(D1 + B1*(D1 - D2)))
        GT = TAU1 + B2*(D1 + D2 + 3.*B1*(D1 - D2))
      ELSEIF (SIG .LE. .5) THEN
!
! 0 < SIG .LE. .5.  Use approximations designed to avoid
!   cancellation error in the hyperbolic functions.
!
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        SINH = SM + SIG
        SINH2 = SM2 + SB2
        E = SIG*SM - CMM - CMM
        F = F1 + AL*TAU1 + A*((CM*SM2-SM*CM2)*(D1+D2) + SIG*&
     &                        (CM*CM2-SINH*SM2)*D1)/(SIG*E)
        GT = TAU1 + ((CM*CM2-SM*SINH2)*(D1+D2) + SIG*&
     &               (CM*SINH2-SINH*CM2)*D1)/E
      ELSE
!
! SIG > .5.  Use negative exponentials in order to avoid
!   overflow.  Note that EMS = EXP(-SIG).
!
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
        F = F1 + AL*S + A*(TM*TM1*TM2*(D1+D2) + SIG*&
     &                     ((E2*TM1*TM1-B1*TS)*D1 +&
     &                      (E1*TM2*TM2-B2*TS)*D2))/(SIG*E)
        TP1 = 1. + E1
        TP2 = 1. + E2
        GT = S + (TM1*(TM*TP2-SIG*E2*TP1)*D1 -&
     &            TM2*(TM*TP1-SIG*E1*TP2)*D2)/E
      ENDIF
!
! Compute GN.
!
      GN = B1*(UN(1)*G1(1) + UN(2)*G1(2) + UN(3)*G1(3)) +&
     &     B2*(UN(1)*G2(1) + UN(2)*G2(2) + UN(3)*G2(3))
!
! Compute G = GT*(UN X P) + GN*UN.
!
      G(1) = GT*(UN(2)*P(3) - UN(3)*P(2)) + GN*UN(1)
      G(2) = GT*(UN(3)*P(1) - UN(1)*P(3)) + GN*UN(2)
      G(3) = GT*(UN(1)*P(2) - UN(2)*P(1)) + GN*UN(3)
      RETURN
!
! P1 X P2 = 0.  Print an error message and terminate
!   processing.
!
    2 WRITE (LUN,100) (P1(I),I=1,3), (P2(I),I=1,3)
  100 FORMAT ('1','ERROR IN ARCINT -- P1 = ',2(F9.6,',  '),&
     &        F9.6/1X,19X,'P2 = ',2(F9.6,',  '),F9.6)
      STOP
      END
      REAL FUNCTION ARCLEN (P,Q)
      REAL P(3), Q(3)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   05/09/92
!
!   This function computes the arc-length (angle in radians)
! between a pair of points on the unit sphere.
!
! On input:
!
!       P,Q = Arrays of length 3 containing the X, Y, and Z
!             coordinates (in that order) of points on the
!             unit sphere.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       ARCLEN = Angle in radians between the unit vectors
!                P and Q.  0 .LE. ARCLEN .LE. PI.
!
! Modules required by ARCLEN:  None
!
! Intrinsic functions called by ARCLEN:  ATAN, SQRT
!
!***********************************************************
!
      INTEGER I
      REAL    D
!
! Local parameters:
!
! D = Euclidean norm squared of P+Q
! I = DO-loop index
!
      D = 0.
      DO 1 I = 1,3
        D = D + (P(I) + Q(I))**2
    1   CONTINUE
      IF (D .EQ. 0.) THEN
!
! P and Q are separated by 180 degrees.
!
        ARCLEN = 4.*ATAN(1.)
      ELSEIF (D .GE. 4.) THEN
!
! P and Q coincide.
!
        ARCLEN = 0.
      ELSE
        ARCLEN = 2.*ATAN(SQRT((4.-D)/D))
      ENDIF
      RETURN
      END
      SUBROUTINE CONSTR (XK,YK,ZK, CX,SX,CY,SY)
      REAL XK, YK, ZK, CX, SX, CY, SY
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   05/09/92
!
!   This subroutine constructs the elements of a 3 by 3
! orthogonal matrix R which rotates a point (XK,YK,ZK) on
! the unit sphere to the north pole, i.e.,
!
!      (XK)     (CY  0 -SY)   (1   0   0)   (XK)     (0)
!  R * (YK)  =  ( 0  1   0) * (0  CX -SX) * (YK)  =  (0)
!      (ZK)     (SY  0  CY)   (0  SX  CX)   (ZK)     (1)
!
! On input:
!
!       XK,YK,ZK = Components of a unit vector to be
!                  rotated to (0,0,1).
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       CX,SX,CY,SY = Elements of R:  CX,SX define a rota-
!                     tion about the X-axis and CY,SY define
!                     a rotation about the Y-axis.
!
! Modules required by CONSTR:  None
!
! Intrinsic function called by CONSTR:  SQRT
!
!***********************************************************
!
      CY = SQRT(YK*YK + ZK*ZK)
      SY = XK
      IF (CY .NE. 0.) THEN
        CX = ZK/CY
        SX = YK/CY
      ELSE
!
! (XK,YK,ZK) lies on the X-axis.
!
        CX = 1.
        SX = 0.
      ENDIF
      RETURN
      END
      REAL FUNCTION FVAL (B1,B2,B3,V1,V2,V3,F1,F2,F3,G1,G2,&
     &                    G3,SIG1,SIG2,SIG3)
      REAL B1, B2, B3, V1(3), V2(3), V3(3), F1, F2, F3,&
     &     G1(3), G2(3), G3(3), SIG1, SIG2, SIG3
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   05/09/92
!
!   Given data values and gradients at the three vertices of
! a spherical triangle containing a point P, this routine
! computes the value of F at P where F interpolates the ver-
! tex data.  Along the triangle sides, the interpolatory
! function F is the Hermite interpolatory tension spline
! defined by the values and tangential gradient components
! at the endpoints, and the gradient component normal to the
! triangle side varies linearly with respect to arc-length
! between the normal gradient components at the endpoints.
! A first-order C-1 blending method is used on the underly-
! ing planar triangle.  Since values and gradients on an arc
! depend only on the vertex data, the method results in C-1
! continuity when used to interpolate over a triangulation.
!
!   The blending method consists of taking F(P) to be a
! weighted sum of the values at PP of three univariate Her-
! mite interpolatory tension splines defined on the line
! segments which join the vertices to the opposite sides and
! pass through PP:  the central projection of P onto the
! underlying planar triangle.  The tension factors for these
! splines are obtained by linear interpolation between the
! pair of tension factors associated with the triangle sides
! which join at the appropriate vertex.
!
!   A tension factor SIGMA associated with a Hermite interp-
! olatory tension spline is a nonnegative parameter which
! determines the curviness of the spline.  SIGMA = 0 results
! in a cubic spline, and the spline approaches the linear
! interpolant as SIGMA increases.
!
! On input:
!
!       B1,B2,B3 = Barycentric coordinates of PP with re-
!                  spect to the (planar) underlying triangle
!                  (V1,V2,V3), where PP is the central
!                  projection of P onto this triangle.
!
!       V1,V2,V3 = Cartesian coordinates of the vertices of
!                  a spherical triangle containing P.  V3
!                  Left V1->V2.
!
!       F1,F2,F3 = Data values associated with the vertices.
!
!       G1,G2,G3 = Gradients associated with the vertices.
!                  Gi is orthogonal to Vi for i = 1,2,3.
!
!       SIG1,SIG2,SIG3 = Tension factors associated with the
!                        triangle sides opposite V1, V2, and
!                        V3, respectively.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       FVAL = Interpolated value at P.
!
! Each vector V above contains X, Y, and Z components in
!   V(1), V(2), and V(3), respectively.
!
! SSRFPACK modules required by FVAL:  ARCINT, ARCLEN, HVAL
!
! Intrinsic function called by FVAL:  SQRT
!
!***********************************************************
!
      REAL    HVAL
      INTEGER I
      REAL    C1, C2, C3, DS, DUM, DV, F, G(3),&
     &        Q1(3), Q2(3), Q3(3), SIG, SUM, S1, S2, S3,&
     &        U1(3), U2(3), U3(3), U1N, U2N, U3N, VAL
!
! Local parameters:
!
! C1,C2,C3 =    Coefficients (weight functions) of partial
!                 interpolants.  C1 = 1 on the edge opposite
!                 V1 and C1 = 0 on the other edges.  Simi-
!                 larly for C2 and C3.  C1+C2+C3 = 1.
! DS =          Directional derivative (scaled by distnace)
!                 at U1, U2, or U3:  DS = (G,U1-V1)/U1N =
!                 -(G,V1)/U1N on side opposite V1, where G/
!                 U1N (plus an orthogonal component) is the
!                 projection of G onto the planar triangle
! DUM =         Dummy variable for calls to ARCINT
! DV =          Directional derivatives (scaled by distance)
!                 at a vertex:  D1 = (G1,U1-V1) = (G1,U1)
! F,G =         Value and gradient at Q1 Q2, or Q3 obtained
!                 by interpolation along one of the arcs of
!                 the spherical triangle
! I =           DO-loop index
! Q1,Q2,Q3 =    Central projections of U1, U2, and U3 onto
!                 the sphere and thus lying on an arc of the
!                 spherical triangle
! SIG =         Tension factor for a side-vertex (partial)
!                 interpolant:  obtained by linear interpo-
!                 lation applied to triangle side tensions
! SUM =         Quantity used to normalize C1, C2, and C3
! S1,S2,S3 =    Sums of pairs of barycentric coordinates:
!                 used to compute U1, U2, U3, and SIG
! U1,U2,U3 =    Points on the boundary of the planar trian-
!                 gle and lying on the lines containing PP
!                 and the vertices.  U1 is opposite V1, etc.
! U1N,U2N,U3N = Quantities used to compute Q1, Q2, and Q3
!                 (magnitudes of U1, U2, and U3)
! VAL =         Local variable used to accumulate the con-
!                 tributions to FVAL
!
!
! Compute weight functions C1, C2, and C3.
!
      C1 = B2*B3
      C2 = B3*B1
      C3 = B1*B2
      SUM = C1 + C2 + C3
      IF (SUM .LE. 0.) THEN
!
! P coincides with a vertex.
!
        FVAL = B1*F1 + B2*F2 + B3*F3
        RETURN
      ENDIF
!
! Normalize C1, C2, and C3.
!
      C1 = C1/SUM
      C2 = C2/SUM
      C3 = C3/SUM
!
! Compute (S1,S2,S3), (U1,U2,U3) and (U1N,U2N,U3N).
!
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
!
! Compute Q1, Q2, and Q3.
!
      U1N = SQRT(U1N)
      U2N = SQRT(U2N)
      U3N = SQRT(U3N)
      DO 2 I = 1,3
        Q1(I) = U1(I)/U1N
        Q2(I) = U2(I)/U2N
        Q3(I) = U3(I)/U3N
    2   CONTINUE
!
! Compute interpolated value (VAL) at P by looping on
!   triangle sides.
!
      VAL = 0.
!
! Contribution from side opposite V1:
!
!   Compute value and gradient at Q1 by interpolating
!     between V2 and V3.
!
      CALL ARCINT (Q1,V2,V3,F2,F3,G2,G3,SIG1, F,G,DUM)
!
!   Add in the contribution.
!
      DV = G1(1)*U1(1) + G1(2)*U1(2) + G1(3)*U1(3)
      DS = -(G(1)*V1(1) + G(2)*V1(2) + G(3)*V1(3))/U1N
      SIG = (B2*SIG3 + B3*SIG2)/S1
      VAL = VAL + C1*HVAL(B1,F1,F,DV,DS,SIG)
!
! Contribution from side opposite V2:
!
!   Compute value and gradient at Q2 by interpolating
!     between V3 and V1.
!
      CALL ARCINT (Q2,V3,V1,F3,F1,G3,G1,SIG2, F,G,DUM)
!
!   Add in the contribution.
!
      DV = G2(1)*U2(1) + G2(2)*U2(2) + G2(3)*U2(3)
      DS = -(G(1)*V2(1) + G(2)*V2(2) + G(3)*V2(3))/U2N
      SIG = (B3*SIG1 + B1*SIG3)/S2
      VAL = VAL + C2*HVAL(B2,F2,F,DV,DS,SIG)
!
! Contribution from side opposite V3:
!
!   Compute interpolated value and gradient at Q3
!     by interpolating between V1 and V2.
!
      CALL ARCINT (Q3,V1,V2,F1,F2,G1,G2,SIG3, F,G,DUM)
!
!   Add in the final contribution.
!
      DV = G3(1)*U3(1) + G3(2)*U3(2) + G3(3)*U3(3)
      DS = -(G(1)*V3(1) + G(2)*V3(2) + G(3)*V3(3))/U3N
      SIG = (B1*SIG2 + B2*SIG1)/S3
      FVAL = VAL + C3*HVAL(B3,F3,F,DV,DS,SIG)
      RETURN
      END
      SUBROUTINE GETSIG (N,X,Y,Z,H,LIST,LPTR,LEND,GRAD,&
     &                   TOL, SIGMA, DSMAX,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IER
      REAL    X(N), Y(N), Z(N), H(N), GRAD(3,N), TOL,&
     &        SIGMA(*), DSMAX
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/21/96
!
!   Given a triangulation of a set of nodes on the unit
! sphere, along with data values H and gradients GRAD at the
! nodes, this subroutine determines, for each triangulation
! arc, the smallest (nonnegative) tension factor SIGMA such
! that the Hermite interpolatory tension spline H(A), de-
! fined by SIGMA and the endpoint values and directional
! derivatives, preserves local shape properties of the data.
! In order to define the shape properties on an arc, it is
! convenient to map the arc to an interval (A1,A2).  Then,
! denoting the endpoint data values by H1,H2 and the deriva-
! tives (tangential gradient components) by HP1,HP2, and
! letting S = (H2-H1)/(A2-A1), the data properties are
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
!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of the nodes.
!
!       H = Array of length N containing data values at the
!           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
!           for I = 1,...,N.
!
!       LIST,LPTR,LEND = Data structure defining the tri-
!                        angulation.  Refer to STRIPACK
!                        Subroutine TRMESH.
!
!       GRAD = Array dimensioned 3 by N whose columns con-
!              tain gradients at the nodes.  GRAD( ,J) must
!              be orthogonal to node J:  GRAD(1,J)*X(J) +
!              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0..  Refer
!              to Subroutines GRADG, GRADL, and SMSURF.
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
!               H(A) preserves the local data properties on
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
! STRIPACK modules required by GETSIG:  LSTPTR, STORE
!
! SSRFPACK modules required by GETSIG:  ARCLEN, SNHCSH
!
! Intrinsic functions called by GETSIG:  ABS, EXP, MAX, MIN,
!                                          SIGN, SQRT
!
!***********************************************************
!
      INTEGER LSTPTR
      REAL    ARCLEN, STORE
      INTEGER ICNT, LP1, LP2, LPL, LUN, N1, N2, NIT, NM1
      REAL    A, AL, C1, C2, COSHM, COSHMM, D0, D1, D1D2,&
     &        D1PD2, D2, DMAX, DSIG, DSM, E, EMS, EMS2, F,&
     &        F0, FMAX, FNEG, FP, FTOL, P1(3), P2(3), RTOL,&
     &        S, S1, S2, SBIG, SCM, SGN, SIG, SIGIN, SINHM,&
     &        SSINH, SSM, STOL, T, T0, T1, T2, TM, TP1,&
     &        UN(3), UNORM
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
  100 FORMAT ('1',13X,'GETSIG -- N =',I4,', TOL = ',E10.3//)
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
! Print a message and compute parameters for the arc:
!   nodal coordinates P1 and P2, arc-length AL,
!   UNORM = magnitude of P1 X P2, and
!   SIGIN = input SIGMA value.
!
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
!
! Compute scaled directional derivatives S1,S2 at the end-
!   points (for the direction N1->N2), first difference S,
!   and second differences D1,D2.
!
        S1 = AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +&
     &               GRAD(3,N1)*P2(3))/UNORM
        S2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +&
     &            GRAD(3,N2)*P1(3))/UNORM
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
! Convexity:  find a zero of F(SIG) = SIG*coshm(SIG)/
!   sinhm(SIG) - TP1.
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
!   Scale sinhm and coshm by 2*exp(-SIG) in order to avoid
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
  120   FORMAT (1X,'CONVEXITY -- SIG = ',E15.8,&
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
  130   FORMAT (1X,'MONOTONICITY -- DSIG = ',E15.8)
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
!   Scale sinhm and coshm by 2*exp(-SIG) in order to avoid
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
! N < 3.
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
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   05/09/92
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
      IF (ABS(AA) .GT. ABS(BB)) THEN
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
      ELSEIF (BB .NE. 0.) THEN
!
! ABS(A) .LE. ABS(B).
!
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
      ELSE
!
! A = B = 0.
!
        C = 1.
        S = 0.
      ENDIF
      RETURN
      END



      SUBROUTINE GRADG (N,X,Y,Z,F,LIST,LPTR,LEND,IFLGS,&
     &                  SIGMA, NIT,DGMAX,GRAD, IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, NIT, IER
      REAL    X(N), Y(N), Z(N), F(N), SIGMA(*), DGMAX,&
     &        GRAD(3,N)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/24/96
!
!   Given a triangulation of N nodes on the unit sphere with
! data values F at the nodes and tension factors SIGMA asso-
! ciated with the arcs, this routine uses a global method
! to compute estimated gradients at the nodes.  The method
! consists of minimizing a quadratic functional Q(G) over
! the N-vector G of gradients, where Q approximates the
! linearized curvature of the restriction to arcs of the
! interpolatory function F defined by Function FVAL.  The
! restriction of F to an arc of the triangulation is the
! Hermite interpolatory tension spline defined by the data
! values and tangential gradient components at the endpoints
! of the arc.  Letting D1F(A) and D2F(A) denote first and
! second derivatives of F with respect to a parameter A var-
! ying along a triangulation arc, Q is the sum of integrals
! over the arcs of D2F(A)**2 + ((SIGMA/L)*(D1F(A)-S))**2,
! where L denotes arc-length, SIGMA is the appropriate ten-
! sion factor, and S is the slope of the linear function of
! A which interpolates the values of F at the endpoints of
! the arc.
!
!   Since the gradient at node K lies in the plane tangent
! to the sphere surface at K, it is effectively defined by
! only two components -- its X and Y components in the coor-
! dinate system obtained by rotating K to the north pole.
! Thus, the minimization problem corresponds to an order-2N
! symmetric positive-definite sparse linear system which is
! solved by a block Gauss-Seidel method with 2 by 2 blocks.
!
!   An alternative method, Subroutine GRADL, computes a
! local approximation to the gradient at a single node and,
! although less efficient when all gradients are needed, was
! found to be generally more accurate (in the case of uni-
! form zero tension) when the nodal distribution is very
! dense, varies greatly, or does not cover the sphere.
! GRADG, on the other hand, was found to be slightly more
! accurate on a uniform distribution of 514 nodes.
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y,Z = Arrays of length N containing Cartesian
!               coordinates of the nodes.  X(I)**2 + Y(I)**2
!               + Z(I)**2 = 1 for I = 1,...,N.
!
!       F = Array of length N containing data values at the
!           nodes.  F(I) is associated with (X(I),Y(I),Z(I))
!           for I = 1,...,N.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
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
!               programs GETSIG, SIG0, SIG1, and SIG2.
!
! The above parameters are not altered by this routine.
!
!       NIT = Maximum number of Gauss-Seidel iterations to
!             be applied.  This maximum will likely be a-
!             chieved if DGMAX is smaller than the machine
!             precision.  NIT .GE. 0.
!
!       DGMAX = Nonnegative convergence criterion.  The
!               method is terminated when the maximum change
!               in a gradient between iterations is at most
!               DGMAX.  The change in a gradient is taken to
!               be the Euclidean norm of the difference (in
!               the rotated coordinate system) relative to 1
!               plus the norm of the old gradient value.
!
!       GRAD = 3 by N array whose columns contain initial
!              solution estimates (zero vectors are suffici-
!              ent).  GRAD(I,J) contains component I of the
!              gradient at node J for I = 1,2,3 (X,Y,Z) and
!              J = 1,...,N.  GRAD( ,J) must be orthogonal to
!              node J -- GRAD(1,J)*X(J) + GRAD(2,J)*Y(J) +
!              GRAD(3,J)*Z(J) = 0.
!
! On output:
!
!       NIT = Number of Gauss-Seidel iterations employed.
!
!       DGMAX = Maximum change in a gradient at the last
!               iteration.
!
!       GRAD = Estimated gradients.  See the description
!              under input parameters.  GRAD is not changed
!              if IER = -1.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and the
!                     convergence criterion was achieved.
!             IER = 1 if no errors were encountered but con-
!                     vergence was not achieved within NIT
!                     iterations.
!             IER = -1 if N or DGMAX is outside its valid
!                      range or NIT .LT. 0 on input.
!             IER = -2 if all nodes are collinear or the
!                      triangulation is invalid.
!             IER = -3 if duplicate nodes were encountered.
!
! SSRFPACK modules required by GRADG:  APLYRT, CONSTR,
!                                        GRCOEF, SNHCSH
!
! Intrinsic functions called by GRADG:  ATAN, MAX, SQRT
!
!***********************************************************
!
      INTEGER IFL, ITER, J, K, LPJ, LPL, MAXIT, NN
      REAL    ALFA, A11, A12, A22, CX, CY, D, DEN, DET,&
     &        DGK(3), DGMX, DG1, DG2, FK, G1, G2, G3, R1,&
     &        R2, SD, SIG, SINAL, SX, SY, T, TOL, XK, YK,&
     &        ZK, XJ, YJ, ZJ, XS, YS
!
! Local parameters:
!
! ALFA =        Arc-length between nodes K and J
! A11,A12,A22 = Matrix components of the 2 by 2 block A*DG
!                 = R where A is symmetric, (DG1,DG2,0) is
!                 the change in the gradient at K, and R is
!                 the residual
! CX,CY =       Components of a rotation mapping K to the
!                 north pole (0,0,1)
! D =           Function of SIG computed by GRCOEF -- factor
!                 in the order-2 system
! DEN =         ALFA*SINAL**2 -- factor in the 2 by 2 system
! DET =         Determinant of the order-2 matrix
! DGK =         Change in GRAD( ,K) from the previous esti-
!                 mate in the original coordinate system
! DGMX =        Maximum change in a gradient between itera-
!                 tions
! DG1,DG2 =     Solution of the 2 by 2 system -- first 2
!                 components of DGK in the rotated coordi-
!                 nate system
! FK =          Data value F(K)
! G1,G2,G3 =    Components of GRAD( ,K)
! IFL =         Local copy of IFLGS
! ITER =        Number of iterations used
! J =           Neighbor of K
! K =           DO-loop and node index
! LPJ =         LIST pointer of node J as a neighbor of K
! LPL =         Pointer to the last neighbor of K
! MAXIT =       Input value of NIT
! NN =          Local copy of N
! R1,R2 =       Components of the residual -- derivatives of
!                 Q with respect to the components of the
!                 gradient at node K
! SD =          Function of SIG computed by GRCOEF -- factor
!                 in the order-2 system
! SIG =         Tension factor associated with ARC K-J
! SINAL =       SIN(ALFA) -- magnitude of the vector cross
!                 product between nodes K and J
! SX,SY =       Components of a rotation mapping K to the
!                 north pole (0,0,1)
! T =           Temporary storage for factors in the system
!                 components
! TOL =         Local copy of DGMAX
! XK,YK,ZK =    Coordinates of node K -- X(K), Y(K), Z(K)
! XJ,YJ,ZJ =    Coordinates of node J in the rotated coor-
!                 dinate system
! XS,YS =       XJ**2, YJ**2
!
      NN = N
      IFL = IFLGS
      MAXIT = NIT
      TOL = DGMAX
!
! Test for errors in input, and initialize iteration count,
!   tension factor, and output value of DGMAX.
!
      IF (NN .LT. 3  .OR.  MAXIT .LT. 0  .OR.  TOL .LT. 0.)&
     &   GO TO 11
      ITER = 0
      SIG = SIGMA(1)
      DGMX = 0.
!
! Top of iteration loop.
!
    1 IF (ITER .EQ. MAXIT) GO TO 4
      DGMX = 0.
!
! Loop on nodes.
!
      DO 3 K = 1,NN
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        FK = F(K)
        G1 = GRAD(1,K)
        G2 = GRAD(2,K)
        G3 = GRAD(3,K)
!
!   Construct the rotation mapping node K to the north pole.
!
        CALL CONSTR (XK,YK,ZK, CX,SX,CY,SY)
!
!   Initialize components of the 2 by 2 system for the
!     change (DG1,DG2,0) in the K-th solution components
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
    2   LPJ = LPTR(LPJ)
          J = ABS(LIST(LPJ))
!
!   Compute the coordinates of J in the rotated system.
!
          T = SX*Y(J) + CX*Z(J)
          YJ = CX*Y(J) - SX*Z(J)
          ZJ = SY*X(J) + CY*T
          XJ = CY*X(J) - SY*T
!
!   Compute arc-length ALFA between J and K, SINAL =
!     SIN(ALFA), and DEN = ALFA*SIN(ALFA)**2.
!
          ALFA = 2.*ATAN(SQRT((1.-ZJ)/(1.+ZJ)))
          XS = XJ*XJ
          YS = YJ*YJ
          SINAL = SQRT(XS+YS)
          DEN = ALFA*(XS+YS)
!
!   Test for coincident nodes and compute functions of SIG:
!     D = SIG*(SIG*COSHM-SINHM)/E and SD = SIG*SINHM/E for
!     E = SIG*SINH-2*COSHM.
!
          IF (DEN .EQ. 0.) GO TO 13
          IF (IFL .GE. 1) SIG = SIGMA(LPJ)
          CALL GRCOEF (SIG, D,SD)
!
!   Update the system components for node J.
!
          T = D/DEN
          A11 = A11 + T*XS
          A12 = A12 + T*XJ*YJ
          A22 = A22 + T*YS
          T = (D+SD)*(FK-F(J))/(ALFA*ALFA*SINAL) +&
     &        ( D*(G1*X(J) + G2*Y(J) + G3*Z(J)) -&
     &          SD*(GRAD(1,J)*XK + GRAD(2,J)*YK +&
     &                 GRAD(3,J)*ZK) )/DEN
          R1 = R1 - T*XJ
          R2 = R2 - T*YJ
!
!   Bottom of loop on neighbors.
!
          IF (LPJ .NE. LPL) GO TO 2
!
!   Solve the 2 by 2 system and update DGMAX.
!
        DET = A11*A22 - A12*A12
        IF (DET .EQ. 0.  .OR.  A11 .EQ. 0.) GO TO 12
        DG2 = (A11*R2 - A12*R1)/DET
        DG1 = (R1 - A12*DG2)/A11
        DGMX = MAX(DGMX,SQRT(DG1*DG1+DG2*DG2)/&
     &             (1.+SQRT(G1*G1+G2*G2+G3*G3)))
!
!   Rotate (DG1,DG2,0) back to the original coordinate
!     system and update GRAD( ,K).
!
        CALL APLYRT (DG1,DG2,CX,SX,CY,SY, DGK)
        GRAD(1,K) = G1 + DGK(1)
        GRAD(2,K) = G2 + DGK(2)
        GRAD(3,K) = G3 + DGK(3)
    3   CONTINUE
!
!   Increment ITER and test for convergence.
!
      ITER = ITER + 1
      IF (DGMX .GT. TOL) GO TO 1
!
! The method converged.
!
      NIT = ITER
      DGMAX = DGMX
      IER = 0
      RETURN
!
! The method failed to converge within NIT iterations.
!
    4 DGMAX = DGMX
      IER = 1
      RETURN
!
! Invalid input parameter.
!
   11 NIT = 0
      DGMAX = 0.
      IER = -1
      RETURN
!
! Node K and its neighbors are collinear.
!
   12 NIT = 0
      DGMAX = DGMX
      IER = -2
      RETURN
!
! Nodes K and J coincide.
!
   13 NIT = 0
      DGMAX = DGMX
      IER = -3
      RETURN
      END
      SUBROUTINE GRADL (N,K,X,Y,Z,W,LIST,LPTR,LEND, G,IER)
      INTEGER N, K, LIST(*), LPTR(*), LEND(N), IER
      REAL    X(N), Y(N), Z(N), W(N), G(3)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/24/96
!
!   Given a triangulation of a set of nodes on the unit
! sphere with their associated data values W, this routine
! estimates a gradient vector at node K as follows:  the
! coordinate system is rotated so that K becomes the north
! pole, node K and a set of nearby nodes are projected
! orthogonally onto the X-Y plane (in the new coordinate
! system), a quadratic is fitted in a weighted least squares
! sense to the data values at the projected nodes such that
! the value (associated with K) at (0,0) is interpolated, X
! and Y-partial derivative estimates DX and DY are computed
! by differentiating the quadratic at (0,0), and the esti-
! mated gradient G is obtained by rotating (DX,DY,0) back to
! the original coordinate system.  Note that G lies in the
! plane tangent to the sphere at node K, i.e., G is orthogo-
! nal to the unit vector represented by node K.  A Marquardt
! stabilization factor is used if necessary to ensure a
! well-conditioned least squares system, and a unique solu-
! tion exists unless the nodes are collinear.
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 7.
!
!       K = Node at which the gradient is sought.  1 .LE. K
!           .LE. N.
!
!       X,Y,Z = Arrays containing the Cartesian coordinates
!               of the nodes.
!
!       W = Array containing the data values at the nodes.
!           W(I) is associated with (X(I),Y(I),Z(I)) for
!           I = 1,...,N.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
!                        Subroutine TRMESH.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       G = X, Y, and Z components (in that order) of the
!           estimated gradient at node K unless IER < 0.
!
!       IER = Error indicator:
!             IER .GE. 6 if no errors were encountered.
!                        IER contains the number of nodes
!                        (including K) used in the least
!                        squares fit.
!             IER = -1 if N or K is outside its valid range.
!             IER = -2 if the least squares system has no
!                      unique solution due to duplicate or
!                      collinear nodes.
!
! STRIPACK module required by GRADL:  GETNP
!
! SSRFPACK modules required by GRADL:  APLYR, APLYRT,
!                                        CONSTR, GIVENS,
!                                        ROTATE, SETUP
!
! Intrinsic functions called by GRADL:  ABS, MIN, REAL, SQRT
!
!***********************************************************
!
      INTEGER   LMN, LMX
      PARAMETER (LMN=10,  LMX=30)
      INTEGER I, IERR, IM1, IP1, J, JP1, KK, L, LM1, LMAX,&
     &        LMIN, LNP, NN, NP, NPTS(LMX)
      REAL    A(6,6), AV, AVSQ, C, CX, CY, DF, DMIN, DTOL,&
     &        DX, DY, RF, RIN, RTOL, S, SF, SUM, SX, SY,&
     &        WK, WT, XP, YP, ZP
!
      DATA    RTOL/1.E-6/, DTOL/.01/, SF/1./
!
! Local parameters:
!
! A =         Transpose of the (upper triangle of the) aug-
!               mented regression matrix
! AV =        Root-mean-square distance (in the rotated
!               coordinate system) between the origin and
!               the nodes (other than K) in the least
!               squares fit.  The first 3 columns of A**T
!               are scaled by 1/AVSQ, the next 2 by 1/AV.
! AVSQ =      AV*AV:  accumulated in SUM
! C,S =       Components of the plane rotation used to
!               triangularize the regression matrix
! CX,SX =     Components of a plane rotation about the X-
!               axis which, together with CY and SY, define
!               a mapping from node K to the north pole
!               (0,0,1)
! CY,SY =     Components of a plane rotation about the Y-
!               axis
! DF =        Negative Z component (in the rotated coordi-
!               nate system) of an element NP of NPTS --
!               increasing function of the angular distance
!               between K and NP.  DF lies in the interval
!               (-1,1).
! DMIN =      Minimum of the magnitudes of the diagonal
!               elements of the triangularized regression
!               matrix
! DTOL =      Tolerance for detecting an ill-conditioned
!               system (DMIN is required to be at least
!               DTOL)
! DX,DY =     X and Y components of the estimated gradient
!               in the rotated coordinate system
! I,J =       Loop indexes
! IERR =      Error flag for calls to GETNP (not checked)
! IM1,IP1 =   I-1, I+1
! JP1 =       J+1
! KK =        Local copy of K
! L =         Number of columns of A**T to which a rotation
!               is applied
! LM1 =       LMIN-1
! LMIN,LMAX = Min(LMN,N), Min(LMX,N)
! LMN,LMX =   Minimum and maximum values of LNP for N
!               sufficiently large.  In most cases LMN-1
!               nodes are used in the fit.  7 .LE. LMN .LE.
!               LMX.
! LNP =       Length of NPTS or LMAX+1
! NN =        Local copy of N
! NP =        Element of NPTS to be added to the system
! NPTS =      Array containing the indexes of a sequence of
!               nodes ordered by angular distance from K.
!               NPTS(1)=K and the first LNP-1 elements of
!               NPTS are used in the least squares fit.
!               unless LNP = LMAX+1, NPTS(LNP) determines R
!               (see RIN).
! RF =        Value of DF associated with NPTS(LNP) unless
!               LNP = LMAX+1 (see RIN)
! RIN =       Inverse of a radius of influence R which
!               enters into WT:  R = 1+RF unless all ele-
!               ments of NPTS are used in the fit (LNP =
!               LMAX+1), in which case R is the distance
!               function associated with some point more
!               distant from K than NPTS(LMAX)
! RTOL =      Tolerance for determining LNP (and hence R):
!               if the increase in DF between two successive
!               elements of NPTS is less than RTOL, they are
!               treated as being the same distance from node
!               K and an additional node is added
! SF =        Marquardt stabilization factor used to damp
!               out the first 3 solution components (second
!               partials of the quadratic) when the system
!               is ill-conditioned.  Increasing SF results
!               in more damping (a more nearly linear fit).
! SUM =       Sum of squared Euclidean distances (in the
!               rotated coordinate system) between the
!               origin and the nodes used in the least
!               squares fit
! WK =        W(K) -- data value at node K
! WT =        Weight for the equation coreesponding to NP:
!               WT = (R-D)/(R*D) = 1/D - RIN, where D = 1-ZP
!               is associated with NP
! XP,YP,ZP =  Coordinates of NP in the rotated coordinate
!               system unless ZP < 0, in which case
!               (XP,YP,0) lies on the equator
!
      NN = N
      KK = K
      WK = W(KK)
!
! Check for errors and initialize LMIN, LMAX.
!
      IF (NN .LT. 7  .OR.  KK .LT. 1  .OR.  KK .GT. NN)&
     &   GO TO 13
      LMIN = MIN(LMN,NN)
      LMAX = MIN(LMX,NN)
!
! Compute NPTS, LNP, AVSQ, AV, and R.
!   Set NPTS to the closest LMIN-1 nodes to K.  DF contains
!   the negative Z component (in the rotated coordinate
!   system) of the new node on return from GETNP.
!
      SUM = 0.
      NPTS(1) = KK
      LM1 = LMIN - 1
      DO 1 LNP = 2,LM1
        CALL GETNP (X,Y,Z,LIST,LPTR,LEND,LNP, NPTS, DF,IERR)
        SUM = SUM + 1. - DF*DF
    1   CONTINUE
!
!   Add additional nodes to NPTS until the increase in
!     R = 1+RF is at least RTOL.
!
      DO 2 LNP = LMIN,LMAX
        CALL GETNP (X,Y,Z,LIST,LPTR,LEND,LNP, NPTS, RF,IERR)
        IF (RF-DF .GE. RTOL) GO TO 3
        SUM = SUM + 1. - RF*RF
    2   CONTINUE
!
!   Use all LMAX nodes in the least squares fit.  R is
!     arbitrarily increased by 5 percent.
!
      RF = 1.05*RF + .05
      LNP = LMAX + 1
!
!   There are LNP-2 equations corresponding to nodes
!     NPTS(2),...,NPTS(LNP-1).
!
    3 AVSQ = SUM/REAL(LNP-2)
      AV = SQRT(AVSQ)
      RIN = 1./(1.+RF)
!
! Construct the rotation.
!
      CALL CONSTR (X(KK),Y(KK),Z(KK), CX,SX,CY,SY)
!
! Set up the first 5 equations of the augmented regression
!   matrix (transposed) as the columns of A, and zero out
!   the lower triangle (upper triangle of A) with Givens
!   rotations.
!
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
!
! Add the additional equations to the system using
!   the last column of A.  I .LE. LNP.
!
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
!
! Test the system for ill-conditioning.
!
    8 DMIN = MIN( ABS(A(1,1)),ABS(A(2,2)),ABS(A(3,3)),&
     &            ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN .GE. DTOL) GO TO 12
      IF (LNP .LE. LMAX) THEN
!
! Add another node to the system and increase R.
!   I = LNP.
!
        LNP = LNP + 1
        IF (LNP .LE. LMAX) CALL GETNP (X,Y,Z,LIST,LPTR,LEND,&
     &                                 LNP,NPTS, RF,IERR)
        RIN = 1./(1.05*(1.+RF))
        GO TO 6
      ENDIF
!
! Stabilize the system by damping second partials.  Add
!   multiples of the first three unit vectors to the first
!   three equations.
!
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
!
! Test the linear portion of the stabilized system for
!   ill-conditioning.
!
      DMIN = MIN( ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN .LT. DTOL) GO TO 14
!
! Solve the 2 by 2 triangular system for the estimated
!   partial derivatives.
!
   12 DY = A(6,5)/A(5,5)
      DX = (A(6,4) - A(5,4)*DY)/A(4,4)/AV
      DY = DY/AV
!
! Rotate the gradient (DX,DY,0) back into the original
!   coordinate system.
!
      CALL APLYRT (DX,DY,CX,SX,CY,SY, G)
      IER = LNP - 1
      RETURN
!
! N or K is outside its valid range.
!
   13 IER = -1
      RETURN
!
! No unique solution due to collinear nodes.
!
   14 IER = -2
      RETURN
      END
      SUBROUTINE GRCOEF (SIGMA, D,SD)
      REAL SIGMA, D, SD
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/21/96
!
!   This subroutine computes factors involved in the linear
! systems solved by Subroutines GRADG and SMSGS.
!
! On input:
!
!       SIGMA = Nonnegative tension factor associated with a
!               triangulation arc.
!
! SIGMA is not altered by this routine.
!
! On output:
!
!       D = Diagonal factor.  D = SIG*(SIG*Coshm(SIG) -
!           Sinhm(SIG))/E where E = SIG*Sinh(SIG) - 2*
!           Coshm(SIG).  D > 0, and D = 4 at SIG = 0.
!
!       SD = Off-diagonal factor.  SD = SIG*Sinhm(SIG)/E.
!            SD > 0, and SD = 2 at SIG = 0.
!
! SSRFPACK module required by GRCOEF:  SNHCSH
!
! Intrinsic function called by GRCOEF:  EXP
!
!***********************************************************
!
      REAL COSHM, COSHMM, E, EMS, SCM, SIG, SINHM, SSINH,&
     &     SSM
      SIG = SIGMA
      IF (SIG .LT. 1.E-9) THEN
!
! Cubic function:
!
        D = 4.
        SD = 2.
      ELSEIF (SIG .LE. .5) THEN
!
! 0 < SIG .LE. .5.
!
! Use approximations designed to avoid cancellation error
!   in the hyperbolic functions when SIGMA is small.
!
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        E = SIG*SINHM - COSHMM - COSHMM
        D = SIG*(SIG*COSHM-SINHM)/E
        SD = SIG*SINHM/E
      ELSE
!
! SIG > .5.
!
! Scale SINHM, COSHM, and E by 2*EXP(-SIG) in order to
!   avoid overflow when SIGMA is large.
!
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
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/21/96
!
!   Given a line segment P1-P2 containing a point P, along
! with values and derivatives at the endpoints, this func-
! tion returns the value H(P), where H is the Hermite inter-
! polatory tension spline defined by the endpoint data.
!
! On input:
!
!       B = Local coordinate of P with respect to P1-P2:
!           P = B*P1 + (1-B)*P2, and thus B = d(P,P2)/
!           d(P1,P2), where d(P1,P2) is the distance between
!           P1 and P2.  B < 0 or B > 1 results in extrapola-
!           tion.
!
!       H1,H2 = Values interpolated at P1 and P2, respec-
!               tively.
!
!       HP1,HP2 = Products of d(P1,P2) with first order der-
!                 ivatives at P1 and P2, respectively.  HP1
!                 may, for example, be the scalar product of
!                 P2-P1 with a gradient at P1.
!
!       SIGMA = Nonnegative tension factor associated with
!               the spline.  SIGMA = 0 corresponds to a
!               cubic spline, and H approaches the linear
!               interpolant of H1 and H2 as SIGMA increases.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       HVAL = Interpolated value H(P).
!
! SSRFPACK module required by HVAL:  SNHCSH
!
! Intrinsic functions called by HVAL:  ABS, EXP
!
!***********************************************************
!
      REAL B1, B2, CM, CM2, CMM, D1, D2, DUMMY, E, E1, E2,&
     &     EMS, S, SB1, SB2, SIG, SM, SM2, TM, TM1, TM2, TS
      B1 = B
      B2 = 1. - B1
!
! Compute slope S and second differences D1 and D2 scaled
!   by the separation between P1 and P2.
!
      S = H2 - H1
      D1 = S - HP1
      D2 = HP2 - S
!
! Test the range of SIGMA.
!
      SIG = ABS(SIGMA)
      IF (SIG .LT. 1.E-9) THEN
!
! Hermite cubic interpolation:
!
        HVAL = H1 + B2*(HP1 + B2*(D1 + B1*(D1 - D2)))
      ELSEIF (SIG .LE. .5) THEN
!
! 0 < SIG .LE. .5.  Use approximations designed to avoid
!   cancellation error in the hyperbolic functions.
!
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        E = SIG*SM - CMM - CMM
        HVAL = H1 + B2*HP1 + ((CM*SM2-SM*CM2)*(D1+D2) + SIG*&
     &                     (CM*CM2-(SM+SIG)*SM2)*D1)/(SIG*E)
      ELSE
!
! SIG > .5.  Use negative exponentials in order to avoid
!   overflow.  Note that EMS = EXP(-SIG).
!
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
        HVAL = H1 + B2*S + (TM*TM1*TM2*(D1+D2) + SIG*&
     &                      ((E2*TM1*TM1-B1*TS)*D1 +&
     &                       (E1*TM2*TM2-B2*TS)*D2))/(SIG*E)
      ENDIF
      RETURN
      END
      SUBROUTINE INTRC0 (N,PLAT,PLON,X,Y,Z,W,LIST,LPTR,&
     &                   LEND, IST, PW,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IST, IER
      REAL    PLAT, PLON, X(N), Y(N), Z(N), W(N), PW
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/24/96
!
!   Given a triangulation of a set of nodes on the unit
! sphere, along with data values at the nodes, this sub-
! routine computes the value at a point P of a continuous
! function which interpolates the data values.  The interp-
! olatory function is linear on each underlying triangle
! (planar triangle with the same vertices as a spherical
! triangle).  If P is not contained in a triangle, an ex-
! trapolated value is taken to be the interpolated value at
! the nearest point of the triangulation boundary.
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       PLAT,PLON = Latitude and longitude of P in radians.
!
!       X,Y,Z = Arrays containing Cartesian coordinates of
!               the nodes.
!
!       W = Array containing data values at the nodes.  W(I)
!           is associated with (X(I),Y(I),Z(I)) for I =
!           1,...,N.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
!                        Subroutine TRMESH.
!
!       IST = Index of the starting node in the search for a
!             triangle containing P.  1 .LE. IST .LE. N.
!             The output value of IST from a previous call
!             may be a good choice.
!
! Input parameters other than IST are not altered by this
!   routine.
!
! On output:
!
!       IST = Index of one of the vertices of the triangle
!             containing P (or nearest P) unless IER = -1
!             or IER = -2.
!
!       PW = Value of the interpolatory function at P if
!            IER .GE. 0.
!
!       IER = Error indicator:
!             IER = 0 if interpolation was performed
!                     successfully.
!             IER = 1 if extrapolation was performed
!                     successfully.
!             IER = -1 if N < 3 or IST is outside its valid
!                      range.
!             IER = -2 if the nodes are collinear.
!             IER = -3 if P is not in a triangle and the
!                      angle between P and the nearest boun-
!                      dary point is at least 90 degrees.
!
! STRIPACK modules required by INTRC0:  JRAND, LSTPTR,
!                                         STORE, TRFIND
!
! Intrinsic functions called by INTRC0:  COS, SIN
!
!***********************************************************
!
      INTEGER I1, I2, I3, LP, N1, N2
      REAL    B1, B2, B3, P(3), PTN1, PTN2, S12, SUM
!
! Local parameters:
!
! B1,B2,B3 = Barycentric coordinates of the central projec-
!              tion of P onto the underlying planar trian-
!              gle, or (B1 and B2) projection of Q onto the
!              underlying line segment N1-N2 when P is
!              exterior.  Unnormalized coordinates are
!              computed by TRFIND when P is in a triangle.
! I1,I2,I3 = Vertex indexes returned by TRFIND
! LP =       LIST pointer to N1 as a neighbor of N2 or N2
!              as a neighbor of N1
! N1,N2 =    Endpoints of a boundary arc which is visible
!              from P when P is not contained in a triangle
! P =        Cartesian coordinates of P
! PTN1 =     Scalar product (P,N1)
! PTN2 =     Scalar product (P,N2)
! S12 =      Scalar product (N1,N2)
! SUM =      Quantity used to normalize the barycentric
!              coordinates
!
      IF (N .LT. 3  .OR.  IST .LT. 1  .OR.  IST .GT. N)&
     &    GO TO 11
!
! Transform (PLAT,PLON) to Cartesian coordinates.
!
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
!
! Find the vertex indexes of a triangle containing P.
!
      CALL TRFIND(IST,P,N,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,&
     &            I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 12
      IST = I1
      IF (I3 .NE. 0) THEN
!
! P is contained in the triangle (I1,I2,I3).  Normalize the
!   barycentric coordinates.
!
        SUM = B1 + B2 + B3
        B1 = B1/SUM
        B2 = B2/SUM
        B3 = B3/SUM
        PW = B1*W(I1) + B2*W(I2) + B3*W(I3)
        IER = 0
        RETURN
      ENDIF
!
! P is exterior to the triangulation, and I1 and I2 are
!   boundary nodes which are visible from P.  Set PW to the
!   interpolated value at Q, where Q is the closest boundary
!   point to P.
!
! Traverse the boundary starting from the rightmost visible
!   node I1.
!
      N1 = I1
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
      IF (I1 .NE. I2) GO TO 2
!
! All boundary nodes are visible from P.  Find a boundary
!   arc N1->N2 such that P Left (N2 X N1)->N1.
!
! Counterclockwise boundary traversal:
!   Set N2 to the first neighbor of N1.
!
    1 LP = LEND(N1)
        LP = LPTR(LP)
        N2 = LIST(LP)
!
! Compute inner products (N1,N2) and (P,N2), and compute
!   B2 = DET(P,N1,N2 X N1).
!
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN2 = P(1)*X(N2) + P(2)*Y(N2) + P(3)*Z(N2)
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 2
!
! P Right (N2 X N1)->N1 -- Iterate.
!
        N1 = N2
        I1 = N1
        PTN1 = PTN2
        GO TO 1
!
! P Left (N2 X N1)->N1, where N2 is the first neighbor of P1.
!   Clockwise boundary traversal:
!
    2 N2 = N1
        PTN2 = PTN1
!
! Set N1 to the last neighbor of N2 and test for
!   termination.
!
        LP = LEND(N2)
        N1 = -LIST(LP)
        IF (N1 .EQ. I1) GO TO 13
!
! Compute inner products (N1,N2) and (P,N1).
!
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
!
! Compute B2 = DET(P,N1,N2 X N1) = DET(Q,N1,N2 X N1)*(P,Q).
!
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 2
!
! Compute B1 = DET(P,N2 X N1,N2) = DET(Q,N2 X N1,N2)*(P,Q).
!
      B1 = PTN1 - S12*PTN2
      IF (B1 .LE. 0.) THEN
!
! Q = N2.
!
        PW = W(N2)
      ELSE
!
! P Strictly Left (N2 X N1)->N2 and P Strictly Left
!   N1->(N2 X N1).  Thus Q lies on the interior of N1->N2.
!   Normalize the coordinates and compute PW.
!
        SUM = B1 + B2
        PW = (B1*W(N1) + B2*W(N2))/SUM
      ENDIF
      IER = 1
      RETURN
!
! N or IST is outside its valid range.
!
   11 IER = -1
      RETURN
!
! Collinear nodes.
!
   12 IER = -2
      RETURN
!
! The angular distance between P and the closest boundary
!   point to P is at least 90 degrees.
!
   13 IER = -3
      RETURN
      END
      SUBROUTINE INTRC1 (N,PLAT,PLON,X,Y,Z,F,LIST,LPTR,LEND,&
     &                   IFLGS,SIGMA,IFLGG,GRAD, IST, FP,&
     &                   IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, IFLGG,&
     &        IST, IER
      REAL    PLAT, PLON, X(N), Y(N), Z(N), F(N), SIGMA(*),&
     &        GRAD(3,N), FP
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/25/96
!
!   Given a triangulation of a set of nodes on the unit
! sphere, along with data values and gradients at the nodes,
! this routine computes a value F(P), where F interpolates
! the nodal data and is once-continuously differentiable
! over the convex hull of the nodes.  Refer to Function FVAL
! for further details.  If P is not contained in a triangle,
! an extrapolated value is computed by extending F beyond
! the boundary in a continuous fashion.
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 3
!           and N .GE. 7 if IFLGG .LE. 0.
!
!       PLAT,PLON = Latitude and longitude in radians of the
!                   point P at which F is to be evaluated.
!
!       X,Y,Z = Arrays of length N containing Cartesian
!               coordinates of the nodes.
!
!       F = Array of length N containing values of F at the
!           nodes:  F(I) = F(X(I),Y(I),Z(I)).
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
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
!               programs FVAL, GETSIG, SIG0, SIG1, and SIG2.
!
!       IFLGG = Gradient option:
!               IFLGG .LE. 0 if INTRC1 is to provide grad-
!                            ient estimates as needed (from
!                            GRADL).
!               IFLGG .GE. 1 if gradients are user-provided
!                            in GRAD.  This is more effici-
!                            ent if INTRC1 is to be called
!                            several times.
!
!       GRAD = 3 by N array whose I-th column contains
!              an estimated gradient at node I if IFLGG .GE.
!              1, or unused dummy parameter if IFLGG .LE. 0.
!              Refer to Subroutines GRADL and GRADG.
!
!       IST = Index of the starting node in the search for a
!             triangle containing P.  The output value of
!             IST from a previous call may be a good choice.
!             1 .LE. IST .LE. N.
!
! Input parameters other than IST are not altered by this
!   routine.
!
! On output:
!
!       IST = Index of one of the vertices of the triangle
!             containing P (or a boundary node if P is not
!             contained in a triangle) unless IER = -1 or
!             IER = -2.
!
!       FP = Value of F at P unless IER < 0, in which case
!            FP is not defined.
!
!       IER = Error indicator and information flag:
!             IER = 0 if no errors were encountered and P is
!                     contained in a triangle.
!             IER = 1 if no errors were encountered and
!                     extrapolation was required.
!             IER = -1 if N or IST is outside its valid
!                      range.
!             IER = -2 if the nodes are collinear.
!             IER = -3 if the angular distance between P and
!                      the nearest point of the triangula-
!                      tion is at least 90 degrees.
!
! STRIPACK modules required by INTRC1: JRAND, LSTPTR, STORE,
!                                        TRFIND
!                    (and optionally)  GETNP if IFLGG .LE. 0
!
! SSRFPACK modules required by INTRC1:  ARCINT, ARCLEN,
!                                         FVAL, HVAL, SNHCSH
!              (and if IFLGG .LE. 0)  APLYR, APLYRT, CONSTR,
!                                       GIVENS, GRADL,
!                                       ROTATE, SETUP
!
! Intrinsic functions called by INTRC1:  COS, SIN, SQRT
!
!***********************************************************
!
      INTEGER LSTPTR
      REAL    ARCLEN, FVAL
      INTEGER I, IERR, I1, I2, I3, LP, N1, N2, NN
      REAL    A, B1, B2, B3, DUM(3), FQ, GQ(3), GQN, G1(3),&
     &        G2(3), G3(3), P(3), P1(3), P2(3), P3(3), PTGQ,&
     &        PTN1, PTN2, Q(3), QNORM, S1, S2, S3, S12, SUM
!
! Local parameters:
!
! A =        Angular separation between P and Q
! B1,B2,B3 = Barycentric coordinates of the central projec-
!              tion of P onto the underlying planar triangle,
!              or (B1 and B2) projection of Q onto the
!              underlying line segment N1-N2 when P is
!              exterior.  Unnormalized coordinates are
!              computed by TRFIND when P is in a triangle.
! DUM =      Dummy parameter for ARCINT
! FQ,GQ =    Interpolated value and gradient at Q
! GQN =      Negative of the component of GQ in the direction
!              Q->P
! G1,G2,G3 = Gradients at I1, I2, and I3, or (G1 and G2) at
!              N1 and N2
! I =        DO-loop index
! IERR =     Error flag for calls to GRADL
! I1,I2,I3 = Vertex indexes returned by TRFIND
! LP =       LIST pointer
! N1,N2 =    Indexes of the endpoints of a boundary arc when
!              P is exterior (not contained in a triangle)
! NN =       Local copy of N
! P =        Cartesian coordinates of P
! P1,P2,P3 = Cartesian coordinates of the vertices I1, I2,
!              and I3, or (P1 and P2) coordinates of N1 and
!              N2 if P is exterior
! PTGQ =     Scalar product (P,GQ) -- factor of the component
!              of GQ in the direction Q->P
! PTN1 =     Scalar product (P,N1) -- factor of B1 and B2
! PTN2 =     Scalar product (P,N2) -- factor of B1 and B2
! Q =        Closest boundary point to P when P is exterior
! QNORM =    Factor used to normalize Q
! S1,S2,S3 = Tension factors associated with the triangle
!              sides opposite I1, I2, and I3, or (S1) the
!              boundary arc N1-N2
! S12 =      Scalar product (N1,N2) -- factor of B1 and B2
! SUM =      Quantity used to normalize the barycentric
!              coordinates
!
      NN = N
      IF (NN .LT. 3  .OR.  (IFLGG .LE. 0  .AND.  NN .LT. 7)&
     &    .OR.  IST .LT. 1  .OR.  IST .GT. NN) GO TO 11
!
! Transform (PLAT,PLON) to Cartesian coordinates.
!
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
!
! Locate P with respect to the triangulation.
!
      CALL TRFIND (IST,P,NN,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,&
     &             I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 12
      IST = I1
      IF (I3 .NE. 0) THEN
!
! P is contained in the triangle (I1,I2,I3).  Store the
!   vertex coordinates, gradients, and tension factors in
!   local variables.
!
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
!
!   Gradients are user-provided.
!
          DO 1 I = 1,3
            G1(I) = GRAD(I,I1)
            G2(I) = GRAD(I,I2)
            G3(I) = GRAD(I,I3)
    1       CONTINUE
        ELSE
!
!   Compute gradient estimates at the vertices.
!
          CALL GRADL (NN,I1,X,Y,Z,F,LIST,LPTR,LEND, G1,IERR)
          IF (IERR .LT. 0) GO TO 12
          CALL GRADL (NN,I2,X,Y,Z,F,LIST,LPTR,LEND, G2,IERR)
          IF (IERR .LT. 0) GO TO 12
          CALL GRADL (NN,I3,X,Y,Z,F,LIST,LPTR,LEND, G3,IERR)
          IF (IERR .LT. 0) GO TO 12
        ENDIF
!
        IF (IFLGS .GT. 0) THEN
!
!   Variable tension:
!
          LP = LSTPTR(LEND(I2),I3,LIST,LPTR)
          S1 = SIGMA(LP)
          LP = LSTPTR(LEND(I3),I1,LIST,LPTR)
          S2 = SIGMA(LP)
          LP = LSTPTR(LEND(I1),I2,LIST,LPTR)
          S3 = SIGMA(LP)
        ELSE
!
!   Uniform tension:
!
          S1 = SIGMA(1)
          S2 = S1
          S3 = S1
        ENDIF
!
! Normalize the coordinates.
!
        SUM = B1 + B2 + B3
        B1 = B1/SUM
        B2 = B2/SUM
        B3 = B3/SUM
        FP = FVAL(B1,B2,B3,P1,P2,P3,F(I1),F(I2),F(I3),G1,&
     &            G2,G3,S1,S2,S3)
        IER = 0
        RETURN
      ENDIF
!
! P is exterior to the triangulation, and I1 and I2 are
!   boundary nodes which are visible from P.  Extrapolate to
!   P by linear (with respect to arc-length) interpolation
!   of the value and directional derivative (gradient comp-
!   onent in the direction Q->P) of the interpolatory
!   surface at Q where Q is the closest boundary point to P.
!
! Determine Q by traversing the boundary starting from I1.
!
      N1 = I1
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
      IF (I1 .NE. I2) GO TO 3
!
! All boundary nodes are visible from P.  Find a boundary
!   arc N1->N2 such that P Left (N2 X N1)->N1.
!
! Counterclockwise boundary traversal:
!   Set N2 to the first neighbor of N1.
!
    2 LP = LEND(N1)
        LP = LPTR(LP)
        N2 = LIST(LP)
!
! Compute inner products (N1,N2) and (P,N2), and compute
!   B2 = Det(P,N1,N2 X N1).
!
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN2 = P(1)*X(N2) + P(2)*Y(N2) + P(3)*Z(N2)
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 3
!
! P Right (N2 X N1)->N1:  iterate.
!
        N1 = N2
        I1 = N1
        PTN1 = PTN2
        GO TO 2
!
! P Left (N2 X N1)->N1 where N2 is the first neighbor of N1.
!   Clockwise boundary traversal:
!
    3 N2 = N1
        PTN2 = PTN1
!
! Set N1 to the last neighbor of N2 and test for
!   termination.
!
        LP = LEND(N2)
        N1 = -LIST(LP)
        IF (N1 .EQ. I1) GO TO 13
!
! Compute inner products (N1,N2) and (P,N1).
!
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
!
! Compute B2 = Det(P,N1,N2 X N1) = Det(Q,N1,N2 X N1)*(P,Q).
!
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 3
!
! Compute B1 = Det(P,N2 X N1,N2) = Det(Q,N2 X N1,N2)*(P,Q).
!
      B1 = PTN1 - S12*PTN2
      IF (B1 .LE. 0.) THEN
!
! Q = N2.  Store value, coordinates, and gradient at Q.
!
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
!
! Extrapolate to P:  FP = FQ + A*(GQ,Q X (PXQ)/SIN(A)),
!   where A is the angular separation between Q and P,
!   and Sin(A) is the magnitude of P X Q.
!
        A = ARCLEN(Q,P)
        PTGQ = P(1)*GQ(1) + P(2)*GQ(2) + P(3)*GQ(3)
        FP = FQ
        IF (A .NE. 0.) FP = FP + PTGQ*A/SIN(A)
        IER = 1
        RETURN
      ENDIF
!
! P Strictly Left (N2 X N1)->N2 and P Strictly Left
!   N1->(N2 X N1).  Thus Q lies on the interior of N1->N2.
!   Store coordinates of N1 and N2 in local variables.
!
      P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
!
! Compute the central projection of Q onto P2-P1 and
!   normalize to obtain Q.
!
      QNORM = 0.
      DO 5 I = 1,3
        Q(I) = B1*P1(I) + B2*P2(I)
        QNORM = QNORM + Q(I)*Q(I)
    5   CONTINUE
      QNORM = SQRT(QNORM)
      DO 6 I = 1,3
        Q(I) = Q(I)/QNORM
    6   CONTINUE
!
! Store gradients at N1 and N2 and tension factor S1.
!
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
!
      IF (IFLGS .LE. 0) S1 = SIGMA(1)
      IF (IFLGS .GE. 1) S1 = SIGMA(LP)
!
! Compute an interpolated value and normal gradient
!   component at Q.
!
      CALL ARCINT (Q,P1,P2,F(N1),F(N2),G1,G2,S1, FQ,DUM,GQN)
!
! Extrapolate to P:  the normal gradient component GQN is
!   the negative of the component in the direction Q->P.
!
      FP = FQ - GQN*ARCLEN(Q,P)
      IER = 1
      RETURN
!
! N or IST is outside its valid range.
!
   11 IER = -1
      RETURN
!
! Collinear nodes encountered.
!
   12 IER = -2
      RETURN
!
! The distance between P and the closest boundary point
!   is at least 90 degrees.
!
   13 IER = -3
      RETURN
      END
      SUBROUTINE ROTATE (N,C,S, X,Y )
      INTEGER N
      REAL    C, S, X(N), Y(N)
!
!***********************************************************
!
!                                              From SSRFPACK
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
      SUBROUTINE SETUP (XI,YI,WI,WK,S1,S2,WT, ROW)
      REAL XI, YI, WI, WK, S1, S2, WT, ROW(6)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   05/09/92
!
!   This subroutine sets up the I-th row of an augmented
! regression matrix for a weighted least squares fit of a
! quadratic function Q(X,Y) to a set of data values Wi,
! where Q(0,0) = Wk.  The first 3 columns (quadratic terms)
! are scaled by 1/S2 and the fourth and fifth columns (lin-
! ear terms) are scaled by 1/S1.
!
! On input:
!
!       XI,YI = Coordinates of node I.
!
!       WI = Data value at node I.
!
!       WK = Data value interpolated by Q at the origin.
!
!       S1,S2 = Inverse scale factors.
!
!       WT = Weight factor corresponding to the I-th
!            equation.
!
!       ROW = Array of length 6.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       ROW = Array containing a row of the augmented re-
!             gression matrix.
!
! Modules required by SETUP:  None
!
!***********************************************************
!
      REAL W1, W2
!
! Local parameters:
!
! W1 = Weighted scale factor for the linear terms
! W2 = Weighted scale factor for the quadratic terms
!
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
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/21/98
!
!   Given a triangulation of a set of nodes on the unit
! sphere, along with an array of tension factors associated
! with the triangulation arcs, this subroutine prints the
! list of arcs (with tension factors) ordered by endpoint
! nodal indexes.  An arc is identified with its smaller
! endpoint index:  N1-N2, where N1 < N2.
!
!   This routine is identical to the similarly named routine
! in SRFPACK.
!
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
!                        gulation.  Refer to STRIPACK
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
! STRIPACK module required by SGPRNT:  LSTPTR
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
      DATA NMAX/9999/,  NLMAX/58/
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
      IF (NB .NE. 0) THEN
        NAT = 3*NM1 - NB
      ELSE
        NAT = 3*N - 6
      ENDIF
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
  100 FORMAT ('1',14X,'TENSION FACTORS,  N =',I5,&
     &        ' NODES'//1X,18X,'N1',5X,'N2',8X,'TENSION'//)
  110 FORMAT (1X,16X,I4,3X,I4,5X,F12.8)
  120 FORMAT (1X,16X,I4,3X,I4,5X,F12.8,3X,F12.8,' *')
  130 FORMAT ('1')
  140 FORMAT (//1X,10X,'NA =',I5,' ARCS')
!
! Error messages:
!
  200 FORMAT (//1X,10X,'*',I5,' ERRORS IN SIGMA')
  210 FORMAT (/1X,10X,'*** ERROR IN TRIANGULATION -- ',&
     &        '3N-NB-3 = ',I5,' ***')
  220 FORMAT (1X,10X,'*** N IS OUT OF RANGE -- NMAX = ',&
     &        I4,' ***')
      END
      REAL FUNCTION SIG0 (N1,N2,N,X,Y,Z,H,LIST,LPTR,LEND,&
     &                    GRAD,IFLGB,HBND,TOL,&
     &                    IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGB,&
     &        IFLGS, IER
      REAL    X(N), Y(N), Z(N), H(N), GRAD(3,N), HBND, TOL,&
     &        SIGMA(*)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/21/96
!
!   Given a triangulation of a set of nodes on the unit
! sphere, along with data values H and gradients GRAD at the
! nodes, this function determines the smallest tension fac-
! tor SIG0 such that the Hermite interpolatory tension
! spline H(A), defined by SIG0 and the endpoint values and
! directional derivatives associated with an arc N1-N2, is
! bounded (either above or below) by HBND for all A in
! (A1,A2), where (A1,A2) denotes an interval corresponding
! to the arc and A is the arc-length.
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
!       X,Y,Z = Arrays of length N containing coordinates of
!               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1.
!
!       H = Array of length N containing data values at the
!           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
!           for I = 1 to N.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
!                        Subroutine TRMESH.
!
!       GRAD = Array dimensioned 3 by N whose columns con-
!              tain gradients at the nodes.  GRAD( ,J) must
!              be orthogonal to node J:  GRAD(1,J)*X(J) +
!              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0.  Refer
!              to Subroutines GRADG, GRADL, and SMSURF.
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
!       SIGMA = Dummy parameter (IFLGS .LE. 0) or array con-
!               taining tension factors associated with arcs
!               in one-to-one correspondence with LIST
!               entries (IFLGS .GE. 1).  Refer to Subroutine
!               GETSIG.
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
!                     = H(A1), and the directional deriva-
!                     tive of H at A1 is negative).
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
! STRIPACK module required by SIG0:  STORE
!
! SSRFPACK modules required by SIG0:  ARCLEN, SNHCSH
!
! Intrinsic functions called by SIG0:  ABS, EXP, LOG, MAX,
!                                        MIN, REAL, SIGN,
!                                        SQRT
!
!***********************************************************
!
      REAL    ARCLEN, STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    A, A0, AA, AL, B, B0, BND, C, C1, C2, COSHM,&
     &        COSHMM, D, D0, D1PD2, D2, DMAX, DSIG, E, EMS,&
     &        F, F0, FMAX, FNEG, FTOL, H1, H2, P1(3), P2(3),&
     &        R, RF, RSIG, RTOL, S, S1, S2, SBIG, SCM, SIG,&
     &        SINHM, SNEG, SSINH, SSM, STOL, T, T0, T1, T2,&
     &        TM, UN(3), UNORM
!
      DATA SBIG/85./,  LUN/-1/
      RF = REAL(IFLGB)
      BND = HBND
!
! Print a heading.
!
      IF (LUN .GE. 0  .AND.  RF .LT. 0.) WRITE (LUN,100) N1,&
     &                                   N2, BND
      IF (LUN .GE. 0  .AND.  RF .GT. 0.) WRITE (LUN,110) N1,&
     &                                   N2, BND
  100 FORMAT (//1X,'SIG0 -- N1 =',I4,', N2 =',I4,&
     &        ', LOWER BOUND = ',E15.8)
  110 FORMAT (//1X,'SIG0 -- N1 =',I4,', N2 =',I4,&
     &        ', UPPER BOUND = ',E15.8)
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
! Store nodal coordinates P1 and P2, compute arc-length AL
!   and unit normal UN = (P1 X P2)/UNORM, and test for
!   coincident nodes.
!
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
!
! Store endpoint data values and test for valid constraint.
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
      S1 = AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +&
     &         GRAD(3,N1)*P2(3))/UNORM
      S2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +&
     &          GRAD(3,N2)*P1(3))/UNORM
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
!     of the Hermite cubic interpolant H0(A) = H2 - (S2*R +
!     B0*R**2 + (A0/3)*R**3), where R(A) = (A2-A)/AL.
!
      S = H2 - H1
      T0 = 3.*S - S1 - S2
      A0 = 3.*(S-T0)
      B0 = T0 - S2
      D0 = T0*T0 - S1*S2
!
!   H0 has local extrema in (A1,A2) iff S1*S2 < 0 or
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
!   HP(R) = (S2 - (C1*sinh(SIG*R) - C2*coshm(SIG*R))/E)/DT
!     = 0 for ESR = (-B +/- sqrt(D))/A = C/(-B -/+ sqrt(D))
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
!   H(R) = H2 - (B*SIG*R + C1 + RF*sqrt(D))/(SIG*E).
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
      REAL FUNCTION SIG1 (N1,N2,N,X,Y,Z,H,LIST,LPTR,LEND,GRAD,&
     &                    IFLGB,HPBND,TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGB,&
     &        IFLGS, IER
      REAL    X(N), Y(N), Z(N), H(N), GRAD(3,N), HPBND, TOL,&
     &        SIGMA(*)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/21/96
!
!   Given a triangulation of a set of nodes on the unit
! sphere, along with data values H and gradients GRAD at the
! nodes, this function determines the smallest tension fac-
! tor SIG1 such that the first derivative HP(A) of the
! Hermite interpolatory tension spline H(A), defined by SIG1
! and the endpoint values and directional derivatives asso-
! ciated with an arc N1-N2, is bounded (either above or
! below) by HPBND for all A in (A1,A2), where (A1,A2) de-
! notes an interval corresponding to the arc and A denotes
! arc-length.
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
!       X,Y,Z = Arrays of length N containing coordinates of
!               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1.
!
!       H = Array of length N containing data values at the
!           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
!           for I = 1 to N.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
!                        Subroutine TRMESH.
!
!       GRAD = Array dimensioned 3 by N whose columns con-
!              gradients at the nodes.  GRAD( ,J) must be
!              orthogonal to node J:  GRAD(1,J)*X(J) +
!              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0.  Refer
!              to Subroutines GRADG, GRADL, and SMSURF.
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
!       SIGMA = Dummy parameter (IFLGS .LE. 0) or array
!               containing tension factors associated with
!               arcs in one-to-one correspondence with LIST
!               entries (IFLGS .GE. 1).  Refer to Subroutine
!               GETSIG.
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
! STRIPACK module required by SIG1:  STORE
!
! SSRFPACK modules required by SIG1:  ARCLEN, SNHCSH
!
! Intrinsic functions called by SIG1:   ABS, EXP, MAX, MIN,
!                                         REAL, SIGN, SQRT
!
!***********************************************************
!
      REAL    ARCLEN, STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    A, A0, AL, B0, BND, C0, C1, C2, COSHM, COSHMM,&
     &        D0, D1, D1PD2, D2, DMAX, DSIG, E, EMS, EMS2,&
     &        F, F0, FMAX, FNEG, FTOL, P1(3), P2(3), RF,&
     &        RTOL, S, S1, S2, SBIG, SIG, SINH, SINHM, STOL,&
     &        T0, T1, T2, TM, UN(3), UNORM
!
      DATA SBIG/85./,  LUN/-1/
      RF = REAL(IFLGB)
      BND = HPBND
!
! Print a heading.
!
      IF (LUN .GE. 0  .AND.  RF .LT. 0.) WRITE (LUN,100) N1,&
     &                                   N2, BND
      IF (LUN .GE. 0  .AND.  RF .GT. 0.) WRITE (LUN,110) N1,&
     &                                   N2, BND
  100 FORMAT (//1X,'SIG1 -- N1 =',I4,', N2 =',I4,&
     &        ', LOWER BOUND = ',E15.8)
  110 FORMAT (//1X,'SIG1 -- N1 =',I4,', N2 =',I4,&
     &        ', UPPER BOUND = ',E15.8)
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
! Store nodal coordinates P1 and P2, compute arc-length AL
!   and unit normal UN = (P1 X P2)/UNORM, and test for
!   coincident nodes.
!
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
!
! Compute first difference S and scaled directional deriva-
!   tives S1,S2 at the endpoints (for the direction N1->N2).
!
      S = H(N2) - H(N1)
      S1 = AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +&
     &         GRAD(3,N1)*P2(3))/UNORM
      S2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +&
     &          GRAD(3,N2)*P1(3))/UNORM
!
! Test for a valid constraint.
!
      IER = -3
      IF ((RF .LT. 0.  .AND.  MIN(S1,S2,S) .LT. BND)  .OR.&
     &    (RF .GT. 0.  .AND.  BND .LT. MAX(S1,S2,S)))&
     &   GO TO 11
!
! Test for infinite tension required.
!
      IER = 1
      SIG = SBIG
      IF (S .EQ. BND  .AND.  (S1 .NE. S  .OR.  S2 .NE. S))&
     &   GO TO 10
!
! Test for SIG = 0 sufficient.  The Hermite cubic interpo-
!   lant H0 has derivative HP0(T) = (S2 + 2*B0*R + A0*R**2)/
!   AL, where R = (T2-T)/AL.
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
      IF (B0*C0 .LE. 0.  .OR.  A0*RF .GT. 0.) GO TO 10
!
!   A0*RF < 0 and HP0(R) = -D0/(DT*A0) at R = -B0/A0.
!
      D0 = T0*T0 - S1*S2
      F0 = (BND + D0/(A0*AL))*RF
      IF (F0 .GE. 0.) GO TO 10
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
!   (AL*(SIG-2.)))*RF -- a value for which F(SIG) .GE. 0 and
!   F(SIG) = 0 for SIG sufficiently large that 2*SIG is in-
!   significant relative to exp(SIG).
!
      FMAX = (BND-S/AL)*RF
      SIG = 2. - A0/(3.*(AL*BND-S))
      IF (LUN .GE. 0) WRITE (LUN,120) F0, FMAX, SIG
  120 FORMAT (1X,9X,'F0 = ',E15.8,', FMAX = ',E15.8/&
     &        1X,8X,'SIG = ',E15.8/)
      IF (STORE(SIG*EXP(-SIG)+.5) .EQ. .5) GO TO 10
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
!   HP(R) = (B+SIGN(A)*SQRT(A*C))/(AL*E) at the critical
!     value of R, where A = C2-C1, B = E*S2-C2, and C = C2 +
!     C1.  Note that RF*A < 0.
!
          F = (BND - (E*S2-C2 - RF*SQRT(T1))/(AL*E))*RF
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
     &    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 10
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
   10 SIG1 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
!
! Error termination.
!
   11 SIG1 = -1.
      RETURN
      END
      REAL FUNCTION SIG2 (N1,N2,N,X,Y,Z,H,LIST,LPTR,LEND,&
     &                    GRAD,TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGS,&
     &        IER
      REAL    X(N), Y(N), Z(N), H(N), GRAD(3,N), TOL,&
     &        SIGMA(*)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/25/96
!
!   Given a triangulation of a set of nodes on the unit
! sphere, along with data values H and gradients GRAD at the
! nodes, this function determines the smallest tension fac-
! tor SIG2 such that the Hermite interpolatory tension
! spline H(A), defined by SIG2 and the endpoint values and
! directional derivatives associated with an arc N1-N2,
! preserves convexity (or concavity) of the data:
!
!   HP1 .LE. S .LE. HP2 implies HPP(A) .GE. 0, and
!   HP1 .GE. S .GE. HP2 implies HPP(A) .LE. 0
!
! for all A in the open interval (A1,A2) corresponding to
! the arc, where HP1 and HP2 are the derivative values of H
! at the endpoints, S is the slope of the linear interpolant
! of the endpoint data values, HPP denotes the second deriv-
! ative of H, and A is arc-length.  Note, however, that
! infinite tension is required if HP1 = S or HP2 = S (unless
! HP1 = HP2 = S).
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
!       X,Y,Z = Arrays of length N containing coordinates of
!               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1.
!
!       H = Array of length N containing data values at the
!           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
!           for I = 1 to N.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
!                        Subroutine TRMESH.
!
!       GRAD = Array dimensioned 3 by N whose columns con-
!              gradients at the nodes.  GRAD( ,J) must be
!              orthogonal to node J:  GRAD(1,J)*X(J) +
!              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0.  Refer
!              to Subroutines GRADG, GRADL, and SMSURF.
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
!       SIGMA = Dummy parameter (IFLGS .LE. 0) or array
!               containing tension factors associated with
!               arcs in one-to-one correspondence with LIST
!               entries (IFLGS .GE. 1).  Refer to Subroutine
!               GETSIG.
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
! STRIPACK module required by SIG2:  STORE
!
! SSRFPACK modules required by SIG2:  ARCLEN, SNHCSH
!
! Intrinsic functions called by SIG2:  ABS, EXP, MAX, MIN,
!                                        SQRT
!
!***********************************************************
!
      REAL    ARCLEN, STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      REAL    AL, COSHM, D1, D1D2, D2, DSIG, DUMMY, EMS, F,&
     &        FP, FTOL, P1(3), P2(3), RTOL, S, SBIG, SIG,&
     &        SINHM, SSM, T, T1, TP1, UN(3), UNORM
!
      DATA SBIG/85./,  LUN/-1/
!
! Print a heading.
!
      IF (LUN .GE. 0) WRITE (LUN,100) N1, N2
  100 FORMAT (//1X,'SIG2 -- N1 =',I4,', N2 =',I4)
!
! Test for errors and set local parameters.
!
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.&
     &    MAX(N1,N2,3) .GT. N) GO TO 11
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
! Store nodal coordinates P1 and P2, compute arc-length AL
!   and unit normal UN = (P1 X P2)/UNORM, and test for
!   coincident nodes.
!
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
!
! Compute first and second differences and test for infinite
!   tension required.
!
      S = H(N2) - H(N1)
      D1 = S - AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +&
     &             GRAD(3,N1)*P2(3))/UNORM
      D2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +&
     &          GRAD(3,N2)*P1(3))/UNORM - S
      D1D2 = D1*D2
      IER = 1
      SIG = SBIG
      IF (D1D2 .EQ. 0.  .AND.  D1 .NE. D2) GO TO 10
!
! Test for a valid constraint.
!
      IER = 2
      SIG = 0.
      IF (D1D2 .LT. 0.) GO TO 10
!
! Test for SIG = 0 sufficient.
!
      IER = 0
      IF (D1D2 .EQ. 0.) GO TO 10
      T = MAX(D1/D2,D2/D1)
      IF (T .LE. 2.) GO TO 10
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
      IF (FP .LE. 0.) GO TO 10
      DSIG = -F/FP
      IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.  .AND.&
     &    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 10
!
!   Bottom of loop:  update SIG.
!
      SIG = SIG + DSIG
      GO TO 6
!
! No errors encountered.
!
   10 SIG2 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
!
! Error termination.
!
   11 SIG2 = -1.
      RETURN
      END
      SUBROUTINE SMSGS (N,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,&
     &                  SIGMA,W,P, NIT,DFMAX,F,GRAD, IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, NIT, IER
      REAL    X(N), Y(N), Z(N), U(N), SIGMA(*), W(N), P,&
     &        DFMAX, F(N), GRAD(3,N)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/25/96
!
!   This subroutine solves the symmetric positive definite
! linear system associated with minimizing the quadratic
! functional Q(F,FX,FY,FZ) described in Subroutine SMSURF.
! Since the gradient at node K lies in the plane tangent to
! the sphere surface at K, it is effectively defined by only
! two components -- its X and Y components in the coordinate
! system obtained by rotating K to the north pole.  Thus,
! the minimization problem corresponds to an order-3N system
! which is solved by the block Gauss-Seidel method with 3 by
! 3 blocks.
!
! On input:
!
!       N,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,SIGMA,W = Parameters
!           as described in Subroutine SMSURF.
!
!       P = Positive smoothing parameter defining Q.
!
! The above parameters are not altered by this routine.
!
!       NIT = Maximum number of iterations to be used.  This
!             maximum will likely be achieved if DFMAX is
!             smaller than the machine precision.  NIT .GE.
!             0.
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
!       GRAD = 3 by N array containing initial estimates of
!              the last 3N solution components (the gradi-
!              ent with FX, FY, and FZ in rows 1, 2, and 3,
!              respectively).
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
!       GRAD = Last 3N solution components -- gradients at
!              the nodes.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and the
!                     convergence criterion was achieved.
!             IER = 1 if no errors were encountered but con-
!                     vergence was not achieved within NIT
!                     iterations.
!             IER = -1 if N, P, NIT, or DFMAX is outside its
!                      valid range on input.  F and GRAD are
!                      not altered in this case.
!             IER = -2 if all nodes are collinear or the
!                      triangulation is invalid.
!             IER = -3 if duplicate nodes were encountered.
!
! SSRFPACK modules required by SMSGS:  APLYRT, CONSTR,
!                                        GRCOEF, SNHCSH
!
! Intrinsic functions called by SMSGS:  ABS, ATAN, MAX, SQRT
!
!***********************************************************
!
      INTEGER IFL, ITER, ITMAX, J, K, LPJ, LPL, NN
      REAL    ALFA, ALFSQ, C11, C12, C13, C22, C23, C33,&
     &        CC22, CC23, CC33, CX, CY, DEN1, DEN2, DET, DF,&
     &        DFMX, DGK(3), DGX, DGY, FK, G1, G2, G3, GJK,&
     &        GKJ, PP, R1, R2, R3, RR2, RR3, SIG, SINAL, SX,&
     &        SY, T, T1, T2, T3, T4, T5, T6, TOL, XJ, XK,&
     &        XS, YJ, YK, YS, ZJ, ZK
!
      NN = N
      IFL = IFLGS
      PP = P
      ITMAX = NIT
      TOL = DFMAX
!
! Test for errors in input and initialize iteration count,
!   tension factor, and output value of DFMAX.
!
      IF (NN .LT. 3  .OR.  PP .LE. 0.  .OR.  ITMAX .LT. 0&
         &.OR.  TOL .LT. 0.) GO TO 5
      ITER = 0
      SIG = SIGMA(1)
      DFMX = 0.
!
! Top of iteration loop.
!
    1 IF (ITER .EQ. ITMAX) GO TO 4
      DFMX = 0.
!
!   Loop on nodes.
!
      DO 3 K = 1,NN
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        FK = F(K)
        G1 = GRAD(1,K)
        G2 = GRAD(2,K)
        G3 = GRAD(3,K)
!
!   Construct the rotation mapping node K to the north pole.
!
        CALL CONSTR (XK,YK,ZK, CX,SX,CY,SY)
!
!   Initialize components of the order-3 system for the
!     change (DF,DGX,DGY) in the K-th solution components.
!
        C11 = PP*W(K)
        C12 = 0.
        C13 = 0.
        C22 = 0.
        C23 = 0.
        C33 = 0.
        R1 = C11*(U(K)-FK)
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
!   Compute the coordinates of J in the rotated system.
!
          T = SX*Y(J) + CX*Z(J)
          YJ = CX*Y(J) - SX*Z(J)
          ZJ = SY*X(J) + CY*T
          XJ = CY*X(J) - SY*T
!
!   Compute arc-length ALFA between K and J, ALFSQ = ALFA*
!     ALFA, SINAL = SIN(ALFA), DEN1 = ALFA*SIN(ALFA)**2, and
!     DEN2 = ALFSQ*SINAL.
!
          ALFA = 2.*ATAN(SQRT((1.-ZJ)/(1.+ZJ)))
          ALFSQ = ALFA*ALFA
          XS = XJ*XJ
          YS = YJ*YJ
          SINAL = SQRT(XS+YS)
          DEN1 = ALFA*(XS+YS)
          DEN2 = ALFSQ*SINAL
!
!   Test for coincident nodes and compute functions of SIG:
!     T1 = SIG*SIG*COSHM/E, T2 = SIG*SINHM/E, and T3 = SIG*
!     (SIG*COSHM-SINHM)/E for E = SIG*SINH - 2*COSHM.
!
          IF (DEN1 .EQ. 0.) GO TO 7
          IF (IFL .GE. 1) SIG = SIGMA(LPJ)
          CALL GRCOEF (SIG, T3,T2)
          T1 = T2 + T3
!
!   Update system components for node J.
!
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
!
!   Bottom of loop on neighbors.
!
          IF (LPJ .NE. LPL) GO TO 2
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
     &      C11 .EQ. 0.) GO TO 6
        DGY = (CC22*RR3 - CC23*RR2)/DET
        DGX = (RR2 - CC23*DGY)/CC22
        DF = (R1 - C12*DGX - C13*DGY)/C11
!
!   Rotate (DGX,DGY,0) back to the original coordinate
!     system, and update GRAD( ,K), F(K), and DFMX.
!
        CALL APLYRT (DGX,DGY,CX,SX,CY,SY, DGK)
        GRAD(1,K) = G1 + DGK(1)
        GRAD(2,K) = G2 + DGK(2)
        GRAD(3,K) = G3 + DGK(3)
        F(K) = FK + DF
        DFMX = MAX(DFMX,ABS(DF)/(1.+ABS(FK)))
    3   CONTINUE
!
!   Increment ITER and test for convergence.
!
      ITER = ITER + 1
      IF (DFMX .GT. TOL) GO TO 1
!
! The method converged.
!
      NIT = ITER
      DFMAX = DFMX
      IER = 0
      RETURN
!
! The method failed to converge within NIT iterations.
!
    4 DFMAX = DFMX
      IER = 1
      RETURN
!
! Invalid input parameter.
!
    5 NIT = 0
      DFMAX = 0.
      IER = -1
      RETURN
!
! Node K and its neighbors are collinear.
!
    6 NIT = 0
      DFMAX = DFMX
      IER = -2
      RETURN
!
! Nodes J and K coincide.
!
    7 NIT = 0
      DFMAX = DFMX
      IER = -3
      RETURN
      END





      SUBROUTINE SMSURF (N,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,&
     &                   SIGMA,W,SM,SMTOL,GSTOL,LPRNT, F,&
     &                   GRAD,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, LPRNT,&
     &        IER
      REAL    X(N), Y(N), Z(N), U(N), SIGMA(*), W(N), SM,&
     &        SMTOL, GSTOL, F(N), GRAD(3,N)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/21/98
!
!   Given a triangulation of N nodes on the unit sphere with
! data values U at the nodes and tension factors SIGMA
! associated with the arcs, this routine determines a set of
! nodal function values F and gradients GRAD = (FX,FY,FZ)
! such that a quadratic functional Q1(F,GRAD) is minimized
! subject to the constraint Q2(F) .LE. SM for Q2(F) =
! (U-F)**T*W*(U-F), where W is a diagonal matrix of positive
! weights.  The functional Q1 is an approximation to the
! linearized curvature over the triangulation of a C-1 fun-
! ction F(V), V a unit vector, which interpolates the nodal
! values and gradients.  Subroutines INTRC1 and UNIF may be
! called to evaluate F at arbitrary points.
!
!   The smoothing procedure is an extension of the method
! for cubic spline smoothing due to C. Reinsch -- Numer.
! Math., 10 (1967) and 16 (1971).  Refer to Function FVAL
! for a further description of the interpolant F.  Letting
! D1F(T) and D2F(T) denote first and second derivatives of F
! with respect to a parameter T varying along a triangula-
! tion arc, Q1 is the sum of integrals over the arcs of
! D2F(T)**2 + ((SIGMA/L)*(D1F(T)-S))**2 where L denotes arc-
! length, SIGMA is the appropriate tension factor, and S is
! the slope of the linear function of T which interpolates
! the values of F at the endpoints of the arc.  Introducing
! a smoothing parameter P, and assuming the constraint is
! active, the problem is equivalent to minimizing Q(P,F,
! GRAD) = Q1(F,GRAD) + P*(Q2(F)-SM).  The secant method is
! used to find a zero of G(P) = 1/SQRT(Q2) - 1/SQRT(SM)
! where F(P) satisfies the order-3N symmetric positive def-
! inite linear system obtained by setting the gradient of Q
! (treated as a function of F and GRAD with GRAD tangent to
! the sphere surface) to zero.  The linear system is solved
! by the block Gauss-Seidel method (refer to SMSGS).
!
!   Note that the method can also be used to select grad-
! ients for the interpolation problem (F = U, SM = 0, and P
! infinite).  This is achieved by a call to Subroutine
! GRADG.
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y,Z = Arrays of length N containing Cartesian
!               coordinates of the nodes.
!
!       U = Array of length N containing data values at the
!           nodes.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
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
!               programs GETSIG, SIG0, SIG1, and SIG2.
!
!       W = Array of length N containing positive weights
!           associated with the data values.  The recommend-
!           ed value of W(I) is 1/DU**2 where DU is the
!           standard deviation associated with U(I).  DU**2
!           is the expected value of the squared error in
!           the measurement of U(I).  (The mean error is
!           assumed to be zero.)
!
!       SM = Positive parameter specifying an upper bound on
!            Q2(F).  Note that F is constant (and Q2(F)
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
!       LPRNT = Logical unit on which diagnostic messages
!               are printed, or negative integer specifying
!               no diagnostics.  For each secant iteration,
!               the following values are printed:  P, G(P),
!               NIT, DFMAX, and DP, where NIT denotes the
!               number of Gauss-Seidel iterations used in
!               the computation of G, DFMAX denotes the max-
!               imum relative change in a solution component
!               in the last Gauss-Seidel iteration, and DP
!               is the change in P computed by linear inter-
!               polation between the current point (P,G) and
!               a previous point.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       F = Array of length N containing nodal function val-
!           ues unless IER < 0.
!
!       GRAD = 3 by N array whose columns contain gradients
!              of F at the nodes unless IER < 0.
!
!       IER = Error indicator and information flag:
!             IER = 0 if no errors were encountered and the
!                     constraint is active -- Q2(F) is ap-
!                     proximately equal to SM.
!             IER = 1 if no errors were encountered but the
!                     constraint is not active -- F and GRAD
!                     are the values and gradients of a con-
!                     stant function which minimizes Q2(F),
!                     and Q1 = 0.
!             IER = 2 if the constraint could not be satis-
!                     fied to within SMTOL due to
!                     ill-conditioned linear systems.
!             IER = -1 if N, W, SM, SMTOL, or GSTOL is out-
!                      side its valid range on input.
!             IER = -2 if all nodes are collinear or the
!                      triangulation is invalid.
!             IER = -3 if duplicate nodes were encountered.
!
! SSRFPACK modules required by SMSURF:  APLYRT, CONSTR,
!                                         GRCOEF, SMSGS,
!                                         SNHCSH
!
! Intrinsic functions called by SMSURF:  ABS, SQRT
!
!***********************************************************
!
      INTEGER I, IERR, ITER, ITMAX, LUN, NIT, NITMAX, NN
      REAL    C, DFMAX, DMAX, DP, G, G0, GNEG, P, Q2, Q2MAX,&
     &        Q2MIN, S, SUMW, TOL, WI
!
! Local parameters:
!
! ITMAX = Maximum number of secant iterations.
! LUN = Local copy of LPRNT.
! NITMAX = Maximum number of Gauss-Seidel iterations for
!          each secant iteration.
! NN = Local copy of N.
! TOL = Local copy of GSTOL.
!
      DATA ITMAX/50/,  NITMAX/40/
!
      NN = N
      TOL = GSTOL
      LUN = LPRNT
      IF (LUN .GT. 99) LUN = -1
!
! Test for errors and initialize F to the weighted least
!   squares fit of a constant function to the data.
!
      IER = -1
      IF (NN .LT. 3  .OR.  SM .LE. 0.  .OR.  SMTOL .LE. 0.&
         &.OR.  SMTOL .GE. 1.  .OR.  TOL .LE. 0.) RETURN
      C = 0.
      SUMW = 0.
      DO 1 I = 1,NN
        WI = W(I)
        IF (WI .LE. 0.) RETURN
        C = C + WI*U(I)
        SUMW = SUMW + WI
    1   CONTINUE
      C = C/SUMW
!
! Compute nodal values and gradients, and accumulate Q2 =
!   (U-F)**T*W*(U-F).
!
      Q2 = 0.
      DO 2 I = 1,NN
        F(I) = C
        GRAD(1,I) = 0.
        GRAD(2,I) = 0.
        GRAD(3,I) = 0.
        Q2 = Q2 + W(I)*(U(I)-F(I))**2
    2   CONTINUE
!
! Compute bounds on Q2 defined by SMTOL, and test for the
!   constraint satisfied by the constant fit.
!
      Q2MIN = SM*(1.-SMTOL)
      Q2MAX = SM*(1.+SMTOL)
      IF (Q2 .LE. Q2MAX) THEN
!
! The constraint is satisfied by a constant function.
!
        IER = 1
        IF (LUN .GE. 0) WRITE (LUN,100)
  100   FORMAT (///1X,'SMSURF -- THE CONSTRAINT IS NOT ',&
     &          'ACTIVE AND THE FITTING FCN IS CONSTANT.')
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
! G(P) is strictly increasing and concave, and G(0) .LT. 0.
!   Initialize parameters for the secant method.  The method
!   uses three points -- (P0,G0), (P,G), and (PNEG,GNEG)
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
    3 NIT = NITMAX
      DFMAX = TOL
      CALL SMSGS (NN,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,SIGMA,W,&
     &            P, NIT,DFMAX,F,GRAD, IERR)
      IF (IERR .LT. 0) IER = IERR
!
!   IERR = -1 in SMSGS could be caused by P = 0 as a result
!     of inaccurate solutions to ill-conditioned systems.
!
      IF (IERR .EQ. -1) IER = 2
      IF (IERR .LT. 0) RETURN
      Q2 = 0.
      DO 4 I = 1,NN
        Q2 = Q2 + W(I)*(U(I)-F(I))**2
    4   CONTINUE
      G = 1./SQRT(Q2) - S
      ITER = ITER + 1
      IF (LUN .GE. 0) WRITE (LUN,120) ITER, P, G, NIT, DFMAX
  120 FORMAT (/1X,I2,' -- P = ',E15.8,', G = ',E15.8,&
     &        ', NIT = ',I2,', DFMAX = ',E12.6)
!
!   Test for convergence.
!
      IF (Q2MIN .LE. Q2  .AND.  Q2 .LE. Q2MAX) RETURN
      IF (ITER .GE. ITMAX) THEN
        IER = 2
        RETURN
      ENDIF
      IF (DMAX .EQ. 0.  .AND.  G .LE. 0.) THEN
!
!   Increase P until G(P) > 0.
!
        P = 10.*P
        DP = P
        GO TO 3
      ENDIF
!
!   A bracketing interval [P0,P] has been found.
!
      IF (G0*G .LE. 0.) THEN
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
    5 DP = -G*DP/(G-G0)
      IF (LUN .GE. 0) WRITE (LUN,130) DP
  130 FORMAT (1X,5X,'DP = ',E15.8)
      IF (ABS(DP) .GT. ABS(DMAX)) THEN
!
!   G0*G .GT. 0 and the new estimate would be outside of the
!     bracketing interval of length ABS(DMAX).  Reset
!     (P0,G0) to (PNEG,GNEG).
!
        DP = DMAX
        G0 = GNEG
        GO TO 5
      ENDIF
!
!   Bottom of loop -- update P, DMAX, and G0.
!
      P = P + DP
      DMAX = DMAX + DP
      G0 = G
      GO TO 3
      END
      SUBROUTINE SNHCSH (X, SINHM,COSHM,COSHMM)
      REAL X, SINHM, COSHM, COSHMM
!
!***********************************************************
!
!                                              From SSRFPACK
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
      SUBROUTINE UNIF (N,X,Y,Z,F,LIST,LPTR,LEND,IFLGS,SIGMA,&
     &                 NROW,NI,NJ,PLAT,PLON,IFLGG, GRAD, FF,&
     &                 IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, NROW, NI,&
     &        NJ, IFLGG, IER
      REAL    X(N), Y(N), Z(N), F(N), SIGMA(*), PLAT(NI),&
     &        PLON(NJ), GRAD(3,N), FF(NROW,NJ)
!
!***********************************************************
!
!                                              From SSRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/25/96
!
!   Given a Delaunay triangulation of a set of nodes on the
! unit sphere, along with data values and tension factors
! associated with the triangulation arcs, this routine
! interpolates the data values to a uniform grid for such
! applications as contouring.  The interpolant is once con-
! tinuously differentiable.  Extrapolation is performed at
! grid points exterior to the triangulation when the nodes
! do not cover the entire sphere.
!
! On input:
!
!       N = Number of nodes.  N .GE. 3 and N .GE. 7 if
!           IFLAG .NE. 1.
!
!       X,Y,Z = Arrays containing Cartesian coordinates of
!               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1
!               for I = 1 to N.
!
!       F = Array containing data values.  F(I) is associ-
!           ated with (X(I),Y(I),Z(I)).
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
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
!               programs GETSIG, SIG0, SIG1, and SIG2.
!
!       NROW = Number of rows in the dimension statement of
!              FF.
!
!       NI,NJ = Number of rows and columns in the uniform
!               grid.  1 .LE. NI .LE. NROW and 1 .LE. NJ.
!
!       PLAT,PLON = Arrays of length NI and NJ, respective-
!                   ly, containing the latitudes and
!                   longitudes of the grid lines.
!
!       IFLGG = Option indicator:
!               IFLGG = 0 if gradient estimates at the ver-
!                         tices of a triangle are to be
!                         recomputed for each grid point in
!                         the triangle and not saved.
!               IFLGG = 1 if gradient estimates are input in
!                         GRAD.
!               IFLGG = 2 if gradient estimates are to be
!                         computed once for each node (by
!                         GRADL) and saved in GRAD.
!
! The above parameters are not altered by this routine.
!
!       GRAD = 3 by N array whose columns contain the X, Y,
!              and Z components (in that order) of the grad-
!              ients at the nodes if IFLGG = 1, array of
!              sufficient size if IFLGG = 2, or dummy para-
!              meter if IFLGG = 0.
!
! Gradient estimates may be computed by Subroutines GRADL or
!   GRADG if IFLGG = 1.
!
!       FF = NROW by NCOL array with NROW .GE. NI and NCOL
!            .GE. NJ.
!
! On output:
!
!       GRAD = Array containing estimated gradients as de-
!              fined above if IFLGG = 2 and IER .GE. 0.
!              GRAD is not altered if IFLGG .NE. 2.
!
!       FF = Interpolated values at the grid points if IER
!            .GE. 0.  FF(I,J) = F(PLAT(I),PLON(J)) for I =
!            1,...,NI and J = 1,...,NJ.
!
!       IER = Error indicator:
!             IER = K if no errors were encountered and K
!                     grid points required extrapolation for
!                     K .GE. 0.
!             IER = -1 if N, NI, NJ, or IFLGG is outside its
!                      valid range.
!             IER = -2 if the nodes are collinear.
!             IER = -3 if extrapolation failed due to the
!                      uniform grid extending too far beyond
!                      the triangulation boundary.
!
! STRIPACK modules required by UNIF:  GETNP, JRAND, LSTPTR,
!                                       STORE, TRFIND
!
! SSRFPACK modules required by UNIF:  APLYR, APLYRT, ARCINT,
!                                       ARCLEN, CONSTR,
!                                       FVAL, GIVENS, GRADL,
!                                       HVAL, INTRC1,
!                                       ROTATE, SETUP,
!                                       SNHCSH
!
!***********************************************************
!
      INTEGER I, J, IERR, IFL, IST, NEX, NN, NST, NX, NY
      DATA    NST/1/
!
! Local parameters:
!
! I,J =   DO-loop indexes
! IERR =  Error flag for calls to GRADL and INTRC1
! IFL =   Local copy of IFLGG
! IST =   Parameter for INTRC1
! NEX =   Number of grid points exterior to the triangula-
!           tion boundary (number of extrapolated values)
! NN =    Local copy of N
! NST =   Initial value for IST
! NX,NY = Local copies of NI and NJ
!
      NN = N
      NX = NI
      NY = NJ
      IFL = IFLGG
      IF (NX .LT. 1  .OR.  NX .GT. NROW  .OR.  NY .LT. 1&
        &.OR.  IFL .LT. 0  .OR.  IFL .GT. 2) GO TO 4
      IST = NST
      IF (IFL .EQ. 2) THEN
!
! Compute gradient estimates at the nodes.
!
        DO 1 I = 1,NN
          CALL GRADL (NN,I,X,Y,Z,F,LIST,LPTR,&
     &                LEND, GRAD(1,I),IERR)
          IF (IERR .LT. 0) GO TO 5
    1     CONTINUE
        IFL = 1
      ENDIF
!
! Compute uniform grid points and interpolated values.
!
      NEX = 0
      DO 3 J = 1,NY
        DO 2 I = 1,NX
          CALL INTRC1 (NN,PLAT(I),PLON(J),X,Y,Z,F,LIST,LPTR,&
     &                 LEND,IFLGS,SIGMA,IFL,&
     &                 GRAD, IST, FF(I,J),IERR)
          IF (IERR .LT. 0) GO TO 5
          NEX = NEX + IERR
    2     CONTINUE
    3   CONTINUE
      IER = NEX
      RETURN
!
! NI, NJ, or IFLGG is outside its valid range.
!
    4 IER = -1
      RETURN
!
! Error in GRADL or INTRC1.
!
    5 IER = IERR
      RETURN
      END
