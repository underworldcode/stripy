      MODULE SWPCOM
      REAL SWTOL
      SAVE SWTOL
      END

      SUBROUTINE ADDCST (NCC,LCC,N,X,Y, LWK,IWK,LIST,LPTR,
     .                   LEND, IER)
      INTEGER NCC, LCC(*), N, LWK, IWK(LWK), LIST(*),
     .        LPTR(*), LEND(N), IER
      REAL    X(N), Y(N)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/12/94
C
C   This subroutine provides for creation of a constrained
C Delaunay triangulation which, in some sense, covers an
C arbitrary connected region R rather than the convex hull
C of the nodes.  This is achieved simply by forcing the
C presence of certain adjacencies (triangulation arcs) cor-
C responding to constraint curves.  The union of triangles
C coincides with the convex hull of the nodes, but triangles
C in R can be distinguished from those outside of R.  The
C only modification required to generalize the definition of
C the Delaunay triangulation is replacement of property 5
C (refer to TRMESH) by the following:
C
C  5')  If a node is contained in the interior of the cir-
C       cumcircle of a triangle, then every interior point
C       of the triangle is separated from the node by a
C       constraint arc.
C
C   In order to be explicit, we make the following defini-
C tions.  A constraint region is the open interior of a
C simple closed positively oriented polygonal curve defined
C by an ordered sequence of three or more distinct nodes
C (constraint nodes) P(1),P(2),...,P(K), such that P(I) is
C adjacent to P(I+1) for I = 1,...,K with P(K+1) = P(1).
C Thus, the constraint region is on the left (and may have
C nonfinite area) as the sequence of constraint nodes is
C traversed in the specified order.  The constraint regions
C must not contain nodes and must not overlap.  The region
C R is the convex hull of the nodes with constraint regions
C excluded.
C
C   Note that the terms boundary node and boundary arc are
C reserved for nodes and arcs on the boundary of the convex
C hull of the nodes.
C
C   The algorithm is as follows:  given a triangulation
C which includes one or more sets of constraint nodes, the
C corresponding adjacencies (constraint arcs) are forced to
C be present (Subroutine EDGE).  Any additional new arcs
C required are chosen to be locally optimal (satisfy the
C modified circumcircle property).
C
C
C On input:
C
C       NCC = Number of constraint curves (constraint re-
C             gions).  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy array of length
C             1 if NCC = 0) containing the index (for X, Y,
C             and LEND) of the first node of constraint I in
C             LCC(I) for I = 1 to NCC.  Thus, constraint I
C             contains K = LCC(I+1) - LCC(I) nodes, K .GE.
C             3, stored in (X,Y) locations LCC(I), ...,
C             LCC(I+1)-1, where LCC(NCC+1) = N+1.
C
C       N = Number of nodes in the triangulation, including
C           constraint nodes.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations, followed by NCC se-
C             quences of constraint nodes.  Only one of
C             these sequences may be specified in clockwise
C             order to represent an exterior constraint
C             curve (a constraint region with nonfinite
C             area).
C
C The above parameters are not altered by this routine.
C
C       LWK = Length of IWK.  This must be at least 2*NI
C             where NI is the maximum number of arcs which
C             intersect a constraint arc to be added.  NI
C             is bounded by N-3.
C
C       IWK = Integer work array of length LWK (used by
C             Subroutine EDGE to add constraint arcs).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C On output:
C
C       LWK = Required length of IWK unless IER = 1 or IER =
C             3.  In the case of IER = 1, LWK is not altered
C             from its input value.
C
C       IWK = Array containing the endpoint indexes of the
C             new arcs which were swapped in by the last
C             call to Subroutine EDGE.
C
C       LIST,LPTR,LEND = Triangulation data structure with
C                        all constraint arcs present unless
C                        IER .NE. 0.  These arrays are not
C                        altered if IER = 1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if NCC, N, or an LCC entry is outside
C                     its valid range, or LWK .LT. 0 on
C                     input.
C             IER = 2 if more space is required in IWK.
C             IER = 3 if the triangulation data structure is
C                     invalid, or failure (in EDGE or OPTIM)
C                     was caused by collinear nodes on the
C                     convex hull boundary.  An error mes-
C                     sage is written to logical unit 6 in
C                     this case.
C             IER = 4 if intersecting constraint arcs were
C                     encountered.
C             IER = 5 if a constraint region contains a
C                     node.
C
C Modules required by ADDCST:  EDGE, LEFT, LSTPTR, OPTIM,
C                                SWAP, SWPTST
C
C Intrinsic functions called by ADDCST:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, IFRST, ILAST, K, KBAK, KFOR, KN, LCCIP1,
     .        LP, LPB, LPF, LPL, LW, LWD2, N1, N2
      LWD2 = LWK/2
C
C Test for errors in input parameters.
C
      IER = 1
      IF (NCC .LT. 0  .OR.  LWK .LT. 0) RETURN
      IF (NCC .EQ. 0) THEN
        IF (N .LT. 3) RETURN
        LWK = 0
        GO TO 9
      ELSE
        LCCIP1 = N+1
        DO 1 I = NCC,1,-1
          IF (LCCIP1 - LCC(I) .LT. 3) RETURN
          LCCIP1 = LCC(I)
    1     CONTINUE
        IF (LCCIP1 .LT. 1) RETURN
      ENDIF
C
C Force the presence of constraint arcs.  The outer loop is
C   on constraints in reverse order.  IFRST and ILAST are
C   the first and last nodes of constraint I.
C
      LWK = 0
      IFRST = N+1
      DO 3 I = NCC,1,-1
        ILAST = IFRST - 1
        IFRST = LCC(I)
C
C   Inner loop on constraint arcs N1-N2 in constraint I.
C
        N1 = ILAST
        DO 2 N2 = IFRST,ILAST
          LW = LWD2
          CALL EDGE (N1,N2,X,Y, LW,IWK,LIST,LPTR,LEND, IER)
          LWK = MAX(LWK,2*LW)
          IF (IER .EQ. 4) IER = 3
          IF (IER .NE. 0) RETURN
          N1 = N2
    2     CONTINUE
    3   CONTINUE
C
C Test for errors.  The outer loop is on constraint I with
C   first and last nodes IFRST and ILAST, and the inner loop
C   is on constraint nodes K with (KBAK,K,KFOR) a subse-
C   quence of constraint I.
C
      IER = 4
      IFRST = N+1
      DO 8 I = NCC,1,-1
        ILAST = IFRST - 1
        IFRST = LCC(I)
        KBAK = ILAST
        DO 7 K = IFRST,ILAST
          KFOR = K + 1
          IF (K .EQ. ILAST) KFOR = IFRST
C
C   Find the LIST pointers LPF and LPB of KFOR and KBAK as
C     neighbors of K.
C
          LPF = 0
          LPB = 0
          LPL = LEND(K)
          LP = LPL
C
    4     LP = LPTR(LP)
            KN = ABS(LIST(LP))
            IF (KN .EQ. KFOR) LPF = LP
            IF (KN .EQ. KBAK) LPB = LP
            IF (LP .NE. LPL) GO TO 4
C
C   A pair of intersecting constraint arcs was encountered
C     if and only if a constraint arc is missing (introduc-
C     tion of the second caused the first to be swapped out).
C
          IF (LPF .EQ. 0  .OR.  LPB .EQ. 0) RETURN
C
C   Loop on neighbors KN of node K which follow KFOR and
C     precede KBAK.  The constraint region contains no nodes
C     if and only if all such nodes KN are in constraint I.
C
          LP = LPF
    5     LP = LPTR(LP)
            IF (LP .EQ. LPB) GO TO 6
            KN = ABS(LIST(LP))
            IF (KN .LT. IFRST  .OR.  KN .GT. ILAST) GO TO 10
            GO TO 5
C
C   Bottom of loop.
C
    6     KBAK = K
    7     CONTINUE
    8   CONTINUE
C
C No errors encountered.
C
    9 IER = 0
      RETURN
C
C A constraint region contains a node.
C
   10 IER = 5
      RETURN
      END
      SUBROUTINE ADDNOD (K,XK,YK,IST,NCC, LCC,N,X,Y,LIST,
     .                   LPTR,LEND,LNEW, IER)
      INTEGER K, IST, NCC, LCC(*), N, LIST(*), LPTR(*),
     .        LEND(*), LNEW, IER
      REAL    XK, YK, X(*), Y(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/27/98
C
C   Given a triangulation of N nodes in the plane created by
C Subroutine TRMESH or TRMSHR, this subroutine updates the
C data structure with the addition of a new node in position
C K.  If node K is inserted into X and Y (K .LE. N) rather
C than appended (K = N+1), then a corresponding insertion
C must be performed in any additional arrays associated
C with the nodes.  For example, an array of data values Z
C must be shifted down to open up position K for the new
C value:  set Z(I+1) to Z(I) for I = N,N-1,...,K.  For
C optimal efficiency, new nodes should be appended whenever
C possible.  Insertion is necessary, however, to add a non-
C constraint node when constraints are present (refer to
C Subroutine ADDCST).
C
C   Note that a constraint node cannot be added by this
C routine.  In order to insert a constraint node, it is
C necessary to add the node with no constraints present
C (call this routine with NCC = 0), update LCC by increment-
C ing the appropriate entries, and then create (or restore)
C the constraints by a call to ADDCST.
C
C   The algorithm consists of the following steps:  node K
C is located relative to the triangulation (TRFIND), its
C index is added to the data structure (INTADD or BDYADD),
C and a sequence of swaps (SWPTST and SWAP) are applied to
C the arcs opposite K so that all arcs incident on node K
C and opposite node K (excluding constraint arcs) are local-
C ly optimal (satisfy the circumcircle test).  Thus, if a
C (constrained) Delaunay triangulation is input, a (con-
C strained) Delaunay triangulation will result.  All indexes
C are incremented as necessary for an insertion.
C
C
C On input:
C
C       K = Nodal index (index for X, Y, and LEND) of the
C           new node to be added.  1 .LE. K .LE. LCC(1).
C           (K .LE. N+1 if NCC=0).
C
C       XK,YK = Cartesian coordinates of the new node (to be
C               stored in X(K) and Y(K)).  The node must not
C               lie in a constraint region.
C
C       IST = Index of a node at which TRFIND begins the
C             search.  Search time depends on the proximity
C             of this node to node K.  1 .LE. IST .LE. N.
C
C       NCC = Number of constraint curves.  NCC .GE. 0.
C
C The above parameters are not altered by this routine.
C
C       LCC = List of constraint curve starting indexes (or
C             dummy array of length 1 if NCC = 0).  Refer to
C             Subroutine ADDCST.
C
C       N = Number of nodes in the triangulation before K is
C           added.  N .GE. 3.  Note that N will be incre-
C           mented following the addition of node K.
C
C       X,Y = Arrays of length at least N+1 containing the
C             Cartesian coordinates of the nodes in the
C             first N positions with non-constraint nodes
C             in the first LCC(1)-1 locations if NCC > 0.
C
C       LIST,LPTR,LEND,LNEW = Data structure associated with
C                             the triangulation of nodes 1
C                             to N.  The arrays must have
C                             sufficient length for N+1
C                             nodes.  Refer to TRMESH.
C
C On output:
C
C       LCC = List of constraint curve starting indexes in-
C             cremented by 1 to reflect the insertion of K
C             unless NCC = 0 or (IER .NE. 0 and IER .NE.
C             -4).
C
C       N = Number of nodes in the triangulation including K
C           unless IER .NE. 0 and IER .NE. -4.  Note that
C           all comments refer to the input value of N.
C
C       X,Y = Arrays updated with the insertion of XK and YK
C             in the K-th positions (node I+1 was node I be-
C             fore the insertion for I = K to N if K .LE. N)
C             unless IER .NE. 0 and IER .NE. -4.
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the addition of node K unless
C                             IER .NE. 0 and IER .NE. -4.
C
C       IER = Error indicator:
C             IER =  0 if no errors were encountered.
C             IER = -1 if K, IST, NCC, N, or an LCC entry is
C                      outside its valid range on input.
C             IER = -2 if all nodes (including K) are col-
C                      linear.
C             IER =  L if nodes L and K coincide for some L.
C             IER = -3 if K lies in a constraint region.
C             IER = -4 if an error flag is returned by SWAP
C                      implying that the triangulation
C                      (geometry) was bad on input.
C
C             The errors conditions are tested in the order
C             specified.
C
C Modules required by ADDNOD:  BDYADD, CRTRI, INDXCC,
C                                INSERT, INTADD, JRAND,
C                                LEFT, LSTPTR, SWAP,
C                                SWPTST, TRFIND
C
C Intrinsic function called by ADDNOD:  ABS
C
C***********************************************************
C
      INTEGER INDXCC, LSTPTR
      INTEGER I, I1, I2, I3, IBK, IO1, IO2, IN1, KK, L,
     .        LCCIP1, LP, LPF, LPO1, NM1
      LOGICAL CRTRI, SWPTST
      KK = K
C
C Test for an invalid input parameter.
C
      IF (KK .LT. 1  .OR.  IST .LT. 1  .OR.  IST .GT. N
     .    .OR.  NCC .LT. 0  .OR.  N .LT. 3) GO TO 7
      LCCIP1 = N+1
      DO 1 I = NCC,1,-1
        IF (LCCIP1-LCC(I) .LT. 3) GO TO 7
        LCCIP1 = LCC(I)
    1   CONTINUE
      IF (KK .GT. LCCIP1) GO TO 7
C
C Find a triangle (I1,I2,I3) containing K or the rightmost
C   (I1) and leftmost (I2) visible boundary nodes as viewed
C   from node K.
C
      CALL TRFIND (IST,XK,YK,N,X,Y,LIST,LPTR,LEND, I1,I2,I3)
C
C Test for collinear nodes, duplicate nodes, and K lying in
C   a constraint region.
C
      IF (I1 .EQ. 0) GO TO 8
      IF (I3 .NE. 0) THEN
        L = I1
        IF (XK .EQ. X(L)  .AND.  YK .EQ. Y(L)) GO TO 9
        L = I2
        IF (XK .EQ. X(L)  .AND.  YK .EQ. Y(L)) GO TO 9
        L = I3
        IF (XK .EQ. X(L)  .AND.  YK .EQ. Y(L)) GO TO 9
        IF (NCC .GT. 0  .AND.  CRTRI(NCC,LCC,I1,I2,I3) )
     .    GO TO 10
      ELSE
C
C   K is outside the convex hull of the nodes and lies in a
C     constraint region iff an exterior constraint curve is
C     present.
C
        IF (NCC .GT. 0  .AND.  INDXCC(NCC,LCC,N,LIST,LEND)
     .      .NE. 0) GO TO 10
      ENDIF
C
C No errors encountered.
C
      IER = 0
      NM1 = N
      N = N + 1
      IF (KK .LT. N) THEN
C
C Open a slot for K in X, Y, and LEND, and increment all
C   nodal indexes which are greater than or equal to K.
C   Note that LIST, LPTR, and LNEW are not yet updated with
C   either the neighbors of K or the edges terminating on K.
C
        DO 2 IBK = NM1,KK,-1
          X(IBK+1) = X(IBK)
          Y(IBK+1) = Y(IBK)
          LEND(IBK+1) = LEND(IBK)
    2     CONTINUE
        DO 3 I = 1,NCC
          LCC(I) = LCC(I) + 1
    3     CONTINUE
        L = LNEW - 1
        DO 4 I = 1,L
          IF (LIST(I) .GE. KK) LIST(I) = LIST(I) + 1
          IF (LIST(I) .LE. -KK) LIST(I) = LIST(I) - 1
    4     CONTINUE
        IF (I1 .GE. KK) I1 = I1 + 1
        IF (I2 .GE. KK) I2 = I2 + 1
        IF (I3 .GE. KK) I3 = I3 + 1
      ENDIF
C
C Insert K into X and Y, and update LIST, LPTR, LEND, and
C   LNEW with the arcs containing node K.
C
      X(KK) = XK
      Y(KK) = YK
      IF (I3 .EQ. 0) THEN
        CALL BDYADD (KK,I1,I2, LIST,LPTR,LEND,LNEW )
      ELSE
        CALL INTADD (KK,I1,I2,I3, LIST,LPTR,LEND,LNEW )
      ENDIF
C
C Initialize variables for optimization of the triangula-
C   tion.
C
      LP = LEND(KK)
      LPF = LPTR(LP)
      IO2 = LIST(LPF)
      LPO1 = LPTR(LPF)
      IO1 = ABS(LIST(LPO1))
C
C Begin loop:  find the node opposite K.
C
    5 LP = LSTPTR(LEND(IO1),IO2,LIST,LPTR)
        IF (LIST(LP) .LT. 0) GO TO 6
        LP = LPTR(LP)
        IN1 = ABS(LIST(LP))
        IF ( CRTRI(NCC,LCC,IO1,IO2,IN1) ) GO TO 6
C
C Swap test:  if a swap occurs, two new arcs are
C             opposite K and must be tested.
C
        IF ( .NOT. SWPTST(IN1,KK,IO1,IO2,X,Y) ) GO TO 6
        CALL SWAP (IN1,KK,IO1,IO2, LIST,LPTR,LEND, LPO1)
        IF (LPO1 .EQ. 0) GO TO 11
        IO1 = IN1
        GO TO 5
C
C No swap occurred.  Test for termination and reset
C   IO2 and IO1.
C
    6   IF (LPO1 .EQ. LPF  .OR.  LIST(LPO1) .LT. 0) RETURN
        IO2 = IO1
        LPO1 = LPTR(LPO1)
        IO1 = ABS(LIST(LPO1))
        GO TO 5
C
C A parameter is outside its valid range on input.
C
    7 IER = -1
      RETURN
C
C All nodes are collinear.
C
    8 IER = -2
      RETURN
C
C Nodes L and K coincide.
C
    9 IER = L
      RETURN
C
C Node K lies in a constraint region.
C
   10 IER = -3
      RETURN
C
C Zero pointer returned by SWAP.
C
   11 IER = -4
      RETURN
      END
      REAL FUNCTION AREAP (X,Y,NB,NODES)
      INTEGER NB, NODES(NB)
      REAL    X(*), Y(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/21/90
C
C   Given a sequence of NB points in the plane, this func-
C tion computes the signed area bounded by the closed poly-
C gonal curve which passes through the points in the
C specified order.  Each simple closed curve is positively
C oriented (bounds positive area) if and only if the points
C are specified in counterclockwise order.  The last point
C of the curve is taken to be the first point specified, and
C this point should therefore not be specified twice.
C
C   The area of a triangulation may be computed by calling
C AREAP with values of NB and NODES determined by Subroutine
C BNODES.
C
C
C On input:
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of a set of points in the plane
C             for some N .GE. NB.
C
C       NB = Length of NODES.
C
C       NODES = Array of length NB containing the ordered
C               sequence of nodal indexes (in the range
C               1 to N) which define the polygonal curve.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       AREAP = Signed area bounded by the polygonal curve,
C              or zero if NB < 3.
C
C Modules required by AREAP:  None
C
C***********************************************************
C
      INTEGER I, ND1, ND2, NNB
      REAL    A
C
C Local parameters:
C
C A =       Partial sum of signed (and doubled) trapezoid
C             areas
C I =       DO-loop and NODES index
C ND1,ND2 = Elements of NODES
C NNB =     Local copy of NB
C
      NNB = NB
      A = 0.
      IF (NNB .LT. 3) GO TO 2
      ND2 = NODES(NNB)
C
C Loop on line segments NODES(I-1) -> NODES(I), where
C   NODES(0) = NODES(NB), adding twice the signed trapezoid
C   areas (integrals of the linear interpolants) to A.
C
      DO 1 I = 1,NNB
        ND1 = ND2
        ND2 = NODES(I)
        A = A + (X(ND2)-X(ND1))*(Y(ND1)+Y(ND2))
    1   CONTINUE
C
C A contains twice the negative signed area of the region.
C
    2 AREAP = -A/2.
      RETURN
      END
      SUBROUTINE BDYADD (KK,I1,I2, LIST,LPTR,LEND,LNEW )
      INTEGER KK, I1, I2, LIST(*), LPTR(*), LEND(*), LNEW
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/22/91
C
C   This subroutine adds a boundary node to a triangulation
C of a set of points in the plane.  The data structure is
C updated with the insertion of node KK, but no optimization
C is performed.
C
C
C On input:
C
C       KK = Index of a node to be connected to the sequence
C            of all visible boundary nodes.  KK .GE. 1 and
C            KK must not be equal to I1 or I2.
C
C       I1 = First (rightmost as viewed from KK) boundary
C            node in the triangulation which is visible from
C            node KK (the line segment KK-I1 intersects no
C            arcs.
C
C       I2 = Last (leftmost) boundary node which is visible
C            from node KK.  I1 and I2 may be determined by
C            Subroutine TRFIND.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Triangulation data structure
C                             created by TRMESH or TRMSHR.
C                             Nodes I1 and I2 must be in-
C                             cluded in the triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the addition of node KK.  Node
C                             KK is connected to I1, I2, and
C                             all boundary nodes in between.
C
C Module required by BDYADD:  INSERT
C
C***********************************************************
C
      INTEGER K, LP, LSAV, N1, N2, NEXT, NSAV
      K = KK
      N1 = I1
      N2 = I2
C
C Add K as the last neighbor of N1.
C
      LP = LEND(N1)
      LSAV = LPTR(LP)
      LPTR(LP) = LNEW
      LIST(LNEW) = -K
      LPTR(LNEW) = LSAV
      LEND(N1) = LNEW
      LNEW = LNEW + 1
      NEXT = -LIST(LP)
      LIST(LP) = NEXT
      NSAV = NEXT
C
C Loop on the remaining boundary nodes between N1 and N2,
C   adding K as the first neighbor.
C
    1 LP = LEND(NEXT)
        CALL INSERT (K,LP,LIST,LPTR,LNEW)
        IF (NEXT .EQ. N2) GO TO 2
        NEXT = -LIST(LP)
        LIST(LP) = NEXT
        GO TO 1
C
C Add the boundary nodes between N1 and N2 as neighbors
C   of node K.
C
    2 LSAV = LNEW
      LIST(LNEW) = N1
      LPTR(LNEW) = LNEW + 1
      LNEW = LNEW + 1
      NEXT = NSAV
C
    3 IF (NEXT .EQ. N2) GO TO 4
        LIST(LNEW) = NEXT
        LPTR(LNEW) = LNEW + 1
        LNEW = LNEW + 1
        LP = LEND(NEXT)
        NEXT = LIST(LP)
        GO TO 3
C
    4 LIST(LNEW) = -N2
      LPTR(LNEW) = LSAV
      LEND(K) = LNEW
      LNEW = LNEW + 1
      RETURN
      END
      SUBROUTINE BNODES (N,LIST,LPTR,LEND, NODES,NB,NA,NT)
      INTEGER N, LIST(*), LPTR(*), LEND(N), NODES(*), NB,
     .        NA, NT
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   Given a triangulation of N points in the plane, this
C subroutine returns an array containing the indexes, in
C counterclockwise order, of the nodes on the boundary of
C the convex hull of the set of points.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C The above parameters are not altered by this routine.
C
C       NODES = Integer array of length at least NB
C               (NB .LE. N).
C
C On output:
C
C       NODES = Ordered sequence of boundary node indexes
C               in the range 1 to N.
C
C       NB = Number of boundary nodes.
C
C       NA,NT = Number of arcs and triangles, respectively,
C               in the triangulation.
C
C Modules required by BNODES:  None
C
C***********************************************************
C
      INTEGER K, LP, N0, NST
C
C Set NST to the first boundary node encountered.
C
      NST = 1
    1 LP = LEND(NST)
        IF (LIST(LP) .LT. 0) GO TO 2
        NST = NST + 1
        GO TO 1
C
C Initialization.
C
    2 NODES(1) = NST
      K = 1
      N0 = NST
C
C Traverse the boundary in counterclockwise order.
C
    3 LP = LEND(N0)
        LP = LPTR(LP)
        N0 = LIST(LP)
        IF (N0 .EQ. NST) GO TO 4
        K = K + 1
        NODES(K) = N0
        GO TO 3
C
C Termination.
C
    4 NB = K
      NT = 2*N - NB - 2
      NA = NT + N - 1
      RETURN
      END
      SUBROUTINE CIRCUM (X1,Y1,X2,Y2,X3,Y3,RATIO, XC,YC,CR,
     .                   SA,AR)
      LOGICAL RATIO
      REAL    X1, Y1, X2, Y2, X3, Y3, XC, YC, CR, SA, AR
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   12/10/96
C
C   Given three vertices defining a triangle, this subrou-
C tine returns the circumcenter, circumradius, signed
C triangle area, and, optionally, the aspect ratio of the
C triangle.
C
C
C On input:
C
C       X1,...,Y3 = Cartesian coordinates of the vertices.
C
C       RATIO = Logical variable with value TRUE if and only
C               if the aspect ratio is to be computed.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       XC,YC = Cartesian coordinates of the circumcenter
C               (center of the circle defined by the three
C               points) unless SA = 0, in which XC and YC
C               are not altered.
C
C       CR = Circumradius (radius of the circle defined by
C            the three points) unless SA = 0 (infinite
C            radius), in which case CR is not altered.
C
C       SA = Signed triangle area with positive value if
C            and only if the vertices are specified in
C            counterclockwise order:  (X3,Y3) is strictly
C            to the left of the directed line from (X1,Y1)
C            toward (X2,Y2).
C
C       AR = Aspect ratio r/CR, where r is the radius of the
C            inscribed circle, unless RATIO = FALSE, in
C            which case AR is not altered.  AR is in the
C            range 0 to .5, with value 0 iff SA = 0 and
C            value .5 iff the vertices define an equilateral
C            triangle.
C
C Modules required by CIRCUM:  None
C
C Intrinsic functions called by CIRCUM:  ABS, SQRT
C
C***********************************************************
C
      INTEGER I
      REAL    DS(3), FX, FY, U(3), V(3)
C
C Set U(K) and V(K) to the x and y components, respectively,
C   of the directed edge opposite vertex K.
C
      U(1) = X3 - X2
      U(2) = X1 - X3
      U(3) = X2 - X1
      V(1) = Y3 - Y2
      V(2) = Y1 - Y3
      V(3) = Y2 - Y1
C
C Set SA to the signed triangle area.
C
      SA = (U(1)*V(2) - U(2)*V(1))/2.
      IF (SA .EQ. 0.) THEN
        IF (RATIO) AR = 0.
        RETURN
      ENDIF
C
C Set DS(K) to the squared distance from the origin to
C   vertex K.
C
      DS(1) = X1*X1 + Y1*Y1
      DS(2) = X2*X2 + Y2*Y2
      DS(3) = X3*X3 + Y3*Y3
C
C Compute factors of XC and YC.
C
      FX = 0.
      FY = 0.
      DO 1 I = 1,3
        FX = FX - DS(I)*V(I)
        FY = FY + DS(I)*U(I)
    1   CONTINUE
      XC = FX/(4.*SA)
      YC = FY/(4.*SA)
      CR = SQRT( (XC-X1)**2 + (YC-Y1)**2 )
      IF (.NOT. RATIO) RETURN
C
C Compute the squared edge lengths and aspect ratio.
C
      DO 2 I = 1,3
        DS(I) = U(I)*U(I) + V(I)*V(I)
    2   CONTINUE
      AR = 2.*ABS(SA)/
     .     ( (SQRT(DS(1)) + SQRT(DS(2)) + SQRT(DS(3)))*CR )
      RETURN
      END
      LOGICAL FUNCTION CRTRI (NCC,LCC,I1,I2,I3)
      INTEGER NCC, LCC(*), I1, I2, I3
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   08/14/91
C
C   This function returns TRUE if and only if triangle (I1,
C I2,I3) lies in a constraint region.
C
C
C On input:
C
C       NCC,LCC = Constraint data structure.  Refer to Sub-
C                 routine ADDCST.
C
C       I1,I2,I3 = Nodal indexes of the counterclockwise-
C                  ordered vertices of a triangle.
C
C Input parameters are altered by this function.
C
C       CRTRI = TRUE iff (I1,I2,I3) is a constraint region
C               triangle.
C
C Note that input parameters are not tested for validity.
C
C Modules required by CRTRI:  None
C
C Intrinsic functions called by CRTRI:  MAX, MIN
C
C***********************************************************
C
      INTEGER I, IMAX, IMIN
      IMAX = MAX(I1,I2,I3)
C
C   Find the index I of the constraint containing IMAX.
C
      I = NCC + 1
    1 I = I - 1
        IF (I .LE. 0) GO TO 2
        IF (IMAX .LT. LCC(I)) GO TO 1
      IMIN = MIN(I1,I2,I3)
C
C P lies in a constraint region iff I1, I2, and I3 are nodes
C   of the same constraint (IMIN >= LCC(I)), and (IMIN,IMAX)
C   is (I1,I3), (I2,I1), or (I3,I2).
C
      CRTRI = IMIN .GE. LCC(I)  .AND.  ((IMIN .EQ. I1 .AND.
     .        IMAX .EQ. I3)  .OR.  (IMIN .EQ. I2  .AND.
     .        IMAX .EQ. I1)  .OR.  (IMIN .EQ. I3  .AND.
     .        IMAX .EQ. I2))
      RETURN
C
C NCC .LE. 0 or all vertices are non-constraint nodes.
C
    2 CRTRI = .FALSE.
      RETURN
      END
      SUBROUTINE DELARC (N,IO1,IO2, LIST,LPTR,LEND,
     .                   LNEW, IER)
      INTEGER N, IO1, IO2, LIST(*), LPTR(*), LEND(N), LNEW,
     .        IER
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/12/94
C
C   This subroutine deletes a boundary arc from a triangula-
C tion.  It may be used to remove a null triangle from the
C convex hull boundary.  Note, however, that if the union of
C triangles is rendered nonconvex, Subroutines DELNOD, EDGE,
C and TRFIND may fail.  Thus, Subroutines ADDCST, ADDNOD,
C DELNOD, EDGE, and NEARND should not be called following
C an arc deletion.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 4.
C
C       IO1,IO2 = Indexes (in the range 1 to N) of a pair of
C                 adjacent boundary nodes defining the arc
C                 to be removed.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Triangulation data structure
C                             created by TRMESH or TRMSHR.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the removal of arc IO1-IO2
C                             unless IER > 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N, IO1, or IO2 is outside its valid
C                     range, or IO1 = IO2.
C             IER = 2 if IO1-IO2 is not a boundary arc.
C             IER = 3 if the node opposite IO1-IO2 is al-
C                     ready a boundary node, and thus IO1
C                     or IO2 has only two neighbors or a
C                     deletion would result in two triangu-
C                     lations sharing a single node.
C             IER = 4 if one of the nodes is a neighbor of
C                     the other, but not vice versa, imply-
C                     ing an invalid triangulation data
C                     structure.
C
C Modules required by DELARC:  DELNB, LSTPTR
C
C Intrinsic function called by DELARC:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER LP, LPH, LPL, N1, N2, N3
      N1 = IO1
      N2 = IO2
C
C Test for errors, and set N1->N2 to the directed boundary
C   edge associated with IO1-IO2:  (N1,N2,N3) is a triangle
C   for some N3.
C
      IF (N .LT. 4  .OR.  N1 .LT. 1  .OR.  N1 .GT. N  .OR.
     .    N2 .LT. 1  .OR.  N2 .GT. N  .OR.  N1 .EQ. N2) THEN
        IER = 1
        RETURN
      ENDIF
C
      LPL = LEND(N2)
      IF (-LIST(LPL) .NE. N1) THEN
        N1 = N2
        N2 = IO1
        LPL = LEND(N2)
        IF (-LIST(LPL) .NE. N1) THEN
          IER = 2
          RETURN
        ENDIF
      ENDIF
C
C Set N3 to the node opposite N1->N2 (the second neighbor
C   of N1), and test for error 3 (N3 already a boundary
C   node).
C
      LPL = LEND(N1)
      LP = LPTR(LPL)
      LP = LPTR(LP)
      N3 = ABS(LIST(LP))
      LPL = LEND(N3)
      IF (LIST(LPL) .LE. 0) THEN
        IER = 3
        RETURN
      ENDIF
C
C Delete N2 as a neighbor of N1, making N3 the first
C   neighbor, and test for error 4 (N2 not a neighbor
C   of N1).  Note that previously computed pointers may
C   no longer be valid following the call to DELNB.
C
      CALL DELNB (N1,N2,N, LIST,LPTR,LEND,LNEW, LPH)
      IF (LPH .LT. 0) THEN
        IER = 4
        RETURN
      ENDIF
C
C Delete N1 as a neighbor of N2, making N3 the new last
C   neighbor.
C
      CALL DELNB (N2,N1,N, LIST,LPTR,LEND,LNEW, LPH)
C
C Make N3 a boundary node with first neighbor N2 and last
C   neighbor N1.
C
      LP = LSTPTR(LEND(N3),N1,LIST,LPTR)
      LEND(N3) = LP
      LIST(LP) = -N1
C
C No errors encountered.
C
      IER = 0
      RETURN
      END
      SUBROUTINE DELNB (N0,NB,N, LIST,LPTR,LEND,LNEW, LPH)
      INTEGER N0, NB, N, LIST(*), LPTR(*), LEND(N), LNEW,
     .        LPH
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/30/98
C
C   This subroutine deletes a neighbor NB from the adjacency
C list of node N0 (but N0 is not deleted from the adjacency
C list of NB) and, if NB is a boundary node, makes N0 a
C boundary node.  For pointer (LIST index) LPH to NB as a
C neighbor of N0, the empty LIST,LPTR location LPH is filled
C in with the values at LNEW-1, pointer LNEW-1 (in LPTR and
C possibly in LEND) is changed to LPH, and LNEW is decremen-
C ted.  This requires a search of LEND and LPTR entailing an
C expected operation count of O(N).
C
C
C On input:
C
C       N0,NB = Indexes, in the range 1 to N, of a pair of
C               nodes such that NB is a neighbor of N0.
C               (N0 need not be a neighbor of NB.)
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the removal of NB from the ad-
C                             jacency list of N0 unless
C                             LPH < 0.
C
C       LPH = List pointer to the hole (NB as a neighbor of
C             N0) filled in by the values at LNEW-1 or error
C             indicator:
C             LPH > 0 if no errors were encountered.
C             LPH = -1 if N0, NB, or N is outside its valid
C                      range.
C             LPH = -2 if NB is not a neighbor of N0.
C
C Modules required by DELNB:  None
C
C Intrinsic function called by DELNB:  ABS
C
C***********************************************************
C
      INTEGER I, LNW, LP, LPB, LPL, LPP, NN
C
C Local parameters:
C
C I =   DO-loop index
C LNW = LNEW-1 (output value of LNEW)
C LP =  LIST pointer of the last neighbor of NB
C LPB = Pointer to NB as a neighbor of N0
C LPL = Pointer to the last neighbor of N0
C LPP = Pointer to the neighbor of N0 that precedes NB
C NN =  Local copy of N
C
      NN = N
C
C Test for error 1.
C
      IF (N0 .LT. 1  .OR.  N0 .GT. NN  .OR.  NB .LT. 1  .OR.
     .    NB .GT. NN  .OR.  NN .LT. 3) THEN
        LPH = -1
        RETURN
      ENDIF
C
C   Find pointers to neighbors of N0:
C
C     LPL points to the last neighbor,
C     LPP points to the neighbor NP preceding NB, and
C     LPB points to NB.
C
      LPL = LEND(N0)
      LPP = LPL
      LPB = LPTR(LPP)
    1 IF (LIST(LPB) .EQ. NB) GO TO 2
        LPP = LPB
        LPB = LPTR(LPP)
        IF (LPB .NE. LPL) GO TO 1
C
C   Test for error 2 (NB not found).
C
      IF (ABS(LIST(LPB)) .NE. NB) THEN
        LPH = -2
        RETURN
      ENDIF
C
C   NB is the last neighbor of N0.  Make NP the new last
C     neighbor and, if NB is a boundary node, then make N0
C     a boundary node.
C
      LEND(N0) = LPP
      LP = LEND(NB)
      IF (LIST(LP) .LT. 0) LIST(LPP) = -LIST(LPP)
      GO TO 3
C
C   NB is not the last neighbor of N0.  If NB is a boundary
C     node and N0 is not, then make N0 a boundary node with
C     last neighbor NP.
C
    2 LP = LEND(NB)
      IF (LIST(LP) .LT. 0  .AND.  LIST(LPL) .GT. 0) THEN
        LEND(N0) = LPP
        LIST(LPP) = -LIST(LPP)
      ENDIF
C
C   Update LPTR so that the neighbor following NB now fol-
C     lows NP, and fill in the hole at location LPB.
C
    3 LPTR(LPP) = LPTR(LPB)
      LNW = LNEW-1
      LIST(LPB) = LIST(LNW)
      LPTR(LPB) = LPTR(LNW)
      DO 4 I = NN,1,-1
        IF (LEND(I) .EQ. LNW) THEN
          LEND(I) = LPB
          GO TO 5
        ENDIF
    4   CONTINUE
C
    5 DO 6 I = 1,LNW-1
        IF (LPTR(I) .EQ. LNW) THEN
          LPTR(I) = LPB
        ENDIF
    6   CONTINUE
C
C No errors encountered.
C
      LNEW = LNW
      LPH = LPB
      RETURN
      END
      SUBROUTINE DELNOD (K,NCC, LCC,N,X,Y,LIST,LPTR,LEND,
     .                   LNEW,LWK,IWK, IER)
      INTEGER K, NCC, LCC(*), N, LIST(*), LPTR(*),
     .        LEND(*), LNEW, LWK, IWK(2,*), IER
      REAL    X(*), Y(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/28/98
C
C   This subroutine deletes node K (along with all arcs
C incident on node K) from a triangulation of N nodes in the
C plane, and inserts arcs as necessary to produce a triangu-
C lation of the remaining N-1 nodes.  If a Delaunay triangu-
C lation is input, a Delaunay triangulation will result, and
C thus, DELNOD reverses the effect of a call to Subroutine
C ADDNOD.
C
C   Note that a constraint node cannot be deleted by this
C routine.  In order to delete a constraint node, it is
C necessary to call this routine with NCC = 0, decrement the
C appropriate LCC entries (LCC(I) such that LCC(I) > K), and
C then create (or restore) the constraints by a call to Sub-
C routine ADDCST.
C
C
C On input:
C
C       K = Index (for X and Y) of the node to be deleted.
C           1 .LE. K .LT. LCC(1).  (K .LE. N if NCC=0).
C
C       NCC = Number of constraint curves.  NCC .GE. 0.
C
C The above parameters are not altered by this routine.
C
C       LCC = List of constraint curve starting indexes (or
C             dummy array of length 1 if NCC = 0).  Refer to
C             Subroutine ADDCST.
C
C       N = Number of nodes in the triangulation on input.
C           N .GE. 4.  Note that N will be decremented
C           following the deletion.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations if NCC > 0.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.  Refer to Sub-
C                             routine TRMESH.
C
C       LWK = Number of columns reserved for IWK.  LWK must
C             be at least NNB-3, where NNB is the number of
C             neighbors of node K, including an extra
C             pseudo-node if K is a boundary node.
C
C       IWK = Integer work array dimensioned 2 by LWK (or
C             array of length .GE. 2*LWK).
C
C On output:
C
C       LCC = List of constraint curve starting indexes de-
C             cremented by 1 to reflect the deletion of K
C             unless NCC = 0 or 1 .LE. IER .LE. 4.
C
C       N = New number of nodes (input value minus one) un-
C           less 1 .LE. IER .LE. 4.
C
C       X,Y = Updated arrays of length N-1 containing nodal
C             coordinates (with elements K+1,...,N shifted
C             up a position and thus overwriting element K)
C             unless 1 .LE. IER .LE. 4.  (N here denotes the
C             input value.)
C
C       LIST,LPTR,LEND,LNEW = Updated triangulation data
C                             structure reflecting the dele-
C                             tion unless IER .NE. 0.  Note
C                             that the data structure may
C                             have been altered if IER .GE.
C                             3.
C
C       LWK = Number of IWK columns required unless IER = 1
C             or IER = 3.
C
C       IWK = Indexes of the endpoints of the new arcs added
C             unless LWK = 0 or 1 .LE. IER .LE. 4.  (Arcs
C             are associated with columns, or pairs of
C             adjacent elements if IWK is declared as a
C             singly-subscripted array.)
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if K, NCC, N, or an LCC entry is out-
C                     side its valid range or LWK < 0 on
C                     input.
C             IER = 2 if more space is required in IWK.
C                     Refer to LWK.
C             IER = 3 if the triangulation data structure is
C                     invalid on input.
C             IER = 4 if K is an interior node with 4 or
C                     more neighbors, and the number of
C                     neighbors could not be reduced to 3
C                     by swaps.  This could be caused by
C                     floating point errors with collinear
C                     nodes or by an invalid data structure.
C             IER = 5 if an error flag was returned by
C                     OPTIM.  An error message is written
C                     to the standard output unit in this
C                     event.
C
C   Note that the deletion may result in all remaining nodes
C being collinear.  This situation is not flagged.
C
C Modules required by DELNOD:  DELNB, LEFT, LSTPTR, NBCNT,
C                                OPTIM, SWAP, SWPTST
C
C Intrinsic function called by DELNOD:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR, NBCNT
      LOGICAL LEFT
      INTEGER I, IERR, IWL, J, LCCIP1, LNW, LP, LP21, LPF,
     .        LPH, LPL, LPL2, LPN, LWKL, N1, N2, NFRST, NIT,
     .        NL, NN, NNB, NR
      LOGICAL BDRY
      REAL    X1, X2, XL, XR, Y1, Y2, YL, YR
C
C Set N1 to K and NNB to the number of neighbors of N1 (plus
C   one if N1 is a boundary node), and test for errors.  LPF
C   and LPL are LIST indexes of the first and last neighbors
C   of N1, IWL is the number of IWK columns containing arcs,
C   and BDRY is TRUE iff N1 is a boundary node.
C
      N1 = K
      NN = N
      IF (NCC .LT. 0  .OR.  N1 .LT. 1  .OR.  NN .LT. 4  .OR.
     .    LWK .LT. 0) GO TO 21
      LCCIP1 = NN+1
      DO 1 I = NCC,1,-1
        IF (LCCIP1-LCC(I) .LT. 3) GO TO 21
        LCCIP1 = LCC(I)
    1   CONTINUE
      IF (N1 .GE. LCCIP1) GO TO 21
      LPL = LEND(N1)
      LPF = LPTR(LPL)
      NNB = NBCNT(LPL,LPTR)
      BDRY = LIST(LPL) .LT. 0
      IF (BDRY) NNB = NNB + 1
      IF (NNB .LT. 3) GO TO 23
      LWKL = LWK
      LWK = NNB - 3
      IF (LWKL .LT. LWK) GO TO 22
      IWL = 0
      IF (NNB .EQ. 3) GO TO 5
C
C Initialize for loop on arcs N1-N2 for neighbors N2 of N1,
C   beginning with the second neighbor.  NR and NL are the
C   neighbors preceding and following N2, respectively, and
C   LP indexes NL.  The loop is exited when all possible
C   swaps have been applied to arcs incident on N1.  If N1
C   is interior, the number of neighbors will be reduced
C   to 3.
C
      X1 = X(N1)
      Y1 = Y(N1)
      NFRST = LIST(LPF)
      NR = NFRST
      XR = X(NR)
      YR = Y(NR)
      LP = LPTR(LPF)
      N2 = LIST(LP)
      X2 = X(N2)
      Y2 = Y(N2)
      LP = LPTR(LP)
C
C Top of loop:  set NL to the neighbor following N2.
C
    2 NL = ABS(LIST(LP))
      IF (NL .EQ. NFRST  .AND.  BDRY) GO TO 5
      XL = X(NL)
      YL = Y(NL)
C
C   Test for a convex quadrilateral.  To avoid an incorrect
C     test caused by collinearity, use the fact that if N1
C     is a boundary node, then N1 LEFT NR->NL and if N2 is
C     a boundary node, then N2 LEFT NL->NR.
C
      LPL2 = LEND(N2)
      IF ( (BDRY  .OR.  LEFT(XR,YR,XL,YL,X1,Y1))  .AND.
     .     (LIST(LPL2) .LT. 0  .OR.
     .      LEFT(XL,YL,XR,YR,X2,Y2)) ) GO TO 3
C
C   Nonconvex quadrilateral -- no swap is possible.
C
      NR = N2
      XR = X2
      YR = Y2
      GO TO 4
C
C   The quadrilateral defined by adjacent triangles
C     (N1,N2,NL) and (N2,N1,NR) is convex.  Swap in
C     NL-NR and store it in IWK.  Indexes larger than N1
C     must be decremented since N1 will be deleted from
C     X and Y.
C
    3 CALL SWAP (NL,NR,N1,N2, LIST,LPTR,LEND, LP21)
      IWL = IWL + 1
      IF (NL .LE. N1) THEN
        IWK(1,IWL) = NL
      ELSE
        IWK(1,IWL) = NL - 1
      ENDIF
      IF (NR .LE. N1) THEN
        IWK(2,IWL) = NR
      ELSE
        IWK(2,IWL) = NR - 1
      ENDIF
C
C   Recompute the LIST indexes LPL,LP and decrement NNB.
C
      LPL = LEND(N1)
      NNB = NNB - 1
      IF (NNB .EQ. 3) GO TO 5
      LP = LSTPTR(LPL,NL,LIST,LPTR)
      IF (NR .EQ. NFRST) GO TO 4
C
C   NR is not the first neighbor of N1.
C     Back up and test N1-NR for a swap again:  Set N2 to
C     NR and NR to the previous neighbor of N1 -- the
C     neighbor of NR which follows N1.  LP21 points to NL
C     as a neighbor of NR.
C
      N2 = NR
      X2 = XR
      Y2 = YR
      LP21 = LPTR(LP21)
      LP21 = LPTR(LP21)
      NR = ABS(LIST(LP21))
      XR = X(NR)
      YR = Y(NR)
      GO TO 2
C
C   Bottom of loop -- test for invalid termination.
C
    4 IF (N2 .EQ. NFRST) GO TO 24
      N2 = NL
      X2 = XL
      Y2 = YL
      LP = LPTR(LP)
      GO TO 2
C
C Delete N1 from the adjacency list of N2 for all neighbors
C   N2 of N1.  LPL points to the last neighbor of N1.
C   LNEW is stored in local variable LNW.
C
    5 LP = LPL
      LNW = LNEW
C
C Loop on neighbors N2 of N1, beginning with the first.
C
    6 LP = LPTR(LP)
        N2 = ABS(LIST(LP))
        CALL DELNB (N2,N1,N, LIST,LPTR,LEND,LNW, LPH)
        IF (LPH .LT. 0) GO TO 23
C
C   LP and LPL may require alteration.
C
        IF (LPL .EQ. LNW) LPL = LPH
        IF (LP .EQ. LNW) LP = LPH
        IF (LP .NE. LPL) GO TO 6
C
C Delete N1 from X, Y, and LEND, and remove its adjacency
C   list from LIST and LPTR.  LIST entries (nodal indexes)
C   which are larger than N1 must be decremented.
C
      NN = NN - 1
      IF (N1 .GT. NN) GO TO 9
      DO 7 I = N1,NN
        X(I) = X(I+1)
        Y(I) = Y(I+1)
        LEND(I) = LEND(I+1)
    7   CONTINUE
C
      DO 8 I = 1,LNW-1
        IF (LIST(I) .GT. N1) LIST(I) = LIST(I) - 1
        IF (LIST(I) .LT. -N1) LIST(I) = LIST(I) + 1
    8   CONTINUE
C
C   For LPN = first to last neighbors of N1, delete the
C     preceding neighbor (indexed by LP).
C
C   Each empty LIST,LPTR location LP is filled in with the
C     values at LNW-1, and LNW is decremented.  All pointers
C     (including those in LPTR and LEND) with value LNW-1
C     must be changed to LP.
C
C  LPL points to the last neighbor of N1.
C
    9 IF (BDRY) NNB = NNB - 1
      LPN = LPL
      DO 13 J = 1,NNB
        LNW = LNW - 1
        LP = LPN
        LPN = LPTR(LP)
        LIST(LP) = LIST(LNW)
        LPTR(LP) = LPTR(LNW)
        IF (LPTR(LPN) .EQ. LNW) LPTR(LPN) = LP
        IF (LPN .EQ. LNW) LPN = LP
        DO 10 I = NN,1,-1
          IF (LEND(I) .EQ. LNW) THEN
            LEND(I) = LP
            GO TO 11
          ENDIF
   10     CONTINUE
C
   11   DO 12 I = LNW-1,1,-1
          IF (LPTR(I) .EQ. LNW) LPTR(I) = LP
   12     CONTINUE
   13   CONTINUE
C
C Decrement LCC entries.
C
      DO 14 I = 1,NCC
        LCC(I) = LCC(I) - 1
   14   CONTINUE
C
C Update N and LNEW, and optimize the patch of triangles
C   containing K (on input) by applying swaps to the arcs
C   in IWK.
C
      N = NN
      LNEW = LNW
      IF (IWL .GT. 0) THEN
        NIT = 4*IWL
        CALL OPTIM (X,Y,IWL, LIST,LPTR,LEND,NIT,IWK, IERR)
        IF (IERR .NE. 0) GO TO 25
      ENDIF
C
C Successful termination.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   21 IER = 1
      RETURN
C
C Insufficient space reserved for IWK.
C
   22 IER = 2
      RETURN
C
C Invalid triangulation data structure.  NNB < 3 on input or
C   N2 is a neighbor of N1 but N1 is not a neighbor of N2.
C
   23 IER = 3
      RETURN
C
C K is an interior node with 4 or more neighbors, but the
C   number of neighbors could not be reduced.
C
   24 IER = 4
      RETURN
C
C Error flag returned by OPTIM.
C
   25 IER = 5
      WRITE (*,100) NIT, IERR
      RETURN
  100 FORMAT (//5X,'*** Error in OPTIM:  NIT = ',I4,
     .        ', IER = ',I1,' ***'/)
      END
      SUBROUTINE EDGE (IN1,IN2,X,Y, LWK,IWK,LIST,LPTR,
     .                 LEND, IER)
      INTEGER IN1, IN2, LWK, IWK(2,*), LIST(*), LPTR(*),
     .        LEND(*), IER
      REAL    X(*), Y(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/23/98
C
C   Given a triangulation of N nodes and a pair of nodal
C indexes IN1 and IN2, this routine swaps arcs as necessary
C to force IN1 and IN2 to be adjacent.  Only arcs which
C intersect IN1-IN2 are swapped out.  If a Delaunay triangu-
C lation is input, the resulting triangulation is as close
C as possible to a Delaunay triangulation in the sense that
C all arcs other than IN1-IN2 are locally optimal.
C
C   A sequence of calls to EDGE may be used to force the
C presence of a set of edges defining the boundary of a non-
C convex and/or multiply connected region (refer to Subrou-
C tine ADDCST), or to introduce barriers into the triangula-
C tion.  Note that Subroutine GETNP will not necessarily
C return closest nodes if the triangulation has been con-
C strained by a call to EDGE.  However, this is appropriate
C in some applications, such as triangle-based interpolation
C on a nonconvex domain.
C
C
C On input:
C
C       IN1,IN2 = Indexes (of X and Y) in the range 1 to N
C                 defining a pair of nodes to be connected
C                 by an arc.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.
C
C The above parameters are not altered by this routine.
C
C       LWK = Number of columns reserved for IWK.  This must
C             be at least NI -- the number of arcs which
C             intersect IN1-IN2.  (NI is bounded by N-3.)
C
C       IWK = Integer work array of length at least 2*LWK.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C On output:
C
C       LWK = Number of arcs which intersect IN1-IN2 (but
C             not more than the input value of LWK) unless
C             IER = 1 or IER = 3.  LWK = 0 if and only if
C             IN1 and IN2 were adjacent (or LWK=0) on input.
C
C       IWK = Array containing the indexes of the endpoints
C             of the new arcs other than IN1-IN2 unless IER
C             .GT. 0 or LWK = 0.  New arcs to the left of
C             IN2-IN1 are stored in the first K-1 columns
C             (left portion of IWK), column K contains
C             zeros, and new arcs to the right of IN2-IN1
C             occupy columns K+1,...,LWK.  (K can be deter-
C             mined by searching IWK for the zeros.)
C
C       LIST,LPTR,LEND = Data structure updated if necessary
C                        to reflect the presence of an arc
C                        connecting IN1 and IN2 unless IER
C                        .NE. 0.  The data structure has
C                        been altered if IER = 4.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if IN1 .LT. 1, IN2 .LT. 1, IN1 = IN2,
C                     or LWK .LT. 0 on input.
C             IER = 2 if more space is required in IWK.
C             IER = 3 if IN1 and IN2 could not be connected
C                     due to either an invalid data struc-
C                     ture or collinear nodes (and floating
C                     point error).
C             IER = 4 if an error flag was returned by
C                     OPTIM.
C
C   An error message is written to the standard output unit
C in the case of IER = 3 or IER = 4.
C
C Modules required by EDGE:  LEFT, LSTPTR, OPTIM, SWAP,
C                              SWPTST
C
C Intrinsic function called by EDGE:  ABS
C
C***********************************************************
C
      LOGICAL LEFT
      INTEGER I, IERR, IWC, IWCP1, IWEND, IWF, IWL, LFT, LP,
     .        LPL, LP21, NEXT, NIT, NL, NR, N0, N1, N2,
     .        N1FRST, N1LST
      REAL    DX, DY, X0, Y0, X1, Y1, X2, Y2
C
C Local parameters:
C
C DX,DY =   Components of arc N1-N2
C I =       DO-loop index and column index for IWK
C IERR =    Error flag returned by Subroutine OPTIM
C IWC =     IWK index between IWF and IWL -- NL->NR is
C             stored in IWK(1,IWC)->IWK(2,IWC)
C IWCP1 =   IWC + 1
C IWEND =   Input or output value of LWK
C IWF =     IWK (column) index of the first (leftmost) arc
C             which intersects IN1->IN2
C IWL =     IWK (column) index of the last (rightmost) are
C             which intersects IN1->IN2
C LFT =     Flag used to determine if a swap results in the
C             new arc intersecting IN1-IN2 -- LFT = 0 iff
C             N0 = IN1, LFT = -1 implies N0 LEFT IN1->IN2,
C             and LFT = 1 implies N0 LEFT IN2->IN1
C LP21 =    Unused parameter returned by SWAP
C LP =      List pointer (index) for LIST and LPTR
C LPL =     Pointer to the last neighbor of IN1 or NL
C N0 =      Neighbor of N1 or node opposite NR->NL
C N1,N2 =   Local copies of IN1 and IN2
C N1FRST =  First neighbor of IN1
C N1LST =   (Signed) last neighbor of IN1
C NEXT =    Node opposite NL->NR
C NIT =     Flag or number of iterations employed by OPTIM
C NL,NR =   Endpoints of an arc which intersects IN1-IN2
C             with NL LEFT IN1->IN2
C X0,Y0 =   Coordinates of N0
C X1,Y1 =   Coordinates of IN1
C X2,Y2 =   Coordinates of IN2
C
C
C Store IN1, IN2, and LWK in local variables and test for
C   errors.
C
      N1 = IN1
      N2 = IN2
      IWEND = LWK
      IF (N1 .LT. 1  .OR.  N2 .LT. 1  .OR.  N1 .EQ. N2  .OR.
     .    IWEND .LT. 0) GO TO 31
C
C Test for N2 as a neighbor of N1.  LPL points to the last
C   neighbor of N1.
C
      LPL = LEND(N1)
      N0 = ABS(LIST(LPL))
      LP = LPL
    1 IF (N0 .EQ. N2) GO TO 30
        LP = LPTR(LP)
        N0 = LIST(LP)
        IF (LP .NE. LPL) GO TO 1
C
C Initialize parameters.
C
      IWL = 0
      NIT = 0
C
C Store the coordinates of N1 and N2.
C
    2 X1 = X(N1)
      Y1 = Y(N1)
      X2 = X(N2)
      Y2 = Y(N2)
C
C Set NR and NL to adjacent neighbors of N1 such that
C   NR LEFT N2->N1 and NL LEFT N1->N2,
C   (NR Forward N1->N2 or NL Forward N1->N2), and
C   (NR Forward N2->N1 or NL Forward N2->N1).
C
C   Initialization:  Set N1FRST and N1LST to the first and
C     (signed) last neighbors of N1, respectively, and
C     initialize NL to N1FRST.
C
      LPL = LEND(N1)
      N1LST = LIST(LPL)
      LP = LPTR(LPL)
      N1FRST = LIST(LP)
      NL = N1FRST
      IF (N1LST .LT. 0) GO TO 4
C
C   N1 is an interior node.  Set NL to the first candidate
C     for NR (NL LEFT N2->N1).
C
    3 IF ( LEFT(X2,Y2,X1,Y1,X(NL),Y(NL)) ) GO TO 4
        LP = LPTR(LP)
        NL = LIST(LP)
        IF (NL .NE. N1FRST) GO TO 3
C
C   All neighbors of N1 are strictly left of N1->N2.
C
      GO TO 5
C
C   NL = LIST(LP) LEFT N2->N1.  Set NR to NL and NL to the
C     following neighbor of N1.
C
    4 NR = NL
        LP = LPTR(LP)
        NL = ABS(LIST(LP))
        IF ( LEFT(X1,Y1,X2,Y2,X(NL),Y(NL)) ) THEN
C
C   NL LEFT N1->N2 and NR LEFT N2->N1.  The Forward tests
C     are employed to avoid an error associated with
C     collinear nodes.
C
          DX = X2-X1
          DY = Y2-Y1
          IF ((DX*(X(NL)-X1)+DY*(Y(NL)-Y1) .GE. 0.  .OR.
     .         DX*(X(NR)-X1)+DY*(Y(NR)-Y1) .GE. 0.)  .AND.
     .        (DX*(X(NL)-X2)+DY*(Y(NL)-Y2) .LE. 0.  .OR.
     .         DX*(X(NR)-X2)+DY*(Y(NR)-Y2) .LE. 0.)) GO TO 6
C
C   NL-NR does not intersect N1-N2.  However, there is
C     another candidate for the first arc if NL lies on
C     the line N1-N2.
C
          IF ( .NOT. LEFT(X2,Y2,X1,Y1,X(NL),Y(NL)) ) GO TO 5
        ENDIF
C
C   Bottom of loop.
C
        IF (NL .NE. N1FRST) GO TO 4
C
C Either the triangulation is invalid or N1-N2 lies on the
C   convex hull boundary and an edge NR->NL (opposite N1 and
C   intersecting N1-N2) was not found due to floating point
C   error.  Try interchanging N1 and N2 -- NIT > 0 iff this
C   has already been done.
C
    5 IF (NIT .GT. 0) GO TO 33
      NIT = 1
      N1 = N2
      N2 = IN1
      GO TO 2
C
C Store the ordered sequence of intersecting edges NL->NR in
C   IWK(1,IWL)->IWK(2,IWL).
C
    6 IWL = IWL + 1
      IF (IWL .GT. IWEND) GO TO 32
      IWK(1,IWL) = NL
      IWK(2,IWL) = NR
C
C   Set NEXT to the neighbor of NL which follows NR.
C
      LPL = LEND(NL)
      LP = LPTR(LPL)
C
C   Find NR as a neighbor of NL.  The search begins with
C     the first neighbor.
C
    7 IF (LIST(LP) .EQ. NR) GO TO 8
        LP = LPTR(LP)
        IF (LP .NE. LPL) GO TO 7
C
C   NR must be the last neighbor, and NL->NR cannot be a
C     boundary edge.
C
      IF (LIST(LP) .NE. NR) GO TO 33
C
C   Set NEXT to the neighbor following NR, and test for
C     termination of the store loop.
C
    8 LP = LPTR(LP)
      NEXT = ABS(LIST(LP))
      IF (NEXT .EQ. N2) GO TO 9
C
C   Set NL or NR to NEXT.
C
      IF ( LEFT(X1,Y1,X2,Y2,X(NEXT),Y(NEXT)) ) THEN
        NL = NEXT
      ELSE
        NR = NEXT
      ENDIF
      GO TO 6
C
C IWL is the number of arcs which intersect N1-N2.
C   Store LWK.
C
    9 LWK = IWL
      IWEND = IWL
C
C Initialize for edge swapping loop -- all possible swaps
C   are applied (even if the new arc again intersects
C   N1-N2), arcs to the left of N1->N2 are stored in the
C   left portion of IWK, and arcs to the right are stored in
C   the right portion.  IWF and IWL index the first and last
C   intersecting arcs.
C
      IWF = 1
C
C Top of loop -- set N0 to N1 and NL->NR to the first edge.
C   IWC points to the arc currently being processed.  LFT
C   .LE. 0 iff N0 LEFT N1->N2.
C
   10 LFT = 0
      N0 = N1
      X0 = X1
      Y0 = Y1
      NL = IWK(1,IWF)
      NR = IWK(2,IWF)
      IWC = IWF
C
C   Set NEXT to the node opposite NL->NR unless IWC is the
C     last arc.
C
   11 IF (IWC .EQ. IWL) GO TO 21
      IWCP1 = IWC + 1
      NEXT = IWK(1,IWCP1)
      IF (NEXT .NE. NL) GO TO 16
      NEXT = IWK(2,IWCP1)
C
C   NEXT RIGHT N1->N2 and IWC .LT. IWL.  Test for a possible
C     swap.
C
      IF ( .NOT. LEFT(X0,Y0,X(NR),Y(NR),X(NEXT),Y(NEXT)) )
     .   GO TO 14
      IF (LFT .GE. 0) GO TO 12
      IF ( .NOT. LEFT(X(NL),Y(NL),X0,Y0,X(NEXT),Y(NEXT)) )
     .   GO TO 14
C
C   Replace NL->NR with N0->NEXT.
C
      CALL SWAP (NEXT,N0,NL,NR, LIST,LPTR,LEND, LP21)
      IWK(1,IWC) = N0
      IWK(2,IWC) = NEXT
      GO TO 15
C
C   Swap NL-NR for N0-NEXT, shift columns IWC+1,...,IWL to
C     the left, and store N0-NEXT in the right portion of
C     IWK.
C
   12 CALL SWAP (NEXT,N0,NL,NR, LIST,LPTR,LEND, LP21)
      DO 13 I = IWCP1,IWL
        IWK(1,I-1) = IWK(1,I)
        IWK(2,I-1) = IWK(2,I)
   13   CONTINUE
      IWK(1,IWL) = N0
      IWK(2,IWL) = NEXT
      IWL = IWL - 1
      NR = NEXT
      GO TO 11
C
C   A swap is not possible.  Set N0 to NR.
C
   14 N0 = NR
      X0 = X(N0)
      Y0 = Y(N0)
      LFT = 1
C
C   Advance to the next arc.
C
   15 NR = NEXT
      IWC = IWC + 1
      GO TO 11
C
C   NEXT LEFT N1->N2, NEXT .NE. N2, and IWC .LT. IWL.
C     Test for a possible swap.
C
   16 IF ( .NOT. LEFT(X(NL),Y(NL),X0,Y0,X(NEXT),Y(NEXT)) )
     .   GO TO 19
      IF (LFT .LE. 0) GO TO 17
      IF ( .NOT. LEFT(X0,Y0,X(NR),Y(NR),X(NEXT),Y(NEXT)) )
     .   GO TO 19
C
C   Replace NL->NR with NEXT->N0.
C
      CALL SWAP (NEXT,N0,NL,NR, LIST,LPTR,LEND, LP21)
      IWK(1,IWC) = NEXT
      IWK(2,IWC) = N0
      GO TO 20
C
C   Swap NL-NR for N0-NEXT, shift columns IWF,...,IWC-1 to
C     the right, and store N0-NEXT in the left portion of
C     IWK.
C
   17 CALL SWAP (NEXT,N0,NL,NR, LIST,LPTR,LEND, LP21)
      DO 18 I = IWC-1,IWF,-1
        IWK(1,I+1) = IWK(1,I)
        IWK(2,I+1) = IWK(2,I)
   18   CONTINUE
      IWK(1,IWF) = N0
      IWK(2,IWF) = NEXT
      IWF = IWF + 1
      GO TO 20
C
C   A swap is not possible.  Set N0 to NL.
C
   19 N0 = NL
      X0 = X(N0)
      Y0 = Y(N0)
      LFT = -1
C
C   Advance to the next arc.
C
   20 NL = NEXT
      IWC = IWC + 1
      GO TO 11
C
C   N2 is opposite NL->NR (IWC = IWL).
C
   21 IF (N0 .EQ. N1) GO TO 24
      IF (LFT .LT. 0) GO TO 22
C
C   N0 RIGHT N1->N2.  Test for a possible swap.
C
      IF ( .NOT. LEFT(X0,Y0,X(NR),Y(NR),X2,Y2) ) GO TO 10
C
C   Swap NL-NR for N0-N2 and store N0-N2 in the right
C     portion of IWK.
C
      CALL SWAP (N2,N0,NL,NR, LIST,LPTR,LEND, LP21)
      IWK(1,IWL) = N0
      IWK(2,IWL) = N2
      IWL = IWL - 1
      GO TO 10
C
C   N0 LEFT N1->N2.  Test for a possible swap.
C
   22 IF ( .NOT. LEFT(X(NL),Y(NL),X0,Y0,X2,Y2) ) GO TO 10
C
C   Swap NL-NR for N0-N2, shift columns IWF,...,IWL-1 to the
C     right, and store N0-N2 in the left portion of IWK.
C
      CALL SWAP (N2,N0,NL,NR, LIST,LPTR,LEND, LP21)
      I = IWL
   23 IWK(1,I) = IWK(1,I-1)
      IWK(2,I) = IWK(2,I-1)
      I = I - 1
      IF (I .GT. IWF) GO TO 23
      IWK(1,IWF) = N0
      IWK(2,IWF) = N2
      IWF = IWF + 1
      GO TO 10
C
C IWF = IWC = IWL.  Swap out the last arc for N1-N2 and
C   store zeros in IWK.
C
   24 CALL SWAP (N2,N1,NL,NR, LIST,LPTR,LEND, LP21)
      IWK(1,IWC) = 0
      IWK(2,IWC) = 0
C
C Optimization procedure --
C
      IF (IWC .GT. 1) THEN
C
C   Optimize the set of new arcs to the left of IN1->IN2.
C
        NIT = 3*(IWC-1)
        CALL OPTIM (X,Y,IWC-1, LIST,LPTR,LEND,NIT,IWK, IERR)
        IF (IERR .NE. 0) GO TO 34
      ENDIF
      IF (IWC .LT. IWEND) THEN
C
C   Optimize the set of new arcs to the right of IN1->IN2.
C
        NIT = 3*(IWEND-IWC)
        CALL OPTIM (X,Y,IWEND-IWC, LIST,LPTR,LEND,NIT,
     .              IWK(1,IWC+1), IERR)
        IF (IERR .NE. 0) GO TO 34
      ENDIF
C
C Successful termination.
C
      IER = 0
      RETURN
C
C IN1 and IN2 were adjacent on input.
C
   30 IER = 0
      RETURN
C
C Invalid input parameter.
C
   31 IER = 1
      RETURN
C
C Insufficient space reserved for IWK.
C
   32 IER = 2
      RETURN
C
C Invalid triangulation data structure or collinear nodes
C   on convex hull boundary.
C
   33 IER = 3
      WRITE (*,130) IN1, IN2
  130 FORMAT (//5X,'*** Error in EDGE:  Invalid triangula',
     .        'tion or null triangles on boundary'/
     .        9X,'IN1 =',I4,', IN2=',I4/)
      RETURN
C
C Error flag returned by OPTIM.
C
   34 IER = 4
      WRITE (*,140) NIT, IERR
  140 FORMAT (//5X,'*** Error in OPTIM:  NIT = ',I4,
     .        ', IER = ',I1,' ***'/)
      RETURN
      END
      SUBROUTINE GETNP (NCC,LCC,N,X,Y,LIST,LPTR,LEND,
     .                  L, NPTS,DS, IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        L, NPTS(L), IER
      REAL    X(N), Y(N), DS(L)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/12/94
C
C   Given a triangulation of N nodes and an array NPTS con-
C taining the indexes of L-1 nodes ordered by distance from
C NPTS(1), this subroutine sets NPTS(L) to the index of the
C next node in the sequence -- the node, other than NPTS(1),
C ...,NPTS(L-1), which is closest to NPTS(1).  Thus, the
C ordered sequence of K closest nodes to N1 (including N1)
C may be determined by K-1 calls to GETNP with NPTS(1) = N1
C and L = 2,3,...,K for K .GE. 2.  Note that NPTS must in-
C clude constraint nodes as well as non-constraint nodes.
C Thus, a sequence of K1 closest non-constraint nodes to N1
C must be obtained as a subset of the closest K2 nodes to N1
C for some K2 .GE. K1.
C
C   The terms closest and distance have special definitions
C when constraint nodes are present in the triangulation.
C Nodes N1 and N2 are said to be visible from each other if
C and only if the line segment N1-N2 intersects no con-
C straint arc (except possibly itself) and is not an interi-
C or constraint arc (arc whose interior lies in a constraint
C region).  A path from N1 to N2 is an ordered sequence of
C nodes, with N1 first and N2 last, such that adjacent path
C elements are visible from each other.  The path length is
C the sum of the Euclidean distances between adjacent path
C nodes.  Finally, the distance from N1 to N2 is defined to
C be the length of the shortest path from N1 to N2.
C
C   The algorithm uses the property of a Delaunay triangula-
C tion that the K-th closest node to N1 is a neighbor of one
C of the K-1 closest nodes to N1.  With the definition of
C distance used here, this property holds when constraints
C are present as long as non-constraint arcs are locally
C optimal.
C
C
C On input:
C
C       NCC = Number of constraints.  NCC .GE. 0.
C
C       LCC = List of constraint curve starting indexes (or
C             dummy array of length 1 if NCC = 0).  Refer to
C             Subroutine ADDCST.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations if NCC > 0.
C
C       LIST,LPTR,LEND = Triangulation data structure.  Re-
C                        fer to Subroutine TRMESH.
C
C       L = Number of nodes in the sequence on output.  2
C           .LE. L .LE. N.
C
C       NPTS = Array of length .GE. L containing the indexes
C              of the L-1 closest nodes to NPTS(1) in the
C              first L-1 locations.
C
C       DS = Array of length .GE. L containing the distance
C            (defined above) between NPTS(1) and NPTS(I) in
C            the I-th position for I = 1,...,L-1.  Thus,
C            DS(1) = 0.
C
C Input parameters other than NPTS(L) and DS(L) are not
C   altered by this routine.
C
C On output:
C
C       NPTS = Array updated with the index of the L-th
C              closest node to NPTS(1) in position L unless
C              IER .NE. 0.
C
C       DS = Array updated with the distance between NPTS(1)
C            and NPTS(L) in position L unless IER .NE. 0.
C
C       IER = Error indicator:
C             IER =  0 if no errors were encountered.
C             IER = -1 if NCC, N, L, or an LCC entry is
C                      outside its valid range on input.
C             IER =  K if NPTS(K) is not a valid index in
C                      the range 1 to N.
C
C Module required by GETNP:  INTSEC
C
C Intrinsic functions called by GETNP:  ABS, MIN, SQRT
C
C***********************************************************
C
      LOGICAL INTSEC
      INTEGER I, IFRST, ILAST, J, K, KM1, LCC1, LM1, LP,
     .        LPCL, LPK, LPKL, N1, NC, NF1, NF2, NJ, NK,
     .        NKBAK, NKFOR, NL, NN
      LOGICAL ISW, VIS, NCF, NJF, SKIP, SKSAV, LFT1, LFT2,
     .        LFT12
      REAL    DC, DL, X1, XC, XJ, XK, Y1, YC, YJ, YK
C
C Store parameters in local variables and test for errors.
C   LCC1 indexes the first constraint node.
C
      IER = -1
      NN = N
      LCC1 = NN+1
      LM1 = L-1
      IF (NCC .LT. 0  .OR.  LM1 .LT. 1  .OR.  LM1 .GE. NN)
     .   RETURN
      IF (NCC .EQ. 0) THEN
        IF (NN .LT. 3) RETURN
      ELSE
        DO 1 I = NCC,1,-1
          IF (LCC1 - LCC(I) .LT. 3) RETURN
          LCC1 = LCC(I)
    1     CONTINUE
        IF (LCC1 .LT. 1) RETURN
      ENDIF
C
C Test for an invalid index in NPTS.
C
      DO 2 K = 1,LM1
        NK = NPTS(K)
        IF (NK .LT. 1  .OR.  NK .GT. NN) THEN
          IER = K
          RETURN
        ENDIF
    2   CONTINUE
C
C Store N1 = NPTS(1) and mark the elements of NPTS.
C
      N1 = NPTS(1)
      X1 = X(N1)
      Y1 = Y(N1)
      DO 3 K = 1,LM1
        NK = NPTS(K)
        LEND(NK) = -LEND(NK)
    3   CONTINUE
C
C Candidates NC for NL = NPTS(L) are the unmarked visible
C   neighbors of nodes NK in NPTS.  ISW is an initialization
C   switch set to .TRUE. when NL and its distance DL from N1
C   have been initialized with the first candidate encount-
C   ered.
C
      ISW = .FALSE.
      DL = 0.
C
C Loop on marked nodes NK = NPTS(K).  LPKL indexes the last
C   neighbor of NK in LIST.
C
      DO 16 K = 1,LM1
        KM1 = K - 1
        NK = NPTS(K)
        XK = X(NK)
        YK = Y(NK)
        LPKL = -LEND(NK)
        NKFOR = 0
        NKBAK = 0
        VIS = .TRUE.
        IF (NK .GE. LCC1) THEN
C
C   NK is a constraint node.  Set NKFOR and NKBAK to the
C     constraint nodes which follow and precede NK.  IFRST
C     and ILAST are set to the first and last nodes in the
C     constraint containing NK.
C
          IFRST = NN + 1
          DO 4 I = NCC,1,-1
            ILAST = IFRST - 1
            IFRST = LCC(I)
            IF (NK .GE. IFRST) GO TO 5
    4       CONTINUE
C
    5     IF (NK .LT. ILAST) THEN
            NKFOR = NK + 1
          ELSE
            NKFOR = IFRST
          ENDIF
          IF (NK .GT. IFRST) THEN
            NKBAK = NK - 1
          ELSE
            NKBAK = ILAST
          ENDIF
C
C   Initialize VIS to TRUE iff NKFOR precedes NKBAK in the
C     adjacency list for NK -- the first neighbor is visi-
C     ble and is not NKBAK.
C
          LPK = LPKL
    6     LPK = LPTR(LPK)
            NC = ABS(LIST(LPK))
            IF (NC .NE. NKFOR  .AND.  NC .NE. NKBAK) GO TO 6
          VIS = NC .EQ. NKFOR
        ENDIF
C
C Loop on neighbors NC of NK, bypassing marked and nonvis-
C   ible neighbors.
C
        LPK = LPKL
    7   LPK = LPTR(LPK)
          NC = ABS(LIST(LPK))
          IF (NC .EQ. NKBAK) VIS = .TRUE.
C
C   VIS = .FALSE. iff NK-NC is an interior constraint arc
C     (NK is a constraint node and NC lies strictly between
C     NKFOR and NKBAK).
C
          IF (.NOT. VIS) GO TO 15
          IF (NC .EQ. NKFOR) VIS = .FALSE.
          IF (LEND(NC) .LT. 0) GO TO 15
C
C Initialize distance DC between N1 and NC to Euclidean
C   distance.
C
          XC = X(NC)
          YC = Y(NC)
          DC = SQRT((XC-X1)*(XC-X1) + (YC-Y1)*(YC-Y1))
          IF (ISW  .AND.  DC .GE. DL) GO TO 15
          IF (K .EQ. 1) GO TO 14
C
C K .GE. 2.  Store the pointer LPCL to the last neighbor
C   of NC.
C
          LPCL = LEND(NC)
C
C Set DC to the length of the shortest path from N1 to NC
C   which has not previously been encountered and which is
C   a viable candidate for the shortest path from N1 to NL.
C   This is Euclidean distance iff NC is visible from N1.
C   Since the shortest path from N1 to NL contains only ele-
C   ments of NPTS which are constraint nodes (in addition to
C   N1 and NL), only these need be considered for the path
C   from N1 to NC.  Thus, for distance function D(A,B) and
C   J = 1,...,K, DC = min(D(N1,NJ) + D(NJ,NC)) over con-
C   straint nodes NJ = NPTS(J) which are visible from NC.
C
          DO 13 J = 1,KM1
            NJ = NPTS(J)
            IF (J .GT. 1  .AND.  NJ .LT. LCC1) GO TO 13
C
C If NC is a visible neighbor of NJ, a path from N1 to NC
C   containing NJ has already been considered.  Thus, NJ may
C   be bypassed if it is adjacent to NC.
C
            LP = LPCL
    8       LP = LPTR(LP)
              IF ( NJ .EQ. ABS(LIST(LP)) ) GO TO 12
              IF (LP .NE. LPCL) GO TO 8
C
C NJ is a constraint node (unless J=1) not adjacent to NC,
C   and is visible from NC iff NJ-NC is not intersected by
C   a constraint arc.  Loop on constraints I in reverse
C   order --
C
            XJ = X(NJ)
            YJ = Y(NJ)
            IFRST = NN+1
            DO 11 I = NCC,1,-1
              ILAST = IFRST - 1
              IFRST = LCC(I)
              NF1 = ILAST
              NCF = NF1 .EQ. NC
              NJF = NF1 .EQ. NJ
              SKIP = NCF  .OR.  NJF
C
C Loop on boundary constraint arcs NF1-NF2 which contain
C   neither NC nor NJ.  NCF and NJF are TRUE iff NC (or NJ)
C   has been encountered in the constraint, and SKIP =
C   .TRUE. iff NF1 = NC or NF1 = NJ.
C
              DO 10 NF2 = IFRST,ILAST
                IF (NF2 .EQ. NC) NCF = .TRUE.
                IF (NF2 .EQ. NJ) NJF = .TRUE.
                SKSAV = SKIP
                SKIP = NF2 .EQ. NC  .OR.  NF2 .EQ. NJ
C
C   The last constraint arc in the constraint need not be
C     tested if none of the arcs have been skipped.
C
                IF ( SKSAV  .OR.  SKIP  .OR.
     .               (NF2 .EQ. ILAST  .AND.
     .               .NOT. NCF  .AND.  .NOT. NJF) ) GO TO 9
                IF ( INTSEC(X(NF1),Y(NF1),X(NF2),Y(NF2),
     .                      XC,YC,XJ,YJ) ) GO TO 12
    9           NF1 = NF2
   10           CONTINUE
              IF (.NOT. NCF  .OR.  .NOT. NJF) GO TO 11
C
C NC and NJ are constraint nodes in the same constraint.
C   NC-NJ is intersected by an interior constraint arc iff
C   1)  NC LEFT NF2->NF1 and (NJ LEFT NF1->NC and NJ LEFT
C         NC->NF2) or
C   2)  NC .NOT. LEFT NF2->NF1 and (NJ LEFT NF1->NC or
C         NJ LEFT NC->NF2),
C   where NF1, NC, NF2 are consecutive constraint nodes.
C
              IF (NC .NE. IFRST) THEN
                NF1 = NC - 1
              ELSE
                NF1 = ILAST
              ENDIF
              IF (NC .NE. ILAST) THEN
                NF2 = NC + 1
              ELSE
                NF2 = IFRST
              ENDIF
              LFT1 = (XC-X(NF1))*(YJ-Y(NF1)) .GE.
     .               (XJ-X(NF1))*(YC-Y(NF1))
              LFT2 = (X(NF2)-XC)*(YJ-YC) .GE.
     .               (XJ-XC)*(Y(NF2)-YC)
              LFT12 = (X(NF1)-X(NF2))*(YC-Y(NF2)) .GE.
     .                (XC-X(NF2))*(Y(NF1)-Y(NF2))
              IF ( (LFT1  .AND.  LFT2)  .OR.  (.NOT. LFT12
     .             .AND.  (LFT1  .OR.  LFT2)) ) GO TO 12
   11         CONTINUE
C
C NJ is visible from NC.  Exit the loop with DC = Euclidean
C   distance if J = 1.
C
            IF (J .EQ. 1) GO TO 14
            DC = MIN(DC,DS(J) + SQRT((XC-XJ)*(XC-XJ) +
     .                  (YC-YJ)*(YC-YJ)))
            GO TO 13
C
C NJ is not visible from NC or is adjacent to NC.  Initial-
C   ize DC with D(N1,NK) + D(NK,NC) if J = 1.
C
   12       IF (J .EQ. 1) DC = DS(K) + SQRT((XC-XK)*(XC-XK)
     .                         + (YC-YK)*(YC-YK))
   13       CONTINUE
C
C Compare DC with DL.
C
          IF (ISW  .AND.  DC .GE. DL) GO TO 15
C
C The first (or a closer) candidate for NL has been
C   encountered.
C
   14     NL = NC
          DL = DC
          ISW = .TRUE.
   15     IF (LPK .NE. LPKL) GO TO 7
   16   CONTINUE
C
C Unmark the elements of NPTS and store NL and DL.
C
      DO 17 K = 1,LM1
        NK = NPTS(K)
        LEND(NK) = -LEND(NK)
   17   CONTINUE
      NPTS(L) = NL
      DS(L) = DL
      IER = 0
      RETURN
      END
      INTEGER FUNCTION INDXCC (NCC,LCC,N,LIST,LEND)
      INTEGER NCC, LCC(*), N, LIST(*), LEND(N)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   08/25/91
C
C   Given a constrained Delaunay triangulation, this func-
C tion returns the index, if any, of an exterior constraint
C curve (an unbounded constraint region).  An exterior con-
C straint curve is assumed to be present if and only if the
C clockwise-ordered sequence of boundary nodes is a subse-
C quence of a constraint node sequence.  The triangulation
C adjacencies corresponding to constraint edges may or may
C not have been forced by a call to ADDCST, and the con-
C straint region may or may not be valid (contain no nodes).
C
C
C On input:
C
C       NCC = Number of constraints.  NCC .GE. 0.
C
C       LCC = List of constraint curve starting indexes (or
C             dummy array of length 1 if NCC = 0).  Refer to
C             Subroutine ADDCST.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       LIST,LEND = Data structure defining the triangula-
C                   tion.  Refer to Subroutine TRMESH.
C
C   Input parameters are not altered by this function.  Note
C that the parameters are not tested for validity.
C
C On output:
C
C       INDXCC = Index of the exterior constraint curve, if
C                present, or 0 otherwise.
C
C Modules required by INDXCC:  None
C
C***********************************************************
C
      INTEGER I, IFRST, ILAST, LP, N0, NST, NXT
      INDXCC = 0
      IF (NCC .LT. 1) RETURN
C
C Set N0 to the boundary node with smallest index.
C
      N0 = 0
    1 N0 = N0 + 1
        LP = LEND(N0)
        IF (LIST(LP) .GT. 0) GO TO 1
C
C Search in reverse order for the constraint I, if any, that
C   contains N0.  IFRST and ILAST index the first and last
C   nodes in constraint I.
C
      I = NCC
      ILAST = N
    2 IFRST = LCC(I)
        IF (N0 .GE. IFRST) GO TO 3
        IF (I .EQ. 1) RETURN
        I = I - 1
        ILAST = IFRST - 1
        GO TO 2
C
C N0 is in constraint I which indexes an exterior constraint
C   curve iff the clockwise-ordered sequence of boundary
C   node indexes beginning with N0 is increasing and bounded
C   above by ILAST.
C
    3 NST = N0
C
    4 NXT = -LIST(LP)
        IF (NXT .EQ. NST) GO TO 5
        IF (NXT .LE. N0  .OR.  NXT .GT. ILAST) RETURN
        N0 = NXT
        LP = LEND(N0)
        GO TO 4
C
C Constraint I contains the boundary node sequence as a
C   subset.
C
    5 INDXCC = I
      RETURN
      END
      SUBROUTINE INSERT (K,LP, LIST,LPTR,LNEW )
      INTEGER K, LP, LIST(*), LPTR(*), LNEW
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This subroutine inserts K as a neighbor of N1 following
C N2, where LP is the LIST pointer of N2 as a neighbor of
C N1.  Note that, if N2 is the last neighbor of N1, K will
C become the first neighbor (even if N1 is a boundary node).
C
C
C On input:
C
C       K = Index of the node to be inserted.
C
C       LP = LIST pointer of N2 as a neighbor of N1.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LNEW = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C On output:
C
C       LIST,LPTR,LNEW = Data structure updated with the
C                        addition of node K.
C
C Modules required by INSERT:  None
C
C***********************************************************
C
      INTEGER LSAV
C
      LSAV = LPTR(LP)
      LPTR(LP) = LNEW
      LIST(LNEW) = K
      LPTR(LNEW) = LSAV
      LNEW = LNEW + 1
      RETURN
      END
      SUBROUTINE INTADD (KK,I1,I2,I3, LIST,LPTR,LEND,LNEW )
      INTEGER KK, I1, I2, I3, LIST(*), LPTR(*), LEND(*),
     .        LNEW
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/22/91
C
C   This subroutine adds an interior node to a triangulation
C of a set of points in the plane.  The data structure is
C updated with the insertion of node KK into the triangle
C whose vertices are I1, I2, and I3.  No optimization of the
C triangulation is performed.
C
C
C On input:
C
C       KK = Index of the node to be inserted.  KK .GE. 1
C            and KK must not be equal to I1, I2, or I3.
C
C       I1,I2,I3 = Indexes of the counterclockwise-ordered
C                  sequence of vertices of a triangle which
C                  contains node KK.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.  Refer to Sub-
C                             routine TRMESH.  Triangle
C                             (I1,I2,I3) must be included
C                             in the triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the addition of node KK.  KK
C                             will be connected to nodes I1,
C                             I2, and I3.
C
C Modules required by INTADD:  INSERT, LSTPTR
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER K, LP, N1, N2, N3
      K = KK
C
C Initialization.
C
      N1 = I1
      N2 = I2
      N3 = I3
C
C Add K as a neighbor of I1, I2, and I3.
C
      LP = LSTPTR(LEND(N1),N2,LIST,LPTR)
      CALL INSERT (K,LP,LIST,LPTR,LNEW)
      LP = LSTPTR(LEND(N2),N3,LIST,LPTR)
      CALL INSERT (K,LP,LIST,LPTR,LNEW)
      LP = LSTPTR(LEND(N3),N1,LIST,LPTR)
      CALL INSERT (K,LP,LIST,LPTR,LNEW)
C
C Add I1, I2, and I3 as neighbors of K.
C
      LIST(LNEW) = N1
      LIST(LNEW+1) = N2
      LIST(LNEW+2) = N3
      LPTR(LNEW) = LNEW + 1
      LPTR(LNEW+1) = LNEW + 2
      LPTR(LNEW+2) = LNEW
      LEND(K) = LNEW + 2
      LNEW = LNEW + 3
      RETURN
      END
      LOGICAL FUNCTION INTSEC (X1,Y1,X2,Y2,X3,Y3,X4,Y4)
      REAL X1, Y1, X2, Y2, X3, Y3, X4, Y4
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   Given a pair of line segments P1-P2 and P3-P4, this
C function returns the value .TRUE. if and only if P1-P2
C shares one or more points with P3-P4.  The line segments
C include their endpoints, and the four points need not be
C distinct.  Thus, either line segment may consist of a
C single point, and the segments may meet in a V (which is
C treated as an intersection).  Note that an incorrect
C decision may result from floating point error if the four
C endpoints are nearly collinear.
C
C
C On input:
C
C       X1,Y1 = Coordinates of P1.
C
C       X2,Y2 = Coordinates of P2.
C
C       X3,Y3 = Coordinates of P3.
C
C       X4,Y4 = Coordinates of P4.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       INTSEC = Logical value defined above.
C
C Modules required by INTSEC:  None
C
C***********************************************************
C
      REAL A, B, D, DX12, DX31, DX34, DY12, DY31, DY34
C
C Test for overlap between the smallest rectangles that
C   contain the line segments and have sides parallel to
C   the axes.
C
      IF ((X1 .LT. X3  .AND.  X1 .LT. X4  .AND.  X2 .LT. X3
     .     .AND.  X2 .LT. X4)  .OR.
     .    (X1 .GT. X3  .AND.  X1 .GT. X4  .AND.  X2 .GT. X3
     .     .AND.  X2 .GT. X4)  .OR.
     .    (Y1 .LT. Y3  .AND.  Y1 .LT. Y4  .AND.  Y2 .LT. Y3
     .     .AND.  Y2 .LT. Y4)  .OR.
     .    (Y1 .GT. Y3  .AND.  Y1 .GT. Y4  .AND.  Y2 .GT. Y3
     .     .AND.  Y2 .GT. Y4)) THEN
        INTSEC = .FALSE.
        RETURN
      ENDIF
C
C Compute A = P4-P3 X P1-P3, B = P2-P1 X P1-P3, and
C   D = P2-P1 X P4-P3 (Z components).
C
      DX12 = X2 - X1
      DY12 = Y2 - Y1
      DX34 = X4 - X3
      DY34 = Y4 - Y3
      DX31 = X1 - X3
      DY31 = Y1 - Y3
      A = DX34*DY31 - DX31*DY34
      B = DX12*DY31 - DX31*DY12
      D = DX12*DY34 - DX34*DY12
      IF (D .EQ. 0.) GO TO 1
C
C D .NE. 0 and the point of intersection of the lines de-
C   fined by the line segments is P = P1 + (A/D)*(P2-P1) =
C   P3 + (B/D)*(P4-P3).
C
      INTSEC = A/D .GE. 0.  .AND.  A/D .LE. 1.  .AND.
     .         B/D .GE. 0.  .AND.  B/D .LE. 1.
      RETURN
C
C D .EQ. 0 and thus either the line segments are parallel,
C   or one (or both) of them is a single point.
C
    1 INTSEC = A .EQ. 0.  .AND.  B .EQ. 0.
      RETURN
      END
      INTEGER FUNCTION JRAND (N, IX,IY,IZ )
      INTEGER N, IX, IY, IZ
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/28/98
C
C   This function returns a uniformly distributed pseudo-
C random integer in the range 1 to N.
C
C
C On input:
C
C       N = Maximum value to be returned.
C
C N is not altered by this function.
C
C       IX,IY,IZ = Integer seeds initialized to values in
C                  the range 1 to 30,000 before the first
C                  call to JRAND, and not altered between
C                  subsequent calls (unless a sequence of
C                  random numbers is to be repeated by
C                  reinitializing the seeds).
C
C On output:
C
C       IX,IY,IZ = Updated integer seeds.
C
C       JRAND = Random integer in the range 1 to N.
C
C Reference:  B. A. Wichmann and I. D. Hill, "An Efficient
C             and Portable Pseudo-random Number Generator",
C             Applied Statistics, Vol. 31, No. 2, 1982,
C             pp. 188-190.
C
C Modules required by JRAND:  None
C
C Intrinsic functions called by JRAND:  INT, MOD, REAL
C
C***********************************************************
C
      REAL U, X
C
C Local parameters:
C
C U = Pseudo-random number uniformly distributed in the
C     interval (0,1).
C X = Pseudo-random number in the range 0 to 3 whose frac-
C       tional part is U.
C
      IX = MOD(171*IX,30269)
      IY = MOD(172*IY,30307)
      IZ = MOD(170*IZ,30323)
      X = (REAL(IX)/30269.) + (REAL(IY)/30307.) +
     .    (REAL(IZ)/30323.)
      U = X - INT(X)
      JRAND = REAL(N)*U + 1.
      RETURN
      END
      LOGICAL FUNCTION LEFT (X1,Y1,X2,Y2,X0,Y0)
      REAL    X1, Y1, X2, Y2, X0, Y0
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This function determines whether node N0 is to the left
C or to the right of the line through N1-N2 as viewed by an
C observer at N1 facing N2.
C
C
C On input:
C
C       X1,Y1 = Coordinates of N1.
C
C       X2,Y2 = Coordinates of N2.
C
C       X0,Y0 = Coordinates of N0.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       LEFT = .TRUE. if and only if (X0,Y0) is on or to the
C              left of the directed line N1->N2.
C
C Modules required by LEFT:  None
C
C***********************************************************
C
      REAL DX1, DY1, DX2, DY2
C
C Local parameters:
C
C DX1,DY1 = X,Y components of the vector N1->N2
C DX2,DY2 = X,Y components of the vector N1->N0
C
      DX1 = X2-X1
      DY1 = Y2-Y1
      DX2 = X0-X1
      DY2 = Y0-Y1
C
C If the sign of the vector cross product of N1->N2 and
C   N1->N0 is positive, then sin(A) > 0, where A is the
C   angle between the vectors, and thus A is in the range
C   (0,180) degrees.
C
      LEFT = DX1*DY2 .GE. DX2*DY1
      RETURN
      END
      INTEGER FUNCTION LSTPTR (LPL,NB,LIST,LPTR)
      INTEGER LPL, NB, LIST(*), LPTR(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This function returns the index (LIST pointer) of NB in
C the adjacency list for N0, where LPL = LEND(N0).
C
C
C On input:
C
C       LPL = LEND(N0)
C
C       NB = Index of the node whose pointer is to be re-
C            turned.  NB must be connected to N0.
C
C       LIST,LPTR = Data structure defining the triangula-
C                   tion.  Refer to Subroutine TRMESH.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       LSTPTR = Pointer such that LIST(LSTPTR) = NB or
C                LIST(LSTPTR) = -NB, unless NB is not a
C                neighbor of N0, in which case LSTPTR = LPL.
C
C Modules required by LSTPTR:  None
C
C***********************************************************
C
      INTEGER LP, ND
C
      LP = LPTR(LPL)
    1 ND = LIST(LP)
        IF (ND .EQ. NB) GO TO 2
        LP = LPTR(LP)
        IF (LP .NE. LPL) GO TO 1
C
    2 LSTPTR = LP
      RETURN
      END
      INTEGER FUNCTION NBCNT (LPL,LPTR)
      INTEGER LPL, LPTR(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This function returns the number of neighbors of a node
C N0 in a triangulation created by Subroutine TRMESH (or
C TRMSHR).
C
C
C On input:
C
C       LPL = LIST pointer to the last neighbor of N0 --
C             LPL = LEND(N0).
C
C       LPTR = Array of pointers associated with LIST.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       NBCNT = Number of neighbors of N0.
C
C Modules required by NBCNT:  None
C
C***********************************************************
C
      INTEGER K, LP
C
      LP = LPL
      K = 1
C
    1 LP = LPTR(LP)
        IF (LP .EQ. LPL) GO TO 2
        K = K + 1
        GO TO 1
C
    2 NBCNT = K
      RETURN
      END
      INTEGER FUNCTION NEARND (XP,YP,IST,N,X,Y,LIST,LPTR,
     .                         LEND, DSQ)
      INTEGER IST, N, LIST(*), LPTR(*), LEND(N)
      REAL    XP, YP, X(N), Y(N), DSQ
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/27/98
C
C   Given a point P in the plane and a Delaunay triangula-
C tion created by Subroutine TRMESH or TRMSHR, this function
C returns the index of the nearest triangulation node to P.
C
C   The algorithm consists of implicitly adding P to the
C triangulation, finding the nearest neighbor to P, and
C implicitly deleting P from the triangulation.  Thus, it
C is based on the fact that, if P is a node in a Delaunay
C triangulation, the nearest node to P is a neighbor of P.
C
C
C On input:
C
C       XP,YP = Cartesian coordinates of the point P to be
C               located relative to the triangulation.
C
C       IST = Index of a node at which TRFIND begins the
C             search.  Search time depends on the proximity
C             of this node to P.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRMESH.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       NEARND = Nodal index of the nearest node to P, or 0
C                if N < 3 or the triangulation data struc-
C                ture is invalid.
C
C       DSQ = Squared distance between P and NEARND unless
C             NEARND = 0.
C
C       Note that the number of candidates for NEARND
C       (neighbors of P) is limited to LMAX defined in
C       the PARAMETER statement below.
C
C Modules required by NEARND:  JRAND, LEFT, LSTPTR, TRFIND
C
C Intrinsic function called by NEARND:  ABS
C
C***********************************************************
C
      INTEGER   LSTPTR
      INTEGER   LMAX
      PARAMETER (LMAX=25)
      INTEGER   I1, I2, I3, L, LISTP(LMAX), LP, LP1, LP2,
     .          LPL, LPTRP(LMAX), N1, N2, N3, NR, NST
      REAL      COS1, COS2, DS1, DSR, DX11, DX12, DX21,
     .          DX22, DY11, DY12, DY21, DY22, SIN1, SIN2
C
C Store local parameters and test for N invalid.
C
      IF (N .LT. 3) GO TO 7
      NST = IST
      IF (NST .LT. 1  .OR.  NST .GT. N) NST = 1
C
C Find a triangle (I1,I2,I3) containing P, or the rightmost
C   (I1) and leftmost (I2) visible boundary nodes as viewed
C   from P.
C
      CALL TRFIND (NST,XP,YP,N,X,Y,LIST,LPTR,LEND, I1,I2,I3)
C
C Test for collinear nodes.
C
      IF (I1 .EQ. 0) GO TO 7
C
C Store the linked list of 'neighbors' of P in LISTP and
C   LPTRP.  I1 is the first neighbor, and 0 is stored as
C   the last neighbor if P is not contained in a triangle.
C   L is the length of LISTP and LPTRP, and is limited to
C   LMAX.
C
      IF (I3 .NE. 0) THEN
        LISTP(1) = I1
        LPTRP(1) = 2
        LISTP(2) = I2
        LPTRP(2) = 3
        LISTP(3) = I3
        LPTRP(3) = 1
        L = 3
      ELSE
        N1 = I1
        L = 1
        LP1 = 2
        LISTP(L) = N1
        LPTRP(L) = LP1
C
C   Loop on the ordered sequence of visible boundary nodes
C     N1 from I1 to I2.
C
    1   LPL = LEND(N1)
          N1 = -LIST(LPL)
          L = LP1
          LP1 = L+1
          LISTP(L) = N1
          LPTRP(L) = LP1
          IF (N1 .NE. I2  .AND.  LP1 .LT. LMAX) GO TO 1
        L = LP1
        LISTP(L) = 0
        LPTRP(L) = 1
      ENDIF
C
C Initialize variables for a loop on arcs N1-N2 opposite P
C   in which new 'neighbors' are 'swapped' in.  N1 follows
C   N2 as a neighbor of P, and LP1 and LP2 are the LISTP
C   indexes of N1 and N2.
C
      LP2 = 1
      N2 = I1
      LP1 = LPTRP(1)
      N1 = LISTP(LP1)
C
C Begin loop:  find the node N3 opposite N1->N2.
C
    2 LP = LSTPTR(LEND(N1),N2,LIST,LPTR)
        IF (LIST(LP) .LT. 0) GO TO 4
        LP = LPTR(LP)
        N3 = ABS(LIST(LP))
C
C Swap test:  Exit the loop if L = LMAX.
C
        IF (L .EQ. LMAX) GO TO 5
        DX11 = X(N1) - X(N3)
        DX12 = X(N2) - X(N3)
        DX22 = X(N2) - XP
        DX21 = X(N1) - XP
C
        DY11 = Y(N1) - Y(N3)
        DY12 = Y(N2) - Y(N3)
        DY22 = Y(N2) - YP
        DY21 = Y(N1) - YP
C
        COS1 = DX11*DX12 + DY11*DY12
        COS2 = DX22*DX21 + DY22*DY21
        IF (COS1 .GE. 0.  .AND.  COS2 .GE. 0.) GO TO 4
        IF (COS1 .LT. 0.  .AND.  COS2 .LT. 0.) GO TO 3
C
        SIN1 = DX11*DY12 - DX12*DY11
        SIN2 = DX22*DY21 - DX21*DY22
        IF (SIN1*COS2 + COS1*SIN2 .GE. 0.) GO TO 4
C
C Swap:  Insert N3 following N2 in the adjacency list for P.
C        The two new arcs opposite P must be tested.
C
    3   L = L+1
        LPTRP(LP2) = L
        LISTP(L) = N3
        LPTRP(L) = LP1
        LP1 = L
        N1 = N3
        GO TO 2
C
C No swap:  Advance to the next arc and test for termination
C           on N1 = I1 (LP1 = 1) or N1 followed by 0.
C
    4   IF (LP1 .EQ. 1) GO TO 5
        LP2 = LP1
        N2 = N1
        LP1 = LPTRP(LP1)
        N1 = LISTP(LP1)
        IF (N1 .EQ. 0) GO TO 5
        GO TO 2
C
C Set NR and DSR to the index of the nearest node to P and
C   its squared distance from P, respectively.
C
    5 NR = I1
      DSR = (X(NR)-XP)**2 + (Y(NR)-YP)**2
      DO 6 LP = 2,L
        N1 = LISTP(LP)
        IF (N1 .EQ. 0) GO TO 6
        DS1 = (X(N1)-XP)**2 + (Y(N1)-YP)**2
        IF (DS1 .LT. DSR) THEN
          NR = N1
          DSR = DS1
        ENDIF
    6   CONTINUE
      DSQ = DSR
      NEARND = NR
      RETURN
C
C Invalid input.
C
    7 NEARND = 0
      RETURN
      END
      SUBROUTINE OPTIM (X,Y,NA, LIST,LPTR,LEND,NIT,IWK, IER)
      INTEGER NA, LIST(*), LPTR(*), LEND(*), NIT, IWK(2,NA),
     .        IER
      REAL    X(*), Y(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/27/98
C
C   Given a set of NA triangulation arcs, this subroutine
C optimizes the portion of the triangulation consisting of
C the quadrilaterals (pairs of adjacent triangles) which
C have the arcs as diagonals by applying the circumcircle
C test and appropriate swaps to the arcs.
C
C   An iteration consists of applying the swap test and
C swaps to all NA arcs in the order in which they are
C stored.  The iteration is repeated until no swap occurs
C or NIT iterations have been performed.  The bound on the
C number of iterations may be necessary to prevent an
C infinite loop caused by cycling (reversing the effect of a
C previous swap) due to floating point inaccuracy when four
C or more nodes are nearly cocircular.
C
C
C On input:
C
C       X,Y = Arrays containing the nodal coordinates.
C
C       NA = Number of arcs in the set.  NA .GE. 0.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C       NIT = Maximum number of iterations to be performed.
C             A reasonable value is 3*NA.  NIT .GE. 1.
C
C       IWK = Integer array dimensioned 2 by NA containing
C             the nodal indexes of the arc endpoints (pairs
C             of endpoints are stored in columns).
C
C On output:
C
C       LIST,LPTR,LEND = Updated triangulation data struc-
C                        ture reflecting the swaps.
C
C       NIT = Number of iterations performed.
C
C       IWK = Endpoint indexes of the new set of arcs
C             reflecting the swaps.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if a swap occurred on the last of
C                     MAXIT iterations, where MAXIT is the
C                     value of NIT on input.  The new set
C                     of arcs in not necessarily optimal
C                     in this case.
C             IER = 2 if NA < 0 or NIT < 1 on input.
C             IER = 3 if IWK(2,I) is not a neighbor of
C                     IWK(1,I) for some I in the range 1
C                     to NA.  A swap may have occurred in
C                     this case.
C             IER = 4 if a zero pointer was returned by
C                     Subroutine SWAP.
C
C Modules required by OPTIM:  LSTPTR, SWAP, SWPTST
C
C Intrinsic function called by OPTIM:  ABS
C
C***********************************************************
C
      LOGICAL SWPTST
      INTEGER I, IO1, IO2, ITER, LP, LP21, LPL, LPP, MAXIT,
     .        N1, N2, NNA
      LOGICAL SWP
C
C Local parameters:
C
C I =       Column index for IWK
C IO1,IO2 = Nodal indexes of the endpoints of an arc in IWK
C ITER =    Iteration count
C LP =      LIST pointer
C LP21 =    Parameter returned by SWAP (not used)
C LPL =     Pointer to the last neighbor of IO1
C LPP =     Pointer to the node preceding IO2 as a neighbor
C             of IO1
C MAXIT =   Input value of NIT
C N1,N2 =   Nodes opposite IO1->IO2 and IO2->IO1,
C             respectively
C NNA =     Local copy of NA
C SWP =     Flag set to TRUE iff a swap occurs in the
C             optimization loop
C
      NNA = NA
      MAXIT = NIT
      IF (NNA .LT. 0  .OR.  MAXIT .LT. 1) GO TO 7
C
C Initialize iteration count ITER and test for NA = 0.
C
      ITER = 0
      IF (NNA .EQ. 0) GO TO 5
C
C Top of loop --
C   SWP = TRUE iff a swap occurred in the current iteration.
C
    1 IF (ITER .EQ. MAXIT) GO TO 6
      ITER = ITER + 1
      SWP = .FALSE.
C
C   Inner loop on arcs IO1-IO2 --
C
      DO 4 I = 1,NNA
        IO1 = IWK(1,I)
        IO2 = IWK(2,I)
C
C   Set N1 and N2 to the nodes opposite IO1->IO2 and
C     IO2->IO1, respectively.  Determine the following:
C
C     LPL = pointer to the last neighbor of IO1,
C     LP = pointer to IO2 as a neighbor of IO1, and
C     LPP = pointer to the node N2 preceding IO2.
C
        LPL = LEND(IO1)
        LPP = LPL
        LP = LPTR(LPP)
    2   IF (LIST(LP) .EQ. IO2) GO TO 3
          LPP = LP
          LP = LPTR(LPP)
          IF (LP .NE. LPL) GO TO 2
C
C   IO2 should be the last neighbor of IO1.  Test for no
C     arc and bypass the swap test if IO1 is a boundary
C     node.
C
        IF (ABS(LIST(LP)) .NE. IO2) GO TO 8
        IF (LIST(LP) .LT. 0) GO TO 4
C
C   Store N1 and N2, or bypass the swap test if IO1 is a
C     boundary node and IO2 is its first neighbor.
C
    3   N2 = LIST(LPP)
        IF (N2 .LT. 0) GO TO 4
        LP = LPTR(LP)
        N1 = ABS(LIST(LP))
C
C   Test IO1-IO2 for a swap, and update IWK if necessary.
C
        IF ( .NOT. SWPTST(N1,N2,IO1,IO2,X,Y) ) GO TO 4
        CALL SWAP (N1,N2,IO1,IO2, LIST,LPTR,LEND, LP21)
        IF (LP21 .EQ. 0) GO TO 9
        SWP = .TRUE.
        IWK(1,I) = N1
        IWK(2,I) = N2
    4   CONTINUE
      IF (SWP) GO TO 1
C
C Successful termination.
C
    5 NIT = ITER
      IER = 0
      RETURN
C
C MAXIT iterations performed without convergence.
C
    6 NIT = MAXIT
      IER = 1
      RETURN
C
C Invalid input parameter.
C
    7 NIT = 0
      IER = 2
      RETURN
C
C IO2 is not a neighbor of IO1.
C
    8 NIT = ITER
      IER = 3
      RETURN
C
C Zero pointer returned by SWAP.
C
    9 NIT = ITER
      IER = 4
      RETURN
      END
      REAL FUNCTION STORE (X)
      REAL X
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   03/18/90
C
C   This function forces its argument X to be stored in a
C memory location, thus providing a means of determining
C floating point number characteristics (such as the machine
C precision) when it is necessary to avoid computation in
C high precision registers.
C
C
C On input:
C
C       X = Value to be stored.
C
C X is not altered by this function.
C
C On output:
C
C       STORE = Value of X after it has been stored and
C               possibly truncated or rounded to the single
C               precision word length.
C
C Modules required by STORE:  None
C
C***********************************************************
C
      REAL Y
      COMMON/STCOM/Y
C
      Y = X
      STORE = Y
      RETURN
      END
      SUBROUTINE SWAP (IN1,IN2,IO1,IO2, LIST,LPTR,
     .                 LEND, LP21)
      INTEGER IN1, IN2, IO1, IO2, LIST(*), LPTR(*), LEND(*),
     .        LP21
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/22/98
C
C   Given a triangulation of a set of points on the unit
C sphere, this subroutine replaces a diagonal arc in a
C strictly convex quadrilateral (defined by a pair of adja-
C cent triangles) with the other diagonal.  Equivalently, a
C pair of adjacent triangles is replaced by another pair
C having the same union.
C
C
C On input:
C
C       IN1,IN2,IO1,IO2 = Nodal indexes of the vertices of
C                         the quadrilateral.  IO1-IO2 is re-
C                         placed by IN1-IN2.  (IO1,IO2,IN1)
C                         and (IO2,IO1,IN2) must be trian-
C                         gles on input.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C On output:
C
C       LIST,LPTR,LEND = Data structure updated with the
C                        swap -- triangles (IO1,IO2,IN1) and
C                        (IO2,IO1,IN2) are replaced by
C                        (IN1,IN2,IO2) and (IN2,IN1,IO1)
C                        unless LP21 = 0.
C
C       LP21 = Index of IN1 as a neighbor of IN2 after the
C              swap is performed unless IN1 and IN2 are
C              adjacent on input, in which case LP21 = 0.
C
C Module required by SWAP:  LSTPTR
C
C Intrinsic function called by SWAP:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER LP, LPH, LPSAV
C
C Local parameters:
C
C LP,LPH,LPSAV = LIST pointers
C
C
C Test for IN1 and IN2 adjacent.
C
      LP = LSTPTR(LEND(IN1),IN2,LIST,LPTR)
      IF (ABS(LIST(LP)) .EQ. IN2) THEN
        LP21 = 0
        RETURN
      ENDIF
C
C Delete IO2 as a neighbor of IO1.
C
      LP = LSTPTR(LEND(IO1),IN2,LIST,LPTR)
      LPH = LPTR(LP)
      LPTR(LP) = LPTR(LPH)
C
C If IO2 is the last neighbor of IO1, make IN2 the
C   last neighbor.
C
      IF (LEND(IO1) .EQ. LPH) LEND(IO1) = LP
C
C Insert IN2 as a neighbor of IN1 following IO1
C   using the hole created above.
C
      LP = LSTPTR(LEND(IN1),IO1,LIST,LPTR)
      LPSAV = LPTR(LP)
      LPTR(LP) = LPH
      LIST(LPH) = IN2
      LPTR(LPH) = LPSAV
C
C Delete IO1 as a neighbor of IO2.
C
      LP = LSTPTR(LEND(IO2),IN1,LIST,LPTR)
      LPH = LPTR(LP)
      LPTR(LP) = LPTR(LPH)
C
C If IO1 is the last neighbor of IO2, make IN1 the
C   last neighbor.
C
      IF (LEND(IO2) .EQ. LPH) LEND(IO2) = LP
C
C Insert IN1 as a neighbor of IN2 following IO2.
C
      LP = LSTPTR(LEND(IN2),IO2,LIST,LPTR)
      LPSAV = LPTR(LP)
      LPTR(LP) = LPH
      LIST(LPH) = IN1
      LPTR(LPH) = LPSAV
      LP21 = LPH
      RETURN
      END
      LOGICAL FUNCTION SWPTST (IN1,IN2,IO1,IO2,X,Y)
      USE SWPCOM, ONLY: SWTOL
      INTEGER IN1, IN2, IO1, IO2
      REAL    X(*), Y(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This function applies the circumcircle test to a quadri-
C lateral defined by a pair of adjacent triangles.  The
C diagonal arc (shared triangle side) should be swapped for
C the other diagonl if and only if the fourth vertex is
C strictly interior to the circumcircle of one of the
C triangles (the decision is independent of the choice of
C triangle).  Equivalently, the diagonal is chosen to maxi-
C mize the smallest of the six interior angles over the two
C pairs of possible triangles (the decision is for no swap
C if the quadrilateral is not strictly convex).
C
C   When the four vertices are nearly cocircular (the
C neutral case), the preferred decision is no swap -- in
C order to avoid unnecessary swaps and, more important, to
C avoid cycling in Subroutine OPTIM which is called by
C DELNOD and EDGE.  Thus, a tolerance SWTOL (stored in
C SWPCOM by TRMESH or TRMSHR) is used to define 'nearness'
C to the neutral case.
C
C
C On input:
C
C       IN1,IN2,IO1,IO2 = Nodal indexes of the vertices of
C                         the quadrilateral.  IO1-IO2 is the
C                         triangulation arc (shared triangle
C                         side) to be replaced by IN1-IN2 if
C                         the decision is to swap.  The
C                         triples (IO1,IO2,IN1) and (IO2,
C                         IO1,IN2) must define triangles (be
C                         in counterclockwise order) on in-
C                         put.
C
C       X,Y = Arrays containing the nodal coordinates.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       SWPTST = .TRUE. if and only if the arc connecting
C                IO1 and IO2 is to be replaced.
C
C Modules required by SWPTST:  None
C
C***********************************************************
C
      REAL DX11, DX12, DX22, DX21, DY11, DY12, DY22, DY21,
     .     SIN1, SIN2, COS1, COS2, SIN12
C
C Tolerance stored by TRMESH or TRMSHR.
C
C
C Local parameters:
C
C DX11,DY11 = X,Y components of the vector IN1->IO1
C DX12,DY12 = X,Y components of the vector IN1->IO2
C DX22,DY22 = X,Y components of the vector IN2->IO2
C DX21,DY21 = X,Y components of the vector IN2->IO1
C SIN1 =      Cross product of the vectors IN1->IO1 and
C               IN1->IO2 -- proportional to sin(T1), where
C               T1 is the angle at IN1 formed by the vectors
C COS1 =      Inner product of the vectors IN1->IO1 and
C               IN1->IO2 -- proportional to cos(T1)
C SIN2 =      Cross product of the vectors IN2->IO2 and
C               IN2->IO1 -- proportional to sin(T2), where
C               T2 is the angle at IN2 formed by the vectors
C COS2 =      Inner product of the vectors IN2->IO2 and
C               IN2->IO1 -- proportional to cos(T2)
C SIN12 =     SIN1*COS2 + COS1*SIN2 -- proportional to
C               sin(T1+T2)
C
C
C Compute the vectors containing the angles T1 and T2.
C
      DX11 = X(IO1) - X(IN1)
      DX12 = X(IO2) - X(IN1)
      DX22 = X(IO2) - X(IN2)
      DX21 = X(IO1) - X(IN2)
C
      DY11 = Y(IO1) - Y(IN1)
      DY12 = Y(IO2) - Y(IN1)
      DY22 = Y(IO2) - Y(IN2)
      DY21 = Y(IO1) - Y(IN2)
C
C Compute inner products.
C
      COS1 = DX11*DX12 + DY11*DY12
      COS2 = DX22*DX21 + DY22*DY21
C
C The diagonals should be swapped iff (T1+T2) > 180
C   degrees.  The following two tests ensure numerical
C   stability:  the decision must be FALSE when both
C   angles are close to 0, and TRUE when both angles
C   are close to 180 degrees.
C
      IF (COS1 .GE. 0.  .AND.  COS2 .GE. 0.) GO TO 2
      IF (COS1 .LT. 0.  .AND.  COS2 .LT. 0.) GO TO 1
C
C Compute vector cross products (Z-components).
C
      SIN1 = DX11*DY12 - DX12*DY11
      SIN2 = DX22*DY21 - DX21*DY22
      SIN12 = SIN1*COS2 + COS1*SIN2
      IF (SIN12 .GE. -SWTOL) GO TO 2
C
C Swap.
C
    1 SWPTST = .TRUE.
      RETURN
C
C No swap.
C
    2 SWPTST = .FALSE.
      RETURN
      END
      SUBROUTINE TRFIND (NST,PX,PY,N,X,Y,LIST,LPTR,LEND, I1,
     .                   I2,I3)
      INTEGER NST, N, LIST(*), LPTR(*), LEND(N), I1, I2, I3
      REAL    PX, PY, X(N), Y(N)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/28/98
C
C   This subroutine locates a point P relative to a triangu-
C lation created by Subroutine TRMESH or TRMSHR.  If P is
C contained in a triangle, the three vertex indexes are
C returned.  Otherwise, the indexes of the rightmost and
C leftmost visible boundary nodes are returned.
C
C
C On input:
C
C       NST = Index of a node at which TRFIND begins the
C             search.  Search time depends on the proximity
C             of this node to P.
C
C       PX,PY = X and y coordinates of the point P to be
C               located.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes in the triangulation.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       I1,I2,I3 = Nodal indexes, in counterclockwise order,
C                  of the vertices of a triangle containing
C                  P if P is contained in a triangle.  If P
C                  is not in the convex hull of the nodes,
C                  I1 indexes the rightmost visible boundary
C                  node, I2 indexes the leftmost visible
C                  boundary node, and I3 = 0.  Rightmost and
C                  leftmost are defined from the perspective
C                  of P, and a pair of points are visible
C                  from each other if and only if the line
C                  segment joining them intersects no trian-
C                  gulation arc.  If P and all of the nodes
C                  lie on a common line, then I1 = I2 = I3 =
C                  0 on output.
C
C Modules required by TRFIND:  JRAND, LEFT, LSTPTR, STORE
C
C Intrinsic function called by TRFIND:  ABS
C
C***********************************************************
C
      INTEGER JRAND, LSTPTR
      LOGICAL LEFT
      REAL    STORE
      INTEGER IX, IY, IZ, LP, N0, N1, N1S, N2, N2S, N3, N4,
     .        NB, NF, NL, NP, NPP
      LOGICAL FRWRD
      REAL    B1, B2, XA, XB, XC, XP, YA, YB, YC, YP
C
      SAVE    IX, IY, IZ
      DATA    IX/1/, IY/2/, IZ/3/
C
C Local parameters:
C
C B1,B2 =    Unnormalized barycentric coordinates of P with
C              respect to (N1,N2,N3)
C IX,IY,IZ = Integer seeds for JRAND
C LP =       LIST pointer
C N0,N1,N2 = Nodes in counterclockwise order defining a
C              cone (with vertex N0) containing P
C N1S,N2S =  Saved values of N1 and N2
C N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
C NB =       Index of a boundary node -- first neighbor of
C              NF or last neighbor of NL in the boundary
C              traversal loops
C NF,NL =    First and last neighbors of N0, or first
C              (rightmost) and last (leftmost) nodes
C              visible from P when P is exterior to the
C              triangulation
C NP,NPP =   Indexes of boundary nodes used in the boundary
C              traversal loops
C XA,XB,XC = Dummy arguments for FRWRD
C YA,YB,YC = Dummy arguments for FRWRD
C XP,YP =    Local variables containing the components of P
C
C Statement function:
C
C FRWRD = TRUE iff C is forward of A->B
C              iff <A->B,A->C> .GE. 0.
C
      FRWRD(XA,YA,XB,YB,XC,YC) = (XB-XA)*(XC-XA) +
     .                           (YB-YA)*(YC-YA) .GE. 0.
C
C Initialize variables.
C
      XP = PX
      YP = PY
      N0 = NST
      IF (N0 .LT. 1  .OR.  N0 .GT. N)
     .  N0 = JRAND(N, IX,IY,IZ )
C
C Set NF and NL to the first and last neighbors of N0, and
C   initialize N1 = NF.
C
    1 LP = LEND(N0)
      NL = LIST(LP)
      LP = LPTR(LP)
      NF = LIST(LP)
      N1 = NF
C
C Find a pair of adjacent neighbors N1,N2 of N0 that define
C   a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.
C
      IF (NL .GT. 0) GO TO 2
C
C   N0 is a boundary node.  Test for P exterior.
C
      NL = -NL
      IF ( .NOT. LEFT(X(N0),Y(N0),X(NF),Y(NF),XP,YP) ) THEN
        NL = N0
        GO TO 9
      ENDIF
      IF ( .NOT. LEFT(X(NL),Y(NL),X(N0),Y(N0),XP,YP) ) THEN
        NB = NF
        NF = N0
        NP = NL
        NPP = N0
        GO TO 11
      ENDIF
      GO TO 3
C
C   N0 is an interior node.  Find N1.
C
    2 IF ( LEFT(X(N0),Y(N0),X(N1),Y(N1),XP,YP) ) GO TO 3
        LP = LPTR(LP)
        N1 = LIST(LP)
        IF (N1 .EQ. NL) GO TO 6
        GO TO 2
C
C   P is to the left of edge N0->N1.  Initialize N2 to the
C     next neighbor of N0.
C
    3 LP = LPTR(LP)
        N2 = ABS(LIST(LP))
        IF ( .NOT. LEFT(X(N0),Y(N0),X(N2),Y(N2),XP,YP) )
     .    GO TO 7
        N1 = N2
        IF (N1 .NE. NL) GO TO 3
      IF ( .NOT. LEFT(X(N0),Y(N0),X(NF),Y(NF),XP,YP) )
     .  GO TO 6
      IF (XP .EQ. X(N0) .AND. YP .EQ. Y(N0)) GO TO 5
C
C   P is left of or on edges N0->NB for all neighbors NB
C     of N0.
C   All points are collinear iff P is left of NB->N0 for
C     all neighbors NB of N0.  Search the neighbors of N0.
C     NOTE -- N1 = NL and LP points to NL.
C
    4 IF ( .NOT. LEFT(X(N1),Y(N1),X(N0),Y(N0),XP,YP) )
     .  GO TO 5
        LP = LPTR(LP)
        N1 = ABS(LIST(LP))
        IF (N1 .EQ. NL) GO TO 17
        GO TO 4
C
C   P is to the right of N1->N0, or P=N0.  Set N0 to N1 and
C     start over.
C
    5 N0 = N1
      GO TO 1
C
C   P is between edges N0->N1 and N0->NF.
C
    6 N2 = NF
C
C P is contained in the wedge defined by line segments
C   N0->N1 and N0->N2, where N1 is adjacent to N2.  Set
C   N3 to the node opposite N1->N2, and save N1 and N2 to
C   test for cycling.
C
    7 N3 = N0
      N1S = N1
      N2S = N2
C
C Top of edge hopping loop.  Test for termination.
C
    8 IF ( LEFT(X(N1),Y(N1),X(N2),Y(N2),XP,YP) ) THEN
C
C   P LEFT N1->N2 and hence P is in (N1,N2,N3) unless an
C     error resulted from floating point inaccuracy and
C     collinearity.  Compute the unnormalized barycentric
C     coordinates of P with respect to (N1,N2,N3).
C
        B1 = (X(N3)-X(N2))*(YP-Y(N2)) -
     .       (XP-X(N2))*(Y(N3)-Y(N2))
        B2 = (X(N1)-X(N3))*(YP-Y(N3)) -
     .       (XP-X(N3))*(Y(N1)-Y(N3))
        IF (STORE(B1+1.) .GE. 1.  .AND.
     .      STORE(B2+1.) .GE. 1.) GO TO 16
C
C   Restart with N0 randomly selected.
C
        N0 = JRAND(N, IX,IY,IZ )
        GO TO 1
      ENDIF
C
C   Set N4 to the neighbor of N2 which follows N1 (node
C     opposite N2->N1) unless N1->N2 is a boundary edge.
C
      LP = LSTPTR(LEND(N2),N1,LIST,LPTR)
      IF (LIST(LP) .LT. 0) THEN
        NF = N2
        NL = N1
        GO TO 9
      ENDIF
      LP = LPTR(LP)
      N4 = ABS(LIST(LP))
C
C   Select the new edge N1->N2 which intersects the line
C     segment N0-P, and set N3 to the node opposite N1->N2.
C
      IF ( LEFT(X(N0),Y(N0),X(N4),Y(N4),XP,YP) ) THEN
        N3 = N1
        N1 = N4
        N2S = N2
        IF (N1 .NE. N1S  .AND.  N1 .NE. N0) GO TO 8
      ELSE
        N3 = N2
        N2 = N4
        N1S = N1
        IF (N2 .NE. N2S  .AND.  N2 .NE. N0) GO TO 8
      ENDIF
C
C   The starting node N0 or edge N1-N2 was encountered
C     again, implying a cycle (infinite loop).  Restart
C     with N0 randomly selected.
C
      N0 = JRAND(N, IX,IY,IZ )
      GO TO 1
C
C Boundary traversal loops.  NL->NF is a boundary edge and
C   P RIGHT NL->NF.  Save NL and NF.

    9 NP = NL
      NPP = NF
C
C Find the first (rightmost) visible boundary node NF.  NB
C   is set to the first neighbor of NF, and NP is the last
C   neighbor.
C
   10 LP = LEND(NF)
      LP = LPTR(LP)
      NB = LIST(LP)
      IF ( .NOT. LEFT(X(NF),Y(NF),X(NB),Y(NB),XP,YP) )
     .  GO TO 12
C
C   P LEFT NF->NB and thus NB is not visible unless an error
C     resulted from floating point inaccuracy and collinear-
C     ity of the 4 points NP, NF, NB, and P.
C
   11 IF ( FRWRD(X(NF),Y(NF),X(NP),Y(NP),XP,YP)  .OR.
     .     FRWRD(X(NF),Y(NF),X(NP),Y(NP),X(NB),Y(NB)) ) THEN
        I1 = NF
        GO TO 13
      ENDIF
C
C   Bottom of loop.
C
   12 NP = NF
      NF = NB
      GO TO 10
C
C Find the last (leftmost) visible boundary node NL.  NB
C   is set to the last neighbor of NL, and NPP is the first
C   neighbor.
C
   13 LP = LEND(NL)
      NB = -LIST(LP)
      IF ( .NOT. LEFT(X(NB),Y(NB),X(NL),Y(NL),XP,YP) )
     .  GO TO 14
C
C   P LEFT NB->NL and thus NB is not visible unless an error
C     resulted from floating point inaccuracy and collinear-
C     ity of the 4 points P, NB, NL, and NPP.
C
      IF ( FRWRD(X(NL),Y(NL),X(NPP),Y(NPP),XP,YP)  .OR.
     .     FRWRD(X(NL),Y(NL),X(NPP),Y(NPP),X(NB),Y(NB)) )
     .  GO TO 15
C
C   Bottom of loop.
C
   14 NPP = NL
      NL = NB
      GO TO 13
C
C NL is the leftmost visible boundary node.
C
   15 I2 = NL
      I3 = 0
      RETURN
C
C P is in the triangle (N1,N2,N3).
C
   16 I1 = N1
      I2 = N2
      I3 = N3
      RETURN
C
C All points are collinear.
C
   17 I1 = 0
      I2 = 0
      I3 = 0
      RETURN
      END
      SUBROUTINE TRLIST (NCC,LCC,N,LIST,LPTR,LEND,NROW, NT,
     .                   LTRI,LCT,IER)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        NROW, NT, LTRI(NROW,*), LCT(*), IER
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   03/22/97
C
C   This subroutine converts a triangulation data structure
C from the linked list created by Subroutine TRMESH or
C TRMSHR to a triangle list.
C
C On input:
C
C       NCC = Number of constraints.  NCC .GE. 0.
C
C       LCC = List of constraint curve starting indexes (or
C             dummy array of length 1 if NCC = 0).  Refer to
C             Subroutine ADDCST.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       LIST,LPTR,LEND = Linked list data structure defin-
C                        ing the triangulation.  Refer to
C                        Subroutine TRMESH.
C
C       NROW = Number of rows (entries per triangle) re-
C              served for the triangle list LTRI.  The value
C              must be 6 if only the vertex indexes and
C              neighboring triangle indexes are to be
C              stored, or 9 if arc indexes are also to be
C              assigned and stored.  Refer to LTRI.
C
C The above parameters are not altered by this routine.
C
C       LTRI = Integer array of length at least NROW*NT,
C              where NT is at most 2N-5.  (A sufficient
C              length is 12N if NROW=6 or 18N if NROW=9.)
C
C       LCT = Integer array of length NCC or dummy array of
C             length 1 if NCC = 0.
C
C On output:
C
C       NT = Number of triangles in the triangulation unless
C            IER .NE. 0, in which case NT = 0.  NT = 2N - NB
C            - 2, where NB is the number of boundary nodes.
C
C       LTRI = NROW by NT array whose J-th column contains
C              the vertex nodal indexes (first three rows),
C              neighboring triangle indexes (second three
C              rows), and, if NROW = 9, arc indexes (last
C              three rows) associated with triangle J for
C              J = 1,...,NT.  The vertices are ordered
C              counterclockwise with the first vertex taken
C              to be the one with smallest index.  Thus,
C              LTRI(2,J) and LTRI(3,J) are larger than
C              LTRI(1,J) and index adjacent neighbors of
C              node LTRI(1,J).  For I = 1,2,3, LTRI(I+3,J)
C              and LTRI(I+6,J) index the triangle and arc,
C              respectively, which are opposite (not shared
C              by) node LTRI(I,J), with LTRI(I+3,J) = 0 if
C              LTRI(I+6,J) indexes a boundary arc.  Vertex
C              indexes range from 1 to N, triangle indexes
C              from 0 to NT, and, if included, arc indexes
C              from 1 to NA = NT+N-1.  The triangles are or-
C              dered on first (smallest) vertex indexes,
C              except that the sets of constraint triangles
C              (triangles contained in the closure of a con-
C              straint region) follow the non-constraint
C              triangles.
C
C       LCT = Array of length NCC containing the triangle
C             index of the first triangle of constraint J in
C             LCT(J).  Thus, the number of non-constraint
C             triangles is LCT(1)-1, and constraint J con-
C             tains LCT(J+1)-LCT(J) triangles, where
C             LCT(NCC+1) = NT+1.
C
C       IER = Error indicator.
C             IER = 0 if no errors were encountered.
C             IER = 1 if NCC, N, NROW, or an LCC entry is
C                     outside its valid range on input.
C             IER = 2 if the triangulation data structure
C                     (LIST,LPTR,LEND) is invalid.  Note,
C                     however, that these arrays are not
C                     completely tested for validity.
C
C Modules required by TRLIST:  None
C
C Intrinsic function called by TRLIST:  ABS
C
C***********************************************************
C
      INTEGER I, I1, I2, I3, ISV, J, JLAST, KA, KN, KT, L,
     .        LCC1, LP, LP2, LPL, LPLN1, N1, N1ST, N2, N3,
     .        NM2, NN
      LOGICAL ARCS, CSTRI, PASS2
C
C Test for invalid input parameters and store the index
C   LCC1 of the first constraint node (if any).
C
      NN = N
      IF (NCC .LT. 0  .OR.  (NROW .NE. 6  .AND.
     .    NROW .NE. 9)) GO TO 12
      LCC1 = NN+1
      IF (NCC .EQ. 0) THEN
        IF (NN .LT. 3) GO TO 12
      ELSE
        DO 1 I = NCC,1,-1
          IF (LCC1-LCC(I) .LT. 3) GO TO 12
          LCC1 = LCC(I)
    1     CONTINUE
        IF (LCC1 .LT. 1) GO TO 12
      ENDIF
C
C Initialize parameters for loop on triangles KT = (N1,N2,
C   N3), where N1 < N2 and N1 < N3.  This requires two
C   passes through the nodes with all non-constraint
C   triangles stored on the first pass, and the constraint
C   triangles stored on the second.
C
C   ARCS = TRUE iff arc indexes are to be stored.
C   KA,KT = Numbers of currently stored arcs and triangles.
C   N1ST = Starting index for the loop on nodes (N1ST = 1 on
C            pass 1, and N1ST = LCC1 on pass 2).
C   NM2 = Upper bound on candidates for N1.
C   PASS2 = TRUE iff constraint triangles are to be stored.
C
      ARCS = NROW .EQ. 9
      KA = 0
      KT = 0
      N1ST = 1
      NM2 = NN-2
      PASS2 = .FALSE.
C
C Loop on nodes N1:  J = constraint containing N1,
C                    JLAST = last node in constraint J.
C
    2 J = 0
      JLAST = LCC1 - 1
      DO 11 N1 = N1ST,NM2
        IF (N1 .GT. JLAST) THEN
C
C N1 is the first node in constraint J+1.  Update J and
C   JLAST, and store the first constraint triangle index
C   if in pass 2.
C
          J = J + 1
          IF (J .LT. NCC) THEN
            JLAST = LCC(J+1) - 1
          ELSE
            JLAST = NN
          ENDIF
          IF (PASS2) LCT(J) = KT + 1
        ENDIF
C
C Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
C   to the last neighbor of N1, and LP2 points to N2.
C
        LPLN1 = LEND(N1)
        LP2 = LPLN1
    3     LP2 = LPTR(LP2)
          N2 = LIST(LP2)
          LP = LPTR(LP2)
          N3 = ABS(LIST(LP))
          IF (N2 .LT. N1  .OR.  N3 .LT. N1) GO TO 10
C
C (N1,N2,N3) is a constraint triangle iff the three nodes
C   are in the same constraint and N2 < N3.  Bypass con-
C   straint triangles on pass 1 and non-constraint triangles
C   on pass 2.
C
          CSTRI = N1 .GE. LCC1  .AND.  N2 .LT. N3  .AND.
     .            N3 .LE. JLAST
          IF ((CSTRI  .AND.  .NOT. PASS2)  .OR.
     .        (.NOT. CSTRI  .AND.  PASS2)) GO TO 10
C
C Add a new triangle KT = (N1,N2,N3).
C
          KT = KT + 1
          LTRI(1,KT) = N1
          LTRI(2,KT) = N2
          LTRI(3,KT) = N3
C
C Loop on triangle sides (I1,I2) with neighboring triangles
C   KN = (I1,I2,I3).
C
          DO 9 I = 1,3
            IF (I .EQ. 1) THEN
              I1 = N3
              I2 = N2
            ELSEIF (I .EQ. 2) THEN
              I1 = N1
              I2 = N3
            ELSE
              I1 = N2
              I2 = N1
            ENDIF
C
C Set I3 to the neighbor of I1 which follows I2 unless
C   I2->I1 is a boundary arc.
C
            LPL = LEND(I1)
            LP = LPTR(LPL)
    4       IF (LIST(LP) .EQ. I2) GO TO 5
              LP = LPTR(LP)
              IF (LP .NE. LPL) GO TO 4
C
C   I2 is the last neighbor of I1 unless the data structure
C     is invalid.  Bypass the search for a neighboring
C     triangle if I2->I1 is a boundary arc.
C
            IF (ABS(LIST(LP)) .NE. I2) GO TO 13
            KN = 0
            IF (LIST(LP) .LT. 0) GO TO 8
C
C   I2->I1 is not a boundary arc, and LP points to I2 as
C     a neighbor of I1.
C
    5       LP = LPTR(LP)
            I3 = ABS(LIST(LP))
C
C Find L such that LTRI(L,KN) = I3 (not used if KN > KT),
C   and permute the vertex indexes of KN so that I1 is
C   smallest.
C
            IF (I1 .LT. I2  .AND.  I1 .LT. I3) THEN
              L = 3
            ELSEIF (I2 .LT. I3) THEN
              L = 2
              ISV = I1
              I1 = I2
              I2 = I3
              I3 = ISV
            ELSE
              L = 1
              ISV = I1
              I1 = I3
              I3 = I2
              I2 = ISV
            ENDIF
C
C Test for KN > KT (triangle index not yet assigned).
C
            IF (I1 .GT. N1  .AND.  .NOT. PASS2) GO TO 9
C
C Find KN, if it exists, by searching the triangle list in
C   reverse order.
C
            DO 6 KN = KT-1,1,-1
              IF (LTRI(1,KN) .EQ. I1  .AND.  LTRI(2,KN) .EQ.
     .            I2  .AND.  LTRI(3,KN) .EQ. I3) GO TO 7
    6         CONTINUE
            GO TO 9
C
C Store KT as a neighbor of KN.
C
    7       LTRI(L+3,KN) = KT
C
C Store KN as a neighbor of KT, and add a new arc KA.
C
    8       LTRI(I+3,KT) = KN
            IF (ARCS) THEN
              KA = KA + 1
              LTRI(I+6,KT) = KA
              IF (KN .NE. 0) LTRI(L+6,KN) = KA
            ENDIF
    9       CONTINUE
C
C Bottom of loop on triangles.
C
   10     IF (LP2 .NE. LPLN1) GO TO 3
   11     CONTINUE
C
C Bottom of loop on nodes.
C
      IF (.NOT. PASS2  .AND.  NCC .GT. 0) THEN
        PASS2 = .TRUE.
        N1ST = LCC1
        GO TO 2
      ENDIF
C
C No errors encountered.
C
      NT = KT
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   12 NT = 0
      IER = 1
      RETURN
C
C Invalid triangulation data structure:  I1 is a neighbor of
C   I2, but I2 is not a neighbor of I1.
C
   13 NT = 0
      IER = 2
      RETURN
      END
      SUBROUTINE TRLPRT (NCC,LCT,N,X,Y,NROW,NT,LTRI,LOUT,
     .                   PRNTX)
      INTEGER NCC, LCT(*), N, NROW, NT, LTRI(NROW,NT),
     .        LOUT
      LOGICAL PRNTX
      REAL    X(N), Y(N)
C
C***********************************************************
C
C                                               From TRLPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/02/98
C
C   Given a triangulation of a set of points in the plane,
C this subroutine prints the triangle list created by
C Subroutine TRLIST and, optionally, the nodal coordinates
C on logical unit LOUT.  The numbers of boundary nodes,
C triangles, and arcs, and the constraint region triangle
C indexes, if any, are also printed.
C
C   All parameters other than LOUT and PRNTX should be
C unaltered from their values on output from TRLIST.
C
C
C On input:
C
C       NCC = Number of constraints.
C
C       LCT = List of constraint triangle starting indexes
C             (or dummy array of length 1 if NCC = 0).
C
C       N = Number of nodes in the triangulation.
C           3 .LE. N .LE. 9999.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes in the triangulation -- not used
C             unless PRNTX = TRUE.
C
C       NROW = Number of rows (entries per triangle) re-
C              served for the triangle list LTRI.  The value
C              must be 6 if only the vertex indexes and
C              neighboring triangle indexes are stored, or 9
C              if arc indexes are also stored.
C
C       NT = Number of triangles in the triangulation.
C            1 .LE. NT .LE. 9999.
C
C       LTRI = NROW by NT array whose J-th column contains
C              the vertex nodal indexes (first three rows),
C              neighboring triangle indexes (second three
C              rows), and, if NROW = 9, arc indexes (last
C              three rows) associated with triangle J for
C              J = 1,...,NT.
C
C       LOUT = Logical unit number for output.  0 .LE. LOUT
C              .LE. 99.  Output is printed on unit 6 if LOUT
C              is outside its valid range on input.
C
C       PRNTX = Logical variable with value TRUE if and only
C               if X and Y are to be printed (to 6 decimal
C               places).
C
C None of the parameters are altered by this routine.
C
C Modules required by TRLPRT:  None
C
C***********************************************************
C
      INTEGER I, K, LUN, NA, NB, NL, NLMAX, NMAX
      DATA    NMAX/9999/,  NLMAX/60/
C
C Local parameters:
C
C   I = DO-loop, nodal index, and row index for LTRI
C   K = DO-loop and triangle index
C   LUN = Logical unit number for output
C   NA = Number of triangulation arcs
C   NB = Number of boundary nodes
C   NL = Number of lines printed on the current page
C   NLMAX = Maximum number of print lines per page
C   NMAX = Maximum value of N and NT (4-digit format)
C
      LUN = LOUT
      IF (LUN .LT. 0  .OR.  LUN .GT. 99) LUN = 6
C
C Print a heading and test for invalid input.
C
      WRITE (LUN,100)
      NL = 1
      IF (N .LT. 3  .OR.  N .GT. NMAX  .OR.
     .    (NROW .NE. 6  .AND.  NROW .NE. 9)  .OR.
     .    NT .LT. 1  .OR.  NT .GT. NMAX) THEN
C
C Print an error message and bypass the loops.
C
        WRITE (LUN,110) N, NROW, NT
        GO TO 3
      ENDIF
      IF (PRNTX) THEN
C
C Print X and Y.
C
        WRITE (LUN,101)
        NL = 6
        DO 1 I = 1,N
          IF (NL .GE. NLMAX) THEN
            WRITE (LUN,106)
            NL = 0
          ENDIF
          WRITE (LUN,102) I, X(I), Y(I)
          NL = NL + 1
    1     CONTINUE
      ENDIF
C
C Print the triangulation LTRI.
C
      IF (NL .GT. NLMAX/2) THEN
        WRITE (LUN,106)
        NL = 0
      ENDIF
      IF (NROW .EQ. 6) THEN
        WRITE (LUN,103)
      ELSE
        WRITE (LUN,104)
      ENDIF
      NL = NL + 5
      DO 2 K = 1,NT
        IF (NL .GE. NLMAX) THEN
          WRITE (LUN,106)
          NL = 0
        ENDIF
        WRITE (LUN,105) K, (LTRI(I,K), I = 1,NROW)
        NL = NL + 1
    2   CONTINUE
C
C Print NB, NA, and NT (boundary nodes, arcs, and
C   triangles).
C
      NB = 2*N - NT - 2
      NA = NT + N - 1
      IF (NL .GT. NLMAX-6) WRITE (LUN,106)
      WRITE (LUN,107) NB, NA, NT
C
C Print NCC and LCT.
C
    3 WRITE (LUN,108) NCC
      IF (NCC .GT. 0) WRITE (LUN,109) (LCT(I), I = 1,NCC)
      RETURN
C
C Print formats:
C
  100 FORMAT (///,24X,'TRIPACK (TRLIST) Output')
  101 FORMAT (//16X,'Node',7X,'X(Node)',10X,'Y(Node)'//)
  102 FORMAT (16X,I4,2E17.6)
  103 FORMAT (//1X,'Triangle',8X,'Vertices',12X,'Neighbors'/
     .        4X,'KT',7X,'N1',5X,'N2',5X,'N3',4X,'KT1',4X,
     .        'KT2',4X,'KT3'/)
  104 FORMAT (//1X,'Triangle',8X,'Vertices',12X,'Neighbors',
     .        14X,'Arcs'/
     .        4X,'KT',7X,'N1',5X,'N2',5X,'N3',4X,'KT1',4X,
     .        'KT2',4X,'KT3',4X,'KA1',4X,'KA2',4X,'KA3'/)
  105 FORMAT (2X,I4,2X,6(3X,I4),3(2X,I5))
  106 FORMAT (///)
  107 FORMAT (/1X,'NB = ',I4,' Boundary Nodes',5X,
     .        'NA = ',I5,' Arcs',5X,'NT = ',I5,
     .        ' Triangles')
  108 FORMAT (/1X,'NCC =',I3,' Constraint Curves')
  109 FORMAT (1X,9X,14I5)
  110 FORMAT (//1X,10X,'*** Invalid Parameter:  N =',I5,
     .        ', NROW =',I5,', NT =',I5,' ***')
      END
      SUBROUTINE TRMESH (N,X,Y, LIST,LPTR,LEND,LNEW,NEAR,
     .                   NEXT,DIST,IER)
      USE SWPCOM, ONLY: SWTOL
      INTEGER N, LIST(*), LPTR(*), LEND(N), LNEW, NEAR(N),
     .        NEXT(N), IER
      REAL    X(N), Y(N), DIST(N)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/28/98
C
C   This subroutine creates a Delaunay triangulation of a
C set of N arbitrarily distributed points in the plane re-
C ferred to as nodes.  The Delaunay triangulation is defined
C as a set of triangles with the following five properties:
C
C  1)  The triangle vertices are nodes.
C  2)  No triangle contains a node other than its vertices.
C  3)  The interiors of the triangles are pairwise disjoint.
C  4)  The union of triangles is the convex hull of the set
C        of nodes (the smallest convex set which contains
C        the nodes).
C  5)  The interior of the circumcircle of each triangle
C        contains no node.
C
C The first four properties define a triangulation, and the
C last property results in a triangulation which is as close
C as possible to equiangular in a certain sense and which is
C uniquely defined unless four or more nodes lie on a common
C circle.  This property makes the triangulation well-suited
C for solving closest point problems and for triangle-based
C interpolation.
C
C   The triangulation can be generalized to a constrained
C Delaunay triangulation by a call to Subroutine ADDCST.
C This allows for user-specified boundaries defining a non-
C convex and/or multiply connected region.
C
C   The algorithm for constructing the triangulation has
C expected time complexity O(N*log(N)) for most nodal dis-
C tributions.  Also, since the algorithm proceeds by adding
C nodes incrementally, the triangulation may be updated with
C the addition (or deletion) of a node very efficiently.
C The adjacency information representing the triangulation
C is stored as a linked list requiring approximately 13N
C storage locations.
C
C
C   The following is a list of the software package modules
C which a user may wish to call directly:
C
C  ADDCST - Generalizes the Delaunay triangulation to allow
C             for user-specified constraints.
C
C  ADDNOD - Updates the triangulation by appending or
C             inserting a new node.
C
C  AREAP  - Computes the area bounded by a closed polygonal
C             curve such as the boundary of the triangula-
C             tion or of a constraint region.
C
C  BNODES - Returns an array containing the indexes of the
C             boundary nodes in counterclockwise order.
C             Counts of boundary nodes, triangles, and arcs
C             are also returned.
C
C  CIRCUM - Computes the area, circumcenter, circumradius,
C             and, optionally, the aspect ratio of a trian-
C             gle defined by user-specified vertices.
C
C  DELARC - Deletes a boundary arc from the triangulation.
C
C  DELNOD - Updates the triangulation with the deletion of a
C             node.
C
C  EDGE   - Forces a pair of nodes to be connected by an arc
C             in the triangulation.
C
C  GETNP  - Determines the ordered sequence of L closest
C             nodes to a given node, along with the associ-
C             ated distances.  The distance between nodes is
C             taken to be the length of the shortest connec-
C             ting path which intersects no constraint
C             region.
C
C  INTSEC - Determines whether or not an arbitrary pair of
C             line segments share a common point.
C
C  JRAND  - Generates a uniformly distributed pseudo-random
C             integer.
C
C  LEFT   - Locates a point relative to a line.
C
C  NEARND - Returns the index of the nearest node to an
C             arbitrary point, along with its squared
C             distance.
C
C  STORE  - Forces a value to be stored in main memory so
C             that the precision of floating point numbers
C             in memory locations rather than registers is
C             computed.
C
C  TRLIST - Converts the triangulation data structure to a
C             triangle list more suitable for use in a fin-
C             ite element code.
C
C  TRLPRT - Prints the triangle list created by Subroutine
C             TRLIST.
C
C  TRMESH - Creates a Delaunay triangulation of a set of
C             nodes.
C
C  TRMSHR - Creates a Delaunay triangulation (more effici-
C             ently than TRMESH) of a set of nodes lying at
C             the vertices of a (possibly skewed) rectangu-
C             lar grid.
C
C  TRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
C             file containing a triangulation plot.
C
C  TRPRNT - Prints the triangulation data structure and,
C             optionally, the nodal coordinates.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.  (X(K),Y(K)) is re-
C             ferred to as node K, and K is referred to as
C             a nodal index.  The first three nodes must not
C             be collinear.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR = Arrays of length at least 6N-12.
C
C       LEND = Array of length at least N.
C
C       NEAR,NEXT,DIST = Work space arrays of length at
C                        least N.  The space is used to
C                        efficiently determine the nearest
C                        triangulation node to each un-
C                        processed node for use by ADDNOD.
C
C On output:
C
C       LIST = Set of nodal indexes which, along with LPTR,
C              LEND, and LNEW, define the triangulation as a
C              set of N adjacency lists -- counterclockwise-
C              ordered sequences of neighboring nodes such
C              that the first and last neighbors of a bound-
C              ary node are boundary nodes (the first neigh-
C              bor of an interior node is arbitrary).  In
C              order to distinguish between interior and
C              boundary nodes, the last neighbor of each
C              boundary node is represented by the negative
C              of its index.
C
C       LPTR = Set of pointers (LIST indexes) in one-to-one
C              correspondence with the elements of LIST.
C              LIST(LPTR(I)) indexes the node which follows
C              LIST(I) in cyclical counterclockwise order
C              (the first neighbor follows the last neigh-
C              bor).
C
C       LEND = Set of pointers to adjacency lists.  LEND(K)
C              points to the last neighbor of node K for
C              K = 1,...,N.  Thus, LIST(LEND(K)) < 0 if and
C              only if K is a boundary node.
C
C       LNEW = Pointer to the first empty location in LIST
C              and LPTR (list length plus one).  LIST, LPTR,
C              LEND, and LNEW are not altered if IER < 0,
C              and are incomplete if IER > 0.
C
C       NEAR,NEXT,DIST = Garbage.
C
C       IER = Error indicator:
C             IER =  0 if no errors were encountered.
C             IER = -1 if N < 3 on input.
C             IER = -2 if the first three nodes are
C                      collinear.
C             IER = -4 if an error flag was returned by a
C                      call to SWAP in ADDNOD.  This is an
C                      internal error and should be reported
C                      to the programmer.
C             IER =  L if nodes L and M coincide for some
C                      M > L.  The linked list represents
C                      a triangulation of nodes 1 to M-1
C                      in this case.
C
C Modules required by TRMESH:  ADDNOD, BDYADD, INSERT,
C                                INTADD, JRAND, LEFT,
C                                LSTPTR, STORE, SWAP,
C                                SWPTST, TRFIND
C
C Intrinsic function called by TRMESH:  ABS
C
C***********************************************************
C
      LOGICAL LEFT
      REAL    STORE
      INTEGER I, I0, J, K, KM1, LCC(1), LP, LPL, NCC, NEXTI,
     .        NN
      REAL    D, D1, D2, D3, EPS
C
C Local parameters:
C
C D =        Squared distance from node K to node I
C D1,D2,D3 = Squared distances from node K to nodes 1, 2,
C              and 3, respectively
C EPS =      Half the machine precision
C I,J =      Nodal indexes
C I0 =       Index of the node preceding I in a sequence of
C              unprocessed nodes:  I = NEXT(I0)
C K =        Index of node to be added and DO-loop index:
C              K > 3
C KM1 =      K-1
C LCC(1) =   Dummy array
C LP =       LIST index (pointer) of a neighbor of K
C LPL =      Pointer to the last neighbor of K
C NCC =      Number of constraint curves
C NEXTI =    NEXT(I)
C NN =       Local copy of N
C SWTOL =    Tolerance for function SWPTST
C
      NN = N
      IF (NN .LT. 3) THEN
        IER = -1
        RETURN
      ENDIF
C
C Compute a tolerance for function SWPTST:  SWTOL = 10*
C   (machine precision)
C
      EPS = 1.
    1 EPS = EPS/2.
        SWTOL = STORE(EPS + 1.)
        IF (SWTOL .GT. 1.) GO TO 1
      SWTOL = EPS*20.
C
C Store the first triangle in the linked list.
C
      IF ( .NOT. LEFT(X(1),Y(1),X(2),Y(2),X(3),Y(3)) ) THEN
C
C   The initial triangle is (3,2,1) = (2,1,3) = (1,3,2).
C
        LIST(1) = 3
        LPTR(1) = 2
        LIST(2) = -2
        LPTR(2) = 1
        LEND(1) = 2
C
        LIST(3) = 1
        LPTR(3) = 4
        LIST(4) = -3
        LPTR(4) = 3
        LEND(2) = 4
C
        LIST(5) = 2
        LPTR(5) = 6
        LIST(6) = -1
        LPTR(6) = 5
        LEND(3) = 6
C
      ELSEIF ( .NOT. LEFT(X(2),Y(2),X(1),Y(1),X(3),Y(3)) )
     .       THEN
C
C   The initial triangle is (1,2,3).
C
        LIST(1) = 2
        LPTR(1) = 2
        LIST(2) = -3
        LPTR(2) = 1
        LEND(1) = 2
C
        LIST(3) = 3
        LPTR(3) = 4
        LIST(4) = -1
        LPTR(4) = 3
        LEND(2) = 4
C
        LIST(5) = 1
        LPTR(5) = 6
        LIST(6) = -2
        LPTR(6) = 5
        LEND(3) = 6
C
      ELSE
C
C   The first three nodes are collinear.
C
        IER = -2
        RETURN
      ENDIF
C
C Initialize LNEW and test for N = 3.
C
      LNEW = 7
      IF (NN .EQ. 3) THEN
        IER = 0
        RETURN
      ENDIF
C
C A nearest-node data structure (NEAR, NEXT, and DIST) is
C   used to obtain an expected-time (N*log(N)) incremental
C   algorithm by enabling constant search time for locating
C   each new node in the triangulation.
C
C For each unprocessed node K, NEAR(K) is the index of the
C   triangulation node closest to K (used as the starting
C   point for the search in Subroutine TRFIND) and DIST(K)
C   is an increasing function of the distance between nodes
C   K and NEAR(K).
C
C Since it is necessary to efficiently find the subset of
C   unprocessed nodes associated with each triangulation
C   node J (those that have J as their NEAR entries), the
C   subsets are stored in NEAR and NEXT as follows:  for
C   each node J in the triangulation, I = NEAR(J) is the
C   first unprocessed node in J's set (with I = 0 if the
C   set is empty), L = NEXT(I) (if I > 0) is the second,
C   NEXT(L) (if L > 0) is the third, etc.  The nodes in each
C   set are initially ordered by increasing indexes (which
C   maximizes efficiency) but that ordering is not main-
C   tained as the data structure is updated.
C
C Initialize the data structure for the single triangle.
C
      NEAR(1) = 0
      NEAR(2) = 0
      NEAR(3) = 0
      DO 2 K = NN,4,-1
        D1 = (X(K)-X(1))**2 + (Y(K)-Y(1))**2
        D2 = (X(K)-X(2))**2 + (Y(K)-Y(2))**2
        D3 = (X(K)-X(3))**2 + (Y(K)-Y(3))**2
        IF (D1 .LE. D2  .AND.  D1 .LE. D3) THEN
          NEAR(K) = 1
          DIST(K) = D1
          NEXT(K) = NEAR(1)
          NEAR(1) = K
        ELSEIF (D2 .LE. D1  .AND.  D2 .LE. D3) THEN
          NEAR(K) = 2
          DIST(K) = D2
          NEXT(K) = NEAR(2)
          NEAR(2) = K
        ELSE
          NEAR(K) = 3
          DIST(K) = D3
          NEXT(K) = NEAR(3)
          NEAR(3) = K
        ENDIF
    2   CONTINUE
C
C Add the remaining nodes.  Parameters for ADDNOD are as
C   follows:
C
C   K = Index of the node to be added.
C   NEAR(K) = Index of the starting node for the search in
C             TRFIND.
C   NCC = Number of constraint curves.
C   LCC = Dummy array (since NCC = 0).
C   KM1 = Number of nodes in the triangulation.
C
      NCC = 0
      DO 7 K = 4,NN
        KM1 = K-1
        CALL ADDNOD (K,X(K),Y(K),NEAR(K),NCC, LCC,KM1,X,Y,
     .               LIST,LPTR,LEND,LNEW, IER)
        IF (IER .NE. 0) RETURN
C
C Remove K from the set of unprocessed nodes associated
C   with NEAR(K).
C
        I = NEAR(K)
        IF (NEAR(I) .EQ. K) THEN
          NEAR(I) = NEXT(K)
        ELSE
          I = NEAR(I)
    3     I0 = I
            I = NEXT(I0)
            IF (I .NE. K) GO TO 3
          NEXT(I0) = NEXT(K)
        ENDIF
        NEAR(K) = 0
C
C Loop on neighbors J of node K.
C
        LPL = LEND(K)
        LP = LPL
    4   LP = LPTR(LP)
          J = ABS(LIST(LP))
C
C Loop on elements I in the sequence of unprocessed nodes
C   associated with J:  K is a candidate for replacing J
C   as the nearest triangulation node to I.  The next value
C   of I in the sequence, NEXT(I), must be saved before I
C   is moved because it is altered by adding I to K's set.
C
          I = NEAR(J)
    5     IF (I .EQ. 0) GO TO 6
          NEXTI = NEXT(I)
C
C Test for the distance from I to K less than the distance
C   from I to J.
C
          D = (X(K)-X(I))**2 + (Y(K)-Y(I))**2
          IF (D .LT. DIST(I)) THEN
C
C Replace J by K as the nearest triangulation node to I:
C   update NEAR(I) and DIST(I), and remove I from J's set
C   of unprocessed nodes and add it to K's set.
C
            NEAR(I) = K
            DIST(I) = D
            IF (I .EQ. NEAR(J)) THEN
              NEAR(J) = NEXTI
            ELSE
              NEXT(I0) = NEXTI
            ENDIF
            NEXT(I) = NEAR(K)
            NEAR(K) = I
          ELSE
            I0 = I
          ENDIF
C
C Bottom of loop on I.
C
          I = NEXTI
          GO TO 5
C
C Bottom of loop on neighbors J.
C
    6     IF (LP .NE. LPL) GO TO 4
    7   CONTINUE
      RETURN
      END
      SUBROUTINE TRMSHR (N,NX,X,Y, NIT, LIST,LPTR,LEND,LNEW,
     .                   IER)
      USE SWPCOM, ONLY: SWTOL
      INTEGER  N, NX, NIT, LIST(*), LPTR(*), LEND(N), LNEW,
     .         IER
      REAL     X(N), Y(N)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/27/98
C
C   This subroutine creates a Delaunay triangulation of a
C set of N nodes in the plane, where the nodes are the vert-
C ices of an NX by NY skewed rectangular grid with the
C natural ordering.  Thus, N = NX*NY, and the nodes are
C ordered from left to right beginning at the top row so
C that adjacent nodes have indexes which differ by 1 in the
C x-direction and by NX in the y-direction.  A skewed rec-
C tangular grid is defined as one in which each grid cell is
C a strictly convex quadrilateral (and is thus the convex
C hull of its four vertices).  Equivalently, any transfor-
C mation from a rectangle to a grid cell which is bilinear
C in both components has an invertible Jacobian.
C
C   If the nodes are not distributed and ordered as defined
C above, Subroutine TRMESH must be called in place of this
C routine.  Refer to Subroutine ADDCST for the treatment of
C constraints.
C
C   The first phase of the algorithm consists of construc-
C ting a triangulation by choosing a diagonal arc in each
C grid cell.  If NIT = 0, all diagonals connect lower left
C to upper right corners and no error checking or additional
C computation is performed.  Otherwise, each diagonal arc is
C chosen to be locally optimal, and boundary arcs are added
C where necessary in order to cover the convex hull of the
C nodes.  (This is the first iteration.)  If NIT > 1 and no
C error was detected, the triangulation is then optimized by
C a sequence of up to NIT-1 iterations in which interior
C arcs of the triangulation are tested and swapped if appro-
C priate.  The algorithm terminates when an iteration
C results in no swaps and/or when the allowable number of
C iterations has been performed.  NIT = 0 is sufficient to
C produce a Delaunay triangulation if the original grid is
C actually rectangular, and NIT = 1 is sufficient if it is
C close to rectangular.  Note, however, that the ordering
C and distribution of nodes is not checked for validity in
C the case NIT = 0, and the triangulation will not be valid
C unless the rectangular grid covers the convex hull of the
C nodes.
C
C
C On input:
C
C       N = Number of nodes in the grid.  N = NX*NY for some
C           NY .GE. 2.
C
C       NX = Number of grid points in the x-direction.  NX
C            .GE. 2.
C
C       X,Y = Arrays of length N containing coordinates of
C             the nodes with the ordering and distribution
C             defined in the header comments above.
C             (X(K),Y(K)) is referred to as node K.
C
C The above parameters are not altered by this routine.
C
C       NIT = Nonnegative integer specifying the maximum
C             number of iterations to be employed.  Refer
C             to the header comments above.
C
C       LIST,LPTR = Arrays of length at least 6N-12.
C
C       LEND = Array of length at least N.
C
C On output:
C
C       NIT = Number of iterations employed.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.  Refer to Sub-
C                             routine TRMESH.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = K if the grid element with upper left
C                     corner at node K is not a strictly
C                     convex quadrilateral.  The algorithm
C                     is terminated when the first such
C                     occurrence is detected.  Note that
C                     this test is not performed if NIT = 0
C                     on input.
C             IER = -1 if N, NX, or NIT is outside its valid
C                      range on input.
C             IER = -2 if NIT > 1 on input, and the optimi-
C                      zation loop failed to converge within
C                      the allowable number of iterations.
C                      The triangulation is valid but not
C                      optimal in this case.
C
C Modules required by TRMSHR:  INSERT, LEFT, LSTPTR, NBCNT,
C                                STORE, SWAP, SWPTST
C
C Intrinsic function called by TRMSHR:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR, NBCNT
      LOGICAL LEFT, SWPTST
      REAL    STORE
      INTEGER I, ITER, J, K, KP1, LP, LPF, LPK, LPL, LPP,
     .        M1, M2, M3, M4, MAXIT, N0, N1, N2, N3, N4, NI,
     .        NJ, NM1, NN, NNB
      LOGICAL TST
      REAL    EPS
C
C Store local variables and test for errors in input
C   parameters.
C
      NI = NX
      NJ = N/NI
      NN = NI*NJ
      MAXIT = NIT
      NIT = 0
      IF (N .NE. NN  .OR.  NJ .LT. 2  .OR.  NI .LT. 2  .OR.
     .    MAXIT .LT. 0) THEN
        IER = -1
        RETURN
      ENDIF
      IER = 0
C
C Compute a tolerance for function SWPTST:  SWTOL = 10*
C   (machine precision)
C
      EPS = 1.
    1 EPS = EPS/2.
        SWTOL = STORE(EPS + 1.)
        IF (SWTOL .GT. 1.) GO TO 1
      SWTOL = EPS*20.
C
C Loop on grid points (I,J) corresponding to nodes K =
C   (J-1)*NI + I.  TST = TRUE iff diagonals are to be
C   chosen by the swap test.  M1, M2, M3, and M4 are the
C   slopes (-1, 0, or 1) of the diagonals in quadrants 1
C   to 4 (counterclockwise beginning with the upper right)
C   for a coordinate system with origin at node K.
C
      TST = MAXIT .GT. 0
      M1 = 0
      M4 = 0
      LP = 0
      KP1 = 1
      DO 6 J = 1,NJ
        DO 5 I = 1,NI
          M2 = M1
          M3 = M4
          K = KP1
          KP1 = K + 1
          LPF = LP + 1
          IF (J .EQ. NJ  .AND.  I .NE. NI) GO TO 2
          IF (I .NE. 1) THEN
            IF (J .NE. 1) THEN
C
C   K is not in the top row, leftmost column, or bottom row
C     (unless K is the lower right corner).  Take the first
C     neighbor to be the node above K.
C
              LP = LP + 1
              LIST(LP) = K - NI
              LPTR(LP) = LP + 1
              IF (M2 .LE. 0) THEN
                LP = LP + 1
                LIST(LP) = K - 1 - NI
                LPTR(LP) = LP + 1
              ENDIF
            ENDIF
C
C   K is not in the leftmost column.  The next (or first)
C     neighbor is to the left of K.
C
            LP = LP + 1
            LIST(LP) = K - 1
            LPTR(LP) = LP + 1
            IF (J .EQ. NJ) GO TO 3
            IF (M3 .GE. 0) THEN
              LP = LP + 1
              LIST(LP) = K - 1 + NI
              LPTR(LP) = LP + 1
            ENDIF
          ENDIF
C
C   K is not in the bottom row.  The next (or first)
C     neighbor is below K.
C
          LP = LP + 1
          LIST(LP) = K + NI
          LPTR(LP) = LP + 1
C
C   Test for a negative diagonal in quadrant 4 unless K is
C     in the rightmost column.  The quadrilateral associated
C     with the quadrant is tested for strict convexity un-
C     less NIT = 0 on input.
C
          IF (I .EQ. NI) GO TO 3
          M4 = 1
          IF (.NOT. TST) GO TO 2
          IF ( LEFT(X(KP1),Y(KP1),X(K+NI),Y(K+NI),X(K),Y(K))
     .         .OR.  LEFT(X(K),Y(K),X(KP1+NI),Y(KP1+NI),
     .                    X(K+NI),Y(K+NI))
     .         .OR.  LEFT(X(K+NI),Y(K+NI),X(KP1),Y(KP1),
     .                    X(KP1+NI),Y(KP1+NI))
     .         .OR.  LEFT(X(KP1+NI),Y(KP1+NI),X(K),Y(K),
     .                    X(KP1),Y(KP1)) )          GO TO 12
          IF ( SWPTST(KP1,K+NI,K,KP1+NI,X,Y) ) GO TO 2
          M4 = -1
          LP = LP + 1
          LIST(LP) = KP1 + NI
          LPTR(LP) = LP + 1
C
C   The next (or first) neighbor is to the right of K.
C
    2     LP = LP + 1
          LIST(LP) = KP1
          LPTR(LP) = LP + 1
C
C   Test for a positive diagonal in quadrant 1 (the neighbor
C     of K-NI which follows K is not K+1) unless K is in the
C     top row.
C
          IF (J .EQ. 1) GO TO 3
          IF (TST) THEN
            M1 = -1
            LPK = LSTPTR(LEND(K-NI),K,LIST,LPTR)
            LPK = LPTR(LPK)
            IF (LIST(LPK) .NE. KP1) THEN
              M1 = 1
              LP = LP + 1
              LIST(LP) = KP1 - NI
              LPTR(LP) = LP + 1
            ENDIF
          ENDIF
C
C   If K is in the leftmost column (and not the top row) or
C     in the bottom row (and not the rightmost column), then
C     the next neighbor is the node above K.
C
          IF (I .NE. 1  .AND.  J .NE. NJ) GO TO 4
          LP = LP + 1
          LIST(LP) = K - NI
          LPTR(LP) = LP + 1
          IF (I .EQ. 1) GO TO 3
C
C   K is on the bottom row (and not the leftmost or right-
C     most column).
C
          IF (M2 .LE. 0) THEN
            LP = LP + 1
            LIST(LP) = K - 1 - NI
            LPTR(LP) = LP + 1
          ENDIF
          LP = LP + 1
          LIST(LP) = K - 1
          LPTR(LP) = LP + 1
C
C   K is a boundary node.
C
    3     LIST(LP) = -LIST(LP)
C
C   Bottom of loop.  Store LEND and correct LPTR(LP).
C     LPF and LP point to the first and last neighbors
C     of K.
C
    4     LEND(K) = LP
          LPTR(LP) = LPF
    5     CONTINUE
    6   CONTINUE
C
C Store LNEW, and terminate the algorithm if NIT = 0 on
C   input.
C
      LNEW = LP + 1
      IF (MAXIT .EQ. 0) RETURN
C
C Add boundary arcs where necessary in order to cover the
C   convex hull of the nodes.  N1, N2, and N3 are consecu-
C   tive boundary nodes in counterclockwise order, and N0
C   is the starting point for each loop around the boundary.
C
      N0 = 1
      N1 = N0
      N2 = NI + 1
C
C   TST is set to TRUE if an arc is added.  The boundary
C     loop is repeated until a traversal results in no
C     added arcs.
C
    7 TST = .FALSE.
C
C   Top of boundary loop.  Set N3 to the first neighbor of
C     N2, and test for N3 LEFT N1 -> N2.
C
    8   LPL = LEND(N2)
          LP = LPTR(LPL)
          N3 = LIST(LP)
          IF ( LEFT(X(N1),Y(N1),X(N2),Y(N2),X(N3),Y(N3)) )
     .       N1 = N2
          IF (N1 .NE. N2) THEN
C
C   Add the boundary arc N1-N3.  If N0 = N2, the starting
C     point is changed to N3, since N2 will be removed from
C     the boundary.  N3 is inserted as the first neighbor of
C     N1, N2 is changed to an interior node, and N1 is
C     inserted as the last neighbor of N3.
C
            TST = .TRUE.
            IF (N2 .EQ. N0) N0 = N3
            LP = LEND(N1)
            CALL INSERT (N3,LP, LIST,LPTR,LNEW )
            LIST(LPL) = -LIST(LPL)
            LP = LEND(N3)
            LIST(LP) = N2
            CALL INSERT (-N1,LP, LIST,LPTR,LNEW )
            LEND(N3) = LNEW - 1
          ENDIF
C
C   Bottom of loops.  Test for termination.
C
          N2 = N3
          IF (N1 .NE. N0) GO TO 8
        IF (TST) GO TO 7
C
C Terminate the algorithm if NIT = 1 on input.
C
      NIT = 1
      IF (MAXIT .EQ. 1) RETURN
C
C Optimize the triangulation by applying the swap test and
C   appropriate swaps to the interior arcs.  The loop is
C   repeated until no swaps are performed or MAXIT itera-
C   tions have been applied.  ITER is the current iteration,
C   and TST is set to TRUE if a swap occurs.
C
      ITER = 1
      NM1 = NN - 1
    9 ITER = ITER + 1
        TST = .FALSE.
C
C   Loop on interior arcs N1-N2, where N2 > N1 and
C     (N1,N2,N3) and (N2,N1,N4) are adjacent triangles.
C
C   Top of loop on nodes N1.
C
        DO 11 N1 = 1,NM1
          LPL = LEND(N1)
          N4 = LIST(LPL)
          LPF = LPTR(LPL)
          N2 = LIST(LPF)
          LP = LPTR(LPF)
          N3 = LIST(LP)
          NNB = NBCNT(LPL,LPTR)
C
C   Top of loop on neighbors N2 of N1.  NNB is the number of
C                                       neighbors of N1.
C
          DO 10 I = 1,NNB
C
C   Bypass the swap test if N1 is a boundary node and N2 is
C     the first neighbor (N4 < 0), N2 < N1, or N1-N2 is a
C     diagonal arc (already locally optimal) when ITER = 2.
C
            IF ( N4 .GT. 0  .AND.  N2 .GT. N1  .AND.
     .          (ITER .NE. 2  .OR.  ABS(N1+NI-N2) .NE. 1) )
     .          THEN
              IF (SWPTST(N3,N4,N1,N2,X,Y) ) THEN
C
C   Swap diagonal N1-N2 for N3-N4, set TST to TRUE, and set
C     N2 to N4 (the neighbor preceding N3).
C
                CALL SWAP (N3,N4,N1,N2, LIST,LPTR,LEND, LPP)
                IF (LPP .NE. 0) THEN
                  TST = .TRUE.
                  N2 = N4
                ENDIF
              ENDIF
            ENDIF
C
C   Bottom of neighbor loop.
C
            IF (LIST(LPL) .EQ. -N3) GO TO 11
            N4 = N2
            N2 = N3
            LP = LSTPTR(LPL,N2,LIST,LPTR)
            LP = LPTR(LP)
            N3 = ABS(LIST(LP))
   10       CONTINUE
   11     CONTINUE
C
C   Test for termination.
C
        IF (TST  .AND.  ITER .LT. MAXIT) GO TO 9
      NIT = ITER
      IF (TST) IER = -2
      RETURN
C
C Invalid grid cell encountered.
C
   12 IER = K
      RETURN
      END
      SUBROUTINE TRPLOT (LUN,PLTSIZ,WX1,WX2,WY1,WY2,NCC,LCC,
     .                   N,X,Y,LIST,LPTR,LEND,TITLE,
     .                   NUMBR, IER)
      CHARACTER*(*) TITLE
      INTEGER LUN, NCC, LCC(*), N, LIST(*), LPTR(*),
     .        LEND(N), IER
      LOGICAL NUMBR
      REAL    PLTSIZ, WX1, WX2, WY1, WY2, X(N), Y(N)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/15/98
C
C   This subroutine creates a level-2 Encapsulated Post-
C script (EPS) file containing a triangulation plot.
C
C
C On input:
C
C       LUN = Logical unit number in the range 0 to 99.
C             The unit should be opened with an appropriate
C             file name before the call to this routine.
C
C       PLTSIZ = Plot size in inches.  The window is mapped,
C                with aspect ratio preserved, to a rectangu-
C                lar viewport with maximum side-length equal
C                to .88*PLTSIZ (leaving room for labels out-
C                side the viewport).  The viewport is
C                centered on the 8.5 by 11 inch page, and
C                its boundary is drawn.  1.0 .LE. PLTSIZ
C                .LE. 8.5.
C
C       WX1,WX2,WY1,WY2 = Parameters defining a rectangular
C                         window against which the triangu-
C                         lation is clipped.  (Only the
C                         portion of the triangulation that
C                         lies in the window is drawn.)
C                         (WX1,WY1) and (WX2,WY2) are the
C                         lower left and upper right cor-
C                         ners, respectively.  WX1 < WX2 and
C                         WY1 < WY2.
C
C       NCC = Number of constraint curves.  Refer to Subrou-
C             tine ADDCST.  NCC .GE. 0.
C
C       LCC = Array of length NCC (or dummy parameter if
C             NCC = 0) containing the index of the first
C             node of constraint I in LCC(I).  For I = 1 to
C             NCC, LCC(I+1)-LCC(I) .GE. 3, where LCC(NCC+1)
C             = N+1.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes with non-constraint nodes in the
C             first LCC(1)-1 locations.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C       TITLE = Type CHARACTER variable or constant contain-
C               ing a string to be centered above the plot.
C               The string must be enclosed in parentheses;
C               i.e., the first and last characters must be
C               '(' and ')', respectively, but these are not
C               displayed.  TITLE may have at most 80 char-
C               acters including the parentheses.
C
C       NUMBR = Option indicator:  If NUMBR = TRUE, the
C               nodal indexes are plotted next to the nodes.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LUN, PLTSIZ, NCC, or N is outside
C                     its valid range.  LCC is not tested
C                     for validity.
C             IER = 2 if WX1 >= WX2 or WY1 >= WY2.
C             IER = 3 if an error was encountered in writing
C                     to unit LUN.
C
C   Various plotting options can be controlled by altering
C the data statement below.
C
C Modules required by TRPLOT:  None
C
C Intrinsic functions called by TRPLOT:  ABS, CHAR, NINT,
C                                          REAL
C
C***********************************************************
C
      INTEGER I, IFRST, IH, ILAST, IPX1, IPX2, IPY1, IPY2,
     .        IW, LP, LPL, N0, N0BAK, N0FOR, N1, NLS
      LOGICAL ANNOT, CNSTR, PASS1
      REAL    DASHL, DX, DY, FSIZN, FSIZT, R, SFX, SFY, T,
     .        TX, TY, X0, Y0
C
      DATA    ANNOT/.TRUE./,  DASHL/4.0/,  FSIZN/10.0/,
     .        FSIZT/16.0/
C
C Local parameters:
C
C ANNOT =     Logical variable with value TRUE iff the plot
C               is to be annotated with the values of WX1,
C               WX2, WY1, and WY2
C CNSTR       Logical variable used to flag constraint arcs:
C               TRUE iff N0-N1 lies in a constraint region
C DASHL =     Length (in points, at 72 points per inch) of
C               dashes and spaces in a dashed line pattern
C               used for drawing constraint arcs
C DX =        Window width WX2-WX1
C DY =        Window height WY2-WY1
C FSIZN =     Font size in points for labeling nodes with
C               their indexes if NUMBR = TRUE
C FSIZT =     Font size in points for the title (and
C               annotation if ANNOT = TRUE)
C I =         Constraint index (1 to NCC)
C IFRST =     Index of the first node in constraint I
C IH =        Height of the viewport in points
C ILAST =     Index of the last node in constraint I
C IPX1,IPY1 = X and y coordinates (in points) of the lower
C               left corner of the bounding box or viewport
C IPX2,IPY2 = X and y coordinates (in points) of the upper
C               right corner of the bounding box or viewport
C IW =        Width of the viewport in points
C LP =        LIST index (pointer)
C LPL =       Pointer to the last neighbor of N0
C N0 =        Nodal index and DO-loop index
C N0BAK =     Predecessor of N0 in a constraint curve
C               (sequence of adjacent constraint nodes)
C N0FOR =     Successor to N0 in a constraint curve
C N1 =        Index of a neighbor of N0
C NLS =       Index of the last non-constraint node
C PASS1 =     Logical variable used to flag the first pass
C               through the constraint nodes
C R =         Aspect ratio DX/DY
C SFX,SFY =   Scale factors for mapping world coordinates
C               (window coordinates in [WX1,WX2] X [WY1,WY2])
C               to viewport coordinates in [IPX1,IPX2] X
C               [IPY1,IPY2]
C T =         Temporary variable
C TX,TY =     Translation vector for mapping world coordi-
C               nates to viewport coordinates
C X0,Y0 =     X(N0),Y(N0) or label location
C
C
C Test for error 1, and set NLS to the last non-constraint
C   node.
C
      IF (LUN .LT. 0  .OR.  LUN .GT. 99  .OR.
     .    PLTSIZ .LT. 1.0  .OR.  PLTSIZ .GT. 8.5  .OR.
     .    NCC .LT. 0  .OR.  N .LT. 3) GO TO 11
      NLS = N
      IF (NCC .GT. 0) NLS = LCC(1)-1
C
C Compute the aspect ratio of the window.
C
      DX = WX2 - WX1
      DY = WY2 - WY1
      IF (DX .LE. 0.0  .OR.  DY .LE. 0.0) GO TO 12
      R = DX/DY
C
C Compute the lower left (IPX1,IPY1) and upper right
C   (IPX2,IPY2) corner coordinates of the bounding box.
C   The coordinates, specified in default user space units
C   (points, at 72 points/inch with origin at the lower
C   left corner of the page), are chosen to preserve the
C   aspect ratio R, and to center the plot on the 8.5 by 11
C   inch page.  The center of the page is (306,396), and
C   T = PLTSIZ/2 in points.
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
      WRITE (LUN,100,ERR=13) IPX1, IPY1, IPX2, IPY2
  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/
     .        '%%BoundingBox:',4I4/
     .        '%%Title:  Triangulation'/
     .        '%%Creator:  TRIPACK'/
     .        '%%EndComments')
C
C Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates
C   of a viewport obtained by shrinking the bounding box by
C   12% in each dimension.
C
      IW = NINT(0.88*REAL(IPX2-IPX1))
      IH = NINT(0.88*REAL(IPY2-IPY1))
      IPX1 = 306 - IW/2
      IPX2 = 306 + IW/2
      IPY1 = 396 - IH/2
      IPY2 = 396 + IH/2
C
C Set the line thickness to 2 points, and draw the
C   viewport boundary.
C
      T = 2.0
      WRITE (LUN,110,ERR=13) T
      WRITE (LUN,120,ERR=13) IPX1, IPY1
      WRITE (LUN,130,ERR=13) IPX1, IPY2
      WRITE (LUN,130,ERR=13) IPX2, IPY2
      WRITE (LUN,130,ERR=13) IPX2, IPY1
      WRITE (LUN,140,ERR=13)
      WRITE (LUN,150,ERR=13)
  110 FORMAT (F12.6,' setlinewidth')
  120 FORMAT (2I4,' moveto')
  130 FORMAT (2I4,' lineto')
  140 FORMAT ('closepath')
  150 FORMAT ('stroke')
C
C Set up a mapping from the window to the viewport.
C
      SFX = REAL(IW)/DX
      SFY = REAL(IH)/DY
      TX = IPX1 - SFX*WX1
      TY = IPY1 - SFY*WY1
      WRITE (LUN,160,ERR=13) TX, TY, SFX, SFY
  160 FORMAT (2F12.6,' translate'/
     .        2F12.6,' scale')
C
C The line thickness (believe it or fucking not) must be
C   changed to reflect the new scaling which is applied to
C   all subsequent output.  Set it to 1.0 point.
C
      T = 2.0/(SFX+SFY)
      WRITE (LUN,110,ERR=13) T
C
C Save the current graphics state, and set the clip path to
C   the boundary of the window.
C
      WRITE (LUN,170,ERR=13)
      WRITE (LUN,180,ERR=13) WX1, WY1
      WRITE (LUN,190,ERR=13) WX2, WY1
      WRITE (LUN,190,ERR=13) WX2, WY2
      WRITE (LUN,190,ERR=13) WX1, WY2
      WRITE (LUN,200,ERR=13)
  170 FORMAT ('gsave')
  180 FORMAT (2F12.6,' moveto')
  190 FORMAT (2F12.6,' lineto')
  200 FORMAT ('closepath clip newpath')
C
C Draw the edges N0->N1, where N1 > N0, beginning with a
C   loop on non-constraint nodes N0.  LPL points to the
C   last neighbor of N0.
C
      DO 3 N0 = 1,NLS
        X0 = X(N0)
        Y0 = Y(N0)
        LPL = LEND(N0)
        LP = LPL
C
C   Loop on neighbors N1 of N0.
C
    2   LP = LPTR(LP)
          N1 = ABS(LIST(LP))
          IF (N1 .GT. N0) THEN
C
C   Add the edge to the path.
C
            WRITE (LUN,210,ERR=13) X0, Y0, X(N1), Y(N1)
  210       FORMAT (2F12.6,' moveto',2F12.6,' lineto')
          ENDIF
          IF (LP .NE. LPL) GO TO 2
    3   CONTINUE
C
C Loop through the constraint nodes twice.  The non-
C   constraint arcs incident on constraint nodes are
C   drawn (with solid lines) on the first pass, and the
C   constraint arcs (both boundary and interior, if any)
C   are drawn (with dashed lines) on the second pass.
C
      PASS1 = .TRUE.
C
C Loop on constraint nodes N0 with (N0BAK,N0,N0FOR) a sub-
C   sequence of constraint I.  The outer loop is on
C   constraints I with first and last nodes IFRST and ILAST.
C
    4 IFRST = N+1
      DO 8 I = NCC,1,-1
        ILAST = IFRST - 1
        IFRST = LCC(I)
        N0BAK = ILAST
        DO 7 N0 = IFRST,ILAST
          N0FOR = N0 + 1
          IF (N0 .EQ. ILAST) N0FOR = IFRST
          LPL = LEND(N0)
          X0 = X(N0)
          Y0 = Y(N0)
          LP = LPL
C
C   Loop on neighbors N1 of N0.  CNSTR = TRUE iff N0-N1 is a
C     constraint arc.
C
C   Initialize CNSTR to TRUE iff the first neighbor of N0
C     strictly follows N0FOR and precedes or coincides with
C     N0BAK (in counterclockwise order).
C
    5     LP = LPTR(LP)
            N1 = ABS(LIST(LP))
            IF (N1 .NE. N0FOR  .AND.  N1 .NE. N0BAK) GO TO 5
          CNSTR = N1 .EQ. N0BAK
          LP = LPL
C
C   Loop on neighbors N1 of N0.  Update CNSTR and test for
C     N1 > N0.
C
    6     LP = LPTR(LP)
            N1 = ABS(LIST(LP))
            IF (N1 .EQ. N0FOR) CNSTR = .TRUE.
            IF (N1 .GT. N0) THEN
C
C   Draw the edge iff (PASS1=TRUE and CNSTR=FALSE) or
C     (PASS1=FALSE and CNSTR=TRUE); i.e., CNSTR and PASS1
C     have opposite values.
C
              IF (CNSTR .NEQV. PASS1)
     .          WRITE (LUN,210,ERR=13) X0, Y0, X(N1), Y(N1)
            ENDIF
            IF (N1 .EQ. N0BAK) CNSTR = .FALSE.
C
C   Bottom of loops.
C
            IF (LP .NE. LPL) GO TO 6
          N0BAK = N0
    7     CONTINUE
    8   CONTINUE
      IF (PASS1) THEN
C
C End of first pass:  paint the path and change to dashed
C   lines for subsequent drawing.  Since the scale factors
C   are applied to everything, the dash length must be
C   specified in world coordinates.
C
        PASS1 = .FALSE.
        WRITE (LUN,150,ERR=13)
        T = DASHL*2.0/(SFX+SFY)
        WRITE (LUN,220,ERR=13) T
  220   FORMAT ('[',F12.6,'] 0 setdash')
        GO TO 4
      ENDIF
C
C Paint the path and restore the saved graphics state (with
C   no clip path).
C
      WRITE (LUN,150,ERR=13)
      WRITE (LUN,230,ERR=13)
  230 FORMAT ('grestore')
      IF (NUMBR) THEN
C
C Nodes in the window are to be labeled with their indexes.
C   Convert FSIZN from points to world coordinates, and
C   output the commands to select a font and scale it.
C
        T = FSIZN*2.0/(SFX+SFY)
        WRITE (LUN,240,ERR=13) T
  240   FORMAT ('/Helvetica findfont'/
     .          F12.6,' scalefont setfont')
C
C   Loop on nodes N0 with coordinates (X0,Y0).
C
        DO 9 N0 = 1,N
          X0 = X(N0)
          Y0 = Y(N0)
          IF (X0 .LT. WX1  .OR.  X0 .GT. WX2  .OR.
     .        Y0 .LT. WY1  .OR.  Y0 .GT. WY2) GO TO 9
C
C   Move to (X0,Y0), and draw the label N0.  The first char-
C     acter will have its lower left corner about one
C     character width to the right of the nodal position.
C
          WRITE (LUN,180,ERR=13) X0, Y0
          WRITE (LUN,250,ERR=13) N0
  250     FORMAT ('(',I3,') show')
    9     CONTINUE
      ENDIF
C
C Convert FSIZT from points to world coordinates, and output
C   the commands to select a font and scale it.
C
      T = FSIZT*2.0/(SFX+SFY)
      WRITE (LUN,240,ERR=13) T
C
C Display TITLE centered above the plot:
C
      Y0 = WY2 + 3.0*T
      WRITE (LUN,260,ERR=13) TITLE, (WX1+WX2)/2.0, Y0
  260 FORMAT (A80/'  stringwidth pop 2 div neg ',F12.6,
     .        ' add ',F12.6,' moveto')
      WRITE (LUN,270,ERR=13) TITLE
  270 FORMAT (A80/'  show')
      IF (ANNOT) THEN
C
C Display the window extrema below the plot.
C
        X0 = WX1
        Y0 = WY1 - 100.0/(SFX+SFY)
        WRITE (LUN,180,ERR=13) X0, Y0
        WRITE (LUN,280,ERR=13) WX1, WX2
        Y0 = Y0 - 2.0*T
        WRITE (LUN,290,ERR=13) X0, Y0, WY1, WY2
  280   FORMAT ('(Window:   WX1 = ',E9.3,',   WX2 = ',E9.3,
     .          ') show')
  290   FORMAT ('(Window:  ) stringwidth pop ',F12.6,' add',
     .          F12.6,' moveto'/
     .          '( WY1 = ',E9.3,',   WY2 = ',E9.3,') show')
      ENDIF
C
C Paint the path and output the showpage command and
C   end-of-file indicator.
C
      WRITE (LUN,300,ERR=13)
  300 FORMAT ('stroke'/
     .        'showpage'/
     .        '%%EOF')
C
C HP's interpreters require a one-byte End-of-PostScript-Job
C   indicator (to eliminate a timeout error message):
C   ASCII 4.
C
      WRITE (LUN,310,ERR=13) CHAR(4)
  310 FORMAT (A1)
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
C DX or DY is not positive.
C
   12 IER = 2
      RETURN
C
C Error writing to unit LUN.
C
   13 IER = 3
      RETURN
      END
      SUBROUTINE TRPRNT (NCC,LCC,N,X,Y,LIST,LPTR,LEND,LOUT,
     .                   PRNTX)
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     .        LOUT
      LOGICAL PRNTX
      REAL    X(N), Y(N)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/30/98
C
C   Given a triangulation of a set of points in the plane,
C this subroutine prints the adjacency lists and, option-
C ally, the nodal coordinates on logical unit LOUT.  The
C list of neighbors of a boundary node is followed by index
C 0.  The numbers of boundary nodes, triangles, and arcs,
C and the constraint curve starting indexes, if any, are
C also printed.
C
C
C On input:
C
C       NCC = Number of constraints.
C
C       LCC = List of constraint curve starting indexes (or
C             dummy array of length 1 if NCC = 0).
C
C       N = Number of nodes in the triangulation.
C           3 .LE. N .LE. 9999.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes in the triangulation -- not used
C             unless PRNTX = TRUE.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C       LOUT = Logical unit number for output.  0 .LE. LOUT
C              .LE. 99.  Output is printed on unit 6 if LOUT
C              is outside its valid range on input.
C
C       PRNTX = Logical variable with value TRUE if and only
C               if X and Y are to be printed (to 6 decimal
C               places).
C
C None of the parameters are altered by this routine.
C
C Modules required by TRPRNT:  None
C
C***********************************************************
C
      INTEGER I, INC, K, LP, LPL, LUN, NA, NABOR(100), NB,
     .        ND, NL, NLMAX, NMAX, NODE, NN, NT
      DATA  NMAX/9999/,  NLMAX/60/
C
      NN = N
      LUN = LOUT
      IF (LUN .LT. 0  .OR.  LUN .GT. 99) LUN = 6
C
C Print a heading and test the range of N.
C
      WRITE (LUN,100) NN
      IF (NN .LT. 3  .OR.  NN .GT. NMAX) THEN
C
C N is outside its valid range.
C
        WRITE (LUN,110)
        GO TO 5
      ENDIF
C
C Initialize NL (the number of lines printed on the current
C   page) and NB (the number of boundary nodes encountered).
C
      NL = 6
      NB = 0
      IF (.NOT. PRNTX) THEN
C
C Print LIST only.  K is the number of neighbors of NODE
C   which are stored in NABOR.
C
        WRITE (LUN,101)
        DO 2 NODE = 1,NN
          LPL = LEND(NODE)
          LP = LPL
          K = 0
C
    1     K = K + 1
            LP = LPTR(LP)
            ND = LIST(LP)
            NABOR(K) = ND
            IF (LP .NE. LPL) GO TO 1
          IF (ND .LE. 0) THEN
C
C   NODE is a boundary node.  Correct the sign of the last
C     neighbor, add 0 to the end of the list, and increment
C     NB.
C
            NABOR(K) = -ND
            K = K + 1
            NABOR(K) = 0
            NB = NB + 1
          ENDIF
C
C   Increment NL and print the list of neighbors.
C
          INC = (K-1)/14 + 2
          NL = NL + INC
          IF (NL .GT. NLMAX) THEN
            WRITE (LUN,106)
            NL = INC
          ENDIF
          WRITE (LUN,103) NODE, (NABOR(I), I = 1,K)
          IF (K .NE. 14) WRITE (LUN,105)
    2     CONTINUE
      ELSE
C
C Print X, Y, and LIST.
C
        WRITE (LUN,102)
        DO 4 NODE = 1,NN
          LPL = LEND(NODE)
          LP = LPL
          K = 0
    3     K = K + 1
            LP = LPTR(LP)
            ND = LIST(LP)
            NABOR(K) = ND
            IF (LP .NE. LPL) GO TO 3
          IF (ND .LE. 0) THEN
C
C   NODE is a boundary node.
C
            NABOR(K) = -ND
            K = K + 1
            NABOR(K) = 0
            NB = NB + 1
          ENDIF
C
C   Increment NL and print X, Y, and NABOR.
C
          INC = (K-1)/8 + 2
          NL = NL + INC
          IF (NL .GT. NLMAX) THEN
            WRITE (LUN,106)
            NL = INC
          ENDIF
          WRITE (LUN,104) NODE, X(NODE), Y(NODE),
     .                    (NABOR(I), I = 1,K)
          IF (K .NE. 8) WRITE (LUN,105)
    4     CONTINUE
      ENDIF
C
C Print NB, NA, and NT (boundary nodes, arcs, and
C   triangles).
C
      NT = 2*NN - NB - 2
      NA = NT + NN - 1
      IF (NL .GT. NLMAX-6) WRITE (LUN,106)
      WRITE (LUN,107) NB, NA, NT
C
C Print NCC and LCC.
C
    5 WRITE (LUN,108) NCC
      IF (NCC .GT. 0) WRITE (LUN,109) (LCC(I), I = 1,NCC)
      RETURN
C
C Print formats:
C
  100 FORMAT (///,26X,'Adjacency Sets,    N = ',I5//)
  101 FORMAT (1X,'Node',32X,'Neighbors of Node'//)
  102 FORMAT (1X,'Node',5X,'X(Node)',8X,'Y(Node)',
     .        20X,'Neighbors of Node'//)
  103 FORMAT (1X,I4,5X,14I5/(1X,9X,14I5))
  104 FORMAT (1X,I4,2E15.6,5X,8I5/(1X,39X,8I5))
  105 FORMAT (1X)
  106 FORMAT (///)
  107 FORMAT (/1X,'NB = ',I4,' Boundary Nodes',5X,
     .        'NA = ',I5,' Arcs',5X,'NT = ',I5,
     .        ' Triangles')
  108 FORMAT (/1X,'NCC =',I3,' Constraint Curves')
  109 FORMAT (1X,9X,14I5)
  110 FORMAT (1X,10X,'*** N is outside its valid',
     .        ' range ***')
      END
