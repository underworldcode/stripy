subroutine addnod ( nst, k, x, y, z, list, lptr, lend, lnew, ier )

!*****************************************************************************80
!
!! ADDNOD adds a node to a triangulation.
!
!  Discussion:
!
!    This subroutine adds node K to a triangulation of the
!    convex hull of nodes 1, ..., K-1, producing a triangulation
!    of the convex hull of nodes 1, ..., K.
!
!    The algorithm consists of the following steps:  node K
!    is located relative to the triangulation (TRFIND), its
!    index is added to the data structure (INTADD or BDYADD),
!    and a sequence of swaps (SWPTST and SWAP) are applied to
!    the arcs opposite K so that all arcs incident on node K
!    and opposite node K are locally optimal (satisfy the circumcircle test).
!
!    Thus, if a Delaunay triangulation of nodes 1 through K-1 is input,
!    a Delaunay triangulation of nodes 1 through K will be output.
!
!  Modified:
!
!    15 May 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, the index of a node at which TRFIND
!    begins its search.  Search time depends on the proximity of this node to
!    K.  If NST < 1, the search is begun at node K-1.
!
!    Input, integer ( kind = 4 ) K, the nodal index (index for X, Y, Z, and
!    LEND) of the new node to be added.  4 <= K.
!
!    Input, real ( kind = 8 ) X(K), Y(K), Z(K), the coordinates of the nodes.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(K),
!    LNEW.  On input, the data structure associated with the triangulation of
!    nodes 1 to K-1.  On output, the data has been updated to include node
!    K.  The array lengths are assumed to be large enough to add node K.
!    Refer to TRMESH.
!
!    Output, integer ( kind = 4 ) IER, error indicator:
!     0 if no errors were encountered.
!    -1 if K is outside its valid range on input.
!    -2 if all nodes (including K) are collinear (lie on a common geodesic).
!     L if nodes L and K coincide for some L < K.
!
!  Local parameters:
!
!    B1,B2,B3 = Unnormalized barycentric coordinates returned by TRFIND.
!    I1,I2,I3 = Vertex indexes of a triangle containing K
!    IN1 =      Vertex opposite K:  first neighbor of IO2
!               that precedes IO1.  IN1,IO1,IO2 are in
!               counterclockwise order.
!    IO1,IO2 =  Adjacent neighbors of K defining an arc to
!               be tested for a swap
!    IST =      Index of node at which TRFIND begins its search
!    KK =       Local copy of K
!    KM1 =      K-1
!    L =        Vertex index (I1, I2, or I3) returned in IER
!               if node K coincides with a vertex
!    LP =       LIST pointer
!    LPF =      LIST pointer to the first neighbor of K
!    LPO1 =     LIST pointer to IO1
!    LPO1S =    Saved value of LPO1
!    P =        Cartesian coordinates of node K
!
  implicit none

  integer ( kind = 4 ) k

  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) ist
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lend(k)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lpo1
  integer ( kind = 4 ) lpo1s
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) nst
  real ( kind = 8 ) p(3)
  logical swptst
  real ( kind = 8 ) x(k)
  real ( kind = 8 ) y(k)
  real ( kind = 8 ) z(k)

  kk = k

  if ( kk < 4 ) then
    ier = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADDNOD - Fatal error!'
    write ( *, '(a)' ) '  K < 4.'
    return
  end if
!
!  Initialization:
!
  km1 = kk - 1
  ist = nst
  if ( ist < 1 ) then
    ist = km1
  end if

  p(1) = x(kk)
  p(2) = y(kk)
  p(3) = z(kk)
!
!  Find a triangle (I1,I2,I3) containing K or the rightmost
!  (I1) and leftmost (I2) visible boundary nodes as viewed
!  from node K.
!
  call trfind ( ist, p, km1, x, y, z, list, lptr, lend, b1, b2, b3, &
    i1, i2, i3 )
!
!  Test for collinear or duplicate nodes.
!
  if ( i1 == 0 ) then
    ier = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADDNOD - Fatal error!'
    write ( *, '(a)' ) '  The nodes are coplanar.'
    return
  end if

  if ( i3 /= 0 ) then

    l = i1

    if ( p(1) == x(l) .and. p(2) == y(l)  .and. p(3) == z(l) ) then
      ier = l
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ADDNOD - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Node ', l, ' is equal to node ', k
      return
    end if

    l = i2

    if ( p(1) == x(l) .and. p(2) == y(l)  .and. p(3) == z(l) ) then
      ier = l
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ADDNOD - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Node ', l, ' is equal to node ', k
      return
    end if

    l = i3
    if ( p(1) == x(l) .and. p(2) == y(l)  .and. p(3) == z(l) ) then
      ier = l
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ADDNOD - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Node ', l, ' is equal to node ', k
      return
    end if

    call intadd ( kk, i1, i2, i3, list, lptr, lend, lnew )

  else

    if ( i1 /= i2 ) then
      call bdyadd ( kk, i1,i2, list, lptr, lend, lnew )
    else
      call covsph ( kk, i1, list, lptr, lend, lnew )
    end if

  end if

  ier = 0
!
!  Initialize variables for optimization of the triangulation.
!
  lp = lend(kk)
  lpf = lptr(lp)
  io2 = list(lpf)
  lpo1 = lptr(lpf)
  io1 = abs ( list(lpo1) )
!
!  Begin loop: find the node opposite K.
!
  do

    lp = lstptr ( lend(io1), io2, list, lptr )

    if ( 0 <= list(lp) ) then

      lp = lptr(lp)
      in1 = abs ( list(lp) )
!
!  Swap test:  if a swap occurs, two new arcs are
!  opposite K and must be tested.
!
      lpo1s = lpo1

      if ( .not. swptst ( in1, kk, io1, io2, x, y, z ) ) then

        if ( lpo1 == lpf .or. list(lpo1) < 0 ) then
          exit
        end if

        io2 = io1
        lpo1 = lptr(lpo1)
        io1 = abs ( list(lpo1) )
        cycle

      end if

      call swap ( in1, kk, io1, io2, list, lptr, lend, lpo1 )
!
!  A swap is not possible because KK and IN1 are already
!  adjacent.  This error in SWPTST only occurs in the
!  neutral case and when there are nearly duplicate nodes.
!
      if ( lpo1 /= 0 ) then
        io1 = in1
        cycle
      end if

      lpo1 = lpo1s

    end if
!
!  No swap occurred.  Test for termination and reset IO2 and IO1.
!
    if ( lpo1 == lpf .or. list(lpo1) < 0 ) then
      exit
    end if

    io2 = io1
    lpo1 = lptr(lpo1)
    io1 = abs ( list(lpo1) )

  end do

  return
end
function arc_cosine ( c )

!*****************************************************************************80
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
!    This routine truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, the argument.
!
!    Output, real ( kind = 8 ) ARC_COSINE, an angle whose cosine is C.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) c
  real ( kind = 8 ) c2

  c2 = c
  c2 = max ( c2, -1.0D+00 )
  c2 = min ( c2, +1.0D+00 )

  arc_cosine = acos ( c2 )

  return
end



function areas ( v1, v2, v3 )

!*****************************************************************************80
!
!! AREAS computes the area of a spherical triangle on the unit sphere.
!
!  Discussion:
!
!    This function returns the area of a spherical triangle
!    on the unit sphere.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the Cartesian coordinates
!    of unit vectors (the three triangle vertices in any order).  These
!    vectors, if nonzero, are implicitly scaled to have length 1.
!
!    Output, real ( kind = 8 ) AREAS, the area of the spherical triangle
!    defined by V1, V2, and V3, in the range 0 to 2*PI (the area of a
!    hemisphere).  AREAS = 0 (or 2*PI) if and only if V1, V2, and V3 lie in (or
!    close to) a plane containing the origin.
!
!  Local parameters:
!
!    A1,A2,A3 =    Interior angles of the spherical triangle.
!
!    CA1,CA2,CA3 = cos(A1), cos(A2), and cos(A3), respectively.
!
!    DV1,DV2,DV3 = copies of V1, V2, and V3.
!
!    I =           DO-loop index and index for Uij.
!
!    S12,S23,S31 = Sum of squared components of U12, U23, U31.
!
!    U12,U23,U31 = Unit normal vectors to the planes defined by
!                 pairs of triangle vertices.
!
  implicit none

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) a3
  real ( kind = 8 ) areas
  real ( kind = 8 ) ca1
  real ( kind = 8 ) ca2
  real ( kind = 8 ) ca3
  real ( kind = 8 ) dv1(3)
  real ( kind = 8 ) dv2(3)
  real ( kind = 8 ) dv3(3)
  real ( kind = 8 ) s12
  real ( kind = 8 ) s23
  real ( kind = 8 ) s31
  real ( kind = 8 ) u12(3)
  real ( kind = 8 ) u23(3)
  real ( kind = 8 ) u31(3)
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)

  dv1(1:3) = v1(1:3)
  dv2(1:3) = v2(1:3)
  dv3(1:3) = v3(1:3)
!
!  Compute cross products Uij = Vi X Vj.
!
  u12(1) = dv1(2) * dv2(3) - dv1(3) * dv2(2)
  u12(2) = dv1(3) * dv2(1) - dv1(1) * dv2(3)
  u12(3) = dv1(1) * dv2(2) - dv1(2) * dv2(1)

  u23(1) = dv2(2) * dv3(3) - dv2(3) * dv3(2)
  u23(2) = dv2(3) * dv3(1) - dv2(1) * dv3(3)
  u23(3) = dv2(1) * dv3(2) - dv2(2) * dv3(1)

  u31(1) = dv3(2) * dv1(3) - dv3(3) * dv1(2)
  u31(2) = dv3(3) * dv1(1) - dv3(1) * dv1(3)
  u31(3) = dv3(1) * dv1(2) - dv3(2) * dv1(1)
!
!  Normalize Uij to unit vectors.
!
  s12 = dot_product ( u12(1:3), u12(1:3) )
  s23 = dot_product ( u23(1:3), u23(1:3) )
  s31 = dot_product ( u31(1:3), u31(1:3) )
!
!  Test for a degenerate triangle associated with collinear vertices.
!
  if ( s12 == 0.0D+00 .or. s23 == 0.0D+00 .or. s31 == 0.0D+00 ) then
    areas = 0.0D+00
    return
  end if

  s12 = sqrt ( s12 )
  s23 = sqrt ( s23 )
  s31 = sqrt ( s31 )

  u12(1:3) = u12(1:3) / s12
  u23(1:3) = u23(1:3) / s23
  u31(1:3) = u31(1:3) / s31
!
!  Compute interior angles Ai as the dihedral angles between planes:
!  CA1 = cos(A1) = -<U12,U31>
!  CA2 = cos(A2) = -<U23,U12>
!  CA3 = cos(A3) = -<U31,U23>
!
  ca1 = - dot_product ( u12(1:3), u31(1:3) )
  ca2 = - dot_product ( u23(1:3), u12(1:3) )
  ca3 = - dot_product ( u31(1:3), u23(1:3) )

  ca1 = max ( ca1, -1.0D+00 )
  ca1 = min ( ca1, +1.0D+00 )
  ca2 = max ( ca2, -1.0D+00 )
  ca2 = min ( ca2, +1.0D+00 )
  ca3 = max ( ca3, -1.0D+00 )
  ca3 = min ( ca3, +1.0D+00 )

  a1 = acos ( ca1 )
  a2 = acos ( ca2 )
  a3 = acos ( ca3 )
!
!  Compute AREAS = A1 + A2 + A3 - PI.
!
  areas = a1 + a2 + a3 - acos ( -1.0D+00 )

  if ( areas < 0.0D+00 ) then
    areas = 0.0D+00
  end if

  return
end



subroutine triareas ( n, x, y, z, nt, triangles, triarea )
!*****************************************************************************
!! TRIAREAS computes the area of every triangle within a triangulation of points
!
! Discussion:
!
!   This subroutine computes the area of each triangle
!   returned by TRMESH and TRLIST. It uses the AREA function
!   to compute the area between 3 points on a sphere.
!
! Parameters:
!
!   Input, integer ( kind = 4 ), n
!   number of points in the triangulation
!
!   Input, real ( kind = 8 ), x(n), y(n), z(n)
!   points that make up the triangulation
!
!   Input, integer ( kind = 4 ), nt
!   number of triangles returned by TRLIST or TRLIST2
!
!   Input, integer ( kind = 8 ), triangles(3,nt)
!   triangles that constitute the triangulation
!
!   Output, real ( kind = 8 ), triarea(nt)
!   area of the each triangle
!
! Local parameters:
!
!   v1, v2, v3 =  x,y,z coordinates for 3 points of a triangle
!   areas      =  area of a single triangle
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) triangles(3,nt)
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
  integer ( kind = 4 ) tri(3)
  real ( kind = 8 ) triarea(nt)
  real ( kind = 8 ) areas
  integer ( kind = 4 ) i

  do i = 1, nt
    tri = triangles(:,i)
    v1(1) = x(tri(1))
    v1(2) = y(tri(1))
    v1(3) = z(tri(1))

    v2(1) = x(tri(2))
    v2(2) = y(tri(2))
    v2(3) = z(tri(2))

    v3(1) = x(tri(3))
    v3(2) = y(tri(3))
    v3(3) = z(tri(3))

    triarea(i) = areas( v1, v2, v3 )
  end do

  return
end


subroutine bdyadd ( kk, i1, i2, list, lptr, lend, lnew )

!*****************************************************************************80
!
!! BDYADD adds a boundary node to a triangulation.
!
!  Discussion:
!
!    This subroutine adds a boundary node to a triangulation
!    of a set of KK-1 points on the unit sphere.  The data
!    structure is updated with the insertion of node KK, but no
!    optimization is performed.
!
!    This routine is identical to the similarly named routine
!    in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KK, the index of a node to be connected to
!    the sequence of all visible boundary nodes.  1 <= KK and
!    KK must not be equal to I1 or I2.
!
!    Input, integer ( kind = 4 ) I1, the first (rightmost as viewed from KK)
!    boundary node in the triangulation that is visible from
!    node KK (the line segment KK-I1 intersects no arcs.
!
!    Input, integer ( kind = 4 ) I2, the last (leftmost) boundary node that
!    is visible from node KK.  I1 and I2 may be determined by TRFIND.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    LNEW, the triangulation data structure created by TRMESH.
!    Nodes I1 and I2 must be included
!    in the triangulation.  On output, the data structure is updated with
!    the addition of node KK.  Node KK is connected to I1, I2, and
!    all boundary nodes in between.
!
!  Local parameters:
!
!    K =     Local copy of KK
!    LP =    LIST pointer
!    LSAV =  LIST pointer
!    N1,N2 = Local copies of I1 and I2, respectively
!    NEXT =  Boundary node visible from K
!    NSAV =  Boundary node visible from K
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lsav
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nsav

  k = kk
  n1 = i1
  n2 = i2
!
!  Add K as the last neighbor of N1.
!
  lp = lend(n1)
  lsav = lptr(lp)
  lptr(lp) = lnew
  list(lnew) = -k
  lptr(lnew) = lsav
  lend(n1) = lnew
  lnew = lnew + 1
  next = -list(lp)
  list(lp) = next
  nsav = next
!
!  Loop on the remaining boundary nodes between N1 and N2,
!  adding K as the first neighbor.
!
  do

    lp = lend(next)
    call insert ( k, lp, list, lptr, lnew )

    if ( next == n2 ) then
      exit
    end if

    next = -list(lp)
    list(lp) = next

  end do
!
!  Add the boundary nodes between N1 and N2 as neighbors of node K.
!
  lsav = lnew
  list(lnew) = n1
  lptr(lnew) = lnew + 1
  lnew = lnew + 1
  next = nsav

  do

    if ( next == n2 ) then
      exit
    end if

    list(lnew) = next
    lptr(lnew) = lnew + 1
    lnew = lnew + 1
    lp = lend(next)
    next = list(lp)

  end do

  list(lnew) = -n2
  lptr(lnew) = lsav
  lend(k) = lnew
  lnew = lnew + 1

  return
end
subroutine bnodes ( n, list, lptr, lend, nodes, nb, na, nt )

!*****************************************************************************80
!
!! BNODES returns the boundary nodes of a triangulation.
!
!  Discussion:
!
!    Given a triangulation of N nodes on the unit sphere created by TRMESH,
!    this subroutine returns an array containing the indexes (if any) of
!    the counterclockwise sequence of boundary nodes, that is, the nodes on
!    the boundary of the convex hull of the set of nodes.  The
!    boundary is empty if the nodes do not lie in a single
!    hemisphere.  The numbers of boundary nodes, arcs, and
!    triangles are also returned.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), the
!    data structure defining the triangulation, created by TRMESH.
!
!    Output, integer ( kind = 4 ) NODES(*), the ordered sequence of NB boundary
!    node indexes in the range 1 to N.  For safety, the dimension of NODES
!    should be N.
!
!    Output, integer ( kind = 4 ) NB, the number of boundary nodes.
!
!    Output, integer ( kind = 4 ) NA, NT, the number of arcs and triangles,
!    respectively, in the triangulation.
!
!  Local parameters:
!
!    K =   NODES index
!    LP =  LIST pointer
!    N0 =  Boundary node to be added to NODES
!    NN =  Local copy of N
!    NST = First element of nodes (arbitrarily chosen to be
!          the one with smallest index)
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nodes(*)
  integer ( kind = 4 ) nst
  integer ( kind = 4 ) nt

  nn = n
!
!  Search for a boundary node.
!
  nst = 0

  do i = 1, nn

    lp = lend(i)

    if ( list(lp) < 0 ) then
      nst = i
      exit
    end if

  end do
!
!  The triangulation contains no boundary nodes.
!
  if ( nst == 0 ) then
    nb = 0
    na = 3 * ( nn - 2 )
    nt = 2 * ( nn - 2 )
    return
  end if
!
!  NST is the first boundary node encountered.
!
!  Initialize for traversal of the boundary.
!
  nodes(1) = nst
  k = 1
  n0 = nst
!
!  Traverse the boundary in counterclockwise order.
!
  do

    lp = lend(n0)
    lp = lptr(lp)
    n0 = list(lp)

    if ( n0 == nst ) then
      exit
    end if

    k = k + 1
    nodes(k) = n0

  end do
!
!  Store the counts.
!
  nb = k
  nt = 2 * n - nb - 2
  na = nt + n - 1

  return
end
subroutine circum ( v1, v2, v3, c, ier )

!*****************************************************************************80
!
!! CIRCUM returns the circumcenter of a spherical triangle.
!
!  Discussion:
!
!    This subroutine returns the circumcenter of a spherical triangle on the
!    unit sphere:  the point on the sphere surface that is equally distant
!    from the three triangle vertices and lies in the same hemisphere, where
!    distance is taken to be arc-length on the sphere surface.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the coordinates of the
!    three triangle vertices (unit vectors) in counter clockwise order.
!
!    Output, real ( kind = 8 ) C(3), the coordinates of the circumcenter unless
!    0 < IER, in which case C is not defined.  C = (V2-V1) X (V3-V1)
!    normalized to a unit vector.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if V1, V2, and V3 lie on a common line:  (V2-V1) X (V3-V1) = 0.
!
!  Local parameters:
!
!    CNORM = Norm of CU:  used to compute C
!    CU =    Scalar multiple of C:  E1 X E2
!    E1,E2 = Edges of the underlying planar triangle:
!            V2-V1 and V3-V1, respectively
!    I =     DO-loop index
!
  implicit none

  real ( kind = 8 ) c(3)
  real ( kind = 8 ) cnorm
  real ( kind = 8 ) cu(3)
  real ( kind = 8 ) e1(3)
  real ( kind = 8 ) e2(3)
  integer ( kind = 4 ) ier
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)

  ier = 0

  e1(1:3) = v2(1:3) - v1(1:3)
  e2(1:3) = v3(1:3) - v1(1:3)
!
!  Compute CU = E1 X E2 and CNORM**2.
!
  cu(1) = e1(2) * e2(3) - e1(3) * e2(2)
  cu(2) = e1(3) * e2(1) - e1(1) * e2(3)
  cu(3) = e1(1) * e2(2) - e1(2) * e2(1)

  cnorm = sqrt ( sum ( cu(1:3)**2 ) )
!
!  The vertices lie on a common line if and only if CU is the zero vector.
!
  if ( cnorm == 0.0D+00 ) then
    ier = 1
    return
  end if

  c(1:3) = cu(1:3) / cnorm

  return
end
subroutine covsph ( kk, n0, list, lptr, lend, lnew )

!*****************************************************************************80
!
!! COVSPH connects an exterior node to boundary nodes, covering the sphere.
!
!  Discussion:
!
!    This subroutine connects an exterior node KK to all
!    boundary nodes of a triangulation of KK-1 points on the
!    unit sphere, producing a triangulation that covers the
!    sphere.  The data structure is updated with the addition
!    of node KK, but no optimization is performed.  All
!    boundary nodes must be visible from node KK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KK = Index of the node to be connected to the
!    set of all boundary nodes.  4 <= KK.
!
!    Input, integer ( kind = 4 ) N0 = Index of a boundary node (in the range
!    1 to KK-1).  N0 may be determined by TRFIND.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    LNEW, the triangulation data structure created by TRMESH.  Node N0 must
!    be included in the triangulation.  On output, updated with the addition
!    of node KK as the last entry.  The updated triangulation contains no
!    boundary nodes.
!
!  Local parameters:
!
!    K =     Local copy of KK
!    LP =    LIST pointer
!    LSAV =  LIST pointer
!    NEXT =  Boundary node visible from K
!    NST =   Local copy of N0
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lsav
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nst

  k = kk
  nst = n0
!
!  Traverse the boundary in clockwise order, inserting K as
!  the first neighbor of each boundary node, and converting
!  the boundary node to an interior node.
!
  next = nst

  do

    lp = lend(next)
    call insert ( k, lp, list, lptr, lnew )
    next = -list(lp)
    list(lp) = next

    if ( next == nst ) then
      exit
    end if

  end do
!
!  Traverse the boundary again, adding each node to K's adjacency list.
!
  lsav = lnew

  do

    lp = lend(next)
    list(lnew) = next
    lptr(lnew) = lnew + 1
    lnew = lnew + 1
    next = list(lp)

    if ( next == nst ) then
      exit
    end if

  end do

  lptr(lnew-1) = lsav
  lend(k) = lnew - 1

  return
end

subroutine crlist ( n, ncol, x, y, z, list, lend, lptr, lnew, &
  ltri, listc, nb, xc, yc, zc, rc, ier )

!*****************************************************************************80
!
!! CRLIST returns triangle circumcenters and other information.
!
!  Discussion:
!
!    Given a Delaunay triangulation of nodes on the surface
!    of the unit sphere, this subroutine returns the set of
!    triangle circumcenters corresponding to Voronoi vertices,
!    along with the circumradii and a list of triangle indexes
!    LISTC stored in one-to-one correspondence with LIST/LPTR
!    entries.
!
!    A triangle circumcenter is the point (unit vector) lying
!    at the same angular distance from the three vertices and
!    contained in the same hemisphere as the vertices.  (Note
!    that the negative of a circumcenter is also equidistant
!    from the vertices.)  If the triangulation covers the
!    surface, the Voronoi vertices are the circumcenters of the
!    triangles in the Delaunay triangulation.  LPTR, LEND, and
!    LNEW are not altered in this case.
!
!    On the other hand, if the nodes are contained in a
!    single hemisphere, the triangulation is implicitly extended
!    to the entire surface by adding pseudo-arcs (of length
!    greater than 180 degrees) between boundary nodes forming
!    pseudo-triangles whose 'circumcenters' are included in the
!    list.  This extension to the triangulation actually
!    consists of a triangulation of the set of boundary nodes in
!    which the swap test is reversed (a non-empty circumcircle
!    test).  The negative circumcenters are stored as the
!    pseudo-triangle 'circumcenters'.  LISTC, LPTR, LEND, and
!    LNEW contain a data structure corresponding to the
!    extended triangulation (Voronoi diagram), but LIST is not
!    altered in this case.  Thus, if it is necessary to retain
!    the original (unextended) triangulation data structure,
!    copies of LPTR and LNEW must be saved before calling this
!    routine.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.  Note that, if N = 3, there are only two Voronoi vertices
!    separated by 180 degrees, and the Voronoi regions are not well defined.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns reserved for LTRI.
!    This must be at least NB-2, where NB is the number of boundary nodes.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of the nodes
!    (unit vectors).
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), the set of adjacency lists.
!    Refer to TRMESH.
!
!    Input, integer ( kind = 4 ) LEND(N), the set of pointers to ends of
!    adjacency lists.  Refer to TRMESH.
!
!    Input/output, integer ( kind = 4 ) LPTR(6*(N-2)), pointers associated
!    with LIST.  Refer to TRMESH.  On output, pointers associated with LISTC.
!    Updated for the addition of pseudo-triangles if the original triangulation
!    contains boundary nodes (0 < NB).
!
!    Input/output, integer ( kind = 4 ) LNEW.  On input, a pointer to the first
!    empty location in LIST and LPTR (list length plus one).  On output,
!    pointer to the first empty location in LISTC and LPTR (list length plus
!    one).  LNEW is not altered if NB = 0.
!
!    Output, integer ( kind = 4 ) LTRI(6,NCOL).  Triangle list whose first NB-2
!    columns contain the indexes of a clockwise-ordered sequence of vertices
!    (first three rows) followed by the LTRI column indexes of the triangles
!    opposite the vertices (or 0 denoting the exterior region) in the last
!    three rows. This array is not generally of any further use outside this
!    routine.
!
!    Output, integer ( kind = 4 ) LISTC(3*NT), where NT = 2*N-4 is the number
!    of triangles in the triangulation (after extending it to cover the entire
!    surface if necessary).  Contains the triangle indexes (indexes to XC, YC,
!    ZC, and RC) stored in 1-1 correspondence with LIST/LPTR entries (or entries
!    that would be stored in LIST for the extended triangulation):  the index
!    of triangle (N1,N2,N3) is stored in LISTC(K), LISTC(L), and LISTC(M),
!    where LIST(K), LIST(L), and LIST(M) are the indexes of N2 as a neighbor
!    of N1, N3 as a neighbor of N2, and N1 as a neighbor of N3.  The Voronoi
!    region associated with a node is defined by the CCW-ordered sequence of
!    circumcenters in one-to-one correspondence with its adjacency
!    list (in the extended triangulation).
!
!    Output, integer ( kind = 4 ) NB, the number of boundary nodes unless
!    IER = 1.
!
!    Output, real ( kind = 8 ) XC(2*N-4), YC(2*N-4), ZC(2*N-4), the coordinates
!    of the triangle circumcenters (Voronoi vertices).  XC(I)**2 + YC(I)**2
!    + ZC(I)**2 = 1.  The first NB-2 entries correspond to pseudo-triangles
!    if 0 < NB.
!
!    Output, real ( kind = 8 ) RC(2*N-4), the circumradii (the arc lengths or
!    angles between the circumcenters and associated triangle vertices) in
!    1-1 correspondence with circumcenters.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if N < 3.
!    2, if NCOL < NB-2.
!    3, if a triangle is degenerate (has vertices lying on a common geodesic).
!
!  Local parameters:
!
!    C =         Circumcenter returned by Subroutine CIRCUM
!    I1,I2,I3 =  Permutation of (1,2,3):  LTRI row indexes
!    I4 =        LTRI row index in the range 1 to 3
!    IERR =      Error flag for calls to CIRCUM
!    KT =        Triangle index
!    KT1,KT2 =   Indexes of a pair of adjacent pseudo-triangles
!    KT11,KT12 = Indexes of the pseudo-triangles opposite N1
!                and N2 as vertices of KT1
!    KT21,KT22 = Indexes of the pseudo-triangles opposite N1
!                and N2 as vertices of KT2
!    LP,LPN =    LIST pointers
!    LPL =       LIST pointer of the last neighbor of N1
!    N0 =        Index of the first boundary node (initial
!                value of N1) in the loop on boundary nodes
!                used to store the pseudo-triangle indexes
!                in LISTC
!    N1,N2,N3 =  Nodal indexes defining a triangle (CCW order)
!                or pseudo-triangle (clockwise order)
!    N4 =        Index of the node opposite N2 -> N1
!    NM2 =       N-2
!    NN =        Local copy of N
!    NT =        Number of pseudo-triangles:  NB-2
!    SWP =       Logical variable set to TRUE in each optimization
!                loop (loop on pseudo-arcs) iff a swap is performed.
!
!    V1,V2,V3 =  Vertices of triangle KT = (N1,N2,N3) sent to subroutine
!                CIRCUM
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncol

  real ( kind = 8 ) c(3)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) kt
  integer ( kind = 4 ) kt1
  integer ( kind = 4 ) kt11
  integer ( kind = 4 ) kt12
  integer ( kind = 4 ) kt2
  integer ( kind = 4 ) kt21
  integer ( kind = 4 ) kt22
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) listc(6*(n-2))
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpn
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) ltri(6,ncol)
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nm2
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nt
  real ( kind = 8 ) rc(2*n-4)
  logical swp
  logical swptst
  real ( kind = 8 ) t
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xc(2*n-4)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yc(2*n-4)
  real ( kind = 8 ) z(n)
  real ( kind = 8 ) zc(2*n-4)

  nn = n
  nb = 0
  nt = 0

  if ( nn < 3 ) then
    ier = 1
    return
  end if
!
!  Search for a boundary node N1.
!
  lp = 0

  do n1 = 1, nn

    if ( list(lend(n1)) < 0 ) then
      lp = lend(n1)
      exit
    end if

  end do
!
!  Does the triangulation already cover the sphere?
!
  if ( lp /= 0 ) then
!
!  There are 3 <= NB boundary nodes.  Add NB-2 pseudo-triangles (N1,N2,N3)
!  by connecting N3 to the NB-3 boundary nodes to which it is not
!  already adjacent.
!
!  Set N3 and N2 to the first and last neighbors,
!  respectively, of N1.
!
    n2 = -list(lp)
    lp = lptr(lp)
    n3 = list(lp)
!
!  Loop on boundary arcs N1 -> N2 in clockwise order,
!  storing triangles (N1,N2,N3) in column NT of LTRI
!  along with the indexes of the triangles opposite
!  the vertices.
!
    do

      nt = nt + 1

      if ( nt <= ncol ) then
        ltri(1,nt) = n1
        ltri(2,nt) = n2
        ltri(3,nt) = n3
        ltri(4,nt) = nt + 1
        ltri(5,nt) = nt - 1
        ltri(6,nt) = 0
      end if

      n1 = n2
      lp = lend(n1)
      n2 = -list(lp)

      if ( n2 == n3 ) then
        exit
      end if

    end do

    nb = nt + 2

    if ( ncol < nt ) then
      ier = 2
      return
    end if

    ltri(4,nt) = 0
!
!  Optimize the exterior triangulation (set of pseudo-
!  triangles) by applying swaps to the pseudo-arcs N1-N2
!  (pairs of adjacent pseudo-triangles KT1 and KT1 < KT2).
!  The loop on pseudo-arcs is repeated until no swaps are
!  performed.
!
    if ( nt /= 1 ) then

      do

        swp = .false.

        do kt1 = 1, nt - 1

          do i3 = 1, 3

            kt2 = ltri(i3+3,kt1)

            if ( kt2 <= kt1 ) then
              cycle
            end if
!
!  The LTRI row indexes (I1,I2,I3) of triangle KT1 =
!  (N1,N2,N3) are a cyclical permutation of (1,2,3).
!
            if ( i3 == 1 ) then
              i1 = 2
              i2 = 3
            else if ( i3 == 2 ) then
              i1 = 3
              i2 = 1
            else
              i1 = 1
              i2 = 2
            end if

            n1 = ltri(i1,kt1)
            n2 = ltri(i2,kt1)
            n3 = ltri(i3,kt1)
!
!  KT2 = (N2,N1,N4) for N4 = LTRI(I,KT2), where LTRI(I+3,KT2) = KT1.
!
            if ( ltri(4,kt2) == kt1 ) then
              i4 = 1
            else if ( ltri(5,kt2 ) == kt1 ) then
              i4 = 2
            else
              i4 = 3
            end if

            n4 = ltri(i4,kt2)
!
!  The empty circumcircle test is reversed for the pseudo-
!  triangles.  The reversal is implicit in the clockwise
!  ordering of the vertices.
!
            if ( .not. swptst ( n1, n2, n3, n4, x, y, z ) ) then
              cycle
            end if
!
!  Swap arc N1-N2 for N3-N4.  KTij is the triangle opposite
!  Nj as a vertex of KTi.
!
            swp = .true.
            kt11 = ltri(i1+3,kt1)
            kt12 = ltri(i2+3,kt1)

            if ( i4 == 1 ) then
              i2 = 2
              i1 = 3
            else if ( i4 == 2 ) then
              i2 = 3
              i1 = 1
            else
              i2 = 1
              i1 = 2
            end if

            kt21 = ltri(i1+3,kt2)
            kt22 = ltri(i2+3,kt2)
            ltri(1,kt1) = n4
            ltri(2,kt1) = n3
            ltri(3,kt1) = n1
            ltri(4,kt1) = kt12
            ltri(5,kt1) = kt22
            ltri(6,kt1) = kt2
            ltri(1,kt2) = n3
            ltri(2,kt2) = n4
            ltri(3,kt2) = n2
            ltri(4,kt2) = kt21
            ltri(5,kt2) = kt11
            ltri(6,kt2) = kt1
!
!  Correct the KT11 and KT22 entries that changed.
!
            if ( kt11 /= 0 ) then
              i4 = 4
              if ( ltri(4,kt11) /= kt1 ) then
                i4 = 5
                if ( ltri(5,kt11) /= kt1 ) i4 = 6
              end if
              ltri(i4,kt11) = kt2
            end if

            if ( kt22 /= 0 ) then
              i4 = 4
              if ( ltri(4,kt22) /= kt2 ) then
                i4 = 5
                if ( ltri(5,kt22) /= kt2 ) then
                  i4 = 6
                end if
              end if
              ltri(i4,kt22) = kt1
            end if

          end do

        end do

        if ( .not. swp ) then
          exit
        end if

      end do

    end if
!
!  Compute and store the negative circumcenters and radii of
!  the pseudo-triangles in the first NT positions.
!
    do kt = 1, nt

      n1 = ltri(1,kt)
      n2 = ltri(2,kt)
      n3 = ltri(3,kt)
      v1(1) = x(n1)
      v1(2) = y(n1)
      v1(3) = z(n1)
      v2(1) = x(n2)
      v2(2) = y(n2)
      v2(3) = z(n2)
      v3(1) = x(n3)
      v3(2) = y(n3)
      v3(3) = z(n3)

      call circum ( v1, v2, v3, c, ierr )

      if ( ierr /= 0 ) then
        ier = 3
        return
      end if
!
!  Store the negative circumcenter and radius (computed from <V1,C>).
!
      xc(kt) = c(1)
      yc(kt) = c(2)
      zc(kt) = c(3)

      t = dot_product ( v1(1:3), c(1:3) )
      t = max ( t, -1.0D+00 )
      t = min ( t, +1.0D+00 )

      rc(kt) = acos ( t )

    end do

  end if
!
!  Compute and store the circumcenters and radii of the
!  actual triangles in positions KT = NT+1, NT+2, ...
!
!  Also, store the triangle indexes KT in the appropriate LISTC positions.
!
  kt = nt
!
!  Loop on nodes N1.
!
  nm2 = nn - 2

  do n1 = 1, nm2

    lpl = lend(n1)
    lp = lpl
    n3 = list(lp)
!
!  Loop on adjacent neighbors N2,N3 of N1 for which N1 < N2 and N1 < N3.
!
    do

      lp = lptr(lp)
      n2 = n3
      n3 = abs ( list(lp) )

      if ( n1 < n2 .and. n1 < n3 ) then

        kt = kt + 1
!
!  Compute the circumcenter C of triangle KT = (N1,N2,N3).
!
        v1(1) = x(n1)
        v1(2) = y(n1)
        v1(3) = z(n1)
        v2(1) = x(n2)
        v2(2) = y(n2)
        v2(3) = z(n2)
        v3(1) = x(n3)
        v3(2) = y(n3)
        v3(3) = z(n3)

        call circum ( v1, v2, v3, c, ierr )

        if ( ierr /= 0 ) then
          ier = 3
          return
        end if
!
!  Store the circumcenter, radius and triangle index.
!
        xc(kt) = c(1)
        yc(kt) = c(2)
        zc(kt) = c(3)

        t = dot_product ( v1(1:3), c(1:3) )
        t = max ( t, -1.0D+00 )
        t = min ( t, +1.0D+00 )

        rc(kt) = acos ( t )
!
!  Store KT in LISTC(LPN), where abs ( LIST(LPN) ) is the
!  index of N2 as a neighbor of N1, N3 as a neighbor
!  of N2, and N1 as a neighbor of N3.
!
        lpn = lstptr ( lpl, n2, list, lptr )
        listc(lpn) = kt
        lpn = lstptr ( lend(n2), n3, list, lptr )
        listc(lpn) = kt
        lpn = lstptr ( lend(n3), n1, list, lptr )
        listc(lpn) = kt

      end if

      if ( lp == lpl ) then
        exit
      end if

    end do

  end do

  if ( nt == 0 ) then
    ier = 0
    return
  end if
!
!  Store the first NT triangle indexes in LISTC.
!
!  Find a boundary triangle KT1 = (N1,N2,N3) with a boundary arc opposite N3.
!
  kt1 = 0

  do

    kt1 = kt1 + 1

    if ( ltri(4,kt1) == 0 ) then
      i1 = 2
      i2 = 3
      i3 = 1
      exit
    else if ( ltri(5,kt1) == 0 ) then
      i1 = 3
      i2 = 1
      i3 = 2
      exit
    else if ( ltri(6,kt1) == 0 ) then
      i1 = 1
      i2 = 2
      i3 = 3
      exit
    end if

  end do

  n1 = ltri(i1,kt1)
  n0 = n1
!
!  Loop on boundary nodes N1 in CCW order, storing the
!  indexes of the clockwise-ordered sequence of triangles
!  that contain N1.  The first triangle overwrites the
!  last neighbor position, and the remaining triangles,
!  if any, are appended to N1's adjacency list.
!
!  A pointer to the first neighbor of N1 is saved in LPN.
!
  do

    lp = lend(n1)
    lpn = lptr(lp)
    listc(lp) = kt1
!
!  Loop on triangles KT2 containing N1.
!
    do

      kt2 = ltri(i2+3,kt1)

      if ( kt2 == 0 ) then
        exit
      end if
!
!  Append KT2 to N1's triangle list.
!
      lptr(lp) = lnew
      lp = lnew
      listc(lp) = kt2
      lnew = lnew + 1
!
!  Set KT1 to KT2 and update (I1,I2,I3) such that LTRI(I1,KT1) = N1.
!
      kt1 = kt2

      if ( ltri(1,kt1) == n1 ) then
        i1 = 1
        i2 = 2
        i3 = 3
      else if ( ltri(2,kt1) == n1 ) then
        i1 = 2
        i2 = 3
        i3 = 1
      else
        i1 = 3
        i2 = 1
        i3 = 2
      end if

    end do
!
!  Store the saved first-triangle pointer in LPTR(LP), set
!  N1 to the next boundary node, test for termination,
!  and permute the indexes:  the last triangle containing
!  a boundary node is the first triangle containing the
!  next boundary node.
!
    lptr(lp) = lpn
    n1 = ltri(i3,kt1)

    if ( n1 == n0 ) then
      exit
    end if

    i4 = i3
    i3 = i2
    i2 = i1
    i1 = i4

  end do

  ier = 0

  return
end
subroutine delarc ( n, io1, io2, list, lptr, lend, lnew, ier )

!*****************************************************************************80
!
!! DELARC deletes a boundary arc from a triangulation.
!
!  Discussion:
!
!    This subroutine deletes a boundary arc from a triangulation
!    It may be used to remove a null triangle from the
!    convex hull boundary.  Note, however, that if the union of
!    triangles is rendered nonconvex, subroutines DELNOD, EDGE,
!    and TRFIND (and hence ADDNOD) may fail.  Also, function
!    NEARND should not be called following an arc deletion.
!
!    This routine is identical to the similarly named routine in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    4 <= N.
!
!    Input, integer ( kind = 4 ) IO1, IO2, indexes (in the range 1 to N) of
!    a pair of adjacent boundary nodes defining the arc to be removed.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    LNEW, the triangulation data structure created by TRMESH.  On output,
!    updated with the removal of arc IO1-IO2 unless 0 < IER.
!
!    Output, integer ( kind = 4 ) IER, error indicator:
!    0, if no errors were encountered.
!    1, if N, IO1, or IO2 is outside its valid range, or IO1 = IO2.
!    2, if IO1-IO2 is not a boundary arc.
!    3, if the node opposite IO1-IO2 is already a boundary node, and thus IO1
!      or IO2 has only two neighbors or a deletion would result in two
!      triangulations sharing a single node.
!    4, if one of the nodes is a neighbor of the other, but not vice versa,
!      implying an invalid triangulation data structure.
!
!  Local parameters:
!
!    LP =       LIST pointer
!    LPH =      LIST pointer or flag returned by DELNB
!    LPL =      Pointer to the last neighbor of N1, N2, or N3
!    N1,N2,N3 = Nodal indexes of a triangle such that N1->N2
!               is the directed boundary edge associated with IO1-IO2
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  n1 = io1
  n2 = io2
!
!  Test for errors.
!
  if ( n < 4 ) then
    ier = 1
    return
  end if

  if ( n1 < 1 ) then
    ier = 1
    return
  end if

  if ( n < n1 ) then
    ier = 1
    return
  end if

  if ( n2 < 1 ) then
    ier = 1
    return
  end if

  if ( n < n2 ) then
    ier = 1
    return
  end if

  if ( n1 == n2 ) then
    ier = 1
    return
  end if
!
!  Set N1->N2 to the directed boundary edge associated with IO1-IO2:
!  (N1,N2,N3) is a triangle for some N3.
!
  lpl = lend(n2)

  if ( -list(lpl) /= n1 ) then
    n1 = n2
    n2 = io1
    lpl = lend(n2)
    if ( -list(lpl) /= n1 ) then
      ier = 2
      return
    end if
  end if
!
!  Set N3 to the node opposite N1->N2 (the second neighbor
!  of N1), and test for error 3 (N3 already a boundary node).
!
  lpl = lend(n1)
  lp = lptr(lpl)
  lp = lptr(lp)
  n3 = abs ( list(lp) )
  lpl = lend(n3)

  if ( list(lpl) <= 0 ) then
    ier = 3
    return
  end if
!
!  Delete N2 as a neighbor of N1, making N3 the first
!  neighbor, and test for error 4 (N2 not a neighbor
!  of N1).  Note that previously computed pointers may
!  no longer be valid following the call to DELNB.
!
  call delnb ( n1, n2, n, list, lptr, lend, lnew, lph )

  if ( lph < 0 ) then
    ier = 4
    return
  end if
!
!  Delete N1 as a neighbor of N2, making N3 the new last neighbor.
!
  call delnb ( n2, n1, n, list, lptr, lend, lnew, lph )
!
!  Make N3 a boundary node with first neighbor N2 and last neighbor N1.
!
  lp = lstptr ( lend(n3), n1, list, lptr )
  lend(n3) = lp
  list(lp) = -n1
!
!  No errors encountered.
!
  ier = 0

  return
end
subroutine delnb ( n0, nb, n, list, lptr, lend, lnew, lph )

!*****************************************************************************80
!
!! DELNB deletes a neighbor from the adjacency list.
!
!  Discussion:
!
!    This subroutine deletes a neighbor NB from the adjacency
!    list of node N0 (but N0 is not deleted from the adjacency
!    list of NB) and, if NB is a boundary node, makes N0 a
!    boundary node.
!
!    For pointer (LIST index) LPH to NB as a neighbor of N0, the empty
!    LIST, LPTR location LPH is filled in with the values at LNEW-1,
!    pointer LNEW-1 (in LPTR and possibly in LEND) is changed to LPH,
!    and LNEW is decremented.
!
!    This requires a search of LEND and LPTR entailing an
!    expected operation count of O(N).
!
!    This routine is identical to the similarly named routine in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N0, NB, indexes, in the range 1 to N, of a
!    pair of nodes such that NB is a neighbor of N0.  (N0 need not be a
!    neighbor of NB.)
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), LNEW,
!    the data structure defining the triangulation.  On output, updated with
!    the removal of NB from the adjacency list of N0 unless LPH < 0.
!
!    Input, integer ( kind = 4 ) LPH, list pointer to the hole (NB as a
!    neighbor of N0) filled in by the values at LNEW-1 or error indicator:
!    >  0, if no errors were encountered.
!    = -1, if N0, NB, or N is outside its valid range.
!    = -2, if NB is not a neighbor of N0.
!
!  Local parameters:
!
!    I =   DO-loop index
!    LNW = LNEW-1 (output value of LNEW)
!    LP =  LIST pointer of the last neighbor of NB
!    LPB = Pointer to NB as a neighbor of N0
!    LPL = Pointer to the last neighbor of N0
!    LPP = Pointer to the neighbor of N0 that precedes NB
!    NN =  Local copy of N
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lnw
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpb
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpp
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nn

  nn = n
!
!  Test for error 1.
!
  if ( n0 < 1 ) then
    lph = -1
    return
  end if

  if ( nn < n0 .or. nb < 1  .or. &
      nn < nb .or. nn < 3 ) then
    lph = -1
    return
  end if
!
!  Find pointers to neighbors of N0:
!
!  LPL points to the last neighbor,
!  LPP points to the neighbor NP preceding NB, and
!  LPB points to NB.
!
  lpl = lend(n0)
  lpp = lpl
  lpb = lptr(lpp)

  do

    if ( list(lpb) == nb ) then
      go to 2
    end if

    lpp = lpb
    lpb = lptr(lpp)

    if ( lpb == lpl ) then
      exit
    end if

  end do
!
!  Test for error 2 (NB not found).
!
  if ( abs ( list(lpb) ) /= nb ) then
    lph = -2
    return
  end if
!
!  NB is the last neighbor of N0.  Make NP the new last
!  neighbor and, if NB is a boundary node, then make N0
!  a boundary node.
!
  lend(n0) = lpp
  lp = lend(nb)

  if ( list(lp) < 0 ) then
    list(lpp) = -list(lpp)
  end if

  go to 3
!
!  NB is not the last neighbor of N0.  If NB is a boundary
!  node and N0 is not, then make N0 a boundary node with
!  last neighbor NP.
!
2 continue

  lp = lend(nb)

  if ( list(lp) < 0 .and. 0 < list(lpl) ) then
    lend(n0) = lpp
    list(lpp) = -list(lpp)
  end if
!
!  Update LPTR so that the neighbor following NB now follows
!  NP, and fill in the hole at location LPB.
!
3 continue

  lptr(lpp) = lptr(lpb)
  lnw = lnew-1
  list(lpb) = list(lnw)
  lptr(lpb) = lptr(lnw)

  do i = nn, 1, -1
    if ( lend(i) == lnw ) then
      lend(i) = lpb
      exit
    end if
  end do

  do i = 1, lnw-1
    if ( lptr(i) == lnw ) then
      lptr(i) = lpb
    end if
  end do
!
!  No errors encountered.
!
  lnew = lnw
  lph = lpb

  return
end
subroutine delnod ( k, n, x, y, z, list, lptr, lend, lnew, lwk, iwk, ier )

!*****************************************************************************80
!
!! DELNOD deletes a node from a triangulation.
!
!  Discussion:
!
!    This subroutine deletes node K (along with all arcs incident on node K)
!    from a triangulation of N nodes on the unit sphere, and inserts arcs as
!    necessary to produce a triangulation of the remaining N-1 nodes.  If a
!    Delaunay triangulation is input, a Delaunay triangulation will result,
!    and thus, DELNOD reverses the effect of a call to ADDNOD.
!
!    Note that the deletion may result in all remaining nodes
!    being collinear.  This situation is not flagged.
!
!  Modified:
!
!    17 June 2002
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, index (for X, Y, and Z) of the node to be
!    deleted.  1 <= K <= N.
!
!    Input/output, integer ( kind = 4 ) N, the number of nodes in the
!    triangulation.  4 <= N.  Note that N will be decremented following the
!    deletion.
!
!    Input/output, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of
!    the nodes in the triangulation.  On output, updated with elements
!    K+1,...,N+1 shifted up one position, thus overwriting element K,
!    unless 1 <= IER <= 4.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    LNEW, the data structure defining the triangulation, created by TRMESH.
!    On output, updated to reflect the deletion unless 1 <= IER <= 4.
!    Note that the data structure may have been altered if 3 < IER.
!
!    Input/output, integer ( kind = 4 ) LWK, the number of columns reserved for
!    IWK.  LWK must be at least NNB-3, where NNB is the number of neighbors of
!    node K, including an extra pseudo-node if K is a boundary node.
!    On output, the number of IWK columns required unless IER = 1 or IER = 3.
!
!    Output, integer ( kind = 4 ) IWK(2,LWK), indexes of the endpoints of the
!    new arcs added unless LWK = 0 or 1 <= IER <= 4.  (Arcs are associated with
!    columns.)
!
!    Output, integer ( kind = 4 ) IER, error indicator:
!    0, if no errors were encountered.
!    1, if K or N is outside its valid range or LWK < 0 on input.
!    2, if more space is required in IWK.  Refer to LWK.
!    3, if the triangulation data structure is invalid on input.
!    4, if K indexes an interior node with four or more neighbors, none of
!      which can be swapped out due to collinearity, and K cannot therefore
!      be deleted.
!    5, if an error flag (other than IER = 1) was returned by OPTIM.  An error
!      message is written to the standard output unit in this case.
!    6, if error flag 1 was returned by OPTIM.  This is not necessarily an
!      error, but the arcs may not be optimal.
!
!  Local parameters:
!
!    BDRY =    Logical variable with value TRUE iff N1 is a boundary node
!    I,J =     DO-loop indexes
!    IERR =    Error flag returned by OPTIM
!    IWL =     Number of IWK columns containing arcs
!    LNW =     Local copy of LNEW
!    LP =      LIST pointer
!    LP21 =    LIST pointer returned by SWAP
!    LPF,LPL = Pointers to the first and last neighbors of N1
!    LPH =     Pointer (or flag) returned by DELNB
!    LPL2 =    Pointer to the last neighbor of N2
!    LPN =     Pointer to a neighbor of N1
!    LWKL =    Input value of LWK
!    N1 =      Local copy of K
!    N2 =      Neighbor of N1
!    NFRST =   First neighbor of N1:  LIST(LPF)
!    NIT =     Number of iterations in OPTIM
!    NR,NL =   Neighbors of N1 preceding (to the right of) and
!              following (to the left of) N2, respectively
!    NN =      Number of nodes in the triangulation
!    NNB =     Number of neighbors of N1 (including a pseudo-
!              node representing the boundary if N1 is a
!              boundary node)
!    X1,Y1,Z1 = Coordinates of N1
!    X2,Y2,Z2 = Coordinates of N2
!    XL,YL,ZL = Coordinates of NL
!    XR,YR,ZR = Coordinates of NR
!
  implicit none

  integer ( kind = 4 ) n

  logical bdry
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iwk(2,*)
  integer ( kind = 4 ) iwl
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical left
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lnw
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpl2
  integer ( kind = 4 ) lpn
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) lwkl
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nbcnt
  integer ( kind = 4 ) nfrst
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nnb
  integer ( kind = 4 ) nr
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) yl
  real ( kind = 8 ) yr
  real ( kind = 8 ) z(*)
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2
  real ( kind = 8 ) zl
  real ( kind = 8 ) zr
!
!  Set N1 to K and NNB to the number of neighbors of N1 (plus
!  one if N1 is a boundary node), and test for errors.  LPF
!  and LPL are LIST indexes of the first and last neighbors
!  of N1, IWL is the number of IWK columns containing arcs,
!  and BDRY is TRUE iff N1 is a boundary node.
!
  n1 = k
  nn = n

  if ( n1 < 1 ) then
    ier = 1
    return
  end if

  if ( nn < n1 ) then
    ier = 1
    return
  end if

  if ( nn < 4 ) then
    ier = 1
    return
  end if

  if ( lwk < 0 ) then
    ier = 1
    return
  end if

  lpl = lend(n1)
  lpf = lptr(lpl)
  nnb = nbcnt(lpl,lptr)
  bdry = list(lpl) < 0

  if ( bdry ) then
    nnb = nnb + 1
  end if

  if ( nnb < 3 ) then
    ier = 3
    return
  end if

  lwkl = lwk
  lwk = nnb - 3

  if ( lwkl < lwk ) then
    ier = 2
    return
  end if

  iwl = 0

  if ( nnb == 3 ) then
    go to 3
  end if
!
!  Initialize for loop on arcs N1-N2 for neighbors N2 of N1,
!  beginning with the second neighbor.  NR and NL are the
!  neighbors preceding and following N2, respectively, and
!  LP indexes NL.  The loop is exited when all possible
!  swaps have been applied to arcs incident on N1.
!
  x1 = x(n1)
  y1 = y(n1)
  z1 = z(n1)
  nfrst = list(lpf)
  nr = nfrst
  xr = x(nr)
  yr = y(nr)
  zr = z(nr)
  lp = lptr(lpf)
  n2 = list(lp)
  x2 = x(n2)
  y2 = y(n2)
  z2 = z(n2)
  lp = lptr(lp)
!
!  Top of loop:  set NL to the neighbor following N2.
!
  do

    nl = abs ( list(lp) )

    if ( nl == nfrst .and. bdry ) then
      exit
    end if

    xl = x(nl)
    yl = y(nl)
    zl = z(nl)
!
!  Test for a convex quadrilateral.  To avoid an incorrect
!  test caused by collinearity, use the fact that if N1
!  is a boundary node, then N1 LEFT NR->NL and if N2 is
!  a boundary node, then N2 LEFT NL->NR.
!
    lpl2 = lend(n2)
!
!  Nonconvex quadrilateral -- no swap is possible.
!
    if ( .not. ((bdry .or. left(xr,yr,zr,xl,yl,zl,x1,y1, &
        z1)) .and. (list(lpl2) < 0  .or. &
        left(xl,yl,zl,xr,yr,zr,x2,y2,z2))) ) then
      nr = n2
      xr = x2
      yr = y2
      zr = z2
      go to 2
    end if
!
!  The quadrilateral defined by adjacent triangles
!  (N1,N2,NL) and (N2,N1,NR) is convex.  Swap in
!  NL-NR and store it in IWK unless NL and NR are
!  already adjacent, in which case the swap is not
!  possible.  Indexes larger than N1 must be decremented
!  since N1 will be deleted from X, Y, and Z.
!
    call swap ( nl, nr, n1, n2, list, lptr, lend, lp21 )

    if ( lp21 == 0 ) then
      nr = n2
      xr = x2
      yr = y2
      zr = z2
      go to 2
    end if

    iwl = iwl + 1

    if ( nl <= n1 ) then
      iwk(1,iwl) = nl
    else
      iwk(1,iwl) = nl - 1
    end if

    if ( nr <= n1 ) then
      iwk(2,iwl) = nr
    else
      iwk(2,iwl) = nr - 1
    end if
!
!  Recompute the LIST indexes and NFRST, and decrement NNB.
!
    lpl = lend(n1)
    nnb = nnb - 1

    if ( nnb == 3 ) then
      exit
    end if

    lpf = lptr(lpl)
    nfrst = list(lpf)
    lp = lstptr ( lpl, nl, list, lptr )
!
!  NR is not the first neighbor of N1.
!  Back up and test N1-NR for a swap again:  Set N2 to
!  NR and NR to the previous neighbor of N1 -- the
!  neighbor of NR which follows N1.  LP21 points to NL
!  as a neighbor of NR.
!
    if ( nr /= nfrst ) then

      n2 = nr
      x2 = xr
      y2 = yr
      z2 = zr
      lp21 = lptr(lp21)
      lp21 = lptr(lp21)
      nr = abs ( list(lp21) )
      xr = x(nr)
      yr = y(nr)
      zr = z(nr)
      cycle

    end if
!
!  Bottom of loop -- test for termination of loop.
!
2   continue

    if ( n2 == nfrst ) then
      exit
    end if

    n2 = nl
    x2 = xl
    y2 = yl
    z2 = zl
    lp = lptr(lp)

  end do
!
!  Delete N1 and all its incident arcs.  If N1 is an interior
!  node and either 3 < NNB or NNB = 3 and N2 LEFT NR->NL,
!  then N1 must be separated from its neighbors by a plane
!  containing the origin -- its removal reverses the effect
!  of a call to COVSPH, and all its neighbors become
!  boundary nodes.  This is achieved by treating it as if
!  it were a boundary node (setting BDRY to TRUE, changing
!  a sign in LIST, and incrementing NNB).
!
3   continue

  if ( .not. bdry ) then

    if ( 3 < nnb ) then
      bdry = .true.
    else
      lpf = lptr(lpl)
      nr = list(lpf)
      lp = lptr(lpf)
      n2 = list(lp)
      nl = list(lpl)
      bdry = left ( x(nr), y(nr), z(nr), x(nl), y(nl), z(nl), &
        x(n2), y(n2), z(n2) )
    end if
!
!  If a boundary node already exists, then N1 and its
!  neighbors cannot be converted to boundary nodes.
!  (They must be collinear.)  This is a problem if 3 < NNB.
!
    if ( bdry ) then

      do i = 1, nn
        if ( list(lend(i)) < 0 ) then
          bdry = .false.
          go to 5
        end if
      end do

      list(lpl) = -list(lpl)
      nnb = nnb + 1

    end if

  end if

5 continue

  if ( .not. bdry .and. 3 < nnb ) then
    ier = 4
    return
  end if
!
!  Initialize for loop on neighbors.  LPL points to the last
!  neighbor of N1.  LNEW is stored in local variable LNW.
!
  lp = lpl
  lnw = lnew
!
!  Loop on neighbors N2 of N1, beginning with the first.
!
6   continue

    lp = lptr(lp)
    n2 = abs ( list(lp) )

    call delnb ( n2, n1, n, list, lptr, lend, lnw, lph )

    if ( lph < 0 ) then
      ier = 3
      return
    end if
!
!  LP and LPL may require alteration.
!
    if ( lpl == lnw ) then
      lpl = lph
    end if

    if ( lp == lnw ) then
      lp = lph
    end if

    if ( lp /= lpl ) then
      go to 6
    end if
!
!  Delete N1 from X, Y, Z, and LEND, and remove its adjacency
!  list from LIST and LPTR.  LIST entries (nodal indexes)
!  which are larger than N1 must be decremented.
!
  nn = nn - 1

  if ( nn < n1 ) then
    go to 9
  end if

  do i = n1, nn
    x(i) = x(i+1)
    y(i) = y(i+1)
    z(i) = z(i+1)
    lend(i) = lend(i+1)
  end do

  do i = 1, lnw-1

    if ( n1 < list(i) ) then
      list(i) = list(i) - 1
    end if

    if ( list(i) < -n1 ) then
      list(i) = list(i) + 1
    end if

  end do
!
!  For LPN = first to last neighbors of N1, delete the
!  preceding neighbor (indexed by LP).
!
!  Each empty LIST,LPTR location LP is filled in with the
!  values at LNW-1, and LNW is decremented.  All pointers
!  (including those in LPTR and LEND) with value LNW-1
!  must be changed to LP.
!
!  LPL points to the last neighbor of N1.
!
9 continue

  if ( bdry ) then
    nnb = nnb - 1
  end if

  lpn = lpl

  do j = 1, nnb

    lnw = lnw - 1
    lp = lpn
    lpn = lptr(lp)
    list(lp) = list(lnw)
    lptr(lp) = lptr(lnw)

    if ( lptr(lpn) == lnw ) then
      lptr(lpn) = lp
    end if

    if ( lpn == lnw ) then
      lpn = lp
    end if

    do i = nn, 1, -1
      if ( lend(i) == lnw ) then
        lend(i) = lp
        exit
      end if
    end do

    do i = lnw-1, 1, -1
      if ( lptr(i) == lnw ) then
        lptr(i) = lp
      end if
    end do

  end do
!
!  Update N and LNEW, and optimize the patch of triangles
!  containing K (on input) by applying swaps to the arcs in IWK.
!
  n = nn
  lnew = lnw

  if ( 0 < iwl ) then

    nit = 4 * iwl

    call optim ( x, y, z, iwl, list, lptr, lend, nit, iwk, ierr )

    if ( ierr /= 0 .and. ierr /= 1 ) then
      ier = 5
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DELNOD - Fatal error.'
      write ( *, '(a)' ) '  OPTIM failed.'
      write ( *, '(a,i8)' ) '  NIT = ', nit
      write ( *, '(a,i8)' ) '  IERR = ', ierr
      return
    end if

    if ( ierr == 1 ) then
      ier = 6
      return
    end if

  end if

  ier = 0

  return
end

subroutine edge ( in1, in2, x, y, z, lwk, iwk, list, lptr, lend, ier )

!*****************************************************************************80
!
!! EDGE swaps arcs to force two nodes to be adjacent.
!
!  Discussion:
!
!    Given a triangulation of N nodes and a pair of nodal
!    indexes IN1 and IN2, this routine swaps arcs as necessary
!    to force IN1 and IN2 to be adjacent.  Only arcs which
!    intersect IN1-IN2 are swapped out.  If a Delaunay triangu-
!    lation is input, the resulting triangulation is as close
!    as possible to a Delaunay triangulation in the sense that
!    all arcs other than IN1-IN2 are locally optimal.
!
!    A sequence of calls to EDGE may be used to force the
!    presence of a set of edges defining the boundary of a
!    non-convex and/or multiply connected region, or to introduce
!    barriers into the triangulation.  Note that
!    GETNP will not necessarily return closest nodes if the
!    triangulation has been constrained by a call to EDGE.
!    However, this is appropriate in some applications, such
!    as triangle-based interpolation on a nonconvex domain.
!
!  Modified:
!
!    17 June 2002
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IN1, IN2, indexes (of X, Y, and Z) in the
!    range 1 to N defining a pair of nodes to be connected by an arc.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Input/output, integer ( kind = 4 ) LWK.  On input, the number of columns
!    reserved for IWK.  This must be at least NI, the number of arcs that
!    intersect IN1-IN2.  (NI is bounded by N-3.)   On output, the number of
!    arcs which intersect IN1-IN2 (but not more than the input value of LWK)
!    unless IER = 1 or IER = 3.  LWK = 0 if and only if IN1 and IN2 were
!    adjacent (or LWK=0) on input.
!
!    Output, integer ( kind = 4 ) IWK(2*LWK), the indexes of the endpoints of
!    the new arcs other than IN1-IN2 unless 0 < IER or LWK = 0.  New arcs to
!    the left of IN1->IN2 are stored in the first K-1 columns (left portion
!    of IWK), column K contains zeros, and new arcs to the right of IN1->IN2
!    occupy columns K+1,...,LWK.  (K can be determined by searching IWK
!    for the zeros.)
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    the data structure defining the triangulation, created by TRMESH.  On
!    output, updated if necessary to reflect the presence of an arc connecting
!    IN1 and IN2 unless 0 < IER.  The data structure has been altered if
!    4 <= IER.
!
!    Output, integer ( kind = 4 ) IER, error indicator:
!    0, if no errors were encountered.
!    1, if IN1 < 1, IN2 < 1, IN1 = IN2, or LWK < 0 on input.
!    2, if more space is required in IWK.  Refer to LWK.
!    3, if IN1 and IN2 could not be connected due to either an invalid
!      data structure or collinear nodes (and floating point error).
!    4, if an error flag other than IER = 1 was returned by OPTIM.
!    5, if error flag 1 was returned by OPTIM.  This is not necessarily
!      an error, but the arcs other than IN1-IN2 may not be optimal.
!
!  Local parameters:
!
!    DPij =     Dot product <Ni,Nj>
!    I =        DO-loop index and column index for IWK
!    IERR =     Error flag returned by Subroutine OPTIM
!    IWC =      IWK index between IWF and IWL -- NL->NR is
!               stored in IWK(1,IWC)->IWK(2,IWC)
!    IWCP1 =    IWC + 1
!    IWEND =    Input or output value of LWK
!    IWF =      IWK (column) index of the first (leftmost) arc
!               which intersects IN1->IN2
!    IWL =      IWK (column) index of the last (rightmost) are
!               which intersects IN1->IN2
!    LFT =      Flag used to determine if a swap results in the
!               new arc intersecting IN1-IN2 -- LFT = 0 iff
!               N0 = IN1, LFT = -1 implies N0 LEFT IN1->IN2,
!               and LFT = 1 implies N0 LEFT IN2->IN1
!    LP =       List pointer (index for LIST and LPTR)
!    LP21 =     Unused parameter returned by SWAP
!    LPL =      Pointer to the last neighbor of IN1 or NL
!    N0 =       Neighbor of N1 or node opposite NR->NL
!    N1,N2 =    Local copies of IN1 and IN2
!    N1FRST =   First neighbor of IN1
!    N1LST =    (Signed) last neighbor of IN1
!    NEXT =     Node opposite NL->NR
!    NIT =      Flag or number of iterations employed by OPTIM
!    NL,NR =    Endpoints of an arc which intersects IN1-IN2
!               with NL LEFT IN1->IN2
!    X0,Y0,Z0 = Coordinates of N0
!    X1,Y1,Z1 = Coordinates of IN1
!    X2,Y2,Z2 = Coordinates of IN2
!
  implicit none

  real ( kind = 8 ) dp12
  real ( kind = 8 ) dp1l
  real ( kind = 8 ) dp1r
  real ( kind = 8 ) dp2l
  real ( kind = 8 ) dp2r
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) iwc
  integer ( kind = 4 ) iwcp1
  integer ( kind = 4 ) iwend
  integer ( kind = 4 ) iwf
  integer ( kind = 4 ) iwk(2,*)
  integer ( kind = 4 ) iwl
  logical left
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) lft
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1frst
  integer ( kind = 4 ) n1lst
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nr
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) z(*)
  real ( kind = 8 ) z0
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2
!
!  Store IN1, IN2, and LWK in local variables and test for errors.
!
  n1 = in1
  n2 = in2
  iwend = lwk

  if ( n1 < 1 .or. n2 < 1 .or. n1 == n2  .or. iwend < 0 ) then
    ier = 1
    return
  end if
!
!  Test for N2 as a neighbor of N1.  LPL points to the last neighbor of N1.
!
  lpl = lend(n1)
  n0 = abs ( list(lpl) )
  lp = lpl

  do

    if ( n0 == n2 ) then
      ier = 0
      return
    end if

    lp = lptr(lp)
    n0 = list(lp)

    if ( lp == lpl ) then
      exit
    end if

  end do
!
!  Initialize parameters.
!
  iwl = 0
  nit = 0
!
!  Store the coordinates of N1 and N2.
!
  do

    x1 = x(n1)
    y1 = y(n1)
    z1 = z(n1)

    x2 = x(n2)
    y2 = y(n2)
    z2 = z(n2)
!
!  Set NR and NL to adjacent neighbors of N1 such that
!  NR LEFT N2->N1 and NL LEFT N1->N2,
!  (NR Forward N1->N2 or NL Forward N1->N2), and
!  (NR Forward N2->N1 or NL Forward N2->N1).
!
!  Initialization:  Set N1FRST and N1LST to the first and
!  (signed) last neighbors of N1, respectively, and
!  initialize NL to N1FRST.
!
    lpl = lend(n1)
    n1lst = list(lpl)
    lp = lptr(lpl)
    n1frst = list(lp)
    nl = n1frst

    if ( n1lst < 0 ) then
      go to 4
    end if
!
!  N1 is an interior node.  Set NL to the first candidate
!  for NR (NL LEFT N2->N1).
!
    do

      if ( left ( x2, y2, z2, x1, y1, z1, x(nl), y(nl), z(nl) ) ) then
        go to 4
      end if

      lp = lptr(lp)
      nl = list(lp)

      if ( nl == n1frst ) then
        exit
      end if

    end do
!
!  All neighbors of N1 are strictly left of N1->N2.
!
    go to 5
!
!  NL = LIST(LP) LEFT N2->N1.  Set NR to NL and NL to the
!  following neighbor of N1.
!
4   continue

    do

      nr = nl
      lp = lptr(lp)
      nl = abs ( list(lp) )
!
!  NL LEFT N1->N2 and NR LEFT N2->N1.  The Forward tests
!  are employed to avoid an error associated with
!  collinear nodes.
!
      if ( left ( x1, y1, z1, x2, y2, z2, x(nl), y(nl), z(nl) ) ) then

        dp12 = x1 * x2 + y1 * y2 + z1 * z2
        dp1l = x1 * x(nl) + y1 * y(nl) + z1 * z(nl)
        dp2l = x2 * x(nl) + y2 * y(nl) + z2 * z(nl)
        dp1r = x1 * x(nr) + y1 * y(nr) + z1 * z(nr)
        dp2r = x2 * x(nr) + y2 * y(nr) + z2 * z(nr)

        if ( ( 0.0D+00 <= dp2l - dp12 * dp1l .or. &
               0.0D+00 <= dp2r - dp12 * dp1r )  .and. &
             ( 0.0D+00 <= dp1l - dp12 * dp2l .or.  &
               0.0D+00 <= dp1r - dp12 * dp2r ) ) then
          go to 6
        end if
!
!  NL-NR does not intersect N1-N2.  However, there is
!  another candidate for the first arc if NL lies on
!  the line N1-N2.
!
        if ( .not. left ( x2, y2, z2, x1, y1, z1, x(nl), y(nl), z(nl) ) ) then
          exit
        end if

      end if
!
!  Bottom of loop.
!
      if ( nl == n1frst ) then
        exit
      end if

    end do
!
!  Either the triangulation is invalid or N1-N2 lies on the
!  convex hull boundary and an edge NR->NL (opposite N1 and
!  intersecting N1-N2) was not found due to floating point
!  error.  Try interchanging N1 and N2 -- NIT > 0 iff this
!  has already been done.
!
5   continue

    if ( 0 < nit ) then
      ier = 3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EDGE - Fatal error!'
      write ( *, '(a)' ) '  Invalid triangulation, or'
      write ( *, '(a)' ) '  null triangles on boundary.'
      write ( *, '(a,i8)' ) '  IN1 = ', in1
      write ( *, '(a,i8)' ) '  IN2 = ', in2
      return
    end if

    nit = 1

    n3 = n1
    n1 = n2
    n2 = n3

  end do
!
!  Store the ordered sequence of intersecting edges NL->NR in
!  IWK(1,IWL)->IWK(2,IWL).
!
6 continue

  iwl = iwl + 1

  if ( iwend < iwl ) then
    ier = 2
    return
  end if

  iwk(1,iwl) = nl
  iwk(2,iwl) = nr
!
!  Set NEXT to the neighbor of NL which follows NR.
!
  lpl = lend(nl)
  lp = lptr(lpl)
!
!  Find NR as a neighbor of NL.  The search begins with the first neighbor.
!
  do

    if ( list(lp) == nr ) then
      go to 8
    end if

    lp = lptr(lp)

    if ( lp == lpl ) then
      exit
    end if

  end do
!
!  NR must be the last neighbor, and NL->NR cannot be a boundary edge.
!
  if ( list(lp) /= nr ) then
    ier = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EDGE - Fatal error!'
    write ( *, '(a)' ) '  Invalid triangulation, or null triangles on boundary.'
    write ( *, '(a,i8)' ) '  IN1 = ', in1
    write ( *, '(a,i8)' ) '  IN2 = ', in2
    return
  end if
!
!  Set NEXT to the neighbor following NR, and test for
!  termination of the store loop.
!
8 continue

  lp = lptr(lp)
  next = abs ( list(lp) )
!
!  Set NL or NR to NEXT.
!
  if ( next /= n2 ) then

    if ( left ( x1, y1, z1, x2, y2, z2, x(next), y(next), z(next) ) ) then
      nl = next
    else
      nr = next
    end if

    go to 6

  end if
!
!  IWL is the number of arcs which intersect N1-N2.
!  Store LWK.
!

  lwk = iwl
  iwend = iwl
!
!  Initialize for edge swapping loop -- all possible swaps
!  are applied (even if the new arc again intersects
!  N1-N2), arcs to the left of N1->N2 are stored in the
!  left portion of IWK, and arcs to the right are stored in
!  the right portion.  IWF and IWL index the first and last
!  intersecting arcs.
!
  iwf = 1
!
!  Top of loop -- set N0 to N1 and NL->NR to the first edge.
!  IWC points to the arc currently being processed.  LFT
!  <= 0 iff N0 LEFT N1->N2.
!
10 continue

  lft = 0
  n0 = n1
  x0 = x1
  y0 = y1
  z0 = z1
  nl = iwk(1,iwf)
  nr = iwk(2,iwf)
  iwc = iwf
!
!  Set NEXT to the node opposite NL->NR unless IWC is the last arc.
!
11 continue

  if (iwc == iwl) then
    go to 21
  end if

  iwcp1 = iwc + 1
  next = iwk(1,iwcp1)

  if ( next /= nl ) then
    go to 16
  end if

  next = iwk(2,iwcp1)
!
!  NEXT RIGHT N1->N2 and IWC < IWL.  Test for a possible swap.
!
  if ( .not. left ( x0, y0, z0, x(nr), y(nr), z(nr), x(next), &
                  y(next), z(next) ) ) then
    go to 14
  end if

  if ( 0 <= lft ) then
    go to 12
  end if

  if ( .not. left ( x(nl), y(nl), z(nl), x0, y0, z0, x(next), &
                  y(next), z(next) ) ) then
    go to 14
  end if
!
!  Replace NL->NR with N0->NEXT.
!
  call swap ( next, n0, nl, nr, list, lptr, lend, lp21 )
  iwk(1,iwc) = n0
  iwk(2,iwc) = next
  go to 15
!
!  Swap NL-NR for N0-NEXT, shift columns IWC+1,...,IWL to
!  the left, and store N0-NEXT in the right portion of IWK.
!
12 continue

  call swap ( next, n0, nl, nr, list, lptr, lend, lp21 )

  do i = iwcp1, iwl
    iwk(1,i-1) = iwk(1,i)
    iwk(2,i-1) = iwk(2,i)
  end do

  iwk(1,iwl) = n0
  iwk(2,iwl) = next
  iwl = iwl - 1
  nr = next
  go to 11
!
!  A swap is not possible.  Set N0 to NR.
!
14 continue

  n0 = nr
  x0 = x(n0)
  y0 = y(n0)
  z0 = z(n0)
  lft = 1
!
!  Advance to the next arc.
!
15 continue

  nr = next
  iwc = iwc + 1
  go to 11
!
!  NEXT LEFT N1->N2, NEXT .NE. N2, and IWC < IWL.
!  Test for a possible swap.
!
16 continue

  if ( .not. &
    left ( x(nl), y(nl), z(nl), x0, y0, z0, x(next), y(next), z(next) ) ) then
    go to 19
  end if

  if ( lft <= 0 ) then
    go to 17
  end if

  if ( .not. &
    left ( x0, y0, z0, x(nr), y(nr), z(nr), x(next), y(next), z(next) ) ) then
    go to 19
  end if
!
!  Replace NL->NR with NEXT->N0.
!
  call swap ( next, n0, nl, nr, list, lptr, lend, lp21 )
  iwk(1,iwc) = next
  iwk(2,iwc) = n0
  go to 20
!
!  Swap NL-NR for N0-NEXT, shift columns IWF,...,IWC-1 to
!  the right, and store N0-NEXT in the left portion of IWK.
!
17 continue

  call swap ( next, n0, nl, nr, list, lptr, lend, lp21 )

  do i = iwc-1, iwf, -1
    iwk(1,i+1) = iwk(1,i)
    iwk(2,i+1) = iwk(2,i)
  end do

  iwk(1,iwf) = n0
  iwk(2,iwf) = next
  iwf = iwf + 1
  go to 20
!
!  A swap is not possible.  Set N0 to NL.
!
19 continue

  n0 = nl
  x0 = x(n0)
  y0 = y(n0)
  z0 = z(n0)
  lft = -1
!
!  Advance to the next arc.
!
20 continue

  nl = next
  iwc = iwc + 1
  go to 11
!
!  N2 is opposite NL->NR (IWC = IWL).
!
21 continue

  if ( n0 == n1 ) then
    go to 24
  end if

  if ( lft < 0 ) then
    go to 22
  end if
!
!  N0 RIGHT N1->N2.  Test for a possible swap.
!
  if ( .not. left ( x0, y0, z0, x(nr), y(nr), z(nr), x2, y2, z2 ) ) then
    go to 10
  end if
!
!  Swap NL-NR for N0-N2 and store N0-N2 in the right portion of IWK.
!
  call swap ( n2, n0, nl, nr, list, lptr, lend, lp21 )
  iwk(1,iwl) = n0
  iwk(2,iwl) = n2
  iwl = iwl - 1
  go to 10
!
!  N0 LEFT N1->N2.  Test for a possible swap.
!
22 continue

  if ( .not. left ( x(nl), y(nl), z(nl), x0, y0, z0, x2, y2, z2 ) ) then
    go to 10
  end if
!
!  Swap NL-NR for N0-N2, shift columns IWF,...,IWL-1 to the
!  right, and store N0-N2 in the left portion of IWK.
!
  call swap ( n2, n0, nl, nr, list, lptr, lend, lp21 )
  i = iwl

  do

    iwk(1,i) = iwk(1,i-1)
    iwk(2,i) = iwk(2,i-1)
    i = i - 1

    if ( i <= iwf ) then
      exit
    end if

  end do

  iwk(1,iwf) = n0
  iwk(2,iwf) = n2
  iwf = iwf + 1
  go to 10
!
!  IWF = IWC = IWL.  Swap out the last arc for N1-N2 and store zeros in IWK.
!
24 continue

  call swap ( n2, n1, nl, nr, list, lptr, lend, lp21 )
  iwk(1,iwc) = 0
  iwk(2,iwc) = 0
!
!  Optimization procedure.
!
!  Optimize the set of new arcs to the left of IN1->IN2.
!
  ier = 0

  if ( 1 < iwc ) then

    nit = 4 * ( iwc - 1 )

    call optim ( x, y, z, iwc-1, list, lptr, lend, nit, iwk, ierr )

    if ( ierr /= 0 .and. ierr /= 1 ) then
      ier = 4
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EDGE - Fatal error!'
      write ( *, '(a)' ) '  An error occurred in OPTIM.'
      write ( *, '(a,i8)' ) '  NIT = ', nit
      write ( *, '(a,i8)' ) '  IER = ', ier
      return
    end if

    if ( ierr == 1 ) then
      ier = 5
    end if

  end if
!
!  Optimize the set of new arcs to the right of IN1->IN2.
!
  if ( iwc < iwend ) then

    nit = 4 * ( iwend - iwc )

    call optim ( x, y, z, iwend-iwc, list, lptr, lend, nit, iwk(1,iwc+1), ierr )

    if ( ierr /= 0 .and. ierr /= 1) then
      ier = 4
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EDGE - Fatal error!'
      write ( *, '(a)' ) '  An error occurred in OPTIM.'
      write ( *, '(a,i8)' ) '  NIT = ', nit
      write ( *, '(a,i8)' ) '  IER = ', ier
      return
    end if

    if ( ierr == 1 ) then
      ier = 5
      return
    end if

  end if

  if ( ier == 5 ) then
    ier = 5
    return
  end if

  return
end
subroutine insert ( k, lp, list, lptr, lnew )

!*****************************************************************************80
!
!! INSERT inserts K as a neighbor of N1.
!
!  Discussion:
!
!    This subroutine inserts K as a neighbor of N1 following
!    N2, where LP is the LIST pointer of N2 as a neighbor of
!    N1.  Note that, if N2 is the last neighbor of N1, K will
!    become the first neighbor (even if N1 is a boundary node).
!
!    This routine is identical to the similarly named routine in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the index of the node to be inserted.
!
!    Input, integer ( kind = 4 ) LP, the LIST pointer of N2 as a neighbor of N1.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LNEW,
!    the data structure defining the triangulation, created by TRMESH.
!    On output, updated with the addition of node K.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lsav

  lsav = lptr(lp)
  lptr(lp) = lnew
  list(lnew) = k
  lptr(lnew) = lsav
  lnew = lnew + 1

  return
end
subroutine intadd ( kk, i1, i2, i3, list, lptr, lend, lnew )

!*****************************************************************************80
!
!! INTADD adds an interior node to a triangulation.
!
!  Discussion:
!
!    This subroutine adds an interior node to a triangulation
!    of a set of points on the unit sphere.  The data structure
!    is updated with the insertion of node KK into the triangle
!    whose vertices are I1, I2, and I3.  No optimization of the
!    triangulation is performed.
!
!    This routine is identical to the similarly named routine in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KK, the index of the node to be inserted.
!    1 <= KK and KK must not be equal to I1, I2, or I3.
!
!    Input, integer ( kind = 4 ) I1, I2, I3, indexes of the
!    counterclockwise-ordered sequence of vertices of a triangle which contains
!    node KK.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), LNEW,
!    the data structure defining the triangulation, created by TRMESH.  Triangle
!    (I1,I2,I3) must be included in the triangulation.
!    On output, updated with the addition of node KK.  KK
!    will be connected to nodes I1, I2, and I3.
!
!  Local parameters:
!
!    K =        Local copy of KK
!    LP =       LIST pointer
!    N1,N2,N3 = Local copies of I1, I2, and I3
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  k = kk
!
!  Initialization.
!
  n1 = i1
  n2 = i2
  n3 = i3
!
!  Add K as a neighbor of I1, I2, and I3.
!
  lp = lstptr ( lend(n1), n2, list, lptr )
  call insert ( k, lp, list, lptr, lnew )

  lp = lstptr ( lend(n2), n3, list, lptr )
  call insert ( k, lp, list, lptr, lnew )

  lp = lstptr ( lend(n3), n1, list, lptr )
  call insert ( k, lp, list, lptr, lnew )
!
!  Add I1, I2, and I3 as neighbors of K.
!
  list(lnew) = n1
  list(lnew+1) = n2
  list(lnew+2) = n3
  lptr(lnew) = lnew + 1
  lptr(lnew+1) = lnew + 2
  lptr(lnew+2) = lnew
  lend(k) = lnew + 2
  lnew = lnew + 3

  return
end
function jrand ( n, ix, iy, iz )

!*****************************************************************************80
!
!! JRAND returns a random integer between 1 and N.
!
!  Discussion:
!
!   This function returns a uniformly distributed pseudorandom integer
!   in the range 1 to N.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Brian Wichmann, David Hill,
!    An Efficient and Portable Pseudo-random Number Generator,
!    Applied Statistics,
!    Volume 31, Number 2, 1982, pages 188-190.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum value to be returned.
!
!    Input/output, integer ( kind = 4 ) IX, IY, IZ = seeds initialized to
!    values in the range 1 to 30,000 before the first call to JRAND, and
!    not altered between subsequent calls (unless a sequence of random
!    numbers is to be repeated by reinitializing the seeds).
!
!    Output, integer ( kind = 4 ) JRAND, a random integer in the range 1 to N.
!
!  Local parameters:
!
!    U = Pseudo-random number uniformly distributed in the interval (0,1).
!    X = Pseudo-random number in the range 0 to 3 whose fractional part is U.
!
  implicit none

  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) jrand
  integer ( kind = 4 ) n
  real ( kind = 8 ) u
  real ( kind = 8 ) x

  ix = mod ( 171 * ix, 30269 )
  iy = mod ( 172 * iy, 30307 )
  iz = mod ( 170 * iz, 30323 )

  x = ( real ( ix, kind = 8 ) / 30269.0D+00 ) &
    + ( real ( iy, kind = 8 ) / 30307.0D+00 ) &
    + ( real ( iz, kind = 8 ) / 30323.0D+00 )

  u = x - int ( x )
  jrand = int( real ( n, kind = 8 ) * u + 1.0D+00 )

  return
end
function left ( x1, y1, z1, x2, y2, z2, x0, y0, z0 )

!*****************************************************************************80
!
!! LEFT determines whether a node is to the left of a plane through the origin.
!
!  Discussion:
!
!    This function determines whether node N0 is in the
!    (closed) left hemisphere defined by the plane containing
!    N1, N2, and the origin, where left is defined relative to
!    an observer at N1 facing N2.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, Z1 = Coordinates of N1.
!
!    Input, real ( kind = 8 ) X2, Y2, Z2 = Coordinates of N2.
!
!    Input, real ( kind = 8 ) X0, Y0, Z0 = Coordinates of N0.
!
!    Output, logical LEFT = TRUE if and only if N0 is in the closed
!    left hemisphere.
!
  implicit none

  logical              left
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) z0
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2
!
!  LEFT = TRUE iff <N0,N1 X N2> = det(N0,N1,N2) >= 0.
!
  left = x0 * ( y1 * z2 - y2 * z1 ) &
       - y0 * ( x1 * z2 - x2 * z1 ) &
       + z0 * ( x1 * y2 - x2 * y1 ) >= 0.0D+00

  return
end
function lstptr ( lpl, nb, list, lptr )

!*****************************************************************************80
!
!! LSTPTR returns the index of NB in the adjacency list.
!
!  Discussion:
!
!    This function returns the index (LIST pointer) of NB in
!    the adjacency list for N0, where LPL = LEND(N0).
!
!    This function is identical to the similarly named function in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LPL, is LEND(N0).
!
!    Input, integer ( kind = 4 ) NB, index of the node whose pointer is to
!    be returned.  NB must be connected to N0.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), the data
!    structure defining the triangulation, created by TRMESH.
!
!    Output, integer ( kind = 4 ) LSTPTR, pointer such that LIST(LSTPTR) = NB or
!    LIST(LSTPTR) = -NB, unless NB is not a neighbor of N0, in which
!    case LSTPTR = LPL.
!
!  Local parameters:
!
!    LP = LIST pointer
!    ND = Nodal index
!
  implicit none

  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nd

  lp = lptr(lpl)

  do

    nd = list(lp)

    if ( nd == nb ) then
      exit
    end if

    lp = lptr(lp)

    if ( lp == lpl ) then
      exit
    end if

  end do

  lstptr = lp

  return
end
function nbcnt ( lpl, lptr )

!*****************************************************************************80
!
!! NBCNT returns the number of neighbors of a node.
!
!  Discussion:
!
!    This function returns the number of neighbors of a node
!    N0 in a triangulation created by TRMESH.
!
!    The number of neighbors also gives the order of the Voronoi
!    polygon containing the point.  Thus, a neighbor count of 6
!    means the node is contained in a 6-sided Voronoi region.
!
!    This function is identical to the similarly named function in TRIPACK.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LPL = LIST pointer to the last neighbor of N0;
!    LPL = LEND(N0).
!
!    Input, integer ( kind = 4 ) LPTR(6*(N-2)), pointers associated with LIST.
!
!    Output, integer ( kind = 4 ) NBCNT, the number of neighbors of N0.
!
!  Local parameters:
!
!    K =  Counter for computing the number of neighbors.
!
!    LP = LIST pointer
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) nbcnt

  lp = lpl
  k = 1

  do

    lp = lptr(lp)

    if ( lp == lpl ) then
      exit
    end if

    k = k + 1

  end do

  nbcnt = k

  return
end

function nearnd ( p, ist, n, x, y, z, list, lptr, lend, al )

!*****************************************************************************80
!
!! NEARND returns the nearest node to a given point.
!
!  Discussion:
!
!    Given a point P on the surface of the unit sphere and a
!    Delaunay triangulation created by TRMESH, this
!    function returns the index of the nearest triangulation
!    node to P.
!
!    The algorithm consists of implicitly adding P to the
!    triangulation, finding the nearest neighbor to P, and
!    implicitly deleting P from the triangulation.  Thus, it
!    is based on the fact that, if P is a node in a Delaunay
!    triangulation, the nearest node to P is a neighbor of P.
!
!    For large values of N, this procedure will be faster than
!    the naive approach of computing the distance from P to every node.
!
!    Note that the number of candidates for NEARND (neighbors of P)
!    is limited to LMAX defined in the PARAMETER statement below.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P(3), the Cartesian coordinates of the point P to
!    be located relative to the triangulation.  It is assumed
!    that P(1)**2 + P(2)**2 + P(3)**2 = 1, that is, that the
!    point lies on the unit sphere.
!
!    Input, integer ( kind = 4 ) IST, the index of the node at which the search
!    is to begin.  The search time depends on the proximity of this
!    node to P.  If no good candidate is known, any value between
!    1 and N will do.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    N must be at least 3.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the Cartesian coordinates of
!    the nodes.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    the data structure defining the triangulation, created by TRMESH.
!
!    Output, real ( kind = 8 ) AL, the arc length between P and node NEARND.
!    Because both points are on the unit sphere, this is also
!    the angular separation in radians.
!
!    Output, integer ( kind = 4 ) NEARND, the index of the nearest node to P.
!    NEARND will be 0 if N < 3 or the triangulation data structure
!    is invalid.
!
!  Local parameters:
!
!    B1,B2,B3 =  Unnormalized barycentric coordinates returned by TRFIND
!    DS1 =       (Negative cosine of the) distance from P to N1
!    DSR =       (Negative cosine of the) distance from P to NR
!    DX1,..DZ3 = Components of vectors used by the swap test
!    I1,I2,I3 =  Nodal indexes of a triangle containing P, or
!                the rightmost (I1) and leftmost (I2) visible
!                boundary nodes as viewed from P
!    L =         Length of LISTP/LPTRP and number of neighbors of P
!    LMAX =      Maximum value of L
!    LISTP =     Indexes of the neighbors of P
!    LPTRP =     Array of pointers in 1-1 correspondence with LISTP elements
!    LP =        LIST pointer to a neighbor of N1 and LISTP pointer
!    LP1,LP2 =   LISTP indexes (pointers)
!    LPL =       Pointer to the last neighbor of N1
!    N1 =        Index of a node visible from P
!    N2 =        Index of an endpoint of an arc opposite P
!    N3 =        Index of the node opposite N1->N2
!    NN =        Local copy of N
!    NR =        Index of a candidate for the nearest node to P
!    NST =       Index of the node at which TRFIND begins the search
!
  implicit none

  integer ( kind = 4 ), parameter :: lmax = 25
  integer ( kind = 4 ) n

  real ( kind = 8 ) al
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) ds1
  real ( kind = 8 ) dsr
  real ( kind = 8 ) dx1
  real ( kind = 8 ) dx2
  real ( kind = 8 ) dx3
  real ( kind = 8 ) dy1
  real ( kind = 8 ) dy2
  real ( kind = 8 ) dy3
  real ( kind = 8 ) dz1
  real ( kind = 8 ) dz2
  real ( kind = 8 ) dz3
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ist
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) listp(lmax)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp1
  integer ( kind = 4 ) lp2
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) lptrp(lmax)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) nearnd
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nst
  real ( kind = 8 ) p(3)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  nearnd = 0
  al = 0.0D+00
!
!  Store local parameters and test for N invalid.
!
  nn = n

  if ( nn < 3 ) then
    return
  end if

  nst = ist

  if ( nst < 1 .or. nn < nst ) then
    nst = 1
  end if
!
!  Find a triangle (I1,I2,I3) containing P, or the rightmost
!  (I1) and leftmost (I2) visible boundary nodes as viewed from P.
!
  call trfind ( nst, p, n, x, y, z, list, lptr, lend, b1, b2, b3, i1, i2, i3 )
!
!  Test for collinear nodes.
!
  if ( i1 == 0 ) then
    return
  end if
!
!  Store the linked list of 'neighbors' of P in LISTP and
!  LPTRP.  I1 is the first neighbor, and 0 is stored as
!  the last neighbor if P is not contained in a triangle.
!  L is the length of LISTP and LPTRP, and is limited to
!  LMAX.
!
  if ( i3 /= 0 ) then

    listp(1) = i1
    lptrp(1) = 2
    listp(2) = i2
    lptrp(2) = 3
    listp(3) = i3
    lptrp(3) = 1
    l = 3

  else

    n1 = i1
    l = 1
    lp1 = 2
    listp(l) = n1
    lptrp(l) = lp1
!
!  Loop on the ordered sequence of visible boundary nodes
!  N1 from I1 to I2.
!
    do

      lpl = lend(n1)
      n1 = -list(lpl)
      l = lp1
      lp1 = l+1
      listp(l) = n1
      lptrp(l) = lp1

      if ( n1 == i2 .or. lmax <= lp1 ) then
        exit
      end if

    end do

    l = lp1
    listp(l) = 0
    lptrp(l) = 1

  end if
!
!  Initialize variables for a loop on arcs N1-N2 opposite P
!  in which new 'neighbors' are 'swapped' in.  N1 follows
!  N2 as a neighbor of P, and LP1 and LP2 are the LISTP
!  indexes of N1 and N2.
!
  lp2 = 1
  n2 = i1
  lp1 = lptrp(1)
  n1 = listp(lp1)
!
!  Begin loop:  find the node N3 opposite N1->N2.
!
  do

    lp = lstptr ( lend(n1), n2, list, lptr )

    if ( 0 <= list(lp) ) then

      lp = lptr(lp)
      n3 = abs ( list(lp) )
!
!  Swap test:  Exit the loop if L = LMAX.
!
      if ( l == lmax ) then
        exit
      end if

      dx1 = x(n1) - p(1)
      dy1 = y(n1) - p(2)
      dz1 = z(n1) - p(3)

      dx2 = x(n2) - p(1)
      dy2 = y(n2) - p(2)
      dz2 = z(n2) - p(3)

      dx3 = x(n3) - p(1)
      dy3 = y(n3) - p(2)
      dz3 = z(n3) - p(3)
!
!  Swap:  Insert N3 following N2 in the adjacency list for P.
!  The two new arcs opposite P must be tested.
!
      if ( dx3 * ( dy2 * dz1 - dy1 * dz2 ) - &
           dy3 * ( dx2 * dz1 - dx1 * dz2 ) + &
           dz3 * ( dx2 * dy1 - dx1 * dy2 ) > 0.0D+00 ) then

        l = l+1
        lptrp(lp2) = l
        listp(l) = n3
        lptrp(l) = lp1
        lp1 = l
        n1 = n3
        cycle

      end if

    end if
!
!  No swap:  Advance to the next arc and test for termination
!  on N1 = I1 (LP1 = 1) or N1 followed by 0.
!
    if ( lp1 == 1 ) then
      exit
    end if

    lp2 = lp1
    n2 = n1
    lp1 = lptrp(lp1)
    n1 = listp(lp1)

    if ( n1 == 0 ) then
      exit
    end if

  end do
!
!  Set NR and DSR to the index of the nearest node to P and
!  an increasing function (negative cosine) of its distance
!  from P, respectively.
!
  nr = i1
  dsr = -( x(nr) * p(1) + y(nr) * p(2) + z(nr) * p(3) )

  do lp = 2, l

    n1 = listp(lp)

    if ( n1 == 0 ) then
      cycle
    end if

    ds1 = -( x(n1) * p(1) + y(n1) * p(2) + z(n1) * p(3) )

    if ( ds1 < dsr ) then
      nr = n1
      dsr = ds1
    end if

  end do

  dsr = -dsr
  dsr = min ( dsr, 1.0D+00 )

  al = acos ( dsr )
  nearnd = nr

  return
end
subroutine optim ( x, y, z, na, list, lptr, lend, nit, iwk, ier )

!*****************************************************************************80
!
!! OPTIM optimizes the quadrilateral portion of a triangulation.
!
!  Discussion:
!
!    Given a set of NA triangulation arcs, this subroutine
!    optimizes the portion of the triangulation consisting of
!    the quadrilaterals (pairs of adjacent triangles) which
!    have the arcs as diagonals by applying the circumcircle
!    test and appropriate swaps to the arcs.
!
!    An iteration consists of applying the swap test and
!    swaps to all NA arcs in the order in which they are
!    stored.  The iteration is repeated until no swap occurs
!    or NIT iterations have been performed.  The bound on the
!    number of iterations may be necessary to prevent an
!    infinite loop caused by cycling (reversing the effect of a
!    previous swap) due to floating point inaccuracy when four
!    or more nodes are nearly cocircular.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(*), Y(*), Z(*), the nodal coordinates.
!
!    Input, integer ( kind = 4 ) NA, the number of arcs in the set.  NA >= 0.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    the data structure defining the triangulation, created by TRMESH.
!    On output, updated to reflect the swaps.
!
!    Input/output, integer ( kind = 4 ) NIT.  On input, the maximum number of
!    iterations to be performed.  NIT = 4*NA should be sufficient.  NIT >= 1.
!    On output, the number of iterations performed.
!
!    Input/output, integer ( kind = 4 ) IWK(2,NA), the nodal indexes of the arc
!    endpoints (pairs of endpoints are stored in columns).  On output, endpoint
!    indexes of the new set of arcs reflecting the swaps.
!
!    Output, integer ( kind = 4 ) IER, error indicator:
!    0, if no errors were encountered.
!    1, if a swap occurred on the last of MAXIT iterations, where MAXIT is the
!      value of NIT on input.  The new set of arcs is not necessarily optimal
!      in this case.
!    2, if NA < 0 or NIT < 1 on input.
!    3, if IWK(2,I) is not a neighbor of IWK(1,I) for some I in the range 1
!      to NA.  A swap may have occurred in this case.
!    4, if a zero pointer was returned by subroutine SWAP.
!
!  Local parameters:
!
!    I =       Column index for IWK
!    IO1,IO2 = Nodal indexes of the endpoints of an arc in IWK
!    ITER =    Iteration count
!    LP =      LIST pointer
!    LP21 =    Parameter returned by SWAP (not used)
!    LPL =     Pointer to the last neighbor of IO1
!    LPP =     Pointer to the node preceding IO2 as a neighbor of IO1
!    MAXIT =   Input value of NIT
!    N1,N2 =   Nodes opposite IO1->IO2 and IO2->IO1, respectively
!    NNA =     Local copy of NA
!    SWP =     Flag set to TRUE iff a swap occurs in the optimization loop
!
  implicit none

  integer ( kind = 4 ) na

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) iwk(2,na)
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) maxit
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nna
  logical              swp
  logical              swptst
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) z(*)

  nna = na
  maxit = nit

  if ( nna < 0 .or. maxit < 1 ) then
    nit = 0
    ier = 2
    return
  end if
!
!  Initialize iteration count ITER and test for NA = 0.
!
  iter = 0

  if ( nna == 0 ) then
    nit = 0
    ier = 0
    return
  end if
!
!  Top of loop.
!  SWP = TRUE iff a swap occurred in the current iteration.
!
  do

    if ( maxit <= iter ) then
      nit = iter
      ier = 1
      return
    end if

    iter = iter + 1
    swp = .false.
!
!  Inner loop on arcs IO1-IO2.
!
    do i = 1, nna

      io1 = iwk(1,i)
      io2 = iwk(2,i)
!
!  Set N1 and N2 to the nodes opposite IO1->IO2 and
!  IO2->IO1, respectively.  Determine the following:
!
!  LPL = pointer to the last neighbor of IO1,
!  LP = pointer to IO2 as a neighbor of IO1, and
!  LPP = pointer to the node N2 preceding IO2.
!
      lpl = lend(io1)
      lpp = lpl
      lp = lptr(lpp)

      do

        if ( list(lp) == io2 ) then
          go to 3
        end if

        lpp = lp
        lp = lptr(lpp)

        if ( lp == lpl ) then
          exit
        end if

      end do
!
!  IO2 should be the last neighbor of IO1.  Test for no
!  arc and bypass the swap test if IO1 is a boundary node.
!
      if ( abs ( list(lp) ) /= io2 ) then
        nit = iter
        ier = 3
        return
      end if

      if ( list(lp) < 0 ) then
        go to 4
      end if
!
!  Store N1 and N2, or bypass the swap test if IO1 is a
!  boundary node and IO2 is its first neighbor.
!
3     continue

      n2 = list(lpp)
!
!  Test IO1-IO2 for a swap, and update IWK if necessary.
!
      if ( 0 <= n2 ) then

        lp = lptr(lp)
        n1 = abs ( list(lp) )

        if ( swptst ( n1, n2, io1, io2, x, y, z ) ) then

          call swap ( n1, n2, io1, io2, list, lptr, lend, lp21 )

          if ( lp21 == 0 ) then
            nit = iter
            ier = 4
            return
          end if

          swp = .true.
          iwk(1,i) = n1
          iwk(2,i) = n2

        end if

      end if

4     continue

    end do

    if ( .not. swp ) then
      exit
    end if

  end do

  nit = iter
  ier = 0

  return
end
function store ( x )

!*****************************************************************************80
!
!! STORE forces its argument to be stored.
!
!  Discussion:
!
!    This function forces its argument X to be stored in a
!    memory location, thus providing a means of determining
!    floating point number characteristics (such as the machine
!    precision) when it is necessary to avoid computation in
!    high precision registers.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value to be stored.
!
!    Output, real ( kind = 8 ) STORE, the value of X after it has been stored
!    and possibly truncated or rounded to the single precision word length.
!
  implicit none

  real ( kind = 8 ) store
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  common /stcom/ y

  y = x
  store = y

  return
end
subroutine swap ( in1, in2, io1, io2, list, lptr, lend, lp21 )

!*****************************************************************************80
!
!! SWAP replaces the diagonal arc of a quadrilateral with the other diagonal.
!
!  Discussion:
!
!    Given a triangulation of a set of points on the unit
!    sphere, this subroutine replaces a diagonal arc in a
!    strictly convex quadrilateral (defined by a pair of adja-
!    cent triangles) with the other diagonal.  Equivalently, a
!    pair of adjacent triangles is replaced by another pair
!    having the same union.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IN1, IN2, IO1, IO2, nodal indexes of the
!    vertices of the quadrilateral.  IO1-IO2 is replaced by IN1-IN2.
!    (IO1,IO2,IN1) and (IO2,IO1,IN2) must be triangles on input.
!
!    Input/output, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    the data structure defining the triangulation, created by TRMESH.
!    On output, updated with the swap; triangles (IO1,IO2,IN1) an (IO2,IO1,IN2)
!    are replaced by (IN1,IN2,IO2) and (IN2,IN1,IO1) unless LP21 = 0.
!
!    Output, integer ( kind = 4 ) LP21, index of IN1 as a neighbor of IN2 after
!    the swap is performed unless IN1 and IN2 are adjacent on input, in which
!    case LP21 = 0.
!
!  Local parameters:
!
!    LP, LPH, LPSAV = LIST pointers
!
  implicit none

  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpsav
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
!
!  Test for IN1 and IN2 adjacent.
!
  lp = lstptr ( lend(in1), in2, list, lptr )

  if ( abs ( list(lp) ) == in2 ) then
    lp21 = 0
    return
  end if
!
!  Delete IO2 as a neighbor of IO1.
!
  lp = lstptr ( lend(io1), in2, list, lptr )
  lph = lptr(lp)
  lptr(lp) = lptr(lph)
!
!  If IO2 is the last neighbor of IO1, make IN2 the last neighbor.
!
  if ( lend(io1) == lph ) then
    lend(io1) = lp
  end if
!
!  Insert IN2 as a neighbor of IN1 following IO1 using the hole created above.
!
  lp = lstptr ( lend(in1), io1, list, lptr )
  lpsav = lptr(lp)
  lptr(lp) = lph
  list(lph) = in2
  lptr(lph) = lpsav
!
!  Delete IO1 as a neighbor of IO2.
!
  lp = lstptr ( lend(io2), in1, list, lptr )
  lph = lptr(lp)
  lptr(lp) = lptr(lph)
!
!  If IO1 is the last neighbor of IO2, make IN1 the last neighbor.
!
  if ( lend(io2) == lph ) then
    lend(io2) = lp
  end if
!
!  Insert IN1 as a neighbor of IN2 following IO2.
!
  lp = lstptr ( lend(in2), io2, list, lptr )
  lpsav = lptr(lp)
  lptr(lp) = lph
  list(lph) = in1
  lptr(lph) = lpsav
  lp21 = lph

  return
end
function swptst ( n1, n2, n3, n4, x, y, z )

!*****************************************************************************80
!
!! SWPTST decides whether to replace a diagonal arc by the other.
!
!  Discussion:
!
!    This function decides whether or not to replace a
!    diagonal arc in a quadrilateral with the other diagonal.
!    The decision will be to swap (SWPTST = TRUE) if and only
!    if N4 lies above the plane (in the half-space not contain-
!    ing the origin) defined by (N1,N2,N3), or equivalently, if
!    the projection of N4 onto this plane is interior to the
!    circumcircle of (N1,N2,N3).  The decision will be for no
!    swap if the quadrilateral is not strictly convex.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, N4, indexes of the four nodes
!    defining the quadrilateral with N1 adjacent to N2, and (N1,N2,N3) in
!    counterclockwise order.  The arc connecting N1 to N2 should be replaced
!    by an arc connecting N3 to N4 if SWPTST = TRUE.  Refer to subroutine SWAP.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Output, logical SWPTST, TRUE if and only if the arc connecting N1
!    and N2 should be swapped for an arc connecting N3 and N4.
!
!  Local parameters:
!
!    DX1,DY1,DZ1 = Coordinates of N4->N1
!    DX2,DY2,DZ2 = Coordinates of N4->N2
!    DX3,DY3,DZ3 = Coordinates of N4->N3
!    X4,Y4,Z4 =    Coordinates of N4
!
  implicit none

  real ( kind = 8 ) dx1
  real ( kind = 8 ) dx2
  real ( kind = 8 ) dx3
  real ( kind = 8 ) dy1
  real ( kind = 8 ) dy2
  real ( kind = 8 ) dy3
  real ( kind = 8 ) dz1
  real ( kind = 8 ) dz2
  real ( kind = 8 ) dz3
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  logical              swptst
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x4
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y4
  real ( kind = 8 ) z(*)
  real ( kind = 8 ) z4

  x4 = x(n4)
  y4 = y(n4)
  z4 = z(n4)
  dx1 = x(n1) - x4
  dx2 = x(n2) - x4
  dx3 = x(n3) - x4
  dy1 = y(n1) - y4
  dy2 = y(n2) - y4
  dy3 = y(n3) - y4
  dz1 = z(n1) - z4
  dz2 = z(n2) - z4
  dz3 = z(n3) - z4
!
!  N4 lies above the plane of (N1,N2,N3) iff N3 lies above
!  the plane of (N2,N1,N4) iff Det(N3-N4,N2-N4,N1-N4) =
!  (N3-N4,N2-N4 X N1-N4) > 0.
!
  swptst =  dx3 * ( dy2 * dz1 - dy1 * dz2 ) &
          - dy3 * ( dx2 * dz1 - dx1 * dz2 ) &
          + dz3 * ( dx2 * dy1 - dx1 * dy2 ) > 0.0D+00

  return
end
subroutine trans ( n, rlat, rlon, x, y, z )

!*****************************************************************************80
!
!! TRANS transforms spherical coordinates to Cartesian coordinates.
!
!  Discussion:
!
!    This subroutine transforms spherical coordinates into
!    Cartesian coordinates on the unit sphere for input to
!    TRMESH.  Storage for X and Y may coincide with
!    storage for RLAT and RLON if the latter need not be saved.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes (points on the unit
!    sphere) whose coordinates are to be transformed.
!
!    Input, real ( kind = 8 ) RLAT(N), latitudes of the nodes in radians.
!
!    Input, real ( kind = 8 ) RLON(N), longitudes of the nodes in radians.
!
!    Output, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates in the
!    range -1 to 1.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1 for I = 1 to N.
!
!  Local parameters:
!
!    COSPHI = cos(PHI)
!    I =      DO-loop index
!    NN =     Local copy of N
!    PHI =    Latitude
!    THETA =  Longitude
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cosphi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nn
  real ( kind = 8 ) phi
  real ( kind = 8 ) rlat(n)
  real ( kind = 8 ) rlon(n)
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  nn = n

  do i = 1, nn
    phi = rlat(i)
    theta = rlon(i)
    cosphi = cos ( phi )
    x(i) = cosphi * cos ( theta )
    y(i) = cosphi * sin ( theta )
    z(i) = sin ( phi )
  end do

  return
end
subroutine trfind ( nst, p, n, x, y, z, list, lptr, lend, b1, b2, b3, i1, &
  i2, i3 )

!*****************************************************************************80
!
!! TRFIND locates a point relative to a triangulation.
!
!  Discussion:
!
!    This subroutine locates a point P relative to a triangulation
!    created by TRMESH.  If P is contained in
!    a triangle, the three vertex indexes and barycentric
!    coordinates are returned.  Otherwise, the indexes of the
!    visible boundary nodes are returned.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, index of a node at which TRFIND begins
!    its search.  Search time depends on the proximity of this node to P.
!
!    Input, real ( kind = 8 ) P(3), the x, y, and z coordinates (in that order)
!    of the point P to be located.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of the
!    triangulation nodes (unit vectors).
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), the
!    data structure defining the triangulation, created by TRMESH.
!
!    Output, real ( kind = 8 ) B1, B2, B3, the unnormalized barycentric
!    coordinates of the central projection of P onto the underlying planar
!    triangle if P is in the convex hull of the nodes.  These parameters
!    are not altered if I1 = 0.
!
!    Output, integer ( kind = 4 ) I1, I2, I3, the counterclockwise-ordered
!    vertex indexes of a triangle containing P if P is contained in a triangle.
!    If P is not in the convex hull of the nodes, I1 and I2 are the rightmost
!    and leftmost (boundary) nodes that are visible from P, and I3 = 0.  (If
!    all boundary nodes are visible from P, then I1 and I2 coincide.)
!    I1 = I2 = I3 = 0 if P and all of the nodes are coplanar (lie on a
!    common great circle.
!
!  Local parameters:
!
!    EPS =      Machine precision
!    IX,IY,IZ = Integer seeds for JRAND
!    LP =       LIST pointer
!    N0,N1,N2 = Nodes in counterclockwise order defining a
!               cone (with vertex N0) containing P, or end-
!               points of a boundary edge such that P Right
!               N1->N2
!    N1S,N2S =  Initially-determined values of N1 and N2
!    N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
!    NEXT =     Candidate for I1 or I2 when P is exterior
!    NF,NL =    First and last neighbors of N0, or first
!               (rightmost) and last (leftmost) nodes
!               visible from P when P is exterior to the
!               triangulation
!    PTN1 =     Scalar product <P,N1>
!    PTN2 =     Scalar product <P,N2>
!    Q =        (N2 X N1) X N2  or  N1 X (N2 X N1) -- used in
!               the boundary traversal when P is exterior
!    S12 =      Scalar product <N1,N2>
!    TOL =      Tolerance (multiple of EPS) defining an upper
!               bound on the magnitude of a negative bary-
!               centric coordinate (B1 or B2) for P in a
!               triangle -- used to avoid an infinite number
!               of restarts with 0 <= B3 < EPS and B1 < 0 or
!               B2 < 0 but small in magnitude
!    XP,YP,ZP = Local variables containing P(1), P(2), and P(3)
!    X0,Y0,Z0 = Dummy arguments for DET
!    X1,Y1,Z1 = Dummy arguments for DET
!    X2,Y2,Z2 = Dummy arguments for DET
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) det
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ), save :: ix = 1
  integer ( kind = 4 ), save :: iy = 2
  integer ( kind = 4 ), save :: iz = 3
  integer ( kind = 4 ) jrand
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1s
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2s
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nst
  real ( kind = 8 ) p(3)
  real ( kind = 8 ) ptn1
  real ( kind = 8 ) ptn2
  real ( kind = 8 ) q(3)
  real ( kind = 8 ) s12
  real ( kind = 8 ) store
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xp
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) yp
  real ( kind = 8 ) z(n)
  real ( kind = 8 ) z0
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2
  real ( kind = 8 ) zp
!
!  Statement function:
!
!  DET(X1,...,Z0) >= 0 if and only if (X0,Y0,Z0) is in the
!  (closed) left hemisphere defined by the plane containing (0,0,0),
!  (X1,Y1,Z1), and (X2,Y2,Z2), where left is defined relative to an
!  observer at (X1,Y1,Z1) facing (X2,Y2,Z2).
!
  det (x1,y1,z1,x2,y2,z2,x0,y0,z0) = x0*(y1*z2-y2*z1) &
       - y0*(x1*z2-x2*z1) + z0*(x1*y2-x2*y1)
!
!  Initialize variables.
!
  xp = p(1)
  yp = p(2)
  zp = p(3)
  n0 = nst

  if ( n0 < 1 .or. n < n0 ) then
    n0 = jrand ( n, ix, iy, iz )
  end if
!
!  Compute the relative machine precision EPS and TOL.
!
  eps = epsilon ( eps )
  tol = 100.0D+00 * eps
!
!  Set NF and NL to the first and last neighbors of N0, and initialize N1 = NF.
!
2 continue

  lp = lend(n0)
  nl = list(lp)
  lp = lptr(lp)
  nf = list(lp)
  n1 = nf
!
!  Find a pair of adjacent neighbors N1,N2 of N0 that define
!  a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.
!
  if ( 0 < nl ) then
!
!  N0 is an interior node.  Find N1.
!
3   continue

    if ( det ( x(n0),y(n0),z(n0),x(n1),y(n1),z(n1),xp,yp,zp ) < 0.0D+00 ) then
      lp = lptr(lp)
      n1 = list(lp)
      if ( n1 == nl ) then
        go to 6
      end if
      go to 3
    end if

  else
!
!  N0 is a boundary node.  Test for P exterior.
!
    nl = -nl
!
!  Is P to the right of the boundary edge N0->NF?
!
    if ( det(x(n0),y(n0),z(n0),x(nf),y(nf),z(nf), xp,yp,zp) < 0.0D+00 ) then
      n1 = n0
      n2 = nf
      go to 9
    end if
!
!  Is P to the right of the boundary edge NL->N0?
!
    if ( det(x(nl),y(nl),z(nl),x(n0),y(n0),z(n0),xp,yp,zp) < 0.0D+00 ) then
      n1 = nl
      n2 = n0
      go to 9
    end if

  end if
!
!  P is to the left of arcs N0->N1 and NL->N0.  Set N2 to the
!  next neighbor of N0 (following N1).
!
4 continue

    lp = lptr(lp)
    n2 = abs ( list(lp) )

    if ( det(x(n0),y(n0),z(n0),x(n2),y(n2),z(n2),xp,yp,zp) < 0.0D+00 ) then
      go to 7
    end if

    n1 = n2

    if ( n1 /= nl ) then
      go to 4
    end if

  if ( det ( x(n0), y(n0), z(n0), x(nf), y(nf), z(nf), xp, yp, zp ) &
    < 0.0D+00 ) then
    go to 6
  end if
!
!  P is left of or on arcs N0->NB for all neighbors NB
!  of N0.  Test for P = +/-N0.
!
  if ( store ( abs ( x(n0 ) * xp + y(n0) * yp + z(n0) * zp) ) &
    < 1.0D+00 - 4.0D+00 * eps ) then
!
!  All points are collinear iff P Left NB->N0 for all
!  neighbors NB of N0.  Search the neighbors of N0.
!  Note:  N1 = NL and LP points to NL.
!
    do

      if ( det(x(n1),y(n1),z(n1),x(n0),y(n0),z(n0),xp,yp,zp) < 0.0D+00 ) then
        exit
      end if

      lp = lptr(lp)
      n1 = abs ( list(lp) )

      if ( n1 == nl ) then
        i1 = 0
        i2 = 0
        i3 = 0
        return
      end if

    end do

  end if
!
!  P is to the right of N1->N0, or P = +/-N0.  Set N0 to N1 and start over.
!
  n0 = n1
  go to 2
!
!  P is between arcs N0->N1 and N0->NF.
!
6 continue

  n2 = nf
!
!  P is contained in a wedge defined by geodesics N0-N1 and
!  N0-N2, where N1 is adjacent to N2.  Save N1 and N2 to
!  test for cycling.
!
7 continue

  n3 = n0
  n1s = n1
  n2s = n2
!
!  Top of edge-hopping loop:
!
8 continue

  b3 = det ( x(n1),y(n1),z(n1),x(n2),y(n2),z(n2),xp,yp,zp )

  if ( b3 < 0.0D+00 ) then
!
!  Set N4 to the first neighbor of N2 following N1 (the
!  node opposite N2->N1) unless N1->N2 is a boundary arc.
!
    lp = lstptr ( lend(n2), n1, list, lptr )

    if ( list(lp) < 0 ) then
      go to 9
    end if

    lp = lptr(lp)
    n4 = abs ( list(lp) )
!
!  Define a new arc N1->N2 which intersects the geodesic N0-P.
!
    if ( det ( x(n0),y(n0),z(n0),x(n4),y(n4),z(n4),xp,yp,zp ) < 0.0D+00 ) then
      n3 = n2
      n2 = n4
      n1s = n1
      if ( n2 /= n2s .and. n2 /= n0 ) then
        go to 8
      end if
    else
      n3 = n1
      n1 = n4
      n2s = n2
      if ( n1 /= n1s .and. n1 /= n0 ) then
        go to 8
      end if
    end if
!
!  The starting node N0 or edge N1-N2 was encountered
!  again, implying a cycle (infinite loop).  Restart
!  with N0 randomly selected.
!
    n0 = jrand ( n, ix, iy, iz )
    go to 2

  end if
!
!  P is in (N1,N2,N3) unless N0, N1, N2, and P are collinear
!  or P is close to -N0.
!
  if ( b3 >= eps ) then
!
!  B3 /= 0.
!
    b1 = det(x(n2),y(n2),z(n2),x(n3),y(n3),z(n3),xp,yp,zp)
    b2 = det(x(n3),y(n3),z(n3),x(n1),y(n1),z(n1),xp,yp,zp)
!
!  Restart with N0 randomly selected.
!
    if ( b1 < -tol .or. b2 < -tol ) then
      n0 = jrand ( n, ix, iy, iz )
      go to 2
    end if

  else
!
!  B3 = 0 and thus P lies on N1->N2. Compute
!  B1 = Det(P,N2 X N1,N2) and B2 = Det(P,N1,N2 X N1).
!
    b3 = 0.0D+00
    s12 = x(n1) * x(n2) + y(n1) * y(n2) + z(n1) * z(n2)
    ptn1 = xp * x(n1) + yp * y(n1) + zp * z(n1)
    ptn2 = xp * x(n2) + yp * y(n2) + zp * z(n2)
    b1 = ptn1 - s12 * ptn2
    b2 = ptn2 - s12 * ptn1
!
!  Restart with N0 randomly selected.
!
    if ( b1 < -tol .or. b2 < -tol ) then
      n0 = jrand ( n, ix, iy, iz )
      go to 2
    end if

  end if
!
!  P is in (N1,N2,N3).
!
  i1 = n1
  i2 = n2
  i3 = n3
  b1 = max ( b1, 0.0D+00 )
  b2 = max ( b2, 0.0D+00 )
  return
!
!  P Right N1->N2, where N1->N2 is a boundary edge.
!  Save N1 and N2, and set NL = 0 to indicate that
!  NL has not yet been found.
!
9 continue

  n1s = n1
  n2s = n2
  nl = 0
!
!  Counterclockwise Boundary Traversal:
!
10 continue

  lp = lend(n2)
  lp = lptr(lp)
  next = list(lp)

  if ( det(x(n2),y(n2),z(n2),x(next),y(next),z(next),xp,yp,zp) >= 0.0D+00 ) then
!
!  N2 is the rightmost visible node if P Forward N2->N1
!  or NEXT Forward N2->N1.  Set Q to (N2 X N1) X N2.
!
    s12 = x(n1) * x(n2) + y(n1) * y(n2) + z(n1) * z(n2)

    q(1) = x(n1) - s12 * x(n2)
    q(2) = y(n1) - s12 * y(n2)
    q(3) = z(n1) - s12 * z(n2)

    if ( xp * q(1) + yp * q(2) + zp * q(3) >= 0.0D+00 ) then
      go to 11
    end if

    if ( x(next) * q(1) + y(next) * q(2) + z(next) * q(3) >= 0.0D+00 ) then
      go to 11
    end if
!
!  N1, N2, NEXT, and P are nearly collinear, and N2 is
!  the leftmost visible node.
!
    nl = n2
  end if
!
!  Bottom of counterclockwise loop:
!
  n1 = n2
  n2 = next

  if ( n2 /= n1s ) then
    go to 10
  end if
!
!  All boundary nodes are visible from P.
!
  i1 = n1s
  i2 = n1s
  i3 = 0
  return
!
!  N2 is the rightmost visible node.
!
11 continue

  nf = n2

  if ( nl == 0 ) then
!
!  Restore initial values of N1 and N2, and begin the search
!  for the leftmost visible node.
!
    n2 = n2s
    n1 = n1s
!
!  Clockwise Boundary Traversal:
!
12  continue

    lp = lend(n1)
    next = -list(lp)

    if ( 0.0D+00 <= &
      det ( x(next), y(next), z(next), x(n1), y(n1), z(n1), xp, yp, zp )  ) then
!
!  N1 is the leftmost visible node if P or NEXT is
!  forward of N1->N2.  Compute Q = N1 X (N2 X N1).
!
      s12 = x(n1) * x(n2) + y(n1) * y(n2) + z(n1) * z(n2)
      q(1) = x(n2) - s12 * x(n1)
      q(2) = y(n2) - s12 * y(n1)
      q(3) = z(n2) - s12 * z(n1)

      if ( xp * q(1) + yp * q(2) + zp * q(3) >= 0.0D+00 ) then
        go to 13
      end if

      if ( x(next) * q(1) + y(next) * q(2) + z(next) * q(3) >= 0.0D+00 ) then
        go to 13
      end if
!
!  P, NEXT, N1, and N2 are nearly collinear and N1 is the rightmost
!  visible node.
!
      nf = n1
    end if
!
!  Bottom of clockwise loop:
!
    n2 = n1
    n1 = next

    if ( n1 /= n1s ) then
      go to 12
    end if
!
!  All boundary nodes are visible from P.
!
    i1 = n1
    i2 = n1
    i3 = 0
    return
!
!  N1 is the leftmost visible node.
!
13   continue

    nl = n1

  end if
!
!  NF and NL have been found.
!
  i1 = nf
  i2 = nl
  i3 = 0

  return
end
subroutine trmesh ( n, x, y, z, list, lptr, lend, ier )

!*****************************************************************************80
!
!! TRMESH creates a Delaunay triangulation on the unit sphere.
!
!  Discussion:
!
!    This subroutine creates a Delaunay triangulation of a
!    set of N arbitrarily distributed points, referred to as
!    nodes, on the surface of the unit sphere.  The Delaunay
!    triangulation is defined as a set of (spherical) triangles
!    with the following five properties:
!
!     1)  The triangle vertices are nodes.
!     2)  No triangle contains a node other than its vertices.
!     3)  The interiors of the triangles are pairwise disjoint.
!     4)  The union of triangles is the convex hull of the set
!           of nodes (the smallest convex set that contains
!           the nodes).  If the nodes are not contained in a
!           single hemisphere, their convex hull is the
!           entire sphere and there are no boundary nodes.
!           Otherwise, there are at least three boundary nodes.
!     5)  The interior of the circumcircle of each triangle
!           contains no node.
!
!    The first four properties define a triangulation, and the
!    last property results in a triangulation which is as close
!    as possible to equiangular in a certain sense and which is
!    uniquely defined unless four or more nodes lie in a common
!    plane.  This property makes the triangulation well-suited
!    for solving closest-point problems and for triangle-based
!    interpolation.
!
!    Provided the nodes are randomly ordered, the algorithm
!    has expected time complexity O(N*log(N)) for most nodal
!    distributions.  Note, however, that the complexity may be
!    as high as O(N**2) if, for example, the nodes are ordered
!    on increasing latitude.
!
!    Spherical coordinates (latitude and longitude) may be
!    converted to Cartesian coordinates by Subroutine TRANS.
!
!    The following is a list of the software package modules
!    which a user may wish to call directly:
!
!    ADDNOD - Updates the triangulation by appending a new node.
!
!    AREAS  - Returns the area of a spherical triangle.
!
!    BNODES - Returns an array containing the indexes of the
!             boundary nodes (if any) in counterclockwise
!             order.  Counts of boundary nodes, triangles,
!             and arcs are also returned.
!
!    CIRCUM - Returns the circumcenter of a spherical triangle.
!
!    CRLIST - Returns the set of triangle circumcenters
!             (Voronoi vertices) and circumradii associated
!             with a triangulation.
!
!    DELARC - Deletes a boundary arc from a triangulation.
!
!    DELNOD - Updates the triangulation with a nodal deletion.
!
!    EDGE   - Forces an arbitrary pair of nodes to be connected
!             by an arc in the triangulation.
!
!    GETNP  - Determines the ordered sequence of L closest nodes
!             to a given node, along with the associated distances.
!
!    INSIDE - Locates a point relative to a polygon on the
!             surface of the sphere.
!
!    INTRSC - Returns the point of intersection between a
!             pair of great circle arcs.
!
!    JRAND  - Generates a uniformly distributed pseudo-random integer.
!
!    LEFT   - Locates a point relative to a great circle.
!
!    NEARND - Returns the index of the nearest node to an
!             arbitrary point, along with its squared
!             distance.
!
!    SCOORD - Converts a point from Cartesian coordinates to
!             spherical coordinates.
!
!    STORE  - Forces a value to be stored in main memory so
!             that the precision of floating point numbers
!             in memory locations rather than registers is
!             computed.
!
!    TRANS  - Transforms spherical coordinates into Cartesian
!             coordinates on the unit sphere for input to
!             Subroutine TRMESH.
!
!    TRLIST - Converts the triangulation data structure to a
!             triangle list more suitable for use in a finite
!             element code.
!
!    TRLPRT - Prints the triangle list created by TRLIST.
!
!    TRMESH - Creates a Delaunay triangulation of a set of
!             nodes.
!
!    TRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
!             file containing a triangulation plot.
!
!    TRPRNT - Prints the triangulation data structure and,
!             optionally, the nodal coordinates.
!
!    VRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
!             file containing a Voronoi diagram plot.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of distinct
!    nodes.  (X(K),Y(K), Z(K)) is referred to as node K, and K is referred
!    to as a nodal index.  It is required that X(K)**2 + Y(K)**2 + Z(K)**2 = 1
!    for all K.  The first three nodes must not be collinear (lie on a
!    common great circle).
!
!    Output, integer ( kind = 4 ) LIST(6*(N-2)), nodal indexes which, along
!    with LPTR, LEND, and LNEW, define the triangulation as a set of N
!    adjacency lists; counterclockwise-ordered sequences of neighboring nodes
!    such that the first and last neighbors of a boundary node are boundary
!    nodes (the first neighbor of an interior node is arbitrary).  In order to
!    distinguish between interior and boundary nodes, the last neighbor of
!    each boundary node is represented by the negative of its index.
!
!    Output, integer ( kind = 4 ) LPTR(6*(N-2)), = Set of pointers (LIST
!    indexes) in one-to-one correspondence with the elements of LIST.
!    LIST(LPTR(I)) indexes the node which follows LIST(I) in cyclical
!    counterclockwise order (the first neighbor follows the last neighbor).
!
!    Output, integer ( kind = 4 ) LEND(N), pointers to adjacency lists.
!    LEND(K) points to the last neighbor of node K.  LIST(LEND(K)) < 0 if and
!    only if K is a boundary node.
!
!    Output, integer ( kind = 4 ) LNEW, pointer to the first empty location
!    in LIST and LPTR (list length plus one).  LIST, LPTR, LEND, and LNEW are
!    not altered if IER < 0, and are incomplete if 0 < IER.
!
!    Workspace, integer ( kind = 4 ) NEAR(N),
!    used to efficiently determine the nearest triangulation node to each
!    unprocessed node for use by ADDNOD.
!
!    Workspace, integer ( kind = 4 ) NEXT(N),
!    used to efficiently determine the nearest triangulation node to each
!    unprocessed node for use by ADDNOD.
!
!    Workspace, real ( kind = 8 ) DIST(N),
!    used to efficiently determine the nearest triangulation node to each
!    unprocessed node for use by ADDNOD.
!
!    Output, integer ( kind = 4 ) IER, error indicator:
!     0, if no errors were encountered.
!    -1, if N < 3 on input.
!    -2, if the first three nodes are collinear.
!     L, if nodes L and M coincide for some L < M.  The data structure
!      represents a triangulation of nodes 1 to M-1 in this case.
!
!  Local parameters:
!
!    D =        (Negative cosine of) distance from node K to node I
!    D1,D2,D3 = Distances from node K to nodes 1, 2, and 3, respectively
!    I,J =      Nodal indexes
!    I0 =       Index of the node preceding I in a sequence of
!               unprocessed nodes:  I = NEXT(I0)
!    K =        Index of node to be added and DO-loop index: 3 < K
!    LP =       LIST index (pointer) of a neighbor of K
!    LPL =      Pointer to the last neighbor of K
!    NEXTI =    NEXT(I)
!    NN =       Local copy of N
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) dist(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical              left
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) near(n)
  integer ( kind = 4 ) next(n)
  integer ( kind = 4 ) nexti
  integer ( kind = 4 ) nn
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  nn = n

  if ( nn < 3 ) then
    ier = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'sTRMESH - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop
  end if
!
!  Store the first triangle in the linked list.
!
  if ( .not. left (x(1),y(1),z(1),x(2),y(2),z(2), &
                   x(3),y(3),z(3) ) ) then
!
!  The first triangle is (3,2,1) = (2,1,3) = (1,3,2).
!
    list(1) = 3
    lptr(1) = 2
    list(2) = -2
    lptr(2) = 1
    lend(1) = 2

    list(3) = 1
    lptr(3) = 4
    list(4) = -3
    lptr(4) = 3
    lend(2) = 4

    list(5) = 2
    lptr(5) = 6
    list(6) = -1
    lptr(6) = 5
    lend(3) = 6

  else if ( .not. left ( x(2),y(2),z(2),x(1),y(1),z(1),x(3),y(3),z(3) ) ) then
!
!  The first triangle is (1,2,3):  3 Strictly Left 1->2,
!  i.e., node 3 lies in the left hemisphere defined by arc 1->2.
!
    list(1) = 2
    lptr(1) = 2
    list(2) = -3
    lptr(2) = 1
    lend(1) = 2

    list(3) = 3
    lptr(3) = 4
    list(4) = -1
    lptr(4) = 3
    lend(2) = 4

    list(5) = 1
    lptr(5) = 6
    list(6) = -2
    lptr(6) = 5
    lend(3) = 6
!
!  The first three nodes are collinear.
!
  else

    ier = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'sTRMESH - Fatal error!'
    write ( *, '(a)' ) '  The first 3 nodes are collinear.'
    write ( *, '(a)' ) '  Try reordering the data.'
    return

  end if
!
!  Initialize LNEW and test for N = 3.
!
  lnew = 7

  if ( nn == 3 ) then
    ier = 0
    return
  end if
!
!  A nearest-node data structure (NEAR, NEXT, and DIST) is
!  used to obtain an expected-time (N*log(N)) incremental
!  algorithm by enabling constant search time for locating
!  each new node in the triangulation.
!
!  For each unprocessed node K, NEAR(K) is the index of the
!  triangulation node closest to K (used as the starting
!  point for the search in Subroutine TRFIND) and DIST(K)
!  is an increasing function of the arc length (angular
!  distance) between nodes K and NEAR(K):  -Cos(a) for arc
!  length a.
!
!  Since it is necessary to efficiently find the subset of
!  unprocessed nodes associated with each triangulation
!  node J (those that have J as their NEAR entries), the
!  subsets are stored in NEAR and NEXT as follows:  for
!  each node J in the triangulation, I = NEAR(J) is the
!  first unprocessed node in J's set (with I = 0 if the
!  set is empty), L = NEXT(I) (if 0 < I) is the second,
!  NEXT(L) (if 0 < L) is the third, etc.  The nodes in each
!  set are initially ordered by increasing indexes (which
!  maximizes efficiency) but that ordering is not main-
!  tained as the data structure is updated.
!
!  Initialize the data structure for the single triangle.
!
  near(1) = 0
  near(2) = 0
  near(3) = 0

  do k = nn, 4, -1

    d1 = -( x(k) * x(1) + y(k) * y(1) + z(k) * z(1) )
    d2 = -( x(k) * x(2) + y(k) * y(2) + z(k) * z(2) )
    d3 = -( x(k) * x(3) + y(k) * y(3) + z(k) * z(3) )

    if ( d1 <= d2 .and. d1 <= d3 ) then
      near(k) = 1
      dist(k) = d1
      next(k) = near(1)
      near(1) = k
    else if ( d2 <= d1 .and. d2 <= d3 ) then
      near(k) = 2
      dist(k) = d2
      next(k) = near(2)
      near(2) = k
    else
      near(k) = 3
      dist(k) = d3
      next(k) = near(3)
      near(3) = k
    end if

  end do
!
!  Add the remaining nodes.
!
  do k = 4, nn

    call addnod ( near(k), k, x, y, z, list, lptr, lend, lnew, ier )

    if ( ier /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'sTRMESH - Fatal error!'
      write ( *, '(a,i8)' ) '  ADDNOD returned error code IER = ', ier
      return
    end if
!
!  Remove K from the set of unprocessed nodes associated with NEAR(K).
!
    i = near(k)

    if ( near(i) == k ) then

      near(i) = next(k)

    else

      i = near(i)

      do

        i0 = i
        i = next(i0)

        if ( i == k ) then
          exit
        end if

      end do

      next(i0) = next(k)

    end if

    near(k) = 0
!
!  Loop on neighbors J of node K.
!
    lpl = lend(k)
    lp = lpl

3   continue

    lp = lptr(lp)
    j = abs ( list(lp) )
!
!  Loop on elements I in the sequence of unprocessed nodes
!  associated with J:  K is a candidate for replacing J
!  as the nearest triangulation node to I.  The next value
!  of I in the sequence, NEXT(I), must be saved before I
!  is moved because it is altered by adding I to K's set.
!
    i = near(j)

    do

      if ( i == 0 ) then
        exit
      end if

      nexti = next(i)
!
!  Test for the distance from I to K less than the distance
!  from I to J.
!
      d = - ( x(i) * x(k) + y(i) * y(k) + z(i) * z(k) )
      if ( d < dist(i) ) then
!
!  Replace J by K as the nearest triangulation node to I:
!  update NEAR(I) and DIST(I), and remove I from J's set
!  of unprocessed nodes and add it to K's set.
!
        near(i) = k
        dist(i) = d

        if ( i == near(j) ) then
          near(j) = nexti
        else
          next(i0) = nexti
        end if

        next(i) = near(k)
        near(k) = i
      else
        i0 = i
      end if

      i = nexti

    end do
!
!  Bottom of loop on neighbors J.
!

    if ( lp /= lpl ) then
      go to 3
    end if

  end do

  return
end


SUBROUTINE INTERP (N,ORDER,PLAT,PLON,X,Y,Z,W,LIST,LPTR,&
                         LEND, IST, PW, IER)
  integer, intent(in) :: order
  INTEGER( kind = 4 ) N, LIST(6*(N-2)),LPTR(6*(N-2)),LEND(N), IST, IER
  REAL( kind = 8 )    PLAT, PLON, X(N), Y(N), Z(N), W(N), PW, GRAD(3,N)
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
!       NNEIGHBOR = logical, if .true. nearest neighbor
!                   interpolation is done
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
  INTEGER( kind = 4 ) I1, I2, I3, LP, N1, N2
  REAL( kind = 8 )    B1, B2, B3, P(3), PTN1, PTN2, S12, SUMM,&
             DIST(3), P1(3), P2(3), P3(3), I(1), ARCLEN
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
! SUMM =      Quantity used to normalize the barycentric
!              coordinates
!
      if (order == 3) then
        CALL INTRC1 (N,PLAT,PLON,X,Y,Z,W,LIST,LPTR,LEND, &
                     0,GRAD, IST, PW,IER)
        RETURN
      endif
      IF (N .LT. 3  .OR.  IST .LT. 1  .OR.  IST .GT. N) GO TO 11
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
                  I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 12
      IST = I1
      IF (I3 .NE. 0) THEN
!
      if (order == 0) then
          P1(1) = X(I1)
          P1(2) = Y(I1)
          P1(3) = Z(I1)
          P2(1) = X(I2)
          P2(2) = Y(I2)
          P2(3) = Z(I2)
          P3(1) = X(I3)
          P3(2) = Y(I3)
          P3(3) = Z(I3)
          dist(1) =  ARCLEN (P,P1)
          dist(2) =  ARCLEN (P,P2)
          dist(3) =  ARCLEN (P,P3)
          I = minloc(dist)
          if (i(1) .eq. 1) pw = w(i1)
          if (i(1) .eq. 2) pw = w(i2)
          if (i(1) .eq. 3) pw = w(i3)
          IER = 0
          RETURN
      endif
!
! P is contained in the triangle (I1,I2,I3).  Normalize the
!   barycentric coordinates.
!
        SUMM = B1 + B2 + B3
        B1 = B1/SUMM
        B2 = B2/SUMM
        B3 = B3/SUMM
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
        SUMM = B1 + B2
        PW = (B1*W(N1) + B2*W(N2))/SUMM
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
      REAL( kind = 8 ) FUNCTION ARCLEN (P,Q)
      REAL( kind = 8 ) P(3), Q(3)
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
      INTEGER( kind = 4 ) I
      REAL( kind = 8 )    D
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
      SUBROUTINE APLYR (X,Y,Z,CX,SX,CY,SY, XP,YP,ZP)
      REAL( kind = 8 ) X, Y, Z, CX, SX, CY, SY, XP, YP, ZP
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE APPLIES THE ROTATION R DEFINED BY SUB-
! ROUTINE CONSTR TO THE UNIT VECTOR (X Y Z)**T, I.E. (X,Y,Z)
! IS ROTATED TO (XP,YP,ZP).  IF (XP,YP,ZP) LIES IN THE
! SOUTHERN HEMISPHERE (ZP .LT. 0), (XP,YP) ARE SET TO THE
! COORDINATES OF THE NEAREST POINT OF THE EQUATOR, ZP RE-
! MAINING UNCHANGED.
!
! INPUT PARAMETERS - X,Y,Z - COORDINATES OF A POINT ON THE
!                            UNIT SPHERE.
!
!              CX,SX,CY,SY - ELEMENTS OF THE ROTATION DE-
!                            FINED BY CONSTR.
!
! INPUT PARAMETERS ARE NOT ALTERED EXCEPT AS NOTED BELOW.
!
! OUTPUT PARAMETERS - XP,YP,ZP - COORDINATES OF THE ROTATED
!                                POINT ON THE SPHERE UNLESS
!                                ZP .LT. 0, IN WHICH CASE
!                                (XP,YP,0) IS THE CLOSEST
!                                POINT OF THE EQUATOR TO THE
!                                ROTATED POINT.  STORAGE FOR
!                                XP, YP, AND ZP MAY COINCIDE
!                                WITH STORAGE FOR X, Y, AND
!                                Z, RESPECTIVELY, IF THE
!                                LATTER NEED NOT BE SAVED.
!
! MODULES REFERENCED BY APLYR - NONE
!
! INTRINSIC FUNCTION CALLED BY APLYR - SQRT
!
!***********************************************************
!
      REAL( kind = 8 ) T
!
! LOCAL PARAMETER -
!
! T = TEMPORARY VARIABLE
!
      T = SX*Y + CX*Z
      YP = CX*Y - SX*Z
      ZP = SY*X + CY*T
      XP = CY*X - SY*T
      IF (ZP .GE. 0.) RETURN
!
! MOVE (XP,YP,ZP) TO THE EQUATOR
!
      T = SQRT(XP*XP + YP*YP)
      IF (T .EQ. 0.) GO TO 1
      XP = XP/T
      YP = YP/T
      RETURN
!
! MOVE THE SOUTH POLE TO AN ARBITRARY POINT OF THE EQUATOR
!
    1 XP = 1.
      YP = 0.
      RETURN
      END
      SUBROUTINE APLYRT (G1P,G2P,CX,SX,CY,SY, G)
      REAL( kind = 8 ) G1P, G2P, CX, SX, CY, SY, G(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE APPLIES THE INVERSE (TRANSPOSE) OF THE
! ROTATION DEFINED BY SUBROUTINE CONSTR TO THE VECTOR
! (G1P G2P 0)**T, I.E. THE GRADIENT (G1P,G2P,0) IN THE ROT-
! ATED COORDINATE SYSTEM IS MAPPED TO (G1,G2,G3) IN THE
! ORIGINAL COORDINATE SYSTEM.
!
! INPUT PARAMETERS - G1P,G2P - X- AND Y-COMPONENTS, RESPECT-
!                              IVELY, OF THE GRADIENT IN THE
!                              ROTATED COORDINATE SYSTEM.
!
!                CX,SX,CY,SY - ELEMENTS OF THE ROTATION R
!                              CONSTRUCTED BY SUBROUTINE
!                              CONSTR.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - G - X-, Y-, AND Z-COMPONENTS (IN THAT
!                         ORDER) OF THE INVERSE ROTATION
!                         APPLIED TO (G1P,G2P,0) -- GRADIENT
!                         IN THE ORIGINAL COORDINATE SYSTEM.
!
! MODULES REFERENCED BY APLYRT - NONE
!
!***********************************************************
!
      REAL( kind = 8 ) T
!
! LOCAL PARAMETERS -
!
! T = TEMPORARY VARIABLE
!
      T = SY*G1P
      G(1) = CY*G1P
      G(2) = CX*G2P - SX*T
      G(3) = -SX*G2P - CX*T
      RETURN
      END
      SUBROUTINE GIVENS ( A,B, C,S)
      REAL( kind = 8 ) A, B, C, S
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE CONSTRUCTS THE GIVENS PLANE ROTATION --
!     ( C  S)
! G = (     ) WHERE C*C + S*S = 1 -- WHICH ZEROS THE SECOND
!     (-S  C)
! ENTRY OF THE 2-VECTOR (A B)-TRANSPOSE.  A CALL TO GIVENS
! IS NORMALLY FOLLOWED BY A CALL TO ROTATE WHICH APPLIES
! THE TRANSFORMATION TO A 2 BY N MATRIX.  THIS ROUTINE WAS
! TAKEN FROM LINPACK.
!
! INPUT PARAMETERS - A,B - COMPONENTS OF THE 2-VECTOR TO BE
!                          ROTATED.
!
! OUTPUT PARAMETERS -  A - OVERWRITTEN BY R = +/-SQRT(A*A
!                          + B*B)
!
!                      B - OVERWRITTEN BY A VALUE Z WHICH
!                          ALLOWS C AND S TO BE RECOVERED
!                          AS FOLLOWS -
!                          C = SQRT(1-Z*Z), S=Z IF ABS(Z)
!                              .LE. 1.
!                          C = 1/Z, S = SQRT(1-C*C) IF
!                              ABS(Z) .GT. 1.
!
!                      C - +/-(A/R)
!
!                      S - +/-(B/R)
!
! MODULES REFERENCED BY GIVENS - NONE
!
! INTRINSIC FUNCTIONS CALLED BY GIVENS - ABS, SQRT
!
!***********************************************************
!
      REAL( kind = 8 ) AA, BB, R, U, V
!
! LOCAL PARAMETERS -
!
! AA,BB = LOCAL COPIES OF A AND B
! R =     C*A + S*B = +/-SQRT(A*A+B*B)
! U,V =   VARIABLES USED TO SCALE A AND B FOR COMPUTING R
!
      AA = A
      BB = B
      IF (ABS(AA) .LE. ABS(BB)) GO TO 1
!
! ABS(A) .GT. ABS(B)
!
      U = AA + AA
      V = BB/U
      R = SQRT(.25 + V*V) * U
      C = AA/R
      S = V * (C + C)
!
! NOTE THAT R HAS THE SIGN OF A, C .GT. 0, AND S HAS
!   SIGN(A)*SIGN(B)
!
      B = S
      A = R
      RETURN
!
! ABS(A) .LE. ABS(B)
!
    1 IF (BB .EQ. 0.) GO TO 2
      U = BB + BB
      V = AA/U
!
! STORE R IN A
!
      A = SQRT(.25 + V*V) * U
      S = BB/A
      C = V * (S + S)
!
! NOTE THAT R HAS THE SIGN OF B, S .GT. 0, AND C HAS
!   SIGN(A)*SIGN(B)
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
      SUBROUTINE ROTATE (N,C,S, X,Y )
      INTEGER( kind = 4) N
      REAL( kind = 8 )    C, S, X(N), Y(N)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!                                            ( C  S)
!   THIS ROUTINE APPLIES THE GIVENS ROTATION (     ) TO THE
!                                            (-S  C)
!               (X(1) ... X(N))
! 2 BY N MATRIX (             ).  THIS ROUTINE WAS TAKEN
!               (Y(1) ... Y(N))
! FROM LINPACK.
!
! INPUT PARAMETERS -   N - NUMBER OF COLUMNS TO BE ROTATED.
!
!                    C,S - ELEMENTS OF THE GIVENS ROTATION.
!                          THESE MAY BE DETERMINED BY
!                          SUBROUTINE GIVENS.
!
!                    X,Y - VECTORS OF LENGTH .GE. N
!                          CONTAINING THE 2-VECTORS TO BE
!                          ROTATED.
!
!   THE PARAMETERS N, C, AND S ARE NOT ALTERED BY THIS
! ROUTINE.
!
! OUTPUT PARAMETERS - X,Y - ROTATED VECTORS
!
! MODULES REFERENCED BY ROTATE - NONE
!
!***********************************************************
!
      INTEGER( kind = 4) I
      REAL( kind = 8 )    XI, YI
!
! LOCAL PARAMETERS -
!
! I =     DO-LOOP INDEX
! XI,YI = X(I), Y(I)
!
      IF (N .LE. 0 .OR. (C .EQ. 1. .AND. S .EQ. 0.)) RETURN
      DO 1 I = 1,N
        XI = X(I)
        YI = Y(I)
        X(I) = C*XI + S*YI
    1   Y(I) = -S*XI + C*YI
      RETURN
      END
      SUBROUTINE SETUP (XI,YI,WI,WK,S1,S2,WT, ROW)
      REAL( kind = 8 ) XI, YI, WI, WK, S1, S2, WT, ROW(6)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE SETS UP THE I-TH ROW OF AN AUGMENTED
! REGRESSION MATRIX FOR A WEIGHTED LEAST SQUARES FIT OF A
! QUADRATIC FUNCTION Q(X,Y) TO A SET OF DATA VALUES WI
! WHERE Q(0,0) = WK.  THE FIRST 3 COLUMNS (QUADRATIC TERMS)
! ARE SCALED BY 1/S2 AND THE FOURTH AND FIFTH COLUMNS (LIN-
! EAR TERMS) ARE SCALED BY 1/S1.
!
! INPUT PARAMETERS - XI,YI - COORDINATES OF NODE I.
!
!                       WI - DATA VALUE AT NODE I.
!
!                       WK - DATA VALUE INTERPOLATED BY Q AT
!                            THE ORIGIN.
!
!                    S1,S2 - INVERSE SCALE FACTORS.
!
!                       WT - WEIGHT FACTOR CORRESPONDING TO
!                            THE I-TH EQUATION.
!
!                      ROW - VECTOR OF LENGTH 6.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER - ROW - VECTOR CONTAINING A ROW OF THE
!                          AUGMENTED REGRESSION MATRIX.
!
! MODULES REFERENCED BY SETUP - NONE
!
!***********************************************************
!
      REAL( kind = 8 ) W1, W2
!
! LOCAL PARAMETERS -
!
! W1 = WEIGHTED SCALE FACTOR FOR THE LINEAR TERMS
! W2 = WEIGHTED SCALE FACTOR FOR THE QUADRATIC TERMS
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
      SUBROUTINE CONSTR (XK,YK,ZK, CX,SX,CY,SY)
      REAL( kind = 8 ) XK, YK, ZK, CX, SX, CY, SY
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE CONSTRUCTS THE ELEMENTS OF A 3 BY 3
! ORTHOGONAL MATRIX R WHICH ROTATES A POINT (XK,YK,ZK) ON
! THE UNIT SPHERE TO THE NORTH POLE, I.E.
!
!      (XK)     (CY  0 -SY)   (1   0   0)   (XK)     (0)
!  R * (YK)  =  ( 0  1   0) * (0  CX -SX) * (YK)  =  (0)
!      (ZK)     (SY  0  CY)   (0  SX  CX)   (ZK)     (1)
!
! INPUT PARAMETERS - XK,YK,ZK - COMPONENTS OF A UNIT VECTOR
!                               TO BE ROTATED TO (0,0,1).
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - CX,SX,CY,SY - ELEMENTS OF R -- CX,SX
!                                   DEFINE A ROTATION ABOUT
!                                   THE X-AXIS AND CY,SY DE-
!                                   FINE A ROTATION ABOUT
!                                   THE Y-AXIS.
!
! MODULES REFERENCED BY CONSTR - NONE
!
! INTRINSIC FUNCTION CALLED BY CONSTR - SQRT
!
!***********************************************************
!
      CY = SQRT(YK*YK + ZK*ZK)
      SY = XK
      IF (CY .EQ. 0.) GO TO 1
      CX = ZK/CY
      SX = YK/CY
      RETURN
!
! (XK,YK,ZK) LIES ON THE X-AXIS
!
    1 CX = 1.
      SX = 0.
      RETURN
      END
      SUBROUTINE ARCINT (P,P1,P2,W1,W2,G1,G2, W,G,GN)
      REAL( kind = 8 )    P(3), P1(3), P2(3), W1, W2, G1(3), G2(3), &
              W, G(3), GN
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN 3 POINTS P, P1, AND P2 LYING ON A COMMON GEODESIC
! OF THE UNIT SPHERE WITH P BETWEEN P1 AND P2, ALONG WITH
! DATA VALUES AND GRADIENTS AT P1 AND P2, THIS SUBROUTINE
! COMPUTES AN INTERPOLATED VALUE W AND A GRADIENT VECTOR G
! AT P.  W IS COMPUTED BY HERMITE CUBIC INTERPOLATION REL-
! ATIVE TO ARC-LENGTH ALONG THE GEODESIC.  THE TANGENTIAL
! COMPONENT OF G IS THE DERIVATIVE (WITH RESPECT TO ARC-
! LENGTH) OF THE CUBIC INTERPOLANT AT P, WHILE THE NORMAL
! COMPONENT OF G IS OBTAINED BY LINEAR INTERPOLATION OF THE
! NORMAL COMPONENTS OF THE GRADIENTS AT P1 AND P2.  THIS
! ALGORITHM IS DUE TO C. L. LAWSON.
!
! INPUT PARAMETERS - P - CARTESIAN COORDINATES OF A POINT
!                        LYING ON THE ARC DEFINED BY P1 AND
!                        P2.
!
!                P1,P2 - COORDINATES OF DISTINCT POINTS ON
!                        THE UNIT SPHERE DEFINING AN ARC
!                        WITH LENGTH LESS THAN 180 DEGREES.
!
!                W1,W2 - DATA VALUES ASSOCIATED WITH P1 AND
!                        P2, RESPECTIVELY.
!
!                G1,G2 - GRADIENT VECTORS ASSOCIATED WITH P1
!                        AND P2.  G1 AND G2 ARE ORTHOGONAL
!                        TO P1 AND P2, RESPECTIVELY.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                    G - ARRAY OF LENGTH 3.
!
! OUTPUT PARAMETERS - W - INTERPOLATED VALUE AT P.
!
!                     G - INTERPOLATED GRADIENT AT P.
!
!                    GN - NORMAL COMPONENT OF G WITH THE
!                         DIRECTION P1 X P2 TAKEN TO BE
!                         POSITIVE.  THE EXTRAPOLATION
!                         PROCEDURE REQUIRES THIS COMPONENT.
!
! FOR EACH VECTOR V, V(1), V(2), AND V(3) CONTAIN X-, Y-,
!   AND Z-COMPONENTS, RESPECTIVELY.
!
! MODULES REFERENCED BY ARCINT - ARCLEN
!
! INTRINSIC FUNCTION CALLED BY ARCINT - SQRT
!
!***********************************************************
!
      INTEGER( kind = 4) I, LUN
      REAL( kind = 8 )    UN(3), UNORM, TAU1, TAU2, A, AL, S, T, GT,&
              ARCLEN
      DATA    LUN/6/
!
! LOCAL PARAMETERS -
!
! I =         DO-LOOP INDEX
! LUN =       LOGICAL UNIT FOR ERROR MESSAGES
! UN =        UNIT NORMAL TO THE PLANE OF P, P1, AND P2
! UNORM =     EUCLIDEAN NORM OF P1 X P2 -- USED TO NORMALIZE
!               UN
! TAU1,TAU2 = TANGENTIAL DERIVATIVES (COMPONENTS OF G1,G2)
!               AT P1 AND P2
! A =         ANGLE IN RADIANS (ARC-LENGTH) BETWEEN P1 AND
!               P2
! AL =        ARC-LENGTH BETWEEN P1 AND P
! S =         NORMALIZED VALUE OF AL -- AS P VARIES FROM P1
!               TO P2, S VARIES FROM 0 TO 1
! T =         1-S -- S AND T ARE BARYCENTRIC COORDINATES OF
!               P WITH RESPECT TO THE ARC FROM P1 TO P2
! GT =        TANGENTIAL COMPONENT OF G -- COMPONENT IN THE
!               DIRECTION UN X P
!
!
! COMPUTE UNIT NORMAL UN
!
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.) GO TO 2
!
! NORMALIZE UN
!
      DO 1 I = 1,3
    1   UN(I) = UN(I)/UNORM
!
! COMPUTE TANGENTIAL DERIVATIVES AT THE ENDPOINTS --
!   TAU1 = (G1,UN X P1) = (G1,P2)/UNORM AND
!   TAU2 = (G2,UN X P2) = -(G2,P1)/UNORM.
!
      TAU1 = (G1(1)*P2(1) + G1(2)*P2(2) + G1(3)*P2(3))/UNORM
      TAU2 =-(G2(1)*P1(1) + G2(2)*P1(2) + G2(3)*P1(3))/UNORM
!
! COMPUTE ARC-LENGTHS A, AL
!
      A = ARCLEN(P1,P2)
      IF (A .EQ. 0.) GO TO 2
      AL = ARCLEN(P1,P)
!
! COMPUTE W BY HERMITE CUBIC INTERPOLATION
!
      S = AL/A
      T = 1. - S
      W = W1*(2.*S+1.)*T*T + W2*(3.-2.*S)*S*S + &
          A*S*T*(TAU1*T - TAU2*S)
!
! COMPUTE TANGENTIAL AND NORMAL DERIVATIVES AT P
!
      GT = 6.*S*T*(W2-W1)/A + &
           TAU1*T*(1.-3.*S) + TAU2*S*(3.*S-2.)
      GN = T*(UN(1)*G1(1) + UN(2)*G1(2) + UN(3)*G1(3)) + &
           S*(UN(1)*G2(1) + UN(2)*G2(2) + UN(3)*G2(3))
!
! COMPUTE G = GT*(UN X P) + GN*UN
!
      G(1) = GT*(UN(2)*P(3) - UN(3)*P(2)) + GN*UN(1)
      G(2) = GT*(UN(3)*P(1) - UN(1)*P(3)) + GN*UN(2)
      G(3) = GT*(UN(1)*P(2) - UN(2)*P(1)) + GN*UN(3)
      RETURN
!
! P1 X P2 = 0.  PRINT AN ERROR MESSAGE AND TERMINATE
!   PROCESSING.
!
    2 WRITE (LUN,100) (P1(I),I=1,3), (P2(I),I=1,3)
  100 FORMAT (1H1,24HERROR IN ARCINT -- P1 = ,2(F9.6,3H,  ), &
              F9.6/1H ,19X,5HP2 = ,2(F9.6,3H,  ),F9.6)
      STOP
      END

SUBROUTINE GETNP (x,y,z,list,lptr,lend,l,npts,df,ier)

! Code converted using TO_F90 by Alan Miller
! Date: 2015-12-20  Time: 09:52:20

REAL( kind = 8 ), INTENT(IN)                         :: x(*)
REAL( kind = 8 ), INTENT(IN)                         :: y(*)
REAL( kind = 8 ), INTENT(IN)                         :: z(*)
INTEGER( kind = 4 ), INTENT(IN OUT)                  :: list(*)
INTEGER( kind = 4 ), INTENT(IN)                      :: lptr(*)
INTEGER( kind = 4 ), INTENT(OUT)                     :: lend(*)
INTEGER( kind = 4 ), INTENT(IN)                      :: l
INTEGER( kind = 4 ), INTENT(IN OUT)                  :: npts(l)
REAL( kind = 8 ), INTENT(OUT)                        :: df
INTEGER( kind = 4 ), INTENT(OUT)                     :: ier



!***********************************************************

!                                              From STRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/28/98

!   Given a Delaunay triangulation of N nodes on the unit
! sphere and an array NPTS containing the indexes of L-1
! nodes ordered by angular distance from NPTS(1), this sub-
! routine sets NPTS(L) to the index of the next node in the
! sequence -- the node, other than NPTS(1),...,NPTS(L-1),
! that is closest to NPTS(1).  Thus, the ordered sequence
! of K closest nodes to N1 (including N1) may be determined
! by K-1 calls to GETNP with NPTS(1) = N1 and L = 2,3,...,K
! for K .GE. 2.

!   The algorithm uses the property of a Delaunay triangula-
! tion that the K-th closest node to N1 is a neighbor of one
! of the K-1 closest nodes to N1.


! On input:

!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of the nodes.

!       LIST,LPTR,LEND = Triangulation data structure.  Re-
!                        fer to Subroutine TRMESH.

!       L = Number of nodes in the sequence on output.  2
!           .LE. L .LE. N.

! The above parameters are not altered by this routine.

!       NPTS = Array of length .GE. L containing the indexes
!              of the L-1 closest nodes to NPTS(1) in the
!              first L-1 locations.

! On output:

!       NPTS = Array updated with the index of the L-th
!              closest node to NPTS(1) in position L unless
!              IER = 1.

!       DF = Value of an increasing function (negative cos-
!            ine) of the angular distance between NPTS(1)
!            and NPTS(L) unless IER = 1.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if L < 2.

! Modules required by GETNP:  None

! Intrinsic function called by GETNP:  ABS

!***********************************************************

INTEGER( kind = 4 ) :: i, lm1, lp, lpl, n1, nb, ni, np
REAL( kind = 8 ) :: dnb, dnp, x1, y1, z1

! Local parameters:

! DNB,DNP =  Negative cosines of the angular distances from
!              N1 to NB and to NP, respectively
! I =        NPTS index and DO-loop index
! LM1 =      L-1
! LP =       LIST pointer of a neighbor of NI
! LPL =      Pointer to the last neighbor of NI
! N1 =       NPTS(1)
! NB =       Neighbor of NI and candidate for NP
! NI =       NPTS(I)
! NP =       Candidate for NPTS(L)
! X1,Y1,Z1 = Coordinates of N1

lm1 = l - 1
IF (lm1 < 1) GO TO 6
ier = 0

! Store N1 = NPTS(1) and mark the elements of NPTS.

n1 = npts(1)
x1 = x(n1)
y1 = y(n1)
z1 = z(n1)
DO  i = 1,lm1
  ni = npts(i)
  lend(ni) = -lend(ni)
END DO

! Candidates for NP = NPTS(L) are the unmarked neighbors
!   of nodes in NPTS.  DNP is initially greater than -cos(PI)
!   (the maximum distance).

dnp = 2.

! Loop on nodes NI in NPTS.

DO  i = 1,lm1
  ni = npts(i)
  lpl = -lend(ni)
  lp = lpl

! Loop on neighbors NB of NI.

  2   nb = ABS(list(lp))
  IF (lend(nb) < 0) GO TO 3

! NB is an unmarked neighbor of NI.  Replace NP if NB is
!   closer to N1.

  dnb = -(x(nb)*x1 + y(nb)*y1 + z(nb)*z1)
  IF (dnb >= dnp) GO TO 3
  np = nb
  dnp = dnb
  3     lp = lptr(lp)
  IF (lp /= lpl) GO TO 2
END DO
npts(l) = np
df = dnp

! Unmark the elements of NPTS.

DO  i = 1,lm1
  ni = npts(i)
  lend(ni) = -lend(ni)
END DO
RETURN

! L is outside its valid range.

6 ier = 1
RETURN
END SUBROUTINE getnp

SUBROUTINE INTRC1 (N,PLAT,PLON,X,Y,Z,W,LIST,LPTR,LEND, &
                         IFLAG,GRAD, IST, PW,IER)
      INTEGER( kind = 4) N,LIST(6*(N-2)),LPTR(6*(N-2)),LEND(N),IFLAG,IST,IER
      REAL( kind = 8 )    PLAT, PLON, X(N), Y(N), Z(N), W(N), &
              GRAD(3,N), PW
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF NODES ON THE UNIT
! SPHERE, ALONG WITH DATA VALUES AT THE NODES, THIS SUB-
! ROUTINE CONSTRUCTS THE VALUE AT A POINT P OF A ONCE CON-
! TINUOUSLY DIFFERENTIABLE FUNCTION WHICH INTERPOLATES THE
! DATA VALUES.  IF THE TRIANGULATION DOES NOT COVER THE
! ENTIRE SPHERE, THE SURFACE IS EXTENDED CONTINUOUSLY BEYOND
! THE BOUNDARY ALLOWING EXTRAPOLATION.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES.  N .GE. 3 AND
!                            N .GE. 7 IF IFLAG = 0.
!
!                PLAT,PLON - LATITUDE AND LONGITUDE OF P IN
!                            RADIANS.
!
!                    X,Y,Z - VECTORS CONTAINING CARTESIAN
!                            COORDINATES OF THE NODES.
!
!                        W - VECTOR CONTAINING DATA VALUES
!                            AT THE NODES.  W(I) IS ASSOCI-
!                            ATED WITH (X(I),Y(I),Z(I)) FOR
!                            I = 1,...,N.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
!                        Subroutine TRMESH.
!
!                    IFLAG - OPTION INDICATOR
!                            IFLAG = 0 IF INTRC1 IS TO PRO-
!                                      VIDE ESTIMATED GRAD-
!                                      IENTS (FROM GRADL).
!                                      N .GE. 7 IN THIS
!                                      CASE.
!                            IFLAG = 1 IF GRADIENTS ARE PRO-
!                                      VIDED IN GRAD.  THIS
!                                      IS MORE EFFICIENT IF
!                                      INTRC1 IS TO BE
!                                      CALLED SEVERAL TIMES.
!
!                     GRAD - ARRAY DIMENSIONED 3 BY N WHOSE
!                            I-TH COLUMN CONTAINS AN ESTI-
!                            MATED GRADIENT AT NODE I IF
!                            IFLAG = 1 (SEE SUBROUTINE
!                            GRADL).  GRAD MAY BE A DUMMY
!                            VARIABLE (NOT USED) IF IFLAG
!                            = 0.
!
!                      IST - INDEX OF THE STARTING NODE IN
!                            THE SEARCH FOR A TRIANGLE CON-
!                            TAINING P.  1 .LE. IST .LE. N.
!                            THE OUTPUT VALUE OF IST FROM A
!                            PREVIOUS CALL MAY BE A GOOD
!                            CHOICE.
!
! INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - IST - INDEX OF ONE OF THE VERTICES OF
!                           THE TRIANGLE CONTAINING P (OR
!                           NEAREST P) UNLESS IER = -1 OR
!                           IER = -2.
!
! OUTPUT PARAMETERS -  PW - INTERPOLATED VALUE AT P IF
!                           IER .GE. 0.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF INTERPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = 1 IF EXTRAPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = -1 IF N, IFLAG, OR IST IS
!                                    OUT OF RANGE.
!                           IER = -2 IF THE NODES ARE COL-
!                                    LINEAR.
!                           IER = -3 IF THE ANGULAR DISTANCE
!                                    BETWEEN P AND THE NEAR-
!                                    EST POINT OF THE TRIAN-
!                                    GULATION IS AT LEAST 90
!                                    DEGREES.
!
! MODULES REFERENCED BY INTRC1 - TRFIND, WVAL, ARCINT,
!                                ARCLEN,
!             (AND OPTIONALLY)   GRADL, GETNP, CONSTR,
!                                APLYR, SETUP, GIVENS,
!                                ROTATE, APLYRT
!
! INTRINSIC FUNCTIONS CALLED BY INTRC1 - COS, SIN, SQRT
!
!***********************************************************
!
      INTEGER( kind = 4) NN, I1, I2, I3, I, IERR, N1, N2
      REAL( kind = 8 )    P(3), P1(3), P2(3), P3(3), W1, W2, W3, G1(3), &
              G2(3), G3(3), B1, B2, B3, SUMM, DUM(3), S12, &
              PTN1, PTN2, Q(3), QNORM, WQ, GQ(3), A, PTGQ, &
              GQN, ARCLEN
!
! LOCAL PARAMETERS -
!
! NN =       LOCAL COPY OF N
! I1,I2,I3 = VERTEX INDICES RETURNED BY TRFIND
! I =        DO-LOOP INDEX
! IERR =     ERROR FLAG FOR CALLS TO GRADL
! N1,N2 =    INDICES OF THE ENDPOINTS OF A BOUNDARY ARC WHEN
!              P IS EXTERIOR (NOT CONTAINED IN A TRIANGLE)
! P =        CARTESIAN COORDINATES OF P
! P1,P2,P3 = CARTESIAN COORDINATES OF THE VERTICES I1, I2,
!              AND I3, OR (P1 AND P2) COORDINATES OF N1 AND
!              N2 IF P IS EXTERIOR
! W1,W2,W3 = DATA VALUES ASSOCIATED WITH I1, I2, AND I3, OR
!              (W1 AND W2) VALUES ASSOCIATED WITH N1 AND
!              N2 IF P IS EXTERIOR
! G1,G2,G3 = GRADIENTS AT I1, I2, AND I3, OR (G1 AND G2) AT
!              N1 AND N2
! B1,B2,B3 = BARYCENTRIC COORDINATES OF THE CENTRAL PROJEC-
!              TION OF P ONTO THE UNDERLYING PLANAR TRIANGLE
!              OR (B1 AND B2) PROJECTION OF Q ONTO THE
!              UNDERLYING LINE SEGMENT N1-N2 WHEN P IS
!              EXTERIOR -- UNNORMALIZED COORDINATES ARE
!              COMPUTED BY TRFIND WHEN P IS IN A TRIANGLE
! SUMM =     QUANTITY USED TO NORMALIZE THE BARYCENTRIC
!             COORDINATES
! DUM =      DUMMY PARAMETER FOR WVAL AND ARCINT
! S12 =      SCALAR PRODUCT (N1,N2) -- FACTOR OF B1 AND B2
! PTN1 =     SCALAR PRODUCT (P,N1) -- FACTOR OF B1 AND B2
! PTN2 =     SCALAR PRODUCT (P,N2) -- FACTOR OF B1 AND B2
! Q =        CLOSEST BOUNDARY POINT TO P WHEN P IS EXTERIOR
! QNORM =    FACTOR USED TO NORMALIZE Q
! WQ,GQ =    INTERPOLATED VALUE AND GRADIENT AT Q
! A =        ANGULAR SEPARATION BETWEEN P AND Q
! PTGQ =     SCALAR PRODUCT (P,GQ) -- FACTOR OF THE COMPONENT
!              OF GQ IN THE DIRECTION Q->P
! GQN =      NEGATIVE OF THE COMPONENT OF GQ IN THE DIRECTION
!              Q->P
!
      NN = N
      IF (NN .LT. 3  .OR.  (IFLAG .NE. 1  .AND.  NN .LT. 7) &
         .OR.  IFLAG .LT. 0  .OR.  IFLAG .GT. 1  .OR. &
         IST .LT. 1  .OR.  IST .GT. NN) GO TO 16
!
! TRANSFORM (PLAT,PLON) TO CARTESIAN COORDINATES
!
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
!
! LOCATE P WITH RESPECT TO THE TRIANGULATION
!
      CALL TRFIND(IST,P,NN,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,&
                  I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 17
      IST = I1
      IF (I3 .EQ. 0) GO TO 4
!
! P IS CONTAINED IN THE TRIANGLE (I1,I2,I3).  STORE THE DATA
!   VALUES AND VERTEX COORDINATES IN LOCAL VARIABLES.
!
      W1 = W(I1)
      W2 = W(I2)
      W3 = W(I3)
      P1(1) = X(I1)
      P1(2) = Y(I1)
      P1(3) = Z(I1)
      P2(1) = X(I2)
      P2(2) = Y(I2)
      P2(3) = Z(I2)
      P3(1) = X(I3)
      P3(2) = Y(I3)
      P3(3) = Z(I3)
      IF (IFLAG .NE. 1) GO TO 2
!
! GRADIENTS ARE USER-PROVIDED
!
      DO 1 I = 1,3
        G1(I) = GRAD(I,I1)
        G2(I) = GRAD(I,I2)
    1   G3(I) = GRAD(I,I3)
      GO TO 3
!
! COMPUTE GRADIENT ESTIMATES AT THE VERTICES
!
    2 CALL GRADL(NN,I1,X,Y,Z,W,LIST,LPTR,LEND, G1,IERR)
      IF (IERR .LT. 0) GO TO 17
      CALL GRADL(NN,I2,X,Y,Z,W,LIST,LPTR,LEND, G2,IERR)
      IF (IERR .LT. 0) GO TO 17
      CALL GRADL(NN,I3,X,Y,Z,W,LIST,LPTR,LEND, G3,IERR)
      IF (IERR .LT. 0) GO TO 17
!
! NORMALIZE THE COORDINATES
!
    3 SUMM = B1 + B2 + B3
      B1 = B1/SUMM
      B2 = B2/SUMM
      B3 = B3/SUMM
      CALL WVAL(B1,B2,B3,P1,P2,P3,W1,W2,W3,G1,G2,G3,0, PW, &
                DUM)
      IER = 0
      RETURN
!
! P IS EXTERIOR TO THE TRIANGULATION, AND I1 AND I2 ARE
!   BOUNDARY NODES WHICH ARE VISIBLE FROM P.  EXTRAPOLATE TO
!   P BY LINEAR (WITH RESPECT TO ARC-LENGTH) INTERPOLATION
!   OF THE VALUE AND DIRECTIONAL DERIVATIVE (GRADIENT COMP-
!   ONENT IN THE DIRECTION Q->P) OF THE INTERPOLATORY
!   SURFACE AT Q WHERE Q IS THE CLOSEST BOUNDARY POINT TO P.
!
! DETERMINE Q BY TRAVERSING THE BOUNDARY STARTING FROM I1
!
    4 N1 = I1
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
      IF (I1 .NE. I2) GO TO 6
!
! ALL BOUNDARY NODES ARE VISIBLE FROM P.  FIND A BOUNDARY
!   ARC N1->N2 SUCH THAT P LEFT (N2 X N1)->N1
!
! COUNTERCLOCKWISE BOUNDARY TRAVERSAL --
!   SET N2 TO THE FIRST NEIGHBOR OF N1
!
    5 lp = lend(n1)
      lp = lptr(lp)
      n2 = list(lp)
!
! COMPUTE INNER PRODUCTS (N1,N2) AND (P,N2), AND COMPUTE
!   B2 = DET(P,N1,N2 X N1)
!
      S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
      PTN2 = P(1)*X(N2) + P(2)*Y(N2) + P(3)*Z(N2)
      B2 = PTN2 - S12*PTN1
      IF (B2 .LE. 0.) GO TO 6
!
! P RIGHT (N2 X N1)->N1 -- ITERATE
!
      N1 = N2
      I1 = N1
      PTN1 = PTN2
      GO TO 5
!
! P LEFT (N2 X N1)->N1 WHERE N2 IS THE FIRST NEIGHBOR OF N1
!   CLOCKWISE BOUNDARY TRAVERSAL --
!
    6 N2 = N1
      PTN2 = PTN1
!
! SET N1 TO THE LAST NEIGHBOR OF N2 AND TEST FOR TERMINATION
!
      lp = lend(n2)
      n1 = -list(lp)
      IF (N1 .EQ. I1) GO TO 18
!
! COMPUTE INNER PRODUCTS (N1,N2) AND (P,N1)
!
      S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
!
! COMPUTE B2 = DET(P,N1,N2 X N1) = DET(Q,N1,N2 X N1)*(P,Q)
!
      B2 = PTN2 - S12*PTN1
      IF (B2 .LE. 0.) GO TO 6
!
! COMPUTE B1 = DET(P,N2 X N1,N2) = DET(Q,N2 X N1,N2)*(P,Q)
!
      B1 = PTN1 - S12*PTN2
      IF (B1 .GT. 0.) GO TO 10
!
! Q = N2.  STORE VALUE, COORDINATES, AND AND GRADIENT AT Q.
!
      WQ = W(N2)
      Q(1) = X(N2)
      Q(2) = Y(N2)
      Q(3) = Z(N2)
      IF (IFLAG .NE. 1) GO TO 8
      DO 7 I = 1,3
    7   GQ(I) = GRAD(I,N2)
      GO TO 9
    8 CALL GRADL(NN,N2,X,Y,Z,W,LIST,LPTR,LEND, GQ,IERR)
      IF (IERR .LT. 0) GO TO 17
!
! EXTRAPOLATE TO P -- PW = WQ + A*(GQ,Q X (PXQ)/SIN(A))
!   WHERE A IS THE ANGULAR SEPARATION BETWEEN Q AND P,
!   AND SIN(A) IS THE MAGNITUDE OF P X Q.
!
    9 A = ARCLEN(Q,P)
      PTGQ = P(1)*GQ(1) + P(2)*GQ(2) + P(3)*GQ(3)
      PW = WQ
      IF (A .NE. 0.) PW = PW + PTGQ*A/SIN(A)
      IER = 1
      RETURN
!
! P STRICTLY LEFT (N2 X N1)->N2 AND P STRICTLY LEFT
!   N1->(N2 X N1).  THUS Q LIES ON THE INTERIOR OF N1->N2.
!   STORE DATA VALUES AND COORDINATES OF N1 AND N2 IN
!   LOCAL VARIABLES
!
   10 W1 = W(N1)
      W2 = W(N2)
      P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
!
! COMPUTE THE CENTRAL PROJECTION OF Q ONTO P2-P1 AND
!   NORMALIZE TO OBTAIN Q
!
      QNORM = 0.
      DO 11 I = 1,3
        Q(I) = B1*P1(I) + B2*P2(I)
   11   QNORM = QNORM + Q(I)*Q(I)
      QNORM = SQRT(QNORM)
      DO 12 I = 1,3
   12   Q(I) = Q(I)/QNORM
!
! STORE OR COMPUTE GRADIENTS AT N1 AND N2
!
      IF (IFLAG .NE. 1) GO TO 14
      DO 13 I = 1,3
        G1(I) = GRAD(I,N1)
   13   G2(I) = GRAD(I,N2)
      GO TO 15
   14 CALL GRADL(NN,N1,X,Y,Z,W,LIST,LPTR,LEND, G1,IERR)
      IF (IERR .LT. 0) GO TO 17
      CALL GRADL(NN,N2,X,Y,Z,W,LIST,LPTR,LEND, G2,IERR)
      IF (IERR .LT. 0) GO TO 17
!
! COMPUTE AN INTERPOLATED VALUE AND NORMAL GRADIENT-
!   COMPONENT AT Q
!
   15 CALL ARCINT(Q,P1,P2,W1,W2,G1,G2, WQ,DUM,GQN)
!
! EXTRAPOLATE TO P -- THE NORMAL GRADIENT COMPONENT GQN IS
!   THE NEGATIVE OF THE COMPONENT IN THE DIRECTION Q->P.
!
      PW = WQ - GQN*ARCLEN(Q,P)
      IER = 1
      RETURN
!
! N, IFLAG, OR IST OUT OF RANGE
!
   16 IER = -1
      RETURN
!
! COLLINEAR NODES
!
   17 IER = -2
      RETURN
!
! THE DISTANCE BETWEEN P AND THE CLOSEST BOUNDARY POINT
!   IS AT LEAST 90 DEGREES
!
   18 IER = -3
      RETURN
      END
      SUBROUTINE SHIFTD (NFRST,NLAST,KK, IARR )
      INTEGER( kind = 4) NFRST, NLAST, KK, IARR(1)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE SHIFTS A SET OF CONTIGUOUS ELEMENTS OF AN
! INTEGER ARRAY KK POSITIONS DOWNWARD (UPWARD IF KK .LT. 0).
! THE LOOPS ARE UNROLLED IN ORDER TO INCREASE EFFICIENCY.
!
! INPUT PARAMETERS - NFRST,NLAST - BOUNDS ON THE PORTION OF
!                                  IARR TO BE SHIFTED.  ALL
!                                  ELEMENTS BETWEEN AND
!                                  INCLUDING THE BOUNDS ARE
!                                  SHIFTED UNLESS NFRST .GT.
!                                  NLAST, IN WHICH CASE NO
!                                  SHIFT OCCURS.
!
!                             KK - NUMBER OF POSITIONS EACH
!                                  ELEMENT IS TO BE SHIFTED.
!                                  IF KK .LT. 0 SHIFT UP.
!                                  IF KK .GT. 0 SHIFT DOWN.
!
!                           IARR - INTEGER ARRAY OF LENGTH
!                                  .GE. NLAST + MAX(KK,0).
!
! NFRST, NLAST, AND KK ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER -        IARR - SHIFTED ARRAY.
!
! MODULES REFERENCED BY SHIFTD - NONE
!
!***********************************************************
!
      INTEGER( kind = 4) INC, K, NF, NL, NLP1, NS, NSL, I, IBAK, INDX, &
              IMAX
      DATA    INC/5/
!
! LOCAL PARAMETERS -
!
! INC =  DO-LOOP INCREMENT (UNROLLING FACTOR) -- IF INC IS
!          CHANGED, STATEMENTS MUST BE ADDED TO OR DELETED
!          FROM THE DO-LOOPS
! K =    LOCAL COPY OF KK
! NF =   LOCAL COPY OF NFRST
! NL =   LOCAL COPY OF NLAST
! NLP1 = NL + 1
! NS =   NUMBER OF SHIFTS
! NSL =  NUMBER OF SHIFTS DONE IN UNROLLED DO-LOOP (MULTIPLE
!          OF INC)
! I =    DO-LOOP INDEX AND INDEX FOR IARR
! IBAK = INDEX FOR DOWNWARD SHIFT OF IARR
! INDX = INDEX FOR IARR
! IMAX = BOUND ON DO-LOOP INDEX
!
      K = KK
      NF = NFRST
      NL = NLAST
      IF (NF .GT. NL  .OR.  K .EQ. 0) RETURN
      NLP1 = NL + 1
      NS = NLP1 - NF
      NSL = INC*(NS/INC)
      IF ( K .LT. 0) GO TO 4
!
! SHIFT DOWNWARD STARTING FROM THE BOTTOM
!
      IF (NSL .LE. 0) GO TO 2
      DO 1 I = 1,NSL,INC
        IBAK = NLP1 - I
        INDX = IBAK + K
        IARR(INDX) = IARR(IBAK)
        IARR(INDX-1) = IARR(IBAK-1)
        IARR(INDX-2) = IARR(IBAK-2)
        IARR(INDX-3) = IARR(IBAK-3)
        IARR(INDX-4) = IARR(IBAK-4)
    1   CONTINUE
!
! PERFORM THE REMAINING NS-NSL SHIFTS ONE AT A TIME
!
    2 IBAK = NLP1 - NSL
    3 IF (IBAK .LE. NF) RETURN
      IBAK = IBAK - 1
      INDX = IBAK + K
      IARR(INDX) = IARR(IBAK)
      GO TO 3
!
! SHIFT UPWARD STARTING FROM THE TOP
!
    4 IF (NSL .LE. 0) GO TO 6
      IMAX = NLP1 - INC
      DO 5 I = NF,IMAX,INC
        INDX = I + K
        IARR(INDX) = IARR(I)
        IARR(INDX+1) = IARR(I+1)
        IARR(INDX+2) = IARR(I+2)
        IARR(INDX+3) = IARR(I+3)
        IARR(INDX+4) = IARR(I+4)
    5   CONTINUE
!
! PERFORM THE REMAINING NS-NSL SHIFTS ONE AT A TIME
!
    6 I = NSL + NF
    7 IF (I .GT. NL) RETURN
      INDX = I + K
      IARR(INDX) = IARR(I)
      I = I + 1
      GO TO 7
      END
      SUBROUTINE GRADL (N,K,X,Y,Z,W,LIST,LPTR,LEND,G,IER)
      INTEGER( kind = 4) N, K,LIST(6*(N-2)),LPTR(6*(N-2)),LEND(N), IER
      REAL( kind = 8 )    X(N), Y(N), Z(N), W(N), G(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF NODES ON THE UNIT
! SPHERE WITH THEIR ASSOCIATED DATA VALUES W, THIS ROUTINE
! ESTIMATES A GRADIENT VECTOR AT NODE K AS FOLLOWS -- THE
! COORDINATE SYSTEM IS ROTATED SO THAT K BECOMES THE NORTH
! POLE, NODE K AND A SET OF NEARBY NODES ARE PROJECTED
! ORTHOGONALLY ONTO THE X-Y PLANE (IN THE NEW COORDINATE
! SYSTEM), A QUADRATIC IS FITTED IN A WEIGHTED LEAST-SQUARES
! SENSE TO THE DATA VALUES AT THE PROJECTED NODES SUCH THAT
! THE VALUE (ASSOCIATED WITH K) AT (0,0) IS INTERPOLATED, X-
! AND Y-PARTIAL DERIVATIVE ESTIMATES DX AND DY ARE COMPUTED
! BY DIFFERENTIATING THE QUADRATIC AT (0,0), AND THE ESTI-
! MATED GRADIENT G IS OBTAINED BY ROTATING (DX,DY,0) BACK TO
! THE ORIGINAL COORDINATE SYSTEM.  NOTE THAT G LIES IN THE
! PLANE TANGENT TO THE SPHERE AT NODE K, I.E. G IS ORTHOGO-
! NAL TO THE UNIT VECTOR REPRESENTED BY NODE K.  A MARQUARDT
! STABILIZATION FACTOR IS USED IF NECESSARY TO ENSURE A
! WELL-CONDITIONED LEAST SQUARES SYSTEM, AND A UNIQUE SOLU-
! TION EXISTS UNLESS THE NODES ARE COLLINEAR.
!
! INPUT PARAMETERS -    N - NUMBER OF NODES IN THE TRIANGU-
!                           LATION.  N .GE. 7.
!
!                       K - NODE AT WHICH THE GRADIENT IS
!                           SOUGHT.  1 .LE. K .LE. N.
!
!                   X,Y,Z - CARTESIAN COORDINATES OF THE
!                           NODES.
!
!                       W - DATA VALUES AT THE NODES.  W(I)
!                           IS ASSOCIATED WITH (X(I),Y(I),
!                           Z(I)) FOR I = 1,...,N.
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to STRIPACK
!                        Subroutine TRMESH.
!
! LIST,LPTR,LEND MAY BE CREATED BY SUBROUTINE TRMESH.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -    G - X-, Y-, AND Z-COMPONENTS (IN
!                            THAT ORDER) OF THE ESTIMATED
!                            GRADIENT AT NODE K UNLESS
!                            IER .LT. 0.
!
!                      IER - ERROR INDICATOR
!                            IER .GE. 6 IF NO ERRORS WERE
!                                       ENCOUNTERED.  IER
!                                       CONTAINS THE NUMBER
!                                       OF NODES (INCLUDING
!                                       K) USED IN THE LEAST
!                                       SQUARES FIT.
!                            IER = -1 IF N OR K IS OUT OF
!                                     RANGE.
!                            IER = -2 IF THE LEAST SQUARES
!                                     SYSTEM HAS NO UNIQUE
!                                     SOLUTION DUE TO DUP-
!                                     LICATE OR COLLINEAR
!                                     NODES.
!
! MODULES REFERENCED BY GRADL - GETNP, CONSTR, APLYR,
!                               SETUP, GIVENS, ROTATE,
!                               APLYRT
!
! INTRINSIC FUNCTIONS CALLED BY GRADL - MIN0, FLOAT, SQRT,
!                                       AMIN1, ABS
!
!***********************************************************
!
      INTEGER( kind = 4) NN, KK, LMN, LMX, LMIN, LMAX, LM1, LNP, &
              NPTS(30), IERR, NP, I, J, IM1, IP1, JP1, L
      REAL( kind = 8 )    WK, SUM, DF, RF, RTOL, AVSQ, AV, RIN, CX, SX, &
              CY, SY, XP, YP, ZP, WT, A(6,6), C, S, DMIN, &
              DTOL, SF, DX, DY
      DATA    LMN/10/
      DATA    LMX/30/, RTOL/1.E-6/, DTOL/.01/, SF/1./
!
! LOCAL PARAMETERS -
!
! NN,KK =     LOCAL COPIES OF N AND K
! LMN,LMX =   MINIMUM AND MAXIMUM VALUES OF LNP FOR N
!               SUFFICIENTLY LARGE.  IN MOST CASES LMN-1
!               NODES ARE USED IN THE FIT.  7 .LE. LMN .LE.
!               LMX.
! LMIN,LMAX = MIN(LMN,N), MIN(LMX,N)
! LM1 =       LMIN-1
! LNP =       LENGTH OF NPTS OR LMAX+1
! NPTS =      ARRAY CONTAINING THE INDICES OF A SEQUENCE OF
!               NODES ORDERED BY ANGULAR DISTANCE FROM K.
!               NPTS(1)=K AND THE FIRST LNP-1 ELEMENTS OF
!               NPTS ARE USED IN THE LEAST SQUARES FIT.
!               UNLESS LNP = LMAX+1, NPTS(LNP) DETERMINES R
!               (SEE RIN).
! IERR =      ERROR FLAG FOR CALLS TO GETNP (NOT CHECKED)
! NP =        ELEMENT OF NPTS TO BE ADDED TO THE SYSTEM
! I,J =       LOOP INDICES
! IM1,IP1 =   I-1, I+1
! JP1 =       J+1
! L =         NUMBER OF COLUMNS OF A**T TO WHICH A ROTATION
!               IS APPLIED
! WK =        W(K) -- DATA VALUE AT NODE K
! SUM =       SUM OF SQUARED EUCLIDEAN DISTANCES (IN THE
!               ROTATED COORDINATE SYSTEM) BETWEEN THE
!               ORIGIN AND THE NODES USED IN THE LEAST
!               SQUARES FIT
! DF =        NEGATIVE Z-COMPONENT (IN THE ROTATED COORDI-
!               NATE SYSTEM) OF AN ELEMENT NP OF NPTS --
!               INCREASING FUNCTION OF THE ANGULAR DISTANCE
!               BETWEEN K AND NP.  DF LIES IN THE INTERVAL
!               (-1,1).
! RF =        VALUE OF DF ASSOCIATED WITH NPTS(LNP) UNLESS
!               LNP = LMAX+1 (SEE RIN)
! RTOL =      TOLERANCE FOR DETERMINING LNP (AND HENCE R) --
!               IF THE INCREASE IN DF BETWEEN TWO SUCCESSIVE
!               ELEMENTS OF NPTS IS LESS THAN RTOL, THEY ARE
!               TREATED AS BEING THE SAME DISTANCE FROM NODE
!               K AND AN ADDITIONAL NODE IS ADDED
! AVSQ =      AV*AV -- ACCUMULATED IN SUM
! AV =        ROOT-MEAN-SQUARE DISTANCE (IN THE ROTATED
!               COORDINATE SYSTEM) BETWEEN THE ORIGIN AND
!               THE NODES (OTHER THAN K) IN THE LEAST
!               SQUARES FIT.  THE FIRST 3 COLUMNS OF A**T
!               ARE SCALED BY 1/AVSQ, THE NEXT 2 BY 1/AV.
! RIN =       INVERSE OF A RADIUS OF INFLUENCE R WHICH
!               ENTERS INTO WT -- R = 1+RF UNLESS ALL ELE-
!               MENTS OF NPTS ARE USED IN THE FIT (LNP =
!               LMAX+1), IN WHICH CASE R IS THE DISTANCE
!               FUNCTION ASSOCIATED WITH SOME POINT MORE
!               DISTANT FROM K THAN NPTS(LMAX)
! CX,SX =     COMPONENTS OF A PLANE ROTATION ABOUT THE X-
!               AXIS WHICH, TOGETHER WITH CY AND SY, DEFINE
!               A MAPPING FROM NODE K TO THE NORTH POLE
!               (0,0,1)
! CY,SY =     COMPONENTS OF A PLANE ROTATION ABOUT THE Y-
!               AXIS
! XP,YP,ZP =  COORDINATES OF NP IN THE ROTATED COORDINATE
!               SYSTEM UNLESS ZP .LT. 0, IN WHICH CASE
!               (XP,YP,0) LIES ON THE EQUATOR
! WT =        WEIGHT FOR THE EQUATION CORRESPONDING TO NP --
!               WT = (R-D)/(R*D) = 1/D - RIN WHERE D = 1-ZP
!               IS ASSOCIATED WITH NP
! A =         TRANSPOSE OF THE (UPPER TRIANGLE OF THE) AUG-
!               MENTED REGRESSION MATRIX
! C,S =       COMPONENTS OF THE PLANE ROTATION USED TO
!               TRIANGULARIZE THE REGRESSION MATRIX
! DMIN =      MINIMUM OF THE MAGNITUDES OF THE DIAGONAL
!               ELEMENTS OF THE TRIANGULARIZED REGRESSION
!               MATRIX
! DTOL =      TOLERANCE FOR DETECTING AN ILL-CONDITIONED
!               SYSTEM -- DMIN IS REQUIRED TO BE AT LEAST
!               DTOL
! SF =        MARQUARDT STABILIZATION FACTOR USED TO DAMP
!               OUT THE FIRST 3 SOLUTION COMPONENTS (SECOND
!               PARTIALS OF THE QUADRATIC) WHEN THE SYSTEM
!               IS ILL-CONDITIONED.  INCREASING SF RESULTS
!               IN MORE DAMPING (A MORE NEARLY LINEAR FIT).
! DX,DY =     X AND Y COMPONENTS OF THE ESTIMATED GRADIENT
!               IN THE ROTATED COORDINATE SYSTEM
!
      NN = N
      KK = K
      WK = W(KK)
!
! CHECK FOR ERRORS AND INITIALIZE LMIN, LMAX
!
      IF (NN .LT. 7  .OR.  KK .LT. 1  .OR.  KK .GT. NN) &
         GO TO 13
      LMIN = MIN0(LMN,NN)
      LMAX = MIN0(LMX,NN)
!
! COMPUTE NPTS, LNP, AVSQ, AV, AND R.
!   SET NPTS TO THE CLOSEST LMIN-1 NODES TO K.  DF CONTAINS
!   THE NEGATIVE Z-COMPONENT (IN THE ROTATED COORDINATE
!   SYSTEM) OF THE NEW NODE ON RETURN FROM GETNP.
!
      SUM = 0.
      NPTS(1) = KK
      LM1 = LMIN - 1
      DO 1 LNP = 2,LM1
        CALL GETNP (X,Y,Z,LIST,LPTR,LEND,LNP, NPTS, DF,IERR)
    1   SUM = SUM + 1. - DF*DF
!
!   ADD ADDITIONAL NODES TO NPTS UNTIL THE INCREASE IN
!     R = 1+RF IS AT LEAST RTOL.
!
      DO 2 LNP = LMIN,LMAX
        CALL GETNP (X,Y,Z,LIST,LPTR,LEND,LNP, NPTS, RF,IERR)
        IF (RF-DF .GE. RTOL) GO TO 3
    2   SUM = SUM + 1. - RF*RF
!
!   USE ALL LMAX NODES IN THE LEAST SQUARES FIT.  R IS
!     ARBITRARILY INCREASED BY 5 PER CENT.
!
      RF = 1.05*RF + .05
      LNP = LMAX + 1
!
!   THERE ARE LNP-2 EQUATIONS CORRESPONDING TO NODES
!     NPTS(2),...,NPTS(LNP-1).
!
    3 AVSQ = SUM/FLOAT(LNP-2)
      AV = SQRT(AVSQ)
      RIN = 1./(1.+RF)
!
! CONSTRUCT THE ROTATION
!
      CALL CONSTR (X(KK),Y(KK),Z(KK), CX,SX,CY,SY)
!
! SET UP THE FIRST 5 EQUATIONS OF THE AUGMENTED REGRESSION
!   MATRIX (TRANSPOSED) AS THE COLUMNS OF A, AND ZERO OUT
!   THE LOWER TRIANGLE (UPPER TRIANGLE OF A) WITH GIVENS
!   ROTATIONS
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
    4     CALL ROTATE (L,C,S, A(JP1,J),A(JP1,I) )
    5   CONTINUE
!
! ADD THE ADDITIONAL EQUATIONS TO THE SYSTEM USING
!   THE LAST COLUMN OF A -- I .LE. LNP.
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
    7     CALL ROTATE (L,C,S, A(JP1,J),A(JP1,6) )
        I = I + 1
        GO TO 6
!
! TEST THE SYSTEM FOR ILL-CONDITIONING
!
    8 DMIN = AMIN1( ABS(A(1,1)),ABS(A(2,2)),ABS(A(3,3)), &
                    ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN .GE. DTOL) GO TO 12
      IF (LNP .GT. LMAX) GO TO 9
!
! ADD ANOTHER NODE TO THE SYSTEM AND INCREASE R --
!   I .EQ. LNP
!
      LNP = LNP + 1
      IF (LNP .LE. LMAX) CALL GETNP (X,Y,Z,LIST,LPTR,LEND,LNP, &
                                     NPTS, RF,IERR)
      RIN = 1./(1.05*(1.+RF))
      GO TO 6
!
! STABILIZE THE SYSTEM BY DAMPING SECOND PARTIALS --ADD
!   MULTIPLES OF THE FIRST THREE UNIT VECTORS TO THE FIRST
!   THREE EQUATIONS.
!
    9 DO 11 I = 1,3
        A(I,6) = SF
        IP1 = I + 1
        DO 10 J = IP1,6
   10     A(J,6) = 0.
        DO 11 J = I,5
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS ( A(J,J),A(J,6), C,S)
   11     CALL ROTATE (L,C,S, A(JP1,J),A(JP1,6) )
!
! TEST THE LINEAR PORTION OF THE STABILIZED SYSTEM FOR
!   ILL-CONDITIONING
!
      DMIN = AMIN1( ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN .LT. DTOL) GO TO 14
!
! SOLVE THE 2 BY 2 TRIANGULAR SYSTEM FOR THE ESTIMATED
!   PARTIAL DERIVATIVES
!
   12 DY = A(6,5)/A(5,5)
      DX = (A(6,4) - A(5,4)*DY)/A(4,4)/AV
      DY = DY/AV
!
! ROTATE THE GRADIENT (DX,DY,0) BACK INTO THE ORIGINAL
!   COORDINATE SYSTEM
!
      CALL APLYRT (DX,DY,CX,SX,CY,SY, G)
      IER = LNP - 1
      RETURN
!
! N OR K IS OUT OF RANGE
!
   13 IER = -1
      RETURN
!
! NO UNIQUE SOLUTION DUE TO COLLINEAR NODES
!
   14 IER = -2
      RETURN
      END

SUBROUTINE WVAL (B1,B2,B3,V1,V2,V3,W1,W2,W3,G1,G2,G3, &
                       IFLAG, PW,PG)
      INTEGER( kind = 4) IFLAG
      REAL( kind = 8 )    B1, B2, B3, V1(3), V2(3), V3(3), W1, W2, W3, &
              G1(3), G2(3), G3(3), PW, PG(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN DATA VALUES AND GRADIENTS AT THE THREE VERTICES OF
! A SPHERICAL TRIANGLE CONTAINING A POINT P, THIS ROUTINE
! COMPUTES THE VALUE AND, OPTIONALLY, THE GRADIENT OF F AT P
! WHERE F INTERPOLATES THE VERTEX DATA.  ALONG THE TRIANGLE
! EDGES, THE INTERPOLATORY FUNCTION F IS THE HERMITE CUBIC
! (WITH RESPECT TO ARC-LENGTH) INTERPOLANT OF THE VALUES AND
! TANGENTIAL GRADIENT COMPONENTS AT THE ENDPOINTS, AND THE
! DERIVATIVE NORMAL TO THE ARC VARIES LINEARLY WITH RESPECT
! TO ARC-LENGTH BETWEEN THE NORMAL GRADIENT COMPONENTS AT
! THE ENDPOINTS.  THUS THE METHOD YIELDS C-1 CONTINUITY WHEN
! USED TO INTERPOLATE OVER A TRIANGULATION.  THE INTERPOLANT
! USES A FIRST-ORDER BLENDING METHOD ON THE UNDERLYING
! PLANAR TRIANGLE.
!
! INPUT PARAMETERS - B1,B2,B3 - BARYCENTRIC COORDINATES OF
!                               PP WITH RESPECT TO THE
!                               (PLANAR) UNDERLYING TRIANGLE
!                               (V1,V2,V3) WHERE PP IS THE
!                               CENTRAL PROJECTION OF P ONTO
!                               THIS TRIANGLE.
!
!                    V1,V2,V3 - CARTESIAN COORDINATES OF THE
!                               VERTICES OF A SPHERICAL TRI-
!                               ANGLE CONTAINING P.  V3 LEFT
!                               V1->V2.
!
!                    W1,W2,W3 - DATA VALUES ASSOCIATED WITH
!                               THE VERTICES.
!
!                    G1,G2,G3 - GRADIENTS ASSOCIATED WITH
!                               THE VERTICES.  GI IS ORTHOG-
!                               ONAL TO VI FOR I = 1,2,3.
!
!                       IFLAG - OPTION INDICATOR
!                               IFLAG = 0 IF ONLY PW IS TO
!                                         BE COMPUTED.
!                               IFLAG = 1 IF BOTH PW AND PG
!                                         ARE DESIRED.
!
!                          PG - VECTOR OF LENGTH 3 IF IFLAG
!                               = 1, NOT USED OTHERWISE.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -      PW - INTERPOLATED VALUE AT P.
!
!                          PG - INTERPOLATED GRADIENT AT P
!                               (ORTHOGONAL TO P) IF IFLAG
!                               = 1.
!
! EACH VECTOR V ABOVE CONTAINS X-, Y-, AND Z-COMPONENTS IN
!   V(1), V(2), AND V(3), RESPECTIVELY.
!
! MODULES REFERENCED BY WVAL - ARCINT, ARCLEN
!
! INTRINSIC FUNCTION CALLED BY WVAL - SQRT
!
!***********************************************************
!
      INTEGER( kind = 4) I
      REAL( kind = 8 )    C1,C2,C3, SUMM, U1(3), U2(3), U3(3), U1N, &
              U2N, U3N, Q1(3), Q2(3), Q3(3), VAL, W, G(3),&
              DUM
!
! LOCAL PARAMETERS -
!
! I =           DO-LOOP INDEX
! C1,C2,C3 =    COEFFICIENTS (WEIGHT FUNCTIONS) OF PARTIAL
!                 INTERPOLANTS.  C1 = 1 ON THE EDGE OPPOSITE
!                 V1 AND C1 = 0 ON THE OTHER EDGES.  SIMI-
!                 LARLY FOR C2 AND C3.  C1+C2+C3 = 1.
! SUMM =         QUANTITY USED TO NORMALIZE C1, C2, AND C3
! U1,U2,U3 =    POINTS ON THE BOUNDARY OF THE PLANAR TRIAN-
!                 GLE AND LYING ON THE LINES CONTAINING PP
!                 AND THE VERTICES.  U1 IS OPPOSITE V1, ETC.
! U1N,U2N,U3N = QUANTITIES USED TO COMPUTE Q1, Q2, AND Q3
!                 (MAGNITUDES OF U1, U2, AND U3)
! Q1,Q2,Q3 =    CENTRAL PROJECTIONS OF U1, U2, AND U3 ONTO
!                 THE SPHERE AND THUS LYING ON AN ARC OF THE
!                 SPHERICAL TRIANGLE
! VAL =         LOCAL VARIABLE USED TO ACCUMULATE THE CON-
!                 TRIBUTIONS TO PW
! W,G =         VALUE AND GRADIENT AT Q1, Q2, OR Q3 OBTAINED
!                 BY INTERPOLATION ALONG ONE OF THE ARCS OF
!                 THE SPHERICAL TRIANGLE
! DUM =         DUMMY VARIABLE FOR THE CALL TO ARCINT
!
!
! COMPUTE WEIGHT FUNCTIONS C1, C2, AND C3
!
      C1 = B2*B3
      C2 = B3*B1
      C3 = B1*B2
      SUMM = C1 + C2 + C3
      IF (SUMM .GT. 0.) GO TO 2
!
! P COINCIDES WITH A VERTEX
!
      PW = B1*W1 + B2*W2 + B3*W3
      IF (IFLAG .NE. 1) RETURN
      DO 1 I = 1,3
    1   PG(I) = B1*G1(I) + B2*G2(I) + B3*G3(I)
      RETURN
!
! NORMALIZE C1, C2, AND C3
!
    2 C1 = C1/SUMM
      C2 = C2/SUMM
      C3 = C3/SUMM
!
! COMPUTE (U1,U2,U3) AND (U1N,U2N,U3N)
!
      U1N = 0.
      U2N = 0.
      U3N = 0.
      DO 3 I = 1,3
        U1(I) = (B2*V2(I) + B3*V3(I))/(B2+B3)
        U2(I) = (B3*V3(I) + B1*V1(I))/(B3+B1)
        U3(I) = (B1*V1(I) + B2*V2(I))/(B1+B2)
        U1N = U1N + U1(I)*U1(I)
        U2N = U2N + U2(I)*U2(I)
    3   U3N = U3N + U3(I)*U3(I)
!
! COMPUTE Q1, Q2, AND Q3
!
      U1N = SQRT(U1N)
      U2N = SQRT(U2N)
      U3N = SQRT(U3N)
      DO 4 I = 1,3
        Q1(I) = U1(I)/U1N
        Q2(I) = U2(I)/U2N
    4   Q3(I) = U3(I)/U3N
!
! COMPUTE INTERPOLATED VALUE (VAL) AT P BY LOOPING ON
!   TRIANGLE SIDES
!
      VAL = 0.
!
! CONTRIBUTION FROM SIDE OPPOSITE V1 --
!
!   COMPUTE VALUE AND GRADIENT AT Q1 BY INTERPOLATING
!     BETWEEN V2 AND V3
!
      CALL ARCINT (Q1,V2,V3,W2,W3,G2,G3, W,G,DUM)
!
!   ADD IN THE CONTRIBUTION
!
      VAL = VAL + C1*( W + B1*B1*(3.-2.*B1)*(W1-W) + &
                  B1*(1.-B1)* &
                  (B1*(G1(1)*U1(1)+G1(2)*U1(2)+G1(3)*U1(3)) &
                   + (1.-B1)* &
                   (G(1)*V1(1)+G(2)*V1(2)+G(3)*V1(3))/U1N) )
!
! CONTRIBUTION FROM SIDE OPPOSITE V2 --
!
!   COMPUTE VALUE AND GRADIENT AT Q2 BY INTERPOLATING
!     BETWEEN V3 AND V1
!
      CALL ARCINT (Q2,V3,V1,W3,W1,G3,G1, W,G,DUM)
!
!   ADD IN THE CONTRIBUTION
!
      VAL = VAL + C2*( W + B2*B2*(3.-2.*B2)*(W2-W) + &
                  B2*(1.-B2)* &
                  (B2*(G2(1)*U2(1)+G2(2)*U2(2)+G2(3)*U2(3)) &
                   + (1.-B2)* &
                   (G(1)*V2(1)+G(2)*V2(2)+G(3)*V2(3))/U2N) )
!
! CONTRIBUTION FROM SIDE OPPOSITE V3 --
!
!   COMPUTE INTERPOLATED VALUE AND GRADIENT AT Q3
!     BY INTERPOLATING BETWEEN V1 AND V2
!
      CALL ARCINT (Q3,V1,V2,W1,W2,G1,G2, W,G,DUM)
!
!   ADD IN THE FINAL CONTRIBUTION
!
      VAL = VAL + C3*( W + B3*B3*(3.-2.*B3)*(W3-W) + &
                  B3*(1.-B3)* &
                  (B3*(G3(1)*U3(1)+G3(2)*U3(2)+G3(3)*U3(3)) &
                   + (1.-B3)* &
                   (G(1)*V3(1)+G(2)*V3(2)+G(3)*V3(3))/U3N) )
      PW = VAL
      IF (IFLAG .NE. 1) RETURN
      RETURN
      END

      subroutine interp_n(npts,nptso,order,olats,olons,x,y,z,datain,lst,&
                           lptr,lend,odata,edata,ierr)

      integer( kind = 4 ), intent(in) :: npts, nptso, order
      integer( kind = 4 ), intent(out) :: ierr
      real( kind = 8 ), intent(in), dimension(nptso) :: olats,olons
      real( kind = 8 ), intent(in), dimension(npts) :: datain,x,y,z
      real( kind = 8 ), intent(out), dimension(nptso) :: odata
      integer( kind = 4 ), intent(in), dimension(npts) :: lend
      integer( kind = 4 ), intent(in), dimension(6*(npts-2)) :: lst,lptr
      integer( kind = 4 ) n,ierr1,ist
      integer( kind = 4 ), intent(out), dimension(nptso) :: edata
      ist = 1
      ierr = 0

      if (order .ne. 0 .and. order .ne. 1 .and. order .ne. 3) then
         print *,'fatal error: interp order must be 0, 1 or 3'
         stop
      endif

      do n=1,nptso
         call interp(npts,order,olats(n),olons(n),x,y,z,datain,lst,lptr,&
                     lend,ist,odata(n),ierr1)

         edata(n) = ierr1

         if (ierr1 .lt. 0) then
           !print *,n,'warning: ierr = ',ierr1,' in interp_n'
           !print *,olats(n), olons(n), npts
           !stop
           ierr = ierr + ierr1

         endif
      enddo

      end subroutine interp_n




      subroutine trlist ( n, list, lptr, lend, nrow, nt, ltri, ier )

!*****************************************************************************80
!
!! TRLIST converts a triangulation data structure to a triangle list.
!
!  Discussion:
!
!    This subroutine converts a triangulation data structure
!    from the linked list created by TRMESH to a triangle list.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), linked
!    list data structure defining the triangulation.  Refer to TRMESH.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows (entries per triangle)
!    reserved for the triangle list LTRI.  The value must be 6 if only the
!    vertex indexes and neighboring triangle indexes are to be stored, or 9
!    if arc indexes are also to be assigned and stored.  Refer to LTRI.
!
!    Output, integer ( kind = 4 ) NT, the number of triangles in the
!    triangulation unless IER /=0, in which case NT = 0.  NT = 2N-NB-2 if
!    NB >= 3 or 2N-4 if NB = 0, where NB is the number of boundary nodes.
!
!    Output, integer ( kind = 4 ) LTRI(NROW,*).  The second dimension of LTRI
!    must be at least NT, where NT will be at most 2*N-4.  The J-th column
!    contains the vertex nodal indexes (first three rows), neighboring triangle
!    indexes (second three rows), and, if NROW = 9, arc indexes (last three
!    rows) associated with triangle J for J = 1,...,NT.  The vertices are
!    ordered counterclockwise with the first vertex taken to be the one
!    with smallest index.  Thus, LTRI(2,J) and LTRI(3,J) are larger than
!    LTRI(1,J) and index adjacent neighbors of node LTRI(1,J).  For
!    I = 1,2,3, LTRI(I+3,J) and LTRI(I+6,J) index the triangle and arc,
!    respectively, which are opposite (not shared by) node LTRI(I,J), with
!    LTRI(I+3,J) = 0 if LTRI(I+6,J) indexes a boundary arc.  Vertex indexes
!    range from 1 to N, triangle indexes from 0 to NT, and, if included,
!    arc indexes from 1 to NA, where NA = 3N-NB-3 if NB >= 3 or 3N-6 if
!    NB = 0.  The triangles are ordered on first (smallest) vertex indexes.
!
!    Output, integer ( kind = 4 ) IER, error indicator.
!    0, if no errors were encountered.
!    1, if N or NROW is outside its valid range on input.
!    2, if the triangulation data structure (LIST,LPTR,LEND) is invalid.
!      Note, however, that these arrays are not completely tested for validity.
!
!  Local parameters:
!
!    ARCS =     Logical variable with value TRUE iff are
!               indexes are to be stored
!    I,J =      LTRI row indexes (1 to 3) associated with
!               triangles KT and KN, respectively
!    I1,I2,I3 = Nodal indexes of triangle KN
!    ISV =      Variable used to permute indexes I1,I2,I3
!    KA =       Arc index and number of currently stored arcs
!    KN =       Index of the triangle that shares arc I1-I2 with KT
!    KT =       Triangle index and number of currently stored triangles
!    LP =       LIST pointer
!    LP2 =      Pointer to N2 as a neighbor of N1
!    LPL =      Pointer to the last neighbor of I1
!    LPLN1 =    Pointer to the last neighbor of N1
!    N1,N2,N3 = Nodal indexes of triangle KT
!    NM2 =      N-2
!
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) nrow

    logical              arcs
    integer ( kind = 4 ) i
    integer ( kind = 4 ) i1
    integer ( kind = 4 ) i2
    integer ( kind = 4 ) i3
    integer ( kind = 4 ) ier
    integer ( kind = 4 ) isv
    integer ( kind = 4 ) j
    integer ( kind = 4 ) ka
    integer ( kind = 4 ) kn
    integer ( kind = 4 ) kt
    integer ( kind = 4 ) lend(n)
    integer ( kind = 4 ) list(6*(n-2))
    integer ( kind = 4 ) lp
    integer ( kind = 4 ) lp2
    integer ( kind = 4 ) lpl
    integer ( kind = 4 ) lpln1
    integer ( kind = 4 ) lptr(6*(n-2))
    integer ( kind = 4 ) ltri(nrow,*)
    integer ( kind = 4 ) n1
    integer ( kind = 4 ) n2
    integer ( kind = 4 ) n3
    integer ( kind = 4 ) nm2
    integer ( kind = 4 ) nt
!
!  Test for invalid input parameters.
!
if ( n < 3 .or. ( nrow /= 6 .and. nrow /= 9 ) ) then
  nt = 0
  ier = 1
  return
end if
!
!  Initialize parameters for loop on triangles KT = (N1,N2,
!  N3), where N1 < N2 and N1 < N3.
!
!  ARCS = TRUE iff arc indexes are to be stored.
!  KA,KT = Numbers of currently stored arcs and triangles.
!  NM2 = Upper bound on candidates for N1.
!
arcs = nrow == 9
ka = 0
kt = 0
nm2 = n-2
!
!  Loop on nodes N1.
!
do n1 = 1, nm2
!
!  Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
!  to the last neighbor of N1, and LP2 points to N2.
!
  lpln1 = lend(n1)
  lp2 = lpln1

1   continue

    lp2 = lptr(lp2)
    n2 = list(lp2)
    lp = lptr(lp2)
    n3 = abs ( list(lp) )

    if ( n2 < n1 .or. n3 < n1 ) then
      go to 8
    end if
!
!  Add a new triangle KT = (N1,N2,N3).
!
    kt = kt + 1
    ltri(1,kt) = n1
    ltri(2,kt) = n2
    ltri(3,kt) = n3
!
!  Loop on triangle sides (I2,I1) with neighboring triangles
!  KN = (I1,I2,I3).
!
    do i = 1, 3

      if ( i == 1 ) then
        i1 = n3
        i2 = n2
      else if ( i == 2 ) then
        i1 = n1
        i2 = n3
      else
        i1 = n2
        i2 = n1
      end if
!
!  Set I3 to the neighbor of I1 that follows I2 unless
!  I2->I1 is a boundary arc.
!
      lpl = lend(i1)
      lp = lptr(lpl)

      do

        if ( list(lp) == i2 ) then
          go to 3
        end if

        lp = lptr(lp)

        if ( lp == lpl ) then
          exit
        end if

      end do
!
!  Invalid triangulation data structure:  I1 is a neighbor of
!  I2, but I2 is not a neighbor of I1.
!
      if ( abs ( list(lp) ) /= i2 ) then
        nt = 0
        ier = 2
        return
      end if
!
!  I2 is the last neighbor of I1.  Bypass the search for a neighboring
!  triangle if I2->I1 is a boundary arc.
!
      kn = 0

      if ( list(lp) < 0 ) then
        go to 6
      end if
!
!  I2->I1 is not a boundary arc, and LP points to I2 as
!  a neighbor of I1.
!
3       continue

      lp = lptr(lp)
      i3 = abs ( list(lp) )
!
!  Find J such that LTRI(J,KN) = I3 (not used if KT < KN),
!  and permute the vertex indexes of KN so that I1 is smallest.
!
      if ( i1 < i2 .and. i1 < i3 ) then
        j = 3
      else if ( i2 < i3 ) then
        j = 2
        isv = i1
        i1 = i2
        i2 = i3
        i3 = isv
      else
        j = 1
        isv = i1
        i1 = i3
        i3 = i2
        i2 = isv
      end if
!
!  Test for KT < KN (triangle index not yet assigned).
!
      if ( n1 < i1 ) then
        cycle
      end if
!
!  Find KN, if it exists, by searching the triangle list in
!  reverse order.
!
      do kn = kt-1, 1, -1
        if ( ltri(1,kn) == i1 .and. &
             ltri(2,kn) == i2 .and. &
             ltri(3,kn) == i3 ) then
          go to 5
        end if
      end do

      cycle
!
!  Store KT as a neighbor of KN.
!
5       continue

      ltri(j+3,kn) = kt
!
!  Store KN as a neighbor of KT, and add a new arc KA.
!
6       continue

      ltri(i+3,kt) = kn

      if ( arcs ) then
        ka = ka + 1
        ltri(i+6,kt) = ka
        if ( kn /= 0 ) then
          ltri(j+6,kn) = ka
        end if
      end if

      end do
!
!  Bottom of loop on triangles.
!
8     continue

    if ( lp2 /= lpln1 ) then
      go to 1
    end if

9     continue

end do

nt = kt
ier = 0

return
end




    subroutine trlist2 ( n, list, lptr, lend, nt, ltri, ier )

!*****************************************************************************80
!
!! TRLIST2 converts a triangulation data structure to a triangle list.
!
!  Discussion:
!
!    This subroutine converts a triangulation data structure
!    from the linked list created by TRMESH to a triangle list.
!
!    It is a version of TRLIST for the special case where the triangle
!    list should only include the nodes that define each triangle.
!
!  Modified:
!
!    21 July 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), linked
!    list data structure defining the triangulation.  Refer to TRMESH.
!
!    Output, integer ( kind = 4 ) NT, the number of triangles in the
!    triangulation unless IER /=0, in which case NT = 0.  NT = 2N-NB-2 if
!    NB >= 3 or 2N-4 if NB = 0, where NB is the number of boundary nodes.
!
!    Output, integer ( kind = 4 ) LTRI(3,*).  The second dimension of LTRI
!    must be at least NT, where NT will be at most 2*N-4.  The J-th column
!    contains the vertex nodal indexes associated with triangle J for
!    J = 1,...,NT.  The vertices are ordered counterclockwise with the first
!    vertex taken to be the one with smallest index.  Thus, LTRI(2,J) and
!    LTRI(3,J) are larger than LTRI(1,J) and index adjacent neighbors of node
!    LTRI(1,J).  The triangles are ordered on first (smallest) vertex indexes.
!
!    Output, integer ( kind = 4 ) IER, error indicator.
!    0, if no errors were encountered.
!    1, if N is outside its valid range on input.
!    2, if the triangulation data structure (LIST,LPTR,LEND) is invalid.
!      Note, however, that these arrays are not completely tested for validity.
!
!  Local parameters:
!
!    I,J =      LTRI row indexes (1 to 3) associated with
!               triangles KT and KN, respectively
!    I1,I2,I3 = Nodal indexes of triangle KN
!    ISV =      Variable used to permute indexes I1,I2,I3
!    KA =       Arc index and number of currently stored arcs
!    KN =       Index of the triangle that shares arc I1-I2 with KT
!    KT =       Triangle index and number of currently stored triangles
!    LP =       LIST pointer
!    LP2 =      Pointer to N2 as a neighbor of N1
!    LPL =      Pointer to the last neighbor of I1
!    LPLN1 =    Pointer to the last neighbor of N1
!    N1,N2,N3 = Nodal indexes of triangle KT
!    NM2 =      N-2
!
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) i
    integer ( kind = 4 ) i1
    integer ( kind = 4 ) i2
    integer ( kind = 4 ) i3
    integer ( kind = 4 ) ier
    integer ( kind = 4 ) isv
    integer ( kind = 4 ) j
    integer ( kind = 4 ) ka
    integer ( kind = 4 ) kn
    integer ( kind = 4 ) kt
    integer ( kind = 4 ) lend(n)
    integer ( kind = 4 ) list(6*(n-2))
    integer ( kind = 4 ) lp
    integer ( kind = 4 ) lp2
    integer ( kind = 4 ) lpl
    integer ( kind = 4 ) lpln1
    integer ( kind = 4 ) lptr(6*(n-2))
    integer ( kind = 4 ) ltri(3,2*n-4)
    integer ( kind = 4 ) n1
    integer ( kind = 4 ) n2
    integer ( kind = 4 ) n3
    integer ( kind = 4 ) nm2
    integer ( kind = 4 ) nt
!
!  Test for invalid input parameters.
!
    if ( n < 3 ) then
      nt = 0
      ier = 1
      return
    end if
!
!  Initialize parameters for loop on triangles KT = (N1,N2,
!  N3), where N1 < N2 and N1 < N3.
!
!  KA,KT = Numbers of currently stored arcs and triangles.
!  NM2 = Upper bound on candidates for N1.
!
    ka = 0
    kt = 0
    nm2 = n-2
!
!  Loop on nodes N1.
!
    do n1 = 1, nm2
!
!  Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
!  to the last neighbor of N1, and LP2 points to N2.
!
      lpln1 = lend(n1)
      lp2 = lpln1

1   continue

    lp2 = lptr(lp2)
    n2 = list(lp2)
    lp = lptr(lp2)
    n3 = abs ( list(lp) )

    if ( n2 < n1 .or. n3 < n1 ) then
      go to 8
    end if

!
!  Add a new triangle KT = (N1,N2,N3).
!
    kt = kt + 1
    ltri(1,kt) = n1
    ltri(2,kt) = n2
    ltri(3,kt) = n3
!
!  Loop on triangle sides (I2,I1) with neighboring triangles
!  KN = (I1,I2,I3).
!
    do i = 1, 3

      if ( i == 1 ) then
        i1 = n3
        i2 = n2
      else if ( i == 2 ) then
        i1 = n1
        i2 = n3
      else
        i1 = n2
        i2 = n1
      end if
!
!  Set I3 to the neighbor of I1 that follows I2 unless
!  I2->I1 is a boundary arc.
!
      lpl = lend(i1)
      lp = lptr(lpl)

      do

        if ( list(lp) == i2 ) then
          go to 3
        end if

        lp = lptr(lp)

        if ( lp == lpl ) then
          exit
        end if

      end do
!
!  Invalid triangulation data structure:  I1 is a neighbor of
!  I2, but I2 is not a neighbor of I1.
!
      if ( abs ( list(lp) ) /= i2 ) then
        nt = 0
        ier = 2
        return
      end if
!
!  I2 is the last neighbor of I1.  Bypass the search for a neighboring
!  triangle if I2->I1 is a boundary arc.
!
      kn = 0

      if ( list(lp) < 0 ) then
        go to 6
      end if
!
!  I2->I1 is not a boundary arc, and LP points to I2 as
!  a neighbor of I1.
!
3       continue

      lp = lptr(lp)
      i3 = abs ( list(lp) )
!
!  Find J such that LTRI(J,KN) = I3 (not used if KT < KN),
!  and permute the vertex indexes of KN so that I1 is smallest.
!
      if ( i1 < i2 .and. i1 < i3 ) then
        j = 3
      else if ( i2 < i3 ) then
        j = 2
        isv = i1
        i1 = i2
        i2 = i3
        i3 = isv
      else
        j = 1
        isv = i1
        i1 = i3
        i3 = i2
        i2 = isv
      end if
!
!  Test for KT < KN (triangle index not yet assigned).
!
      if ( n1 < i1 ) then
        cycle
      end if
!
!  Find KN, if it exists, by searching the triangle list in
!  reverse order.
!
      do kn = kt-1, 1, -1
        if ( ltri(1,kn) == i1 .and. &
             ltri(2,kn) == i2 .and. &
             ltri(3,kn) == i3 ) then
          go to 5
        end if
      end do

      cycle

5       continue

6       continue

      end do
!
!  Bottom of loop on triangles.
!
8     continue

    if ( lp2 /= lpln1 ) then
      go to 1
    end if

9     continue

    end do

    nt = kt
    ier = 0

    return
    end
