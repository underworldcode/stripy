! Copyright 2016-2017 Louis Moresi, Ben Mather, Romain Beucher
!
! This file is part of Quagmire.
!
! Quagmire is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or any later version.
!
! Quagmire is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with Quagmire.  If not, see <http://www.gnu.org/licenses/>.
!

subroutine ntriw ( n, x, y, nt, ltri, area, weight )
!*****************************************************************************
!! NTRIW computes the pointwise area to calculate local areas on a mesh
!
! Parameters:
!
!   Input, integer ( kind = 4 ), n
!   number of points in the triangulation
!
!   Input, real ( kind = 8 ), x(n), y(n)
!   x and y coordinates that make up the triangulation
!
!   Input, integer ( kind = 4 ), n
!   number of points in the triangulation
!
!   Input, integer ( kind = 4 ), nt
!   number of triangles in the triangulation
!
!   Input, integer ( kind = 4 ), ltri(nt)
!   list of triangles in the triangulation
!
!   Ouput, real ( kind = 8 ), area(n), weight(n)
!   areas and weights for each point

  implicit none

  integer ( kind = 4 ) n,nt,i
  real ( kind = 8 ) x(n),y(n),area(n)
  integer ( kind = 4 ) weight(n)
  integer ( kind = 4 ) ltri(3,nt)
  real ( kind = 8 ) v1x,v1y,v2x,v2y
  integer ( kind = 4 ) tri(3)

!
! Get 2 sides of triangle
!
  do i = 1, nt
    tri = ltri(:,i)
    v1x = x(tri(2)) - x(tri(1))
    v1y = y(tri(2)) - y(tri(1))
    v2x = x(tri(1)) - x(tri(3))
    v2y = y(tri(1)) - y(tri(3))

    area(tri) = area(tri) + abs(v1x*v2y - v1y*v2x)
    weight(tri) = weight(tri) + 1
  end do
!
! Now we divide each element by 6
!
  area = area/6
  return
end

subroutine lambert_equal_area( lon, lat, n, lons, lats, x, y )
! Project lons and lats relative to maintain a consistent area
! relative to lon and lat
! 
! The Lambert azimuthal equal-area projection is a particular mapping
! from a sphere to a disk. It accurately represents area in all regions
! of the sphere, but it does not accurately represent angles.
! 
! https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection

  implicit none

! INPUT VARIABLES
  integer ( kind = 4 ) n
  real ( kind = 8 ) lon, lat, lons(n), lats(n)
! WORK VARIABLES
  real ( kind = 8 ) cost1, sint1, sint(n), cost(n), k_prime(n)
  real ( kind = 8 ) costcosll0(n)
! OUTPUT VARIABLES
  real ( kind = 8 ) x(n), y(n)

  cost1 = cos(lat)
  sint1 = sin(lat)
  cost = cos(lats)
  sint = sin(lats)
  costcosll0 = cost*cos(lons - lon)
  
  k_prime = sqrt(2.0/(1.0 + sint1*sint + cost1*costcosll0))
  x = k_prime*cost*sin(lons - lon)
  y = k_prime*cost1*sint - sint1*costcosll0
  return
end subroutine

function geocentric_radius( lat, r1, r2 )
! Calculate the radius of an oblate spheroid (like the earth)
!
! Parameters
! ----------
! lat : array of floats
!     latitudinal coordinates in radians
! r1 : float
!     radius at the equator (in metres)
! r2 : float
!     radius at the poles (in metres)
!
! Returns
! -------
! r : array of floats
!     radius at provided latitudes `lat` in metres

  implicit none

  real ( kind = 8 ), intent(in) :: lat, r1, r2
  real ( kind = 8 ) coslat, sinlat, num, den
  real ( kind = 8 ) geocentric_radius

  coslat = cos(lat)
  sinlat = sin(lat)
  num = ((r1**2)*coslat)**2 + ((r2**2)*sinlat)**2
  den = (r1*coslat)**2 + (r2*sinlat)**2
  geocentric_radius = sqrt(num/den)
  return
end function

subroutine ntriw_proj( n, lon, lat, nt, ltri, r1, r2, area, weight)
!*****************************************************************************
!! NTRIW_S computes the pointwise area to calculate local areas
!! on a spherical mesh
!
! Parameters:
!
!   Input, integer ( kind = 4 ), n
!   number of points in the triangulation
!
!   Input, real ( kind = 8 ), x(n), y(n)
!   lon and lat coordinates that make up the triangulation
!
!   Input, integer ( kind = 4 ), n
!   number of points in the triangulation
!
!   Input, integer ( kind = 4 ), nt
!   number of triangles in the triangulation
!
!   Input, integer ( kind = 4 ), ltri(nt)
!   list of triangles in the triangulation
!
!   Ouput, real ( kind = 8 ), area(n), weight(n)
!   areas and weights for each point

  implicit none

  integer ( kind = 4 ) n, nt, i
  real ( kind = 8 ) lon(n)
  real ( kind = 8 ) lat(n)
  real ( kind = 8 ) area(n)
  real ( kind = 8 ) v1x, v1y, v2x, v2y, r, r1, r2
  real ( kind = 8 ) x(3), y(3), ilon(3), ilat(3)
  real ( kind = 8 ) geocentric_radius
  integer ( kind = 4 ) weight(n)
  integer ( kind = 4 ) ltri(3,nt)
  integer ( kind = 4 ) tri(3)

  area(:) = 0
  weight(:) = 0

  do i = 1, nt
    tri = ltri(:,i)
    ilon = lon(tri)
    ilat = lat(tri)

    call lambert_equal_area(ilon(1), ilat(1), 3, ilon, ilat, x, y)
    r = geocentric_radius(ilat(1), r1, r2)
    x = x * r
    y = y * r

    v1x = x(2) - x(1)
    v1y = y(2) - y(1)
    v2x = x(1) - x(3)
    v2y = y(1) - y(3)

    area(tri) = area(tri) + abs(v1x*v2y - v1y*v2x)
    weight(tri) = weight(tri) + 1
  end do

  area = area / 6
  return
end subroutine

subroutine ntriw_s( n, nt, ltri, tri_area, area, weight)
!*****************************************************************************
!! NTRIW_S computes the pointwise area to calculate local areas
!! on a spherical mesh
!
! Parameters:
!
!   Input, integer ( kind = 4 ), n
!   number of points in the triangulation
!
!   Input, real ( kind = 8 ), x(n), y(n)
!   lon and lat coordinates that make up the triangulation
!
!   Input, integer ( kind = 4 ), n
!   number of points in the triangulation
!
!   Input, integer ( kind = 4 ), nt
!   number of triangles in the triangulation
!
!   Input, integer ( kind = 4 ), ltri(nt)
!   list of triangles in the triangulation
!
!   Ouput, real ( kind = 8 ), area(n), weight(n)
!   areas and weights for each point

  implicit none

  integer ( kind = 4 ) n, nt, i
  real ( kind = 8 ) area(n)
  real ( kind = 8 ) tri_area(nt)
  integer ( kind = 4 ) weight(n)
  integer ( kind = 4 ) ltri(3,nt)
  integer ( kind = 4 ) tri(3)

  area(:) = 0
  weight(:) = 0

  do i = 1, nt
    tri = ltri(:,i)
    area(tri) = area(tri) + tri_area(i)
    weight(tri) = weight(tri) + 1
  end do
  area = area / 3
  return
end subroutine

subroutine fill_mask_to_idx( rows, cols, mask, idx )

  implicit none

  integer ( kind = 4 ) rows, cols, i
  integer ( kind = 4 ) idx(rows)
  logical ( kind = 4 ) mask(rows,cols)

  mask(:,:) = .false.

  do i = 1, rows
    mask(i,1:idx(i)) = .true.
  end do
  return
end subroutine

subroutine add_pt ( pt, array, n )
!*****************************************************************************
! ADD_PT adds a point to an integer array if it does not already exist
! in this way it mimics a set in Python

  implicit none

  integer ( kind = 4 ) pt, n, i
  integer ( kind = 4 ) array(n)

  do i = 1, n
    if (array(i) .eq. 0) then
      exit
    else if (array(i) .eq. pt) then
      return
    end if
  end do

  array(i) = pt
  return
end subroutine

subroutine ncloud ( nt, ltri, n, nnz, ecloud )
!*****************************************************************************
! NCLOUD finds all neighbours and extended neighbours for every point
! in a triangulation

  implicit none

  integer ( kind = 4 ) nt, nnz, n
  integer ( kind = 4 ) ltri(3,nt)
  integer ( kind = 4 ) ecloud(n,nnz*nnz/2)
  integer ( kind = 4 ) cloud(n,nnz)
  integer ( kind = 4 ) tri(3)
  integer ( kind = 4 ) neighbours(nnz), eneighbours(nnz)
  integer ( kind = 4 ) i, t, ncol, np, enp, pt, ept

  ncol = nnz*nnz/2

  ecloud(:,:) = 0
  do t = 1, nt
    tri = ltri(:,t)

    call add_pt(tri(2), ecloud(tri(1),:), ncol)
    call add_pt(tri(3), ecloud(tri(1),:), ncol)
    call add_pt(tri(1), ecloud(tri(2),:), ncol)
    call add_pt(tri(3), ecloud(tri(2),:), ncol)
    call add_pt(tri(1), ecloud(tri(3),:), ncol)
    call add_pt(tri(2), ecloud(tri(3),:), ncol)

  end do

  cloud = ecloud(:,1:nnz)

  do i = 1, n
    neighbours = cloud(i,:)
    do np = 1, nnz
      ! Get neighbours
      pt = neighbours(np)
      if (pt .eq. 0) then
        exit
      endif
      eneighbours = cloud(pt,:)
      do enp = 1, nnz
        ! Get extended neighbours
        ept = eneighbours(enp)
        if (ept .eq. 0) then
          exit
        endif
        call add_pt(ept, ecloud(i,:), ncol)
      end do
    end do
  end do

  return
end subroutine
