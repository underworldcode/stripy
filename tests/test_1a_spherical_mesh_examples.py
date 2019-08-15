import pytest
import stripy

# These tests just test the mesh construction doesn't fail at the moment
# ToDo: add meaningful assertions

# It's hard to parameterise these as the calls are a little
# bit different in each case

def test_kd_tree_building():
    mesh = stripy.spherical_meshes.icosahedral_mesh(include_face_points=True, refinement_levels=2, tree=True)


def test_triangulated_cube_mesh():

    mesh = stripy.spherical_meshes.triangulated_cube_mesh(refinement_levels=0)
    mesh = stripy.spherical_meshes.triangulated_cube_mesh(refinement_levels=2)


def test_octahedral_mesh():

    mesh = stripy.spherical_meshes.octahedral_mesh(include_face_points=False, refinement_levels=0)
    mesh = stripy.spherical_meshes.octahedral_mesh(include_face_points=False, refinement_levels=2)
    mesh = stripy.spherical_meshes.octahedral_mesh(include_face_points=True,  refinement_levels=0)
    mesh = stripy.spherical_meshes.octahedral_mesh(include_face_points=True,  refinement_levels=2)

def test_icosahedral_mesh():

    mesh = stripy.spherical_meshes.icosahedral_mesh(include_face_points=False, refinement_levels=0)
    mesh = stripy.spherical_meshes.icosahedral_mesh(include_face_points=False, refinement_levels=2)
    mesh = stripy.spherical_meshes.icosahedral_mesh(include_face_points=True,  refinement_levels=0)
    mesh = stripy.spherical_meshes.icosahedral_mesh(include_face_points=True,  refinement_levels=2)


def test_triangulated_soccerball_mesh():

    mesh = stripy.spherical_meshes.triangulated_soccerball_mesh(refinement_levels=0)
    mesh = stripy.spherical_meshes.triangulated_soccerball_mesh(refinement_levels=2)

def test_random_mesh():

    mesh = stripy.spherical_meshes.random_mesh(number_of_points=10)
    mesh = stripy.spherical_meshes.random_mesh(number_of_points=10000)


def test_uniform_ring_mesh():

    mesh = stripy.spherical_meshes.uniform_ring_mesh(resolution=9, refinement_levels=0)
    mesh = stripy.spherical_meshes.uniform_ring_mesh(resolution=3, refinement_levels=3)


## Cartesian meshes

def test_cart_square_border():
    mesh = stripy.cartesian_meshes.square_border((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=0)
    mesh = stripy.cartesian_meshes.square_border((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=2)
    mesh = stripy.cartesian_meshes.square_border((-1.0,1.0,-1.0,1.0), 0.1, 0.2, refinement_levels=0)


def test_cart_square_mesh():
    mesh = stripy.cartesian_meshes.square_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=0)
    mesh = stripy.cartesian_meshes.square_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=2)
    mesh = stripy.cartesian_meshes.square_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.2, refinement_levels=0)


def test_elliptical_mesh():
    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=0)
    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=2)
    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.2, refinement_levels=0)
    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.1, random_scale=0.01, refinement_levels=0)
    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.2, random_scale=0.01, refinement_levels=0)


def test_cart_random_mesh():
    mesh = stripy.cartesian_meshes.random_mesh((-1.0,1.0,-1.0,1.0), 5000)
    mesh = stripy.cartesian_meshes.random_mesh((-1.0,1.0,-1.0,1.0), 10000)
    mesh = stripy.cartesian_meshes.random_mesh((-1.0,1.0,-1.0,1.0), 20000)
