import pytest
import stripy

# These tests just test the mesh construction doesn't fail at the moment
# ToDo: add meaningful assertions

# It's hard to parameterise these as the calls are a little
# bit different in each case

def test_kd_tree_building_cartesian():
    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.2,
                                                    random_scale=0.01, refinement_levels=0, tree=True)

def test_cart_square_border():
    mesh = stripy.cartesian_meshes.square_border((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=0)
    mesh = stripy.cartesian_meshes.square_border((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=2)
    mesh = stripy.cartesian_meshes.square_border((-1.0,1.0,-1.0,1.0), 0.1, 0.2, refinement_levels=0)


def test_elliptical_mesh():

    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=0)
    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=2)
    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.2, refinement_levels=0)
    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.1, random_scale=0.01, refinement_levels=0)
    mesh = stripy.cartesian_meshes.elliptical_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.2, random_scale=0.01, refinement_levels=0)


def test_cart_square_mesh():
    mesh = stripy.cartesian_meshes.square_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=0)
    mesh = stripy.cartesian_meshes.square_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.1, refinement_levels=2)
    mesh = stripy.cartesian_meshes.square_mesh((-1.0,1.0,-1.0,1.0), 0.1, 0.2, refinement_levels=0)


def test_cart_random_mesh():
    mesh = stripy.cartesian_meshes.random_mesh((-1.0,1.0,-1.0,1.0), 5000)
    mesh = stripy.cartesian_meshes.random_mesh((-1.0,1.0,-1.0,1.0), 10000)
    mesh = stripy.cartesian_meshes.random_mesh((-1.0,1.0,-1.0,1.0), 20000)
