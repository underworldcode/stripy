import pytest
import stripy

# These tests just test the mesh construction doesn't fail at the moment
# ToDo: add meaningful assertions

# It's hard to parameterise these as the calls are a little
# bit different in each case


def test_warped_xy_mesh_sphere():

    mesh = stripy.hybrid_st_meshes.warped_xy_mesh_sphere(res_lon=93, res_lat=28, epsilon=0.001)

    r_array = mesh.XX**2 + mesh.YY**2 + mesh.ZZ**2

    assert r_array.min() > 0.99999
    assert r_array.max() < 1.00001
