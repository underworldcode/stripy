import pytest
import stripy

# These tests just test the mesh construction doesn't fail at the moment
# ToDo: add meaningful assertions

# It's hard to parameterise these as the calls are a little
# bit different in each case


def test_hybrid_latlon_mesh():

    mesh = stripy.cartesian_meshes.hybrid_latlon_sphere(res_lon=93, res_lat=28, epsilon=0.001)

    
