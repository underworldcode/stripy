"""
Copyright 2017 Louis Moresi

This file is part of Litho1pt.

Stripy is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

Stripy is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Litho1pt0.  If not, see <http://www.gnu.org/licenses/>.
"""

from . import documentation

import os.path as path
import stripy as stripy
import numpy as np
import subprocess
import shutil as shutil
import pkg_resources

from collections import OrderedDict as _OrderedDict

DATA_PATH  = pkg_resources.resource_filename('litho1pt0', 'data/')
L1pt0_FILE = pkg_resources.resource_filename('litho1pt0', 'data/litho_data.npz')
C1pt0_FILE = pkg_resources.resource_filename('litho1pt0', 'data/Crust1pt0-regionalisation.npz')

##
## These are the data structures needed to decode the litho1.0 and crust1.0 files
##

l1_layer_decode = _OrderedDict()
l1_layer_decode["ASTHENO-TOP"   ] =  0
l1_layer_decode["LID-BOTTOM"    ] =  1
l1_layer_decode["LID-TOP"       ] =  2
l1_layer_decode["CRUST3-BOTTOM" ] =  3
l1_layer_decode["CRUST3-TOP"    ] =  4
l1_layer_decode["CRUST2-BOTTOM" ] =  5
l1_layer_decode["CRUST2-TOP"    ] =  6
l1_layer_decode["CRUST1-BOTTOM" ] =  7
l1_layer_decode["CRUST1-TOP"    ] =  8
l1_layer_decode["SEDS3-BOTTOM"  ] =  9
l1_layer_decode["SEDS3-TOP"     ] =  10
l1_layer_decode["SEDS2-BOTTOM"  ] =  11
l1_layer_decode["SEDS2-TOP"     ] =  12
l1_layer_decode["SEDS1-BOTTOM"  ] =  13
l1_layer_decode["SEDS1-TOP"     ] =  14
l1_layer_decode["WATER-BOTTOM"  ] =  15
l1_layer_decode["WATER-TOP"     ] =  16
l1_layer_decode["ICE-BOTTOM"    ] =  17
l1_layer_decode["ICE-TOP"       ] =  18


l1_data_decode = _OrderedDict()
l1_data_decode[ "DEPTH"    ] = 0   # (m)
l1_data_decode[ "DENSITY"  ] = 1   # (kg/m^3)
l1_data_decode[ "VP"       ] = 2   # (m/s)
l1_data_decode[ "VS"       ] = 3   # (m/s)
l1_data_decode[ "QKAPPA"   ] = 4   #   -
l1_data_decode[ "QMU"      ] = 5   #   -
l1_data_decode[ "VP2"      ] = 6   # (m/s)
l1_data_decode[ "VS2"      ] = 7   # (m/s)
l1_data_decode[ "ETA"      ] = 8   #   -

c1_region_descriptor = []
c1_region_descriptor.append("Platform")
c1_region_descriptor.append("Slow, thin Platform")
c1_region_descriptor.append("Archean (Antarctica)")
c1_region_descriptor.append("Early Archean")
c1_region_descriptor.append("Late Archean")
c1_region_descriptor.append("Early/mid  Proter.,")
c1_region_descriptor.append("Early/mid  Proter. (Antarctica, slow)")
c1_region_descriptor.append("Late Proter.")
c1_region_descriptor.append("Slow late Proter.")
c1_region_descriptor.append("Island arc")
c1_region_descriptor.append("Forearc")
c1_region_descriptor.append("Continental arc")
c1_region_descriptor.append("Slow continental arc")
c1_region_descriptor.append("Extended crust")
c1_region_descriptor.append("Fast extended crust (Antarctica)")
c1_region_descriptor.append("Orogen (Antarctica), thick upper, thin lower crust")
c1_region_descriptor.append("Orogen, thick upper crust, very thin lower crust")
c1_region_descriptor.append("Orogen, thick upper crust, fast middle crust")
c1_region_descriptor.append("Orogen with slow lower crust (Andes)")
c1_region_descriptor.append("Slow orogen (Himalaya)")
c1_region_descriptor.append("Margin-continent/shield  transition")
c1_region_descriptor.append("Slow Margin/Shield (Antarctica)")
c1_region_descriptor.append("Rift")
c1_region_descriptor.append("Phanerozoic")
c1_region_descriptor.append("Fast Phanerozoic (E. Australia, S. Africa, N. Siberia)")
c1_region_descriptor.append("Normal oceanic")
c1_region_descriptor.append("Oceans 3 Myrs and younger")
c1_region_descriptor.append("Melt affected o.c. and oceanic plateaus")
c1_region_descriptor.append("Continental shelf")
c1_region_descriptor.append("Continental slope, margin, transition")
c1_region_descriptor.append("Inactive ridge, Alpha Ridge")
c1_region_descriptor.append("Thinned cont. crust, Red Sea")
c1_region_descriptor.append("Oceanic plateau with cont. crust")
c1_region_descriptor.append("Caspian depression")
c1_region_descriptor.append("Intermed. cont./oc. crust, Black Sea")
c1_region_descriptor.append("Caspian Sea oceanic")

# This is only useful when converting the original data
_c1_region_decode = _OrderedDict()
_c1_region_decode["D-"]= 0
_c1_region_decode["E-"]= 1
_c1_region_decode["F-"]= 2
_c1_region_decode["G1"]= 3
_c1_region_decode["G2"]= 4
_c1_region_decode["H1"]= 5
_c1_region_decode["H2"]= 6
_c1_region_decode["I1"]= 7
_c1_region_decode["I2"]= 8
_c1_region_decode["J-"]= 9
_c1_region_decode["K-"]= 10
_c1_region_decode["L1"]= 11
_c1_region_decode["L2"]= 12
_c1_region_decode["M-"]= 13
_c1_region_decode["N-"]= 14
_c1_region_decode["O-"]= 15
_c1_region_decode["P-"]= 16
_c1_region_decode["Q-"]= 17
_c1_region_decode["R1"]= 18
_c1_region_decode["R2"]= 19
_c1_region_decode["T-"]= 20
_c1_region_decode["U-"]= 21
_c1_region_decode["X-"]= 22
_c1_region_decode["Z1"]= 23
_c1_region_decode["Z2"]= 24
_c1_region_decode["A1"]= 25
_c1_region_decode["A0"]= 26
_c1_region_decode["B-"]= 27
_c1_region_decode["C-"]= 28
_c1_region_decode["S-"]= 29
_c1_region_decode["V1"]= 30
_c1_region_decode["V2"]= 31
_c1_region_decode["W-"]= 32
_c1_region_decode["Y1"]= 33
_c1_region_decode["Y2"]= 34
_c1_region_decode["Y3"]= 35


###
###  Initialise the module
###

_l1_data = np.load(L1pt0_FILE)
_litho_data = _l1_data["litho1_all_data"]
_mesh_coords = _l1_data["litho1_mesh_coords"]
_interpolator = stripy.sTriangulation(np.radians(_mesh_coords.T[2]), np.radians(_mesh_coords.T[0]))

_c1_data = np.load(C1pt0_FILE)
_c1_crust_type_lat_lon = _c1_data['latlonDescriptor']

def layer_depth( lat, lon, layerID="LID-BOTTOM"):
    """Returns layer depth at lat / lon (degrees)
    where lat/lon may be arrays (of equal size).
    Depths are returned in metres.
    """

    ## Must wrap longitude from 0 to 360 ...

    lon1 = np.array(lon)%360.0
    lat1 = np.array(lat)

    # ## Must wrap longitude from -180 to 180 ...
    #
    # lon1[np.where(lon1 > 180.0)] = 360.0 - lon1[np.where(lon1 > 180.0)]
    #
    data, err = _interpolator.interpolate( np.radians(lon1), np.radians(lat1),
                                      _litho_data[l1_layer_decode[layerID], l1_data_decode["DEPTH"]], order=1 )

    return data


def crust_type_at(lat=None, lon=None):
    """
    lat, lon (degrees)
    """
    # Get lon into appropriate format

    lats = np.array(lat)
    lons = np.array(lon%360)

    iVals = ((90.0-lats)%180).astype(np.int)
    jVals = (lons%360.0).astype(int)

    # i = int((-lat+90.0)%180)
    # j = int(lon)

    t = _c1_crust_type_lat_lon[iVals,jVals]

    # t = _c1_crust_type_lat_lon[i,j]
    # des = litho.c1_region_descriptor[t]

    return t





    return t


def property_at_lat_lon_depth_points(lat, lon, depth, quantity_ID="DENSITY"):
    """
    Lat / Lon are in degrees
    Depth in km
    quantity_ID needs to match those in the litho1 model

    Points that are not found are given the out-of-range value of -99999
    """

    nlayers = len(l1_layer_decode)
    shape = np.array(lon).shape

    lon1   = np.array(lon).reshape(-1)
    lat1   = np.array(lat).reshape(-1)
    depth1 = np.array(depth).reshape(-1)
    point_properties = np.empty_like(depth1)

    layer_depths     = np.empty((nlayers, lat1.shape[0]))
    layer_properties = np.ones((nlayers+1, lat1.shape[0])) * -99999.0  # if the point is not found it will end up in the overshoot !

    # should assert here that the three arrays are equal size

    for i in range(0, nlayers, 1 ):
        layer_depths[i], err = _interpolator.interpolate( lon1 * np.pi / 180.0, lat1 * np.pi / 180.0,
                                      _litho_data[i,l1_data_decode["DEPTH"]], order=1)
        layer_properties[i], err = _interpolator.interpolate( lon1 * np.pi / 180.0, lat1 * np.pi / 180.0,
                                      _litho_data[i,l1_data_decode[quantity_ID]], order=1)


    A = -layer_depths
    B = -depth1 * 1000.0
    C = divmod(np.searchsorted(A.ravel(), B), A.shape[1])[0] # YEP - this seems to be the best way !!

    # point_properties = np.diag(layer_properties[C[:],:])

    point_properties = np.empty_like(depth1)
    for i,layer in enumerate(C):

        point_properties[i] = layer_properties[layer,i]

    return C, point_properties.reshape(shape)

def property_on_depth_profile(lat, lon, depths, quantity_ID="DENSITY"):
    """
    Lat / Lon are in degrees
    Depth in km
    quantity_ID needs to match those in the litho1 model

    Points that are not found are given the out-of-range value of -99999

    """

    lon1 = np.array((lon,))
    lat1 = np.array((lat,))
    depths1 = np.array(depths)

    nlayers = len(l1_layer_decode)
    point_properties = np.empty_like(depths)

    layer_depths     = np.empty((nlayers))
    layer_properties = np.ones((nlayers+1)) * -99999.0  # if the point is not found it will end up in the overshoot !

    # should assert here that the three arrays are equal size

    for i in range(0, nlayers, 1 ):
        layer_depths[i], err = _interpolator.interpolate( np.radians(lon1), np.radians(lat1),
                                      _litho_data[i,l1_data_decode["DEPTH"]], order=1)
        layer_properties[i], err = _interpolator.interpolate( np.radians(lon1), np.radians(lat1),
                                      _litho_data[i,l1_data_decode[quantity_ID]], order=1)


    A = -layer_depths
    B = -depths1 * 1000.0
    C = np.searchsorted(A, B)

    point_properties = np.empty_like(depths)
    for i,layer in enumerate(C):

        point_properties[i] = layer_properties[layer]

    return C, point_properties



## This stuff makes a lot less sense for the package ... should be in a setup script.


def preprocess_raw_litho1_data(model_path, truncated_model_path):

    print ("Running text-file processing script (this may take a while)")
    truncation_script = path.join(path.dirname(__file__),"scripts","truncate_litho1_model_files.sh")
    output = subprocess.check_output([truncation_script, "-r", model_path, "-t", truncated_model_path])

    # now copy the coordinate files as well
    print ("Copying node coordinate data")

    shutil.copy(path.join(model_path,"Icosahedron_Level7_LatLon_mod.txt" ),
                path.join(truncated_model_path,"Icosahedron_Level7_LatLon_mod.txt" ) )

    print ("Done ! ")
    print ("Processed model files can be found in truncated_model_path as *.model_tr ")

    return


def process_raw_litho1_data(model_path):
    """
    Don't forget to strip the model data first
    truncate_raw_litho1_data(model path, truncated_model_path)
    """

    npoints = 40962
    nlayers = len(l1_layer_decode)
    nentries = 9

    litho_data = np.ones((nlayers, nentries, npoints)) * -99999.0

    for model in range(0,npoints):
        model_name = "node"+str(model+1)+".model_tr"
        file = path.join(model_path, model_name)

        thisnodedata  = np.loadtxt(file, usecols=np.arange(0,9), comments="nlayer")
        thisnodeident = np.loadtxt(file, usecols=(9,), dtype=str, comments="nlayer")

        if (model % 1000 == 0):
            print ("Reading node ", model)

        for i, ident in enumerate(thisnodeident):
            litho_data[l1_layer_decode[ident], :,  model] = thisnodedata[i,:]

        ## Post process - not all layers in the model are populated (gives -99999 which is always illegal)
        ## For depth, it makes more sense to have the layer simply have no thickness but appear in the default order

        for layer in range(9,nlayers):
            missing_entries = np.where(litho_data[layer, l1_data_decode["DEPTH"]] == -99999)
            litho_data[layer, l1_data_decode["DEPTH"], missing_entries] = litho_data[layer-1, l1_data_decode["DEPTH"], missing_entries]

    grid_points_location = path.join(model_path,"Icosahedron_Level7_LatLon_mod.txt")
    litho_points = np.loadtxt( grid_points_location )

    return litho_data, litho_points

def write_processed_litho_data(filename, litho_data, litho_points):
    """
    Ensures that the data is stored in a format which is valid for initialising the class
    """

    np.savez_compressed(filename, litho1_all_data=litho_data, litho1_mesh_coords=litho_points)

    return
