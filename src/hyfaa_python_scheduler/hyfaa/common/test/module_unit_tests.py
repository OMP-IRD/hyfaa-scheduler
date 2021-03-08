#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################################
#  This code has been developped by Magellium SAS
#
#  Licensing:
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>
#######################################################################

from hyfaa.common.common_functions import *
from hyfaa.common.test.common_test_functions import *


def test_lonlat_selection():
    
    geo_sel = {'lonmin': -12., 'lonmax': 16., 'latmin': 4., 'latmax': 25.}
    nlon_in, nlat_in = 360, 181
    lon_in, lat_in = np.linspace(0., 360.-360./nlon_in, nlon_in), np.linspace(-90., 90., nlat_in)
    data_in = np.zeros((nlat_in, nlon_in), dtype=np.float64)
    data_in[lat_in<geo_sel['latmin'],:] = 1.
    data_in[lat_in>geo_sel['latmax'],:] = 1.
    data_in[:,np.logical_and(lon_in>geo_sel['lonmax'], lon_in<geo_sel['lonmin']+360.)] = 1.
    lonlat_selection, lat_selection, lon_selection = compute_lonlat_selection(lon_in, lat_in, geo_sel)
    data_sel = grid_select_lonlat(data_in, lonlat_selection)
    
    check_condition_utest('Check that grid selection is a perfect match for exact bounds', all([np.prod(np.shape(data_sel))==np.sum(data_in==0.), \
        np.max(data_sel) == 0.]))
    
    geo_sel2 = {'lonmin': -12.1, 'lonmax': 16., 'latmin': 4., 'latmax': 25.}
    data_in2 = np.zeros((nlat_in, nlon_in), dtype=np.float64)
    data_in2[lat_in<geo_sel2['latmin'],:] = 1.
    data_in2[lat_in>geo_sel2['latmax'],:] = 1.
    data_in2[:,np.logical_and(lon_in>geo_sel2['lonmax'], lon_in<geo_sel2['lonmin']+360.)] = 1.
    lonlat_selection, lat_selection, lon_selection = compute_lonlat_selection(lon_in, lat_in, geo_sel2)
    data_sel2 = grid_select_lonlat(data_in2, lonlat_selection)
    
    check_condition_utest('Check that grid selection is larger for inexact bounds', np.prod(np.shape(data_sel2))==np.prod(np.shape(data_sel))+np.shape(data_sel)[0])
    

if __name__ == '__main__':
    
    test_lonlat_selection()

