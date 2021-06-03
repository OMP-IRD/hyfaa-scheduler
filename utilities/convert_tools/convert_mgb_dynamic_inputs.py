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

import os, sys, shutil
from datetime import datetime, timedelta
import numpy as np
import netCDF4
from pdb import set_trace
import time

        
def MGB_rain_to_netcdf(pbi_file, output_nc, start_date, nc, dt, nt=None):
    
    if hasattr(start_date,'split'):
        #then it is a string => convert to datetime
        start_date = datetime.strptime(start_date, '%Y-%m-%dT%H:%M:%S') 
    
    #create netCDF file from PBI file
    with open(pbi_file, mode='rb') as pbi_ds, netCDF4.Dataset(output_nc, mode='w') as ds:
        
        #check some file dimensions
        pbi_ds.seek(0,2)
        nentries = pbi_ds.tell()
        
        if nentries%(4*nc) != 0:
            raise Exception('input file is does not match a nc=%d * nt=? float32 file'%nc)
        ntmax = nentries//(4*nc)
        if nt is None:
            nt = ntmax
        else:
            if nt > ntmax:
                raise Exception('User requested nt size is %d but file only contains %d lines'%(nt, ntmax))
        pbi_ds.seek(0,0)
        ar_loc = np.empty((nc, nt), dtype='f4')
        for ii in range(nt):
            ar_loc[:,ii] = np.frombuffer(pbi_ds.read(4*nc), dtype='f4')
        
        #convert
        ds.createDimension('n_meshes', nc)
        ds.createDimension('nt', nt)
        var = ds.createVariable('dates', 'f8', ('nt', ), zlib=True, complevel=4, shuffle=True)
        var[:] = np.arange(nt, dtype=np.float64)*dt + (start_date - datetime(1950,1,1)).total_seconds()/86400.0
        var.long_name = 'dates in CNES julian days (days since 01/01/1950)'
        var = ds.createVariable('rain', 'f4', ('n_meshes', 'nt'), zlib=True, complevel=4, shuffle=True)
        var[:] = ar_loc
        var.long_name = 'rain data'

        
    print('PBI file %s successfully converted to netCDF\n  => %s'%(pbi_file, output_nc))



if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to convert MGB code input files to the new (netCDF) format')
    parser.add_argument("-input", type=str, help="input folder containing MGB static files in old MGB input format")
    parser.add_argument("-output", type=str, help="output folder => will contain files in new format")
    parser.add_argument("-start_date", type=str, help="First date of input file in YYYY-MM-DDTHH:MM:SS format")
    parser.add_argument("-nc", type=int, help="number of cells in input file")
    parser.add_argument("-dt", type=float, help="time step in days")
    parser.add_argument("--nt", type=int, help="Optional: number of times in input file")
    args = parser.parse_args()

    MGB_rain_to_netcdf(args.input, args.output, args.start_date, args.nc, args.dt, nt=args.nt)


