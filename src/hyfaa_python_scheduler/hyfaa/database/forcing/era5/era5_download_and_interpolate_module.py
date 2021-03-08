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
from calendar import monthrange
from scipy.interpolate import RegularGridInterpolator
import sqlite3
import pandas
import cdsapi
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

from hyfaa.common.test.common_test_functions import *


    


def create_daily_files(filename, daily_dates, temp_folder_base=None, geo_selection=None, verbose=0):
    
    variable_name_dict={'lon': 'longitude', 'lat': 'latitude', 'grid': 'tp'}
    dim_names_lonlat={'lon': 'longitude', 'lat': 'latitude'}
    
    if temp_folder_base is None:
        if 'TMPDIR' in os.environ:
            temp_folder_base = os.environ['TMPDIR']
        else:
            temp_folder_base = os.path.abspath('.')
    temp_folder = tempfile.mkdtemp(dir=temp_folder_base, prefix='crdaily')
    os.makedirs(temp_folder, exist_ok=True)

        
    with netCDF4.Dataset(filename) as ds:
          
        times_dt = [datetime(1900,1,1)+timedelta(time_loc/24.) for time_loc in ds.variables['time'][:]]
        times_index = {el: ii for ii, el in enumerate(times_dt)}
        lon = ds.variables[variable_name_dict['lon']][:]
        lat = ds.variables[variable_name_dict['lat']][:]
        if lat[1] < lat[0]:
            reverse_lat = True
            lat = lat[::-1]
        else:
            reverse_lat = False
        dims = ds.variables[variable_name_dict['grid']].dimensions
        if dims[-2] != dim_names_lonlat['lat'] or dims[-1] != dim_names_lonlat['lon']:
            raise Exception('expecting last dimensions of grid to be (lat, lon)')
        
        #get spatial selection grid parameters (lonlat_selection dict) and lon,lat vectors reduced to selection
        if geo_selection is not None:
            lonlat_selection, lat_selection, lon_selection = compute_lonlat_selection(lon, lat, geo_selection)
        else:
            lat_selection, lon_selection = lat, lon

        for date_loc in daily_dates:
            
            filename_loc = os.path.join(temp_folder, 'forcing_grid_%s_%s.nc'%(date2strtag(date_loc), date2strtag(datetime.now())))
            
            grid_out = None
            for hour in range(24):
                time_loc = date_loc + timedelta(hour/24.)
                assert time_loc in times_index, 'missing date %s'%time_loc
                it = times_index[time_loc]
            
                #load grid
                if len(dims) == 4 and ds.variables[variable_name_dict['grid']].dimensions[1] == 'number':
                    grid_in = np.mean(ds.variables[variable_name_dict['grid']][it,:,:,:], axis=0)
                elif len(dims) == 3:
                    grid_in = ds.variables[variable_name_dict['grid']][it,:,:]
                else:
                    raise Exception('only dimension 3 grid (time,lat,lon) or dimension 4 grid (time,number,lat,lon) accepted')
                grid_in *= 1000. #to mm
                    
                #if lat is reversed (from 90 to -90) then reorder grid along the latitude dimension
                if reverse_lat:
                    grid_in = grid_in[::-1,:]
                    
                #geographical selection
                if geo_selection is not None:
                    grid_in = grid_select_lonlat(grid_in, lonlat_selection)
                    
                if grid_out is None:
                    grid_out = grid_in.copy()
                else:
                    grid_out += grid_in
                
            #write forcing grid    
            write_forcing_grid(lon_selection, lat_selection, grid_out, filename_loc, {'time': date2str(date_loc)}, verbose=verbose)
        
            yield {'file_path': filename_loc, 'date_data': date_loc, 'product_type': 'analysis', 'grid_status': 'complete'}


    
                
                
    
def retrieve_era5(product_type, year, month, file_out, days=None, overwrite=False):

    if os.path.exists(file_out):
        if not overwrite:
            return 'existed'
        else:
            os.unlink(file_out)
            print('File %s removed...'%file_out)
            
    if days is None:
        days = list(range(1, monthrange(year, month)[1]+1))
            
    retrieve_dict = {'product_type': product_type, 'format': 'netcdf', 'variable': 'total_precipitation', 'year': '%s'%year, 'month': '%02d'%month, \
        'day': ['%02d'%el for el in days], 'time': ['%02d:00'%ii for ii in range(24)]}
    
    make_necessary_directories(file_out)

    c = cdsapi.Client(url='https://cds.climate.copernicus.eu/api/v2', key='4785:84084f93-49e0-4050-87ab-27c8b9591461', verify=0)
    c.retrieve("reanalysis-era5-single-levels", retrieve_dict, file_out)

    print('Downloaded %s ...'%file_out)
    return 'downloaded'
    

                
def interpolate_file(file_in, file_out, lon, lat):
    geo_selection = {'lonmin': np.min(lon), 'lonmax': np.max(lon), 'latmin': np.min(lat), 'latmax': np.max(lat)}
    with netCDF4.Dataset(file_in, mode='r') as ds_in:
        lonlat_selection, lat_selection, lon_selection = compute_lonlat_selection(ds_in.variables['longitude'][:], ds_in.variables['latitude'][:], geo_selection)
        grid_selection = grid_select_lonlat(ds_in.variables['rain'][:], lonlat_selection)
    f_interp = RegularGridInterpolator((lat_selection, lon_selection), grid_selection, method='nearest')
    make_necessary_directories(file_out)
    with netCDF4.Dataset(file_out, mode='w') as ds_out:
        ds_out.createDimension('n_meshes', len(lon))
        var = ds_out.createVariable('rain', 'f4', ('n_meshes', ), zlib=True, complevel=4, shuffle=True)
        var[:] = f_interp([(lat_loc, lon_loc) for lat_loc, lon_loc in zip(lat, lon)])
          
        
        
                
############################
#retrieve_forcing_data
def retrieve_forcing_data(missing_dates, geo_selection=None, nprocs=1, temp_folder_base=None, verbose=0):
    
    if temp_folder_base is None:
        temp_folder_base = '.'
    os.makedirs(temp_folder_base, exist_ok=True)
    temp_folder = tempfile.mkdtemp(dir=temp_folder_base, prefix='rfd_')
    
    missing_dates = sorted(list(set(missing_dates))) #unicity of dates
    for el in missing_dates:
        assert el == datetime(el.year, el.month, el.day) #make sure everything is at 00:00:00 hour
        
    try:
        for year, month in sorted(list(set([(el.year, el.month) for el in missing_dates]))):
            filename_loc = '%s/era5_%d_%d.nc'%(temp_folder, year, month)
            if verbose > 2:
                print('Retrieving file %s ...'%filename_loc)
            dates_loc = sorted([el for el in missing_dates if (el.year, el.month) == (year, month)])
            status_out = retrieve_era5('reanalysis', year, month, filename_loc, days=[el.day for el in dates_loc])
            if status_out == 'failed':
                raise Exception('Retrieval of file for month %02d/%04d failed'%(month, year))
            if verbose > 2:
                print('Creating daily files from %s ...'%filename_loc)
            dates_loc_retrieved = set()
            for file_info in create_daily_files(filename_loc, dates_loc, temp_folder_base=temp_folder, geo_selection=geo_selection, verbose=0):
                dates_loc_retrieved.add(file_info['date_data'])
                yield file_info
            for date_not_retrieved in sorted(list(set(dates_loc) - dates_loc_retrieved)):
                print('%s not retrieved'%(date_not_retrieved.strftime('%Y-%m-%d')))
    finally:
        shutil.rmtree(temp_folder)
        
        
        
############################





############################
#interpolate_forcing_data
def interpolate_forcing_data(files_info_in, lon, lat, nprocs=1, temp_folder_base=None, verbose=0):
    
    if temp_folder_base is None:
        temp_folder_base = '.'
    os.makedirs(temp_folder_base, exist_ok=True)
    temp_folder = tempfile.mkdtemp(dir=temp_folder_base, prefix='ifd_')
        
    try:
        for file_info_loc in files_info_in:
            file_info_new = copy.deepcopy(file_info_loc)
            file_info_new['file_path'] = '%s/data_%s_%s.nc'%(temp_folder, file_info_new['date_data'].isoformat(), date2strtag(datetime.now()))
            if verbose > 2:
                print('Interpolating to mesh: %s -> %s'%(file_info_loc['file_path'], file_info_new['file_path']))
            interpolate_file(file_info_loc['file_path'], file_info_new['file_path'], lon, lat)
            if verbose > 1:
                print('Successfully interpolated to mesh: %s -> %s'%(file_info_loc['file_path'], file_info_new['file_path']))
            yield file_info_new
    finally:
        shutil.rmtree(temp_folder)




if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="retrieve era5 forcing data", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--date_min", type=str, required=True, help="date min in Y-m-d format")
    parser.add_argument("--date_max", type=str, required=True, help="date max in Y-m-d format")
    parser.add_argument("--output_folder", type=str, required=True, help="output folder")
    parser.add_argument("--temp_dir", type=str, help="temp_dir")
    args = parser.parse_args()
    
    args.date_min = datetime.strptime(args.date_min, '%Y-%m-%d')
    args.date_max = datetime.strptime(args.date_max, '%Y-%m-%d')
    required_dates = [args.date_min + timedelta(el) for el in range((args.date_max-args.date_min).days+1)]
    os.makedirs(args.output_folder, exist_ok=True)
    for file_info in retrieve_forcing_data(required_dates, geo_selection=None, nprocs=1, temp_folder_base=args.temp_dir, verbose=0):
        shutil.move(file_info['file_path'], os.path.join(args.output_folder, 'rain_era5_%s.nc'%file_info['date_data'].strftime('%Y%m%dT%H%M%S')))
    
    
    
