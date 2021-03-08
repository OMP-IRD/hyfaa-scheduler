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
from scipy.interpolate import RegularGridInterpolator
from calendar import monthrange
import sqlite3
import pandas
import ftputil

from hyfaa.common.test.common_test_functions import *


        

def extract_gsmap_data_to_netcdf_grid(gsmap_file, output_nc, geo_selection=None, temp_dir=None):
    
    if gsmap_file.split('.')[-1] == 'gz':
        file_extracted = True
        if temp_dir is None:
            if 'TMPDIR' in os.environ:
                temp_dir = os.environ['TMPDIR']
            else:
                temp_dir = os.path.abspath('.')
        os.makedirs(temp_dir, exist_ok=True)
        temp_dir_session = tempfile.mkdtemp(dir=temp_dir, prefix='gsmapex')
        gsmap_unzipped_file = os.path.join(temp_dir_session, 'gsmap.dat')
        subprocess.check_call('gzip -dck %s > %s'%(gsmap_file, gsmap_unzipped_file), shell=True)
    else:
        file_extracted = False
        gsmap_unzipped_file = gsmap_file
    
    date_file = datetime.strptime(os.path.basename(gsmap_file).split('.')[1], '%Y%m%d')
    nlat, nlon = 1200, 3600
    dlon, dlat = 0.1, 0.1
    lat = np.linspace(-60.+dlat/2., 60.-dlat/2., nlat)
    lon = np.linspace(dlon/2., 360.-dlon/2., nlon)
    
    #get spatial selection grid parameters (lonlat_selection dict) and lon,lat vectors reduced to selection
    if geo_selection is not None:
        lonlat_selection, lat_selection, lon_selection = compute_lonlat_selection(lon, lat, geo_selection)
    else:
        lat_selection, lon_selection = lat, lon
    
    #create netCDF file from PBI file
    with open(gsmap_unzipped_file, mode='rb') as gsmap_ds:
        
        #check some file dimensions
        gsmap_ds.seek(0,2)
        nentries = gsmap_ds.tell()
        assert nentries == 4*nlon*nlat, 'wrong number of bytes in gsmap file %s (unzipped %s):  %d, expected %d'%(gsmap_file, gsmap_unzipped_file, nentries, 4*nlon*nlat)
        gsmap_ds.seek(0,0)
        rain_array = np.frombuffer(gsmap_ds.read(4*nlon*nlat), dtype='f4')
    rain_array = rain_array.reshape((nlat, nlon))[::-1,:]
    if geo_selection is not None:
        rain_array = grid_select_lonlat(rain_array, lonlat_selection)
        
    rain_array *= 24 #mm/hr to mm/day
        
    with netCDF4.Dataset(output_nc, mode='w') as ds:
        ds.createDimension('nlon', len(lon_selection))
        ds.createDimension('nlat', len(lat_selection))
        var = ds.createVariable('longitude', 'f8', ('nlon', ), zlib=True, complevel=4, shuffle=True)
        var[:] = lon_selection
        var = ds.createVariable('latitude', 'f8', ('nlat', ), zlib=True, complevel=4, shuffle=True)
        var[:] = lat_selection
        var = ds.createVariable('rain', 'f4', ('nlat', 'nlon'), zlib=True, complevel=4, shuffle=True)
        var[:] = rain_array
        var.unit = 'mm/day'
        var.long_name = 'rain data'
        ds.setncattr('date', date_file.strftime('%Y-%m-%d'))
        ds.setncattr('information', 'this file contains GSMAP rain data for a full day')
        
    if file_extracted:
        shutil.rmtree(temp_dir_session)
        

class gsmapFTP():
    
    __ftp_adress = 'hokusai.eorc.jaxa.jp'
    __datapath = '/realtime_ver/v7/daily0.1_G/00Z-23Z'
    __default_username = 'rainmap'
    __default_password= 'Niskur+1404'
    __host = None
    __tries_if_failed = 5
    
    def __init__(self, username=None, password=None):
        if username is None:
            self.__username = self.__default_username
            self.__password = self.__default_password
        else:
            assert password is not None
            self.__username = username
            self.__password = password
        self.__connect__()
        
            
    def __connect__(self):
        try:
            self.__host.close()
        except:
            pass
        self.__host = ftputil.FTPHost(self.__ftp_adress, self.__username, self.__password)
        
        
    def listdir(self, path):
        try:
            list_out = self.__host.listdir(path)
        except:
            self.__connect__()
            list_out = self.__host.listdir(path)
        if list_out is None:
            list_out = []
        return list_out
        
    def listmonths(self):
        return self.listmonths_from_rawlist(self.listdir(self.__datapath))
        
    @staticmethod
    def listmonths_from_rawlist(rawlist):
        return sorted([(int(el[0:4]), int(el[4:6])) for el in rawlist])
        
    def listfiles_in_month(self, year, month):
        fol_month = os.path.join(self.__datapath, '%04d%02d'%(year, month))
        return self.listfiles_in_month_from_rawlist(self.listdir(fol_month), fol_month, year, month)
        
    @staticmethod
    def listfiles_in_month_from_rawlist(rawlist, fol_month, year, month):
        dico_out = dict()
        for filename in rawlist:
            datetime_loc = datetime.strptime(filename.split('.')[1], '%Y%m%d')
            assert (datetime_loc.year, datetime_loc.month) == (year, month)
            dico_out[datetime_loc] = os.path.join(fol_month, filename)
        return dico_out
        
        
    def download_and_extract(self, filepath_remote, filepath_local, geo_selection=None, temp_dir=None, tries_if_failed=None):
        
        if tries_if_failed is None:
            tries_if_failed = self.__tries_if_failed
        assert tries_if_failed >= 1
        
        if temp_dir is None:
            if 'TMPDIR' in os.environ:
                temp_dir = os.environ['TMPDIR']
            else:
                temp_dir = os.path.abspath('.')
        os.makedirs(temp_dir, exist_ok=True)
        
        succeeded = False
        itries = 0
        while(True):
            itries += 1
            temp_dir_session = tempfile.mkdtemp(dir=temp_dir, prefix='gsmapdl')
            dl_path = os.path.join(temp_dir_session, os.path.basename(filepath_remote))
            try:
                #download
                self.__host.download(filepath_remote, dl_path)
                #extract
                extract_gsmap_data_to_netcdf_grid(dl_path, filepath_local, geo_selection=geo_selection, temp_dir=temp_dir_session)
                succeeded = True
                shutil.rmtree(temp_dir_session)
                break
            except:
                shutil.rmtree(temp_dir_session)
                print('Download or extraction failed for file %s'%filepath_remote)
                self.__connect__()
            if itries >= tries_if_failed:
                break
        return succeeded
        
        
        
   
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
def retrieve_forcing_data(missing_dates, gsmap_folder_local=None, geo_selection=None, nprocs=1, temp_folder_base=None, verbose=0):
    
    if temp_folder_base is None:
        temp_folder_base = os.path.abspath('.')
    os.makedirs(temp_folder_base, exist_ok=True)
    temp_folder = tempfile.mkdtemp(dir=temp_folder_base, prefix='rfd_')
    
    try:
        missing_dates_set = set(missing_dates)
        year_months_days = set([(el.year, el.month, el.day) for el in missing_dates])
        year_months = set([(el[0], el[1]) for el in year_months_days])
        
        if gsmap_folder_local is None:
            gsmap_ftp = gsmapFTP()
            year_months_available = set(gsmap_ftp.listmonths())
        else:
            year_months_available = set(gsmapFTP.listmonths_from_rawlist(os.listdir(gsmap_folder_local)))
        
        files_info = []
        for year, month in sorted(list(year_months & year_months_available)):
            if gsmap_folder_local is None:
                dico_files_ftp = gsmap_ftp.listfiles_in_month(year, month)
            else:
                dico_files_ftp = gsmapFTP.listfiles_in_month_from_rawlist(os.listdir(os.path.join(gsmap_folder_local, '%04d%02d'%(year, month))), \
                    os.path.join(gsmap_folder_local, '%04d%02d'%(year, month)), year, month)
            for day in range(1, monthrange(year, month)[1]+1):
                datetime_loc = datetime(year, month, day)
                if (year, month, day) not in year_months_days:
                    continue
                if datetime_loc not in dico_files_ftp:
                    continue
                #download daily file to temp folder
                output_nc = os.path.join(temp_folder, 'gsmap_%04d%02d%02d.nc'%(year, month, day))
                if gsmap_folder_local is None:
                    gsmap_ftp.download_and_extract(dico_files_ftp[datetime_loc], output_nc, geo_selection=geo_selection, temp_dir=temp_folder)
                else:
                    extract_gsmap_data_to_netcdf_grid(dico_files_ftp[datetime_loc], output_nc, geo_selection=geo_selection, temp_dir=temp_folder)
                if verbose >= 2:
                    print('GSMAP grid data for date %s successfully extracted to %s'%(datetime_loc.strftime('%Y-%m-%dT%H:%M:%S'), output_nc))
                yield {'file_path': output_nc, 'date_data': datetime_loc, 'product_type': 'analysis', 'grid_status': 'complete'}
    finally:
        shutil.rmtree(temp_folder)
############################





############################
#interpolate_forcing_data
def interpolate_forcing_data(files_info_in, lon, lat, nprocs=1, temp_folder_base=None, verbose=0):
    
    if temp_folder_base is None:
        temp_folder_base = os.path.abspath('.')
    os.makedirs(temp_folder_base, exist_ok=True)
    temp_folder = tempfile.mkdtemp(dir=temp_folder_base, prefix='ifd_')
        
    try:
        for file_info_loc in files_info_in:
            file_info_new = copy.deepcopy(file_info_loc)
            file_info_new['file_path'] = os.path.join(temp_folder, 'data_%s_%s.nc'%(date2strtag(file_info_new['date_data']), date2strtag(datetime.now())))
            if verbose >= 2:
                print('Interpolating to mesh: %s -> %s'%(file_info_loc['file_path'], file_info_new['file_path']))
            interpolate_file(file_info_loc['file_path'], file_info_new['file_path'], lon, lat)
            if verbose >= 2:
                print('Successfully interpolated to mesh: %s -> %s'%(file_info_loc['file_path'], file_info_new['file_path']))
            yield file_info_new
    finally:
        shutil.rmtree(temp_folder)
############################





    
    
