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

import os, sys, shutil, re, subprocess
if sys.version_info.major < 3:
    raise Exception('sorry, only python 3 supported')
import time
import numpy as np
from datetime import datetime,timedelta
import netCDF4
import tempfile
import copy

#we tolerate up to 0.01 seconds date difference when comparing dates
dt_date_tolerance_seconds = 1.e-2
dt_date_tolerance_days = dt_date_tolerance_seconds/(24.*3600.)



def compute_size_du(path, human_readable=False):
    """disk usage"""
    if human_readable:
        return subprocess.check_output(['du','-hs', path]).split()[0].decode('utf-8')
    else:
        return subprocess.check_output(['du','-s', path]).split()[0].decode('utf-8')
        

def nc_compress(file_in, file_out=None, complevel=None, verbose=2):
    
    if complevel is None:
        complevel = 4
    complevel = int(complevel)
    assert 0 <= complevel <= 9
    
    inplace = False
    if file_out is None:
        inplace = True
        temp_dir_loc = tempfile.mkdtemp(dir=os.path.abspath('.'), prefix='nccomp')
        file_out = os.path.join(temp_dir_loc, 'nccomp.nc')
        
    if verbose >= 2:
        size_original = float(compute_size_du(file_in))
    subprocess.check_call(['nccopy', '-d', '%d'%complevel, file_in, file_out])
    if verbose >= 2:
        size_new = float(compute_size_du(file_out))
        print('Compression level %d: size changed from %s to %s (%s%% of original)'%(complevel, size_original, size_new, 100.*size_new/size_original))
    
    if inplace:
        shutil.move(file_out, file_in)
        shutil.rmtree(temp_dir_loc)
    
    

def matmul_bigdata(mat1, mat2):
    
    max_size = 1.e5
    
    shp1 = np.shape(mat1)
    shp2 = np.shape(mat2)
    assert shp1[1] == shp2[0]
    assert mat1.dtype == mat2.dtype
    
    #it the matrix is small enough, simply perform the basic matrix operation
   # if shp1[0]*shp1[1]*shp2[1] < max_size:
   #     return np.matmul(mat1, mat2)
        
    #if the matrix is too big, simply split it in smaller part along the biggest dimension
    mat_out = np.empty((shp1[0], shp2[1]), dtype=mat1.dtype)
    if shp1[0] > shp2[1]:
        block_size = max(1, int(np.floor(max_size/(1.*shp1[1]*shp2[1]))))
        it0 = 0
        while(True):
            it1 = it0 + block_size
            if it1 > shp1[0]:
                it1 = shp1[0]
            mat_out[it0:it1,:] = np.matmul(mat1[it0:it1,:], mat2)
            if it1 >= shp1[0]:
                break
            it0 = it1
    else:
        block_size = max(1, int(np.floor(max_size/(1.*shp1[1]*shp1[0]))))
        it0 = 0
        while(True):
            it1 = it0 + block_size
            if it1 > shp2[1]:
                it1 = shp2[1]
            mat_out[:,it0:shp2[1]] = np.matmul(mat1, mat2[:,it0:shp2[1]])
            if it1 >= shp2[1]:
                break
            it0 = it1
    return mat_out


def multiply_bigdata(mat1, mat2):
    max_size = 1.e5

    shp1 = np.shape(mat1)
    shp2 = np.shape(mat2)
    assert shp1 == shp2
    assert mat1.dtype == mat2.dtype

    # it the matrix is small enough, simply perform the basic matrix operation
    if shp1[0] * shp1[1] * shp2[1] < max_size:
        return np.multiply(mat1, mat2)

    # if the matrix is too big, simply split it in smaller part along the biggest dimension
    mat_out = np.empty((shp1), dtype=mat1.dtype)
    if shp1[0] > shp1[1]:
        block_size = max(1, int(np.floor(max_size / (1. * shp1[0] * shp1[1]))))
        it0 = 0
        while (True):
            it1 = it0 + block_size
            if it1 > shp1[0]:
                it1 = shp1[0]
            mat_out[it0:it1, :] = np.multiply(mat1[it0:it1, :], mat2[it0:it1,:])
            if it1 >= shp1[0]:
                break
            it0 = it1
    else:
        block_size = max(1, int(np.floor(max_size / (1. * shp1[1] * shp1[0]))))
        it0 = 0
        while (True):
            it1 = it0 + block_size
            if it1 > shp1[1]:
                it1 = shp1[1]
            mat_out[:, it0:shp1[1]] = np.multiply(mat1[:, it0:shp1[1]], mat2[:, it0:shp1[1]])
            if it1 >= shp1[1]:
                break
            it0 = it1
    return mat_out

def internet_available(hostname="www.google.com"):
    internet_test_dir = 'internet_test_dir'
    os.makedirs(internet_test_dir, exist_ok=True)
    try:
        subprocess.check_call('wget %s -q -P %s -T 2'%(hostname, internet_test_dir), shell=True)
        return True
    except:
        pass
    if os.path.exists(internet_test_dir):
        os.system('rm -R -f %s'%internet_test_dir)
    return False

#############################################
#basic functions
def list_form(list_in):
    if type(list_in) == list:
        list_out = list_in
    else:
        list_out = [list_in]
    return list_out


def check_dict(dico, keys, check_none=True, prefix=None):
    none_keys = []
    not_present_keys = []
    for key in keys:
        if key not in dico:
            not_present_keys.append(key)
        else:
            if dico[key] is None:
                none_keys.append(key)
    msg = ''
    if len(not_present_keys) > 0:
        msg += 'Missing keys: %s'%(', '.join(not_present_keys))
    if (len(none_keys) > 0) and check_none:
        if len(msg) > 0:
            msg += '\n'
        msg += 'None keys: %s'%(', '.join(none_keys))
    if len(msg) > 0:
        if prefix is not None:
            msg = '%s%s'%(prefix, msg)
        raise Exception(msg)


def stack_generator(list_generator, stack_size):
    stack = []
    for elem in list_generator:
        stack.append(elem)
        if len(stack) >= stack_size:
            yield stack
            stack = []
    if len(stack) > 0:
        yield stack
        
        

def list_argsort(li):
    return sorted(range(len(li)), key=lambda k: li[k])
    
def is_sorted(li):
    if len(li) <= 1:
        return True
    return all([li[ii] >= li[ii-1] for ii in range(1,len(li))])


def classic_mask(var, v_min=None, v_max=None):
    if v_min is None:
        v_min = -1e18
    if v_max is None:
        v_max = 1e18
    var = np.ma.masked_invalid(var)
    var = np.ma.masked_outside(var, v_min, v_max)
    return var
    
def nan_array(tu):
    ar = np.empty(tu)
    ar[:] = np.nan
    return ar

    
def mask_logic(l, operation):
    """Input list of masks"""
    if operation not in ['or', 'and']:
        raise Exception('operation parameter must be in [or,and]')
    n = len(l)
    if n == 1:
        return np.copy(l[0])
    elif n > 1:
        if np.any(np.array([np.shape(el) for el in l]) != np.shape(l[0])):
            raise Exception('shapes mismatch')
        mask = np.copy(l[0])
        if operation == 'or':
            for ii in range(1, n):
                mask = np.logical_or(mask, l[ii])
        elif operation == 'and':
            for ii in range(1, n):
                mask = np.logical_and(mask, l[ii])
        return mask
    else:
        raise Exception('Empty list')
        
        
def pad_txt(txt_in, len_aim):
    ncomplete = len_aim-len(txt_in)
    if ncomplete < 0:
        raise Exception('padding problem, text size would be reduced...')
    return ' '*ncomplete + txt_in


###########
#era5 hours
def datetime_to_era5hours(datetime_obj):
    return (datetime_obj-datetime(1900,1,1)).total_seconds()/3600.

def era5hours_to_datetime(era5hours):
    return datetime(1900,1,1)+timedelta(era5hours/24.)

def datetime_to_era5hours_array(datetime_ar):
    return np.array([(datetime_obj-datetime(1900,1,1)).total_seconds()/3600. for datetime_obj in datetime_ar])

def era5hours_to_datetime_array(era5hours_ar):
    return datetime(1900,1,1)+np.array([timedelta(era5hours/24.) for era5hours in era5hours_ar])
###########


###########
#cnes julian days
def datetime_to_julianday(datetime_obj):
    return (datetime_obj-datetime(1950,1,1)).total_seconds()/(24.*3600.)

def julianday_to_datetime(jday):
    return datetime(1950,1,1)+timedelta(jday)

def datetime_to_julianday_array(datetime_ar):
    return np.array([(datetime_obj-datetime(1950,1,1)).total_seconds()/(24.*3600.) for datetime_obj in datetime_ar])

def julianday_to_datetime_array(jday_ar):
    return datetime(1950,1,1) + np.array([timedelta(jday) for jday in jday_ar])
###########

def multiple_split(txt, split_list):
    out = [txt]
    for split_item in split_list:
        out_new = []
        for el in out:
            out_new += el.split(split_item)
        out = [el for el in out_new if len(el) > 0]
    return out
    
##########


def str2date(datetime_str):
    return datetime.strptime(datetime_str, '%Y%m%dT%H%M%SU%f')
    
def date2str(datetime_obj):
    return datetime_obj.strftime('%Y%m%dT%H%M%SU%f')
    
def date2strtag(datetime_obj, mode='microsecond'):
    if mode == 'microsecond':
        return datetime_obj.strftime('%Y%m%dT%H%M%SU%f')
    elif mode == 'second':
        return datetime_obj.strftime('%Y%m%dT%H%M%S')
    elif mode == 'minute':
        return datetime_obj.strftime('%Y%m%dT%H%M')
    elif mode == 'hour':
        return datetime_obj.strftime('%Y%m%dT%H')
    elif mode == 'day':
        return datetime_obj.strftime('%Y%m%d')
    else:
        raise Exception('mode %s unknown'%mode)


def make_necessary_directories(filename):
    dir_filename = os.path.dirname(filename)
    if len(dir_filename) > 0:
        os.makedirs(dir_filename, exist_ok=True)


#lonlat selection
def minmax_id_container(ar, min_val, max_val):
    if min_val < ar[0] or max_val > ar[-1]:
        print(min_val, ar[0], max_val, ar[-1])
        raise Exception('min_val < ar[0] or max_val > ar[-1]')
    return np.where(ar <= min_val)[0][-1], np.where(ar >= max_val)[0][0]+1
    
        
            
def compute_lonlat_selection(lon, lat, geo_selection):
    lat_ids = minmax_id_container(lat, geo_selection['latmin'], geo_selection['latmax'])
    lat_selection = lat[lat_ids[0]:lat_ids[1]]
    dico_out = {'lat': {'idmin': lat_ids[0], 'idmax': lat_ids[1]}}
    lon_ids = minmax_id_container(np.concatenate((lon-360., lon, lon[0:1]+360.), axis=0), geo_selection['lonmin'], geo_selection['lonmax'])
    nlon = len(lon)
    if lon_ids[0] < nlon:
        if lon_ids[1] <= nlon:
            lon_selection = lon[lon_ids[0]:lon_ids[1]]-360.
            dico_out['lon'] = [{'idmin': lon_ids[0], 'idmax': lon_ids[1]}]
        else:
            lon_selection = np.concatenate((lon[lon_ids[0]:nlon]-360., lon[0:lon_ids[1]-nlon]), axis=0)
            dico_out['lon'] = [{'idmin': lon_ids[0], 'idmax': nlon}, {'idmin': 0, 'idmax': lon_ids[1]-nlon}]
    else:
        lon_selection = lon[lon_ids[0]-nlon:lon_ids[1]-nlon]
        dico_out['lon'] = [{'idmin': lon_ids[0]-nlon, 'idmax': lon_ids[1]-nlon}]
    return dico_out, lat_selection, lon_selection
    

def grid_select_lonlat(grid_in, lonlat_selection):
    return np.concatenate(tuple([grid_in[lonlat_selection['lat']['idmin']:lonlat_selection['lat']['idmax'],idlon['idmin']:idlon['idmax']] for idlon in lonlat_selection['lon']]), axis=1)
    


def write_forcing_grid(lon, lat, grid_in, filename, info_dict, complevel=4, verbose=0):
    if 'time' not in info_dict:
        raise Exception('time attribute must be present in info_dict parameter')
    with netCDF4.Dataset(filename, mode='w') as ds:
        #dimensions
        ds.createDimension('longitude', len(lon))
        ds.createDimension('latitude', len(lat))
        for el in info_dict:
            ds.setncattr(el, info_dict[el])
        #variables
        var_loc = ds.createVariable('longitude', 'f8', ('longitude', ), zlib=True, complevel=complevel, shuffle=True)
        var_loc[:] = lon
        var_loc = ds.createVariable('latitude', 'f8', ('latitude', ), zlib=True, complevel=complevel, shuffle=True)
        var_loc[:] = lat
        var_loc = ds.createVariable('rain', 'f4', ('latitude', 'longitude'), zlib=True, complevel=complevel, shuffle=True)
        var_loc.setncattr('unit', 'mm/day')
        var_loc[:] = grid_in
        for el in info_dict:
            ds.setncattr(el, info_dict[el])
    if verbose > 0:
        print('%s written...'%filename)
            
            

