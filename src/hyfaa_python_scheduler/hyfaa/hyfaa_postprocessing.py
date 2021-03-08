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
import copy
from datetime import datetime, timedelta
import numpy as np
import netCDF4
import time

from hyfaa.common.test.common_test_functions import *
from hyfaa.check_input_parameters.check_input_parameters import *
from hyfaa.common.yaml.yaml_parser import load_yaml
from hyfaa.database.hydrostates.hydrostates_db import HydroStates_DBManager


portal_default_round = 3
vars_info_portal = {'water_elevation_catchment': {'info_dict': {'units': 'cm', \
    'long_name': 'Water elevation in catchment relative to river bottom.', \
    'standard_name': 'water_elevation_catchment', \
    'comments': 'Water elevation in MGB model should not be interpreted as real water surface height.'}, 'round': 3}, \
    'streamflow_catchment': {'info_dict': {'units': 'm^3/s', \
    'long_name': 'Flowrate in catchment', \
    'standard_name': 'streamflow_catchment', \
    'comments': 'Flowrate in catchment.'}, 'round': 3}}


    

def hyfaa_postprocessing(yaml_file_or_dict, verbose=None):
    """main scheduler processing function
    
    :param yaml_file_or_dict: python dict or yaml file
    """
    
    time_now = datetime.utcnow()
    
    #check if input is a python dict or a yaml file
    valid_obj = False
    if isinstance(yaml_file_or_dict, dict):
        valid_obj = True
        dico = yaml_file_or_dict
    else:
        if hasattr(yaml_file_or_dict, 'replace'): #pythonic way of checking if this is a python3 string, python2 string or python2 unicode
            if os.path.exists(yaml_file_or_dict):
                valid_obj = True
                dico = load_yaml(yaml_file_or_dict, env_vars=True)
    if not valid_obj:
        raise Exception('input must be a python dict or a path to an existing yaml file')
        
    #check parameters
    dico = check_parameters(dico)
    exec_time = dico['exec_time']
    if verbose is None:
        if 'verbose' in dico:
            verbose = dico['verbose']
    if verbose is None:
        verbose = 0
        
    if verbose >= 1:
        print('Launching HYFAA post-processing...')
        
    #get some dimensions from static_data.nc file
    with netCDF4.Dataset(dico['mgb']['static_data_file']) as ds:
        n_cells = ds.dimensions['n_cells'].size
        mesh_lon = ds.variables['longitude_center'][:]
        mesh_lat = ds.variables['latitude_center'][:]
    
    
    #post-processing science file
    
    #create output directory if it doesn't exist
    if not os.path.exists(os.path.dirname(dico['post_processing']['science_file'])):
        os.makedirs(os.path.dirname(dico['post_processing']['science_file']), exist_ok=True)
        
    #get previous data from science file to avoid reprocessing everything if possible. WARNING: changes in variables requested will trigger full reprocessing
    dico_science_previous, old_science_file = None, None
    vars_expected = {'longitude', 'latitude', 'time'}
    for el in ['control', 'analysis']:
        vars_expected |= {'%s_%s'%(el1, el) for el1 in dico['post_processing']['variables'] + ['time_added_to_hydb']}
    if os.path.exists(dico['post_processing']['science_file']):
        #read old file and then move it to a temporary location (to be erased after new file is created)
        try:
            old_science_file = dico['post_processing']['science_file'] + '_old%s'%time_now.strftime('%Y%m%d%H%M%S')
            shutil.move(dico['post_processing']['science_file'], old_science_file)
            with netCDF4.Dataset(old_science_file, mode='r') as ds_science_old:
                assert ds_science_old.dimensions['n_cells'].size == n_cells, 'n_cells mismatch with old science file: %d vs %d'%(ds.dimensions['n_cells'].size, n_cells)
                assert ds_science_old.dimensions['n_ensemble'].size == dico['n_ensemble'], 'n_ensemble mismatch with old science file: %d vs %d'%(ds.dimensions['n_ensemble'].size, dico['n_ensemble'])
                assert np.all(ds_science_old.variables['longitude'][:] == mesh_lon) and np.all(ds_science_old.variables['latitude'][:] == mesh_lat), 'lonlat mismatch with old science file'
            ds_science_old = netCDF4.Dataset(old_science_file, mode='r')
        except:
            raise
            #if input file is unreadable or variables don't match, reprocess everything
            os.unlink(old_science_file)
            old_science_file = None
        
    #get main information from hydrological state database
    with HydroStates_DBManager(dico['hydrological_states_database_directory'], mode='r', verbose=verbose) as hydrostates_db:
        
        file_info = dict()
        for index, row in hydrostates_db.read_as_pandas_dataframe("SELECT * FROM FILEINFO WHERE file_status=?", params=['added']).iterrows():
            date_loc = hydrostates_db.str2date(row['date_data'])
            if date_loc not in file_info:
                file_info[date_loc] = dict()
            file_info[date_loc][row['type']] = {'file_path': hydrostates_db.get_full_path(row), 'date_added_to_db': hydrostates_db.str2date(row['date_added_to_db'])}
        time_new = sorted(list(file_info.keys()))
        n_time = len(time_new)
            
    #pre-allocate array that will load each variable for each time step
    data_loc = np.ma.masked_invalid(np.zeros((n_cells, dico['n_ensemble']), dtype=np.float64))
            
    #open science and portal file and add information to them either from the old science file (if not modifications since last time), or from hydrological state database
    with netCDF4.Dataset(dico['post_processing']['science_file'], mode='w') as ds_science, netCDF4.Dataset(dico['post_processing']['portal_file'], mode='w') as ds_portal:
            
        #get existing dates and date_added_to_db in old science file
        if old_science_file is not None:
            dict_times_old = {'time_added_to_hydb_%s'%type_loc: ds_science_old.variables['time_added_to_hydb_%s'%type_loc][:] for type_loc in ['control', 'analysis']}
            dict_times_old['time'] = ds_science_old.variables['time'][:]
            
        #initialize science file variables
        ds_science.createDimension('n_time', n_time)
        ds_science.createDimension('n_cells', n_cells)
        ds_science.createDimension('n_ensemble', dico['n_ensemble'])
        ds_science.createVariable('time', np.float64, ['n_time'])
        ds_science.variables['time'][:] = np.ma.masked_invalid(np.array([datetime_to_julianday(el) for el in time_new], dtype=np.float64))
        ds_science.createVariable('longitude', np.float64, ['n_cells'])
        ds_science.variables['longitude'][:] = np.ma.masked_invalid(mesh_lon)
        ds_science.createVariable('latitude', np.float64, ['n_cells'])
        ds_science.variables['latitude'][:] = np.ma.masked_invalid(mesh_lat)
        for type_loc in ['control', 'analysis']:
            var_name = 'time_added_to_hydb_%s'%type_loc
            ds_science.createVariable(var_name, np.float64, ['n_time'])
            ds_science.variables[var_name].setncattr('longname', var_name)
            for var_name in ['%s_%s'%(elem, type_loc) for elem in dico['post_processing']['variables']]:
                ds_science.createVariable(var_name, np.float64, ['n_time', 'n_cells', 'n_ensemble'])
                ds_science.variables[var_name].setncattr('longname', var_name)
        ds_science.setncattr('date_created', time_now.strftime('%Y-%m-%dT%H:%M:%S'))
        
        
        #initialize portal file variables
        ds_portal.createDimension('n_time', n_time)
        ds_portal.createDimension('n_cells', n_cells)
        ds_portal.createVariable('time', np.float64, ['n_time'])
        ds_portal.variables['time'][:] = np.ma.masked_invalid(np.array([datetime_to_julianday(el) for el in time_new], dtype=np.float64))
        ds_portal.variables['time'].setncattr('units', 'days since 1950-01-01 00:00:00.0')
        ds_portal.variables['time'].setncattr('long_name', 'time (days since 1950-01-01)')
        ds_portal.variables['time'].setncattr('standard_name', 'time')
        ds_portal.variables['time'].setncattr('calendar', 'gregorian')
        ds_portal.createVariable('longitude', np.float64, ['n_cells'])
        ds_portal.variables['longitude'][:] = np.ma.masked_invalid(mesh_lon)
        ds_portal.variables['longitude'].setncattr('units', "degrees_east")
        ds_portal.variables['longitude'].setncattr('long_name', "longitude")
        ds_portal.variables['longitude'].setncattr('standard_name', "longitude")
        ds_portal.variables['longitude'].setncattr('comments', "East longitude relative to Greenwich meridian")
        ds_portal.createVariable('latitude', np.float64, ['n_cells'])
        ds_portal.variables['latitude'][:] = np.ma.masked_invalid(mesh_lat)
        ds_portal.variables['latitude'].setncattr('units', "degrees_north")
        ds_portal.variables['latitude'].setncattr('long_name', "latitude")
        ds_portal.variables['latitude'].setncattr('standard_name', "latitude")
        ds_portal.variables['latitude'].setncattr('comments', "Positive latitude is North latitude, negative latitude is South latitude.")
        var_name = 'time_added_to_hydb'
        ds_portal.createVariable(var_name, np.float64, ['n_time'])
        ds_portal.variables[var_name].setncattr('units', 'days since 1950-01-01 00:00:00.0')
        ds_portal.variables[var_name].setncattr('long_name', 'time added to hydrological database (days since 1950-01-01)')
        ds_portal.variables[var_name].setncattr('standard_name', var_name)
        ds_portal.variables[var_name].setncattr('calendar', 'gregorian')
        
        for var_name in dico['post_processing']['variables']:
            for stat_type in ['mean', 'median', 'std']:
                ds_portal.createVariable(var_name + '_' + stat_type, np.float64, ['n_time', 'n_cells'])
                if var_name in vars_info_portal:
                    ds_portal.variables[var_name + '_' + stat_type].setncattr('units', vars_info_portal[var_name]['info_dict']['units'])
                    ds_portal.variables[var_name + '_' + stat_type].setncattr('long_name', vars_info_portal[var_name]['info_dict']['long_name'] + ' (%s)'%stat_type)
                    ds_portal.variables[var_name + '_' + stat_type].setncattr('standard_name', vars_info_portal[var_name]['info_dict']['standard_name'] + '_%s'%stat_type)
                    ds_portal.variables[var_name + '_' + stat_type].setncattr('comments', vars_info_portal[var_name]['info_dict']['comments'] + ' (%s)'%stat_type)
                else:
                    ds_portal.variables[var_name + '_' + stat_type].setncattr('long_name', var_name + ' (%s)'%stat_type)
                    ds_portal.variables[var_name + '_' + stat_type].setncattr('standard_name', var_name + '_%s'%stat_type)
        ds_portal.setncattr('date_created', time_now.strftime('%Y-%m-%dT%H:%M:%S'))
        ds_portal.setncattr('n_ensemble', dico['n_ensemble'])
        
        #iterate over each time step and read information either from the old science file (if unchanged since last time), or from hydrological state database
        it_decile = 1
        tstart_decile = datetime.utcnow()
        for it, time_loc in enumerate(time_new):
            
            ratio_done = (it+1)*1./(len(time_new)*1.)
            if ratio_done*10. >= it_decile:
                print('%d%% done after %s, %s time remaining (estimation)'%(ratio_done*100., datetime.utcnow()-tstart_decile, (datetime.utcnow()-tstart_decile)*(1.-ratio_done)/ratio_done))
                it_decile = max(it_decile+1, int(np.floor(ratio_done*10.)))
                
            time_loc_jday = datetime_to_julianday(time_loc)
            
            best_type_loc = 'control'
            if 'analysis' in file_info[time_loc]:
                best_type_loc = 'analysis'
            
            for type_loc in ['control', 'analysis']:

                #time_loc must be at least in control results, if it is not in analysis step then skip (analysis data for this time step will be masked in science and portal files)
                if type_loc == 'control':
                    assert type_loc in file_info[time_loc]
                elif type_loc not in file_info[time_loc]:
                    continue
                
                date_added_loc = datetime_to_julianday(file_info[time_loc][type_loc]['date_added_to_db'])
                
                #get matching id in old science file
                if old_science_file is not None:
                    ids_match = np.where(np.logical_and(dict_times_old['time'] == time_loc_jday, dict_times_old['time_added_to_hydb_%s'%type_loc] == date_added_loc))[0]
                    if len(ids_match) == 0:
                        it0 = None
                    elif len(ids_match) == 1:
                        it0 = ids_match[0]
                    else:
                        raise Exception('%d match found in old science file, should not happen'%len(ids_match))
                else:
                    it0 = None
                
                if it0 is None:
                    ds_in = netCDF4.Dataset(file_info[time_loc][type_loc]['file_path']) #open database file
                for elem in dico['post_processing']['variables']:
                    var_name = '%s_%s'%(elem, type_loc)
                    if it0 is None:
                        for i_ensemble in range(dico['n_ensemble']):
                            data_loc[:,i_ensemble] = np.ma.masked_invalid(ds_in.variables['%s_%d'%(elem, i_ensemble)][:])
                    else:
                        data_loc = ds_science_old.variables[var_name]['data'][it0,:,:]
                    ds_science.variables[var_name][it,:,:] = data_loc
                    if type_loc == best_type_loc:
                        if elem in vars_info_portal:
                            round_value_loc = vars_info_portal[elem]['round']
                        else:
                            round_value_loc = portal_default_round
                        ds_portal.variables[elem + '_mean'][it,:] = np.round(np.ma.mean(data_loc, axis=1), round_value_loc)
                        ds_portal.variables[elem + '_median'][it,:] = np.round(np.ma.median(data_loc, axis=1), round_value_loc)
                        ds_portal.variables[elem + '_std'][it,:] = np.round(np.ma.std(data_loc, axis=1), round_value_loc)
                    
                if it0 is None:
                    ds_in.close()
                
                ds_science.variables['time_added_to_hydb_%s'%type_loc][it] = date_added_loc
                if type_loc == best_type_loc: 
                    ds_portal.variables['time_added_to_hydb'][it] = date_added_loc
    
    #compress portal file.
    #NB: no need compressing the science file, it will not yield any significant gains in size
    nc_compress(dico['post_processing']['portal_file'])
    
    if old_science_file is not None:
        try:
            ds_science_old.close()
        except:
            pass
        os.unlink(old_science_file)
        
    
    if verbose >= 1:
        print('  Post processing complete in %s'%(datetime.utcnow()-time_now))
    
    
    
    
    
def test_main_program():
    
    print_utest_message('Test post processing')
    

    
    
if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="HYFAA post-processing chain", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--test", action='store_true', help="launch unit tests")
    parser.add_argument("--input_yaml_file", type=str, help="path to input yaml file")  
    parser.add_argument("--verbose", type=int, help="verbose level overload, default to the one contained in the yaml file")  
    args = parser.parse_args()

    if args.test:
        test_main_program()
    else:
        assert args.input_yaml_file is not None
        hyfaa_postprocessing(args.input_yaml_file, verbose=args.verbose)

