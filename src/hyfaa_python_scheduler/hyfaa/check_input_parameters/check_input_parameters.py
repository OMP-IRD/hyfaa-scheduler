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
from hyfaa.common.yaml.yaml_parser import load_yaml, write_simple_1level_dict_to_yaml_file



def get_lonlat_minmax_from_mgb_static_file(static_file_path):
    if not os.path.exists(static_file_path):
        raise Exception('file %s does not exist'%static_file_path)
    with netCDF4.Dataset(static_file_path) as ds:
        longitudes = ds.variables['longitude_center'][:]
        latitudes = ds.variables['latitude_center'][:]
    return {'latmin': min(latitudes), 'latmax': max(latitudes), 'lonmin': min(longitudes), 'lonmax': max(longitudes)}
        


def check_parameters(dico):
    #1st level
    check_dict(dico, ['init_conditions', 'scheduler_time_step', 'forecast_time_span', 'n_ensemble', 'operational_mode', 'retreatment_time_span', \
        'hydrological_states_database_directory', 'perturb_static_data', 'forcing_source', 'forcing_grid_database_directory', 'forcing_onmesh_database_directory', \
        'forcing_dates_dt_max', 'rain_uncertainty', 'assimilation_database_directory', 'assimilation_sources', 'assim_params_file', \
        'mgb', 'model_min_time_step', 'post_processing'], check_none=True)
    check_dict(dico, ['exec_time', 'forcing_grid_geo_selection', 'nprocs', 'temporary_files_directory'], check_none=False, prefix=None)
    if dico['temporary_files_directory'] is None:
        if 'TMPDIR' in os.environ:
            dico['temporary_files_directory'] = os.path.abspath(os.environ['TMPDIR'])
        else:
            dico['temporary_files_directory'] = os.path.abspath('.')
    for el in ['scheduler_time_step', 'retreatment_time_span', 'forecast_time_span', 'forcing_dates_dt_max', 'model_min_time_step']:
        dico[el] = float(dico[el])
        if dico[el] <= 0.:
            raise Exception('parameter %s must be > 0'%el)
    dico['n_ensemble'] = int(dico['n_ensemble'])
    if dico['n_ensemble'] < 1:
        raise Exception('n_ensemble must be >= 1')
    if dico['nprocs'] == None:
        dico['nprocs'] = max([1, cpu_count()-1])
    if 'verbose' not in dico:
        dico['verbose'] = 1
    if dico['verbose'] is None:
        dico['verbose'] = 1
        
    if dico['assimilation_sources'] is None:
        dico['assimilation_sources'] = []
    dico['assimilation_sources'] = list_form(dico['assimilation_sources'])
    
    #mgb executable and paths
    check_dict(dico['mgb'], ['executable', 'static_data_file', 'input_model'], check_none=True, prefix='in mgb:')
    
    subdict = dico['forcing_grid_geo_selection']
    if subdict == 'auto':
        #read mini.gtp file and guess minmax
        dico['forcing_grid_geo_selection'] = get_lonlat_minmax_from_mgb_static_file(dico['mgb']['static_data_file'])
    if subdict is not None:
        check_dict(subdict, ['latmin', 'latmax', 'lonmin', 'lonmax'], check_none=False, prefix='in forcing_grid_geo_selection:')
        for el in ['latmin', 'latmax', 'lonmin', 'lonmax']:
            if subdict[el] is not None:
                subdict[el] = float(subdict[el])
        if subdict['latmin'] is None:
            subdict['latmin'] = -90.
        if subdict['latmax'] is None:
            subdict['latmax'] = 90.
        if (subdict['lonmin'] is None) and (subdict['lonmax'] is None):
            subdict['lonmin'], subdict['lonmax'] = 0.,360.
        elif (subdict['lonmin'] is None) or (subdict['lonmax'] is None):
            raise Exception('in forcing_grid_geo_selection: lonmin and lonmax must either both be None or both filled')
    
    #2nd level
    #init conditions
    check_dict(dico['init_conditions'], ['time_start'], check_none=True, prefix='in init_conditions: ')
    check_dict(dico['init_conditions'], ['hydrological_state_start'], check_none=False, prefix='in init_conditions: ')
    if not hasattr(dico['init_conditions']['time_start'], 'strftime'):
        dico['init_conditions']['time_start'] = datetime.strptime(dico['init_conditions']['time_start'], '%Y-%m-%dT%H:%M:%S.%f')
    
    if dico['exec_time'] is None:
        dico['exec_time'] = datetime.utcnow()
    elif not hasattr(dico['exec_time'], 'strftime'):
        dico['exec_time'] = datetime.strptime(dico['exec_time'], '%Y-%m-%dT%H:%M:%S.%f')
        
    #perturb static_data
    check_dict(dico['perturb_static_data'], ['activate', 'folder_store', 'varying_parameters','type','mode'], check_none=True, prefix='in perturb_static_data: ')
    if not isinstance(dico['perturb_static_data']['activate'], bool):
        raise Exception('perturb_static_data:activate must be a boolean')
    if (dico['n_ensemble'] < 2) and (dico['perturb_static_data']['activate']):
        print('Cannot perturb static data with an sensemble size of 1 => setting perturb_static_data:activate to false')
        dico['perturb_static_data']['activate'] = False
    if dico['perturb_static_data']['activate']:
        dico['perturb_static_data']['varying_parameters'] = check_perturb_static_data_varying_parameters(dico['perturb_static_data']['varying_parameters'])
        var_params_save_path = '%s/varying_parameters.yaml'%dico['perturb_static_data']['folder_store']
        if os.path.exists(var_params_save_path):
            var_params_saved = check_perturb_static_data_varying_parameters(load_yaml(var_params_save_path), comparison_list=dico['perturb_static_data']['varying_parameters'])
            if len(os.listdir(dico['perturb_static_data']['folder_store'])) != dico['n_ensemble']+1:
                print(dico['n_ensemble'])
                print(os.listdir(dico['perturb_static_data']['folder_store']))
                raise Exception('saved ensemble size mismatch with current ensemble')
        if dico['perturb_static_data']['type'] not in ['saltelli', 'normal']:
            raise Exception('type must be in [saltelli, normal], type %s unknown'%dico['perturb_static_data']['type'])
        if dico['perturb_static_data']['mode'] not in ['per_cell', 'per_variable']:
            raise Exception('mode must be in [per_cell, per_variable], mode %s unknown'%dico['perturb_static_data']['mode'])   
            
    #forcing source
    assert dico['forcing_source'] in ['gsmap', 'era5']
            
    #post_processing
    check_dict(dico['post_processing'], ['science_file', 'portal_file', 'variables'], check_none=True, prefix='in post_processing: ')
    if not isinstance(dico['post_processing']['variables'], list):
        assert isinstance(dico['post_processing']['variables'], str)
        dico['post_processing']['variables'] = [dico['post_processing']['variables']]
    
    return dico
    
    
    
    
def check_perturb_static_data_varying_parameters(list_in, comparison_list=None):
    
    list_in = list_form(list_in)
    if len(list_in) < 1:
        raise Exception('there must be at least one varying parameter')
    for elem in list_in:
        check_dict(elem, ['name', 'error'], check_none=True, prefix='in perturb_static_data:varying_parameters: ')
        elem['error'] = float(elem['error'])
        if elem['error']<0.0:
            raise Exception('Error: %s parameter is negative.'%elem['name'])
            
    if comparison_list is not None:
        dico_in = {el['name']: el['error'] for el in list_in}
        dico_comp = {el['name']: el['error'] for el in comparison_list}
        if dico_in != dico_comp:
            raise Exception('saved ensemble dictionnary mismatch')
        
    return list_in




def get_scheduler_basic_main_parameters_dict():
    dico_input = load_yaml("""
#scheduler main processing input parameters
#main processing = operational analysis and forecasting using MGB-IPH code + assimilation processes
nprocs: 1
verbose: 1

#init conditions (for first init i.e. if hydrological_states_database_directory is empty)
#these conditions are stored into the hydrological states database => if the database exists and the conditions don't match, the program returns an error
init_conditions:
  time_start: '2011-01-01T00:00:00.000000'
  #to start model with an real initial hydrological state => if None than no hydrological state is loaded
  #default configuration is to use the result of a 1200 days simulation starting 01/01/2011 and ending 15/04/2014, this helps initialize the model faster
  #hydrological_state_start: ${hyfaa_run_dir}/data_niger/default_start_hydrostate_niger_20141504.nc
  hydrological_state_start: 
scheduler_time_step: 1.
forecast_time_span: 3.
#choose an ensemble size of 1 if you do not wish to perturb rain inputs
n_ensemble: 1
exec_time: 

#operational_mode: do not activate, it is only to be used when hyfaa scheduler_processing_main is launched by an operational routine and deletes 
#the last retreatment_time_span days to incorporate assimilationn data
operational_mode: false
retreatment_time_span: 1.


#hydrological state database : used as an input and an output to see which date was treated first
hydrological_states_database_directory: ${hyfaa_run_dir}/data_hydro/hydrological_states_db
  
#perturb parameters within static_data file (only those f(i_cells)) : n_ensemble will be generated
perturb_static_data:
  activate: false
  type: normal
  mode: per_variable
  folder_store: ${hyfaa_run_dir}/data_hydro/mini_gtp_ensemble
  varying_parameters:
  -
    name: river_width
    error: 0.3
    min: 2.4
    max: 
  -
    name: river_depth
    error: 0.3
    min: 0.2
    max: 
  -
    name: manning_coefficient
    error: 0.3
    min: 0.03
    max: 0.25
  -
    name: main_river_slope
    error: 0.3  
    min: 0.001
    max: 

#input databases
#forcing source: gsmap, era5, (imerg, arpege : not implemented yet)
forcing_source: gsmap
forcing_grid_geo_selection:
  lonmin: -12.5
  lonmax: 16.5
  latmin: 4.
  latmax: 25.
forcing_grid_database_directory: ${hyfaa_run_dir}/data_niger/forcing_grid_db
forcing_onmesh_database_directory: ${hyfaa_run_dir}/data_niger/forcing_onmesh_db
forcing_dates_dt_max: 1.e-2
rain_uncertainty: 0.5

assimilation_database_directory: ${hyfaa_run_dir}/data_niger/assimilation_db
assimilation_sources:
- ${hyfaa_run_dir}/cmd/hysope_svs.yaml
assim_params_file: ${hyfaa_run_dir}/cmd/assimilation_parameters.yaml



#executable for MGB-IPH simulation
mgb:
  executable: mgb_iph
  static_data_file: ${hyfaa_run_dir}/data_niger/static_data_laetitia.nc
  input_model: ${hyfaa_run_dir}/data_niger/mgb_iph_model_input.yaml
model_min_time_step: 1.
  

#temporary files directory
temporary_files_directory: ${hyfaa_run_dir}/temp



#post-processing
post_processing:
  science_file: ${hyfaa_run_dir}/data_hydro/post_processing_science.nc
  portal_file: ${hyfaa_run_dir}/data_hydro/post_processing_portal.nc
  variables:
    - water_elevation_catchment
    - streamflow_catchment

""", env_vars=True, string_input=True)
    dico_input['assim_params_file'] = get_scheduler_standard_assimilation_dict()
    return dico_input


def get_scheduler_standard_assimilation_dict():
    dico_assim = load_yaml("""
statvars:
  soil_moisture: true
  water_baseflow_reservoir: true
  water_subsurface_reservoir: true
  water_surface_volume: true
  temperature: true
  water_canopy_interception_volume: true
  runoff_generation: true
  streamflow_catchment: true
  water_depth_catchment: true
  previous_water_storage_catchment: true
  flooded_area_catchment: true
  updated_flow_interconnections: true
parameters:
  activate: true
  param_list:
    manning_coefficient: true
    river_depth: true
    river_width: true
    main_river_slope: false
localization: 
  activate: false
  method: velocity
  velocity_file: ${hyfaa_run_dir}/data/mean_velocity.csv
  lengthscale: 10
""", env_vars=False, string_input=True)
    return dico_assim
    
    
    
    
