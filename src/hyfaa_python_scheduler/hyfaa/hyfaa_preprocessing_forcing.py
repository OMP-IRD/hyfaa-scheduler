#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

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
from multiprocessing import cpu_count
import glob

from hyfaa.common.test.common_test_functions import *
from hyfaa.check_input_parameters.check_input_parameters import *
from hyfaa.common.yaml.yaml_parser import load_yaml
from hyfaa.common.parallel_processing.easy_parallel import simple_parallel_run, Job
from hyfaa.database.forcing.forcing_grid_db import ForcingGrid_DBManager
from hyfaa.database.forcing.forcing_onmesh_db import ForcingOnMesh_DBManager
from prometheus_client import CollectorRegistry, Gauge, Summary, write_to_textfile, push_to_gateway

prometheus_registry = CollectorRegistry()

prom_job_duration = Summary("job_duration_seconds", "Duration of a job run",
                           registry=prometheus_registry, labelnames=["app", "instance"])
prom_job_last_success = Gauge(
            "job_last_success_unixtime",
            "Last time a batch job successfully finished",
            registry=prometheus_registry,
            labelnames=["app", "instance"]
        )

prom_mgb_files_most_recent = Gauge(
    "mgb_files_most_recent_unixtime",
    "Most recent file generated during the hyfaa processes. Type can be forcing, onmesh, hydrostate. Source can be gsmap, era5, ecmwf, mgb",
    registry=prometheus_registry,
    labelnames=["app", "instance", "type", "source"]
)

def get_mesh_cell_centers_from_static_data_file(static_data_file):
    with netCDF4.Dataset(static_data_file) as ds:
        longitudes = ds.variables['longitude_center'][:]
        latitudes = ds.variables['latitude_center'][:]
    return longitudes, latitudes


@prom_job_duration.labels(app="hyfaa", instance="guyane").time()
def hyfaa_preprocessing_forcing(yaml_file_or_dict, gsmap_folder_local=None, verbose=None):
    """main scheduler processing function
    
    :param yaml_file_or_dict: python dict or yaml file
    """
    
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
        verbose = 1
        
        
    #load module for forcing data download
    if dico['forcing_source'] == 'gsmap':
        from hyfaa.database.forcing.gsmap.gsmap_download_and_interpolate_module import retrieve_forcing_data, interpolate_forcing_data
    elif dico['forcing_source'] == 'era5':
        from hyfaa.database.forcing.era5.era5_download_and_interpolate_module import retrieve_forcing_data, interpolate_forcing_data
    elif dico['forcing_source'] in ['imerg', 'arpege']:
        raise NotImplementedError('retrieval from forcing data source %s has not been implemented yet'%dico['forcing_source'])
    else:
        raise Exception('forcing data source %s unknown'%dico['forcing_source'])
    
    mesh_lon, mesh_lat = get_mesh_cell_centers_from_static_data_file(dico['mgb']['static_data_file'])
    
    if verbose >= 2:
        print('Forcing pre-proc: Parameters and mesh coordinates loaded...')

    
    #make temporary folder and hydro state database folders if they do not exist
    for fol in [dico['forcing_grid_database_directory'], dico['forcing_onmesh_database_directory'], dico['temporary_files_directory']]:
        os.makedirs(fol, exist_ok=True)
    main_temp_dir = tempfile.mkdtemp(dir=dico['temporary_files_directory'], prefix='temphyfaaprefor')
    
    min_time = dico['init_conditions']['time_start']
    min_day = datetime(min_time.year, min_time.month, min_time.day)
    max_time = dico['exec_time'] + timedelta(dico['forecast_time_span']+1.)
    max_day = datetime(max_time.year, max_time.month, max_time.day) + timedelta(1)
    ndays = (max_day-min_day).days
    expected_forcing_dates = min_time + np.array([timedelta(day_loc) for day_loc in range(ndays+1)])
    
    if verbose >= 2:
        print('Forcing pre-proc: Opening forcing grid database...')

    #update forcing grid DB
    with ForcingGrid_DBManager(dico['forcing_grid_database_directory'], mode='w', verbose=0) as gr_db:
        
        gr_db_dates, _ = gr_db.get_dates(date_min=expected_forcing_dates[0], date_max=expected_forcing_dates[-1])
        missing_dates = set(expected_forcing_dates) - set(gr_db_dates)
        
        if verbose >= 1:
            print('Forcing pre-proc: Downloading data for forcing grid database')
        for file_info in retrieve_forcing_data(missing_dates, geo_selection=dico['forcing_grid_geo_selection'], gsmap_folder_local=gsmap_folder_local, \
            nprocs=dico['nprocs'], temp_folder_base=main_temp_dir, verbose=verbose):
            file_info['data_type'] = 'rain'
            gr_db.add(file_info)
            if verbose >= 1:
                print('Forcing pre-processing: Added date %s to grid DB'%file_info['date_data'].strftime('%Y-%m-%dT%H:%M:%S'))
        gr_db_dates, gr_files_info = gr_db.get_dates(date_min=expected_forcing_dates[0], date_max=expected_forcing_dates[-1])
        gr_db_data_path = gr_db.get_data_path()
        # Set metric
        if gr_db_dates:
            latest_time_seconds = gr_db_dates[-1].timestamp()
            prom_mgb_files_most_recent.labels(app="hyfaa", instance="guyane", type="forcing",
                                              source=dico['forcing_source']).set(latest_time_seconds)
    
    
    if verbose >= 2:
        print('Forcing pre-proc: Opening forcing on-mesh database...')

    #update forcing mesh DB
    with ForcingOnMesh_DBManager(dico['forcing_onmesh_database_directory'], mode='w', verbose=0) as mesh_db:
        
        mesh_db_dates, mesh_db_info = mesh_db.get_dates(date_min=expected_forcing_dates[0], date_max=expected_forcing_dates[-1])
        msh_db_dates_set = set(mesh_db_dates)
        mesh_db_dates = np.array(mesh_db_dates)
        gr_db_dates = np.array(gr_db_dates)
        
        if verbose >= 1:
            print('Forcing pre-proc: Analysing tasks to complete on-mesh DB from grid DB...')
        
        
        files_info_interp_mesh = []
        for ii in range(len(gr_db_dates)):
            file_info_loc = {'file_path': os.path.join(gr_db_data_path, gr_files_info['file_path'][ii]), 'data_type': 'rain', 'date_data': gr_db_dates[ii], \
                    'product_type': gr_files_info['product_type'][ii], 'grid_status': gr_files_info['grid_status'][ii]}
            if gr_db_dates[ii] not in msh_db_dates_set:
                files_info_interp_mesh.append(file_info_loc)
                continue
            ids_match = np.where(mesh_db_dates==gr_db_dates[ii])[0]
            assert len(ids_match) >= 1, 'date not found when it should have been found'
            exact_match = False
            for i0 in ids_match:
                if all([file_info_loc[el] == mesh_db_info[el][i0] for el in ['data_type', 'product_type', 'grid_status']]):
                    exact_match = True
                    break
            if not exact_match:
                files_info_interp_mesh.append(file_info_loc)
                
        if verbose >= 1:
            print('Forcing pre-proc: Interpolating %d grid files to mesh files...'%(len(files_info_interp_mesh)))
    
        for file_info in interpolate_forcing_data(files_info_interp_mesh, mesh_lon, mesh_lat, nprocs=dico['nprocs'], \
            temp_folder_base=main_temp_dir, verbose=verbose):
            
            mesh_db.add(file_info)
            if verbose >= 1:
                print('Forcing pre-proc: Added file %s to on-mesh DB'%file_info['file_path'])
        # Set metric
        if files_info_interp_mesh:
            latest_time_seconds = sorted([ i.get("date_data") for i in files_info_interp_mesh])[-1].timestamp()
            prom_mgb_files_most_recent.labels(app="hyfaa", instance="guyane", type="onmesh",
                                              source=dico['forcing_source']).set(latest_time_seconds)

    
    
    

    
    
def test_main_program():
    
    test_use_web = False
    if 'test_use_web' in os.environ:
        if os.environ['test_use_web'].lower() == 'true':
            test_use_web = internet_available()
    
    if test_use_web:
        print_utest_message('Test full retrieval process for forcing data: success if runs without errors')
        
        main_test_dir = 'temp_%s'%(date2strtag(datetime.now()))
        os.makedirs(main_test_dir, exist_ok=True)
        
        dico_input = get_scheduler_basic_main_parameters_dict()

        dico_input['init_conditions']['time_start'] = datetime(2000,1,1)
        dico_input['exec_time'] = datetime(2000,1,2)
        dico_input['hydrological_states_database_directory'] = '%s/hydrological_states_database_directory'%main_test_dir
        dico_input['forcing_source'] = 'gsmap'
        dico_input['forcing_grid_geo_selection'] = {'lonmin': -12., 'lonmax': 16., 'latmin': 4., 'latmax': 25.}
        dico_input['forcing_grid_database_directory'] = '%s/forcing_grid_database_directory'%main_test_dir
        dico_input['forcing_onmesh_database_directory'] = '%s/forcing_onmesh_database_directory'%main_test_dir
        dico_input['assimilation_database_directory'] = '%s/assimilation_database_directory'%main_test_dir
        dico_input['assimilation_database_directory'] = '%s/assimilation_database_directory'%main_test_dir
        dico_input['assimilation_sources'] = [{'function': 'get_hysope_data', \
            'function_info': [{'sv_name': 'L_baikal', 'mesh_id': 10}, {'sv_name': 'R_amz_nap_env_0422_01', 'mesh_id': 12}]}]
        dico_input['temporary_files_directory'] = '%s/temporary_files_directory'%main_test_dir
        dico_input['output_directory']['directory'] = '%s/output_directory'%main_test_dir
        dico_input['assim_params_file'] = 'None'

        try:
            scheduler_pre_processing_forcing(dico_input, verbose=2)
            success_message_utest()
        except:
            fail_message_utest()
            
        shutil.rmtree(main_test_dir)


def write_prometheus_metrics(registry, metrics_filepath=None, pushgateway_url=None):
    if pushgateway_url:
        push_to_gateway(pushgateway_url, job='preprocessing_forcing', registry=registry)
    if metrics_filepath:
        write_to_textfile(metrics_filepath, registry)


if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="HYFAA pre-processing forcing chain", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--test", action='store_true', help="launch unit tests")
    parser.add_argument("--input_yaml_file", type=str, help="path to input yaml file")
    parser.add_argument("--gsmap_folder_local", type=str, help="Use a local folder that contains GSMAP data yearmonth folders. This is useful for tests or simulation especially when FTP is not working.")
    parser.add_argument("--verbose", type=int, help="verbose level overload, default to the one contained in the yaml file")
    args = parser.parse_args()

    if args.test:
        test_main_program()
    else:
        assert args.input_yaml_file is not None
        hyfaa_preprocessing_forcing(args.input_yaml_file, gsmap_folder_local=args.gsmap_folder_local, verbose=args.verbose)

    # prom_gauge.set_to_current_time(process="preprocessing_forcing")
    prom_job_last_success.labels(app="hyfaa", instance="guyane").set_to_current_time()
    prom_pushgateway_url = os.getenv("PROM_PUSHGATEWAY_URL")
    prom_metrics_textfile = os.getenv("PROM_METRICS_TEXTFILE")
    write_prometheus_metrics(prometheus_registry,
                             metrics_filepath=prom_metrics_textfile,
                             pushgateway_url=prom_pushgateway_url)
