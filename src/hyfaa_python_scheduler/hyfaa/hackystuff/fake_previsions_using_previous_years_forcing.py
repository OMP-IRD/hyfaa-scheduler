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
from hyfaa.common.yaml.yaml_parser import load_yaml, write_simple_1level_dict_to_yaml_file, write_dict_to_yaml_file
from hyfaa.common.timer.easy_timer import Timer
from hyfaa.check_input_parameters.check_input_parameters import *
from hyfaa.database.forcing.forcing_onmesh_db import ForcingOnMesh_DBManager
from hyfaa.database.hydrostates.hydrostates_db import HydroStates_DBManager, generate_ensemble_files, extract_ensemble_member
from hyfaa.simulation_manager.simulation_manager import SimulationTasks
from hyfaa.hyfaa_postprocessing import hyfaa_postprocessing



def fake_previsions_using_previous_years_forcing(yaml_file_or_dict, ndays, output_dir=None, nyearsmax=None, nprocs=None, verbose=None):
    """main scheduler processing function
    
    :param yaml_file_or_dict: python dict or yaml file
    """
    
    tim = Timer()
    
    print('Starting fake_previsions_using_previous_years_forcing ...')
    
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
    
    
    #make temporary folder and hydro state database folders if they do not exist
    for fol in [dico['hydrological_states_database_directory'], dico['temporary_files_directory']]:
        os.makedirs(fol, exist_ok=True)
    main_temp_dir = tempfile.mkdtemp(dir=dico['temporary_files_directory'], prefix='temphyfaahacky')

    #outputs
    if output_dir is None:
        assert dico['post_processing']['portal_file'] is not None
        output_dir = os.path.join(os.path.dirname(dico['post_processing']['portal_file']), 'prevision_using_previous_years')
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    start_mode = 'initial_solution'
    hydrological_state_start_path = os.path.join(output_dir, 'hydrostate_start.nc')
    output_hydrostate_db_dir = os.path.join(output_dir, 'hydrological_states_db')
    os.makedirs(output_hydrostate_db_dir)
    output_science_file = os.path.join(output_dir, 'post_processing_science.nc')
    output_portal_file = os.path.join(output_dir, 'post_processing_portal.nc')
    
    #get mgb iph input model:
    mgb_iph_input_model = load_yaml(dico['mgb']['input_model'], env_vars=True)
    
    
    #get some dimensions from static_data.nc file
    with netCDF4.Dataset(dico['mgb']['static_data_file']) as ds:
        n_cells = ds.dimensions['n_cells'].size
        n_soil_types = ds.dimensions['n_soil_types'].size
        river_length = ds.variables['longest_river_length'][:]
        mesh_lon = ds.variables['longitude_center'][:]
        mesh_lat = ds.variables['latitude_center'][:]
        
    
    #get file containing last hydrological state in HYFAA hydrological state database
    with HydroStates_DBManager(dico['hydrological_states_database_directory'], mode='r', verbose=verbose) as hydrostates_db_in:
        _, last_date_hydrostate = hydrostates_db_in.get_minmax_dates()
        _, files_match = hydrostates_db_in.get_paths_matching_dates([last_date_hydrostate], dt_max=0.)
    hydrological_state_ens_start_path = files_match[0][0]
    streamflow_ar = np.ma.masked_invalid(np.zeros((dico['n_ensemble'], n_cells), dtype=np.float64))
    with netCDF4.Dataset(hydrological_state_ens_start_path) as ds:
        for ii_ens in range(dico['n_ensemble']):
            streamflow_ar[ii_ens, :] = np.ma.masked_invalid(ds.variables['streamflow_catchment_%d'%ii_ens][:])
    median_streamflow = np.ma.median(streamflow_ar, axis=0)
    ii_ens_selection = np.argmin(np.sum((streamflow_ar-median_streamflow)**2, axis=0))
    extract_ensemble_member(hydrological_state_ens_start_path, hydrological_state_start_path, ii_ens_selection, compression_level=4)
        
    
    #point to the right static data file
    if dico['perturb_static_data']['activate']:
        static_data_file = os.path.join(dico['perturb_static_data']['folder_store'], 'static_data_%d.nc'%ii_ens_select)
        assert os.path.exists(static_data_file), 'static data file %s not found'%static_data_file
    else:
        static_data_file = dico['mgb']['static_data_file']

    #dates to compute
    dates_compute = [last_date_hydrostate + timedelta(ii) for ii in range(1,ndays+1)]

            

    #open forcing on_mesh DB and get years over which fake previsions can be computed
    forcing_db = ForcingOnMesh_DBManager(dico['forcing_onmesh_database_directory'], mode='r', verbose=verbose)
    dates_forcing_available, _ = forcing_db.get_dates()
    nyears_use = int(np.floor((dates_compute[0] - min(dates_forcing_available)).total_seconds()/(3600.*24.*365.25)))
    if nyears_use < 1:
        raise Exception('0 previous years can be used for prevision, exiting...')
    if nyearsmax:
        nyears_use = min(nyears_use, nyearsmax)
    years_use_ar = dates_compute[0].year-np.arange(1,nyears_use+1, dtype=np.int32)
    print('  -> Will compute previsions using forcing data from %d previous years'%nyears_use)
            
        

    
    #open hydrological state database and simulation task manager
    with HydroStates_DBManager(output_hydrostate_db_dir, mode='w', verbose=verbose) as hydrostates_db, \
        SimulationTasks(dico['nprocs'], dico['mgb']['executable'], verbose=0) as sim_tasks:
            
            
        #split ensemble files
        last_hydrological_state_ensemble_paths, temp_input_files_created = generate_ensemble_files(hydrological_state_start_path, nyears_use, main_temp_dir)
					
        
        ############################
        #LOOP ON SCHEDULER TIME STEPS
        if verbose >= 1:
            print('HYFAA initialized: %s'%tim)
        
        last_date = last_date_hydrostate
        for i_compute, date_compute in enumerate(dates_compute):
            tim.level_up()
            common_temp_dir = tempfile.mkdtemp(dir=main_temp_dir, prefix='temphyfaa_%s'%date2strtag(date_compute))
            
            print('Computing hydrological state %s -> %s (step %d/%d):'%(date2str(last_date), date2str(date_compute), i_compute+1, len(dates_compute)))
            

            if verbose >= 1:
                print('  0) Writing MGB parameter files and dynamic forcing files')
            #configure common (date) parameters            
            param_dict = copy.deepcopy(mgb_iph_input_model)
            param_dict['forcing_dates_dt_max'] = 0.
            param_dict['day'] = last_date.day
            param_dict['month'] = last_date.month
            param_dict['year'] = last_date.year
            param_dict['hour'] = last_date.hour
            dt_full = (date_compute - last_date).total_seconds()/(24.*3600.)
            param_dict['nt'] = int(np.ceil(dt_full/dico['model_min_time_step']))
            dt_days = dt_full / (1.*param_dict['nt'])
            param_dict['dt'] = dt_days*24.*3600.
            
            #forcing dates
            dates_forcing_in = [last_date+timedelta(ii*dt_days) for ii in range(param_dict['nt'])]
            overide_dates_forcing = [datetime(param_dict['year'], param_dict['month'], param_dict['day'], param_dict['hour'])+timedelta(ii*dt_days) for ii in range(param_dict['nt'])]
            

            #forcing data from previous years
            for ii_year in range(nyears_use):
                dt_years = timedelta(int(np.floor(365.25*(ii_year+1))))
                if ii_year == 0:
                    _, rain_data = forcing_db.build_forcing_data([el-dt_years for el in dates_forcing_in], data_type='rain', search_delta_before=dico['forcing_dates_dt_max'])
                    perturbed_rain_vectors = np.empty((nyears_use, n_cells, len(dates_forcing_in)), dtype=rain_data.dtype)
                    perturbed_rain_vectors[ii_year,:,:] = rain_data
                else:
                    _, perturbed_rain_vectors[ii_year,:,:] = forcing_db.build_forcing_data([el-dt_years for el in dates_forcing_in], data_type='rain', search_delta_before=dico['forcing_dates_dt_max'])


                
            #############################################
            #BUILD MGB JOBS
            if verbose >= 1:
                print('  2) Preparing %d MGB simulations...'%(nyears_use))
            mgb_output_files_ordered = []
            job_params = []
            

            for i_ensemble in range(nyears_use):
                #configure
                local_temp_dir = '%s/ENS%d'%(common_temp_dir, i_ensemble)
                os.system('mkdir -p %s/output'%local_temp_dir)
                param_dict_loc = copy.deepcopy(param_dict)
                if dico['perturb_static_data']['activate']:
                    param_dict_loc['static_data_file'] = static_data_ensemble_files[i_ensemble]
                else:
                    param_dict_loc['static_data_file'] = dico['mgb']['static_data_file']
                param_dict_loc['output_directory'] = '%s/output/'%local_temp_dir
                param_dict_loc['hydrological_state_read_file'] = last_hydrological_state_ensemble_paths[i_ensemble]
                param_dict_loc['hydrological_state_write_file'] = os.path.join(main_temp_dir, 'hydrostate_out_%d.nc'%i_ensemble)
                mgb_output_files_ordered.append(param_dict_loc['hydrological_state_write_file'])
                param_dict_loc['forcing_file'] = '%s/dynamic_forcing.nc'%local_temp_dir
                #write param file
                param_file = '%s/mgb_iph_input.yaml'%local_temp_dir
                
                write_simple_1level_dict_to_yaml_file(param_dict_loc, param_file)
                
                #make precipitation file
                forcing_db.write_forcing_file(dates_forcing_in, param_dict_loc['forcing_file'], data_type='rain', \
                    forcing_data=perturbed_rain_vectors[i_ensemble,:,:], overide_dates=overide_dates_forcing)
                
                job_params.append([param_file])
            if verbose >= 1:
                print('    => MGB simulations prepared successfuly: %s'%tim)
            
            
   
            
            
            ############################
            #LAUNCH MGB SIMULATION
            if verbose >= 1:
                print('  3) Launching %d MGB simulations...'%(len(job_params)))
                
            job_output_dict = sim_tasks.run(job_params) #MGB simulation !!!
            
            missing_output_files = [el for el in mgb_output_files_ordered if not os.path.exists(el)]
            if len(missing_output_files) > 0:
                print('    => %d MGB simulations failed: %s'%(len(missing_output_files), tim))
                print('MGB calculations failed, simulation output files missing:\n%s'%('\n'.join([' - %s'%el for el in missing_output_files])))
                error_log_file = os.path.join(dico['temporary_files_directory'],'error_jobs_yaml.log')
                print('Dumping job dict to %s'%error_log_file)
                write_dict_to_yaml_file(job_output_dict, error_log_file)
                sys.exit(1)
            if verbose >= 1:
                print('    => MGB simulations ended successfuly: %s'%tim)
            ########################
            
            
                    
                
            #make sure that all ensemble files in the process of being joined where indeed joined before removing them
            hydrostates_db.update_joined_ensemble_files_from_worker()
            if temp_input_files_created:
                for filepath in last_hydrological_state_ensemble_paths:
                    os.unlink(filepath)
            #replace output file paths by input file paths for next iteration
            last_hydrological_state_ensemble_paths = []
            for i_ensemble, filepath in enumerate(mgb_output_files_ordered):
                filepath_new = os.path.join(main_temp_dir, 'hydrostate_in_%d.nc'%i_ensemble)
                shutil.move(filepath, filepath_new)
                last_hydrological_state_ensemble_paths.append(filepath_new)
            temp_input_files_created = True
            ############################

                

            if verbose >= 1:
                print('  4) Adding last_rain_data_loaded vector to MGB simulation outputs and adding resulting hydrological state to database...')
            
            #add last_rain_data_loaded vector to MGB simulation outputs (and compute inertial time step condition to raise weird cases)
            for i_ensemble, filepath in enumerate(last_hydrological_state_ensemble_paths):
                with netCDF4.Dataset(filepath, mode='a') as ds:
                    var_data = ds.createVariable('last_rain_data_loaded', 'f4', ('n_meshes', ), zlib=True, complevel=4, shuffle=True)
                    var_data[:] = perturbed_rain_vectors[i_ensemble,:,-1].astype(np.float32)
            
            #add to database
            hydrostates_db.add({'date_data': date_compute, 'forcing_confidence_coefficient': 1.0}, last_hydrological_state_ensemble_paths)
            
            
            last_date = date_compute
            if verbose >= 1:
                print('    => %s'%tim)

                
            shutil.rmtree(common_temp_dir)
            tim_dico = tim.get_full_info()
            tim.level_down()
            print('  => step complete : %.2f seconds (%.2f from beginning at %s)'%(tim_dico['level'], tim_dico['start'], tim_dico['start_date']))
            
    forcing_db._close_()
    sim_tasks.close()
    
    
    #post-processing
    dico_post = copy.deepcopy(dico)
    dico_post['n_ensemble'] = nyears_use
    dico_post['hydrological_states_database_directory'] = output_hydrostate_db_dir
    dico_post['post_processing']['science_file'] = output_science_file
    dico_post['post_processing']['portal_file'] = output_portal_file
    hyfaa_postprocessing(dico_post)
    for ncfile in [output_science_file, output_portal_file]:
        with netCDF4.Dataset(ncfile, mode='a') as ds:
            ds.createVariable('years_used', np.int32, ['n_ensemble'], zlib=True, complevel=4, shuffle=True)
            ds.variables['years_used'][:] = years_use_ar
            ds.setncattr('fake_prevision_info', "This file was generated using previous years's forcing data, constituting an ensemble.")

    
    
if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="HYFAA hackystuff : make fake previsions based on standard MGB simulation using previous years forcing", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--input_yaml_file", type=str, required=True, help="path to HYFAA input yaml file") 
    parser.add_argument("--ndays", type=int, required=True, help="simulation over ndays")
    parser.add_argument("--output_dir", type=str, help="output directory path. By default it will be in a 'prevision_using_previous_years' folder in the same directory as the portal file. WARNING: will delete contents.")
    parser.add_argument("--verbose", type=int, help="verbose level overload, default to the one contained in the yaml file")
    parser.add_argument("--nyearsmax", type=int, help="maximum number of years to use for fake forcing prevision (will use most recent)")
    parser.add_argument("--nprocs", type=int, help="nprocs overload, default to the one contained in the yaml file")
    args = parser.parse_args()

    assert args.input_yaml_file is not None
    fake_previsions_using_previous_years_forcing(args.input_yaml_file, args.ndays, output_dir=args.output_dir, nyearsmax=args.nyearsmax, nprocs=args.nprocs, verbose=args.verbose)
    

