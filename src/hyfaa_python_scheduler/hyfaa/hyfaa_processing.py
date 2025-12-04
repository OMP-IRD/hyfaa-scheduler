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
from multiprocessing import cpu_count

from hyfaa.common.test.common_test_functions import *
from hyfaa.common.yaml.yaml_parser import load_yaml, write_simple_1level_dict_to_yaml_file, write_dict_to_yaml_file
from hyfaa.common.timer.easy_timer import Timer
from hyfaa.check_input_parameters.check_input_parameters import *
from hyfaa.database.forcing.forcing_onmesh_db import ForcingOnMesh_DBManager
from hyfaa.database.hydrostates.hydrostates_db import HydroStates_DBManager, generate_ensemble_files
from hyfaa.database.assimilation.assimilation_db import Assimilation_Database
from hyfaa.simulation_manager.simulation_manager import SimulationTasks
import hyfaa.assimilation.EnKF_filter as EnKF_filter
import hyfaa.assimilation.assim_tools as assim_tools

from prometheus_client import CollectorRegistry, Gauge, Summary

from hyfaa.utils.monitoring import write_prometheus_metrics

prometheus_registry = CollectorRegistry()


prom_job_duration = Summary("job_duration_seconds", "Duration of a job run",
                           registry=prometheus_registry, labelnames=["app", "instance", "process"])
prom_job_last_success = Gauge(
            "job_last_success_unixtime",
            "Last time a batch job successfully finished",
            registry=prometheus_registry,
            labelnames=["app", "instance", "process"]
        )


def check_assimilation_static_vars_are_in_hydrological_states(last_hydrological_state_ensemble_paths, list_contr_params):
    for ii, hydro_file in enumerate(last_hydrological_state_ensemble_paths):
        with netCDF4.Dataset(hydro_file) as ds:
            for elem in list_contr_params:
                if elem not in ds.variables:
                    #raise Exception('Static variable %s should not be missing from hydrological state file %s'%(elem, hydro_file) + \
                    #    '\n -> Maybe you tried to activate parameter assimilation but launched from a previous calculation where this option was not activated.')
                    return False
                else:
                    return True
                    
                        


def add_assimilation_static_vars_to_hydrological_states(last_hydrological_state_ensemble_paths, static_data_ensemble_files, list_contr_params, verbose=1):
    for ii, (hydro_file, param_file) in enumerate(zip(last_hydrological_state_ensemble_paths, static_data_ensemble_files)):
        with netCDF4.Dataset(hydro_file, mode='a') as ds_w, netCDF4.Dataset(param_file) as ds_r:
            for elem in list_contr_params:
                if elem in ds_w.variables:
                    raise Exception('variable %s should not already be in file %s'%(elem, hydro_file))
                assert elem in ds_r.variables, 'Variable %s is not contained in parameter file' % elem
                # check if dimension exist
                for dim in ds_r.variables[elem].dimensions:
                    if dim not in ds_w.dimensions:
                        # create dimension
                        ds_w.createDimension('%s' % dim, ds_r.dimensions[dim].size)
                    else:
                        assert ds_r.dimensions[dim].size == ds_w.dimensions[dim].size, 'dimensions mismatch between static file %s and hydrological state file %s'%(param_file, hydro_file)
                dspar = ds_w.createVariable(elem, ds_r.variables[elem].dtype, [el for el in ds_r.variables[elem].dimensions])
                # upload data to new netCDF file
                dspar[:] = ds_r.variables[elem][:]
        if verbose > 0:
            print('Loaded static data from param file %s\n -> written to hydrological state file %s'%(param_file, hydro_file))


@prom_job_duration.labels(app="hyfaa", instance="guyane", process="processing").time()
def hyfaa_processing(yaml_file_or_dict, verbose=None):
    """main scheduler processing function
    
    :param yaml_file_or_dict: python dict or yaml file
    """
    
    tim = Timer()
    
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
        
    #load assimilation parameters
    if isinstance(dico['assim_params_file'], dict):
        dico_assim = dico['assim_params_file']
    elif os.path.exists(dico['assim_params_file']):
        dico_assim = load_yaml(dico['assim_params_file'], env_vars=True)
    else:
        raise Exception('There is no file for DA parameters')
    #be sure to remove parameters from control list if their perturbation is not activated
    if not dico['perturb_static_data']['activate']:
        dico_assim['parameters']['activate'] = False
    
    
    #make temporary folder and hydro state database folders if they do not exist
    for fol in [dico['hydrological_states_database_directory'], dico['temporary_files_directory']]:
        os.makedirs(fol, exist_ok=True)
    main_temp_dir = tempfile.mkdtemp(dir=dico['temporary_files_directory'], prefix='temphyfaamain')

    #rain_uncertainty
    prec_error=dico['rain_uncertainty']
    
    #get mgb iph input model:
    mgb_iph_input_model = load_yaml(dico['mgb']['input_model'], env_vars=True)
    alpha_inertial = float(mgb_iph_input_model['alpha'])
    
    
    #get some dimensions from static_data.nc file
    with netCDF4.Dataset(dico['mgb']['static_data_file']) as ds:
        n_cells = ds.dimensions['n_cells'].size
        n_soil_types = ds.dimensions['n_soil_types'].size
        river_length = ds.variables['longest_river_length'][:]
        mesh_lon = ds.variables['longitude_center'][:]
        mesh_lat = ds.variables['latitude_center'][:]
        
    
    if dico['n_ensemble'] < 2:
        print('Ensemble size of 1 => setting to false all options to perturb parameters or input forcing, as well as assimilation.')
        dico['activate_assimilation'] = False
        dico['perturb_static_data']['activate'] = False
        dico_assim['parameters']['activate'] = False
    if dico['perturb_static_data']['activate']:
        static_data_ensemble_files = ['%s/static_data_%d.nc'%(dico['perturb_static_data']['folder_store'], ii) for ii in range(dico['n_ensemble'])]
        files_present = [el for el in static_data_ensemble_files if os.path.exists(el)]
        if len(files_present) == 0:
            print('Generating perturbed static_data files ...')
            assim_tools.generate_perturbed_static_data_files(dico['mgb']['static_data_file'], dico['perturb_static_data']['folder_store'], \
                dico['perturb_static_data']['varying_parameters'], dico['n_ensemble'],dico['perturb_static_data']['type'],dico['perturb_static_data']['mode'])
        elif len(files_present) != dico['n_ensemble']:
            raise Exception('some ensemble files missing')
    else:
        static_data_ensemble_files = [dico['mgb']['static_data_file']]

    #get list of control parameters
    list_contr_params=[]
    if dico_assim['parameters']['activate']:
        print('control parameters exist and will be added to hydrological states files')
        for elem, activated in dico_assim['parameters']['param_list'].items():
            if activated:
                list_contr_params.append(elem)
    if len(list_contr_params)==0:
        print('No parameters are listed, control parameters option is deactivated')
            

    #open forcing om-mesh database
    forcing_db = ForcingOnMesh_DBManager(dico['forcing_onmesh_database_directory'], mode='r', verbose=verbose)
    
    ##create assim filter object
    Assimilation_filter = EnKF_filter.EnKFobj(dico_assim, dico['n_ensemble'], n_cells, n_soil_types, mesh_lon, mesh_lat, dico['model_min_time_step'])
    
    #open assimilation database
    Assim_db = Assimilation_Database(dico['assimilation_database_directory'], mode='r', verbose=verbose)
    #open hydrological state database and simulation task manager
    with HydroStates_DBManager(dico['hydrological_states_database_directory'], mode='w', verbose=verbose) as hydrostates_db, \
        SimulationTasks(dico['nprocs'], dico['mgb']['executable'], verbose=0) as sim_tasks:
            
        if dico['operational_mode']:
            #remove all dates after exec_time - retreatment_time_span => so it takes into account new assimilation data
            hydrostates_db.remove_after_date(exec_time -timedelta(dico['retreatment_time_span']))
         #get last date (if it exists)
        _, last_date_hydrostate = hydrostates_db.get_minmax_dates()
        
        
        ###############
        #check start modes : if it restarts from an existing calculation, if it is initialized from an initial solution or nothing
        if last_date_hydrostate is None:
            start_mode = 'null'
            last_date_hydrostate = dico['init_conditions']['time_start']
            hydrological_state_start_path = dico['init_conditions']['hydrological_state_start']
            if hydrological_state_start_path is not None:
                start_mode = 'initial_solution'
        else:
            start_mode = 'restart'
            _, files_match = hydrostates_db.get_paths_matching_dates([last_date_hydrostate], dt_max=0., type_request='control')
            files_match = files_match[0]
            #pick last date
            hydrological_state_start_path = files_match[0]
            #add parameters to hydrological states
            
        #split ensemble files
        last_hydrological_state_ensemble_paths, temp_input_files_created = generate_ensemble_files(hydrological_state_start_path, dico['n_ensemble'], main_temp_dir)
        
        #if last_hydrological_state_ensemble_paths contains real file paths and 
        if dico['n_ensemble'] > 1 and len(list_contr_params) > 0:
            if start_mode == 'initial_solution':
                add_assimilation_static_vars_to_hydrological_states(last_hydrological_state_ensemble_paths, static_data_ensemble_files, list_contr_params, verbose=1)
            elif start_mode == 'restart':
                condinit_param=check_assimilation_static_vars_are_in_hydrological_states(last_hydrological_state_ensemble_paths, list_contr_params)
                if not condinit_param:
                    add_assimilation_static_vars_to_hydrological_states(last_hydrological_state_ensemble_paths, static_data_ensemble_files, list_contr_params, verbose=1)
					

        #dates to compute
        #from beginning
        dates_compute = np.arange(datetime_to_julianday(dico['init_conditions']['time_start']), \
            datetime_to_julianday(exec_time) + dico['forecast_time_span']+dico['scheduler_time_step'], dico['scheduler_time_step'])
        #select only those > last_date_hydrostate + model_min_time_step/3
        dates_compute = dates_compute[dates_compute > datetime_to_julianday(last_date_hydrostate)+dico['model_min_time_step']/3.]
        #select only those <= last_date_forcing + model_min_time_step
        dates_forcing_available, _ = forcing_db.get_dates()
        dates_forcing_available = set(dates_forcing_available)
        last_date_forcing = max(dates_forcing_available)
        dates_compute = dates_compute[dates_compute <= datetime_to_julianday(last_date_forcing)+dico['model_min_time_step']]
        #convert to datetimes
        dates_compute = [julianday_to_datetime(el) for el in dates_compute]
        last_date = last_date_hydrostate

        ############################
        #LOOP ON SCHEDULER TIME STEPS
        if verbose >= 1:
            print('HYFAA initialized: %s'%tim)
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
            

            #forcing perturbations
            _, rain_data = forcing_db.build_forcing_data(dates_forcing_in, data_type='rain', search_delta_before=dico['forcing_dates_dt_max'])
            
            #############################################
            #get fields of multiplication factors
            perturbed_rain_vectors = np.empty((dico['n_ensemble'], np.shape(rain_data)[0], np.shape(rain_data)[1]), dtype=rain_data.dtype)
            
            if dico['n_ensemble'] == 1:
                print('    => ensemble size = 1, no rain perturbation')
                for ii in range(dico['n_ensemble']):
                    perturbed_rain_vectors[ii,:,:] = rain_data
            else:
                perturbed_multfact_vector = assim_tools.build_gaussian_error_fields(n_cells, dico['n_ensemble'], dico['forcing_grid_geo_selection'], mesh_lon, mesh_lat)
                #build fields of perturbed precipitation
                
                for ii in range(np.shape(rain_data)[1]):
                    perturbed_rain_vectors[:,:,ii] = assim_tools.rain_perturbation(perturbed_multfact_vector, rain_data[:,ii],prec_error)
                    
                    #perturbed_rain_vectors[:,:,ii] =
            perturbed_rain_vectors[perturbed_rain_vectors<0.] = 0.
            #############################################
            


            #############################################
            #perform analysis
            if verbose >= 1:
                print('    => %s'%tim)
                print('  1) Analysis step')
            if dico['activate_assimilation']:
                Assim_obs = Assim_db.get_values_between_dates(date_start=last_date, date_end=date_compute, dt_max=0.0, start_strict=False, end_strict=True)
                if len(Assim_obs) > 0 and (start_mode == 'start_mode' or (start_mode != 'start_mode' and i_compute>0)):
                    if verbose >= 2:
                        print('   => Observations data avalaible ===> Performing analysis...')
                    #make sure that all ensemble files in the process of being joined where indeed joined before the assimilation process modifies them
                    hydrostates_db.update_joined_ensemble_files_from_worker()
                    #perform analysis
                    Assimilation_filter.perform_analysis(Assim_obs, last_hydrological_state_ensemble_paths, inplace=True, timestep_analysis=None)
                    #add corrected hydro states to DB
                    hydrostates_db.add({'date_data': last_date, 'forcing_confidence_coefficient': 1.0, 'number_obs_used': len(Assim_obs),'type': 'analysis'}, \
                        last_hydrological_state_ensemble_paths)
                    if verbose >= 1:
                        print('    => Assimilated data: %s'%tim)
                else:
                    if verbose >= 1:
                        print('    => No data assimilated: %s'%tim)
            else:
                if verbose >= 1:
                    print('    => Assimilation not activated, no data assimilation: %s'%tim)
                
                
                
                
            #############################################
            #BUILD MGB JOBS
            if verbose >= 1:
                print('  2) Preparing %d MGB simulations...'%(dico['n_ensemble']))
            mgb_output_files_ordered = []
            job_params = []
            

            for i_ensemble in range(dico['n_ensemble']):
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
            # ~ print(job_output_dict[0]['output'])
            
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
            
            
            
            
            
            
            
            
                
                
            #in the case where the scheduler was started with no initial hydrological state file and static parameters are marked up for assimilation, after the first MGB job, 
            #hydrological state files do not contain static parameters, so we add them here
            if dico['n_ensemble'] > 1 and len(list_contr_params) > 0 and i_compute==0 and start_mode == 'null':
                add_assimilation_static_vars_to_hydrological_states(mgb_output_files_ordered, static_data_ensemble_files, list_contr_params, verbose=1)
                    
                
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
    Assim_db._close_()
    sim_tasks.close()
    
    if temp_input_files_created:
        for filepath in last_hydrological_state_ensemble_paths:
            os.unlink(filepath)
    
    
    
def test_main_program(main_test_dir=None):
    
    print_utest_message('Test main scheduler process: success if runs without errors')
        
    main_test_dir = tempfile.mkdtemp(prefix='temp_spm_', dir='.')
    
    dico_input = get_scheduler_basic_main_parameters_dict()
      
    dico_input['init_conditions']['time_start'] = datetime(2000,1,1)
    dico_input['exec_time'] = datetime(2000,1,3)
    dico_input['hydrological_states_database_directory'] = '%s/hydrological_states_database_directory'%main_test_dir
    dico_input['forcing_grid_database_directory'] = '%s/forcing_grid_database_directory'%main_test_dir
    dico_input['forcing_onmesh_database_directory'] = '%s/forcing_onmesh_database_directory'%main_test_dir
    dico_input['assimilation_database_directory'] = '%s/assimilation_database_directory'%main_test_dir
    dico_input['temporary_files_directory'] = '%s/temporary_files_directory'%main_test_dir
    dico_input['perturb_static_data']['folder_store'] = '%s/perturbed_static_data_store'%main_test_dir
    
    with netCDF4.Dataset(dico_input['mgb']['static_data_file']) as ds:
        n_cells = ds.dimensions['n_cells'].size
    
    ntimes = int(np.ceil((dico_input['exec_time'] + timedelta(dico_input['forecast_time_span']) - dico_input['init_conditions']['time_start']).total_seconds() / \
        (24.*3600.*dico_input['scheduler_time_step']))) + 1
    model_times = dico_input['init_conditions']['time_start']+np.array([timedelta(ii*dico_input['scheduler_time_step']) for ii in range(ntimes)])
    
    #create test forcing database
    os.makedirs(dico_input['forcing_onmesh_database_directory'])
    with ForcingOnMesh_DBManager(dico_input['forcing_onmesh_database_directory'], mode='w') as db:
        for ii in range(ntimes):
            product_type = 'analysis'
            if model_times[ii] >= dico_input['exec_time']:
                product_type = 'forecast'
            loc_file = '%s/forcing_%d.nc'%(main_test_dir, ii)
            with netCDF4.Dataset(loc_file, mode='w') as ds:
                ds.createDimension('n_meshes', n_cells)
                var = ds.createVariable('rain', 'f4', ('n_meshes', ), zlib=True, complevel=4, shuffle=True)
                var[:] = np.random.rand(n_cells).astype(np.float32)*10.
                ds.setncattr('date', date2str(model_times[ii]))
                ds.setncattr('date_julianday', datetime_to_julianday(model_times[ii]))
            db.add({'file_path': loc_file, 'data_type': 'rain', 'date_data': model_times[ii], 'product_type': product_type, 'grid_status': 'complete'})
    
    #create test assimilation database
    n_cells_assim_per_day = 10
    meas_info = []
    for icell in range(n_cells_assim_per_day):
        mesh_id, sv_name, lon, lat = np.random.randint(1, n_cells), 'SV_%d'%icell, np.random.rand()*360., np.random.rand()*180.-90.
        for ii in range(ntimes):
            meas_info.append({'mesh_id': mesh_id, 'value': np.random.rand()*10., 'uncertainty': np.random.rand()+1., \
                'date_data': model_times[ii]-timedelta(np.random.rand()*0.5), 'lon': lon, 'lat': lat, \
                'sv_name': sv_name, 'info': ''})
    os.makedirs(dico_input['assimilation_database_directory'])
    with Assimilation_Database(dico_input['assimilation_database_directory'], mode='w') as db:
        db.add(meas_info)

    try:
        hyfaa_processing(dico_input)
        success_message_utest()
    except:
        fail_message_utest()
        
    shutil.rmtree(main_test_dir)
    
    
if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="HYFAA main processing chain", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--test", action='store_true', help="launch unit tests")
    parser.add_argument("--input_yaml_file", type=str, help="path to input yaml file") 
    parser.add_argument("--verbose", type=int, help="verbose level overload, default to the one contained in the yaml file")
    args = parser.parse_args()

    if args.test:
        test_main_program()
    else:
        assert args.input_yaml_file is not None
        hyfaa_processing(args.input_yaml_file, verbose=args.verbose)

    prom_job_last_success.labels(app="hyfaa", instance="guyane", process="processing").set_to_current_time()
    prom_pushgateway_url = os.getenv("PROM_PUSHGATEWAY_URL")
    prom_metrics_textfile = os.getenv("PROM_METRICS_TEXTFILE")
    write_prometheus_metrics(prometheus_registry,
                             metrics_filepath=prom_metrics_textfile,
                             pushgateway_url=prom_pushgateway_url)


