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
import sqlite3
import pandas

from hyfaa.common.test.common_test_functions import *
from hyfaa.database.forcing.forcing_db import Forcing_DBManager



class ForcingOnMesh_DBManager(Forcing_DBManager):
    """This class handles the forcing on-mesh database.
    """
    
    def __init__(self, db_sql_path, mode='r', verbose=1):
        super().__init__(db_sql_path, mode=mode, verbose=verbose)
        self._priority_coefficients = {'product_type': {'analysis': 1.2, 'forecast': 1.0}, \
            'grid_status': {'complete': 1.3, 'incomplete': 0.6, 'replaced': 0.2}}
    
    def __enter__(self):
        return super().__enter__()
    
    def __exit__(self, type, value, traceback):
        return super().__exit__(type, value, traceback)            
            
    def build_forcing_data(self, dates_in, data_type='rain', search_delta_before=None, search_delta_after=None):
        """Builds forcing data n_cells*nt matrix by retrieving data from database matching dates_in
        
        :param data_type: rain, ...
        :param dates_in: list of dates (datetimes)
        :param dt_max: positive float
        """
        
        dates_of_files_chosen, files_chosen = self.get_best_paths_matching_dates_and_data_type(dates_in, data_type=data_type, \
            search_delta_before=search_delta_before, search_delta_after=search_delta_after, error_if_none_found=True)
        
        nt = len(files_chosen)
        
        #open first file to get dimensions
        with netCDF4.Dataset(files_chosen[0]) as ds:
            n_cells = ds.dimensions['n_meshes'].size
            data = np.empty((n_cells, nt), dtype=np.float32)
            data[:,0] = ds.variables[data_type][:]
        for ii in range(1,len(files_chosen)):
            with netCDF4.Dataset(files_chosen[ii]) as ds_in:
                if n_cells != ds_in.dimensions['n_meshes'].size:
                    raise Exception('File %s has different mesh size than file %s'%(files_chosen[0], files_chosen[ii]))
                data[:,ii] = ds_in.variables[data_type][:]
                
        return dates_of_files_chosen, data
        
        
    def write_forcing_file(self, dates_in, file_out, data_type='rain', dt_max=0.01, forcing_data=None, overide_dates=None):
        """Builds forcing file containing n_cells*nt matrix by retrieving data from database matching dates_in
        
        :param data_type: rain, ...
        :param dates_in: list of dates (datetimes)
        :param file_out: output netCDF file
        :param dt_max: positive float
        """
        
        if overide_dates is not None:
            if len(overide_dates) != len(dates_in):
                raise Exception('overide dates (%d) must have the same length as dates_in (%d)'%(len(overide_dates), len(dates_in)))
        
        if forcing_data is None:
            dates_of_files_chosen, forcing_data = self.build_forcing_data(dates_in, data_type=data_type, dt_max=dt_max)
        else:
            dates_of_files_chosen = dates_in
                
        with netCDF4.Dataset(file_out, mode='w') as ds:
            ds.createDimension('n_meshes', np.shape(forcing_data)[0])
            ds.createDimension('nt', len(dates_in))
            var_dates = ds.createVariable('dates', 'f8', ('nt', ), zlib=True, complevel=4, shuffle=True)
            var_dates_original = ds.createVariable('dates_original', 'f8', ('nt', ), zlib=True, complevel=4, shuffle=True)
            var_data = ds.createVariable(data_type, 'f4', ('n_meshes', 'nt'), zlib=True, complevel=4, shuffle=True)
            if overide_dates is None:
                var_dates[:] = np.array([datetime_to_julianday(el) for el in dates_in], dtype=np.float64)
            else:
                var_dates[:] = np.array([datetime_to_julianday(overide_dates[ii]) for ii in range(len(overide_dates))], dtype=np.float64)
            var_dates_original[:] = np.array([datetime_to_julianday(el) for el in dates_of_files_chosen], dtype=np.float64)
            var_data[:] = forcing_data


        
        
    def build_forcing_file_chunked(self, dates_in, file_out, data_type='rain', max_chunksize=1000, dt_max=0.01, overide_dates = None):
        """Builds forcing file containing n_cells*nt matrix by retrieving data from database matching dates_in
        
        :param data_type: rain, ...
        :param dates_in: list of dates (datetimes)
        :param file_out: output netCDF file
        :param dt_max: positive float
        """
        
        if overide_dates is not None:
            if len(overide_dates) != len(dates_in):
                raise Exception('overide dates must have the same length as dates_in')
        
        initialized = False
        i0 = 0
        with netCDF4.Dataset(file_out, mode='w') as ds:
            for ids_chunk in stack_generator(range(len(dates_in)), max_chunksize):
                dates_in_chunk = [dates_in[ii] for ii in ids_chunk]
                nchunk = len(dates_in_chunk)
                dates_of_files_chosen_chunk, var_data_chunk = self.build_forcing_data(dates_in_chunk, data_type=data_type, dt_max=dt_max)
                if not initialized:
                    ds.createDimension('n_meshes', np.shape(var_data_chunk)[0])
                    ds.createDimension('nt', len(dates_in))
                    var_dates = ds.createVariable('dates', 'f8', ('nt', ), zlib=True, complevel=4, shuffle=True)
                    var_dates_original = ds.createVariable('dates_original', 'f8', ('nt', ), zlib=True, complevel=4, shuffle=True)
                    var_data = ds.createVariable(data_type, 'f4', ('n_meshes', 'nt'), zlib=True, complevel=4, shuffle=True)
                    initialized = True
                if overide_dates is None:
                    var_dates[i0:i0+nchunk] = np.array([datetime_to_julianday(el) for el in dates_in_chunk], dtype=np.float64)
                else:
                    var_dates[i0:i0+nchunk] = np.array([datetime_to_julianday(overide_dates[ii]) for ii in ids_chunk], dtype=np.float64)
                var_dates_original[i0:i0+nchunk] = np.array([datetime_to_julianday(el) for el in dates_of_files_chosen_chunk], dtype=np.float64)
                var_data[:,i0:i0+nchunk] = var_data_chunk
                i0 += nchunk
        





def test_db():
    """main test for database object
    """
    
    print_utest_message('\n\nRunning main forcing on-mesh database manager diagnostics:\n')
    
    main_test_dir = 'temp_%s'%(date2strtag(datetime.now()))
    if os.path.exists(main_test_dir):
        shutil.rmtree(main_test_dir)
    sql_db_dir = os.path.join(main_test_dir, 'forcing_database')
    fake_data_dir = os.path.join(main_test_dir, 'fake_data')
        
    #create directories
    for fol in [sql_db_dir, fake_data_dir]:
        os.makedirs(fol, exist_ok=True)
            
    #create fake data
    nt_data = 10
    ncells_data = 100
    dt_data = 1.
    time_files = []
    filenames = []
    for ii in range(nt_data):
        time_files.append(datetime(2011,1,27)+timedelta(dt_data*ii))
        filenames.append('%s/rain_proper_%s'%(fake_data_dir, ForcingOnMesh_DBManager.date2str(time_files[-1])))
        with netCDF4.Dataset(filenames[-1], mode='w') as ds:
            ds.createDimension('n_meshes', ncells_data)
            var = ds.createVariable('rain', 'f4', ('n_meshes', ), zlib=True, complevel=4, shuffle=True)
            var[:] = np.ones(ncells_data, dtype=np.float32)*ii
    time_improper_file = time_files[-1]+timedelta(dt_data*5.)
    improper_file = '%s/rain_improper_%s'%(fake_data_dir, ForcingOnMesh_DBManager.date2str(time_improper_file))
    with netCDF4.Dataset(improper_file, mode='w') as ds:
        ds.createDimension('n_meshes', ncells_data+9)
        var = ds.createVariable('rain', 'f4', ('n_meshes', ), zlib=True, complevel=4, shuffle=True)
        var[:] = np.ones(ncells_data+9, dtype=np.float32)*ii
            
    
    #1 : test opening empty database in read mode
    print_utest_message('Test that opening empty database in read mode fails')
    try:
        with ForcingOnMesh_DBManager(sql_db_dir, mode='r', verbose=0) as db:
            fail_message_utest()
    except:
        success_message_utest()
    check_condition_utest('Test that opening empty database in read mode does not generate file creation', len(os.listdir(sql_db_dir)) == 0)

    
    #2 : test opening empty database in write mode
    print_utest_message('Test that opening empty database in write mode succeeds')
    try:
        with ForcingOnMesh_DBManager(sql_db_dir, mode='w', verbose=0) as db:
            success_message_utest()
            check_condition_utest('Test that opening empty database in write mode creates necessary files', \
                all([os.path.exists(el) for el in ['%s/database_manager.sql'%sql_db_dir, '%s/database_manager.sql_lock'%sql_db_dir, '%s/data_store'%sql_db_dir]]))
    except:
        shutil.rmtree(main_test_dir)
        fail_message_utest()
    check_condition_utest('Test that lock file is removed upon database closure', not os.path.exists('%s/database_manager.sql_lock'%sql_db_dir))
    
    
    
    #3 : check that opening database without context manager succeeds in read mode
    print_utest_message('Test that opening database in read mode without context manager succeeds')
    try:
        db = ForcingOnMesh_DBManager(sql_db_dir, mode='r', verbose=0)
        db._close_()
        success_message_utest()
    except:
        fail_message_utest()
    #4 : check that it fails in write mode
    print_utest_message('Test that opening database in write mode without context manager fails')
    try:
        db = ForcingOnMesh_DBManager(sql_db_dir, mode='w', verbose=0)
        db._close_()
        fail_message_utest()
    except:
        success_message_utest()
    
    
    
    #5 : manually create lock file and test database opening (should be ok in read mode and fail in write mode)
    shutil.copy('%s/database_manager.sql'%sql_db_dir, '%s/database_manager.sql_lock'%sql_db_dir)
    #read mode
    print_utest_message('Test reading while lock is active')
    try:
        with ForcingOnMesh_DBManager(sql_db_dir, mode='r', verbose=0) as db:
            success_message_utest()
    except:
        shutil.rmtree(main_test_dir)
        fail_message_utest()
    #write mode
    print_utest_message('Test if write is prevented while lock is active')
    try:
        with ForcingOnMesh_DBManager(sql_db_dir, mode='w', verbose=0) as db:
            fail_message_utest()
    except:
        success_message_utest()
    #check that lock file still exists after ForcingOnMesh_DBManager exit in write mode
    check_condition_utest('Test that lock file still exists when ForcingOnMesh_DBManager exits because another instance of the manager is already opened in write mode', \
        os.path.exists('%s/database_manager.sql_lock'%sql_db_dir))
    os.unlink('%s/database_manager.sql_lock'%sql_db_dir)
    
        
    #6 test writing
    print_utest_message('Test writing in database')
    try:
        with ForcingOnMesh_DBManager(sql_db_dir, mode='w', verbose=0) as db:
            for ii in range(len(filenames)):
                db.add({'file_path': filenames[ii], 'data_type': 'rain', 'date_data': time_files[ii], 'product_type': 'analysis', 'grid_status': 'complete'})
            db.add({'file_path': improper_file, 'data_type': 'rain', 'date_data': time_improper_file, 'product_type': 'analysis', 'grid_status': 'complete'})
        success_message_utest()
    except:
        shutil.rmtree(main_test_dir)
        fail_message_utest()


    #7 test reading
    db = ForcingOnMesh_DBManager(sql_db_dir, mode='r', verbose=0)
    
    print_utest_message('Test get_best_paths_matching_dates_and_data_type for exact dates')
    times_data, file_paths = db.get_best_paths_matching_dates_and_data_type(time_files, data_type='rain', dt_max=0.)
    success_message_utest()
    check_condition_utest('Test if dates retrieved match dates added', all([abs((times_data[ii]-time_files[ii]).total_seconds())<=dt_date_tolerance_seconds for ii in range(len(time_files))]))
        
    print_utest_message('Test get_best_paths_matching_dates_and_data_type for near dates')
    times_data, file_paths = db.get_best_paths_matching_dates_and_data_type([el+timedelta((np.random.rand(1)[0]-0.5)*0.99) for el in time_files], data_type='rain', dt_max=0.5)
    success_message_utest()
    check_condition_utest('Test if dates retrieved match dates added', all([abs((times_data[ii]-time_files[ii]).total_seconds())<=dt_date_tolerance_seconds for ii in range(len(time_files))]))
        
    print_utest_message('Test get_best_paths_matching_dates_and_data_type for near dates with large dt_max')
    times_data, file_paths = db.get_best_paths_matching_dates_and_data_type([el+timedelta((np.random.rand(1)[0]-0.5)*0.99) for el in time_files], data_type='rain', dt_max=2.5)
    success_message_utest()
    check_condition_utest('Test if dates retrieved match dates added', all([abs((times_data[ii]-time_files[ii]).total_seconds())<=dt_date_tolerance_seconds for ii in range(len(time_files))]))
        
    print_utest_message('Test build_forcing_data')
    times_data, data_retrieved = db.build_forcing_data(time_files, data_type='rain', dt_max=0.)
    success_message_utest()
    check_condition_utest('Test if dates retrieved match dates added', all([abs((times_data[ii]-time_files[ii]).total_seconds())<=dt_date_tolerance_seconds for ii in range(len(time_files))]))
    check_condition_utest('Test if data retrieved is correct', all([np.any(data_retrieved[:,ii]==np.float32(ii)) for ii in range(len(time_files))]))
        
    #build forcing file
    print_utest_message('Test build_forcing_file_chunked')
    file_out = os.path.join(fake_data_dir, 'test_build_forcing_file_chunked.nc')
    db.build_forcing_file_chunked(time_files, file_out, data_type='rain', max_chunksize=2*len(time_files), dt_max=0.)
    success_message_utest()
    check_condition_utest('Test if file has been written', os.path.exists(file_out))
        
    #build forcing file by chunks
    print_utest_message('Test build_forcing_file_chunked by chunks')
    file_out_chunk = os.path.join(fake_data_dir, 'test_build_forcing_file_chunked_chunk.nc')
    db.build_forcing_file_chunked(time_files, file_out_chunk, data_type='rain', max_chunksize=len(time_files)//3, dt_max=0.)
    success_message_utest()
    check_condition_utest('Test if file has been written', os.path.exists(file_out_chunk))
    
    #check that the two forcing files are identical
    are_identical = True
    ds_in = netCDF4.Dataset(file_out)
    ds_in_chunk = netCDF4.Dataset(file_out_chunk)
    for el in ds_in.variables:
        var1 = ds_in.variables[el][:]
        var2 = ds_in_chunk.variables[el][:]
        if 2.*np.mean(np.abs(var1-var2))/np.mean(np.abs(var1+var2)) > 1.e-6:
            print('Variable %s from forcing file written by chunks %s not identical to the one written in one shot %s'%(el, file_out, file_out_chunk))
            are_identical = False
    ds_in.close()
    ds_in_chunk.close()
    check_condition_utest('Test if file written by chunks and file written in one shot are identical', are_identical)
    
    
    print_utest_message('test if written file can be read')
    ds = netCDF4.Dataset(file_out)
    times_data, data_retrieved = ds.variables['dates_original'][:], ds.variables['rain'][:]
    ds.close()
    success_message_utest()
    
    check_condition_utest('Test if dates retrieved match dates added', all([abs((julianday_to_datetime(times_data[ii])-time_files[ii]).total_seconds())<=dt_date_tolerance_seconds for ii in range(len(time_files))]))
    check_condition_utest('Test if data retrieved is correct', all([np.any(data_retrieved[:,ii]==np.float32(ii)) for ii in range(len(time_files))]))
    db._close_()
        

    #8 Test if build_forcing_file_chunked with inconsistent mesh size fails (as it should)
    print_utest_message('Test opening database in read mode')
    with ForcingOnMesh_DBManager(sql_db_dir, mode='r', verbose=0) as db:
        success_message_utest()
        
        print_utest_message('Test if build_forcing_file_chunked with inconsistent mesh size fails (as it should)')
        file_out = os.path.join(fake_data_dir, 'test_build_forcing_file_chunked2.nc')
        try:
            db.build_forcing_file_chunked(time_files + [time_improper_file], file_out, data_type='rain', dt_max=0.)
            fail_message_utest()
        except:
            success_message_utest()

    shutil.rmtree(main_test_dir)
    



