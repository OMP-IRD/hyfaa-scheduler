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

import os, sys, shutil, re
if sys.version_info.major < 3:
    raise Exception('sorry, only python 3 supported')
from datetime import datetime, timedelta
import numpy as np
import netCDF4
import sqlite3
import pandas

from hyfaa.common.common_functions import *
from hyfaa.common.test.common_test_functions import *
from hyfaa.database.basic_filehandler_db import BasicFileHandler_DBManager

import multiprocessing
from queue import Empty as queue_empty





class HydroStates_DBManager(BasicFileHandler_DBManager):
    """This class handles the hydrological state database. 
    Instead of a single file path per entry, the file path is here instead considered as a file prefix 
    for multiple files to allow for storage of ensembles.
    """

    def __init__(self, db_sql_path, mode='r', verbose=1, async_writes=False):
        super().__init__(db_sql_path, mode=mode, verbose=verbose)
        self.async_writes = async_writes
        if self.async_writes:
            self.workQueue = multiprocessing.Queue(0)
            self.returnQueue = multiprocessing.Queue(0)
            self.addingList = []
            self.worker = JoinEnsembleWorker(self.workQueue, self.returnQueue)
            self.worker.start()

    def __enter__(self):
        return super().__enter__()
    
    def __exit__(self, type, value, traceback):
        self.update_joined_ensemble_files_from_worker()
        if self.async_writes:
            self.worker.join()
        return super().__exit__(type, value, traceback)

    def check_column_values(self, values):
        """check if column values are ok for database inputs"""
        
        none_keys = sorted(list(self._necessary_input_columns.intersection(set([elem for elem in self._columns if values[self.column_id[elem]] in [None, 'None']]))))
        if len(none_keys) > 0:
            raise Exception('missing_keys in HydroStates_DBManager add function parameter file_info:\n%s\n'%('\n'.join([' - %s'%elem for elem in none_keys])))
        if values[self.column_id['file_status']] not in ['adding', 'added', 'removing', 'removed']:
            raise Exception('file_status must be in [adding,added,removing,removed]')

    def initialize_columns(self):
        """initialize columns"""
        self._columns = dict([('file_path', 'str'), ('file_status', 'str'), ('date_data', 'str'), ('date_added_to_db', 'str'), \
            ('ensemble_size', 'int'), ('forcing_confidence_coefficient', 'float'), ('type', 'str')])
        self._necessary_input_columns = set(['file_path', 'file_status', 'date_data', 'date_added_to_db', 'ensemble_size', 'forcing_confidence_coefficient', 'type'])
        
        
    def add(self, file_info, ordered_ensemble_files):
        """Add files to SQL database => overloads function from parent class"""
        self._check_writable_()
        
        #add file to db with status adding
        if 'type' not in file_info.keys():
            file_info['type'] = 'control'
        file_info['file_status'] = 'adding'
        file_info['date_added_to_db'] = datetime.now()
        file_info['file_path'] = self.filename_from_fileinfo(file_info)
        file_info['ensemble_size'] = len(ordered_ensemble_files)
        list_write = [file_info[el] if el in file_info else None for el in self._columns.keys()]
        #check for proper inputs
        self.check_column_values(list_write)
        
        
        #add to db under status adding
        self._cursor.execute('INSERT INTO FILEINFO VALUES (%s)'%(','.join(['?' for el in self._columns.keys()])), tuple(self.convert_column_dates2str(list_write)))
        self._conn.commit()
        
        #fuse ensemble file into a single file
        if not os.path.exists(os.path.join(self._data_path,file_info['type'])):
            os.makedirs(os.path.join(self._data_path,file_info['type']), exist_ok=True)
        self.join_ensemble_members(ordered_ensemble_files, os.path.join(self._data_path, file_info['type'], file_info['file_path']), compression_level=4)


    def update_joined_ensemble_files_from_worker(self):
        if not self.async_writes:
            return
        while True:
            if len(self.addingList) == 0:
                return
            try:
                output_ensemble_file = self.returnQueue.get(block=True,timeout=0.01)
                self._cursor.execute("UPDATE FILEINFO SET file_status='added' WHERE file_path=?", (os.path.basename(output_ensemble_file), ))
                self._conn.commit()
                self.addingList.remove(output_ensemble_file)
            except queue_empty:
                pass
            except:
                raise


    def join_ensemble_members(self, ordered_ensemble_files, output_ensemble_file, compression_level=4):
        if self.async_writes:
            self.workQueue.put((ordered_ensemble_files, output_ensemble_file, compression_level))
            self.addingList.append(output_ensemble_file)
        else:
            join_ensemble_members(ordered_ensemble_files, output_ensemble_file, compression_level=compression_level)
            self._cursor.execute("UPDATE FILEINFO SET file_status='added' WHERE file_path=?", (os.path.basename(output_ensemble_file), ))
            self._conn.commit()
            
            
    def get_full_path(self, dico_in):
        return os.path.join(self._data_path, dico_in['type'], dico_in['file_path'])
            
    
    def get_paths_matching_dates(self, dates_in, type_request=None, dt_max=0.0):
        """Retrieves filenames matching dates_in.
        dt_max is the maximum date difference between data date and dates_in allowed => to retrieve exact dates, input dt_max=0.
            
        :param dates_in: list of dates (datetimes)
        :param dt_max: positive float
        """
        
        if not is_sorted(dates_in):
            raise Exception('input dates must be sorted')

        if dt_max < 0.:
            raise Exception('dt_max must be > 0')
        elif dt_max < dt_date_tolerance_days:
            dt_max = dt_date_tolerance_days
            
        if type_request is None:
            info_db = self.read_as_pandas_dataframe("SELECT * FROM FILEINFO WHERE file_status=? AND datetime(date_data) BETWEEN datetime(?) AND datetime(?)", \
                    params=['added', self.date2str(dates_in[0]-timedelta(dt_max)), self.date2str(dates_in[-1]+timedelta(dt_max))])
        else:
            info_db = self.read_as_pandas_dataframe("SELECT * FROM FILEINFO WHERE file_status=? AND type=? AND datetime(date_data) BETWEEN datetime(?) AND datetime(?)", \
                    params=['added', type_request, self.date2str(dates_in[0]-timedelta(dt_max)), self.date2str(dates_in[-1]+timedelta(dt_max))])
        dates_data = [self.str2date(info_db['date_data'][ii]) for ii in range(len(info_db['date_data']))]
        
        dates_of_files_chosen = []
        files_chosen = []
        
        for date_loc in dates_in:
            date_min_ok = date_loc-timedelta(dt_max)
            date_max_ok = date_loc+timedelta(dt_max)
            index_match = [ii for ii in range(len(dates_data)) if dates_data[ii]>=date_min_ok and dates_data[ii]<=date_max_ok]
            #return date and files that match in order of closest date and then of added last (i.e. most recent) to added first
            #it also ensures that the analysis result is returned first
            sort_index = list_argsort([(abs((dates_data[ii]-date_loc).total_seconds()), -datetime_to_julianday(self.str2date(info_db['date_added_to_db'][ii]))) for ii in index_match])
            dates_of_files_chosen.append([dates_data[index_match[ii]] for ii in sort_index])
            files_chosen.append([self.get_full_path(info_db.loc[index_match[ii]]) for ii in sort_index])
    
        return dates_of_files_chosen, files_chosen
        
        
    def remove_files(self, filenames):
        """Remove files from database"""
        raise NotImplementedError('not implemented yet')
        
        
    def remove_after_date(self, datetime_last):
        info_db = self.read_as_pandas_dataframe("SELECT date_data,type,file_path FROM FILEINFO WHERE datetime(date_data) > datetime(?)", params=[self.date2str(datetime_last)])
        for i0 in range(len(info_db['date_data'])):
            print('Removed date %s'%info_db['date_data'][i0])
            os.unlink(self.get_full_path(info_db.loc[i0]))
            self._cursor.execute('DELETE FROM FILEINFO WHERE type=? AND file_path=?', (info_db['type'][i0], info_db['file_path'][i0], ))
        self._conn.commit()



class JoinEnsembleWorker(multiprocessing.Process):
    
    def __init__(self, workQueue, returnQueue):
        super(JoinEnsembleWorker, self).__init__()
        self.workQueue = workQueue
        self.returnQueue = returnQueue
        
    def run(self):
        while True:
            try:
                tuple_seq = self.workQueue.get(block=True,timeout=0.01)
                if tuple_seq is None:
                    break
                ordered_ensemble_files, output_ensemble_file, compression_level = tuple_seq
                join_ensemble_members(ordered_ensemble_files, output_ensemble_file, compression_level=compression_level)
                self.returnQueue.put(output_ensemble_file)
            except queue_empty:
                pass
            except:
                raise
                return
        
        
        
        

def join_ensemble_members(ordered_ensemble_files, output_ensemble_file, compression_level=4):
    """takes hydrological state files (1 for each ensemble member) and joins them into a single ensemble hydrological state file"""
    
    common_dims = ['n_meshes', 'n_meshes_add1', 'n_hrus']
    separate_dims = ['n_flat_muskingum']
    par_dims = ['n_cells', 'n_soil_types', 'n_basins', 'n_months_year', 'n_regions', 'n_flood_plain_points_max',\
                'n_faces', 'n_deltas', 'n_specific_outlets']

    dsnew = netCDF4.Dataset(output_ensemble_file, mode='w')
    dsnew.createDimension('n_ensemble', len(ordered_ensemble_files))
    
    #copy dimensions
    #iterate on input files
    for i0, input_file in enumerate(ordered_ensemble_files):
        with netCDF4.Dataset(input_file) as ds:
            if i0 == 0:
                for el in common_dims+par_dims:
                    if el in ds.dimensions:
                        dsnew.createDimension(el, ds.dimensions[el].size)
                #copy general attributes
                for el in set(ds.ncattrs()) - set(['ensemble_member_id']):
                    dsnew.setncattr(el, ds.getncattr(el))
            for el in separate_dims:
                dsnew.createDimension('%s_%d'%(el, i0), ds.dimensions[el].size)

        
    #iterate on input files
    for i0, input_file in enumerate(ordered_ensemble_files):
        with netCDF4.Dataset(input_file) as ds:
            #copy variables
            for elem in ds.variables.keys():
                #create variable
                dsvar = dsnew.createVariable(elem + '_%d'%i0, ds.variables[elem].dtype, [el if el in common_dims+par_dims else '%s_%d'%(el, i0) for el in ds.variables[elem].dimensions], \
                    zlib=True, complevel=compression_level, shuffle=True)
                #copy variable attributes
                ncattrs_already_added = dsvar.ncattrs()
                for el in ds.variables[elem].ncattrs():
                    if el not in ncattrs_already_added:
                        dsvar.setncattr(el, ds.variables[elem].getncattr(el))
                #upload data to new netCDF file
                dsvar[:] = ds.variables[elem][:]

    dsnew.close()

        
        
        

def extract_ensemble_member(input_ensemble_file, output_file, ensemble_member, compression_level=4):
    """takes 1 ensemble hydrological state file and extracts an hydrological state file corresponding to a single ensemble member"""
    
    with netCDF4.Dataset(input_ensemble_file) as ds, netCDF4.Dataset(output_file, mode='w') as dsnew:
        
        n_ensemble = ds.dimensions['n_ensemble'].size
        assert ensemble_member < n_ensemble
        par_dims=['n_cells','n_soil_types','n_basins','n_months_year','n_regions','n_flood_plain_points_max','n_faces''n_deltas','n_specific_outlets']
        #copy dimensions
        dim_converter = dict()
        for dim in ['n_meshes', 'n_meshes_add1', 'n_hrus', 'n_flat_muskingum']:
            if dim == 'n_flat_muskingum':
                dim_converter[dim] = dim + '_%d'%ensemble_member
            else:
                dim_converter[dim] = dim
            dsnew.createDimension(dim, ds.dimensions[dim_converter[dim]].size)
        for dim in par_dims:
            if dim in ds.dimensions:
                dim_converter[dim] = dim
                dsnew.createDimension(dim, ds.dimensions[dim].size)

        dim_converter_rev = {value: key for key,value in dim_converter.items()}
        
        #copy variables
        for elem in ds.variables.keys():
            try:
                elem_out = '_'.join(elem.split('_')[0:-1])
                ens_member_loc = int(elem.split('_')[-1])
                assert ens_member_loc < n_ensemble
            except:
                raise Exception('variable name %s does not end by an ensemble member ID')
            if ens_member_loc != ensemble_member:
                continue
            #create variable
            dsvar = dsnew.createVariable(elem_out, ds.variables[elem].dtype, [dim_converter_rev[dim] for dim in ds.variables[elem].dimensions], \
                zlib=True, complevel=compression_level, shuffle=True)
            #copy variable attributes
            ncattrs_already_added = dsvar.ncattrs()
            for el in ds.variables[elem].ncattrs():
                if el not in ncattrs_already_added:
                    dsvar.setncattr(el, ds.variables[elem].getncattr(el))
            #upload data to new netCDF file
            dsvar[:] = ds.variables[elem][:]
            
        #copy general attributes
        for el in ds.ncattrs():
            dsnew.setncattr(el, ds.getncattr(el))
        dsnew.setncattr('ensemble_member_id', ensemble_member)
    


def generate_ensemble_files(hydrostate_file, n_ensemble, empty_folder_store_hydrostates):
    
    if hydrostate_file is None:
        #then we simply return None for all ensemble members => MGB-IPH code will make the initialization
        return [None]*n_ensemble, False

    with netCDF4.Dataset(hydrostate_file) as ds:
        if 'n_ensemble' in ds.dimensions.keys():
            assert ds.dimensions['n_ensemble'].size == n_ensemble, 'ensemble size %d from file %s does not match ensemble size requested %d'%(ds.dimensions['n_ensemble'].size, hydrostate_file, n_ensemble)
        else:
            #then we guess it is not an ensemble file we simply use the same input file for all ensemble members
            for ii in range(n_ensemble):
                shutil.copy(hydrostate_file, os.path.join(empty_folder_store_hydrostates,'hydrostate_in_%d.nc'%ii))
            return [os.path.join(empty_folder_store_hydrostates,'hydrostate_in_%d.nc'%ii) for ii in range(n_ensemble)], True

    #then we extract each ensemble member to a separate file
    hydrostates_ensemble_files = []
    for ii in range(n_ensemble):
        output_file = os.path.join(empty_folder_store_hydrostates, 'hydrostate_in_%d.nc'%ii)
        extract_ensemble_member(hydrostate_file, output_file, ii)
        hydrostates_ensemble_files.append(output_file)
    return hydrostates_ensemble_files, True



def test_db():
    """main test for database object
    """
    
    print_utest_message('\n\nRunning main hydrostate database manager diagnostics:\n')
    
    main_test_dir = 'temp_%s'%(date2strtag(datetime.now()))
    if os.path.exists(main_test_dir):
        os.system('rm -R -f %s'%main_test_dir)
    sql_db_dir = '%s/forcing_database'%main_test_dir
    fake_data_dir = '%s/fake_data'%main_test_dir
        
    #create directories
    for fol in [sql_db_dir, fake_data_dir]:
        if not os.path.exists(fol):
            os.system('mkdir -p %s'%fol)
            
    #create test data
    ensemble_groups = []
    n_ensemble = 10
    n_groups = 7
    dt_data = 1.
    time_files = []
    for ii in range(n_groups):
        time_files.append(datetime(2011,1,27)+timedelta(dt_data*ii))
        loc_file = '%s/hydrostate_test_%d_mean.nc'%(fake_data_dir, ii)
        make_hydrostate(loc_file)
        loc_dico = {'input_ensemble_files_ordered': [], 'nonensemble_input_files': {'mean': loc_file}}
        for i0 in range(n_ensemble):
            loc_file = '%s/hydrostate_test_%d_%d.nc'%(fake_data_dir, ii, i0)
            make_hydrostate(loc_file)
            loc_dico['input_ensemble_files_ordered'].append(loc_file)
        ensemble_groups.append(loc_dico)
            
    
    #1 : test opening empty database in read mode
    print_utest_message('Test that opening empty database in read mode fails')
    try:
        with HydroStates_DBManager(sql_db_dir, mode='r', verbose=0) as db:
            fail_message_utest()
    except:
        success_message_utest()
    check_condition_utest('Test that opening empty database in read mode does not generate file creation', len(os.listdir(sql_db_dir)) == 0)

    
    #2 : test opening empty database in write mode
    print_utest_message('Test that opening empty database in write mode succeeds')
    try:
        with HydroStates_DBManager(sql_db_dir, mode='w', verbose=0) as db:
            success_message_utest()
            check_condition_utest('Test that opening empty database in write mode creates necessary files', \
                all([os.path.exists(el) for el in ['%s/database_manager.sql'%sql_db_dir, '%s/database_manager.sql_lock'%sql_db_dir, '%s/data_store'%sql_db_dir]]))
    except:
        os.system('rm -R -f %s'%main_test_dir)
        fail_message_utest()
    check_condition_utest('Test that lock file is removed upon database closure', not os.path.exists('%s/database_manager.sql_lock'%sql_db_dir))
    
    
    
    #3 : check that opening database without context manager succeeds in read mode
    print_utest_message('Test that opening database in read mode without context manager succeeds')
    try:
        db = HydroStates_DBManager(sql_db_dir, mode='r', verbose=0)
        db._close_()
        success_message_utest()
    except:
        fail_message_utest()
    #4 : check that it fails in write mode
    print_utest_message('Test that opening database in write mode without context manager fails')
    try:
        db = HydroStates_DBManager(sql_db_dir, mode='w', verbose=0)
        db._close_()
        fail_message_utest()
    except:
        success_message_utest()
    
    
    
    #5 : manually create lock file and test database opening (should be ok in read mode and fail in write mode)
    shutil.copy('%s/database_manager.sql'%sql_db_dir, '%s/database_manager.sql_lock'%sql_db_dir)
    #read mode
    print_utest_message('Test reading while lock is active')
    try:
        with HydroStates_DBManager(sql_db_dir, mode='r', verbose=0) as db:
            success_message_utest()
    except:
        os.system('rm -R -f %s'%main_test_dir)
        fail_message_utest()
    #write mode
    print_utest_message('Test if write is prevented while lock is active')
    try:
        with HydroStates_DBManager(sql_db_dir, mode='w', verbose=0) as db:
            fail_message_utest()
    except:
        success_message_utest()
    #check that lock file still exists after HydroStates_DBManager exit in write mode
    check_condition_utest('Test that lock file still exists when HydroStates_DBManager exits because another instance of the manager is already opened in write mode', \
        os.path.exists('%s/database_manager.sql_lock'%sql_db_dir))
    os.unlink('%s/database_manager.sql_lock'%sql_db_dir)
    
        
    #6 test writing
    print_utest_message('Test writing in database')
    try:
        with HydroStates_DBManager(sql_db_dir, mode='w', verbose=0) as db:
            for ii in range(n_groups):
                db.add({'date_data': time_files[ii], 'forcing_confidence_coefficient': 1.0, 'number_obs_used': 10, 'type': 'analysis'}, \
                    ensemble_groups[ii]['input_ensemble_files_ordered'])
        success_message_utest()
    except:
        os.system('rm -R -f %s'%main_test_dir)
        fail_message_utest()


    #7 test reading
    db = HydroStates_DBManager(sql_db_dir, mode='r', verbose=0)
    
    print_utest_message('Test get_paths_matching_dates for exact dates')
    times_data, file_paths = db.get_paths_matching_dates(time_files, dt_max=0.)
    success_message_utest()
    check_condition_utest('Test that a single date is retrieved for all dates added', all([len(times_data[ii]) == 1 for ii in range(len(times_data))]))
    check_condition_utest('Test if dates retrieved match dates added', all([abs((times_data[ii][0]-time_files[ii]).total_seconds())<=dt_date_tolerance_seconds for ii in range(len(time_files))]))
        
    print_utest_message('Test get_paths_matching_dates for near dates')
    times_data, file_paths = db.get_paths_matching_dates([el+timedelta((np.random.rand(1)[0]-0.5)*0.99) for el in time_files], dt_max=0.5)
    success_message_utest()
    check_condition_utest('Test that a single date is retrieved for all dates added', all([len(times_data[ii]) == 1 for ii in range(len(times_data))]))
    check_condition_utest('Test if dates retrieved match dates added', all([abs((times_data[ii][0]-time_files[ii]).total_seconds())<=dt_date_tolerance_seconds for ii in range(len(time_files))]))
        
    print_utest_message('Test get_paths_matching_dates for near dates with large dt_max')
    times_data, file_paths = db.get_paths_matching_dates([el+timedelta((np.random.rand(1)[0]-0.5)*0.99) for el in time_files], dt_max=2.5)
    success_message_utest()
    check_condition_utest('Test if dates retrieved match dates added', all([abs((times_data[ii][0]-time_files[ii]).total_seconds())<=dt_date_tolerance_seconds for ii in range(len(time_files))]))
    
    db._close_()
        


    os.system('rm -R -f %s'%main_test_dir)
    
    
    
if __name__ == '__main__':
    
    test_db()



