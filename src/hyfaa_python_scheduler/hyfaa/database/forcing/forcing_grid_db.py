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
import sqlite3
import pandas

from hyfaa.common.test.common_test_functions import *
from hyfaa.database.forcing.forcing_db import Forcing_DBManager
from hyfaa.common.timer.easy_timer import Timer



class ForcingGrid_DBManager(Forcing_DBManager):
    """This class handles the forcing on-mesh database.
    """
    
    def __init__(self, db_sql_path, mode='r', verbose=0):
        super().__init__(db_sql_path, mode=mode, verbose=verbose)
    
    def __enter__(self):
        return super().__enter__()
    
    def __exit__(self, type, value, traceback):
        return super().__exit__(type, value, traceback)


  

        
def test_db():
    """main test for database object
    """
    
    print_utest_message('\n\nRunning main forcing grid database manager diagnostics:\n')
    
    main_test_dir = 'temp_%s'%(date2strtag(datetime.now()))
    
    input_data_dir = '%s/era5-like_data'%main_test_dir
    sql_db_dir = '%s/forcing_grid_db'%main_test_dir
    
    for fol in [input_data_dir, sql_db_dir]:
        if not os.path.exists(fol):
            os.system('mkdir -p %s'%fol)
            
            
    #make input data
    filenames = dict()
    nlon_in, nlat_in = 200, 121
    lon_in, lat_in = np.linspace(0., 360.-360./nlon_in, nlon_in), np.linspace(90., -90., nlat_in)
    for year in range(2000,2001):
        for month in range(1,3):
            for day in range(1,monthrange(year,month)[1]+1):
                for hour in range(0,24,3):
                    filenames[(year,month,day,hour)] = '%s/era5like_%d_%d_%d_%d.nc'%(input_data_dir, year, month, day, hour)
                    write_forcing_grid(lon_in, lat_in, np.random.rand(len(lon_in)*len(lat_in)).reshape((len(lat_in), len(lon_in))), \
                        filenames[(year,month,day,hour)], {'time': date2str(datetime(year,month,day,hour))}, verbose=0)
                
                
    #1 : test opening empty database in read mode
    print_utest_message('Test that opening empty database in read mode fails')
    try:
        with ForcingGrid_DBManager(sql_db_dir, mode='r', verbose=0) as db:
            fail_message_utest()
    except:
        success_message_utest()
    check_condition_utest('Test that opening empty database in read mode does not generate file creation', len(os.listdir(sql_db_dir)) == 0)

    
    #2 : test opening empty database in write mode
    print_utest_message('Test that opening empty database in write mode succeeds')
    try:
        with ForcingGrid_DBManager(sql_db_dir, mode='w', verbose=0) as db:
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
        db = ForcingGrid_DBManager(sql_db_dir, mode='r', verbose=0)
        db._close_()
        success_message_utest()
    except:
        fail_message_utest()
    #4 : check that it fails in write mode
    print_utest_message('Test that opening database in write mode without context manager fails')
    try:
        db = ForcingGrid_DBManager(sql_db_dir, mode='w', verbose=0)
        db._close_()
        fail_message_utest()
    except:
        success_message_utest()
    
    
    
    #5 : manually create lock file and test database opening (should be ok in read mode and fail in write mode)
    shutil.copy('%s/database_manager.sql'%sql_db_dir, '%s/database_manager.sql_lock'%sql_db_dir)
    #read mode
    print_utest_message('Test reading while lock is active')
    try:
        with ForcingGrid_DBManager(sql_db_dir, mode='r', verbose=0) as db:
            success_message_utest()
    except:
        shutil.rmtree(main_test_dir)
        fail_message_utest()
    #write mode
    print_utest_message('Test if write is prevented while lock is active')
    try:
        with ForcingGrid_DBManager(sql_db_dir, mode='w', verbose=0) as db:
            fail_message_utest()
    except:
        success_message_utest()
    #check that lock file still exists after ForcingGrid_DBManager exit in write mode
    check_condition_utest('Test that lock file still exists when ForcingGrid_DBManager exits because another instance of the manager is already opened in write mode', \
        os.path.exists('%s/database_manager.sql_lock'%sql_db_dir))
    os.unlink('%s/database_manager.sql_lock'%sql_db_dir)
    
                
                
    #6 write test create DB
    print_utest_message('Test creating database')
    dates_check = []
    try:
        with ForcingGrid_DBManager(sql_db_dir, mode='w') as db:
            for key in sorted(filenames.keys()):
                date_loc = datetime(key[0],key[1],key[2],key[3])
                dates_check.append(date_loc)
                db.add({'file_path': filenames[key], 'data_type': 'rain', 'date_data': date_loc, 'product_type': 'analysis', 'grid_status': 'complete'})
        success_message_utest()
    except:
        fail_message_utest()
        
        
    #7 read test
    print_utest_message('Test reading database')
    try:
        db = ForcingGrid_DBManager(sql_db_dir)
        dates_of_files_chosen, files_chosen = db.get_best_paths_matching_dates_and_data_type(dates_check, data_type='rain', dt_max=0.01, error_if_none_found=False)
        db._close_()
        success_message_utest()
    except:
        fail_message_utest()
        
    check_condition_utest('Check that dates recovered are an exact match', \
        all([abs((dates_check[ii]-dates_of_files_chosen[ii]).total_seconds())<0.01 for ii in range(len(dates_check))]))
    
    shutil.rmtree(main_test_dir)
        


