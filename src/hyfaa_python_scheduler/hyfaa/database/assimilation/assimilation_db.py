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



class Assimilation_Database(object):
    """This class handles the forcing on-mesh database.
    """
    
    def __init__(self, db_sql_path, mode='r', verbose=1):
        """Basic init function ; most initialization operations are deferred to the initialization function triggered by __enter__
        i.e. entering context manager
        
        :param file_path: sql databse path (sqlite3 file)
        :param mode: string 'r' for read_only or 'w' for write
        :param verbose: int representing verbose level
        """
        
        self._exec_date = datetime.now()
        if not os.path.exists(db_sql_path):
            raise Exception('path %s does not exist'%db_sql_path)
        if not os.path.isdir(db_sql_path):
            raise Exception('path %s is not a directory'%db_sql_path)
        db_sql_path = os.path.abspath(db_sql_path)
        self._path = os.path.join(db_sql_path, 'database.sql')
        if mode not in ['r', 'w']:
            raise Exception('Mode must be either r: read only or w: read and write')
        self._mode = mode
        self._verbose = verbose
        self.__within_context = False
        self.__special_exit = None
        self.__initialized = False
        if self._mode == 'r':
            self.__initialization__()
        
    def __enter__(self):
        """Enter context and launch initialization.
        Entering context is required for this object to function properly with write operations 
        (see paragraph about SQL database security in object header)."""
        self.__within_context = True
        if not self.__initialized:
            self.__initialization__()
        return self
            
            
    def __exit__(self, type, value, traceback):
        """Removes _lock file if it was created and exits context."""
        self._conn.close()
        if self._mode == 'w' and self.__special_exit != 'lock exists':
            if os.path.exists('%s_lock'%self._path):
                os.unlink('%s_lock'%self._path)
            else:
                print('Exiting write mode but no lock file exists => should not happen !')
                return 1
        return 0
        
    
    def _close_(self):
        if self._mode == 'r':
            self._conn.close()
        else:
            raise Exception('close function only allowed in read mode')
        
        
    def _check_within_context_(self):
        """Checks if object was opened with context manager : required for all write operations"""
        if not self.__within_context:
            raise Exception('Object must be within a context')
            
            
    def _check_writable_(self):
        """Checks if mode is write"""
        self._check_within_context_()
        if self._mode != 'w':
            raise Exception('Cannot update database: read only mode')
            
    def check_column_values(self, values):
        """check if column values are ok for database inputs"""
        none_keys = sorted(list(self._necessary_input_columns.intersection(set([elem for elem in self._columns if values[self.column_id[elem]] in [None, 'None']]))))
        if len(none_keys) > 0:
            raise Exception('missing_keys in ForcingOnMesh_DBManager add function parameter file_info:\n%s\n'%('\n'.join([' - %s'%elem for elem in none_keys])))

    def initialize_columns(self):
        """initialize columns"""
        self._columns = dict([('mesh_id', 'int'), ('value', 'float'), ('uncertainty', 'float'), ('date_data', 'str'), \
            ('lon', 'float'), ('lat', 'float'), ('date_added_to_db', 'str'), ('sv_name', 'str'), ('info', 'str')])
        self._necessary_input_columns = set(['mesh_id', 'sv_name', 'date_data', 'value'])
        
        
    def __initialization__(self):
        """The real initialization function is separated from __init__ and called from __enter__ to require entering the context manager.
            - _lock file is created if self._mode is write
            - If database file does not exist, it is created
            - SQL database coherence with initialization parameters is checked:
                - FILEINFO columns must be indentical to self.columns
        """
        
        if self._mode != 'r':
            self._check_within_context_()
        
        if os.path.exists('%s_lock'%self._path) and self._mode == 'w':
            self.__special_exit = 'lock exists'
            raise Exception('It seems that database %s is currently being modified i.e. %s_lock file exists'%(self._path, self._path))
                
            
        #initialize self.columns
        self.initialize_columns()
        self.column_id = {elem: ii for ii,elem in enumerate(self._columns)}
        self.column_is_date = dict()
        for el in self._columns.keys():
            if 'date_' in el:
                self.column_is_date[el] = True
            else:
                self.column_is_date[el] = False

        
        #create db if it does not exist, otherwise simply open it
        if not os.path.exists(self._path):
            #this function creates the table => requires write priviledges
            self._check_writable_()
            if self._verbose > 0:
                print('TableManager database does not exist, creating it:\n  => %s'%self._path)
            self._conn = sqlite3.connect(self._path)
            self._cursor = self._conn.cursor()
            self._cursor.execute('CREATE TABLE FILEINFO (%s)'%(', '.join(['%s %s'%(el, self._columns[el]) for el in self._columns.keys()])))
            self._conn.commit()
            #lock database for writes (we already know mode == 'w')
            shutil.copy(self._path, '%s_lock'%self._path)
        else:
            self._conn = sqlite3.connect(self._path)
            self._cursor = self._conn.cursor()
            if self._mode == 'w':
                #create lock if writing is intended
                shutil.copy(self._path, '%s_lock'%self._path)
            #check columns in db vs columns in tablemanager
            self._cursor.execute('SELECT * FROM FILEINFO')
            columns_db = [el[0] for el in self._cursor.description]
            if set(columns_db) != set(self._columns.keys()):
                raise Exception('columns in db not consistent with columns in tablemanager')

    
    def read_as_pandas_dataframe(self, sql_query, params=None):
        """returns SQL read as a panda dataframe (much more convenient)"""
        return pandas.read_sql_query(sql_query, self._conn, params=params)
        
        
    def convert_column_dates2str(self, info_in, output='list'):
        """converts datetime objects in column entries to datetime strings"""
        if hasattr(info_in, 'keys'):
            items = [(el, el) for el in self._columns.keys()]
        elif hasattr(info_in, '__getitem__'):
            items = [(ii, el) for ii,el in enumerate(self._columns.keys())]
        else:
            raise Exception('Only accepts dict, dict or list')
            
        if output == 'dict':
            return dict([(el1, self.date2str(info_in[el0])) if self.column_is_date[el1] else (el1, info_in[el0]) for el0, el1 in items])
        elif output == 'list':
            return [self.date2str(info_in[el0]) if self.column_is_date[el1] else info_in[el0] for el0, el1 in items]
        else:
            raise Exception('output type %s unkown'%output)
        
        
    def convert_column_str2dates(self, info_in, output='list'):
        """converts datetime strings in column entries to datetime objects"""
        if hasattr(info_in, 'keys'):
            items = [(el, el) for el in self._columns.keys()]
        elif hasattr(info_in, '__getitem__'):
            items = [(ii, el) for el in enumerate(self._columns.keys())]
        else:
            raise Exception('Only accepts dict, dict or list')
        
        if output == 'dict':
            return dict([(el1, self.str2date(info_in[el0])) if self.column_is_date[el1] else (el1, info_in[el0]) for el0, el1 in items])
        elif output == 'list':
            return [self.str2date(info_in[el0]) if self.column_is_date[el1] else info_in[el0] for el0, el1 in items]
        else:
            raise Exception('output type %s unkown'%output)
        
    
    @staticmethod
    def date2str(datetime_object):
        """converts datetime object to datetime string"""
        if datetime_object is None:
            return 'None'
        return datetime_object.strftime('%Y-%m-%dT%H:%M:%S.%f')[0:-3]
    
    
    @staticmethod
    def str2date(string):
        """converts datetime string to datetime object"""
        if string is 'None':
            return None
        return datetime.strptime(string + '000', '%Y-%m-%dT%H:%M:%S.%f')
        

    def add(self, file_infos):
        """Add a file to SQL database"""
        self._check_writable_()
        
        for file_info in file_infos:
            #columns = mesh_id, value, date_data, lon, lat, date_added_to_db, sv_name, info
            #add file to db with status adding
            file_info['date_added_to_db'] = datetime.now()
            list_write = [file_info[el] if el in file_info else None for el in self._columns.keys()]
            #check for proper inputs
            self.check_column_values(list_write)
            
            #add to db
            self._cursor.execute('INSERT INTO FILEINFO VALUES (%s)'%(','.join(['?' for el in self._columns.keys()])), tuple(self.convert_column_dates2str(list_write)))
        self._conn.commit()
        
        
    def remove(self, mesh_sv_date_tuples):
        for mesh_id, sv_name, date_data in mesh_sv_date_tuples:
            self._cursor.execute('DELETE FROM FILEINFO WHERE mesh_id=? AND sv_name=? and date_data=?', (mesh_id, sv_name, self.date2str(date_data)))
        self._conn.commit()
        
    def get_values_between_dates(self, date_start=None, date_end=None, dt_max=0.0, start_strict=False, end_strict=True):
        """Retrieves values between 2 dates (sorted by date).
        dt_max is the maximum date difference between data date and dates_in allowed => to retrieve exact dates, input dt_max=0.
            
        :param date_start: start date (datetime object)
        :param date_end: end date (datetime object)
        :param dt_max: positive float
        """
        
        if start_strict:
            start_diff_operator = '>'
        else:
            start_diff_operator = '>='
        if end_strict:
            end_diff_operator = '<'
        else:
            end_diff_operator = '<='
            
        if dt_max < 0.:
            raise Exception('dt_max must be > 0')
            
        if (date_start is not None) and (date_end is not None):
            return self.read_as_pandas_dataframe("SELECT * FROM FILEINFO WHERE datetime(date_data) %s datetime(?) AND datetime(date_data) %s datetime(?) ORDER BY datetime(date_data)"%(start_diff_operator, end_diff_operator), \
                params=[self.date2str(date_start-timedelta(dt_max)), self.date2str(date_end+timedelta(dt_max))])
        elif (date_start is not None):
            return self.read_as_pandas_dataframe("SELECT * FROM FILEINFO WHERE datetime(date_data) %s datetime(?) ORDER BY datetime(date_data)"%start_diff_operator, \
                params=[self.date2str(date_start-timedelta(dt_max))])
        elif (date_end is not None):
            return self.read_as_pandas_dataframe("SELECT * FROM FILEINFO WHERE datetime(date_data) %s datetime(?) ORDER BY datetime(date_data)"%end_diff_operator, \
                params=[self.date2str(date_end+timedelta(dt_max))])
        else:
            return self.read_as_pandas_dataframe("SELECT * FROM FILEINFO ORDER BY datetime(date_data)")

    




def test_db():
    """main test for database object
    """
    
    print_utest_message('\n\nRunning main assimilation database diagnostics:\n')
    
    main_test_dir = tempfile.mkdtemp(dir=os.path.abspath('.'), prefix='testaassimdb_')
    if os.path.exists(main_test_dir):
        shutil.rmtree(main_test_dir)
    sql_db_dir = '%s/assimilation_database'%main_test_dir
        
    #create directories
    for fol in [sql_db_dir]:
        os.makedirs(fol, exist_ok=True)
            
    
    #1 : test opening empty database in read mode
    print_utest_message('Test that opening empty database in read mode fails')
    try:
        with Assimilation_Database(sql_db_dir, mode='r', verbose=0) as db:
            fail_message_utest()
    except:
        success_message_utest()
    check_condition_utest('Test that opening empty database in read mode does not generate file creation', len(os.listdir(sql_db_dir)) == 0)

    
    #2 : test opening empty database in write mode
    print_utest_message('Test that opening empty database in write mode succeeds')
    try:
        with Assimilation_Database(sql_db_dir, mode='w', verbose=0) as db:
            success_message_utest()
            check_condition_utest('Test that opening empty database in write mode creates necessary files', \
                all([os.path.exists(el) for el in ['%s/database.sql'%sql_db_dir, '%s/database.sql_lock'%sql_db_dir]]))
    except:
        shutil.rmtree(main_test_dir)
        fail_message_utest()
    check_condition_utest('Test that lock file is removed upon database closure', not os.path.exists('%s/database.sql_lock'%sql_db_dir))
    
    
    
    #3 : check that opening database without context manager succeeds in read mode
    print_utest_message('Test that opening database in read mode without context manager succeeds')
    try:
        db = Assimilation_Database(sql_db_dir, mode='r', verbose=0)
        db._close_()
        success_message_utest()
    except:
        fail_message_utest()
    #4 : check that it fails in write mode
    print_utest_message('Test that opening database in write mode without context manager fails')
    try:
        db = Assimilation_Database(sql_db_dir, mode='w', verbose=0)
        db._close_()
        fail_message_utest()
    except:
        success_message_utest()
    
    
    #5 : manually create lock file and test database opening (should be ok in read mode and fail in write mode)
    shutil.copy('%s/database.sql'%sql_db_dir, '%s/database.sql_lock'%sql_db_dir)
    #read mode
    print_utest_message('Test reading while lock is active')
    try:
        with Assimilation_Database(sql_db_dir, mode='r', verbose=0) as db:
            success_message_utest()
    except:
        shutil.rmtree(main_test_dir)
        fail_message_utest()

    #write mode
    print_utest_message('Test if write is prevented while lock is active')
    try:
        with Assimilation_Database(sql_db_dir, mode='w', verbose=0) as db:
            fail_message_utest()
    except:
        success_message_utest()

    #check that lock file still exists after Assimilation_Database exits in write mode
    check_condition_utest('Test that lock file still exists when Assimilation_Database exits because another instance of the manager is already opened in write mode', \
        os.path.exists('%s/database.sql_lock'%sql_db_dir))
    os.unlink('%s/database.sql_lock'%sql_db_dir)
    
        
    #6 test writing
    print_utest_message('Test writing in database')
    meas_info = []
    n_data = 10
    for ii in range(n_data):
        meas_info.append({'mesh_id': np.random.randint(10000), 'value': np.random.rand(), 'uncertainty': np.random.rand()+0.1, \
            'date_data': datetime(1993,1,1)+timedelta(ii+np.random.rand()), 'lon': np.random.rand()*360., 'lat': np.random.rand()*180.-90., \
            'sv_name': 'atlantis_south', 'info': ''})
    try:
        with Assimilation_Database(sql_db_dir, mode='w', verbose=0) as db:
            db.add(meas_info)
        success_message_utest()
    except:
        fail_message_utest()

    meas_info = [meas_info[ii] for ii in np.argsort(np.array([el['date_data'] for el in meas_info]))]

    #7 test reading exact match
    print_utest_message('Test that data is correctly recovered')
    db = Assimilation_Database(sql_db_dir, mode='r', verbose=0)
    try:
        data = db.get_values_between_dates(date_start=min([el['date_data'] for el in meas_info]), date_end=max([el['date_data'] for el in meas_info]), dt_max=0.0)
        for ii in range(len(data['date_data'])):
            date_keys = set(['date_data', 'date_added_to_db'])
            for el in data.keys():
                if el in date_keys:
                    if abs((db.str2date(data[el][ii]) - meas_info[ii][el]).total_seconds()) > dt_date_tolerance_seconds:
                        print(db.str2date(data[el][ii]), meas_info[ii][el])
                        raise Exception('%s mismatch'%el)
                else:
                    if data[el][ii] != meas_info[ii][el]:
                        print(data[el][ii], meas_info[ii][el])
                        raise Exception('%s mismatch'%el)
        success_message_utest()
    except:
        fail_message_utest()
    db._close_()


    #8 test removing data
    n_data_remove = 2
    print_utest_message('Test removing data')
    try:
        with Assimilation_Database(sql_db_dir, mode='w', verbose=0) as db:
            db.remove([(meas_info[ii]['mesh_id'], meas_info[ii]['sv_name'], meas_info[ii]['date_data']) for ii in range(n_data_remove)])
        success_message_utest()
    except:
        fail_message_utest()
        
    #9 test that the right data was removed
    print_utest_message('Test that the right data was removed')
    db = Assimilation_Database(sql_db_dir, mode='r', verbose=0)
    try:
        data = db.get_values_between_dates(date_start=min([el['date_data'] for el in meas_info]), date_end=max([el['date_data'] for el in meas_info]), dt_max=0.0)
        for ii in range(len(data['date_data'])):
            date_keys = set(['date_data', 'date_added_to_db'])
            for el in data.keys():
                if el in date_keys:
                    if abs((db.str2date(data[el][ii]) - meas_info[ii+n_data_remove][el]).total_seconds()) > dt_date_tolerance_seconds:
                        print(db.str2date(data[el][ii]), meas_info[ii+n_data_remove][el])
                        raise Exception('%s mismatch'%el)
                else:
                    if data[el][ii] != meas_info[ii+n_data_remove][el]:
                        print(data[el][ii], meas_info[ii+n_data_remove][el])
                        raise Exception('%s mismatch'%el)
        success_message_utest()
    except:
        fail_message_utest()
    db._close_()
        
        
    shutil.rmtree(main_test_dir)
    
    
        
        
        
if __name__ == '__main__':
    
    
    test_db()
    

