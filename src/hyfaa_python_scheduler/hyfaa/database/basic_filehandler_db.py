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



class BasicFileHandler_DBManager(object):
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
        self._path = '%s/database_manager.sql'%db_sql_path
        self._data_path = '%s/data_store'%db_sql_path
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
            if not self.__within_context:
                print('Database has not been initialized => this must be done in write mode...')
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
            
            if self._necessary_input_columns>set(columns_db):
                raise Exception('columns in db not consistent with columns in tablemanager')
                
        #create data folder
        if not os.path.exists(self._data_path):
            os.system('mkdir -p %s'%self._data_path)
            
        #remove all entries/files with status adding => means it failed and that there may be corrupted files taking disk space
        info_db = self.read_as_pandas_dataframe("SELECT file_path FROM FILEINFO WHERE file_status=?", params=['adding'])
        if len(info_db['file_path']) > 0:
            self.remove_files(info_db['file_path'])
            
    def get_data_path(self):
        return self._data_path
    

    
    def read_as_pandas_dataframe(self, sql_query, params=None):
        """returns SQL read as a panda dataframe (much more convenient)"""
        return pandas.read_sql_query(sql_query, self._conn, params=params)
        
        
        
    def get_minmax_dates(self, output_mode='datetime_object'):
        if output_mode not in ['datetime_object', 'datetime_string']:
            raise Exception('output_mode %s unknown'%output_mode)
        minmax_dates = [self.read_as_pandas_dataframe("SELECT MIN(date_data) FROM FILEINFO WHERE file_status=?", params=['added'])['MIN(date_data)'][0], \
            self.read_as_pandas_dataframe("SELECT MAX(date_data) FROM FILEINFO WHERE file_status=?", params=['added'])['MAX(date_data)'][0]]
        if output_mode == 'datetime_string':
            return minmax_dates
        for ii in range(len(minmax_dates)):
            if minmax_dates[ii] is not None:
                minmax_dates[ii] = self.str2date(minmax_dates[ii])
        return minmax_dates
        
    
    def remove_files(self, filenames):
        """Remove files from database"""
        for filename in filenames:
            full_path = os.path.join(self._data_path, filename)
            os.unlink(full_path)
            self._cursor.execute('DELETE FROM FILEINFO WHERE file_path=?', (filename, ))
        self._conn.commit()

        
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
        

    def add(self, file_info, store_file_option='move'):
        """Add a file to SQL database"""
        self._check_writable_()
        
        #add file to db with status adding
        input_file = file_info['file_path']
        file_info['file_status'] = 'adding'
        file_info['date_added_to_db'] = datetime.now()
        file_info['file_path'] = self.filename_from_fileinfo(file_info)
        list_write = [file_info[el] if el in file_info else None for el in self._columns.keys()]
        #check for proper inputs
        self.check_column_values(list_write)
        
        #add to db under status adding
        self._cursor.execute('INSERT INTO FILEINFO VALUES (%s)'%(','.join(['?' for el in self._columns.keys()])), tuple(self.convert_column_dates2str(list_write)))
        self._conn.commit()
        
        #move (or copy) file
        output_file = '%s/%s'%(self._data_path, file_info['file_path'])
        if store_file_option == 'move':
            shutil.move(input_file, output_file)
        elif store_file_option == 'copy':
            shutil.copy(input_file, output_file)
        else:
            raise Exception('store_file_option must be in [copy,move]')
            
        #set status to added
        self._cursor.execute("UPDATE FILEINFO SET file_status='added' WHERE file_path=?", (file_info['file_path'], ))
        self._conn.commit()
        return output_file
        
        
    def get_paths_matching_dates(self, dates_in, dt_max=0.0):
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
            
        info_db = self.read_as_pandas_dataframe("SELECT * FROM FILEINFO WHERE file_status=? AND " + \
                "datetime(date_data) BETWEEN datetime(?) AND datetime(?)", params=['added', self.date2str(dates_in[0]-timedelta(dt_max)), self.date2str(dates_in[-1]+timedelta(dt_max))])
        dates_data = [self.str2date(info_db['date_data'][ii]) for ii in range(len(info_db['date_data']))]
        
        dates_of_files_chosen = []
        files_chosen = []
        for date_loc in dates_in:
            date_min_ok = date_loc-timedelta(dt_max)
            date_max_ok = date_loc+timedelta(dt_max)
            index_match = [ii for ii in range(len(dates_data)) if dates_data[ii]>=date_min_ok and dates_data[ii]<=date_max_ok]
            #return date and files that match in order of closest date and then of added last (i.e. most recent) to added first
            sort_index = list_argsort([(abs((dates_data[ii]-date_loc).total_seconds()), -datetime_to_julianday(self.str2date(info_db['date_added_to_db'][ii]))) for ii in index_match])
            dates_of_files_chosen.append([dates_data[index_match[ii]] for ii in sort_index])
            files_chosen.append(['%s/%s'%(self._data_path, info_db['file_path'][index_match[ii]]) for ii in sort_index])
                
        return dates_of_files_chosen, files_chosen
        
        
    def filename_from_fileinfo(self, file_info):
        """filename used for storage into database : the only important thing is that it be unique so we use date of data, add_date and a random number"""
        return 'DATA_%s_%s_%s.nc'%(file_info['date_data'].strftime('%Y%m%dT%H%M%S%f'), file_info['date_added_to_db'].strftime('%Y%m%dT%H%M%S%f'), '%04d'%np.random.randint(0,9999))
        
    def check_column_values(self, values):
        """empty function meant to be overloaded"""
        raise Exception('this function should be oveloaded')

    def initialize_columns(self):
        """empty function meant to be overloaded"""
        raise Exception('this function should be oveloaded')

        
