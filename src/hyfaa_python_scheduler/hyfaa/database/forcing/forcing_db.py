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
from hyfaa.database.basic_filehandler_db import BasicFileHandler_DBManager



class Forcing_DBManager(BasicFileHandler_DBManager):
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

    def check_column_values(self, values):
        """check if column values are ok for database inputs"""
        
        none_keys = sorted(list(self._necessary_input_columns.intersection(set([elem for elem in self._columns if values[self.column_id[elem]] in [None, 'None']]))))
        if len(none_keys) > 0:
            raise Exception('missing_keys in ForcingOnMesh_DBManager add function parameter file_info:\n%s\n'%('\n'.join([' - %s'%elem for elem in none_keys])))
        if values[self.column_id['product_type']] not in ['analysis', 'forecast']:
            raise Exception('product_type must be in [analysis,forecast]')
        if values[self.column_id['file_status']] not in ['adding', 'added', 'removing', 'removed']:
            raise Exception('file_status must be in [adding,added,removing,removed]')
        if values[self.column_id['grid_status']] not in ['complete', 'incomplete', 'replaced']:
            raise Exception('grid_status must be in [complete,incomplete,replaced]')

    def initialize_columns(self):
        """initialize columns"""
        self._columns = dict([('file_path', 'str'), ('data_type', 'str'), ('date_data', 'str'), ('date_added_to_db', 'str'), \
            ('date_created', 'str'), ('product_type', 'str'), ('file_status', 'str'), ('grid_status', 'str')])
        self._necessary_input_columns = set(['file_path', 'data_type', 'date_data', 'product_type', 'file_status', 'grid_status'])
        
        
    def get_priority_coefficient(self, product_type, grid_status, date_difference):
        coeff = self._priority_coefficients['product_type'][product_type] * self._priority_coefficients['grid_status'][grid_status]
        coeff *= np.exp(-2.*date_difference)
        return coeff
        
    @staticmethod
    def check_and_correct_delta(delta_in):
        if delta_in is None:
            delta_in = timedelta(0.001)
        if not isinstance(delta_in, timedelta):
            delta_in = timedelta(delta_in)
        if delta_in < timedelta(dt_date_tolerance_days):
            delta_in = timedelta(dt_date_tolerance_days)
        if delta_in < timedelta(0.):
            raise Exception('delta_in must be > 0')
        return delta_in
        
    
    def get_dates(self, date_min=None, date_max=None, data_type='rain', search_delta_before=None, search_delta_after=None):
        
        search_delta_before = self.check_and_correct_delta(search_delta_before)
        search_delta_after = self.check_and_correct_delta(search_delta_after)
            
        if (date_min is not None) and (date_max is not None):
            sql_request = "SELECT * FROM FILEINFO WHERE file_status=? AND data_type=? AND " + \
                    "datetime(date_data) BETWEEN datetime(?) AND datetime(?) ORDER BY datetime(date_data)"
            params = ['added', data_type, self.date2str(date_min-search_delta_before), self.date2str(date_max+search_delta_after)]
        elif (date_min is not None):
            sql_request = "SELECT * FROM FILEINFO WHERE file_status=? AND data_type=? AND " + \
                    "datetime(date_data) > datetime(?) ORDER BY datetime(date_data)"
            params = ['added', data_type, self.date2str(date_min-search_delta_before)]
        elif (date_max is not None):
            sql_request = "SELECT * FROM FILEINFO WHERE file_status=? AND data_type=? AND " + \
                    "datetime(date_data) < datetime(?) ORDER BY datetime(date_data)"
            params = ['added', data_type, self.date2str(date_max+search_delta_after)]
        else:
            sql_request = "SELECT * FROM FILEINFO WHERE file_status=? AND data_type=? ORDER BY datetime(date_data)"
            params = ['added', data_type]
        info_db = self.read_as_pandas_dataframe(sql_request, params=params)
        return [self.str2date(info_db['date_data'][ii]) for ii in range(len(info_db['date_data']))], info_db
        
    
    def get_best_paths_matching_dates_and_data_type(self, dates_in, data_type='rain', search_delta_before=None, search_delta_after=None, error_if_none_found=False):
        """Retrieves filenames matching dates_in.
        The "best" file will be chosen according to :
            - product_type (analysis is supposed to be better than a forecast)
            - grid_status (complete > incomplete >> replaced (means that it contains climatology instead of real data) )
            - data date difference from dates_in => after 0.1 days, the priority coefficient decreases as 0.1*np.sqrt(1./date_diff) to reach 0.1 for a 1 day difference.
        dt_max is the maximum date difference between data date and dates_in allowed => to retrieve exact dates, input dt_max=0.
            
        :param data_type: rain, ...
        :param dates_in: list of dates (datetimes)
        :param dt_max: positive float
        """
            
        if not is_sorted(dates_in):
            raise Exception('input dates must be sorted')

        search_delta_before = self.check_and_correct_delta(search_delta_before)
        search_delta_after = self.check_and_correct_delta(search_delta_after)
            
        dates_data, info_db = self.get_dates(date_min=dates_in[0], date_max=dates_in[-1], search_delta_before=search_delta_before, search_delta_after=search_delta_after)
        
        dates_of_files_chosen = []
        files_chosen = []
        for date_loc in dates_in:
            date_min_ok = date_loc-search_delta_before
            date_max_ok = date_loc+search_delta_after
            index_match = [ii for ii in range(len(dates_data)) if dates_data[ii]>=date_min_ok and dates_data[ii]<=date_max_ok]
            if len(index_match) == 0:
                if error_if_none_found:
                    raise Exception('No data found in database for %s at date %s'%(data_type, self.date2str(date_loc)))
                else:
                    dates_of_files_chosen.append(None)
                    files_chosen.append(None)
                    continue
            ii_best = None
            coeff = 0.
            last_coeff = 0.
            for ii in index_match:
                date_diff = abs((self.str2date(info_db['date_data'][ii])-date_loc).total_seconds()/(24.*3600.))
                coeff = self.get_priority_coefficient(info_db['product_type'][ii], info_db['grid_status'][ii], date_diff)
                is_best_choice = False
                if coeff > last_coeff:
                    is_best_choice = True
                elif coeff == last_coeff and info_db['date_created'][ii] != 'None':
                    if info_db['date_created'][ii_best] == 'None':
                        is_best_choice = True
                    else:
                        if self.str2date(info_db['date_created'][ii]) > self.str2date(info_db['date_created'][ii_best]):
                            is_best_choice = True
                if is_best_choice:
                    ii_best = ii
                    last_coeff = coeff
            dates_of_files_chosen.append(dates_data[ii_best])
            files_chosen.append('%s/%s'%(self._data_path, info_db['file_path'][ii_best]))
                
        return dates_of_files_chosen, files_chosen
            




