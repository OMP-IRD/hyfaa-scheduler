#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys, shutil
from datetime import datetime, timedelta
import sqlite3
import pandas

    
def convert_db_using_colon_datetags_filenames_to_new_standard(db_path):
    sql_file_path = os.path.join(db_path, 'database_manager.sql')
    data_store_path = os.path.join(db_path, 'data_store')
    assert os.path.exists(sql_file_path)
    assert os.path.exists(data_store_path)
    assert os.path.isdir(data_store_path)
    
    conn = sqlite3.connect(sql_file_path)
    cursor = conn.cursor()
    info_db = pandas.read_sql_query("SELECT file_path FROM FILEINFO WHERE file_status=?", conn, params=['added'])
    for relative_path in info_db['file_path']:
        assert os.path.exists(os.path.join(data_store_path, relative_path))
        filename_split = os.path.basename(relative_path).split('_')
        assert len(filename_split) == 4
        for ii in [1,2]:
            filename_split[ii] = datetime.strptime(filename_split[ii] + '000', '%Y-%m-%dT%H:%M:%S.%f').strftime('%Y%m%dT%H%M%S%f')
        new_relative_path = relative_path.replace(os.path.basename(relative_path), '_'.join(filename_split))
        shutil.move(os.path.join(data_store_path, relative_path), os.path.join(data_store_path, new_relative_path))
        cursor.execute("UPDATE FILEINFO SET file_path='%s' WHERE file_path=?"%new_relative_path, (relative_path, ))
        conn.commit()
        print('%s -> %s OK'%(relative_path, new_relative_path))
        
    print('success converting db !')

if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to convert all database that used colons in date tags in filenames to new standard so it can run on windows file systems.')
    parser.add_argument("db_path", type=str, help="database path")
    args = parser.parse_args()
    
    
    
    convert_db_using_colon_datetags_filenames_to_new_standard(args.db_path)
    


    
    
    
    
