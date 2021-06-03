#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys, shutil
from datetime import datetime, timedelta
import numpy as np
import netCDF4
from forcing_grid_db import ForcingGrid_DBManager

    

def rain_netcdf_to_forcing_database(rain_file, database_directory):
    """Convertit un fichier pbi en une base de donnee de pluie utilisable par le sequenceur
    
    :param rain_file: path to the file containing rain(n_meshes, dates)
    :param database_directory: directory for database (must be empty beforehand)
    """
    
    if not os.path.exists(database_directory):
        os.system('mkdir -p %s'%database_directory)
    if not os.path.isdir(database_directory):
        raise Exception('found a file instead of directory at path %s'%database_directory)
    if len(os.listdir(database_directory)) > 0:
        raise Exception('database_directory %s is not empty'%database_directory)
        
    temp_dir = 'temp_create_meshdb_%s'%datetime.now().strftime('%Y%m%dT%H%M%SU%f')
    os.system('mkdir -p %s'%temp_dir)
        
    with netCDF4.Dataset(rain_file) as ds, ForcingGrid_DBManager(database_directory, mode='w', verbose=1) as mesh_db:
        n_cells = ds.dimensions['n_meshes'].size
        dates_juliandays = ds.variables['dates'][:]
        for ii in range(len(dates_juliandays)):
            datetime_loc = datetime(1950,1,1)+timedelta(dates_juliandays[ii])
            file_out = '%s/rain_%d.nc'%(temp_dir, ii)
            with netCDF4.Dataset(file_out, mode='w') as ds_out:
                ds_out.createDimension('n_meshes', n_cells)
                var = ds_out.createVariable('rain', 'f4', ('n_meshes', ), zlib=True, complevel=4, shuffle=True)
                var[:] = ds.variables['rain'][:,ii]
                ds_out.setncattr('date_measurement', datetime_loc.strftime('%Y-%m-%dT%H:%M:%S.%f'))
                ds_out.setncattr('date_measurement_CNES_jday', dates_juliandays[ii])
                ds_out.setncattr('date_created', datetime.now().strftime('%Y-%m-%dT%H:%M:%S.%f'))
            mesh_db.add({'file_path': file_out, 'data_type': 'rain', 'date_data': datetime_loc, \
                'product_type': 'analysis', 'grid_status': 'complete'})
                
    os.system('rm -R -f %s'%temp_dir)
    

if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to create a forcing on-mesh database from a rain_file containing rain(n_meshes, dates)')
    parser.add_argument("rain_netcdf_file", type=str, help="rain_netcdf_file")
    parser.add_argument("db_folder", type=str, help="db_folder: path to a non-exisitng folder where database will be created")
    args = parser.parse_args()
    
    fol = os.path.realpath(args.db_folder)
    
    #Exemple de base de donnée remplie avec des nombres aléatoires
    rain_netcdf_to_forcing_database(args.rain_netcdf_file, fol)


    
    
    
    
