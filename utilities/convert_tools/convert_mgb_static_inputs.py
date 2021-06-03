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

import os, sys, shutil
from datetime import datetime, timedelta
import numpy as np
import netCDF4


def search_str(str_list, search_list, case_sensible=True, single=False, error_none=False):
    
    if not isinstance(search_list, list):
        search_list = [search_list]
        
    if not case_sensible:
        str_match_list = [el for el in str_list if all([search_str.lower() in el.lower() for search_str in search_list])]
    else:
        str_match_list = [el for el in str_list if all([search_str in el for search_str in search_list])]
        
    if error_none and len(str_match_list) == 0:
        raise Exception('no match for search_list [%s]'%(', '.join(search_list)))
        
    if single:
        if len(str_match_list) > 1:
            raise Exception('multiple match for search_list [%s]'%(', '.join(search_list)))
        elif len(str_match_list) == 0:
            return None
        else:
            return str_match_list[0]
            
    return str_match_list
    
    
    
###############################################################
#read main mesh information file (mini.gtp...)
def get_info_mgb_data(input_file, verbose=1):
    """Read mini.gtp file and translate each variable into a dictionnary entry containing 
    data_type,longname,dimensions,data,etc... information so it can be written in a netCDF file"""
    
    if verbose > 0:
        print('Reading main mesh information (mini.gtp) file %s'%input_file)
    
    with open(input_file) as ds:
        lines = ds.readlines()
    header = lines[0].split()
    nitems = len(header)
    lines = [line.split() for line in lines[1:]]
    wrong_size_line_indexes = [ii for ii in range(len(lines)) if len(lines[ii])!=nitems]
    if wrong_size_line_indexes:
        print(wrong_size_line_indexes)
        raise Exception('Above line indexes have number of element different from header...')
        
    dico = {'attributes': dict(), 'variables': dict()}
    dico['variables'] = {"cell_id_original": {'column_name': 'CatID', 'data_type': 'i4', 'longname': 'cell/mini-basin ID in the original build'}, \
    "longitude_center": {'column_name': 'Xcen', 'data_type': 'f4', 'longname': 'longitude of the cell center', 'unit': 'degrees east (from greenwich)'}, \
    "latitude_center": {'column_name': 'Ycen', 'data_type': 'f4', 'longname': 'latitude of the cell center', 'unit': 'degrees north (from equator)'}, \
    "cell_basin_id": {'column_name': 'Sub', 'data_type': 'i4', 'longname': 'cell basin ID : gives the ID of the basin to which the cell/mini-basin belongs'}, \
    "cell_area": {'column_name': 'Area_(km2)', 'data_type': 'f4', 'longname': 'Cell/mini-basin area', 'unit': 'km^2'}, \
    "total_drained_area": {'column_name': 'AreaM_(km2)', 'data_type': 'f4', 'longname': 'Total area drained by cell/mini-basin', 'unit': 'km^2'}, \
    "main_river_length": {'column_name': 'Ltr_(km)', 'data_type': 'f4', 'longname': 'Length of the main river within cell/mini-basin', 'unit': 'km'}, \
    "main_river_slope": {'column_name': 'Str_(m/km)', 'data_type': 'f4', 'longname': 'Slope of the main river within cell/mini-basin', 'unit': 'km'}, \
    "longest_river_length": {'column_name': 'Lrl_(km)', 'data_type': 'f4', 'longname': 'Length of the longest river within cell/mini-basin', 'unit': 'km'}, \
    "longest_river_slope": {'column_name': 'Srl_(m/km)', 'data_type': 'f4', 'longname': 'Slope of the longest river within cell/mini-basin', 'unit': 'km'}, \
    "cell_id_upstream": {'column_name': 'MiniJus', 'data_type': 'i4', 'longname': 'Upstream cell/mini-basin ID'}, \
    "hydrological_order": {'column_name': 'Ordem', 'data_type': 'i4', 'longname': 'Hydrological order (Strahler)'}, \
    "hydrodynamic_module": {'column_name': 'Hdr', 'data_type': 'i4', 'longname': 'Hydrological module'}, \
    "river_width": {'column_name': 'Width', 'data_type': 'f4', 'longname': 'River width'}, \
    "river_depth": {'column_name': 'Depth', 'data_type': 'f8', 'longname': 'River depth'}, \
    "manning_coefficient": {'column_name': 'NMAN', 'data_type': 'f4', 'longname': 'Manning roughness coefficient'}}
    
    for key, dico_loc in dico['variables'].items():
        index = header.index(dico_loc['column_name'])
        dico_loc['dimensions'] = ['n_cells']
        dico_loc['data'] = np.array([line[index] for line in lines], dtype=dico_loc['data_type'])
        
    nhrus = len([el for el in header if 'BLC_' in el])
    indexes = [header.index('BLC_%02d'%ii) for ii in range(1,nhrus+1)]
    dico['variables']['hydrological_response_unit_coefficients'] = {'column_name': 'BCL_*', 'data_type': 'f4', 'dimensions': ['n_soil_types', 'n_cells'], \
        'longname': 'Hydrological response unit (HRU) coefficients : represents the percentage of each soil type within the cell, described by the hydrological_response_unit_climatology_* variables', \
        'data': np.array([[line[index] for index in indexes] for line in lines], dtype='f4').T}
        
    basin_ids_set = set(dico['variables']['cell_basin_id']['data'])
    dico['variables']['basin_id'] = {'data_type': 'i4', 'longname': 'basin IDs', 'dimensions': ['n_basins'], 'data': np.arange(1, max(basin_ids_set)+1, dtype='i4')}
    if set(dico['variables']['basin_id']['data']) != basin_ids_set:
        print('All basin IDS: %s'%(','.join(['%d'%ii for ii in sorted(list(basin_ids_set))])))
        print('WARNING: some basins are not represented by any cells: %s'%(', '.join(['%d'%el for el in sorted(list(set(dico['variables']['basin_id']['data']) - basin_ids_set))])))
        
    if verbose > 0:
        print('  -> Main mesh information (mini.gtp) file successfully read')
    
    return dico
###############################################################
    
    
    
    
    
###############################################################
#read vegetation parameter climatology file (ALBIAF...)
def get_data_from_hru_region_lines(lines, verbose=1):
    """Read a block of lines extracted from a hydrological response unit file containing albedo, 
    leaf_area_index, trees_height and superficial_resistance climatology information for a single climatic region, 
    as well as soil_type names. Translate it into dictionnary entries."""
    dico = {'soil_types': None, 'data': dict()}
    ii_line = 0
    while('albedo' not in lines[ii_line]):
        ii_line += 1
    while(ii_line<len(lines)-1):
        line = lines[ii_line]
        if 'albedo' in line:
            key = 'albedo'
        elif 'LAI  !leaf area index' in line:
            key = 'leaf_area_index'
        elif 'Z !trees heigth' in line:
            key = 'trees_height'
        elif 'SR  !superficial resistance' in line:
            key = 'superficial_resistance'
        else:
            break
        lines_block = []
        while(ii_line<len(lines)-1):
            ii_line += 1
            line_splitted = lines[ii_line].split()
            if len(line_splitted) == 13:
                lines_block.append(line_splitted)
            else:
                break
        if len(lines_block) < 2:
            raise Exception('retrieved block with less than 2 lines')
        if lines_block[0][0] != 'use':
            raise Exception('First line of block should contain "use" in first column')
        lines_block = lines_block[1:]
        if dico['soil_types'] is None:
            dico['soil_types'] = np.array([line[0] for line in lines_block])
        else:
            assert all(dico['soil_types'] == np.array([line[0] for line in lines_block]))
        dico['data'][key] = np.array([line[1:] for line in lines_block], dtype='f4')
        
    missing_keys = sorted(list(set(['albedo', 'leaf_area_index', 'trees_height', 'superficial_resistance']) - set(list(dico['data'].keys()))))
    if len(missing_keys) > 0:
        raise Exception('Could not recover the following vegetation variables: %s'%(', '.join(missing_keys)))
            
    return dico
    
    
    
def read_region_file(region_file, verbose=1):
    """Read region file and translate region_id into dictionnary entries"""
    if verbose > 0:
        print('Reading region information file %s'%region_file)
    with open(region_file) as ds:
        lines = ds.readlines()
    lines = [line.replace('\n','').split() for line in lines[1:] if len(line.replace('\n','').split()) == 2]
    assert all(np.array([int(line[0]) for line in lines], dtype='i4') == np.arange(len(lines), dtype='i4')+1), 'cell ids from region file are not monotonous fr 1 to n_cells'
    if verbose > 0:
        print('  -> Region information file successfully read')
    return np.array([int(line[1]) for line in lines], dtype='i4')



def get_hru_data(hru_file, cell_basin_ids, region_file=None, congo_region_mode=False, verbose=1):
    """Read hydrological response unit file (and region file if it exists) and translate each variable into a dictionnary entry containing 
    data_type,longname,dimensions,data,etc... information so it can be written in a netCDF file"""
    
    n_cells = len(cell_basin_ids)
    
    if verbose > 0:
        print('Reading Hydrological Response Unit information file %s'%hru_file)
    
    #read hru file
    with open(hru_file) as ds:
        lines = ds.readlines()
    lines = [line.replace('\n','') for line in lines if len(line.replace(' ','').replace('\n','')) > 0]
    
    #handle case when HRU file is composed of different regions
    if region_file is not None:
        #read region file
        region_ids = read_region_file(region_file, verbose=verbose)
        assert len(region_ids) == n_cells, 'length of region ids do not match number of cells'
        
        #if congo region mode, force region ids in ARUWIMI region
        if congo_region_mode:
            print('Applying congo region mode')
            for ii in range(n_cells):
                if cell_basin_ids[ii] == 14:
                    region_ids[ii] = 4
                elif cell_basin_ids[ii] == 13:
                    region_ids[ii] = 6
                    
    else:
        assert not congo_region_mode
        region_ids = None
    if not any(['region' in line.lower() for line in lines]):
        lines.insert(0, 'region 1')
        if region_ids is not None:
            if any(region_ids != 1):
                raise Exception('single region encountered : region_id values in region file should be all "1" but this is not the case')
        else:
            region_ids = np.ones(n_cells, dtype='i4')
    else:
        if region_ids is None:
            raise Exception('multiple regions encountered in hru_file : region_file must be present')
    
    n_regions = sum(['region' in line.lower() for line in lines])
    assert set(list(range(1,n_regions+1))) == set(region_ids), 'some regions are not represented by any cells'
    
    dico = {'attributes': dict(), 'variables': dict()}
    dico['variables'] = {"cell_region_id": {'data_type': 'i4', 'longname': 'cell region id : to which region the celle belongs', 'dimensions': ['n_cells'], 'data': region_ids}, \
        "albedo_climatology": {'data_type': 'f4', 'longname': 'Albedo function of month and soil type (vegetation varies with seasons)', 'dimensions': ['n_months_year', 'n_soil_types', 'n_regions']}, \
        "leaf_area_index_climatology": {'data_type': 'f4', 'longname': 'Leaf area index function of month and soil type (vegetation varies with seasons)', 'dimensions': ['n_months_year', 'n_soil_types', 'n_regions']}, \
        "trees_height_climatology": {'data_type': 'f4', 'longname': 'Trees height function of month and soil type (vegetation varies with seasons)', 'dimensions': ['n_months_year', 'n_soil_types', 'n_regions']}, \
        "superficial_resistance_climatology": {'data_type': 'f4', 'longname': 'Superficial resistance function of month and soil type (vegetation varies with seasons)', 'dimensions': ['n_months_year', 'n_soil_types', 'n_regions']}}

    id_lines_region = [ii for ii in range(len(lines)) if 'region' in lines[ii].lower()]
    for index_line, ii_line in enumerate(id_lines_region):
        line_splitted = lines[ii_line].split()
        
        #check that region ids gor from 1 to n_regions in HRU file
        ii = 0
        while(line_splitted[ii].lower() != 'region'):
            ii += 1
        id_region = int(line_splitted[ii+1])
        assert id_region == index_line+1
        
        id_min_block = ii_line+1
        if index_line < len(id_lines_region)-1:
            id_max_block = id_lines_region[index_line+1]
        else:
            id_max_block = len(lines)
        dico_loc = get_data_from_hru_region_lines(lines[id_min_block:id_max_block], verbose=verbose)
        if index_line == 0:
            n_soil_types = len(dico_loc['soil_types'])
            soil_types = dico_loc['soil_types']
            for el in dico_loc['data'].keys():
                dico['variables']['%s_climatology'%el]['data'] = np.empty((12, n_soil_types, n_regions), dtype='f4')
        assert all(soil_types == dico_loc['soil_types'])
        for el in dico_loc['data'].keys():
            assert np.shape(dico_loc['data'][el]) == (n_soil_types, 12)
            dico['variables']['%s_climatology'%el]['data'][:,:,index_line] = dico_loc['data'][el].T
            
    dico['attributes'] = {'soil_types': ','.join(soil_types)}
    if verbose > 0:
        print('  -> Hydrological Response Unit information file successfully read')
            
    return dico
###############################################################




###############################################################
#read climatology forcing data
vars_climato_dict = {'temperature': 'ic_temp', 'relative_humidity': 'ic_urel', 'sunshine': 'ic_sunp', 'wind_speed': 'ic_wind', 'pressure': 'ic_patm'}

def read_climatology_file(climatology_file, n_cells, verbose=1):
    if verbose > 0:
        print('Reading forcing climatology information file %s'%climatology_file)
    with open(climatology_file) as ds:
        lines = ds.readlines()
    lines = [line.replace('\n','').split() for line in lines]
    lines = [line for line in lines if len(line) == 14]
    dico = dict()
    for key, value in vars_climato_dict.items():
        dico[key] = [line[1:] for line in lines if line[0].lower() == value]
    del lines
    for key, lines in dico.items():
        if np.any(np.array([line[0] for line in lines], dtype='i4') != np.arange(n_cells, dtype='i4')+1):
            raise Exception('%s lines do not match cell ids'%key)
        dico[key] = np.array([[float(el) for el in line[1:]] for line in dico[key]], dtype='f4').T
    return dico



#get all climatology information
def get_climatology_data(climato_list, n_cells, verbose=1):
    """Read climatology forcing information file and translate each variable into a dictionnary entry containing 
    data_type,longname,dimensions,data,etc... information so it can be written in a netCDF file"""

    dico_climato_files = dict()
    for el in climato_list:
        els = el.split(',')
        assert len(els) == 2
        if els[0] != 'default':
            els[0] = int(els[0])
        assert els[0] not in dico_climato_files, '%s key encountered multiple times in climato_list'%els[0]
        dico_climato_files[els[0]] = els[1]
    assert 'default' in dico_climato_files
    n_climatology = len(dico_climato_files)
    
    order_climato = ['default'] + sorted([el for el in dico_climato_files if el != 'default'])
    specific_climatology_years = np.array([-1] + order_climato[1:], dtype='i4')
    dico = {'attributes': dict(), 'variables': dict()}
    dico['variables']['specific_climatology_years'] = {'data_type': 'i4', 'longname': 'specific climatology years', 'dimensions': ['n_climatology'], 'data': specific_climatology_years}
    for key in vars_climato_dict.keys():
        dico['variables']["%s_climatology"%key] = {'data_type': 'f4', 'longname': '%s climatology on each cell'%key, 'dimensions': ['n_months_year', 'n_cells', 'n_climatology'], \
            'data': np.zeros((12, n_cells, n_climatology), dtype='f4')}
    
    for i0, el in enumerate(order_climato):
        dico_loc = read_climatology_file(dico_climato_files[el], n_cells, verbose=verbose)
        for key in vars_climato_dict.keys():
            dico['variables']["%s_climatology"%key]['data'][:,:,i0] = dico_loc[key]


    if verbose > 0:
        print('  -> Forcing climatology information file successfully read')
    return dico
###############################################################
    
    
    
    

###############################################################
#read flood plain information
def get_flp_cota_area_data(input_file, n_cells, verbose=1):
    """Read flood plain cota_area file and translate each variable into a dictionnary entry containing 
    data_type,longname,dimensions,data,etc... information so it can be written in a netCDF file"""
    
    if verbose > 0:
        print('Reading flood plain cota_area information file %s'%input_file)
    with open(input_file) as ds:
        lines = ds.readlines()
    lines = [line.replace('\n','').split() for line in lines[1:] if len(line.split()) == 4]
    dico = {'attributes': dict(), 'variables': dict()}
    dico['variables'] = {"flood_plain_points_number": {'data_type': 'i4', 'longname': 'number of vertical points to represent flood plain', 'dimensions': ['n_cells'], 'data': np.zeros(n_cells, dtype='i4')}, \
        "flood_plain_points_lowest_altitude": {'data_type': 'f8', 'longname': 'lowest altitude of flood plain representation', 'dimensions': ['n_cells'], 'data': np.zeros(n_cells, dtype='f8')}, \
        "flood_plain_points_altitude": {'data_type': 'f8', 'longname': 'lowest altitude of flood plain representation', 'dimensions': ['n_cells', 'n_flood_plain_points_max']}, \
        "flood_plain_points_area_flooded": {'data_type': 'f8', 'longname': 'lowest altitude of flood plain representation', 'dimensions': ['n_cells', 'n_flood_plain_points_max']}}
    
    #get number of cota_area points per cell
    old_index = -10
    for line in lines:
        index = int(line[0])-1
        dico['variables']['flood_plain_points_number']['data'][index] += 1
        if index != old_index:
            dico['variables']['flood_plain_points_lowest_altitude']['data'][index] = float(line[1])
        else:
            assert dico['variables']['flood_plain_points_lowest_altitude']['data'][index] == float(line[1])
        old_index = index
    
    n_flood_plain_points_max = max(dico['variables']['flood_plain_points_number']['data'])
    dico['variables']['flood_plain_points_altitude']['data'] = np.zeros((n_cells, n_flood_plain_points_max), dtype='f8')
    dico['variables']['flood_plain_points_area_flooded']['data'] = np.zeros((n_cells, n_flood_plain_points_max), dtype='f8')
    old_index = -10
    for line in lines:
        index = int(line[0])-1
        if index != old_index:
            ii = 0
        else:
            ii += 1
        dico['variables']['flood_plain_points_altitude']['data'][index,ii] = float(line[2])
        dico['variables']['flood_plain_points_area_flooded']['data'][index,ii] = float(line[3])
        old_index = index
    
    if verbose > 0:
        print('  -> Flood plain cota_area information file successfully read')
    return dico


def get_flp_face_data(input_file, verbose=1):
    """Read flood plain face.con file and translate each variable into a dictionnary entry containing 
    data_type,longname,dimensions,data,etc... information so it can be written in a netCDF file"""
    if verbose > 0:
        print('Reading flood plain face.con information file %s'%input_file)

    with open(input_file) as ds:
        lines = ds.readlines()
    n_columns = np.median([len(line.replace('\n','').split()) for line in lines])
    lines = [line.replace('\n','').split() for line in lines if len(line.replace('\n','').split()) == n_columns]
    n_faces = len(lines)
    assert all(np.array([line[0] for line in lines], dtype='i4') == np.arange(n_faces, dtype='i4')+1), 'face ids shoud be monotonously increasing in file %s'%input_file
    dico = {'attributes': dict(), 'variables': dict()}
    if n_columns == 4:
        dico['variables'] = {"face_adjacent_cell_1": {'data_type': 'i4', 'longname': 'face adjacent cell 1', 'dimensions': ['n_faces'], 'data': np.array([line[1] for line in lines], dtype='i4')}, \
                "face_adjacent_cell_2": {'data_type': 'i4', 'longname': 'face adjacent cell 2', 'dimensions': ['n_faces'], 'data': np.array([line[2] for line in lines], dtype='i4')}, \
                "face_dx": {'data_type': 'f8', 'longname': 'face dx', 'dimensions': ['n_faces'], 'data': np.array([line[3] for line in lines], dtype='f8')}}
    elif n_columns == 6:
        dico['variables'] = {"face_adjacent_cell_1": {'data_type': 'i4', 'longname': 'face adjacent cell 1', 'dimensions': ['n_faces'], 'data': np.array([line[1] for line in lines], dtype='i4')}, \
                "face_ymin_1": {'data_type': 'f8', 'longname': 'face ymin', 'dimensions': ['n_faces'], 'data': np.array([line[2] for line in lines], dtype='f8')}, \
                "face_adjacent_cell_2": {'data_type': 'i4', 'longname': 'face adjacent cell 2', 'dimensions': ['n_faces'], 'data': np.array([line[3] for line in lines], dtype='i4')}, \
                "face_ymin_2": {'data_type': 'f8', 'longname': 'face ymin', 'dimensions': ['n_faces'], 'data': np.array([line[4] for line in lines], dtype='f8')}, \
                "face_dx": {'data_type': 'f8', 'longname': 'face dx', 'dimensions': ['n_faces'], 'data': np.array([line[5] for line in lines], dtype='f8')}}
    elif n_columns == 7:
        dico['variables'] = {"face_adjacent_cell_1": {'data_type': 'i4', 'longname': 'face adjacent cell 1', 'dimensions': ['n_faces'], 'data': np.array([line[1] for line in lines], dtype='i4')}, \
                "face_ymin_1": {'data_type': 'f8', 'longname': 'face ymin', 'dimensions': ['n_faces'], 'data': np.array([line[2] for line in lines], dtype='f8')}, \
                "face_adjacent_cell_2": {'data_type': 'i4', 'longname': 'face adjacent cell 2', 'dimensions': ['n_faces'], 'data': np.array([line[3] for line in lines], dtype='i4')}, \
                "face_ymin_2": {'data_type': 'f8', 'longname': 'face ymin', 'dimensions': ['n_faces'], 'data': np.array([line[4] for line in lines], dtype='f8')}, \
                "face_dx": {'data_type': 'f8', 'longname': 'face dx', 'dimensions': ['n_faces'], 'data': np.array([line[5] for line in lines], dtype='f8')}, \
                "face_bflow": {'data_type': 'f8', 'longname': 'face bflow', 'dimensions': ['n_faces'], 'data': np.array([line[6] for line in lines], dtype='f8')}}
    else:
        raise Exception('number of columns in file %s must be in [4,6,7]'%input_file)
    
    if verbose > 0:
        print('  -> Flood plain face.con information file successfully read')
    return dico
    

def get_delta_data(delta_list, n_cells, verbose=1):
    """Read flood plain delta files and translate each variable into a dictionnary entry containing 
    data_type,longname,dimensions,data,etc... information so it can be written in a netCDF file"""
    if verbose > 0:
        print('Reading flood plain delta information')
    n_deltas = len(delta_list)
    dico = {'attributes': dict(), 'variables': dict()}
    dico['attributes'] = {'delta_names': []}
    dico['variables'] = {"cell_delta_id": {'data_type': 'i4', 'longname': 'cell delta id or 0 if cell does not belong to a delta.', 'dimensions': ['n_cells', 'n_deltas'], 'data': np.zeros((n_cells, n_deltas), dtype='i4')}, \
        'delta_outlet_value': {'data_type': 'f4', 'longname': 'delta outlet value', 'dimensions': ['n_deltas'], 'data': -np.ones(n_deltas, dtype='f4')}}
        
    for i_delta, delta_info in enumerate(delta_list):
        dico['attributes']['delta_names'].append(os.path.basename(delta_info['file_path']).split('.')[0].lower().replace('flag','').replace('delta','').rstrip('_').lstrip('_'))
        with open(delta_info['file_path']) as ds:
            lines = ds.readlines()
        ids_lines = [ii for ii in range(len(lines)) if len(lines[ii].replace('\n','').replace(' ','')) > 0]
        assert ids_lines == list(range(len(ids_lines))), 'there are some holes in delta file %s'%delta_info['file_path']
        for ii in ids_lines:
            if int(lines[ii]) == 1:
                dico['variables']['cell_delta_id']['data'][ii,i_delta] = 1
        if 'outlet_value' in delta_info:
            if delta_info['outlet_value'] is not None:
                dico['variables']['delta_outlet_value']['data'][i_delta] = delta_info['outlet_value']
    dico['attributes']['delta_names'] = ','.join(dico['attributes']['delta_names'])
    if verbose > 0:
        print('  -> Flood plain delta information successfully read')
    return dico


def get_outlet_data(outlet_list, n_faces, verbose=1):
    """Get outlet data and translate each variable into a dictionnary entry containing 
    data_type,longname,dimensions,data,etc... information so it can be written in a netCDF file"""
    if verbose > 0:
        print('Getting flood plain specific outlet information')
    n_outlets = len(outlet_list)

    dico = {'attributes': dict(), 'variables': dict()}
    dico['attributes'] = {'specific_outlet_names': []}
    dico['variables'] = {'specific_outlet_face_index': {'data_type': 'i4', 'longname': 'delta outlet face index', 'dimensions': ['n_specific_outlets'], 'data': np.zeros(n_outlets, dtype='i4')}, \
        'specific_outlet_value': {'data_type': 'f4', 'longname': 'delta outlet value', 'dimensions': ['n_specific_outlets'], 'data': np.zeros(n_outlets, dtype='f4')}}
        
    for i_outlet, outlet_info in enumerate(outlet_list):
        dico['attributes']['specific_outlet_names'].append(outlet_info['outlet_name'])
        dico['variables']['specific_outlet_face_index']['data'][i_outlet] = outlet_info['outlet_face_index']
        dico['variables']['specific_outlet_value']['data'][i_outlet] = outlet_info['outlet_value']
    dico['attributes']['specific_outlet_names'] = ','.join(dico['attributes']['specific_outlet_names'])
    if verbose > 0:
        print('  -> Flood plain specific outlet information successfully read')
    return dico
###############################################################



###############################################################
#read calibration file
def get_calibration_information(calibration_file, n_basins, n_soil_types, verbose=1):
    """Read calibration (paruso.cal) file and translate each variable into a dictionnary entry containing 
    data_type,longname,dimensions,data,etc... information so it can be written in a netCDF file"""
    if verbose > 0:
        print('Reading calibration (paruso.cal) information file %s'%calibration_file)

    with open(calibration_file) as ds:
        lines = ds.readlines()
    lines = [line.replace('\n','').split() for line in lines if len(line.replace(' ','').replace('\n','')) > 0]
    
    dico = {'attributes': dict(), 'variables': dict()}
    dico['variables'] = {"calibration_soil_max_water_capacity": {'data_type': 'f4', 'longname': 'Calibration: soil max water capacity', 'dimensions': ['n_soil_types', 'n_basins'], 'data': np.zeros((n_soil_types, n_basins), dtype='f4')}, \
        "calibration_soil_arno_model_b": {'data_type': 'f4', 'longname': 'Calibration: parameter B for ARNO model', 'dimensions': ['n_soil_types', 'n_basins'], 'data': np.zeros((n_soil_types, n_basins), dtype='f4')}, \
        "calibration_soil_groundwater_flow_coefficient": {'data_type': 'f4', 'longname': 'Calibration: groundwater flow coefficient', 'dimensions': ['n_soil_types', 'n_basins'], 'data': np.zeros((n_soil_types, n_basins), dtype='f4')}, \
        "calibration_soil_subsurface_flow_coefficient": {'data_type': 'f4', 'longname': 'Calibration: subsurface flow coefficient', 'dimensions': ['n_soil_types', 'n_basins'], 'data': np.zeros((n_soil_types, n_basins), dtype='f4')}, \
        "calibration_soil_arno_model_lambda": {'data_type': 'f4', 'longname': 'Calibration: parameter LAMBDA for ARNO model', 'dimensions': ['n_soil_types', 'n_basins'], 'data': np.zeros((n_soil_types, n_basins), dtype='f4')}, \
        "calibration_soil_capillarity_depth": {'data_type': 'f4', 'longname': 'Calibration: water depth of water ascension by capillarity', 'dimensions': ['n_soil_types', 'n_basins'], 'data': np.zeros((n_soil_types, n_basins), dtype='f4')}, \
        "calibration_soil_capillarity_limit": {'data_type': 'f4', 'longname': 'Calibration: limit to start ascending capilar flux', 'dimensions': ['n_soil_types', 'n_basins'], 'data': np.zeros((n_soil_types, n_basins), dtype='f4')}, \
        "calibration_routing_cs": {'data_type': 'f4', 'longname': 'Calibration routing parameter CS', 'dimensions': ['n_basins'], 'data': np.zeros(n_basins, dtype='f4')}, \
        "calibration_routing_ci": {'data_type': 'f4', 'longname': 'Calibration routing parameter CI', 'dimensions': ['n_basins'], 'data': np.zeros(n_basins, dtype='f4')}, \
        "calibration_routing_cb": {'data_type': 'f4', 'longname': 'Calibration routing parameter CB', 'dimensions': ['n_basins'], 'data': np.zeros(n_basins, dtype='f4')}, \
        "calibration_routing_qesp": {'data_type': 'f4', 'longname': 'Calibration routing parameter QESP', 'dimensions': ['n_basins'], 'data': np.zeros(n_basins, dtype='f4')}}

    ii = 0
    for i_basin in range(n_basins):
        ii += 2
        for i_soil_type in range(n_soil_types):
            dico['variables']['calibration_soil_max_water_capacity']['data'][i_soil_type, i_basin] = float(lines[ii][1])
            dico['variables']['calibration_soil_arno_model_b']['data'][i_soil_type, i_basin] = float(lines[ii][2])
            dico['variables']['calibration_soil_groundwater_flow_coefficient']['data'][i_soil_type, i_basin] = float(lines[ii][3])
            dico['variables']['calibration_soil_subsurface_flow_coefficient']['data'][i_soil_type, i_basin] = float(lines[ii][4])
            dico['variables']['calibration_soil_arno_model_lambda']['data'][i_soil_type, i_basin] = float(lines[ii][5])
            dico['variables']['calibration_soil_capillarity_depth']['data'][i_soil_type, i_basin] = float(lines[ii][6])
            dico['variables']['calibration_soil_capillarity_limit']['data'][i_soil_type, i_basin] = float(lines[ii][7])
            ii += 1
        dico['variables']['calibration_routing_cs']['data'][i_basin] = float(lines[ii][1])
        ii += 1
        dico['variables']['calibration_routing_ci']['data'][i_basin] = float(lines[ii][1])
        ii += 1
        dico['variables']['calibration_routing_cb']['data'][i_basin] = float(lines[ii][1])
        ii += 1
        dico['variables']['calibration_routing_qesp']['data'][i_basin] = float(lines[ii][1])
        ii += 1

    if verbose > 0:
        print('  -> Calibration (paruso.cal) information file successfully read')
    return dico
###############################################################



    
###############################################################
#overall converter
def convert_cell_info(input_folder, output_file, climato_list=None, delta_list=None, outlet_list=None, congo_region_mode=False, complevel=4, shuffle=None, verbose=1):
    
    filenames = [os.path.join(input_folder, filename) for filename in os.listdir(input_folder)]
    
    if complevel > 0:
        if shuffle is None:
            shuffle = True
        elif shuffle in [True, 1]:
            shuffle = True
        elif shuffle in [False, 0]:
            shuffle = False
        else:
            raise Exception('unknown input for shuffle: %s'%shuffle)
        zlib = True
    else:
        if shuffle is None:
            shuffle = False
        elif shuffle in [True, 1]:
            raise Exception('compression shuffle option incompatible with complevel=0')
        elif shuffle in [False, 0]:
            shuffle = False
        else:
            raise Exception('unknown input for shuffle: %s'%shuffle)
        zlib = False
    
    #main mesh information
    dico = get_info_mgb_data(search_str(filenames, 'mini.gtp', case_sensible=False, single=True, error_none=True), verbose=verbose)
    n_cells, n_basins = len(dico['variables']['cell_id_original']['data']), len(dico['variables']['basin_id']['data'])
    
    #hydrological response unit information (vegetation parameters)
    dico_loc = get_hru_data(search_str(filenames, 'albiaf', case_sensible=False, single=True, error_none=True), dico['variables']['cell_basin_id']['data'], \
        region_file=search_str(filenames, 'region', case_sensible=False, single=True, error_none=False), congo_region_mode=congo_region_mode, verbose=verbose)
    for key in ['attributes', 'variables']:
        dico[key].update(dico_loc[key])
    n_regions, n_soil_types, _ = np.shape(dico['variables']['albedo_climatology']['data'])
    
    #read climatology file
    if climato_list is None:
        climato_list = ['default,%s'%search_str(filenames, 'medias', case_sensible=False, single=True, error_none=True)]
    dico_loc = get_climatology_data(climato_list, n_cells, verbose=verbose)
    for key in ['attributes', 'variables']:
        dico[key].update(dico_loc[key])
    
    #read flood plane cota_area file
    dico_loc = get_flp_cota_area_data(search_str(filenames, 'cota_area', case_sensible=False, single=True, error_none=True), n_cells, verbose=verbose)
    for key in ['attributes', 'variables']:
        dico[key].update(dico_loc[key])
    
    #read flood plane face.con file
    dico_loc = get_flp_face_data(search_str(filenames, 'face.con', case_sensible=False, single=True, error_none=True), verbose=verbose)
    for key in ['attributes', 'variables']:
        dico[key].update(dico_loc[key])
    n_faces = len(dico['variables']['face_dx']['data'])
    
    #read delta files
    if delta_list is not None:
        if len(delta_list) > 0:
            dico_loc = get_delta_data(delta_list, n_cells, verbose=verbose)
            for key in ['attributes', 'variables']:
                dico[key].update(dico_loc[key])
    
    #specific outlets
    if outlet_list is not None:
        if len(outlet_list) > 0:
            dico_loc = get_outlet_data(outlet_list, n_faces, verbose=verbose)
            for key in ['attributes', 'variables']:
                dico[key].update(dico_loc[key])
    #read calibration information
    dico_loc = get_calibration_information(search_str(filenames, 'paruso', case_sensible=False, single=True, error_none=True), n_basins, n_soil_types, verbose=verbose)
    for key in ['attributes', 'variables']:
        dico[key].update(dico_loc[key])
        
    #write output netCDF file
    if verbose > 0:
        print('Data retrieval successful => writing netCDF file...')
    with netCDF4.Dataset(output_file, mode='w') as ds:
        #dimensions
        dimension_dict = dict()
        for key, dico_var in dico['variables'].items():
            for dimension, dimension_size in zip(dico_var['dimensions'], np.shape(dico_var['data'])):
                if dimension in dimension_dict:
                    assert dimension_dict[dimension] == dimension_size
                else:
                    dimension_dict[dimension] = dimension_size
        for dimension, dimension_size in dimension_dict.items():
            ds.createDimension(dimension, dimension_size)
            
        #variables
        for key, dico_var in dico['variables'].items():
            var = ds.createVariable(key, dico_var['data_type'], tuple(dico_var['dimensions']), zlib=zlib, complevel=complevel, shuffle=shuffle)
            var[:] = dico_var['data']
            for param_key in sorted(list(set(list(dico_var.keys())) - set(['data', 'data_type', 'dimensions']))):
                var.setncattr(param_key, dico_var[param_key])
                
        #attributes
        for key, value in dico['attributes'].items():
            ds.setncattr(key, value)
        
        

        



if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to convert MGB code input files to the new (single netCDF) format')
    parser.add_argument("input_folder", type=str, help="input folder containing MGB static data files in old MGB input format")
    parser.add_argument("output_file", type=str, help="output file => will contain MGB static data in new format")
    parser.add_argument("--climato", type=str, action='append', help="Climatology. default,filename for the default climatology file. year,filename for year specific climatology file. " + \
        "If only the default, it is not necessary to specify the file, the program will search for it under 'medias*' name.")
    parser.add_argument("--delta", type=str, action='append', default=[], help="Add deltas using 'file_path(,outlet_value)'. " + \
        "Example for Niger northern delta : flag_northernDelta.txt,500")
    parser.add_argument("--specific_outlet", type=str, action='append', help="Add specific outlets using 'outlet_name,outlet_value'. " + \
        "Example for diaka distribuary : diaka_distribuary,600")
    parser.add_argument("--congo_region_mode", action='store_true', help="activate congo region mode : it will overwrite input regions and set region 4 for basin 14 and region 6 for basin 13 (ARUWIMI zone).")
    parser.add_argument("--verbose", type=int, help="verbosity level, default=1", default=1)
    parser.add_argument("--compression_level", type=int, help="compression level [0:9]. 0 is not compression, 9 is most compressed but longer to read. default is 4", default=4)
    parser.add_argument("--compression_shuffle", type=int, help="1 to activate, 0 to deactivate. Default is 1 if compression_level > 0, 0 if compression_level = 0.")
    args = parser.parse_args()
    
    delta_list = []
    for elem in args.delta:
        elem_s = elem.split(',')
        if len(elem_s) == 1:
            delta_list.append({'file_path': elem})
        elif len(elem_s) == 2:
            delta_list.append({'file_path': elem_s[0], 'outlet_value': float(elem_s[1])})
        else:
            raise Exception('--delta must have syntax file_path(,index_outlet,outlet_value) and therefore either 1 or 2 coma-splitted elements, got %d for entry : %s'%(len(elem_s), elem))
    outlet_list = []
    if args.specific_outlet!=None:
        for elem in args.specific_outlet:
            elem_s = elem.split(',')
            if len(elem_s) == 3:
                outlet_list.append({'outlet_name': elem_s[0], 'outlet_face_index': int(elem_s[1]), 'outlet_value': float(elem_s[2])})
            else:
                raise Exception('--specific_outlet must have syntax outlet_name,outlet_face_index,outlet_value and therefore 3 coma-splitted elements, got %d for entry : %s'%(len(elem_s), elem))
    
    #convert static data : mini.gtp, delta flags, face.con, cota_area.flp as well as climatology files ALBIAF.fix and medias.cru
    convert_cell_info(args.input_folder, args.output_file, climato_list=args.climato, \
        delta_list=delta_list, outlet_list=outlet_list, congo_region_mode=args.congo_region_mode, \
        complevel=args.compression_level, shuffle=args.compression_shuffle, verbose=args.verbose)


