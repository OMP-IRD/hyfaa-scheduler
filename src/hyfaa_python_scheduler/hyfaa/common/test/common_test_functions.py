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
from hyfaa.check_input_parameters.check_input_parameters import get_scheduler_basic_main_parameters_dict, get_scheduler_standard_assimilation_dict


#for tests
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    

def print_utest_message(message):
    print(bcolors.OKBLUE + bcolors.BOLD + message + bcolors.ENDC + bcolors.ENDC)
    
def check_condition_utest(message, condition):
    print_utest_message(message)
    if condition:
        success_message_utest()
    else:
        fail_message_utest()
    
def success_message_utest():
    print(bcolors.OKGREEN + bcolors.BOLD + ' -> Succeeded' + bcolors.ENDC + bcolors.ENDC)
    
    
def fail_message_utest():
    print(bcolors.FAIL + bcolors.BOLD + ' -> Failed' + bcolors.ENDC + bcolors.ENDC)
    raise


def make_hydrostate(output_file, static_assimilation_params=None):
    
    dico_main = get_scheduler_basic_main_parameters_dict()
    dico_assim = get_scheduler_standard_assimilation_dict()
    shutil.copy(dico_main['init_conditions']['hydrological_state_start'], output_file)
    if static_assimilation_params is None:
        static_assimilation_params = []
    if len(static_assimilation_params) > 0:
        with netCDF4.Dataset(output_file, mode='a') as ds_w, netCDF4.Dataset(dico_main['mgb']['static_data_file']) as ds_r:
            for elem in static_assimilation_params:
                if elem in ds_w.variables:
                    raise Exception('variable %s is already present in dummy hydrostate file %s'%(elem, dico_main['init_conditions']['hydrological_state_start']))
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
                
    with netCDF4.Dataset(output_file, mode='a') as ds:
        for elem in ds.variables:
            data_loc = ds.variables[elem][:]
            ds.variables[elem][:] = data_loc*(1.+0.2*(np.random.rand(np.prod(np.shape(data_loc))).reshape(np.shape(data_loc))-0.5))



        

