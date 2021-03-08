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

""" 
    Author: rjugier
    This module is about "tuning" yaml I/O
"""

import os, sys, re
if sys.version_info.major < 3:
    raise Exception('sorry, only python 3 supported')
import yaml
if int(yaml.__version__.split('.')[0]) >= 5:
    from yaml import load, add_implicit_resolver, add_constructor, Loader, CLoader, dump, CDumper
else:
    from yaml import load, add_implicit_resolver, add_constructor, Loader, dump, Dumper
    CLoader, CDumper = Loader, Dumper

pattern = re.compile( r'^(.*)\$\{(.*)\}(.*)$' )



def pathex_constructor(loader,node):
    precedingPath = loader.construct_scalar(node)
    output_strings = []
    while ('${' in precedingPath) and ('}' in precedingPath):
        precedingPath, envVar, remainingPath = pattern.match(precedingPath).groups()
        output_strings.append(remainingPath)
        output_strings.append(os.environ[envVar])
    output_strings.append(precedingPath)
    return ''.join(output_strings[::-1])


add_implicit_resolver("!pathex", pattern)
add_constructor('!pathex', pathex_constructor)


def load_yaml(filename_or_string_input, env_vars=False, complicated_input=False, string_input=False):
    
    if string_input:
        input_string = filename_or_string_input
    else:
        with open(filename_or_string_input) as descr:
            input_string = descr.read()
    if env_vars or complicated_input:
        dico = load(input_string, Loader=Loader)
    else:
        dico = load(input_string, Loader=CLoader)

    if 'include' in dico:
        if isinstance(dico['include'], list):
            for el in dico['include']:
                if os.path.exists(el):
                    dico.update(load_yaml(el, env_vars=env_vars, complicated_input=complicated_input))
                else:
                    raise Exception('file %s does not exist'%el)
        elif isinstance(dico['include'], dict):
            for el in dico['include']:
                if os.path.exists(dico['include'][el]):
                    dico.update({el: load_yaml(dico['include'][el], env_vars=env_vars, complicated_input=complicated_input)})
                else:
                    raise Exception('file %s does not exist'%dico['include'][el])
        else:
            if os.path.exists(dico['include']):
                dico.update(load_yaml(dico['include'], env_vars=env_vars, complicated_input=complicated_input))
            else:
                raise Exception('file %s does not exist'%dico['include'])
    return dico


def write_dict_to_yaml_file(yaml_dict, filename):
    with open(filename, mode='w') as ds:
        dump(yaml_dict, ds, default_flow_style=False, Dumper=CDumper)
        
        


def write_simple_1level_dict_to_yaml_file(yaml_dict, filename):
    with open(filename, mode='w') as ds:
        for key in yaml_dict:
            value = yaml_dict[key]
            if value is None:
                ds.write('%s: \n'%key)
            elif hasattr(value, 'replace'):
                ds.write("%s: %s\n"%(key, value))
            elif hasattr(value, '__len__'):
                ds.write('%s: \n'%key)
                for el in value:
                    ds.write("  - %s\n"%el)
            else:
                ds.write('%s: %s\n'%(key, value))

                    
                    
                    
