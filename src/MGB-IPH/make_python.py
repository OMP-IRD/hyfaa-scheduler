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

import os, sys, subprocess

python_version = sys.version_info.major

def multiple_split(txt, split_list):
    out = [txt]
    for split_item in split_list:
        out_new = []
        for el in out:
            out_new += el.split(split_item)
        out = [el for el in out_new if len(el) > 0]
    return out
    
    
def exclude_between(txt, id_first, id_last):
    ntxt, n1, n2 = len(txt), len(id_first), len(id_last)
    if n1 < 1 or n2 < 1:
        raise Exception('len(id_first) < 1 or len(id_last) < 1')
    if (id_first not in txt) or (id_last not in txt):
        return txt
    ids_valid = [True for el in txt]
    i0 = None
    ii = 0
    within = False
    while(True):
        if not within:
            if ntxt-ii <= n1:
                break
            if txt[ii:ii+n1] == id_first:
                ii = ii+n1
                i0 = ii
                within = True
            else:
                ii += 1
        else:
            if ntxt-ii < n2:
                break
            if txt[ii:ii+n2] == id_first:
                for i2 in range(i0,ii+n2):
                    ids_valid[i2] = False
                ii = ii+n2
                i0 = None
                within = False
            else:
                ii+= 1
    return ''.join([txt[ii] for ii in range(len(txt)) if ids_valid[ii]])
    
    
    
    

def is_subroutine_module(fortran_file):
    fortran_file_uni = fortran_file.lower().split('.')[0]
    if python_version == 3:
        with open(fortran_file, encoding="ISO-8859-1") as ds:
            lines = ds.readlines()
    else:
        with open(fortran_file) as ds:
            lines = ds.readlines()
    lines = [line.replace('\n','').lower().split() for line in lines if len(line.replace('\n','').replace(' ','')) > 5]
    lines = [line for line in lines if line[0] in ['subroutine', 'function', 'module']]
    if (len(lines) == 0) or (len(lines) > 1):
        return False
    if lines[0] == 'module':
        return False
    return lines[0][1].split('(')[0] == fortran_file_uni
    
    
    

def get_necessary_modules_list(fortran_file, fortran_subroutine_modules_uni):
    
    fortran_file_uni = fortran_file.lower().split('.')[0]
    
    if python_version == 3:
        with open(fortran_file, encoding="ISO-8859-1") as ds:
            lines = ds.readlines()
    else:
        with open(fortran_file) as ds:
            lines = ds.readlines()
    lines = [line.lower().split('!')[0] for line in lines]
    
    module_files = []
    for line in lines:
        if len(line.replace(' ','').replace('\n','')) < 5:
            continue
        line_split = line.split()
        #modules and subroutine modules
        if len(line_split) > 1:
            if line_split[0] in ['use', 'call']:
                module_file = line_split[1].split(',')[0].split('!')[0].split('(')[0]
                if line_split[0] == 'call':
                    if module_file in fortran_subroutine_modules_uni:
                        module_files.append(module_file)
                else:
                    module_files.append(module_file)
        #function modules
        line_split = multiple_split(exclude_between(exclude_between(line, "'", "'"), '"', '"'), ['(',',','+','-','*','/','*','=',' '])
        line_split = [el for el in line_split if len(el) > 2]
        for el in line_split:
            if el == fortran_file_uni:
                continue
            if el in fortran_subroutine_modules_uni:
                module_files.append(el)
    
    module_files = sorted(list(set(module_files)))
    
    return module_files
    
    

    
    
    
def get_fortran_compile_list(compilation_list, add_file, modules_available, modules_infinite_loop_uni, fortran_subroutine_modules_uni, exceptions=['datetime_module', 'netcdf']):
    
    add_file_uni = add_file.lower().split('.')[0]
    if add_file_uni in modules_infinite_loop_uni:
        print(add_file_uni, modules_infinite_loop_uni)
        raise Exception('use module loop detected')
    
    compilation_list_uni = [el.lower().split('.')[0] for el in compilation_list]
    if add_file_uni in compilation_list_uni:
        return compilation_list
        
    modules_available_uni_dict = {el.lower().split('.')[0]: el for el in modules_available}
    necessary_modules_uni = get_necessary_modules_list(add_file, fortran_subroutine_modules_uni)
    for el in necessary_modules_uni:
        if el in exceptions:
            continue
        if el not in modules_available_uni_dict.keys():
            raise Exception('module %s required by file %s is not available'%(el, add_file))
        compilation_list = get_fortran_compile_list(compilation_list, modules_available_uni_dict[el], modules_available, \
            modules_infinite_loop_uni + [add_file_uni], fortran_subroutine_modules_uni, exceptions=exceptions)

    
    compilation_list.insert(0, add_file)
    
    return compilation_list
    
    
 
def main(debug_mode=False):
    
    make_clean()
    
    main_program_file = '1main.f90'
    sufixes = ['.f90', '.F90']
    
    for fol in ['obj', 'bin']:
        if not os.path.exists(fol):
            os.system('mkdir -p %s'%fol)
        if not os.path.exists('%s/.gitignore'%fol):
            with open('%s/.gitignore'%fol, mode='w') as ds:
                ds.write('# Ignore everything in this directory\n*\n# Except this file\n!.gitignore\n')
    
    include_dirs = [el for el in ['/usr/local/include', '/usr/include', '%s/build/include'%os.environ['datetime_fortran_path'], \
        '/softs/rh7/netcdf/4.4.1/include'] if os.path.exists(el)]
    lib_dirs = [el for el in ['%s/build/lib'%os.environ['datetime_fortran_path'], '/softs/rh7/netcdf/4.4.1/lib'] if os.path.exists(el)]
    
    compiler = 'gfortran'
    if debug_mode:
        cflags = '-g -fcheck=all -fbacktrace -DUNIX -DGFORTRAN -ffree-line-length-none'
    else:
        cflags = '-O3 -ffree-line-length-none'
        
    fortran_files = os.listdir('.')
    fortran_files = [el for el in fortran_files if any([sufix in el for sufix in sufixes])]
    fortran_subroutine_modules_uni = [el.lower().split('.')[0] for el in fortran_files if is_subroutine_module(el)]

    if main_program_file not in fortran_files:
        raise Exception('main program file %s missing'%main_program_file)
    compilation_list = get_fortran_compile_list([], main_program_file, fortran_files, [], fortran_subroutine_modules_uni)[::-1]

    
    lines = ['FC=%s'%compiler,'CFLAGS=%s'%cflags, 'SRC=%s'%(' \\\n'.join(compilation_list))]
    lines += ['OBJ=$(addprefix obj/, $(patsubst %.F, %.o, $(patsubst %.F90, %.o, $(patsubst %.f90, %.o, $(SRC)))))']
    for sufix in sufixes:
        lines += ['obj/%%.o: %%%s\n\t$(FC) -cpp $(CFLAGS) -c $< -o $@ %s'%(sufix, ' '.join(['-I%s'%el for el in include_dirs]))]
    lines += ['bin/mgb-iph: $(OBJ)\n\t$(FC) -o bin/mgb_iph $(OBJ) %s -ldatetime -lnetcdff'%(' '.join(['-L%s'%el for el in lib_dirs]))]
    lines += ['clean:\n\trm obj/*.o *.mod bin/mgb_iph']

        
    with open('Makefile', mode='w') as ds:
        ds.write('\n%s\n'%('\n\n'.join(lines)))
        
    subprocess.check_call('make', shell=True)
    if not debug_mode:
        os.unlink('Makefile')


def make_clean():
    os.system('rm -f Makefile; rm -f bin/mgb_iph; rm -f *.mod; rm -f obj/*.o')
    
    
if __name__ == '__main__':
    
    debug_mode = False
    clean_mode = False
    if len(sys.argv) == 2:
        if sys.argv[1] == 'clean':
            clean_mode = True
        if sys.argv[1] == 'debug':
            debug_mode = True
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    if clean_mode:
        make_clean()
    else:
        main(debug_mode=debug_mode)

    
