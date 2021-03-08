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

import os, sys


def main(mgb_folder):
    
    if not os.path.exists('%s/Makefile'%mgb_folder):
        raise Exception('folder %s does not contain a Makefile file'%mgb_folder)
    
    fortran_endings = ['.f90', '.F90']
    fortran_files = set([el for el in os.listdir(mgb_folder) if any([ending in el for ending in fortran_endings])])
    
    with open('%s/Makefile'%mgb_folder) as ds:
        lines = ds.readlines()
    id_line_src = [ii for ii in range(len(lines)) if 'SRC=' in lines[ii]][0]
    id_line_obj = [ii for ii in range(len(lines)) if 'OBJ=' in lines[ii]][0]
    lines = lines[id_line_src:id_line_obj]
    fortran_files_makefile = set()
    for line in lines:
        filename = line.replace('SRC=','').replace('\\','').replace('\t','').replace(' ','').replace('\n','')
        print(filename)
        if len(filename) > 0:
            fortran_files_makefile.add(filename)
    
    missing_files = fortran_files_makefile - fortran_files
    unnecessary_files = fortran_files - fortran_files_makefile
   
    if len(missing_files)+len(unnecessary_files) == 0:
        print('Fortran files in folder %s match Makefile perfectly'%mgb_folder)
    
    if len(missing_files) > 0:
        print('Missing files :\n%s\n'%('\n'.join(sorted(list(missing_files)))))
    
    if len(unnecessary_files) > 0:
        print('Unnecessary files :\n%s\n'%('\n'.join(sorted(list(unnecessary_files)))))


    
if __name__ == '__main__':
    
    if len(sys.argv) != 2:
        raise Exception('program + mgb_folder')
        
    main(sys.argv[1])
    
