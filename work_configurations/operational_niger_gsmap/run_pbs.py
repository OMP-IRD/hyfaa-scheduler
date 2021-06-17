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


import os, sys, shutil, subprocess
import yaml


def run_pbs(pbs_name=None):
    
    cwd = os.getcwd()

    if pbs_name is None:
        pbs_name = 'hyfaa'
    try:
        shutil.which('qsub')
    except:
        raise Exception('Cannot launch on PBS, qsub command not accessible')
    nprocs = max([yaml.load(open(os.path.join('config', el)))['nprocs'] for el in os.listdir('config') if 'input_' in el and '.yaml' in el])
    memory = max(8000, nprocs*4000)
    txt = 'qsub -N %s -v hyfaa_workdir=%s -l select=1:ncpus=%d:mem=%dmb:os=rh7 -l walltime=100:00:00 run.sh'%(pbs_name, cwd, nprocs, memory)
    print(txt)
    subprocess.check_call(txt, shell=True)


    
if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="run pbs", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--pbs_name", type=str, help="name of PBS job. Using this option activates --pbs option.")
    args = parser.parse_args()
    
    run_pbs(pbs_name=args.pbs_name)


