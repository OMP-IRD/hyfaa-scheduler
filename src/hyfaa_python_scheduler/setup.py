#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os, sys, subprocess, shutil
import tempfile
from setuptools import setup, find_packages





def find_shell(dir_in, expr=None, is_dir=False, is_file=False):
    if not os.path.exists(dir_in):
        raise Exception('directory %s not found'%dir_in)
    dir_in = os.path.realpath(dir_in)
    cmd = ['find', dir_in]
    if is_dir and is_file:
        raise Exception('is_dir and is_file cannot both be true')
    elif is_dir:
        cmd += ['-type', 'd']
    elif is_file:
        cmd += ['-type', 'f']
    if expr is not None:
        cmd += ['-iname', expr]
    list_find = subprocess.check_output(cmd).decode("utf-8").split('\n')[0:-1]
    for el in list_find:
        if not os.path.exists(el):
            raise Exception('find returned non existing path %s'%el)
        if is_dir and (not os.path.isdir(el)):
            raise Exception('find returned non directory path %s'%el)
        elif (not is_dir) and os.path.isdir(el):
            raise Exception('find returned non file path %s'%el)
    return list_find
    
    
def is_program(filename):
    with open(filename) as ds:
        lines = ds.readlines()
    return any([('if __name__ == "__main__":' in line) or ("if __name__ == '__main__':" in line) for line in lines])



if __name__ == '__main__':
    
    package_name='hyfaa'
    version = '4.0.0'

    hyfaa_src_dir = os.path.abspath(os.environ['hyfaa_src_dir'])
    scripts = [el.replace(hyfaa_src_dir + 'hyfaa_python_scheduler/', '') for el in find_shell(package_name, expr='*.py', is_file=True) if is_program(el)]
    print('Scripts:\n%s\n'%('\n'.join(sorted([os.path.basename(el) for el in scripts]))))

    setup(name=package_name,
        version=version,
        description='HYFAA scheduler by Magellium',
        long_description=open(os.path.join(os.path.dirname(hyfaa_src_dir), 'README.md')).read(),
        classifiers=[
            'Development Status :: 1',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3.7',
            'Topic :: Software Development :: Libraries :: Python Modules',
            'License :: GNU'
            ],
        keywords='hydrology',
        author='Magellium',
        author_email='xxxxx@magellium.fr',
        packages=find_packages(),
        include_package_data=False,
        install_requires=[],
        dependency_links=[],
        scripts= scripts,
        zip_safe=False
    )



