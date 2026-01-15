#!/bin/bash
# If you already have MGB installed as a standard installation, you just need to configure your
# python interpreter to use this source code when running hyfaa
# This scripts does just this

# This environment variable needs to be exported, since it is used in the setup.py file that does the installation
export hyfaa_src_dir=$(pwd)/src
cd ${hyfaa_src_dir}/hyfaa_python_scheduler
pip install .
cd -