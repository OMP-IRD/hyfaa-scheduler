#install script for HYFAA

set -e

if [ $# -eq 0 ]
then
    standard_install=1
    mgb_iph_install_dir='/usr/local/bin'
else
    standard_install=0
    mgb_iph_install_dir=$(realpath $1)
    if [ -d "${mgb_iph_install_dir}" ]
    then
        echo 'target install directory '${mgb_iph_install_dir}' must not exist prior to running this script'
        exit
    fi
fi


export hyfaa_src_dir=$(realpath $(dirname "$0")/src)

#create a temp directory to compile sources for datetime-fortran
temp_dir=${hyfaa_src_dir}/temp
if [ -d "$temp_dir" ]; then rm -Rf ${temp_dir}; fi
mkdir ${temp_dir}
cd ${temp_dir}

#compile datetime-fortran
export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)
tar -zxf ${hyfaa_src_dir}/install_libs/datetime-fortran.tar.gz
cd datetime-fortran
mkdir build
cd build
cmake ..
make -j 4
ctest
export datetime_fortran_path=${temp_dir}/datetime-fortran
        
#compile MGB-IPH code
${hyfaa_src_dir}/MGB-IPH/make_python.py
#remove temp dir for datetime-fortran
rm -Rf ${temp_dir}

if [ ${standard_install} -eq 0 ]
then
    if ! [ -d "${mgb_iph_install_dir}" ]; then mkdir -p ${mgb_iph_install_dir}; fi
    mv ${hyfaa_src_dir}/MGB-IPH/bin/mgb_iph ${mgb_iph_install_dir}/
else
    echo 'sudo priviledges needed to copy mgb_iph binary to '${mgb_iph_install_dir}
    sudo mv ${hyfaa_src_dir}/MGB-IPH/bin/mgb_iph ${mgb_iph_install_dir}/
fi
#clean MGB-IPH source directory
${hyfaa_src_dir}/MGB-IPH/make_python.py clean

#pip install hyfaa python scheduler
cd ${hyfaa_src_dir}/hyfaa_python_scheduler
if [ ${standard_install} -eq 1 ]
then
    pip install .
else
    pip install --prefix=${mgb_iph_install_dir} .
    echo ''
    echo 'IMPORTANT : HYFAA and mgb_iph have been installed to non standard directory '${mgb_iph_install_dir}
    echo 'The following paths must be added to PATH and PYTHONPATH for HYFAA to work:'
    echo 'export PATH='${mgb_iph_install_dir}':'${mgb_iph_install_dir}'/bin:$PATH'
    echo "export PYTHONPATH=$(find ${mgb_iph_install_dir} -type d -iname 'site-packages')"':$PYTHONPATH'
    echo ''
fi
