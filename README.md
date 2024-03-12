
# HYFAA V4 : Global Hydrological Operational Forecast and Assimilated System version 4

Install markdown reader browser extension, enable access to urls, and open this file to view it with a nicer rendering.

## About HYFAA

HYFAA is a python scheduler for operational hydrological forecasting. To achieve this, it combines :
- a hydrological simulation model : MGB-IPH
- retrieval routines that gather the necessary forcing and assimilation data into organised databases
- a main operational scheduler that configures MGB-IPH simulations and handles assimilation routines
- a post-processing routine to make easy-to-use files containing main variables

HYFAA is composed of :
- python HYFAA code provided in the __src/scheduler_python__ folder
- fortran MGB-IPH code provided in the __src/MGB-IPH__ folder called by the python HYFAA code

Hardware requirements :
- CPU : no minimum requirements (at least 1 ^^)
- RAM : 8GB minimum, or more for large ensemble calculations (post-processing step requires it)
- Disk space : it is recommended to have 1GB + 1MB * (ensemble_size * number_of_days) i.e. ~ 185GB space for a 10 year simulation with an ensemble size of 50



## Install HYFAA

### With Docker

Make docker image with `make_docker.sh`

### On a linux PC without docker

1. `apt-get install build-essential gfortran cmake libnetcdf-dev libnetcdff-dev`
2. `pip install numpy numba scipy netCDF4 pyyaml progress pandas geopandas pytest requests SALib ftputil`

#### Standard install

3. Run `./install.sh`

NB:
- This will store `hyfaa` python modules in your default python site-package which is already in your import paths, so no action required.
- `mgb_iph` script will be stored in `/usr/local/bin` (requires root privileges) which should be in your $PATH, so no action should be necessary.

#### Local install
Use for instance if you do not have root privileges on your machine

3. Run `./install.sh ${mgb_iph_install_dir}`, `${mgb_iph_install_dir}` being any directory (must not exist prior to installation).
4. You will need to add paths to $PATH and $PYTHONPATH to directories within `${mgb_iph_install_dir}` that contain `hyfaa` and `mgb_iph` executables:
```
export PATH=${mgb_iph_install_dir}:${mgb_iph_install_dir}/bin:$PATH
export PYTHONPATH=$(find ${mgb_iph_install_dir} -type d -iname 'site-packages'):$PYTHONPATH
```

NB:
- The export commands necessary will be shown at the end of the install script.
- __To avoid entering those lines everytime you open a new terminal, simply add them to your ~/.bashrc__

### Standard linux install on CNES cluster

#### Method 1
This is the preferred method as it installs an independent python environment and therefore allows more flexibility to add other libraries in the future and prevents problems with changes on HAL python modules.
- `module purge` to avoid module conflicts
- `module load cmake netcdf/4.4.1 conda`
- `conda create -n hyfaa_env python=3.7`
- `conda activate hyfaa_env`
- Follow __Local install__ from linux PC without docker installation from step 2

#### Method 2
This solely relies on HAL modules for python environment, rendering installation easier but changes on HAL python module may cause installation to fail in the future, and SALib library is not available.

- `module purge` to avoid module conflicts
- `module load cmake netcdf/4.4.1 python`
- Follow __Local install__ from linux PC without docker installation from step 3


### Windows or Mac

Use docker (best option), or use a unix virtual machine.


## Use HYFAA

- choose a configuration folder in `work_configurations`
- follow the `README.md` inside of the configuration folder to download input_data folder (hydrological static data configuration) and initialized `databases` (may be optional)


#### With Docker

- edit `run_docker.sh` to
  - mount the configuration folder chosen as `/work`,
  - adjust the hydroweb credentials
and launch `run_docker.sh`

#### Without docker

- go to the chosen configuration folder and launch `run.sh`

#### Using PBS on CNES cluster

- go to the chosen configuration folder and launch `./run_pbs.py`

__WARNING:__ So that modules and paths, pythonpaths are set on the node, you must either add them to your ~/.bashrc, or to the `run.sh` script in the configuration folder.

__NB:__ You can use the `--pbs_name` option to set your job name; `hyfaa` by default.
