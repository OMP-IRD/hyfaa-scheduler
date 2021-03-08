
To run MGB simulation:
0) Copy (including subdirectories) the fortran code (git: in src/MGB/code_MGB_fortran) to a new directory (let's name it $code_dir)
    => to avoid compliling/running the program in git tracked directories because it may leave tempory files that we don't want to track
1) Go to directory and compile using command 'make' ('make clean' command cleans directory of temporary compiled files as well as program)
2) create output directory at any adress (let's name it $output_dir)
3) Fill main parameter file (let's call it $info_mgb_file e.g. infoMGB.yaml), point to necessary static files :
    => niger configuration in src/MGB/Inputs/MGB_inputs_config_niger
    => amazon configuration in src/MGB/Inputs/MGB_inputs_config_amazon
4) Copy / build precipitation file as well as obsevation and substitute stations file if needed
    => those files are not included in the git repository because they are too large
    => small test versions (only for reduced number of iterations) of those files will be included later in the git repository for automated tests

5) Run program :
    => $code_dir/bin/mgb_iph.x $info_mgb_file $output_dir

WARNING:
- The MGB code must be compiled with gfortran, otherwise record length of binary files may not be handled properly => to be corrected in the future by using the netCDF standard

Some librairies are necessary:
- datetime-fortran (datetime_module) => installation instructions on https://github.com/wavebitscientific/datetime-fortran
  => 
    git clone https://github.com/wavebitscientific/datetime-fortran
      => alread contained in this folder
    cd datetime-fortran
    mkdir build
    cd build
    cmake ..
    make
    ctest
  =>
    then change paths ...datetime-fortran/build... from the MGB code Makefile to the path you put the library in
- lapack,blas,netcdf
  => 
    simply use "sudo apt install liblapack-dev", "sudo apt install libblas-dev", "sudo apt install libnetcdf-dev", "sudo apt install libnetcdff-dev" 
    
