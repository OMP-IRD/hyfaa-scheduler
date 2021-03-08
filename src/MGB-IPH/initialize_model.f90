!*********************************************************************************
!
!  SUBROUTINE initialize_model reads yaml input file and calls load_static_data routine
!
!---------------------------------------------------------------------------------
!  Discussion:
! 
!    This routine reads simulation setup patameters
!
!
!	 Allocates and saves global variable: 
!
!
!  	Usage:
!
!    * no subroutine is called in this subroutine
!
!    where
!
!    * no arguments are passed in this subroutine
!
!    uses modules and functions
!
!    * module     vars_main   in      vars_main.f90
!    * subroutine load_static_data in load_static_data.f90
!
!	 opens
!
!    * Opens yaml input file containing simulation setup parameters.
!    * Opens netCDF input file (specified in yaml file) containing static data (mesh information, HRUs, climatology)
!
!    reads
!
!    * Reads yaml input file  containing simulation setup parameters
!    * Reads netCDF input file (specified in yaml file) containing static data (mesh information, HRUs, climatology)
!
!    creates
!
!    * Does not create files
!    
!
!---------------------------------------------------------------------------------
!  Licensing:
!
!    This code is distributed under the ...
!
!  Version/Modified: 
!
!    2015.21.06 - 21 June 2015 (By: Fernando Mainardi Fan) 
!    2019.04.09 => conversion to yaml inputs by Rémi Jugier (Magellium)   
!
!  Authors:
!
!    Original fortran version by Walter Collischonn
!    Present fortran version by:
!    * Walter Collischonn
!    * Rodrigo Cauduro Dias de Paiva
!    * Diogo da Costa Buarque
!    * Paulo Pontes RÃ³genes
!    * Mino  Viana Sorribas
!    * Fernando Mainardi Fan
!    * Juan Martin Bravo 
!    * Rémi Jugier
!
!  Main Reference:
!
!    Walter Collischonn,
!    Modelo de Grandes Bacias - Thesis
!    Porto Alegre, 2001
!    ISBN: XXXXXXXXXXX,
!
!---------------------------------------------------------------------------------
!  Variables and Parameters:
!
!   *Variables declarations and routines calls are all commented below.
!	* All variables are global!?
!
!---------------------------------------------------------------------------------		


subroutine initialize_model(yaml_input_file)
	use global_variables
    use yaml
    use julianday_mod
    use static_data_loader
    
	implicit none
    character(len=400), intent(in):: yaml_input_file
    character(len=200):: hardcoded_configuration
	integer K, IOP
	
    
    class(type_node),pointer :: yaml_dict

    !load yaml parameters (should be a dict)write_hydrological_state
    yaml_dict => yaml_read(trim(adjustl(yaml_input_file)))
        
    !~     call yaml_dump(yaml_dict, 'check.yaml')
    
    !hardcoded configuration : some coefficients/parameters are hardcoded into the MGB code and are specific to the basin
    hardcoded_configuration = dictget_strvalue(yaml_dict, 'hardcoded_configuration')
    if (trim(adjustl(hardcoded_configuration))=='niger') then
        is_congo = .false.
    elseif (trim(adjustl(hardcoded_configuration))=='congo') then
        is_congo = .true.
    else
        print*, 'hardcoded_configuration must be niger or congo'
        call abort()
    end if
    
    ! 0: simulation, 1: calibration
    ! must be 0 for use with scheduler
    icalib = dictget_intvalue(yaml_dict, 'run_mode')
    if (icalib == 0) then
        param_calib_file = dictget_strvalue(yaml_dict, 'calibration_mode_input_file')
    end if
    
    ! hydrodynamic_module
    hydrodynamic_module = dictget_strvalue(yaml_dict, 'hydrodynamic_module')
    
    !main output directory
    output_directory = dictget_strvalue(yaml_dict, 'output_directory')
    if (trim(adjustl(output_directory)) == 'none') then
        print*, 'output directory parameter must be filled i.e. != none'
        call abort()
    end if
    trash_str = adjustr(output_directory)
    if (trash_str(1024:1024) == '/') then
        output_directory = trim(adjustl(trash_str(1:1023)))
    end if
    call system('mkdir -p '//trim(adjustl(output_directory)))
    
    ! start date, number of steps, dt of step, etc...
    IDIA = dictget_intvalue(yaml_dict, 'day')
    IMES = dictget_intvalue(yaml_dict, 'month')
    IANO = dictget_intvalue(yaml_dict, 'year')
    IDINI=floor(julianday(IANO,IMES,IDIA)) 
    HORAINI = dictget_intvalue(yaml_dict, 'hour')
    NT = dictget_intvalue(yaml_dict, 'nt')
    DTP = dictget_realvalue(yaml_dict, 'dt')
    alpha = dictget_dpvalue(yaml_dict, 'alpha')
    simulation_start_date = julianday(IANO,IMES,IDIA,HORAINI)
    
    !dynamic forcing information
    forcing_data_file = dictget_strvalue(yaml_dict, 'forcing_file')
    forcing_data_dates_dt_max = dictget_dpvalue(yaml_dict, 'forcing_dates_dt_max')
    forcing_use_dynamic_temperature = dictget_intvalue(yaml_dict, 'forcing_use_dynamic_temperature')
    forcing_use_dynamic_relative_humidity = dictget_intvalue(yaml_dict, 'forcing_use_dynamic_relative_humidity')
    forcing_use_dynamic_wind_speed = dictget_intvalue(yaml_dict, 'forcing_use_dynamic_wind_speed')
    forcing_use_dynamic_sunshine = dictget_intvalue(yaml_dict, 'forcing_use_dynamic_sunshine')
    forcing_use_dynamic_pressure = dictget_intvalue(yaml_dict, 'forcing_use_dynamic_pressure')
    !open dynamic file
    call check_ncrequest( nf90_open(trim(adjustl(forcing_data_file)), NF90_NOwrite, FILPLU) )
    nc = nc_dimension(FILPLU, "n_meshes")
    !get dates in dynamic forcing file
    call check_ncrequest( nf90_inq_varid(FILPLU, "dates", forcing_dates_id) )
    call check_ncrequest( nf90_inq_dimid(FILPLU, 'nt', trash_int) )
    call check_ncrequest( nf90_inquire_dimension(FILPLU, trash_int, len = n_forcing_file_dates) )
    allocate(forcing_file_julianday_dates(n_forcing_file_dates))
    call check_ncrequest( nf90_get_var(FILPLU, forcing_dates_id, forcing_file_julianday_dates) )
    !get variables id for dynamic forcing file
    call check_ncrequest( nf90_inq_varid(FILPLU, "rain", forcing_rain_id) )
    if (forcing_use_dynamic_temperature == 1) then
        call check_ncrequest( nf90_inq_varid(FILPLU, "temperature", forcing_varid_dynamic_temperature) )
    end if
    if (forcing_use_dynamic_relative_humidity == 1) then
        call check_ncrequest( nf90_inq_varid(FILPLU, "relative_humidity", forcing_varid_dynamic_relative_humidity) )
    end if
    if (forcing_use_dynamic_wind_speed == 1) then
        call check_ncrequest( nf90_inq_varid(FILPLU, "wind_speed", forcing_varid_dynamic_wind_speed) )
    end if
    if (forcing_use_dynamic_sunshine == 1) then
        call check_ncrequest( nf90_inq_varid(FILPLU, "sunshine", forcing_varid_dynamic_sunshine) )
    end if
    if (forcing_use_dynamic_pressure == 1) then
        call check_ncrequest( nf90_inq_varid(FILPLU, "pressure", forcing_varid_dynamic_pressure) )
    end if

    
    !Reads Boundary condition (Water surface slope)
!~     WSslope = dictget_dpvalue(yaml_dict, 'water_surface_slope')*0.001 !converting from m/km to m/m
    
    !files for hydrological state I/O => for hotstart
    hydrological_state_read_file = dictget_strvalue(yaml_dict, 'hydrological_state_read_file')
    hydrological_state_write_prefix_each_step = dictget_strvalue(yaml_dict, 'hydrological_state_write_prefix_each_step')
    hydrological_state_write_file = dictget_strvalue(yaml_dict, 'hydrological_state_write_file')    
    
    flagaclimed = 1
    cru_index = 1

    !observations
    NOBS = dictget_listlength(yaml_dict, 'obs_cell_indexes')
    if (NOBS .gt. 0) then
        ARQOBS = dictget_strvalue(yaml_dict, 'obs_file')
        allocate(IQOBS(NOBS))
        IQOBS = dictget_intarray(yaml_dict, 'obs_cell_indexes', NOBS)
        call LEQOBS	!Subroutine for reading the file observed streamflow data (File with extension .QOB)
    else
        ARQOBS = 'none'
    end if
    
    !substitute stations
    NUMSUBST = dictget_listlength(yaml_dict, 'subst_cell_indexes')
    if (NUMSUBST .gt. 0) then
        ARQSUBST = dictget_strvalue(yaml_dict, 'subst_file')
        allocate(ISUBST(NUMSUBST),ISUBSTAUX(NUMSUBST))
        ISUBST = dictget_intarray(yaml_dict, 'subst_cell_indexes', NUMSUBST)
        call LESUBST !Then calls the subroutine for the reading the file with to be substituded streamflow data (File: QSUBST.QSB) 
    else
        ARQSUBST = 'none'
    end if

    !Load static data
    call load_static_data(dictget_strvalue(yaml_dict, 'static_data_file'))
    
    
end subroutine initialize_model
