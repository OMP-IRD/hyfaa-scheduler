    !*********************************************************************************
    !
    !  SUBROUTINE SIMULA is the main routine for hydrological simulation
	!                    (wo calibration)
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine has the main loop of the MGB-IPH, the time loop, from iT=1,nT
    !     where iT is the time interval and nT is the number of time steps.
	!
	!	 For each time interval date info is set, then rainfall and climate data are
	!	   loaded through subroutine read_forcing. Afterwards catchment flow
	!	   generation is done using subroutine CELULA. Finally, river routing is 
	!	   achieved by calling (i) REDE, for Muskingum-Cunge Method or
	!	   (ii) flood_inertial, for a 1D Hydrodynamic Model (w/ inertia and pressure)
	!	
	!	At the end discharge time series are stored in :
	!		QRB: calculated discharge time series in all subbasins
	!		QR:  calculated discharge time series in subbasin outlet w observed data
	!	Those are recorded in files at the end when returns to SIMULA
	!
	!
	!	 * iT value should not be changed inside subroutines!
	!
	!    * AUX_MOD module from full hydrodynamic is deactivated (commented)!
	!    * Hidrodinamico2 subroutine from full hydrodynamic is deactivated (commented)!
    !
    !  Usage:
    !
    !    call SIMULA
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
	!	 * module     IFPORT  			from (visual studio 2008)
    !    * module     vars_main   		in    vars_main.f90
	!    * module     vars_calib   		in    vars_calib.f90
    !    * module     vars_inerc  		in    VARSINERC.f90  !Quais? 
 	!	 * subroutine CONDINIC	  		in	  CONDINIC.f90
	!	 * subroutine MODELO	  		in	  MODELO.f90
	!	 * subroutine FOBJ	  			in	  FOBJ.f90	
	!
    !
    !	 opens
    !
    !    * Opens QPROP.PRP  	   output  ascii  file with calculated discharge by origin (i.e. surface, subsurface, gw)
	!	 * Opens VAZAO.QCL 		   output  ascii  file with calculated discharge in subbasins defined in setup file
	!	 * Opens Qmes90.TXT 	   output  ascii  file with calculated Q90 in all catchments
	!	 * Opens Qmes95.TXT 	   output  ascii  file with calculated Q95 in all catchments
	!	 * Opens RESUMO_SIAQUA.TXT output  ascii  file for SIAQUA	
    !
    !    reads
    !
    !    * Does not read files
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
    ! 2015.06.21 - 21 June 2015 (By: Fernando Mainardi Fan)   
    !
    !  Authors:
    !
    !    Original fortran version by Walter Collischonn
    !    Present fortran version by:
    !    * Walter Collischonn
    !    * Rodrigo Cauduro Dias de Paiva
    !    * Diogo da Costa Buarque
    !    * Paulo Pontes Rógenes
    !    * Mino  Viana Sorribas
    !    * Fernando Mainardi Fan
    !    * Juan Martin Bravo 
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
    !   * Variables declarations and routines calls are all commented below.
    !
    !---------------------------------------------------------------------------------

	SUBROUTINE SIMULA()	
	use global_variables
	use hydrological_state_restart_io !hydrological state restart inputs/outputs subroutines
	
	implicit none
 
    integer:: time_meas1, time_meas2, clock_rate, clock_max
   
    
    call system_clock(time_meas1, clock_rate, clock_max)
    
    ! Main Routines
    print*, 'Is there somebody here?'
    print*, 'Calling initial conditions...'
    call CONDINIC  ! Initial Conditions
    if (trim(adjustl(hydrological_state_read_file)).ne.'none') then
        print*, 'Loading hydrological state from file '//trim(adjustl(hydrological_state_read_file))//' ...'
        call read_hydrological_state(trim(adjustl(hydrological_state_read_file)))
    end if
    print*, 'HS read...'
    call system_clock(time_meas2, clock_rate, clock_max)
    print*, '    TIMING: simula init = ', real(time_meas2-time_meas1)/real(clock_rate)
    call system_clock(time_meas1, clock_rate, clock_max)



    print*, 'Running main model routine...'
    call MODELO    ! Time Loop
    call system_clock(time_meas2, clock_rate, clock_max)
    print*, '    TIMING: MODELO = ', real(time_meas2-time_meas1)/real(clock_rate)
    call system_clock(time_meas1, clock_rate, clock_max)
    

    !dump hydrological state if requested
    if (trim(adjustl(hydrological_state_write_file)) .ne. 'none') then
        call write_hydrological_state(trim(adjustl(hydrological_state_write_file)))
        call system_clock(time_meas2, clock_rate, clock_max)
        print*, '    TIMING: write hydro state = ', real(time_meas2-time_meas1)/real(clock_rate)
        call system_clock(time_meas1, clock_rate, clock_max)
    end if



	end
