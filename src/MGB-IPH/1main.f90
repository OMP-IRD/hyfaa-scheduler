    PROGRAM MGB_Inercial
!*********************************************************************************
!***************MODELO HIDROLÓGICO DE GRANDES BACIAS (MGB-IPH)********************
!*********************************************************************************
!*		Large Scale Hydrology Research Group          	      *
!*	      Instituto de Pesquisas Hidráulicas (IPH)      	  *
!*	 Universidade Federal do Rio Grande do Sul (UFRGS)        *
!*                  	   June 2015                          *
!*                  	  Version 2.0                         *
!*                                                            *
!*********************************************************************************
!*********************************************************************************
!
!  PROGRAM MGB_Inercial (1main.f90) is the main routine from MGB-IPH.
!  Modelo de Grandes Bacias - Instituto de Pesquisas Hidraulicas
!
!---------------------------------------------------------------------------------
!  Opening:
!
!    This model and routines were developed by Walter Collischonn
!    during his research thesis in early 2001, in Brazil.
!
!    From there on MGB-IPH was modified for research and project purposes.
!    The current package is the official version of MGB-IPH.
!
!   As references we list below some publications where the model was used:
!
!   MELLER, A. ; COLLISCHONN, W. ; FAN, F. M. ; BUARQUE, D. C. ; PAIVA, R. C. D. ; DIAS, P. ; MOREIRA, D. . Previsão de Cheias por Conjunto em Curto Prazo. Revista Brasileira de Recursos Hídricos, v. 19, p. 33-49, 2014.
!   FAN, F. M.; COLLISCHONN, W. ; MELLER, A. ; BOTELHO, L. C. M. Ensemble streamflow forecasting experiments in a tropical basin: The São Francisco river case study. Journal of Hydrology (Amsterdam), v. SI, p. 10.1016/j.jhydr, 2014.
!	FAN, F. M.; COLLISCHONN, W. .Integração do Modelo MGB-IPH com Sistema de Informação Geográfica. Revista Brasileira de Recursos Hídricos, v. 19, p. 243-254, 2014.
!	Paiva, R. C. D. ; Collischonn, W. ; Buarque, D. C. . Validation of a full hydrodynamic model for large scale hydrologic modelling in the Amazon. Hydrological Processes, 27, p. 333–346. DOI: 10.1002/hyp.8425, 2013.
!	Paiva, R. C. D. ; Paiva, R. C. D. ; Collischonn, W. ; Bonnet, M.-P. ; De Gonçalves, L. G. G. ; Calmant, S. ; Getirana, A. ; Santos Da Silva, J. . Assimilating in situ and radar altimetry data into a large-scale hydrologic-hydrodynamic model for streamflow forecast in the Amazon. Hydrology and Earth System Sciences Discussions (Online), v. 10, p. 2879-2925, 2013.
!	Paiva, R. C. D.; Buarque, D. C. ; Collischonn, W. ; Bonnet, M.-P. ; Frappart, F.; Calmant, S.; Bulhões Mendes, C. A.. Large-scale hydrologic and hydrodynamic modeling of the Amazon River basin. Water Resources Research, v. 49, p. 1226-1243, 2013.
!	Collischonn, W. ; Meller, A. ; Fan, F. ; Moreira, D.S. ; Silva Dias, P.L. ; Buarque, D.; Bravo, J. M. . Short-term Ensemble Flood Forecasting Experiments in Brazil. Geophysical Research Abstracts, v. 15, p. 11910, 2013.
!	Bravo, J. M. ; Allasia, D. ; Paz, A. R. ; Collischonn, W. ; Tucci, C. E. M. . Coupled Hydrologic-Hydraulic Modeling of the Upper Paraguay River Basin. Journal of Hydrologic Engineering, v. 17, p. 635, 2012.
!	Paiva, R. C. D.; Collischonn, W.; Bonnet, M. P.; De Gonçalves, L. G. G. . On the sources of hydrological prediction uncertainty in the Amazon. Hydrology and Earth System Sciences, v. 16, p. 3127-3137, 2012.
!	Paiva, R. C. D.; Collischonn, W.; Buarque, D. C.. Validation of a full hydrodynamic model for large-scale hydrologic modelling in the Amazon. Hydrological Processes (Print), v. -, p. n/a-n/a, 2012.
!	Pontes, P. R. M. ; Collischonn, W. . Conservação de Volume em Modelos Simplificados de Propagação de Vazão. Revista Brasileira de Recursos Hídricos, v. 17, p. 83-96, 2012.
!	Sorribas, M. V.; Collischonn, W.; Marques, D. M.; Fragoso Jr., C. R.; Castro, N. M. R.; Souza, R. S. . Modelagem Distribuída do Carbono em Bacias Hidrográficas. Revista Brasileira de Recursos Hídricos, v. 17, p. 225-240, 2012.
!	Meller, A. ; Bravo, J. M. ; Collischonn, W. . Assimilação de Dados de Vazão na Previsão de Cheias em Tempo-Real com o Modelo Hidrológico MGB-IPH. Revista Brasileira de Recursos Hídricos, v. 17, p. 209-224, 2012.
!	Buarque, D. C. ; Collischonn, W. ; Paiva, R. C. D. . Coupling a basin erosion and river sediment transport model into a large scale hydrological model: an application in the Amazon basin. Geophysical Research Abstracts, v. 14, p. 11935, 2012.
!	Nóbrega, M. T. ; Collischonn, W. ; Tucci, C. E. M. ; Paz, A. R. . Uncertainty in climate change impacts on water resources in the Rio Grande Basin, Brazil. Hydrology and Earth System Sciences, v. 15, p. 585-595, 2011.
!	Paz, Adriano Rolim da ; Collischonn, Walter ; Tucci, Carlos E. M. ; Padovani, Carlos R. . Large-scale modelling of channel flow and floodplain inundation dynamics and its application to the Pantanal (Brazil). Hydrological Processes (Print), v. 25, p. 1498-1516, 2011.
!	Paiva, Rodrigo C.D. ; Collischonn, Walter ; Tucci, Carlos E.M. . Large scale hydrologic and hydrodynamic modeling using limited data and a GIS based approach. Journal of Hydrology (Amsterdam), v. 406, p. 170-181, 2011.
!	Collischonn, B. ; Paiva, R. C. D. ; Meirelles, F. S. C. ; Collischonn, W. ; Fan, F. M. ; Camano, E. . Modelagem Hidrológica de Uma Bacia com Uso Intensivo de Água: Caso do Rio Quaraí-RS. Revista Brasileira de Recursos Hídricos, v. 16, p. 119-133, 2011.
!	Paiva, R. C. D. ; Buarque, D. C. ; Collischonn, W. ; Sorribas, M. ; Allasia, D. G. ; Mendes, C. A. B. ; Tucci, C. E. M. ; Bonnet, M. P. . Hydrologic and Hydrodynamic Modelling of the Amazon Basin using TRMM Rainfall Estimates. Geophysical Research Abstracts, v. 13, p. 12666, 2011.
!	Paiva, R. ; Buarque, Diogo ; Collischonn, W. ; Sorribas, M. ; Allasia, D. G. P. ; Mendes, C. A. B. ; Tucci, C. E. M. ; Tucci, C. E. M. ; Bonnet, M. . Using TRMM rainfall estimates in hydrological and hydrodynamic modelling of the Amazon Basin. IAHS-AISH Publication, v. 343, p. 72-77, 2011.
!	Paz, A. R. ; Bravo, J. M. ; Allasia, D. ; Collischonn, W. ; Tucci, C. E. M. . Large-Scale Hydrodynamic Modeling of a Complex River Network and Floodplains. Journal of Hydrologic Engineering, v. 15, p. 152-165, 2010.
!	Paz, A. R. ; Collischonn, W. ; Tucci, C. E. M. . Simulação hidrológica de rios com grandes planícies de inundação. Revista Brasileira de Recursos Hídricos, v. 15, p. 31-43, 2010.
!	Paz, A. R. ; Collischonn, W. . Derivação de rede de drenagem a partir de dados do SRTM. Revista Geográfica Acadêmica, v. 2, p. 84-95, 2008.
!	Larentis, D. G. ; Collischonn, W. ; Tucci, C. E. M. . Simulação da qualidade de água em grandes bacias: Rio Taquari-Antas, RS. Revista Brasileira de Recursos Hídricos, v. 13, p. 5-22, 2008.
!	Tucci, Carlos E. M. ; Collischonn, W. ; Clarke, R. T. ; Clarke, Robin T. ; Paz, Adriano R. ; Allasia, Daniel . Short- and long-term flow forecasting in the Rio Grande watershed (Brazil). Atmospheric Science Letters, v. 9, p. 53-56, 2008.
!	Collischonn, W. ; Tucci, C. E. M. ; Clarke, R. T. ; Chou, S. C. ; Guilhon, L. G. ; Cataldi, M. ; Allasia, D. G. . Medium-range reservoir inflow predictions based on quantitative precipitation forecasts. Journal of Hydrology, v. 344, p. 112-122, 2007.
!	Collischonn, W. ; Allasia, D. G. ; Silva, B. C. ; Tucci, C. E. M. . The MGB-IPH model for large-scale rainfall-runoff modelling. Hydrological Sciences Journal, v. 52, p. 878-895, 2007.
!	Collischonn, W. ; Tucci, C. E. M. ; Clarke, R. T. ; Delgado, M. C. ; Silva, B. C. ; Collischonn, B. ; Allasia, D. G. ; Paz, A. R. . Modelo hidrológico distribuído para previsão de vazões incrementais na bacia do rio Paranaíba entre Itumbiara e São Simão. Revista Brasileira de Recursos Hídricos, v. 12, p. 43-56, 2007.
!	Silva, B. C. ; Collischonn, W. ; Tucci, C. E. M. ; Clarke, R. T. ; Delgado, M. C. . Previsão hidroclimática de vazão de curto prazo na bacia do rio São Francisco. Revista Brasileira de Recursos Hídricos, v. 12, p. 31-42, 2007.
!	Collischonn, W. ; Silva, B. C. ; Tucci, C. E. M. ; Allasia, D. G. . Large basin simulation experience in South America. IAHS Pubblication n. 303, v. 303, p. 360-370, 2006.
!	Silva, B. C. ; Tucci, C. E. M. ; Collischonn, W. . Previsão de vazão com modelos hidroclimáticos. Revista Brasileira de Recursos Hídricos, v. 11, p. 15-30, 2006.
!	Andreolli, I. ; Collischonn, W. ; Tucci, C. E. M. ; Haas, R. ; Regina, J. V. M. . Previsão de vazão afluente a um reservatório utilizando previsão quantitativa de chuva. Revista Brasileira de Recursos Hídricos, v. 11, p. 55-70, 2006.
!	Ribeiro Neto, A. ; Collischonn, W. ; Silva, R. C. V. ; Tucci, C. E. . Hydrological modelling in Amazonia use of the MGB-IPH model and alternative databases. IAHS-AISH Publication, v. 303, p. 246-254, 2006.
!	Collischonn, W. ; Tucci, C. E. M. ; Haas, R. ; Andreolli, I. . Forecasting river Uruguay flow using rainfall forecasts from a regional weather-prediction model. Journal of Hydrology, v. 305, p. 87-98, 2005.
!	Collischonn, W. ; Tucci, C. E. M. . Previsão Sazonal de vazão na bacia do rio Uruguai 1: Ajuste e verificação do modelo hidrológico distribuído. Revista Brasileira de Recursos Hídricos, v. 10, n.4, p. 43-59, 2005.
!	Collischonn, W. ; Tucci, C. E. M. ; Clarke, R. T. ; Dias, P. L. S. ; Sampaio, G. O. . Previsão sazonal de vazão na bacia do rio Uruguai 2: Previsão Climática-Hidrológica. Revista Brasileira de Recursos Hídricos, v. 10, n.4, p. 60-72, 2005.
!	Bravo, J. M. ; Collischonn, W. ; Pilar, J. V. ; Silva, B. C. ; Tucci, C. E. M. . Operação de um reservatório com múltiplos usos com base na previsão de curto prazo. Revista Brasileira de Energia, v. 11, p. 85-110, 2005.
!	Collischonn, W. ; Tucci, C. E. M. ; Clarke, R. T. . Variabilidade temporal no regime hidrológico da bacia do rio Paraguai. Revista Brasileira de Recursos Hídricos, Porto Alegre RS, v. 8, n.1, p. 201-211, 2003.
!	Tucci, C. E. M. ; Dias, P. L. S. ; Clarke, R. T. ; Sampaio, G. O. ; Collischonn, W. . Long-term flow forecasts based on climate and hydrologic modeling: Uruguay river basin. Water Resources Research, New York, v. 39, n.7, p. 1-2, 2003.
!	Collischonn, W. ; Tucci, C. E. M. . Ajuste multiobjetivo dos parâmetros de um modelo hidrológico. Revista Brasileira de Recursos Hídricos, Porto Alegre, v. 8, n.3, p. 27-39, 2003.
!	Collischonn, W. ; Tucci, C. E. M. . Simulação hidrológica de grandes bacias. Revista Brasileira de Recursos Hídricos, v. 6, n.2, 2001.
!   COLLISCHONN, W. Modelagem de Grandes Bacias - ph.d. Thesis. 2001 
!
!    Older source codes have been distributed by some people in
!    past research, thus any other Large Scale Model with similar source
!    codes are non-official and doesn't respect copyrights and patents.
!
!    This distribution is a version of the model in English.
!
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!  1main.f90
!---------------------------------------------------------------------------------
!  Discussion:
! 
!    This is the main routine of the model.
!    This routines starts the model, open some input files and calls
!    model simulation modes (i.e. simulation, calibration) etc...
!
!  Usage:
!
!    PROGRAM MGB_Inercial
!
!    uses modules, functions, and subroutines
!
!	*Module	PORTLIB     !library to calculate the processing time
!	*Module	IFPORT      !library to calculate the processing time
!	*Module	global_variables   !module containing the main variables used in the model
!	*Function	TEMPO(0,ISEED)  ! INICIA CONTAGEM DO TEMPO
!	*Subroutine	initialize_model(yaml_file)  ! initialize main parameters and launch load_static_data which loads all model information and climatology
!	*Subroutine	ARQCLISUB	  !Subroutine for reading the files with daily/hourly meteorological information (Files with extension .CLI)
!	*Subroutine	LEQOBS	!Subroutine for reading the file observed streamflow data (File with extension .QOB)
!	*Subroutine	LESUBST !Then calls the subroutine for the reading the file with to be substituded streamflow data (File: QSUBST.QSB) 
!	*Subroutine	PARCEL !Calculates parameters related to MicroBains and rivers
!	*Subroutine	PARCUNGE !Calcualtes Muskingum-Cunge parameters
!	*Subroutine	flood_TOPO !Creates a matrix with topology information that is used in the Local Inertial routing method
!	*Subroutine	flood_TAB !Creates a table with the volume of water in the floodplain from the table of water depth vs. area obtained by DEM preprocessing (File with extension .FLP)
!	*Subroutine	SIMULA !Subroutine that calls a standard model simulation run
!	*Subroutine	LeCalib !Subroutine that reads automatic calibration input files
!	*Subroutine	ALLOCA_CALIB(0) !Subroutine that allocates automatic calibration input files
!	*Subroutine	CALIBRA !Subroutine that alls a model automatic calibration run
!	*Subroutine	ALLOCA_CALIB(1) !Subroutine that cleans automatic calibration input files
!	*Subroutine	PREVISAO !Subroutine that calls a forecasting run
!	*Function	TEMPO(1,0)  !Finish time count
!	*Subroutine	ALLOCA_VARS(1) !Cleans main variables
!
!	 opens
!
!    * CHUVAbin.pbi, the precipitation binary file
!    * NOSOLO.tx, information about soil water at one minibasin
!    * AJUSTE.fob, model performance coeficcients
!
!    reads
!
!    * no files are read in this routine
!
!    creates
!
!    * NOSOLO.tx, information about soil water at one minibasin
!    * AJUSTE.fob, model performance coeficcients
!
!---------------------------------------------------------------------------------
!  Licensing:
!
!    This code is distributed under the...
!
!  Version/Modified: 
!
!    21 June 2015
!    By: Fernando Mainardi Fan
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
!  Main References:
!  COLLISCHONN, W. ; ALLASIA, D. G. ; SILVA, B. C. ; TUCCI, C. E. M. The MGB-IPH model for large-scale rainfall-runoff modelling. Hydrological Sciences Journal, v. 52, p. 878-895, 2007.
!  COLLISCHONN, W., TUCCI, C. E. M. Simulação hidrológica de grandes bacias. Revista Brasileira de Recursos Hídricos, v. 6, n. 2, 2001.
!  COLLISCHONN, W. Modelagem de Grandes Bacias - ph.d. Thesis. 2001 
!
!---------------------------------------------------------------------------------
!  Variables and Parameters:
!  *Variables declarations and routines calls are all commented below.
!---------------------------------------------------------------------------------
! End of header
!---------------------------------------------------------------------------------

!  Calls used libraries and modules 

    use global_variables
    use basic_functions
    use julianday_mod
    use netcdf
    use netcdf_addons
    
!  Variables declaration 

    implicit none !States that non declared variables cannot be used
    real TIME_SECONDS !Variable to calculate the processing time
    integer::nargs_command_line
    logical::use_init_nc
    character(len=400)::init_nc_path, yaml_input_file
    integer:: time_meas1, time_meas2, clock_rate, clock_max
!  Statements used to calculate the processing time


    call TEMPO(0,ISEED)  ! INICIA CONTAGEM DO TEMPO
    TIME_SECONDS = SECNDS(0.0)	
    
    !   main input file
    nargs_command_line = command_argument_count()
    if (nargs_command_line.ne.1) then
        print*, 'Expecting only one command line argument (path to yaml input file)'
        call abort()
    end if
    call get_command_argument(1, yaml_input_file, trash_int, status_int)
    if (status_int .ne. 0) then
        print*, 'Could not retrieve INFO_MGB_FILE from command line, command must be: PROGRAM INFO_MGB_FILE'
        call abort()
    end if
    call system_clock(time_meas1, clock_rate, clock_max)
    call initialize_model(yaml_input_file)  !Subroutine that reads the model main input files with the running main information 
    call system_clock(time_meas2, clock_rate, clock_max)
    print*, '    TIMING: Total static data loading time = ', real(time_meas2-time_meas1)/real(clock_rate)


!  End of openning MGB-IPH main Input and Output files that are read and written over the simulation

!  Ask for verification
!  Informs the model general setup


    print*,'Input data reading process ended'
    print*,''
    write(*,'(A30,I8)')'Number of catchments:',NC
    write(*,'(A30,I8)')'Number of subbasins:',NB
    write(*,'(A30,I8)')'Number of HRUs:',NU
    write(*,'(A30,4I5)')'Start date (dd mm yyyy hh):',IDIA,IMES,IANO,HORAINI 
    write(*,'(A30,F16.2)')'Time step (s):',DTP
    write(*,'(A30,I8)')'Number of time steps:',NT
    write(*,'(A30,I8)')'Number of observation gauges:',NOBS
    write(*,'(A30,I8)')'Number of substitution gauges:',NUMSUBST
    write(*,'(A30,I8)')'Number of climate gauges:',NC
    print*,'ICALIB (simulation=0; autocalibration=1; forecasting=2)'
    write(*,'(A30,I8)')'ICALIB :',ICALIB
    print*,'PLEASE VERIFY INPUT INFORMATION AND PRESS ENTER TO START RUN...'
    print*,''
    !PAuse


!  Select the mode of MGB-IPH run: Simulation (ICALIB=0), Automatic Calibration (ICALIB=1), Forecasting (ICALIB=2)

	CALIBRA_case: select case (ICALIB) !Verify the mode
	case (0) !Standard model simulation run
		print*,'SIMULATION MODE' !Warns about the chosen mode
		call SIMULA() !Subroutine that calls a standard model simulation run
        call system_clock(time_meas2, clock_rate, clock_max)
        print*, '    TIMING: SIMULA = ', real(time_meas2-time_meas1)/real(clock_rate)
        call system_clock(time_meas1, clock_rate, clock_max)
	case (1) !Automatic Calibration run
		print*,'AUTOMATIC CALIBRATION MODE' !Warns about the chosen mode
		call lecalib !Subroutine that reads automatic calibration input files
		call calibra !Subroutine that alls a model automatic calibration run
	case default
		 print*, ' ERROR: UNKNOWM ICALIB!!!' !Warns about unknown ICALIB's
		call abort()
	end select CALIBRA_case

!  End of MGB-IPH

!  After running the model, finish the time count, close all files, and clean all variables 

    call TEMPO(1,0)  !Finish time count
    deallocate(forcing_file_julianday_dates)
    call check_ncrequest( nf90_close(FILPLU) )
    
! End of closings and cleanings

! Inform the user that the model run is over and informations about total processing time

    print*,'MGB-IPH endED SUCCESSFULLY'
    print *, char(7)
    TIME_SECONDS = SECNDS(TIME_SECONDS)
    print*,'TOTAL TIME:',TIME_SECONDS/60.,'MINUTES'
    print*,'PRESS ENTER TO close...'
    !PAuse
    !READ(*,*)
    
end         
