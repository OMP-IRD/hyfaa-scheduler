module global_variables

	!Variables declaration

	implicit none
	save


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !main variables

	!Minibasins and muskingum-cunge routing
	integer,PARAMETER:: NUMUSK=150 !Maximum number of stretchs for MUSKINGUM CUNGE routing
	integer,PARAMETER:: NFP=3 !Number of objective functions for each gauge
	integer,PARAMETER:: NTMUP=20 !Number of points in the table of non-linear muskingun Cunge
	
	!Numbers
	integer NC,NU !minibasins, uses
	integer NT !time intervals
	integer NB !number of sub-basins
    integer::n_regions
	integer NOBS !Number of observed streamflow series, number of minibasins with observed data
	integer NUMSUBST !Number of minibasins where streamflow calculated will be replaced by the information read from a file

    
	integer::trash_int, trash_int2, trash_int3, trash_int4, trash_int5, status_int
	real::trash_real, trash_real2
	double precision::trash_dp, trash_dp2
    logical::is_congo
    character(len=100):: hydrodynamic_module
	character(len=1024):: trash_str, temp_str, temp_str1, temp_str2
	CHARACTER(len=1024):: output_directory, forcing_data_file, param_calib_file
	CHARACTER(len=1024):: hydrological_state_read_file, hydrological_state_write_file, hydrological_state_write_prefix_each_step
	CHARACTER(len=1024):: old_water_balance_state_write_file, old_water_inertial_state_write_file
	CHARACTER(len=1024):: initial_population_file, geomorphological_relation_file, par_calib_file
	CHARACTER(len=10) CABE(1000)
	integer:: FILPLU, forcing_dates_id, forcing_rain_id,  n_forcing_file_dates
    integer::forcing_use_dynamic_temperature, forcing_use_dynamic_relative_humidity, forcing_use_dynamic_wind_speed, forcing_use_dynamic_sunshine, forcing_use_dynamic_pressure
    integer::forcing_varid_dynamic_temperature, forcing_varid_dynamic_relative_humidity, forcing_varid_dynamic_wind_speed, forcing_varid_dynamic_sunshine, forcing_varid_dynamic_pressure
	double precision, dimension(:), allocatable::forcing_file_julianday_dates
	double precision:: simulation_start_date
	double precision::forcing_data_dates_dt_max
	integer,PARAMETER:: filtemp=775
	integer,PARAMETER:: filtemp2=776
	integer, dimension(:), allocatable::IQOBS,ISUBST,ISUBSTAUX

	!Counters
	integer IT !Time Interval

	!Variables with numbers of files
	!Input Files
	integer,PARAMETER:: FILFIX=131 !File of fixed parameters
	integer,PARAMETER:: FILVAR=132 !File of calibrated parameters
	integer,PARAMETER:: FILUSO=133 !File of parameters associated with the HRU's
	integer,PARAMETER:: FILHIG=134 !File of minibasins informations
	integer,PARAMETER:: FILMED=136 !File of monthly mean climatological data
	integer,PARAMETER:: FILOBS=137 !File of observed streamflows
	integer,PARAMETER:: FILLIM=138 !File of parameters limits
	integer,PARAMETER:: FILSUBS=140 !Hydrographs of file read to substitute the calculated ones
	integer,PARAMETER:: FILOUTROS=141 !Hydrograph of file observed on the axis of san francisco
	integer,PARAMETER:: FILCALIB=142 !File with the description of parameters for calibration
    integer,PARAMETER:: FILRG=143 !File of geomorphological relationships

	!Output Files
	integer,PARAMETER:: FILPRP=145 !Hydrograph file with proportions
	integer,PARAMETER:: FILAJU=146 !Model Perfomance file
	integer,PARAMETER:: FILBAC=147 !Output file of averaged data to the basins
	integer,PARAMETER:: FILSOL=148 !Soil file (one HRU)
	integer,PARAMETER:: FILSOL2=149 !Soil file (one HRU)
	integer,PARAMETER:: FILEVO=150 !EVOLUTION parameters file evolution
	integer,PARAMETER:: FILBOM=151 !File containing the best parameters yet
	integer,PARAMETER:: FILTUD=152 !File containing all mini-basins calculated hydrographs

!	!Calibration
	integer ICALIB
	real,ALLOCATABLE:: R2(:),ERRV(:),R2L(:) !coef. nash, error in volumes, coef nash logarithms
	real VFO(NFP) !Vector with values of the objective functions

	!Variables with dimensions used
	real,ALLOCATABLE:: QR(:,:),QOBS(:,:) !Streamflow in the outlets of all sub-basins
	real,ALLOCATABLE:: QLIDO(:,:) !Streamflow at some point to replace the calculated values
	real,ALLOCATABLE:: QMANTEIGA(:),QSAOROM(:),QSAOFRA(:),QPMCRUZ(:),QMANGA(:),QCARI(:)
	real,ALLOCATABLE:: QBJLAPA(:),QGAMEL(:),QMORP(:)
	integer,ALLOCATABLE:: KCB(:)!Nnumber of cells in each sub-basin
	
	!Mean vales at each Sub-basin
	real,ALLOCATABLE:: DSUM(:,:),DINM(:,:),DBAM(:,:)
	real,ALLOCATABLE:: SIM(:,:),EM(:,:),WBM(:,:),PM2(:,:)
	
	!Streamflow in accordance with the origin for each one of the time steps
	real,ALLOCATABLE:: QB(:),QBI(:),QBIS(:)
	
	integer IDIA,IMES,IANO,IDINI,IHORA !Date and day of the Julian 
	real,ALLOCATABLE:: QESP(:)!Especific baseflow (M3/S/KM2)
	
	!Albedo, Leaf Area Index, average height, and surface resistance
    integer,allocatable,dimension(:)::cell_region_id
	real,ALLOCATABLE:: ALB(:,:,:),RIAF(:,:,:),Z(:,:,:),RS(:,:,:)
	
	!Number of stretchs and time intervals used in the calculation  of muskingum-cunge
	integer,ALLOCATABLE:: NSUBT(:)
	real,ALLOCATABLE:: DT(:)
	real,ALLOCATABLE:: CEL(:),TIND(:)!River celerity and concentration time of the minibasin
	
	
	!File names of monthly climatological data
	CHARACTER (20) ACLIMED
	
	!File name of the streamflows observed data
	CHARACTER (20) ARQOBS
	
	!Name of the file that have the streamflows to override the flows calculated in some minibasins
	CHARACTER (20) ARQSUBST
	
	!Parameters related to the soil
	real,ALLOCATABLE:: WM(:,:),B(:,:),KINS(:,:),KBAS(:,:)
	real,ALLOCATABLE:: WMOLD(:,:),BOLD(:,:),KBOLD(:,:),KIOLD(:,:)
	real,ALLOCATABLE:: PLAM(:,:),CAP(:,:),WC(:,:)
	real,ALLOCATABLE:: CI(:),CB(:),CS(:)!Propagation parameters of the minibasin
	real,ALLOCATABLE:: CIOLD(:),CSOLD(:)!Copy of the parameters of propagation in the minibasin
	
	!Parameters related to minibasins
	integer,ALLOCATABLE:: IEXUT(:)!Indicates the outlet minibasins
	real,ALLOCATABLE:: BRIO(:)!River width
	real BRX !River width (auxiliary variable)
	real,ALLOCATABLE:: PM(:) !Average rainfall
	integer,ALLOCATABLE:: IBAC(:),CELJUS(:)!Watershed, downstream minibasin
	real,ALLOCATABLE::  HCEL(:),LCEL(:) ! Declividade e comprimento do rio mais longo
	real,ALLOCATABLE:: X(:),Y(:)!Coordinates of the minibasin center
	!Area of the minibasin, drainage area in number of minibasins, river length and slope of the river
	real,ALLOCATABLE:: ACEL(:),ACUR(:),SRIO(:),DECL(:)
	!Proportion of HRUs in the minibasins
	real,ALLOCATABLE:: PUSO(:,:)
	integer,ALLOCATABLE:: ICBOM(:)!Number of climatological stations in the cell
	real,ALLOCATABLE:: QREF(:) !Reference streamflow
	real QRX !Reference streamflow (auxiliary variable)
	! temp., humidity, wind, insolation, pressure, daily at each station
	real,ALLOCATABLE:: TD(:,:),UD(:,:),VD(:,:),SOLD(:,:),PAD(:,:) 
	real,ALLOCATABLE:: PMB(:)!By basin average rainfall
	integer,ALLOCATABLE:: KPM(:) !Auxiliar variable
	! temp., humidity, wind, insolation, pressure, monthly values
    integer:: n_climatology
    integer, allocatable:: specific_climatology_years(:)
	real,ALLOCATABLE:: TAMM(:,:,:),URMM(:,:,:),VVMM(:,:,:),PAMM(:,:,:),SOLMM(:,:,:)
	real,ALLOCATABLE:: XYC(:,:)!x and y coordinates of climatological stations

	!Counters
	integer IC
	

	!More variables
	real,ALLOCATABLE:: SI(:,:) !Amount of water intercepted in the surface
	real XLADO,DH,ARX !Length of the minibasin side, difference in altitude, area (auxiliary variable)
	real,ALLOCATABLE:: QRB(:,:)	!Stores hydrographs at the outlets
	real,ALLOCATABLE:: QM1(:),QJ1(:),QM2(:),QJ2(:) !Streamflows upstream and downstream in the cell i at times 1 and 2
	real,ALLOCATABLE:: QCEL1(:),QCEL2(:) !Streamflows originated in the cell at time t and t + 1 in cell i
	real,ALLOCATABLE:: PMB2(:),PMI2(:),PMS2(:),PJB2(:),PJI2(:),PJS2(:) !Proportions of origin of the flows in the river
	real,ALLOCATABLE:: VRB(:),VRI(:),VRS(:) !Volumes of proportions in the stretch
	real,ALLOCATABLE:: ET(:,:) 	!Total evapotraspiration
	real,ALLOCATABLE:: CAF(:,:)		!Capillary flow upward
	real,ALLOCATABLE:: W(:,:) 	!mount of of water over the soil
	real,ALLOCATABLE:: QBAS(:),QINT(:),QSUP(:)		!Flows in the minibasins
	real,ALLOCATABLE:: VBAS(:),VINT(:),VSUP(:)		!Volumes in the minibasins
	real,ALLOCATABLE:: TONTEM(:)		!Temperature the day before
	real,ALLOCATABLE:: P(:) !Rainfall at the interval in the cell
	real,ALLOCATABLE:: TA(:),UR(:),VV(:),SOL(:),PA(:)	!Temperature, humidity, wind, insolation, pressure
	integer HORAINI !Time that the simulation starts
	integer JDIA !Julian day
	CHARACTER (40) TITULO !Variable alpha  to read titles in the input files
	real DTP !Time interval of main calculation (in seconds given in the input file)
	real PEVAPO !Parameter that multiplies the evaporation
	!Variables related to radiation
	real GLAT !Latitude in decimal degrees (- south)
	real SOLX,T1,TAR2 !Insolation (hours of sun) and temperature (°c)
	real ALBX !ALBEDO
	real RLX !Resulting net radiation  (MJ/m2/dia)
	real RDIA !Julian day (real)
	real URX
	real,PARAMETER:: MESP=1000.0 !Water specific mass (KG/M3)
	double precision,PARAMETER:: PI=3.14159265359 !pi
	real,PARAMETER:: STEBOL=4.903E-9 !CONST. STEFAN BOLTZMANN
	!Internal variables of the subroutine
	real CLATE !latent heat of vaporization
	real SDECL,SANG !solar declination, angle of birth
	real HIM !maximum duration of sunshine (hours)
	real DR ! Relative distance Earth - Sun
	real GS,STO !Heat flux to the soil, Rariadion in the top of the atmosphere
	real SSUP,SN !Heat flux
	real ED,ES !Vapor pressures (real, Saturation)
	real SLONG !Longwave radiation

    integer CRU_INDEX !Indicates whether the climatological averages come from CRU or INMET data format    
	!Variables from Subroutine CELULA
	!Auxiliaries
	real WX,WMX,PAX,TAX,VVX,BX,KIX,KBX,RIAFX,ZX,EX,PX,SIX,VVX2
	real XL,T2,RSX,CAPX,ETX,WCX,VBX,VIX,VSX,PPU,WXB,SIB,ETB,DCAP

	!Variables from Subroutine EVAPO
	real D !Vapor pressure deficit
	real DEL !Derivative of the function ES x T
	real GAMA !Psychrometric constant
	real MSPA !Density of air
	real RUG !Roughness
	real SF !IDEM
	real SIL !Maxima intercepted amount of water
	real EIX !Amount of water evaporated from interception (POTENTIAL)
	real EIXP !Amount of water evaporated from the water surface  (POTENTIAL)
	real,ALLOCATABLE:: EVQ(:) !Direct evaporation from liquid surfaces (Water URH) from the minibasin in M3/S
	real REIX !Amount of water evaporated from interception  (real)
	real FDE !Fraction of available evaporative demand
	real,PARAMETER:: CSPA=0.001013 	!Heat especific from the humid air  MJ/kg/C
	real RA !Aerodynamic resistance
	real RAEP !Aerodynamic resistance to potential evapotranspiration
	real WL !Limit soil moisture above which no soil moisture restricts evapotranspiration (Shuttleworth)
	real WPM !Wilting point  (SHUTTLEWORTH)
	real F4 !Correction factor for surface resistence as a function of soil moisture deficit

	integer IB,IU !Counter of Sub-basins and HRUs
	
	integer ITCHUVA

	real,ALLOCATABLE:: QCONTORM(:,:),QCONTORJ(:,:) !Flow of boundary condition of MUSKINGUM CUNGE
	real,ALLOCATABLE:: QRIOINI(:,:) !Initial condition for MUSKINGUM CUNGE in the river stretch

	real,ALLOCATABLE:: BPLAN(:) !Width of the floodplain (includes own river)
	real,ALLOCATABLE:: HCALHA1(:),HCALHA2(:) !Depht that starts and that the floodplain is totally flooded
	real,ALLOCATABLE:: QMUP(:,:),AMUP(:,:),BMUP(:,:),CMUP(:,:) !Table muskingum non-linear
	integer,ALLOCATABLE:: ICODMUSK(:) !Code that indicates whether linear (0) or nonlinear (1)

	integer,allocatable:: OD(:) !RP
	real,allocatable:: CalibFLAG(:,:)
	real,allocatable:: CBOLD(:),PLAMOLD(:,:),CAPOLD(:,:),WCOLD(:,:)
	integer,allocatable::	sbtFLAG(:)	! RP

	! Variables to calculate indicators of maximum evapotranspiration:
	real,allocatable:: E0agua(:),E0TOPO(:),E0SUP(:)
	integer :: flagaclimed   !To indicate that we are working with CRU
	integer,ALLOCATABLE:: DIAH(:),MESH(:),ANOH(:),HORAH(:) !Date and time corresponding to the time interval counter

    !Forecasting module extra variables
    !FMF 21/06/2015
    integer,ALLOCATABLE:: SUBJUS(:),ATIVABACIA(:) !Sub-basin downstream from the actual one, flag to activate sub-basin
    integer,PARAMETER:: FILTUD2=1329 !forecasting module discharge results
	integer,PARAMETER:: FILTUD3=1330 !incremental discharges from the forecasting module

	
	!MANNING !FMF 09/09/2015 
	real,ALLOCATABLE:: RUGMAN(:) !Manning coefficient
    
    real ACELX
    
    real WMMINI
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !inertial variables
    integer, dimension(:,:), allocatable::cell_delta_info !if cell belongs to a delta, it will contain the delta id (>0), if not it will contain 0
    integer::n_deltas, n_specific_outlets, n_flood_plain_points_max
    integer, dimension(:), allocatable::specific_outlet_face_index
    real, dimension(:), allocatable::delta_outlet_value, specific_outlet_value

    
    !real*8,PARAMETER:: nMan=0.030             !Rugosidade considerada no modelo Inercial !FMF 09/09/2015 
    real*8:: alpha              !Coeficiente alpha utilizado para calculo do dt ideal do modelo Inercial
    real*8,PARAMETER:: g=9.81                 !Coeficiente alpha utilizado para calculo do dt ideal do modelo Inercial
    
    real*8,ALLOCATABLE:: nMan(:)            !Rugosidade considerada no modelo Inercial !FMF 09/09/2015 
    
    integer,ALLOCATABLE:: HDFLAG(:)             !Código do modelo hidrodinâmico ou inercial
    integer hdFLAG0                             !Flag do modelo hidrinâmico ou inercial
    integer,ALLOCATABLE:: MINIMONT(:,:)         !MATRIZ COM RELAÇÕES TOPOLÓGICAS DE CADA MINI-BACIA
    integer,ALLOCATABLE:: NPFL(:)               !NUMERO DE PONTOS DA TABELA COTA-AREA EM CADA MINI-BACIA
    real*8,ALLOCATABLE:: ZFUNDOFL(:)           !COTA DO FUNDO DA PLANICIE DA TABELA COTA-ÁREA
    real*8,ALLOCATABLE:: ZFL(:,:)              !COTA DA PLANICIE PARA TABELA COTA-AREA
    real*8,ALLOCATABLE:: AFL(:,:)                 !ÁREA DA PLANICIE PARA TABELA COTA-AREA
    real*8,ALLOCATABLE:: HRIO(:)                  !PROFUNDIDADE DE CALHA CHEIA DO RIO
    real*8 HRX                                    !PROFUNDIDADE DE CALHA CHEIA DO RIO (VARIÁVEL AUXILIAR) 
    real*8,ALLOCATABLE:: ZTAB(:,:)                !COTA PARA TABELA COTA-VOLUME DE CADA MINI-BACIA
    real*8,ALLOCATABLE:: VTAB(:,:)                !VOLUME PARA TABELA COTA-VOLUME DE CADA MINI-BACIA
    real*8,ALLOCATABLE:: ATAB(:,:)
    real*8 dtfloodmax,dtflood,dtflood0,tflood     !Variáveis relacionadas ao intervalo de tempo do modelo inercial
    real*8,ALLOCATABLE:: dtfloodIC(:)
    real*8 hmaxfl                                 !Variável que recebe a profundidade máxima do vetor Hfl
    real*8,ALLOCATABLE:: Q2fl(:),Vel2fl(:)        !Vazão e velocidade calculada pelo modelo inercial em cada minibacia
    real*8,ALLOCATABLE:: Qmont(:),Vol2(:),Vol1(:) !Vazão a montante e Volumes no tempo t e t+1 em uma determinada minibacia
    real*8:: SumQup
    real*8,ALLOCATABLE:: Area2(:)
    real*8,ALLOCATABLE:: Hfl(:),Yfl(:)            !Profundidade e Nível de água em cada minibacia
    real*8 nfroude                                !Numero de Froude para testes de regime supercritico
    real,ALLOCATABLE:: wwm_mean(:) !average soil moisture saturation in the InlandDelta
                                                
    !Variáveis da rotina discharge
    real*8 z1,y1,z2,y2,Sflow,hflow
	real*8 dxflow,bflow,q0,q, xMan
	integer iCJus
    
    !Variáveis da rotina continuity
    integer Nentradas, Kent, Jent,itab1,itab2,imeio
    real*8 y2_fl 
    
    integer :: nface,iFACE,KCAT,KCAT2                             !NUMERO DE PONTOS DA TABELA DE FACES
    integer,allocatable::nFACECAT1(:),nFACECAT2(:)
    real*8,ALLOCATABLE::Q2face(:),nFACEDX(:),Q2viz(:)                     !Vazão nas faces
    integer,allocatable:: jtab(:)
    
    integer,ALLOCATABLE:: INLAND_DELTA(:) !Flag se minibacia está no inland delta 
    integer,ALLOCATABLE:: NORTHERN_DELTA(:) !Flag se minibacia está no inland delta 
    integer,ALLOCATABLE:: SOUTHERN_DELTA(:) !Flag se minibacia está no inland delta 
    
    real,ALLOCATABLE:: DINFILT(:) !INFILTRATION FROM FLOODPLAIN TO SOIL
    real DINFILTX !INFILTRATION FROM FLOODPLAIN TO SOIL
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calibration variables
	integer ISEED
	integer KOUNTF,ISHUFFLE,ISOMA,IDOMIN,IDD,NPC
	real FMIN,SPARET,RMAX,XRAND,PINI,DIFPAR
	integer IREJECT,LRUIM
	real,ALLOCATABLE:: FO(:,:)
	integer, ALLOCATABLE:: IPARET(:),IRUIM(:)
	real, ALLOCATABLE:: PPAR(:,:) 
	real, ALLOCATABLE:: PAR(:,:),PARX(:)
	real, ALLOCATABLE:: PMIN(:),PMAX(:)
	real, ALLOCATABLE:: FOLD(:)
	real, ALLOCATABLE:: MEDIA(:)
	real,ALLOCATABLE::XPAR(:)
	real,ALLOCATABLE:: SPAR(:,:,:),FPLEX(:,:,:)
	real,ALLOCATABLE:: SOMAPAR(:),REFLEX(:),CONTRA(:)
	real,ALLOCATABLE:: PROB(:)
	real,ALLOCATABLE:: VMIN(:) !Minimum value of the objective functions
	integer:: NPAR !Number of parameters of the function to be optimized
	integer:: NS !Number of points in the initial sample
	integer:: NF !Number of objective functions to be optimized
	
	integer,ALLOCATABLE:: p_Calib(:,:),iFO(:),IBCONGEL(:)
	integer:: iMaxGen,NCONGEL
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! variables added for assimilation with python scheduler
    ! insertion of parameters in control variables 
    ! V. Pedinotti
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::count_params_to_assim
    character(len=200),dimension(:),allocatable::list_params_in_control
    
end module global_variables
