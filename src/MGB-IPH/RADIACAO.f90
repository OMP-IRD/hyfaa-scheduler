    !*********************************************************************************
    !
    !  SUBROUTINE RADIACAO is the routine that calculates net radiation at terr. surface
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine calculates net radiation based on energy equations using meterological
	!		data (i.e. air temperature, humidity, atm pressure, insolation, wind speed),
	!		soil/vegetation (i.e. canopy, surface roughness, albedo),
	!		date (e.g. julian day) and location
	!
	!	 RADIACAO is called inside subroutine CELULA.
	!
	!	 Saves global variables used for vertical water balance and evap indicators in CELULA.
	! 		RLX: net radiation (MJ/m2.dia)
	! 		STO: top-of-atmosphere radiation (MJ/m2.dia)
	! 		SSUP: incident radiation at terr. surface (MJ/m2.dia) - only insolation discounts.
	!
	!	* Cloudness correction based on Amazon Basin assumption. talk to Rodrigo.
	!
	!
    !
    !  	Usage:
    !
    !    call RADIACAO
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
    !    * module     vars_main   in      vars_main.f90
    !
    !	 opens
    !
    !    * Does not open files
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
    !    2014.25.11 - 25 November 2014 (By: Mino V. Sorribas)    
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
    !   *Variables declarations and routines calls are all commented below.
	!	* All variables are global!?
    !
    !---------------------------------------------------------------------------------		
	
	SUBROUTINE RADIACAO
	use global_variables
	implicit none
	!real GLAT 				! latitude in decimal degrees (negative for south)
	!real SOLX,T1,T2 		! insolation (number of hours), and air temperature (celsius)
	!real ALBX 				! albedo
	!real RLX 				! net radiation (MJ/m2/dia)
	!integer JDIA 			! julian day - integer
	!real RDIA 				! julian day - float
	!real URX
	!real,PARAMETER:: MESP=1000.0 ! water specific weigth (kg/m3)
	!real,PARAMETER:: PI=3.141592 ! pi
	!real,PARAMETER:: STEBOL=4.903E-9 ! stefan-boltzmann constant
	
	!internal variables
	!real CLATE 			!latent heat (vap.)
	!real SDECL,SANG 		!solar declination (angle at sunset)
	!real HIM 				!maximum duration of insolation (hours)
	!real DR 				! relative distance [earth - sun] 
	!real GS,STO 			! ground heat flux and top-of-atmosphere radiation
	!real SSUP,SN 			! heat fluxes
	!real ED,ES 			! air vapor pressure (real and saturation)
	!real SLONG 			! short-wave radiation

	! latent heat of vaporization (MJ/Kg) - ch4. Handbook of Hydrology
	CLATE=(2.501-0.002361*T2) !usa temperatura em C

	! Saturation Vapor pressure
	! ES=0.6108*EXP((12.27*T2)/(237.3+T2)) !in KPa
	ES=0.6108*EXP((17.27*T2)/(237.3+T2)) !in KPa  ! test modif ETP Cécile (from Ayan)
	
	! Real Vapor Pressure
	ED=ES*URX/100.0 		!in KPa
	
	! - ch4. - Handbook of Hydrology
	! Solar Declination (rad)
	RDIA=JDIA
	SDECL=0.4093*SIN((2*PI/365.0)*RDIA-1.405)
	
	! Sunset hour angle (rad)
	SANG=ACOS(-TAN(GLAT*PI/180.)*TAN(SDECL)) !intrinsic function TAND is for degree arguments
	
	! maximum insolation (hours)
	HIM=24*SANG/PI
	! - ch4. - !
	
	! Distance Earth-Sun
	DR=1.0+0.033*COS(2*PI*RDIA/365.0)
	
	! Top-of-atmosphere Radiation (mm/day)
	STO=15.392*DR*(SANG*SIN(GLAT*PI/180.)*SIN(SDECL)+COS(GLAT*PI/180.)*COS(SDECL)*SIN(SANG))
	
	! Convert to MJ/m2.day
	STO=STO/1000.0 		! to meters/day
	STO=STO*MESP*CLATE  ! to MJ/m2.day

	!Calculates incoming radiation to terr. surface		
	IF (flagaclimed == 1) then			! from insolation
        if (CRU_INDEX==1) then
            SSUP=STO*(0.25+0.5*(SOLX)/100.)  ! Sunshine data in % for CRU data
        else
            SSUP=STO*(0.25+0.5*SOLX/HIM) 	! (eq. 4.2.6 from Handbook Hydrology)
        endif
	ELSE		!from incident shortwave on terr surface data
        if (is_congo) then
            SSUP=STO*(0.25+0.5*SOLX/HIM) !Climate daily data
        else
            SSUP=SOLX	
        end if
    endIF

    ! Discount reflection (albedo)
	SN=SSUP*(1-ALBX)
	
	! Long-wave losses
	SLONG=STEBOL*(T2+273.2)**4.0 	  ! Temperature is converted to Kelvin
	SLONG=SLONG*(0.34-0.14*(ED)**0.5) ! Emissivity correction
	
	! Cloudness
	! SLONG=SLONG*(0.9*SOLX/HIM+0.1)	
	! SOLX is in terms of solar radiation, then uses Handbook of Hydrology eq.
	! Parameters from wetlands to address Amazon basin
	SLONG=SLONG*(1.0*SSUP/STO+0.0) 	  ! Cloudness correction


	! Ground heat flux !MJ/m2.dia
	GS=0.38*(T2-T1) 				  ! eq. 4.2.18 (Handbook of Hydrology)
	
	! Net Radiation (Energy Available) (MJ/m2.dia) - discount longwave and groun heatflux	
	RLX=SN-SLONG-GS
	

	
	RETURN
	end
