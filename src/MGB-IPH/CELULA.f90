    !*********************************************************************************
    !
    !  SUBROUTINE CELULA is the routine for catchment water balance in MGB-IPH
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    It calculates volumes and discharges for each local catchment.
    !    The routine has an outer loop for each catchment (iC=1,nC) and an inner loop for each URH (iU=1,NU)
	!    For each URH the subroutine SOLO is called and returns in-soil water fluxes	
	!
	!    It also uses subroutines RADIACAO and EVAPO to update vertical balance in each URH.
	!	 RADIACAO updates incoming radiation used for evap indicators at end of this routine.
	!
	!    Catchment volumes are then updated from vertical balance and discharge is calculated
    !    using linear reservoirs (i.e. surface, subsurface and groundwater. Final volume
    !    is updated from linear reservoir discharge.
	!
	!	 Main variables calculated here are:
	!		- catchment discharge/lateral inflow: QSUP, QINT and QBAS
	!		- linear reservoir volumes: VSUP,VINT and VBAS
	!
	!	 * iC and iU should not be changed inside subroutines!
	!
	!    *This routine limits min(VSUP)=0. from Parana River Basin work
    !
    !  Usage:
    !
    !    call CELULA
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
    !    * module     vars_main   in      vars_main.f90
    !    * module     vars_inerc  in      VARSINERC.f90  !Quais?  
    !    * subroutine SOLO        in      SOLO.f90 	
	!    * subroutine RADIACAO    in      RADIACAO.f90 	
	!    * subroutine EVAPO       in      EVAPO.f90 	
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
    !  2015.06.21 - 21 June 2015 (By: Fernando Mainardi Fan)    
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
    !
    !---------------------------------------------------------------------------------
	SUBROUTINE CELULA	       	                 
	
	use global_variables
	use julianday_mod
	
	implicit none
	
	real TKB(NC),TKI(NC),TKS(NC)    ! Delay and concentration times	
	real DSUP,DINT,DBAS             ! Drainage water depths from vertical balance (i.e. surface, subsurf,gw)
	real DSB,DIB,DBB				! Catchment contribution to sub-basin water by source (i.e. surface, subsurf, gw)
	!Auxiliar variables
	!	real WX,WMX,PAX,TAX,VVX,URX,SOLX,BX,KIX,KBX,RIAFX,ALBX,ZX,EX,PX,RLX,SIX
	!	real XL,T1,T2,RSX,CAPX,ETX,WCX,VBX,VIX,VSX,GLAT,PPU,WXB,SIB,ETB,DCAP
	!OUTRAS
	!	real CLATE,ED,ES
	!integer IB,IU
	integer JULDAY
    
	! Returns month, day and year to the actual julian day - will use month for climatology
	call CALDAT(JDIA,IMES,IDIA,IANO)
	JDIA=JDIA-floor(julianday(IANO,1,1))+1


	! Catchments loop
	DO IC=1,NC	 
        
        !SIMULATE ONLY RED FLOODS
        !IF(IBAC(IC)<5.AND.IBAC(IC)>7) cycle
        !IF(IBAC(IC).NE.4) cycle
        !IF(IBAC(IC)>4) cycle
        
        !IF((IBAC(IC)>4).AND.(IBAC(IC).NE.9)) cycle
        
        !IF(IBAC(IC)>9) cycle
        !IF(IBAC(IC).NE.9) cycle
        !IF(IBAC(IC).NE.5) cycle
        !IF(IBAC(IC)<10.OR.IBAC(IC)>14) cycle   !modif cécile
  
		if (sbtFLAG(iC)==1) cycle ! if catchment uses data substitution (e.g. dam operation), cycle.

		IB=IBAC(IC)       
	
		EVQ(IC)=0.0 !Direct evaporation from open waters
        
		
		! URH (GRU) Loop
		DO IU=1,NU

			IF(PUSO(IC,IU).LT.0.0001)THEN   ! iu-th URH not exist in ic-th catchment
				CYCLE
			endIF
            
            !Temporary variables
 			PX=P(IC)              ! Loads rain
			WX=W(IC,IU)			  ! Soil Water storage
			SIX=SI(IC,IU)		  ! Interceptation
			WMX=WM(IB,IU)		  ! Maximum soil water storage
			PAX=PA(IC)            ! air pressure
			TAX=TA(IC)			  ! air temperature
			VVX=VV(IC)			  ! wind speed
            ! relative humidity
            if (is_congo) then
                URX=UR(IC)/1.2
            else
                URX=UR(IC)
            end if
			SOLX=SOL(IC)		  ! solar radiation (or number of hours-insolation)
			XL=PLAM(IB,IU)		  ! lambda parameter from ARNO model
			T1=TONTEM(IC)         ! day before air temperature
			T2=TA(IC)             ! air temperature
			RIAFX=RIAF(cell_region_id(ic),IU,IMES)   ! leaf-area index
			ALBX=ALB(cell_region_id(ic),IU,IMES)     ! albedo
			ZX=Z(cell_region_id(ic),IU,IMES)		  ! vegetation height
			RSX=RS(cell_region_id(ic),IU,IMES)		  ! surface resistnace
			BX=B(IB,IU)	          ! b parameter from ARNO model
			KIX=KINS(IB,IU)       ! Subsurface water flow
			KBX=KBAS(IB,IU)       ! Groundwater flow
			CAPX=CAP(IB,IU)	      ! Capilarity
			WCX=WC(IB,IU)	      ! Soil Storage Water limit for capilar ascension
            

            

		!	IF(PAX.GT.130.0.OR.PAX.LT.70.0)STOP 'ERROR PRESSURE SUB CELULA'
			IF(TAX.LT.-50.0.OR.TAX.GT.50.0)STOP 'ERROR AIR TEMP SUB CELULA'
			IF(VVX.LT.0.0.OR.VVX.GT.100.)STOP 'ERROR WIND SUB CELULA'
			IF(URX.LT.0.0.OR.URX.GT.100.01)STOP 'ERROR UMIDITY SUB CELULA'
		!	IF(SOLX.LT.0.0.OR.SOLX.GT.24.)STOP 'ERRO INSOLATION SUB CELULA'

			! Calculates Net Radiation			
			GLAT=Y(IC)      !Loads catchment latitude
			call RADIACAO   !Returns RLX - Net radiation
            



			! Calculates evapotranspiration and interceptation (updates PX)
			call EVAPO
            


!			ET(IC,IU)=ETX

!Linha para testar modelo com evapo zero:
!			EIX=0.0
!			EX=0.0
!			ET(IC,IU)=0.0
!			PX=P(IC)

            ACELX=ACEL(IC)
			! Soil water balance 
			IF(WMX.GT.0.001) THEN
				call SOLO(PX,EX,WX,WMX,BX,KIX,KBX,XL,DSUP,DINT,DBAS,CAPX,WCX,DTP,hdflag(ic),area2(ic),is_congo,ACELX,dinfiltx)   !
				ETX=EX+REIX
			ELSE
				IF (HDFLAG(iC)==0) THEN            !if open water (Wm=0.0), surface flow equals interceptation
					DSUP=PX
					DINT=0.0
					DBAS=0.0
					CAPX=0.0
					EVQ(IC)=(EIXP*1000.*ACEL(IC)*(PUSO(IC,IU)/100.))/DTP !Direct evaporation from open waters in m3/s
				ELSE	
					!*******
					! Considers evaporation in floodplain water !RP
					!DSUP=PX 
					DINT=0.0
					DBAS=0.0
					CAPX=0.0
					!EVQ(IC)=(EIX*1000.*ACEL(IC)*(PUSO(IC,IU)/100.))/DTP
					EVQ(IC)=EIXP
					!*******
				endIF
				ETX=EIX
			endIF

			!IF (HDFLAG(iC)==1) EVQ(IC)=EIXP   !PARA GARANTIR QUE TODAS AS MINIBACIAS TEM EVQ(IC)
			PPU=PUSO(IC,IU)/100.
			DINFILT(IC)=DINFILT(IC)+DINFILTX*PPU   !mm/dtp
                
			ET(IC,IU)=ETX

			W(IC,IU)=WX
			SI(IC,IU)=SIX
			CAF(IC,IU)=CAPX


			 	!HRU area  fraction
			WXB=WXB+WX*PPU 			!Average Water Soil-Storage in catchment
			SIB=SIB+(P(IC)-PX)*PPU  !Average Interceptation in catchment
			ETB=ETB+ETX*PPU 		!Average Evapotranspiration in catchment
			WMMINI=WMMINI+WMX*PPU   !Average Wm in the catchment
            
			!Fix units; multiplies by area and HRU fraction
			! Converts DSUP, DINT and DBAS from mm/dtp to m3/s			
			!DSUP=((DSUP*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(DTP)	!DSUP=((DSUP*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			!DINT=((DINT*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(DTP)	!DINT=((DINT*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			!DBAS=((DBAS*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(DTP)	!DBAS=((DBAS*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			if (HDFLAG(iC)==0) then
			DSUP=((DSUP*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(DTP)	!DSUP=((DSUP*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			DINT=((DINT*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(DTP)	!DINT=((DINT*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			DBAS=((DBAS*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(DTP)	!DBAS=((DBAS*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			else
			DSUP=((DSUP*1000.)*max((ACEL(IC)-Area2(IC)),0.0)*(PUSO(IC,IU)/100.))/(DTP)	!DSUP=((DSUP*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			DINT=((DINT*1000.)*max((ACEL(IC)-Area2(IC)),0.0)*(PUSO(IC,IU)/100.))/(DTP)	!DINT=((DINT*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			DBAS=((DBAS*1000.)*max((ACEL(IC)-Area2(IC)),0.0)*(PUSO(IC,IU)/100.))/(DTP)	!DBAS=((DBAS*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
		   


            endif
			!Capilar flux - converts from mm/dtp to m3/s
			DCAP=DCAP+(PUSO(IC,IU)/100.)*(CAPX*1000.*ACEL(IC))/(DTP) !DCAP=DCAP+(PUSO(IC,IU)/100.)*(CAPX*1000.*ACEL(IC))/(3600.*24.) 
			
			!Updates linear reservoir volumes(m3)
			VBAS(IC)=VBAS(IC)+DBAS*DTP   !VBAS(IC)=VBAS(IC)+DBAS*3600.*24.
			VINT(IC)=VINT(IC)+DINT*DTP   !VINT(IC)=VINT(IC)+DINT*3600.*24.
			VSUP(IC)=VSUP(IC)+DSUP*DTP   !VSUP(IC)=VSUP(IC)+DSUP*3600.*24.

		endDO !end HRU loop
	  
		! Takes water from groundwater reservoir by capilar flux
		!VBAS(IC)=VBAS(IC)-DCAP*DTP !VBAS(IC)=VBAS(IC)-DCAP*3600.*24. !ATUALIZA VOLUME
		VBAS(IC)=MAX(VBAS(IC)-DCAP*DTP,0.0) 					      !Update volumes. Fernando & Paulo fixed in 25/01/2013 - limit negative volume in linear reservoir. FMF
	
		DCAP=0.0
		
		! Calculates catchment discharges
		TKB(IC)=CB(IB)*3600. 	!CB is defined in hours in parameter file ! convert to seconds
		TKI(IC)=CI(IB)*TIND(IC) !CI is defined in parameter file
		TKS(IC)=CS(IB)*TIND(IC) !CS is defined in parameter file
		QBAS(IC)=VBAS(IC)/TKB(IC)
		QINT(IC)=VINT(IC)/TKI(IC)
		QSUP(IC)=VSUP(IC)/TKS(IC)
		
		! Updates Volumes (m3)
		VBX=VBAS(IC)-QBAS(IC)*DTP       !VBX=VBAS(IC)-QBAS(IC)*3600.*24.
		VIX=VINT(IC)-QINT(IC)*DTP       !VIX=VINT(IC)-QINT(IC)*3600.*24.
		VSX=VSUP(IC)-QSUP(IC)*DTP       !VSX=VSUP(IC)-QSUP(IC)*3600.*24.
		
		
		! Verifies drying linear reservoir
		! Surface:
		! MGB-IPH original code
!		IF(VSX.LT.0.0)THEN 		!surface
!			QSUP(IC)=VSUP(IC)/(DTP)            !QSUP(IC)=VSUP(IC)/(3600.*24.)
!			VSUP(IC)=0.0
!		ELSE
!			VSUP(IC)=VSX
!		endIF

		! Assume that surface reservoir can have negative volumes to address floodplain evaporation !RP		
		if (hdFLAG(iC)==0) then
			IF(VSX.LT.0.0)THEN 		!surface
				QSUP(IC)=VSUP(IC)/(DTP)            !QSUP(IC)=VSUP(IC)/(3600.*24.)
				VSUP(IC)=0.0
			ELSE
				VSUP(IC)=VSX
			endIF
		else

			IF(VSUP(IC)<0.0)THEN 		
				! Negative discharge is instantly taken in next time interval, without phase shift
				QSUP(IC)=VSUP(IC)/(DTP)  !QSUP(IC)=VSUP(IC)/(3600.*24.)
				VSUP(IC)=0.0
			ELSEIF(VSX<0.0)THEN
				QSUP(IC)=VSUP(IC)/(DTP)            !QSUP(IC)=VSUP(IC)/(3600.*24.)
				VSUP(IC)=0.0
			ELSE
				VSUP(IC)=VSX
			endIF		
		
		endif	

		IF(VIX.LT.0.0)THEN !sub-surface
			QINT(IC)=VINT(IC)/(DTP)             !QINT(IC)=VINT(IC)/(3600.*24.)
			VINT(IC)=0.0
		ELSE
			VINT(IC)=VIX
		endIF
		IF(VBX.LT.0.0)THEN !groundwater
			QBAS(IC)=VBAS(IC)/(DTP)             !QBAS(IC)=VBAS(IC)/(3600.*24.)
			VBAS(IC)=0.0
		ELSE
			VBAS(IC)=VBX
		endIF
		!Ends drying linear reservoirs verification

        
		!IB is the sub-basin number where catchment IC is located
		PM2(IB,IT)=PM2(IB,IT)+P(IC)/KCB(IB) 	 !update sub-basin average precipitation
		
		DSB=(QSUP(IC)*3600.*24.)/(1000.*ACEL(IC))
		DIB=(QINT(IC)*3600.*24.)/(1000.*ACEL(IC))
		DBB=(QBAS(IC)*3600.*24.)/(1000.*ACEL(IC))
		
		SIM(IB,IT) =SIM(IB,IT)+SIB/KCB(IB)  	!update sub-basin average interception
		EM(IB,IT)  =EM(IB,IT)+ETB/KCB(IB)  		!update sub-basin average evapotranspiration
		WBM(IB,IT) =WBM(IB,IT)+WXB/KCB(IB)  	!update sub-basin average in-soil water
		DSUM(IB,IT)=DSUM(IB,IT)+DSB/KCB(IB) 	!update sub-basin average surface flow
		DINM(IB,IT)=DINM(IB,IT)+DIB/KCB(IB) 	!update sub-basin average subsurface flow
		DBAM(IB,IT)=DBAM(IB,IT)+DBB/KCB(IB) 	!update sub-basin average gw flow
		
		!clear temporary variables
		SIB=0.0
		ETB=0.0
		WXB=0.0
		DSB=0.0
		DIB=0.0
		DBB=0.0
        WXB=0.0
        WMMINI=0.0
	

		! Calculates evapotranspiration indicators
		! Evaporation from open water
		E0agua(iC)=EIXP !mm/dtp
		! Evaporação using all energy on top of atmosphere (based on solar constant calculations)
		E0TOPO(iC)=1000.*STO/(MESP*CLATE)
		! Evaporação using shortwave radiation available at terrestrial surface (usind input data)
		E0SUP(iC)=1000.*SSUP/(MESP*CLATE)


	endDO !End catchment loop

		

79	FORMAT(I5,12F6.1)
	RETURN
	end
