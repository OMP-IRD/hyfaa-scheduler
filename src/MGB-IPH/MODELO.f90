    !*********************************************************************************
    !
    !  SUBROUTINE MODELO controls the main time loop in MGB-IPH and call routines
	!					that realize catchment and river routing
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine has the main loop of the MGB-IPH, the time loop, from iT=1,nT
    !     where iT is the time interval and nT is the number of time steps.
	!
	!	 For each time interval date info is set, then rainfall and climate data are
	!	   loaded through subroutines read_forcing. Afterwards catchment flow
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
    !    call MODELO
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
	!	 * function CALDAT	      		in	  caldat.f90
    !    * module     vars_main   		in    vars_main.f90
    !    * module     vars_inerc  		in    VARSINERC.f90  !Quais? 
 	!	 * subroutine read_forcing	    in	  read_forcing.f90
	!	 * subroutine CELULA	  		in	  CELULA.f90
	!	 * subroutine REDE	      		in	  REDE.f90	
	!	 * subroutine flood_inercial	in	  flood_inercial.f90
	!
	!    Deactivated:	
	!	 * module	  AUX_MOD        	in	  HD_AUX_MOD.f90
	!	 * subroutine Hidrodinamico2	in	  Hidrodinamico2.f90
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
    !   * Variables declarations and routines calls are all commented below.
    !
    !---------------------------------------------------------------------------------
	SUBROUTINE MODELO

	use global_variables
    use hydrological_state_restart_io
	!use AUX_MOD
	implicit none
	integer K				!indexes and counters
	integer KC,JB,JC,KB,JCB,MWP 	!... same
    integer count_daniel,dt_daniel
	integer iTwrite
    

	! Initialize
	IT=0
    Dt_daniel=INT((NT-IT)/10)
 
    !CALCULA PREC E ET MENSAL PARA INLAND DELTA
    
     ! Time Loop     
    DO WHILE (IT<NT)
		
		IT=IT+1
		
       if(mod(it,100)==0)  print*,100*IT/NT,' %'
          


		if (icalib==0.and.itWrite<iT) then
			itWrite=itWrite+1
        endif
        
		if(it==count_daniel)then
			701		FORMAT(A2,$)
			count_daniel=Count_daniel+Dt_daniel
        endif


		JDIA=IDINI+INT((IT+HORAINI-1)/(86400./DTP)) 	! Gets calendar day (big number)
		call CALDAT(JDIA,IMES,IDIA,IANO)

		DIAH(IT)=IDIA !stores day corresponding to iT time interval
		MESH(IT)=IMES !stores month corresponding to iT time interval
		ANOH(IT)=IANO !stores year corresponding to iT time interval
		HORAH(IT)=MOD(IT,24)+HORAINI 	!hour of day corresponding to iT time interval


		!Subroutine that reads forcing data (from forcing_file or climatology)
		ITCHUVA=IT
		TONTEM=TA
		call read_forcing
        
        
	    
		DO KC=1,NC
			PM(KC)=PM(KC)+P(KC) !Cumulative rainfall (used later for average)
		endDO		
	
		! Calls catchment routing/discharge for lateral inflow
		DINFILT=0.0
		call CELULA	

	
		! Saves detailed info for soil state variables in file NOSOLO.HIG.
		! Uses JC index for catchment and JB index for URH of interest.
		IF(ICALIB.EQ.0)THEN !only if not calibrating.
			JB=1
			JC=1
			!JC=57
			write(FILSOL,75)IT,P(JC),W(JC,JB),SI(JC,JB),ET(JC,JB),CAF(JC,JB),QBAS(JC),QINT(JC),QSUP(JC)
			JB=2
			JC=2
			write(FILSOL2,75)IT,P(JC),W(JC,JB),SI(JC,JB),ET(JC,JB),CAF(JC,JB),QBAS(JC),QINT(JC),QSUP(JC)			

!			write(971,66) (E0agua(iC),iC=1,nC)
!			write(972,66) (E0topo(iC),iC=1,nC)
!			write(973,66) (E0sup(iC),iC=1,nC)

		endIF

	
		! Call main river routing routine using Muskingum-Cunge Method
		if(hdFLAG0==0)then
		    call REDE
		endif
		

        
		! Calculates lateral inflow 
        do ic=1,nc
		    QCEL2(IC)=QBAS(IC)+QINT(IC)+QSUP(IC) !sums surface, subsurface and 
		enddo
        
        ! Calls river routing routine using Inertial Method
        if(hdFLAG0>0)then
            call flood_inercial
        endif   


		! Stores discharges for file ouput when in simulation model
		IF(ICALIB.EQ.0)THEN 	! only if it is not calibration
			DO KB=1,NB				! store discharge in sub-basin
				JCB=IEXUT(KB) 		! outlet catchment of KB-th subbasin
                if (jcb.gt.0) then
                    QRB(KB,IT)=QJ2(JCB)
                end if
			endDO
	
			! Store discharge by water origin (i.e. surface, subsurface, groundwater) in a specified catchment
			MWP=1 						!catchment (this is manual)
			QB(IT)=QJ2(MWP)*PJB2(MWP)
			QBI(IT)=QB(IT)+QJ2(MWP)*PJI2(MWP)
			QBIS(IT)=QBI(IT)+QJ2(MWP)*PJS2(MWP)
		endIF
        
        
        !dump hydrological state if requested
        if (trim(adjustl(hydrological_state_write_prefix_each_step)) .ne. 'none') then
            write(trash_str, *) IT
            trash_str = trim(adjustl(hydrological_state_write_prefix_each_step)) // '_' // trim(adjustl(trash_str)) // '.nc'
            print*, 'Hydrological state written to ' // trim(adjustl(trash_str))
            call write_hydrological_state(trim(adjustl(trash_str)))
        end if
        
	
	endDO !End time loop
	

	

75	FORMAT(I6,8F10.4)

66	FORMAT(F10.4)

     
	RETURN
	end
