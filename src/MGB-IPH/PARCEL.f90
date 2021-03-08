    !*********************************************************************************
    !
    !  SUBROUTINE PARCEL is the routine that calculates additional catchment and river
	!			parameters necessary for routing and other operations
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine sets catchment additional parameters, mainly width and depth given by
	!		geomorphological relations. Also reference discharge for Muskingum-Cunge method
	!		and concentration time for catchment routing are calculated.
	!
	!	 PARCEL is called in MGB-IPH Main routine.
	!	
    !	 It calls subroutine INTECLIM to assing nearest climate station for each catchment
	!	 	also it calls subroutine REGION to get hydraulic parameters.
	!
    !
    !  	Usage:
    !
    !    call PARCEL
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
    !    * module     vars_main   in      vars_main.f90
    !    * module     vars_inerc  in      vars_inerc.f90  
	!    * subroutine INTECLIM    in      INTECLIM.f90 
	!    * subroutine REGION      in      REGION.f90 
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
    !
    !---------------------------------------------------------------------------------	
	SUBROUTINE PARCEL

	use global_variables
	implicit none
	integer IW
	integer AUXRG

	! Calculates Time of Concentration (i.e. time of travel) for Catchment
	DO IW=1,NC
		!XLADO=ACEL(IW)**0.5 	! catchment length estimative (old)
		XLADO=LCEL(IW) 			! main river length in catchment (in km) !RP
		!DH=HCEL(IW)-LCEL(IW)
		DH=HCEL(IW)*LCEL(IW) 	! altitude variation in main river (in meters) !RP

		TIND(IW)=3600.*((0.868*XLADO**3.0)/DH)**0.385	!Concentration Time by Kirpich eq. (ni seconds)
	endDO			

	!--------------------------------------------------------------------------------------- !PR 15SET14
	! Code for river width and depth based on Geomorphological relations assigned in REGION subroutine	
	if(1==0)then
	! Calculates width, depth and reference discharge for MC-Method
	DO IC=1,NC
		ARX=ACUR(IC)  	! Cumulative Drainage Area
		call REGION		! Subroutine calculates river width and full bank depth
		BRIO(IC)=BRX	 ! store river width
		HRIO(IC)=HRX !	 ! store full bank depth
		QREF(IC)=QRX 	 ! store reference discharge for Muskingum-Cunge calcs.
	endDO	
	endif
	
	!--------------------------------------------------------------------------------------- !PR 15SET14
	! Code for river width and depth based on Geomorphological relations read in external file
	if(1==0)then	
	open(FILRG,FILE=trim(adjustl(geomorphological_relation_file)),STATUS='OLD')
	DO IW=1,NC
        READ(FILRG,'(I12,3F16.3)')AUXRG,BRIO(IW),HRIO(IW),QREF(IW)
        !print*,IW,AUXRG
    endDO
	close(FILRG)
	!--------------------------------------------------------------------------------------- !PR 15SET14
	endif

	RETURN
	end
