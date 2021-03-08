    !*********************************************************************************
    !
    !  SUBROUTINE CalibParam is the routine passes MOCOM-UA multipliers in calibration
	!                        to update MGB-IPH parameters
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine gets PARX, a vector with parameters multipliers and applies in
	!    	MGB-IPH parameters for new simulation during calibration process
    !
	!	 (i) 7 parameters are URH/Sub-Basin based:
	!	 	 WM,B,KINS,KBAS,PLAM,CAP,WC
	!
	!	 (ii) 3 parameters are Sub-basin based:
	!		 CS,CI,CB
	!
	!	Multipliers are transfered from HRU flagged with -[number_of_urh] in input file
	!
	!	Some sub-basin can have parameters "frozen" so they are not affected by calibration.
	!		This is achieved by setting Frozen sub-basin in calibration input file.
	!		
    !
    !
    !  Usage:
    !
    !    call CALIBPARAM
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
    !    * module     vars_main   in      vars_main.f90
    !    * module     vars_calib  in      vars_calib.f90  
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
	!	 2015.02.02 - 
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
	
	SUBROUTINE CalibParam

	use global_variables
	implicit none

	integer:: i1,i2,i3,L,i,j,k			! loop counters
	
	! Load reference values for calibration parameters
	WM  =WMOLD
	B   =BOLD
	KINS=KIOLD
	KBAS=KBOLD
	PLAM=PLAMOLD
	CAP =CAPOLD
	WC  =WCOLD
	CS  =CSOLD
	CI  =CIOLD
	CB  =CBOLD
	L=0
	
	! Check and modify parameters
	do i1=1,NU
		do i2=1,7
			if (p_Calib(i1,i2)==0) cycle	!wotn change this parameter

			if (p_Calib(i1,i2)==1) then
				L=L+1
				selectcase (i2)
				case(1)
					WM(:,i1)=WMOLD(:,i1)*PARX(L)
				case(2)
					B(:,i1)=BOLD(:,i1)*PARX(L)
				case(3)
					KBAS(:,i1)=KBOLD(:,i1)*PARX(L)
					!KINS(:,i1)=KIOLD(:,i1)*PARX(L)	!someone changed the order from older version
				case(4)
					KINS(:,i1)=KIOLD(:,i1)*PARX(L)
					!KBAS(:,i1)=KBOLD(:,i1)*PARX(L)
				case(5)
					PLAM(:,i1)=PLAMOLD(:,i1)*PARX(L)
				case(6)
					CAP(:,i1)=CAPOLD(:,i1)*PARX(L)
				case(7)
					WC(:,i1)=WCOLD(:,i1)*PARX(L)
				endselect
			elseif (p_Calib(i1,i2)<0) then
				i3=-p_Calib(i1,i2)
				selectcase (i2)
				case(1)
					WM(:,i1)=WMOLD(:,i1)*WM(:,i3)/WMOLD(:,i3)
				case(2)
					B(:,i1)=BOLD(:,i1)*B(:,i3)/BOLD(:,i3)
				case(3)
					KBAS(:,i1)=KBOLD(:,i1)*KBAS(:,i3)/KBOLD(:,i3)
					!KINS(:,i1)=KIOLD(:,i1)*KINS(:,i3)/KIOLD(:,i3)
				case(4)
					KINS(:,i1)=KIOLD(:,i1)*KINS(:,i3)/KIOLD(:,i3)
					!KBAS(:,i1)=KBOLD(:,i1)*KBAS(:,i3)/KBOLD(:,i3)

				case(5)
					PLAM(:,i1)=PLAMOLD(:,i1)*PLAM(:,i3)/PLAMOLD(:,i3)
				case(6)
					CAP(:,i1)=CAPOLD(:,i1)*CAP(:,i3)/CAPOLD(:,i3)
				case(7)
					WC(:,i1)=WCOLD(:,i1)*WC(:,i3)/WCOLD(:,i3)
				endselect
			endif
		enddo
	enddo
	if (p_Calib(NU+1,1)==1) then
		L=L+1
		CS=CSOLD*PARX(L)
	endif
	if (p_Calib(NU+2,1)==1) then
		L=L+1
		CI=CIOLD*PARX(L)
	endif
	if (p_Calib(NU+3,1)==1) then
		L=L+1
		CB=CBOLD*PARX(L)
	endif

	! Sub-basins where parameters are not going to be modified ("frozen")
	if (NCONGEL>0) then
		do i1=1,NCONGEL
			IB=IBCONGEL(i1)
			WM(IB,:)=WMOLD(IB,:)
			B(IB,:)=BOLD(IB,:)
			KINS(IB,:)=KIOLD(IB,:)
			KBAS(IB,:)=KBOLD(IB,:)
			PLAM(IB,:)=PLAMOLD(IB,:)
			CAP(IB,:)=CAPOLD(IB,:)
			WC(IB,:)=WCOLD(IB,:)
			CS(IB)=CSOLD(IB)
			CI(IB)=CIOLD(IB)
			CB(IB)=CBOLD(IB)
		enddo
	endif
		

	RETURN
	end
