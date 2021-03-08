    !---------------------------------------------------------------------------------
    !  Discussion:
    !
    !    This sub-routine generates initial conditions to MGB-IPH model.
    !
    !    uses modules and functions
    !
    !    * module vars_main in VARSMAIN.f90
    !    * module vars_inerc (Inertial Version)
    !    
    !    reads
    !
    !    * no files are read in this routine
    !
    !    creates
    !
    !    * no files are created in this routine
    !
    !---------------------------------------------------------------------------------
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. - VER ISSO.
    !
    !  Version/Modified: 
    !
    !    2014.09.001 - 09 September 2014
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
    ! Variables and Parameters:
    ! *Variables declarations and routines calls are all commented below.
    !---------------------------------------------------------------------------------
    ! End of header
    !---------------------------------------------------------------------------------
    
	subroutine CondInic

    !  Variables and Parameters:
	use global_variables

	implicit none

	integer JUS,i,NN,TK
	!-------------------------------------------------------------------------------------

	print*, 'Starting CONDINIC'
	!Initializing variables
	PM=0.0
	PM2=0.0
	SIM=0.0  !Intercepted water
	EM=0.0   !Evapotranspiration water
	WBM=0.0  !mean stored water
	DSUM=0.0 !Runoff in the basin
	DINM=0.0 !Subsurface flow in the basin
	DBAM=0.0 !underground flow in the basin
	Vol1 = 0.
			
	!---------------------------------------------------------------------------------------------
	!Initial conditions
	
	! Intercepted volume:
	SI=0.0
	! Soil moisture:
	DO IC=1,NC
		IB=IBAC(IC)
		DO IU=1,NU
			W(IC,IU)=WM(IB,IU)*0.15 !Soil has 15% of water at the beginning of simulation 			
		endDO
	endDO
	! Volumes in Reservoirs:
	DO IC=1,NC
		IB=IBAC(IC)
		VBAS(IC)=QESP(IB)*ACEL(IC)*CB(IB)*3600.
		VINT(IC)=0.0
		VSUP(IC)=0.0

	endDO

	! Flow in reservoirs:
	DO IC=1,NC
		! Sub-catchment
		IB=IBAC(IC)
		
		! Underground flow:
		TK=CB(IB)*3600.
		QBAS(IC)=VBAS(IC)/TK
		! Subsurface flow:
		TK=CI(IB)*TIND(IC)
		QINT(IC)=VINT(IC)/TK
		! Runoff:
		TK=CS(IB)*TIND(IC)
		QSUP(IC)=VSUP(IC)/TK
	endDO

	!Temperature in the last day:
	TA=20.0
	
	! Flows Inicial Conditions
	! OBS RP Melhorar condicoes iniciais do MGB, usar QREF no lugar de vazao zero nos rios ********************

	QM2=0.0
	PMB2=0.0
	PMI2=0.0
	PMS2=0.0
	QM1=0.0
	QJ1=0.0
	QJ2=0.0
	QCEL1=0.0
	QCEL2=0.0
	QRIOINI=0.0
	QCONTORM=0.0
	QCONTORJ=0.0
	

	DO IC=1,NC
	    !Sum the flows generated in each unit-catchment
		QCEL2(IC)=QBAS(IC)+QINT(IC)+QSUP(IC)
		QCEL1(IC)=QCEL2(IC)

		IF(NSUBT(IC).GT.0)THEN	!unit-catchment with river

			!Accumulates flow
			QM2(IC)=QM2(IC)+QCEL2(IC)
			
			!Constant flow in river
			QJ2(IC)=QM2(IC)	

			!Continuity 
			JUS=CELJUS(IC)
			if (JUS>0) QM2(JUS)=QM2(JUS)+QJ2(IC) 

		ELSE !unit-catchment without river
			QJ2(IC)=QCEL2(IC)
			QM2(IC)=QCEL2(IC)

			JUS=CELJUS(IC)
			
			if (JUS>0) QM2(JUS)=QM2(JUS)+QCEL2(IC)
		endIF


		! Verifies if there is flow replacement
		do i=1,NUMSUBST 
			if (iC==ISUBST(i)) then
				QM2(JUS)=QLIDO(i,1)
				exit			
			endif
		enddo
	endDO
	QJ1=QJ2
	QM1=QM2

	!Muskingum-Cunge boundary conditions
	do iC=1,nC	
		QRIOINI(IC,:)=QM2(IC)
		QCONTORM(IC,:)=QM2(IC)
		QCONTORJ(IC,:)=QJ2(IC)
	enddo

	! Inercial boundary conditinos (only to Inercial version)
    IF (hdFLAG0>0) then
        call flood_ini
    endIF 
    !******************************************************************

	return
	end subroutine
