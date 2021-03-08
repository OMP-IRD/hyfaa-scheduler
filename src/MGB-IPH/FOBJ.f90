    !*********************************************************************************
    !
    !  SUBROUTINE FOBJ evaluates objective functions for model assesment and calibration
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
	!    This routine calculates error functions for observed x calculated discharges
	!       time series. There are three error functions:
	!		- Nash-Sutcliff effiency
	!		- Nash-Sutcliff effiency of log-tranformed values
	!		- Bias on Volume (i.e. sum[Vcalc(i)]/sum[Vobs(i)]-1. )
	!
	!	For calibration, each error function can be weighted by a a number of subbasin(or stations)
	!		as defined in calibration input file.
	!
	!	Errors are stored in variable VFO
	!
    !  Usage:
    !
    !    call FOBJ
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
	
	SUBROUTINE FOBJ
	use global_variables
	implicit none
	
	integer NANO
	integer KB				! loop counter
	real SR,ERRVT,SRL		! for summation of errors
	real POND(NOBS,3)		! sub-basin weight for each objective-function (i.e. 3 function)
!	real POND(NOBS,nf)

	integer INIAJU 			! time interval for start of comparisons 
	real SOMAOBS 			! sum of observed values
	real SOMACAL 			! sum of calculated values
	real SOMALOGS 			! sum of log-tranformed observed values
	real SOMALCAL 			! sum of log-tranformed calculated values
	integer ITVAL 			! counts number of intervals w observed data
	real XMOBS,XLOGS 		! average of observed and log-transformed observed
	real SQDES(NOBS) 		! sum of squared deviations for observed
	real SOMQ 				! sum of squared deviations between observed and its average
	real SQLOG 				! sum of squared deviations for log-tranformed observed
	real SOLOG 				! sum of squared deviations between log-transformed observed and its average
	integer i,j

    iniaju=360		!starts only after 360 intervals (1st year for daily data)
	SQDES=0.0
	R2=-100.0
	R2L=-100.0
	ERRV=1000.0

	DO KB=1,NOBS
		IB=KB
		
		SOMAOBS=0.0
		SOMACAL=0.0
		SOMALOGS=0.0
		SOMALCAL=0.0
		SOMQ=0.0
		SQLOG=0.0
		SOLOG=0.0
		ITVAL=0
		IF(INIAJU.LT.NT)THEN
			DO IT=INIAJU,NT
				IF(QOBS(IB,IT).GE.0.0) THEN
					SOMAOBS=SOMAOBS+QOBS(IB,IT)
					SOMACAL=SOMACAL+QR(IB,IT)
					SOMALOGS=SOMALOGS+LOG(QOBS(IB,IT)+0.0001) !Sums a small value to avoid undefinition in log(0.0)
                    SOMALCAL=SOMALCAL+LOG(QR(IB,IT)+0.0001)
					ITVAL=ITVAL+1
				endIF
			endDO
            
			XMOBS=SOMAOBS/ITVAL
			XLOGS=SOMALOGS/ITVAL
			IMES=0
			NANO=0
			DO IT=INIAJU,NT
				IF(QOBS(IB,IT).GE.0.0) THEN
					SQDES(IB)=SQDES(IB)+(QOBS(IB,IT)-QR(IB,IT))**2.0
					SOMQ=SOMQ+(QOBS(IB,IT)-XMOBS)**2.0
					SQLOG=SQLOG+(LOG(QOBS(IB,IT)+0.0001)-LOG(QR(IB,IT)+0.0001))**2.0
					SOLOG=SOLOG+(LOG(QOBS(IB,IT)+0.0001)-XLOGS)**2.0
				endIF
			endDO
			
			if (ITVAL>0) then
				R2(IB)=1.-SQDES(IB)/SOMQ
				R2L(IB)=1.-SQLOG/SOLOG
				ERRV(IB)=(SOMACAL/SOMAOBS-1.0)*100.
			endif
		endIF
	endDO



	
	SR=0.0
	SRL=0.0
	ERRVT=0.0




	! Considers objective function as weighted-average from observed sub-basin data (or station) that
	! were selected in calibration. Only for calibration!
    if (icalib==0) then
		calibFLAG=1.0
	endif
	
	POND=0.0
	do kb=1,nobs
		do j=1,nf
			if (calibFLAG(kb,j)>0.0) POND(kb,j)=calibFLAG(kb,j)/sum(calibFLAG(:,j))
			    print*, '(calibFLAG(kb,j)>0.0) POND(kb,j)=calibFLAG(kb,j)/sum(calibFLAG(:,j))'
			    call abort()
		enddo
    enddo
    
	VFO=0.0
	do kb=1,nobs
		if (POND(kb,1)>0.0)  VFO(1)=VFO(1)+(1.-R2(kb))*POND(kb,1)
		if (POND(kb,2)>0.0) VFO(2)=VFO(2)+(1.-R2L(kb))*POND(kb,2)
		if (POND(kb,3)>0.0) VFO(3)=VFO(3)+ABS(ERRV(kb))*POND(kb,3)
	enddo


	RETURN
	end
