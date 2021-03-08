    !*********************************************************************************
    !
    !  SUBROUTINE LEQOBS reads the file observed streamflow data (File with extension .QOB)
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine reads the file observed streamflow data (File with extension .QOB)

	!
	!	 LEQOBS is called inside 1main.
	!
	!	 Saves discharge time series: 
	!     	QOBS(.,.)
	!
	!
	!
    !
    !  	Usage:
    !
    !    * no subroutine is called in this subroutine
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
    !    * File with extension .QOB containing time series of observed stream flow. 
    !
    !    reads
    !
    !    * File with extension .QOB containing time series of observed stream flow
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
    !    2014.26.11 - 25 November 2014 (By: Rodrigo Paiva)    
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

	SUBROUTINE LEQOBS
		use global_variables
        implicit none
        integer I,J,K,L
        if (trim(adjustl(ARQOBS)) .ne. 'none') then
            open(FILOBS,FILE=trim(adjustl(ARQOBS)),STATUS='OLD')
            READ(FILOBS,701)(CABE(K),K=1,NOBS)
            DO IT=1,NT
                READ(FILOBS,*)I,J,L,(QOBS(K,IT),K=1,NOBS)
            endDO
            close (FILOBS)
        else
            print*, 'No observation file ARQOBS defined...'
        end if
701	FORMAT(1000A10)
        RETURN
	end
