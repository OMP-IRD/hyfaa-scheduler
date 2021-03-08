    !*********************************************************************************
    !
    !  SUBROUTINE LESUBST reads  the file with time series of stream flow to substitute simulation at selected points.  (File: QSUBST.QSB)
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine reads the file with time series of stream flow to substitute simulation at selected points.  (File: QSUBST.QSB)

	!
	!	 LESUBST is called inside 1main.
	!
	!	 Saves discharge time series: 
	!     	QLIDO(.,.)
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

	SUBROUTINE LESUBST
        !Reads discharge time series to be substituted by simulation at selected points. 
        use global_variables
        implicit none
        integer I,J,K,L
        if (trim(adjustl(ARQSUBST)) .ne. 'none') then
            open(FILSUBS,FILE=trim(adjustl(ARQSUBST)),STATUS='OLD')
            READ(FILSUBS,701)(CABE(K),K=1,NUMSUBST)
            DO IT=1,NT
            READ(FILSUBS,*)I,J,L,(QLIDO(K,IT),K=1,NUMSUBST)
            endDO
            close (FILSUBS)
        else
            print*, 'No subtitute station file ARQSUBST defined...'
            call abort()
        end if
701	FORMAT(20A10)
        RETURN
	end
