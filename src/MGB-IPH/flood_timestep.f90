	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine calculates the time-step of inertial model  (Inertial version).
    !
    !
    ! Usage:
    !
    ! *
    !
    ! uses modules, functions, and subroutines
    !
    ! * use global_variables
    ! * use vars_inerc (only to Inertial version)
    !
    ! opens
    !
    ! * no files are created in this routine
    !
    ! reads
    !
    ! * no files are created in this routine
    !
    ! creates
    !
    ! * no files are created in this routine
    !
    !---------------------------------------------------------------------------------
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. - VER ISSO.
    !
    !  Version/Modified: 
    !
    ! 2015.07.06 - 07 July 2015 (by Paulo Pontes)
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
	
	subroutine flood_timestep
	!--------------------------------------------------------------------------------------
	! Variables and parameters:
	use global_variables
	implicit none
	
	real*8 dtminIC
	real*8 dtteste
	!-------------------------------------------------------------------------------------
    
    dtteste=DBLE(DTP)
    

    DO IC=1,NC  !Catchment loop
        ! MAX H to calculates the time-step of Inertial model:
        hmaxfl=(Hfl(IC))
        ! Correct if equal to zero:
        hmaxfl=max(hmaxfl,0.001)
        ! Compute time interval:
        dtflood=alpha*SRIO(IC)*1000./(g*hmaxfl)**0.5
        if(dtflood<dtteste) dtteste=dtflood
    endDO

	
    ! MAX time-step of MGB input data (e.g. dialy):
    dtfloodmax=DBLE(DTP)
    dtminIC=dtteste
    dtflood=min(dtfloodmax,dtminIC)   !dt for all catchments calculated as function of H. 
 	
 	endsubroutine
