	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine is the main subroutine of Inertial model  (Inertial version).
    !
    !
    ! Usage: flood_timestep, flood_continuity and flood_discharge
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

    SUBROUTINE flood_inercial
    
    ! Variables and parameters
    use global_variables
    
    implicit none
    real(8) elapsed_time 
    
    !-------------------------------------------------------------------------------------
    tflood=0.0

    do while (tflood<dtflood0)
    
        call flood_timestep
        dtflood=min(dtflood,dtflood0-tflood)
        tflood=tflood+dtflood
        !print*, 'Flood inundation iT=',iT, 100.*tflood/dtflood0,'%',', dt=',dtflood,' s'
        ! Inertial equation:
        call flood_discharge
        

        ! Continuity equation
        call flood_continuity    


    enddo

    
	return
	endsubroutine
    
    
    
