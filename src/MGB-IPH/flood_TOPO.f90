    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine identifies the catchments which draining to a catchment IC analysed  (Inertial version).
    !    Generates a matrix where each IC has a upstream catchment number and the catchments identifier
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
    
    SUBROUTINE flood_TOPO

    !Parameters and variables
    use global_variables
    
    MINIMONT=0 
 
    DO IC=1,NC 
        J=CELJUS(IC) !Identifying the downstream catchemnt of IC 
        IF(J<0)CYCLE !This catchment is a outlet
        MINIMONT(J,1)=MINIMONT(J,1)+1 
        K=MINIMONT(J,1) 
        MINIMONT(J,1+K)=IC
    endDO
        
    RETURN
    end
