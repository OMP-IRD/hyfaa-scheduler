    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine uses the computer clock to generate the seed (random process)
    !
    ! Usage:
    !
    !   * no modules, functions, or subroutines are used in this funcation
    !
    !    uses modules and functions
    !
    !    * module PORTLIB
    !
    ! opens
    !
    !   * no files are opened in this routine
    !
    ! reads
    !
    !   * no files are read in this routine
    !
    ! creates
    !
    !   * no files are created in this routine
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
    
	SUBROUTINE SEMENTE(ISEED)

    ! Variables and parameters
	integer ISEED
	ISEED=TIME()
	RETURN
	end
