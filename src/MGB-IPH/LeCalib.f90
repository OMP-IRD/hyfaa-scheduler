    !*********************************************************************************
    !
    !  SUBROUTINE LECALIB Reads information file for automatic calibration ParCalib.mgb
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine reads information file for automatic calibration ParCalib.mgb

	!
	!	 LECALIB is called inside 1main.
	!
	!	 Saves all variables related to automatic calibration: 
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
	!    * module     vars_calib   in      vars_calib.f90 
	!
    !	 opens
    !
    !    * file for automatic calibration ParCalib.mgb. 
    !
    !    reads
    !
    !    * file for automatic calibration ParCalib.mgb
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
    !    2015.27.05 - 27 MAY 2015 (By: Rodrigo Paiva)    
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
	subroutine LeCalib
	! Reads information file for automatic calibration ParCalib.mgb
	use global_variables
	implicit none
	! Local Vars:
	! Variáveis locais:
	integer:: i,j,k
	character(60):: trashTex
	character(10):: trashTex2
	!----------------------------------------------------------------------------
	print*,'Read calibration parameters:'
	open(FILCALIB,FILE=trim(adjustl(param_calib_file)),STATUS='OLD',ACTION='READ')
	READ(FILCALIB,75) trashTex
	print*,trashTex
	READ(FILCALIB,75) trashTex
	print*,trashTex
	READ(FILCALIB,*) NS ! Population size MOCOM-UA
	print*,NS
	read(*,*)
	READ(FILCALIB,75) trashTex
	print*,trashTex
	READ(FILCALIB,*) NF ! Number of objective functions
	print*, NF
	read(*,*)
	allocate(iFO(nF))
	do i=1,4
		READ(FILCALIB,75) trashTex
		print*,trashTex
	enddo	
	READ(FILCALIB,*) (iFO(i),i=1,NF) ! Objective function type
	print*, (iFO(i),i=1,NF)
	read(*,*)
	READ(FILCALIB,75) trashTex
	print*,trashTex
	READ(FILCALIB,*) iMaxGen
	print*,iMaxGen
	READ(FILCALIB,75) trashTex
	print*,trashTex
	READ(FILCALIB,*) NCONGEL ! Number of subbasins
	print*, NCONGEL
	allocate(IBCONGEL(NCONGEL))
	READ(FILCALIB,75) trashTex
	print*,trashTex
	READ(FILCALIB,*) (IBCONGEL(i),i=1,NCONGEL) ! Subbasin codes
	print*, (IBCONGEL(i),i=1,NCONGEL)
	read(*,*)
	allocate(p_Calib(NU+3,7))
	read(*,*)
	p_Calib=0
	do i=1,5
		READ(FILCALIB,75) trashTex
		print*,trashTex
	enddo
	do i=1,NU
		READ(FILCALIB,73) trashTex2,(p_Calib(i,j),j=1,7)
		write(*,73) trashTex2,(p_Calib(i,j),j=1,7)
	enddo
	do i=NU+1,NU+3
		READ(FILCALIB,73) trashTex2,p_Calib(i,1)
		write(*,73) trashTex2,p_Calib(i,1)
	enddo
	NPAR=0
	do i=1,NU+3
		do j=1,7
			if (p_Calib(i,j)==1) NPAR=NPAR+1
		enddo
	enddo
	print*,'Number of parameters:'
	print*, NPAR
	read(*,*)
	READ(FILCALIB,75) trashTex
	print*,trashTex
	READ(FILCALIB,75) trashTex
	print*,trashTex
	ALLOCATE (PMIN(NPAR),PMAX(NPAR)) ! Limits of parameters.
	do i=1,NPAR
		READ(FILCALIB,*) PMIN(i),PMAX(i),trashTex2
		print*, PMIN(i),PMAX(i),trashTex2
	enddo
	allocate (CalibFLAG(NOBS,nf)) ! RP
	READ(FILCALIB,75) trashTex
	print*,trashTex
	do i=1,NOBS
		READ(FILCALIB,*) (calibFLAG(i,j),j=1,NF)
		print*, (calibFLAG(i,j),j=1,NF)
	enddo
	close (FILCALIB)
71	FORMAT(6I10)
72	FORMAT(5A10)
73	FORMAT(A10,7I6)
74	FORMAT(A20)
75	FORMAT(A60)
76	FORMAT(I10,F10.1)
77	FORMAT(A20)
78	FORMAT(A10,1I6)
!79 FORMAT(A10,2F6)
	RETURN
	end
