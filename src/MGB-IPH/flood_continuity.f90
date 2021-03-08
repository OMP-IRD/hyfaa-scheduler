	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine calculates the level and depth for each catchment from Continuity equation   (Inertial version).
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

	subroutine flood_continuity

	!--------------------------------------------------------------------------------------
	! Variables and parameters:
	use global_variables
	implicit none
	real*8 areajtab, yjtab
    
    !$OMP PARALLEL NUM_THREADS(8)
    !$OMP DO PRIVATE(Nentradas,SumQup,Areajtab,yjtab)
    !-------------------------------------------------------------------------------------
    
        
    do iC=1,nC
    
        IB=IBAC(IC) 
    
        !SIMULATE ONLY RED FLOODS
        !IF(IBAC(IC)<5.AND.IBAC(IC)>7) cycle
        !IF(IBAC(IC).NE.4) cycle
        !IF(IBAC(IC)>4) cycle
        !IF((IBAC(IC)>4).AND.(IBAC(IC).NE.9)) cycle
        !IF(IBAC(IC)>9) cycle
        !IF(IBAC(IC).NE.9) cycle
        !IF(IBAC(IC).NE.5) cycle
        !IF(IBAC(IC)<10.OR.IBAC(IC)>14) cycle   !modif cécile
  
        
            !Number of upstream catchments of IC
            Nentradas = MINIMONT(iC,1)
            SumQup=0.0
            !Sum of downstream IC flows
            if(Nentradas==0)then
                SumQup=0.0
            else
                SumQup = SUM(Q2fl(MINIMONT(iC,2:1+Nentradas)))
            endif
            
            
            !>>>>>>>>>>>>>>>>>>>>>FOR FLOW SUBSTITUTION (E.G. ANSONGO)
            if (NUMSUBST .gt. 0) then
                if(iC==ISUBST(1)) then
                    Q2fl(iC)=QLIDO(1,IT)
                    Vol2(iC) = Vol1(iC)+dtflood*(QCEL2(iC)-Q2fl(iC)+Q2viz(IC))-(E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))
                else
                    Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC))-((E0agua(IC)+DINFILT(IC))*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) 
                end if
            else
                Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC))-((E0agua(IC)+DINFILT(IC))*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))    
            endif


            Vol2(iC) = max(Vol2(iC),0.0)
            
            
            !Interpolates the Area and Level from Volume
            call hunt(VTAB(:,iC),Vol2(iC),jtab(iC),ATAB(:,iC),Areajtab,ZTAB(:,IC),yjtab,NPFL(IC)+2)
            Area2(iC) = min(Areajtab,ACEL(iC))
            

            
            
            !Updates variables:
            y2_fl=max(yjtab,ZTAB(1,IC))
            y2_fl=y2_fl+0.001
            
            !Calculates depth
            Hfl(iC)=y2_fl-ZTAB(1,IC)
            Yfl(iC)=y2_fl       
            
            !Updates the volume
            Vol1(iC)=Vol2(iC)      
            
    
    enddo
    

    
    !$OMP end DO
    !$OMP end PARALLEL
	
	endsubroutine
