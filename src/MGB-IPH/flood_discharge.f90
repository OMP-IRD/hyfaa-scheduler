	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine calculates the flow for each catchment from Inertial equation   (Inertial version).
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

	subroutine flood_discharge

	!--------------------------------------------------------------------------------------
	! Variables and parameters:
	use global_variables
	implicit none
    integer K ,i_delta, i_outlet					!indexes and counters
    logical::is_delta_outlet
    real*8,PARAMETER:: dxart=5.               !dx artificial caso o comprimento do rio de uma minibacia seja muito pequeno Ver flood_timestep => used for congo and not for niger
    
    !$OMP PARALLEL NUM_THREADS(8)
    !$OMP DO PRIVATE(y1,y2,z1,z2,iCJus,hflow,dxflow,bflow,xMan,q0,Sflow,q)

	
    !-------------------------------------------------------------------------------------
    do iC=1,nC

        IB=IBAC(IC) 

            ! Bottom level and water level of IC catchment:
            z1=ZTAB(1,iC)
            y1=Hfl(iC)+z1

            ! Bottom level and water level of downstream IC catchment:
            iCJus = CELJUS(iC)
            
            if(iCJus == -1)then
                
                if (is_congo) then
                    z2=z1-0.0005*dxart*1000.
                    y2=y1-0.0005*dxart*1000.
                else
                    z2=z1-0.0001*SRIO(IC)*1000. !assuming a 0.0001 m/m slope at the downstream boundary
                    y2=y1-0.0001*SRIO(IC)*1000.
                end if
                
                !boundary condition for outflow in ocean
                !z2=z1
                !y2=0.0
            else
                z2=ZTAB(1,iCJus)
                y2=Hfl(iCJus)+z2
            endif

            
            ! Calculates the hflow variable:
            hflow=max(y2,y1)-max(z2,z1)
            hflow=max(hflow,0.0)
            
           
            if(iCJus /= -1)then
                dxflow=DBLE(SRIO(IC)*1000.) + DBLE(SRIO(iCJus)*1000.)       
                dxflow=dxflow/2.
            else
                if (is_congo) then
                    dxflow=DBLE(SRIO(IC)*1000.) + DBLE(dxart*1000.)
                    dxflow=dxflow/2.
                else
                    dxflow=DBLE(SRIO(IC)*1000.)
                end if
            endif
            

            !River width
            bflow=DBLE(BRIO(iC))
                      
            !River manning coefficient
            xMan=nMan(iC)
            
            ! Flow in the last time-step
            q0=Q2fl(iC)/bflow ! in m2/s
            
                    
            ! Water table slope:
            Sflow=-(y1-y2)/dxflow
            
            ! Calculates flow from Inertial equation (m2/s) for each IC:
            if (hflow>0.0) then
                q=(q0-(g*dtflood*hflow*Sflow))/(1.+g*dtflood*hflow*xMan*xMan*abs(q0)/(hflow**(10.0/3.0)))
                q=q*bflow !(m3/s)
            else
                q=0.0;
            endif
            
            ! Updates variable flow inertial:
            Q2fl(iC)=q ! in m3/s
            QJ2(iC)=Q2fl(iC) ! in m3/s
            

            
           
    enddo

    
    !$OMP end DO
    !$OMP end PARALLEL
    
        !Loop que calcula a vazão nas interconexões entre minibacias.
    if (.not.is_congo) then !liga conexões
    Q2viz=0.0
    !Q2face=0.0
    !$OMP PARALLEL DO NUM_THREADS(8) default (SHARED)
    do iFACE=1,nface
    
    
            is_delta_outlet = .false.
            KCAT=nFACECAT1(iFACE)
            KCAT2=nFACECAT2(iFACE)
            
            ! Nível de Fundo e Nível da Água da minibacia iC:
            z1=ZTAB(1,KCAT)
            y1=Hfl(KCAT)+z1
            z2=ZTAB(1,KCAT2)
            y2=Hfl(KCAT2)+z2

            ! Cálculo da profundidade de escoamento:
            hflow=max(y2,y1)-max(z2,z1)
            !Limitador de fundo
!            hflow=hflow-1.0
            !Correção de valores negativos
            hflow=max(hflow,0.0)
            
            !A rotina DBLE transforma a variável de entrada em um real*8
            !Média dos dx de IC e ICJUS
            dxflow=DBLE(nFACEDX(iFACE))       !Verificar se precisa de um limitador do dx
            bflow=100.0
            !bflow=DBLE(BRIO(KCAT))

            
            do i_delta=1,n_deltas
                if (cell_delta_info(i_delta,kcat)==1 .or. cell_delta_info(i_delta,kcat2)==1) then
                    if (delta_outlet_value(i_delta) > 0.) then
                        bflow = delta_outlet_value(i_delta)
                        is_delta_outlet = .true.
                    end if
                end if
            end do
            
            
            !WIDTH FOR ESPECIFIC CONNECTIONS (E.G. RIVER DEFLUENCES)
            
            !specific outlets (example on Niger: diaka distribuary at face 230 and bflow value 600)
            do i_outlet=1,n_specific_outlets
                if (iface == specific_outlet_face_index(i_outlet)) then
                    bflow = specific_outlet_value(i_outlet)
                end if
            end do
            

            xMan=nMan(KCAT)  
            !xMan=0.055 !modif Cécile - Manning pour le delta intérieur
               
            ! Vazão no tempo anterior:
            q0=Q2face(iFACE)/bflow ! em m2/s
                    
            ! Declividade da linha de água:
            Sflow=-(y1-y2)/dxflow
            
                
            ! Cálculo da vazão Inercial (por unidade de largura do rio) na face de jusante da minibacia iC:
            if (hflow>0.0) then
                q=(q0-(g*dtflood*hflow*Sflow))/(1.+g*dtflood*hflow*xMan*xMan*abs(q0)/(hflow**(10.0/3.0)))
                q=q*bflow
            else
                q=0.0;
            endif
            
            ! Calcula a nova vazão no próximo intervalo de tempo: EXCLUI FACES QUE NÃO ESTÃO LOCALIZADAS NO INLAND DELTA
            IF(is_delta_outlet) THEN
                Q2face(iFACE)=q ! em m3/s
                Q2viz(KCAT)=Q2viz(KCAT)- Q2face(iFACE)
                Q2viz(KCAT2)=Q2viz(KCAT2)+ Q2face(iFACE)
            ELSE
                Q2face(iFACE)=0
                Q2viz(KCAT)=0
                Q2viz(KCAT2)=0
            endIF
            

            
    enddo 
    !$OMP end PARALLEL DO
    endif

	endsubroutine
