module allocate_global_variables
	use global_variables

contains

    subroutine allocate_main_variables
        implicit none	
        save

		allocate (DIAH(NT),MESH(NT),ANOH(NT),HORAH(NT)) !Date and time corresponding to the time interval counter
		allocate (ICODMUSK(NC)) !Code that indicates linear or non-linear
		allocate (BPLAN(NC)) !Width of the floodplain (includes own river)
		allocate (HCALHA1(NC),HCALHA2(NC)) !Depth that starts the floodplain and that this floodplain is totally flooded
		allocate (QMUP(NC,NTMUP),AMUP(NC,NTMUP),BMUP(NC,NTMUP),CMUP(NC,NTMUP)) !Non-linear Muskingum table
		allocate (EVQ(NC)) !Direct evaporation of the liquid surface of the cell in m3 / s
		allocate (QRIOINI(NC,NUMUSK+1)) !Initial condition of muskingum Cunge rounting
		allocate (QCONTORM(NC+1,25),QCONTORJ(NC+1,25)) !Streamflow boundary condition for Muskingum Cunge
		allocate (SI(NC,NU))
        allocate (specific_climatology_years(n_climatology))
        allocate (TAMM(n_climatology,NC,12),URMM(n_climatology,NC,12),VVMM(n_climatology,NC,12))
        allocate (PAMM(n_climatology,NC,12),SOLMM(n_climatology,NC,12))
		allocate (PMB(NB)) !Average rainfall at each sub-basin
		allocate (KPM(NB)) !Auxiliar variable
        
        allocate(ICBOM(NC))
		
		 allocate (TD(NC,NT),UD(NC,NT),VD(NC,NT),SOLD(NC,NT))
		 allocate (PAD(NC,NT))
		!end if

		allocate (QREF(NC)) !Reference flow
	 	allocate (PUSO(NC,NU))!Proportion of uses in minibasin
	 	allocate (ACEL(NC),ACUR(NC),SRIO(NC),DECL(NC))
		allocate (X(NC),Y(NC))!Coordinates of the cell center
		allocate (IBAC(NC),HCEL(NC),LCEL(NC),CELJUS(NC))!SubBasin,HMAX,HMIN,Downstream minibasin
	 	allocate (PM(NC)) !Average rainfall and rainfall at the timestep in the minibasin
		allocate (BRIO(NC))!River width
		allocate (IEXUT(NB)) !Indicate the outlets of the basin
		allocate (CI(NB),CB(NB),CS(NB)) !Pparameters of propagation in cell
		allocate (CIOLD(NB),CSOLD(NB)) !Pparameters of propagation in cell
		allocate (PLAM(NB,NU),CAP(NB,NU),WC(NB,NU))
		allocate (WMOLD(NB,NU),BOLD(NB,NU),KBOLD(NB,NU),KIOLD(NB,NU))
		allocate (WM(NB,NU),B(NB,NU),KINS(NB,NU),KBAS(NB,NU))
		allocate (NSUBT(NC)) !Number of stretchs for MUSKINGUM CUNGE routing (diferent at each minibasin)
		allocate (DT(NC)) !Time steps for MUSKINGUM CUNGE (diferent at each minibasin)
		allocate (CEL(NC),TIND(NC))
        allocate(cell_region_id(nc))
		allocate (ALB(n_regions,NU,12),RIAF(n_regions,NU,12),Z(n_regions,NU,12),RS(n_regions,NU,12))
		allocate (QESP(NB))	 !Specific baseflow (M3/S/KM2)
		allocate (QB(NT),QBI(NT),QBIS(NT))!Streamflow according to its origin
		allocate (DSUM(NB,NT),DINM(NB,NT),DBAM(NB,NT))!Averages over the SubBasins
		allocate (SIM(NB,NT),EM(NB,NT),WBM(NB,NT),PM2(NB,NT))!Averages over the SubBasins
		allocate (KCB(NB)) !Number of minibasins inside each subbasin
		allocate (QR(NOBS,NT),QOBS(NOBS,NT))!Calculated and Observed streamflow in the outlet of each microbasin
		allocate (QLIDO(NUMSUBST,NT)) !Flow that is read to replace calculated
		allocate (QMANTEIGA(NT),QSAOROM(NT),QSAOFRA(NT),QPMCRUZ(NT),QMANGA(NT),QCARI(NT))
		allocate (QBJLAPA(NT),QGAMEL(NT),QMORP(NT))
		allocate (QRB(NB,NT)) !Hydrographs of sub-basins
		allocate (QM1(NC+1),QJ1(NC),QM2(NC+1),QJ2(NC)) !Streamflows upstream and downstream in i
		allocate (QCEL1(NC),QCEL2(NC)) !Original streamflows at minibasins at the instants t E t+1 at minibasin i
		allocate (PMB2(NC+1),PMI2(NC+1),PMS2(NC+1),PJB2(NC+1),PJI2(NC+1),PJS2(NC+1)) !Proportions source of streamflow to the river
		allocate (VRB(NC),VRI(NC),VRS(NC)) !Volumes of proportions source of streamflow to the river
		allocate (ET(NC,NU)) 	!Total evapotranspiration
		allocate (CAF(NC,NU))		!Upward capillary flow
		allocate (W(NC,NU)) 	!Amount of water in the soil
		allocate (QBAS(NC),QINT(NC),QSUP(NC))		!Flows at the minibasins
		allocate (VBAS(NC),VINT(NC),VSUP(NC))		!Volume at the minibasin
		allocate (TONTEM(NC))		!Previous day temperature
		allocate (P(NC)) !Rainfall at the interval in the cell
		allocate (TA(NC),UR(NC),VV(NC),SOL(NC),PA(NC))	!Temperature, Umidity, Wind, Insolation, Pressure
		allocate (R2(NOBS),ERRV(NOBS),R2L(NOBS)) !COEF. NASH-SUTCLIFFE, VOLUME ERROR, COEF. NASH LOGARITHMS

		!allocate (OD(NC),hdFLAG(NC)) !RP
		allocate (OD(NC)) !Inertial modification
 
		allocate (CBOLD(NB),PLAMOLD(NB,NU),CAPOLD(NB,NU),WCOLD(NB,NU)) ! RP
		allocate (sbtFLAG(nC))

		! Variables to calculate indicators of maximum evapotranspiration:
		allocate (E0agua(nC),E0TOPO(nC),E0SUP(nC))

        !FMF 21/06/2015
        allocate (ATIVABACIA(NB)) !flag to activate sub-basins or not
		
		!FMF 09/09/2015 
		allocate (RUGMAN(NC)) !Manning
        
        allocate (nFACECAT1(nface),nFACECAT2(nface),nFACEDX(nface),Q2face(nface),Q2viz(nc))
        
        
    end subroutine allocate_main_variables


    subroutine allocate_inertial_variables
        implicit none	
        save

        allocate(cell_delta_info(n_deltas,nc), delta_outlet_value(n_deltas))
        allocate(specific_outlet_face_index(n_specific_outlets), specific_outlet_value(n_specific_outlets))
        allocate (HDFLAG(NC))               !Código para rodar hidrodinâmico ou inercial
        allocate (MINIMONT(NC,10))           !MATRIZ DE LIGAÇÕES TOPOLOGICAS PARA PROPAGAÇÃO MODELO INERCIAL
        allocate (NPFL(NC),ZFUNDOFL(NC))    !Número de pontos do arquivo Cota-Area e Nível de Fundo da planície
        allocate (ZFL(n_flood_plain_points_max,NC))               !COTAS PARA TABELA COTA-AREA DA PLANICIE
        allocate (AFL(n_flood_plain_points_max,NC))               !ÁREAS PARA TABELA COTA-AREA DA PLANICIE
        allocate (HRIO(NC))                 !PROFUNDAIDE DE CALHA CHEIA DO RIO
        allocate (ZTAB(n_flood_plain_points_max+2,NC),VTAB(n_flood_plain_points_max+2,NC))  !COTAS E VOLUMES PARA TABELA COTA-VOLUME
        allocate (ATAB(n_flood_plain_points_max+2,NC))
        allocate (dtfloodIC(NC))            !Intervalo de tempo em cada minibacia
        allocate (Q2fl(NC),Vel2fl(NC))      !Vazão e velocidade calculada pelo modelo inercial em cada minibacia
        allocate (Qmont(3))                 !Vazão a montante de uma determinada minibacia
        allocate (Vol2(NC),Vol1(NC))        !Volumes no tempo t e t+1 em uma determinada minibacia
        allocate (Area2(nc))
        allocate (Hfl(NC),Yfl(NC))          !Profundidade e Nível de água em cada minibacia
        allocate (jtab(NC))
        !FMF 09/09/2015 
        allocate (nMan(NC))
        allocate(DINFILT(NC)) !INFILTRATION FROM FLOODPLAIN TO SOIL	
        allocate(wwm_mean(NT)) !W/WM AVERAGE FOR EACH FLOW GAUGE
      
    end subroutine allocate_inertial_variables
    
    
    subroutine allocate_calibration_variables
    
		allocate (IRUIM(NS))
		allocate (FO(NS,NF)) !Objective functions
		allocate (IPARET(NS))
		allocate (XPAR(NPAR))
		allocate (SOMAPAR(NPAR),REFLEX(NPAR),CONTRA(NPAR)) !Sums of COORD., reflection point, point of contraction
		allocate (MEDIA(NS)) !Average of objective functions of the population
		allocate (PPAR(NPAR,NS)) !relative value of the parameter at the validity range
		allocate (PAR(NPAR,NS)) !parameter value 
!		allocate (PMIN(NPAR),PMAX(NPAR)) !validity range of the parameters
		allocate (FOLD(NS)) !former value of the function at each point of the sample
		allocate (PROB(NS)) !probability of choice of point 1..NS
		allocate (VMIN(NF)) !minimum value of objective functions
		allocate (PARX(NPAR)) !parameter value
        
    end subroutine allocate_calibration_variables
    
    

end module allocate_global_variables
