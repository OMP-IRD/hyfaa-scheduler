module static_data_loader

	use global_variables
    use allocate_global_variables
    use netcdf
    use netcdf_addons


contains

    subroutine load_static_data(static_data_file)
        implicit none
        character(len=*), intent(in)::static_data_file
        integer::var_id, ncid
        integer::ii,nc_rep
        logical::specific_climatology_year_mode
        real, dimension(1,11,12)::test1
        real, dimension(12,11,1)::test2

        print*, 'Opening static_data ncfile'
        call check_ncrequest( nf90_open(trim(adjustl(static_data_file)), NF90_NOwrite, ncid) )

        
        !check if mesh dimension matches
        print*, 'Reading main dimensions'
        call check_dimension_match(ncid, "n_cells", NC)
        nb = nc_dimension(ncid, "n_basins")
        nu = nc_dimension(ncid, "n_soil_types")
        n_regions = nc_dimension(ncid, "n_regions")
        nface = nc_dimension(ncid, "n_faces")
        n_flood_plain_points_max = nc_dimension(ncid, "n_flood_plain_points_max")
        
        nc_rep = nf90_inq_dimid(ncid, "n_deltas", trash_int)
        if (nc_rep /= nf90_noerr) then
            n_deltas = 0 !if this dimension does not exist, it means that it is 0
        else
            n_deltas = nc_dimension(ncid, "n_deltas")
        end if
        
        nc_rep = nf90_inq_dimid(ncid, "n_specific_outlets", trash_int)
        if (nc_rep /= nf90_noerr) then
            n_specific_outlets = 0 !if this dimension does not exist, it means that it is 0
        else
            n_specific_outlets = nc_dimension(ncid, "n_specific_outlets")
        end if
        
        nc_rep = nf90_inq_dimid(ncid, "n_climatology", trash_int)
        if (nc_rep /= nf90_noerr) then
            specific_climatology_year_mode = .false.
            n_climatology = 1 !if this dimension does not exist, it means that it is 0
        else
            specific_climatology_year_mode = .true.
            n_climatology = nc_dimension(ncid, "n_climatology")
        end if
        
        print*, 'Allocating variables'
        call allocate_main_variables !Allocate the model main variables
        call allocate_inertial_variables !Allocate the variables of the model Local Inertial flow routing model
        
        !climatological stations used to be defined on a separate mesh, this is not the case anymore => you need to interpolate data on mesh in pre-processing
        do ii=1,nc
            icbom(ii) = ii
        end do
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LEVAR
        !read hydrological response unit variables that used to be in ALBIAF.FIX file
        print*, 'Reading hydrological response unit variables'
        call check_ncrequest( nf90_inq_varid(ncid, "cell_region_id", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, cell_region_id) )
        call check_ncrequest( nf90_inq_varid(ncid, "albedo_climatology", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, alb) )
        call check_ncrequest( nf90_inq_varid(ncid, "leaf_area_index_climatology", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, riaf) )
        call check_ncrequest( nf90_inq_varid(ncid, "trees_height_climatology", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, z) )
        call check_ncrequest( nf90_inq_varid(ncid, "superficial_resistance_climatology", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, rs) )



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LEUSO
        !read calibration variables that used to be in PARUSO.CAL file
        print*, 'Reading calibration variables'
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_soil_max_water_capacity", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, wm) )
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_soil_arno_model_b", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, b) )
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_soil_groundwater_flow_coefficient", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, kbas) )
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_soil_subsurface_flow_coefficient", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, kins) )
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_soil_arno_model_lambda", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, plam) )
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_soil_capillarity_depth", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, cap) )
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_soil_capillarity_limit", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, wc) )
        

        
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_routing_cs", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, cs) )
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_routing_ci", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, ci) )
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_routing_cb", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, cb) )
        call check_ncrequest( nf90_inq_varid(ncid, "calibration_routing_qesp", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, qesp) )

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LECELL
        !read variables that used to be in mini.gtp file
        print*, 'Reading main hydrological variables (previously named mini.gtp)'
        call check_ncrequest( nf90_inq_varid(ncid, "longitude_center", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, X) )
        call check_ncrequest( nf90_inq_varid(ncid, "latitude_center", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, Y) )
        call check_ncrequest( nf90_inq_varid(ncid, "cell_basin_id", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, IBAC) )
        call check_ncrequest( nf90_inq_varid(ncid, "cell_area", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, ACEL) )
        call check_ncrequest( nf90_inq_varid(ncid, "total_drained_area", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, ACUR) )
        call check_ncrequest( nf90_inq_varid(ncid, "main_river_length", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, SRIO) )
        call check_ncrequest( nf90_inq_varid(ncid, "main_river_slope", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, DECL) )
        call check_ncrequest( nf90_inq_varid(ncid, "longest_river_length", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, LCEL) )
        call check_ncrequest( nf90_inq_varid(ncid, "longest_river_slope", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, HCEL) )
        call check_ncrequest( nf90_inq_varid(ncid, "cell_id_upstream", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, CELJUS) )
        call check_ncrequest( nf90_inq_varid(ncid, "hydrological_order", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, OD) )
        call check_ncrequest( nf90_inq_varid(ncid, "hydrodynamic_module", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, hdFLAG) )
        call check_ncrequest( nf90_inq_varid(ncid, "river_width", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, BRIO) )
        call check_ncrequest( nf90_inq_varid(ncid, "river_depth", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, HRIO) )
        call check_ncrequest( nf90_inq_varid(ncid, "manning_coefficient", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, RUGMAN) )
        call check_ncrequest( nf90_inq_varid(ncid, "hydrological_response_unit_coefficients", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, PUSO) )

        print*, 'Fill and check additional variables'
        !fill complementary variables (used to be in LECELL subroutine)
        nman = rugman !Fill inertial model subroutines manning 
        qref = 0.333*acur**0.7 !Calculates reference initial discharge
        
        !check that hydrological response unit percentage sum up to a complete 100%, and that there are no negative values
        do ii=1,nc
            if (minval(PUSO(ii,:))<0.0) then
                print*, 'Negative values encountered in the % of HRUs in the catchments',ii
                call abort()
            end if
            trash_dp = sum(PUSO(ii,:))
            if ((trash_dp>100.01) .or. (trash_dp<99.99)) then
                print*, 'Values encountered in the % of HRUs in the catchments',ii
                call abort()
            endif
            !set it at an exact 100% => activate after validation
!~             PUSO(ii,:) = 100.*PUSO(ii,:)/trash_dp
        end do
        
        ! Define outlet of sub-basins: Outlet equals catchment of larger ID code.
        IEXUT=-999999
        do ii=1,nc
            if (IEXUT(IBAC(ii))<ii) then
                IEXUT(IBAC(ii)) = ii
            end if
        end do
        ! Correct length of main river. Minimum length = 2.0 km:   
        do ii=1,nc
            srio(ii)=max(srio(ii),2.0)
        enddo
        ! Change units of river slope to m/m
        decl=decl/1000.0 !RP
        ! Minimun tributarie lenght set to 0.001 km
        do ic=1,nc
            lcel(ic)=max(lcel(ic),0.001)
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LECLIMED
        print*, 'Reading climatology data'
        !read climatology for forcing variables that used to be in medias.cru file
        if (specific_climatology_year_mode) then
            call check_ncrequest( nf90_inq_varid(ncid, "specific_climatology_years", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, specific_climatology_years) )
        else
            specific_climatology_years = -1
        end if
        
        call check_ncrequest( nf90_inq_varid(ncid, "temperature_climatology", var_id) )
        !check that dimensions match specific_climatology_year_mode
        if (specific_climatology_year_mode.and.(nc_ndims(ncid, var_id).ne.3)) then
            print*, 'temperature_climatology dimensions should be 3 in specific_climatology_year_mode'
            call abort()
        else if ((.not.specific_climatology_year_mode).and.(nc_ndims(ncid, var_id).ne.2)) then
            print*, 'temperature_climatology dimensions should be 2 when not in specific_climatology_year_mode'
            call abort()
        end if
        if (specific_climatology_year_mode) then
            call check_ncrequest( nf90_get_var(ncid, var_id, TAMM) )
            call check_ncrequest( nf90_inq_varid(ncid, "relative_humidity_climatology", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, URMM) )
            call check_ncrequest( nf90_inq_varid(ncid, "sunshine_climatology", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, SOLMM) )
            call check_ncrequest( nf90_inq_varid(ncid, "wind_speed_climatology", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, VVMM) )
            call check_ncrequest( nf90_inq_varid(ncid, "pressure_climatology", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, PAMM) )
        else
            call check_ncrequest( nf90_get_var(ncid, var_id, TAMM(1,:,:)) )
            call check_ncrequest( nf90_inq_varid(ncid, "relative_humidity_climatology", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, URMM(1,:,:)) )
            call check_ncrequest( nf90_inq_varid(ncid, "sunshine_climatology", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, SOLMM(1,:,:)) )
            call check_ncrequest( nf90_inq_varid(ncid, "wind_speed_climatology", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, VVMM(1,:,:)) )
            call check_ncrequest( nf90_inq_varid(ncid, "pressure_climatology", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, PAMM(1,:,:)) )
        end if

        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LECLIMED
        !read flood plain related static data that used to be in cota_area.flp
        print*, 'Reading flood plain data'
        if (trim(adjustl(hydrodynamic_module)) == '0') then
            hdFLAG=0
        elseif (trim(adjustl(hydrodynamic_module)) == '1') then
            hdFLAG=1
            print*, 'hdflag = 1 !!!'
        end if
        hdFLAG0=sum(hdFLAG)
        !load parameters that used to be in cota_area.flp
        call check_ncrequest( nf90_inq_varid(ncid, "flood_plain_points_number", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, npfl) )
        call check_ncrequest( nf90_inq_varid(ncid, "flood_plain_points_lowest_altitude", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, zfundofl) )
        call check_ncrequest( nf90_inq_varid(ncid, "flood_plain_points_altitude", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, zfl) )
        call check_ncrequest( nf90_inq_varid(ncid, "flood_plain_points_area_flooded", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, afl) )
        

        !load parameters that used to be in face.con
        print*, 'Reading cell data (used to be in face.con)'
        call check_ncrequest( nf90_inq_varid(ncid, "face_adjacent_cell_1", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, nfacecat1) )
        call check_ncrequest( nf90_inq_varid(ncid, "face_adjacent_cell_2", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, nfacecat2) )
        call check_ncrequest( nf90_inq_varid(ncid, "face_dx", var_id) )
        call check_ncrequest( nf90_get_var(ncid, var_id, nfacedx) )
        
        !get delta information => to be generalized
        if (n_deltas.gt.0) then
            print*, 'Reading delta information'
            call check_ncrequest( nf90_inq_varid(ncid, "cell_delta_id", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, cell_delta_info) )
            call check_ncrequest( nf90_inq_varid(ncid, "delta_outlet_value", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, delta_outlet_value) )
        end if
        !get specific outlet information
        if (n_specific_outlets.gt.0) then
            print*, 'Reading specific outlet information'
            call check_ncrequest( nf90_inq_varid(ncid, "specific_outlet_face_index", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, specific_outlet_face_index) )
            call check_ncrequest( nf90_inq_varid(ncid, "specific_outlet_value", var_id) )
            call check_ncrequest( nf90_get_var(ncid, var_id, specific_outlet_value) )
        end if


        call PARCUNGE !Calculates Muskingum-Cunge parameters
        call PARCEL !Calculates parameters related to MicroBains and rivers
        call flood_TOPO !Creates a matrix with topology information that is used in the Local Inertial routing method
        call flood_TAB !Creates a table with the volume of water in the floodplain from the table of water depth vs. area obtained by DEM preprocessing (File with extension .FLP)


    end subroutine load_static_data
    
    
    
    

end module static_data_loader
