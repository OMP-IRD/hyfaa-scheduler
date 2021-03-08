module hydrological_state_restart_io

	use global_variables
    use netcdf
    use netcdf_addons
    use julianday_mod

contains
    
  subroutine write_hydrological_state(hydro_state_file)
    implicit none
    character(len=*),intent(in)::hydro_state_file
    integer::ncid,nu_id,n_meshes_id,n_flat_muskingum,n_flat_muskingum_id,var_id,ii,ii0,ii1,n_meshes_id_add1,ip
    integer:: n_cells_id
    real, dimension(:), allocatable::real_vector
    double precision::decimal_julian_day_model

    
    
    
    decimal_julian_day_model = simulation_start_date+DTP*IT*1.d0/(24.d0*3600.d0)

    call check_ncrequest( nf90_create(trim(adjustl(hydro_state_file)), NF90_HDF5, ncid) )
    
    !save time at which hydrological state is saved
    call check_ncrequest( nf90_put_att(ncid, NF90_GLOBAL, 'date', datestr_from_julianday(decimal_julian_day_model)) )
    call check_ncrequest( nf90_put_att(ncid, NF90_GLOBAL, 'date_julianday', decimal_julian_day_model) )
    
    !check if mesh dimension matches
    call check_ncrequest( nf90_def_dim(ncid, "n_meshes", NC, n_meshes_id) )
    call check_ncrequest( nf90_def_dim(ncid, "n_meshes_add1", NC+1, n_meshes_id_add1) )
    call check_ncrequest( nf90_def_dim(ncid, "n_hrus", NU, nu_id) )
    n_flat_muskingum = sum(NSUBT)+NC
    call check_ncrequest( nf90_def_dim(ncid, "n_flat_muskingum", n_flat_muskingum, n_flat_muskingum_id) )
    print*,  "count_params_to_assim=",count_params_to_assim
    print*,  "params_to_assim=",list_params_in_control
    if (count_params_to_assim>0) then
        ip=1
        if (trim(list_params_in_control(ip))=="manning_coefficient" .or. &
        &trim(list_params_in_control(ip))=="river_width" .or. &
        &trim(list_params_in_control(ip))=="river_depth" .or. &
        &trim(list_params_in_control(ip))=="main_river_slope") then
            call check_ncrequest( nf90_def_dim(ncid, "n_cells", NC, n_cells_id) )
        end if
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !defining standard variables
    
    !soil moisture => W
    call check_ncrequest( nf90_def_var(ncid, "soil_moisture", NF90_FLOAT, (/ n_meshes_id, nu_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )

    !water in baseflow reservoir => VBAS
    call check_ncrequest( nf90_def_var(ncid, "water_baseflow_reservoir", NF90_FLOAT, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )

    !water in subsurface reservoir => VINT
    call check_ncrequest( nf90_def_var(ncid, "water_subsurface_reservoir", NF90_FLOAT, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )

    !water surface volume => VSUP
    call check_ncrequest( nf90_def_var(ncid, "water_surface_volume", NF90_FLOAT, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )

    !temperature => TA
    call check_ncrequest( nf90_def_var(ncid, "temperature", NF90_FLOAT, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )

    !flow entering catchment => QM2
    call check_ncrequest( nf90_def_var(ncid, "flow_entering_catchment", NF90_FLOAT, (/ n_meshes_id_add1/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )

    !water canopy interception volume => SI
    call check_ncrequest( nf90_def_var(ncid, "water_canopy_interception_volume", NF90_FLOAT, (/ n_meshes_id, nu_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )

    !number of stretches for the muskingum cunge routing => NSUBT
    call check_ncrequest( nf90_def_var(ncid, "stretchs_muskingum_cunge_routing", NF90_INT, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )

    !initial condition of the muskingum cunge routing => QRIOINI
    call check_ncrequest( nf90_def_var(ncid, "initial_condition_muskingum_cunge_routing", NF90_FLOAT, (/ n_flat_muskingum_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )

    !flow entering catchment => QCEL2
    call check_ncrequest( nf90_def_var(ncid, "runoff_generation", NF90_FLOAT, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !defining inertial variables
    
    !streamflow at catchment => Q2fl
    call check_ncrequest( nf90_def_var(ncid, "streamflow_catchment", NF90_DOUBLE, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
    
    !water depth at catchment => Hfl
    call check_ncrequest( nf90_def_var(ncid, "water_depth_catchment", NF90_DOUBLE, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
    
    !water elevation at catchment => Yfl
    call check_ncrequest( nf90_def_var(ncid, "water_elevation_catchment", NF90_DOUBLE, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
    
    !previous water storage at catchment => Vol1
    call check_ncrequest( nf90_def_var(ncid, "previous_water_storage_catchment", NF90_DOUBLE, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
    
    !water storage at catchment => Vol2
    call check_ncrequest( nf90_def_var(ncid, "water_storage_catchment", NF90_DOUBLE, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
    
    !flooded area at catchment => Area2
    call check_ncrequest( nf90_def_var(ncid, "flooded_area_catchment", NF90_DOUBLE, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
    
    !streamflow at catchment => jtab
    call check_ncrequest( nf90_def_var(ncid, "index_volume_level_catchment", NF90_INT, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
    
    !streamflow coming (or leaving) the catchment through connections => Q2face
    call check_ncrequest( nf90_def_var(ncid, "streamflow_io_connections_catchment", NF90_DOUBLE, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
    
    !streamflow at catchment => Q2viz
    call check_ncrequest( nf90_def_var(ncid, "updated_flow_interconnections", NF90_DOUBLE, (/ n_meshes_id/), var_id) )
    call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
    
    !!manning coefficient!!
    do ip=1,count_params_to_assim
        if (trim(list_params_in_control(ip))=="manning_coefficient")then
            call check_ncrequest( nf90_def_var(ncid, "manning_coefficient", NF90_FLOAT, (/ n_cells_id/), var_id) )
            call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
        end if
        if (trim(list_params_in_control(ip))=="river_depth")then
            call check_ncrequest( nf90_def_var(ncid, "river_depth", NF90_DOUBLE, (/ n_cells_id/), var_id) )
            call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
        end if
        if (trim(list_params_in_control(ip))=="river_width")then
            call check_ncrequest( nf90_def_var(ncid, "river_width", NF90_FLOAT, (/ n_cells_id/), var_id) )
            call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
        end if
        if (trim(list_params_in_control(ip))=="main_river_slope")then
            call check_ncrequest( nf90_def_var(ncid, "main_river_slope", NF90_FLOAT, (/ n_cells_id/), var_id) )
            call check_ncrequest( nf90_def_var_deflate(ncid, var_id, 1, 1, 4) )
        end if
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !ending definition mode
    call check_ncrequest( nf90_enddef(ncid) )
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !writing standard variables
    
    !soil moisture => W
    call check_ncrequest( nf90_inq_varid(ncid, "soil_moisture", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, W) )
    
    !water in baseflow reservoir => VBAS
    call check_ncrequest( nf90_inq_varid(ncid, "water_baseflow_reservoir", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, VBAS) )

    !water in subsurface reservoir => VINT
    call check_ncrequest( nf90_inq_varid(ncid, "water_subsurface_reservoir", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, VINT) )

    !water surface volume => VSUP
    call check_ncrequest( nf90_inq_varid(ncid, "water_surface_volume", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, VSUP) )

    !temperature => TA
    call check_ncrequest( nf90_inq_varid(ncid, "temperature", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, TA) )

    !flow entering catchment => QM2
    call check_ncrequest( nf90_inq_varid(ncid, "flow_entering_catchment", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, QM2) )

    !water canopy interception volume => SI
    call check_ncrequest( nf90_inq_varid(ncid, "water_canopy_interception_volume", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, SI) )
    
    !number of stretches for the muskingum cunge routing => NSUBT
    call check_ncrequest( nf90_inq_varid(ncid, "stretchs_muskingum_cunge_routing", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, NSUBT) )
    
    !initial condition of the muskingum cunge routing => QRIOINI
    allocate(real_vector(n_flat_muskingum))
    ii = 0
    do ii0=1,NC
      do ii1=1,NSUBT(ii0)+1
        ii = ii +1
        real_vector(ii) = QRIOINI(ii0,ii1)
      end do
    end do
    call check_ncrequest( nf90_inq_varid(ncid, "initial_condition_muskingum_cunge_routing", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, real_vector) )
    deallocate(real_vector)
    
    !flow entering catchment => QCEL2
    call check_ncrequest( nf90_inq_varid(ncid, "runoff_generation", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, QCEL2) )
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !defining inertial variables
    
    !streamflow at catchment => Q2fl
    call check_ncrequest( nf90_inq_varid(ncid, "streamflow_catchment", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, Q2fl) )
    
    !water depth at catchment => Hfl
    call check_ncrequest( nf90_inq_varid(ncid, "water_depth_catchment", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, Hfl) )
    
    !water elevation at catchment => Yfl
    call check_ncrequest( nf90_inq_varid(ncid, "water_elevation_catchment", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, Yfl) )
    
    !previous water storage at catchment => Vol1
    call check_ncrequest( nf90_inq_varid(ncid, "previous_water_storage_catchment", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, Vol1) )
    
    !water storage at catchment => Vol2
    call check_ncrequest( nf90_inq_varid(ncid, "water_storage_catchment", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, Vol2) )
    
    !flooded area at catchment => Area2
    call check_ncrequest( nf90_inq_varid(ncid, "flooded_area_catchment", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, Area2) )
    
    !streamflow at catchment => jtab
    call check_ncrequest( nf90_inq_varid(ncid, "index_volume_level_catchment", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, jtab) )
    
    !streamflow coming (or leaving) the catchment through connections => Q2face
    call check_ncrequest( nf90_inq_varid(ncid, "streamflow_io_connections_catchment", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, Q2face) )
    
    !streamflow at catchment => Q2viz
    call check_ncrequest( nf90_inq_varid(ncid, "updated_flow_interconnections", var_id) )
    call check_ncrequest( nf90_put_var(ncid, var_id, Q2viz) )
    
    !!!!!!! addition for parameters!!!!!!
   
   !!manning coefficient!!
   do ip=1,count_params_to_assim
        if (list_params_in_control(ip)=="manning_coefficient")then
            
            call check_ncrequest( nf90_inq_varid(ncid, "manning_coefficient", var_id) )
            call check_ncrequest( nf90_put_var(ncid, var_id, RUGMAN) )
        end if
        if (list_params_in_control(ip)=="river_depth")then

            call check_ncrequest( nf90_inq_varid(ncid, "river_depth", var_id) )
            call check_ncrequest( nf90_put_var(ncid, var_id, HRIO) )
        end if
        if (list_params_in_control(ip)=="river_width")then

            call check_ncrequest( nf90_inq_varid(ncid, "river_width", var_id) )
            call check_ncrequest( nf90_put_var(ncid, var_id, BRIO) )
        end if
        if (list_params_in_control(ip)=="main_river_slope")then

            call check_ncrequest( nf90_inq_varid(ncid, "main_river_slope", var_id) )
            call check_ncrequest( nf90_put_var(ncid, var_id, SRIO) )
        end if
   end do
    
    call check_ncrequest( nf90_close(ncid) )

  end subroutine write_hydrological_state
  
  
  
  
  
  
  subroutine read_hydrological_state(hydro_state_file)
    implicit none
    character(len=*),intent(in)::hydro_state_file
    integer::ncid,var_id,ii,ii0,ii1,cc
    real, dimension(:), allocatable::real_vector
    logical :: param_in_hydro
    
    print*,  "Entering read routine"
    call check_ncrequest( nf90_open(trim(adjustl(hydro_state_file)), NF90_NOwrite, ncid) )
    !check if mesh dimension matches
    call check_dimension_match(ncid, "n_meshes", NC)
    call check_dimension_match(ncid, "n_meshes_add1", NC+1)
    call check_dimension_match(ncid, "n_hrus", NU)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !reading standard variables
    

    !soil moisture => W
    call check_ncrequest( nf90_inq_varid(ncid, "soil_moisture", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, W) )
    

    !water in baseflow reservoir => VBAS
    call check_ncrequest( nf90_inq_varid(ncid, "water_baseflow_reservoir", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, VBAS) )
    

    !water in subsurface reservoir => VINT
    call check_ncrequest( nf90_inq_varid(ncid, "water_subsurface_reservoir", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, VINT) )
    

    !water surface volume => VSUP
    call check_ncrequest( nf90_inq_varid(ncid, "water_surface_volume", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, VSUP) )
    

    !temperature => TA
    call check_ncrequest( nf90_inq_varid(ncid, "temperature", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, TA) )
    

    !flow entering catchment => QM2
    call check_ncrequest( nf90_inq_varid(ncid, "flow_entering_catchment", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, QM2) )
    

    !water canopy interception volume => SI
    call check_ncrequest( nf90_inq_varid(ncid, "water_canopy_interception_volume", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, SI) )
    

    !initial condition of the muskingum cunge routing => QRIOINI
    !check number of stretchs for muskingum cunge routing
    call check_ncrequest( nf90_inq_varid(ncid, "stretchs_muskingum_cunge_routing", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, nsubt) )
    allocate(real_vector(sum(NSUBT)+NC))
    call check_ncrequest( nf90_inq_varid(ncid, "initial_condition_muskingum_cunge_routing", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, real_vector) )
    QRIOINI=0.
    ii = 0
    do ii0=1,NC
      do ii1=1,NSUBT(ii0)+1
        ii = ii +1
        QRIOINI(ii0,ii1) = real_vector(ii)
      end do
    end do
    deallocate(real_vector)
    

    !flow entering catchment => QCEL2
    call check_ncrequest( nf90_inq_varid(ncid, "runoff_generation", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, QCEL2) )
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !reading inertial variables

    !streamflow at catchment => Q2fl
    call check_ncrequest( nf90_inq_varid(ncid, "streamflow_catchment", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, Q2fl) )
    

    !water depth at catchment => Hfl
    call check_ncrequest( nf90_inq_varid(ncid, "water_depth_catchment", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, Hfl) )
    

    !water elevation at catchment => Yfl
    call check_ncrequest( nf90_inq_varid(ncid, "water_elevation_catchment", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, Yfl) )
    

    !previous water storage at catchment => Vol1
    call check_ncrequest( nf90_inq_varid(ncid, "previous_water_storage_catchment", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, Vol1) )
    

    !water storage at catchment => Vol2
    call check_ncrequest( nf90_inq_varid(ncid, "water_storage_catchment", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, Vol2) )
    

    !flooded area at catchment => Area2
    call check_ncrequest( nf90_inq_varid(ncid, "flooded_area_catchment", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, Area2) )
    

    !streamflow at catchment => jtab
    call check_ncrequest( nf90_inq_varid(ncid, "index_volume_level_catchment", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, jtab) )
    

    !streamflow coming (or leaving) the catchment through connections => Q2face
    call check_ncrequest( nf90_inq_varid(ncid, "streamflow_io_connections_catchment", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, Q2face) )
    

    !streamflow at catchment => Q2viz
    call check_ncrequest( nf90_inq_varid(ncid, "updated_flow_interconnections", var_id) )
    call check_ncrequest( nf90_get_var(ncid, var_id, Q2viz) )
    print*,  "State variables read, counting control parameters"
!!!!!!! addition for parameters!!!!!!
    count_params_to_assim=0
    call check_ncrequest_opt( nf90_inq_varid(ncid, "manning_coefficient", var_id),param_in_hydro )
    if (param_in_hydro) then
        count_params_to_assim=count_params_to_assim+1
    end if
    print*,  "counting params..." ,count_params_to_assim
    call check_ncrequest_opt( nf90_inq_varid(ncid, "river_depth", var_id),param_in_hydro )
    if (param_in_hydro) then
        count_params_to_assim=count_params_to_assim+1
    end if
    print*,  "counting params..." ,count_params_to_assim
    call check_ncrequest_opt( nf90_inq_varid(ncid, "river_width", var_id),param_in_hydro )
    if (param_in_hydro) then
        count_params_to_assim=count_params_to_assim+1
    end if
    print*,  "counting params..." ,count_params_to_assim
    call check_ncrequest_opt( nf90_inq_varid(ncid, "main_river_slope", var_id),param_in_hydro )
    if (param_in_hydro) then
        count_params_to_assim=count_params_to_assim+1
    end if
    print*,  "counting params..." ,count_params_to_assim

    print*,  "number of assimilated parameters=",count_params_to_assim
    allocate(list_params_in_control(count_params_to_assim))    
    
    cc=0
    call check_ncrequest_opt( nf90_inq_varid(ncid, "manning_coefficient", var_id),param_in_hydro )
    if (param_in_hydro) then
        call check_dimension_match(ncid, "n_cells", NC)
        call check_ncrequest( nf90_get_var(ncid, var_id, RUGMAN) )
        list_params_in_control(cc+1)="manning_coefficient"
        cc=cc+1
        print*,  "RUGMAN=",SUM(RUGMAN)
    end if
    call check_ncrequest_opt( nf90_inq_varid(ncid, "river_width", var_id),param_in_hydro )
    if (param_in_hydro) then
        call check_ncrequest( nf90_get_var(ncid, var_id, BRIO) )
        list_params_in_control(cc+1)="river_width"
        cc=cc+1
    end if
    call check_ncrequest_opt( nf90_inq_varid(ncid, "river_depth", var_id),param_in_hydro )
    if (param_in_hydro) then
        call check_ncrequest( nf90_get_var(ncid, var_id, HRIO) )
        list_params_in_control(cc+1)="river_depth"
        cc=cc+1
        print*,  "river depth=",SUM(HRIO)
    end if
    call check_ncrequest_opt( nf90_inq_varid(ncid, "main_river_slope", var_id),param_in_hydro )
    if (param_in_hydro) then
        call check_ncrequest( nf90_get_var(ncid, var_id, SRIO) )
        list_params_in_control(cc+1)="main_river_slope"
        cc=cc+1
    end if
    print*,  "list of assimilated parameters=",list_params_in_control
    call check_ncrequest( nf90_close(ncid) )


    if (count_params_to_assim>0) then
      call PARCUNGE !Calculates Muskingum-Cunge parameters
      call PARCEL !Calculates parameters related to MicroBains and rivers
      call flood_TOPO !Creates a matrix with topology information that is used in the Local Inertial routing method
      call flood_TAB !Creates a table with the volume of water in the floodplain from the table of water depth vs. area obtained by DEM preprocessing (File with extension .FLP)
       print*,  "Hydrodynamics parameters recalculated"
    end if
  end subroutine read_hydrological_state
  
  
end module hydrological_state_restart_io
