
subroutine read_forcing

  !RAIN input must be in mm and match model time step. For 1 day time step, input must be in mm/day for example.
  !routine will look for available forcing data prior to current model date, with a date difference, in days, < forcing_data_dates_dt_max

  use global_variables
  use julianday_mod
  use netcdf
  use netcdf_addons

  implicit none
  double precision::decimal_julian_day_model,dt_loc
  integer::closest_date_id,mm,iyyy,id,specific_climatology_years_id,it_loc
  
  call caldat(jdia,mm,id,iyyy)
  
  decimal_julian_day_model = simulation_start_date+DTP*(ITCHUVA-1)*1.d0/(24.d0*3600.d0)
  
  !get best forcing date match for before or equal model date within forcing_data_dates_dt_max time range
  closest_date_id = -1
  it_loc = 1
  do while (.true.)
      dt_loc = decimal_julian_day_model-forcing_file_julianday_dates(it_loc)
      if (dt_loc .lt. 0.d0) then
          exit
      end if
      if (dt_loc .le. forcing_data_dates_dt_max) then
          if (closest_date_id==-1) then
              closest_date_id = it_loc
          else
              if (dt_loc .lt. (decimal_julian_day_model-forcing_file_julianday_dates(it_loc))) then
                  closest_date_id = it_loc
              end if
          end if
      end if
      it_loc = it_loc + 1
      if (it_loc .gt. n_forcing_file_dates) then
          exit
      end if
  end do
  if (closest_date_id==-1) then
      print*, 'could not find matching forcing date within allowed time range'
      call abort()
  end if
  
  write(temp_str,*) ITCHUVA
  write(temp_str2,*) closest_date_id
  temp_str = "iteration "//trim(adjustl(temp_str))//' , date '// &
    & trim(adjustl(datestr_from_julianday(decimal_julian_day_model)))//' :'//&
    & ' loading precipitations from date '// &
    & trim(adjustl(datestr_from_julianday(forcing_file_julianday_dates(closest_date_id))))// &
    & ' (index='//trim(adjustl(temp_str2))//')'
  print*, trim(adjustl(temp_str))
  
  
  !load precipitation data from netCDF rain variable (forcing_rain_id obtained in 1main.f90)
  call check_ncrequest( nf90_get_var(FILPLU, forcing_rain_id, p, start = (/closest_date_id,1/), count = (/1,NC/)) )
  
  !load other forcing variables from dynamic forcing or climatology
  
    !find out if this year is a specific climatology year
    specific_climatology_years_id = 1
    if (n_climatology.gt.1) then
        do it_loc=2,n_climatology
            if (iyyy==specific_climatology_years(it_loc)) then
                if (specific_climatology_years_id.gt.1) then
                    print*, 'multiple match for climatology year, this means that the specific_climatology_years vector has duplicate years'
                    call abort()
                end if
                specific_climatology_years_id = it_loc
            end if
        end do
    end if
    
    
    !temperature
    if (forcing_use_dynamic_temperature==1) then
        call check_ncrequest( nf90_get_var(FILPLU, forcing_varid_dynamic_temperature, ta, start = (/closest_date_id,1/), count = (/1,NC/)) )
    else
        do ic=1,nc
            ta(ic)=tamm(specific_climatology_years_id,ic,mm)
        end do
    end if
    
    !relative humidity
    if (forcing_use_dynamic_relative_humidity==1) then
        call check_ncrequest( nf90_get_var(FILPLU, forcing_varid_dynamic_relative_humidity, ta, start = (/closest_date_id,1/), count = (/1,NC/)) )
    else
        do ic=1,nc
            ur(ic)=urmm(specific_climatology_years_id,ic,mm)
        end do
    end if
    
    !wind speed
    if (forcing_use_dynamic_wind_speed==1) then
        call check_ncrequest( nf90_get_var(FILPLU, forcing_varid_dynamic_wind_speed, ta, start = (/closest_date_id,1/), count = (/1,NC/)) )
    else
        do ic=1,nc
            vv(ic)=vvmm(specific_climatology_years_id,ic,mm)
        end do
    end if
    
    !sunhine
    if (forcing_use_dynamic_sunshine==1) then
        call check_ncrequest( nf90_get_var(FILPLU, forcing_varid_dynamic_sunshine, ta, start = (/closest_date_id,1/), count = (/1,NC/)) )
    else
        do ic=1,nc
            sol(ic)=solmm(specific_climatology_years_id,ic,mm)
        end do
    end if
    
    !pressure
    if (forcing_use_dynamic_pressure==1) then
        call check_ncrequest( nf90_get_var(FILPLU, forcing_varid_dynamic_pressure, ta, start = (/closest_date_id,1/), count = (/1,NC/)) )
    else
        do ic=1,nc
            pa(ic)=pamm(specific_climatology_years_id,ic,mm)
        end do
    end if

end subroutine read_forcing

