module netcdf_addons
use basic_functions
use netcdf

contains

  
  
  
  subroutine check_ncrequest(status)
    implicit none
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check_ncrequest
  
  
  function nc_dimension(ncid, dim_name)
    implicit none
    integer, intent(in)::ncid
    character(len=*),intent(in)::dim_name
    integer::nc_dimension,dim_id
    call check_ncrequest( nf90_inq_dimid(ncid, trim(adjustl(dim_name)), dim_id) )
    call check_ncrequest( nf90_inquire_dimension(ncid, dim_id, len = nc_dimension) )
  end function nc_dimension
  
  
  subroutine check_dimension_match(ncid, dim_name, dim_value_comp)
    implicit none
    integer, intent(in)::ncid,dim_value_comp
    character(len=*),intent(in)::dim_name
    character(len=50)::print_str1,print_str2
    integer::dim_value
    dim_value = nc_dimension(ncid,dim_name)
    if (dim_value .ne. dim_value_comp) then
      write(print_str1,*) dim_value
      write(print_str2,*) dim_value_comp
      print*, 'number of mesh elements ('//trim(adjustl(print_str1))//') '// &
        & 'different from declared number NC='//trim(adjustl(print_str2))
        call abort()
    end if
  end subroutine check_dimension_match

  subroutine check_ncrequest_opt(status,paramvar)
    implicit none
    integer, intent ( in) :: status
    logical, intent(out) :: paramvar
    
    if(status /= nf90_noerr) then 
      paramvar=.false.
    else
      paramvar=.true.
    end if
  end subroutine check_ncrequest_opt

  function nc_ndims(ncid, varid)
    implicit none
    integer, intent(in)::ncid,varid
    integer::nc_ndims,trash_int,trash_int2
    character(50)::trash_str
    integer, dimension(10)::dimids
    call check_ncrequest( nf90_inquire_variable(ncid, varid, trash_str, trash_int, nc_ndims, dimids, trash_int2) )
  end function nc_ndims
  
end module netcdf_addons


