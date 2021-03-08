module julianday_mod

  use datetime_module, only:datetime,timedelta
  implicit none
  
  
contains

    function julianday(year,month,day,hour_in,minute_in,second_in,millisecond_in)
      implicit none
      integer, intent(in)::year,month,day
      integer, optional, intent(in)::hour_in,minute_in,second_in,millisecond_in
      integer::hour,minute,second,millisecond
      double precision::julianday
      type(timedelta) :: timedelta_obj
      
      if (present(hour_in)) then
        hour = hour_in
      else
        hour = 0
      end if
      if (present(minute_in)) then
        minute = minute_in
      else
        minute = 0
      end if
      if (present(second_in)) then
        second = second_in
      else
        second = 0
      end if
      if (present(millisecond_in)) then
        millisecond = millisecond_in
      else
        millisecond = 0
      end if
      timedelta_obj=datetime(year,month,day,hour,minute,second,millisecond)-datetime(1950,1,1,0,0,0)
      julianday = timedelta_obj%total_seconds()*1.d0/(24.d0*3600.d0)

    end function julianday
    
    
    function datetime_from_julianday(decimal_julian_day) result(datetime_obj)
      implicit none
      double precision,intent(in)::decimal_julian_day
      double precision::decimal_second_residue
      integer::days,seconds,milliseconds
      type(datetime) :: datetime_obj
      days = floor(decimal_julian_day)
      decimal_second_residue = (decimal_julian_day-days*1.d0)*24.d0*3600.d0
      seconds = floor(decimal_second_residue)
      milliseconds = nint((decimal_second_residue-seconds)*1.d3)
      datetime_obj=datetime(1950,1,1,0,0,0)+timedelta(days=days,seconds=seconds,milliseconds=milliseconds)
    end function datetime_from_julianday
    
    
    subroutine yearmonthday_from_julianday(decimal_julian_day,year,month,day)
      implicit none
      double precision,intent(in)::decimal_julian_day
      integer, intent(out)::year,month,day
      type(datetime) :: datetime_obj
      datetime_obj = datetime_from_julianday(decimal_julian_day)
      year = datetime_obj%getYear()
      month = datetime_obj%getMonth()
      day = datetime_obj%getDay()
    end subroutine yearmonthday_from_julianday
    
    
    function datestr_from_julianday(decimal_julian_day, separator_str_in) result(datestr)
      implicit none
      double precision,intent(in)::decimal_julian_day
      character(len=1),optional,intent(in)::separator_str_in
      character(len=1)::separator_str
      type(datetime)::datetime_obj
      character(len=23)::datestr
      separator_str = 'T'
      datetime_obj = datetime_from_julianday(decimal_julian_day)
      write(datestr,'(I0.3)') datetime_obj%getMillisecond()
      datestr = trim(adjustl(datetime_obj%strftime("%Y-%m-%d"//separator_str//"%H:%M:%S.")))//trim(adjustl(datestr))
    end function datestr_from_julianday
    

end module julianday_mod
