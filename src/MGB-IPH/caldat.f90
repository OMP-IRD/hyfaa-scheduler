subroutine caldat(julian,mm,id,iyyy)
  use julianday_mod
  implicit none
  integer,intent(in)::julian
  integer,intent(out)::mm,id,iyyy
  call yearmonthday_from_julianday(julian*1.d0,iyyy,mm,id)
end subroutine
      
