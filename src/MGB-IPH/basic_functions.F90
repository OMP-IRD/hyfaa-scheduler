module basic_functions

    implicit none

contains

  elemental function str_to_real(x) result(y)
      character(len=*), intent(in) :: x
      real :: y
      read(x,*) y
  end function



  elemental function str_to_dp(x) result(y)
      character(len=*), intent(in) :: x
      double precision :: y
      read(x,*) y
  end function



  elemental function str_to_int(x) result(y)
      character(len=*), intent(in) :: x
      integer :: y
      read(x,*) y
  end function
    
  function coma_split(input_string, ndimensions, char_len) result(string_array)
    implicit none
    integer, intent(in)::ndimensions,char_len
    character(len=*), intent(in)::input_string
    character(len=char_len), dimension(ndimensions)::string_array
    integer::nchar,ii,last_id,string_id
    string_array = ''
    nchar = len(input_string)
    if (nchar == 0) then
        return
    else if (input_string(1:1) == ',' .or. input_string(nchar:nchar) == ',') then
        print*, 'input string in coma_split cannot start or end with a coma (,)'
        call abort()
    end if
    last_id = 1
    ii = 2
    string_id = 1
    do while(ii.le.nchar)
        if (input_string(ii:ii) == ',') then
            string_array(string_id) = input_string(last_id:ii-1)
            last_id = ii+1
            string_id = string_id+1
        end if
        ii = ii+1
    end do
    string_array(string_id) = input_string(last_id:ii-1)
  end function
  
  
  function count_string_occurrence(input_string, search_string) result(occurrence)
    implicit none
    character(len=*), intent(in)::input_string,search_string
    integer::occurrence,nchar_input,nchar_search,ii
    
    nchar_input = len(input_string)
    nchar_search = len(search_string)
    occurrence = 0
    do ii=1,nchar_input
        if (ii+nchar_search-1>nchar_input) then
          return
        end if
        if (input_string(ii:ii+nchar_search-1) == search_string) then
          occurrence = occurrence + 1
        end if
    end do
  end function
  
  

  
  
  subroutine find_integer_match_id_in_array(int_in,len_array_in,array_in,int_id)
    implicit none
    integer,intent(in)::int_in,len_array_in
    integer,dimension(len_array_in)::array_in
    integer,intent(out)::int_id
    integer::ii,count_match
    count_match = 0
    do ii=1,len_array_in
      if (array_in(ii) == int_in) then
        int_id = ii
        count_match = count_match+1
      end if
    end do
    if (count_match<1) then
      print*, 'No match found in integer array for specific integer'
      call abort()
    else if (count_match>1) then
      print*, 'Mutliple matches found in integer array for specific integer'
      call abort()
    end if
  end subroutine
  

    
end module basic_functions
