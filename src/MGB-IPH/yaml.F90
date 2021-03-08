! -----------------------------------------------------------------------------
! This file is part of Fortran-YAML: a lightweight YAML parser written in
! object-oriented Fortran.
!
! Official repository: https://github.com/BoldingBruggeman/fortran-yaml
!
! Copyright 2013-2016 Bolding & Bruggeman ApS.
!
! This is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation (https://www.gnu.org/licenses/gpl.html). It is distributed in the
! hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! A copy of the license is provided in the COPYING file.
! -----------------------------------------------------------------------------

module yaml
   use basic_functions

   implicit none
   
   private

   public string_length,type_node,type_scalar,type_null,type_error,real_kind
   public type_dictionary,type_key_value_pair
   public type_list,type_list_item
   public parse,error_length,yaml_read,yaml_dump
   public dictget_strvalue,dictget_realvalue,dictget_dpvalue,dictget_intvalue
   public dictget_strarray,dictget_realarray,dictget_dparray,dictget_intarray
   public get_node_stringvalue,get_valuenode_from_node_and_key,str_to_int,str_to_dp,str_to_real,dictget_listlength
   
   integer,parameter :: string_length = 1024
   integer,parameter :: real_kind = 8
   integer,parameter :: line_length  = 2048
   integer,parameter :: error_length = 2048
   
   type type_file
      integer                 :: unit   = -1
      character(line_length)  :: line   = ''
      integer                 :: indent = 0
      logical                 :: eof    = .false.
      integer                 :: iline  = 0
      character(error_length) :: error_message = ''
      logical                 :: has_error     = .false.
   contains
      procedure :: next_line
      procedure :: set_error
   end type

   type,abstract :: type_node
      character(len=string_length) :: path = ''
   contains
      procedure (node_dump),deferred :: dump
      procedure                      :: set_path => node_set_path
      procedure                      :: finalize => node_finalize
   end type

   abstract interface
      subroutine node_dump(self,unit,indent)
         import type_node
         class (type_node),intent(in) :: self
         integer,intent(in) :: unit,indent
      end subroutine
   end interface

   type,extends(type_node) :: type_scalar
      character(len=string_length) :: string = ''
   contains
      procedure :: dump       => value_dump
      procedure :: to_logical => scalar_to_logical
      procedure :: to_integer => scalar_to_integer
      procedure :: to_real    => scalar_to_real
   end type

   type,extends(type_node) :: type_null
   contains
      procedure :: dump => null_dump
   end type

   type type_key_value_pair
      character(len=string_length)       :: key   = ''
      class (type_node),         pointer :: value => null()
      logical                            :: accessed = .false.
      type (type_key_value_pair),pointer :: next  => null()
   end type

   type,extends(type_node) :: type_dictionary
      type (type_key_value_pair),pointer :: first => null()
   contains
      procedure :: get            => dictionary_get
      procedure :: get_scalar     => dictionary_get_scalar
      procedure :: get_dictionary => dictionary_get_dictionary
      procedure :: get_list       => dictionary_get_list
      procedure :: get_string     => dictionary_get_string
      procedure :: get_logical    => dictionary_get_logical
      procedure :: get_integer    => dictionary_get_integer
      procedure :: get_real       => dictionary_get_real
      procedure :: set            => dictionary_set
      procedure :: set_string     => dictionary_set_string
      procedure :: dump           => dictionary_dump
      procedure :: flatten        => dictionary_flatten
      procedure :: reset_accessed => dictionary_reset_accessed
      procedure :: set_path       => dictionary_set_path
      procedure :: finalize       => dictionary_finalize
   end type

   type type_list_item
      class (type_node),    pointer :: node => null()
      type (type_list_item),pointer :: next => null()
   end type

   type,extends(type_node) :: type_list
      type (type_list_item),pointer :: first => null()
   contains
      procedure :: append   => list_append
      procedure :: dump     => list_dump
      procedure :: set_path => list_set_path
   end type

   type type_error
      character(len=string_length) :: message
   end type

contains

    subroutine yaml_dump(root_in, path, unit_open_in)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::path
        integer, optional, intent(in)::unit_open_in
        integer::unit_open
        if (present(unit_open_in)) then
            unit_open = unit_open_in
        else
            unit_open = 17567
        end if
        open(unit=unit_open,file=path,action='write')
        call root_in%dump(unit=unit_open,indent=0)
        close(unit_open)
    end subroutine yaml_dump

    function yaml_read(path, unit_open_in) result(output_node)
        character(len=*), intent(in)::path
        integer, optional, intent(in)::unit_open_in
        class (type_node), pointer::output_node
        integer::unit_open
        character(len=error_length)::error
        
        if (trim(adjustl(path))=='') then
          print*, 'ERROR: path to YAML file not provided.'
          call abort()
        end if
        if (present(unit_open_in)) then
            unit_open = unit_open_in
        else
            unit_open = 17567
        end if
        output_node => parse(path,unit=unit_open,error=error)
        if (error/='') then
          print*, 'PARSE ERROR: '//trim(error)
          call abort()
        end if
    end function yaml_read

    function dictget_realvalue(root_in, key) result(output_value)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        real::output_value
        output_value = str_to_real(dictget_strvalue(root_in, key))
    end function dictget_realvalue
    
    function dictget_dpvalue(root_in, key) result(output_value)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        double precision::output_value
        output_value = str_to_dp(dictget_strvalue(root_in, key))
    end function dictget_dpvalue
    
    function dictget_intvalue(root_in, key) result(output_value)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        integer::output_value
        output_value = str_to_int(dictget_strvalue(root_in, key))
    end function dictget_intvalue

    function dictget_strvalue(root_in, key) result(output_value)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        character(string_length)::output_value
        class (type_node), pointer::temp_obj
        
        if (value_is_none(root_in, key)) then
            output_value = 'none'
            return
        end if
        
        ! get value node
        temp_obj => get_valuenode_from_node_and_key(root_in, key)
        
        ! get the string value from the node (expecting node to be of type scalar)
        output_value = get_node_stringvalue(temp_obj)
        
    end function dictget_strvalue


    
    function dictget_realarray(root_in, key, output_len) result(output_array)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        integer, intent(in)::output_len
        real, dimension(output_len)::output_array
        output_array = str_to_real(dictget_strarray(root_in, key, output_len, 50))
    end function dictget_realarray
    
    
    function dictget_dparray(root_in, key, output_len) result(output_array)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        integer, intent(in)::output_len
        double precision, dimension(output_len)::output_array
        output_array = str_to_dp(dictget_strarray(root_in, key, output_len, 50))
    end function dictget_dparray
    
    
    function dictget_intarray(root_in, key, output_len) result(output_array)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        integer, intent(in)::output_len
        integer, dimension(output_len)::output_array
        output_array = str_to_int(dictget_strarray(root_in, key, output_len, 50))
    end function dictget_intarray
    
    
    function dictget_strarray(root_in, key, output_len, output_str_len) result(output_array)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        integer, intent(in)::output_len, output_str_len
        character(output_str_len), dimension(output_len)::output_array
        class (type_node), pointer::temp_obj
        type (type_list_item), pointer :: temp_item
        logical::go_on
        integer::ii
        
        ! get value node
        temp_obj => get_valuenode_from_node_and_key(root_in, key)
        
        select type (temp_obj)
           class is (type_list)
              if (.not.check_associated_list_item(temp_obj%first)) then
                  print*, 'list is empty'
                  call abort()
              end if
              
              temp_item => temp_obj%first
              go_on = .TRUE.
              ii = 1
              do while (go_on)
              
                  output_array(ii) = trim(adjustl(get_node_stringvalue(temp_item%node)))
                  
                  if (check_associated_list_item(temp_item%next)) then
                      temp_item => temp_item%next
                      ii = ii+1
                  else
                      go_on = .FALSE.
                  end if
              end do

           class default
              print*, 'temp_obj not type_list'
              call abort()
        end select

    end function dictget_strarray
    
    
    function value_is_none(root_in, key)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        class (type_node), pointer::temp_node
        logical::value_is_none
        
        value_is_none = .FALSE.
        temp_node => get_valuenode_from_node_and_key(root_in, key)
        if (.not.associated(temp_node)) then
            value_is_none = .TRUE.
        end if
        select type (temp_node)
            class is (type_scalar)
                value_is_none = len(trim(adjustl(temp_node%string))) == 0
                if (.not.value_is_none) then
                  if ((trim(adjustl(temp_node%string)) == 'None').or.(trim(adjustl(temp_node%string)) == 'none')) then
                     value_is_none = .TRUE.
                  end if
                end if
            class is (type_null)
                value_is_none = .TRUE.
        end select
    end function value_is_none
    
    
    
    function dictget_listlength(root_in, key) result(output_len)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        integer::output_len
        class (type_node), pointer::temp_obj
        type (type_list_item), pointer :: temp_item
        logical::go_on
        
        ! get value node
        temp_obj => get_valuenode_from_node_and_key(root_in, key)
        
        output_len = 0
        select type (temp_obj)
           class is (type_list)
              if (check_associated_list_item(temp_obj%first)) then
                  temp_item => temp_obj%first
                  output_len = output_len+1
                  go_on = .TRUE.
                  do while (go_on)
                      if (check_associated_list_item(temp_item%next)) then
                          temp_item => temp_item%next
                          output_len = output_len+1
                      else
                          go_on = .FALSE.
                      end if
                  end do
              end if
        end select
        
    end function dictget_listlength
    
    
    function check_associated_list_item(list_item) result(valid)
        type (type_list_item), pointer :: list_item
        class (type_node), pointer::temp_obj
        logical::valid
        
        valid = .FALSE.
        if (.not.associated(list_item)) then
            return
        end if
        if (.not.associated(list_item%node)) then
            return
        end if
        temp_obj => list_item%node
        select type (temp_obj)
            class is (type_scalar)
                valid = .TRUE.
        end select
        
    end function check_associated_list_item

    
    function get_valuenode_from_node_and_key(root_in, key) result(temp_obj)
        class (type_node), pointer, intent(in)::root_in
        character(len=*), intent(in)::key
        class (type_node), pointer::temp_obj
        
        !check if input dict pointer is associated
        if (.not.associated(root_in)) then
            print*, 'input root_in not associated'
            call abort()
        end if

        ! check if input dict is a dict
        ! make temp_obj type_node pointer point towards the right (key,node)
        select type (root_in)
           class is (type_dictionary)
              temp_obj => root_in%get(trim(adjustl(key)))
           class default
              print*, 'input root_in not a dictionnary'
              call abort()
        end select
        
    end function get_valuenode_from_node_and_key


    function get_node_stringvalue(temp_obj) result(output_str)
        class (type_node), pointer, intent(in)::temp_obj
        character(string_length)::output_str
        
        !check if input type_node pointer is associated
        if (.not.associated(temp_obj)) then
            print*, 'input temp_obj pointer not associated'
            call abort()
        end if
        
        ! check if input type_node pointer is of type scalar
        ! load string from type_scalar
        select type (temp_obj)
           class is (type_scalar)
              output_str = temp_obj%string
           class default
              print*, 'temp_obj not type_scalar (i.e. not a string, but a more complex type)'
              call abort()
        end select

    end function get_node_stringvalue


























!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   function parse(path,unit,error) result(root)
      integer,                intent(in)  :: unit
      character(len=*),       intent(in)  :: path
      character(error_length),intent(out) :: error
      class (type_node),pointer           :: root

      type (type_file) :: file
      logical          :: already_open

      nullify(root)
      error = ''

      inquire(unit=unit, opened=already_open)
      if (.not.already_open) open(unit=unit,file=path,status='old',action='read',err=90)
      file%unit = unit
      file%eof = .false.
      call file%next_line()
      if (.not.file%has_error) root => read_value(file)
      if (.not.already_open) close(file%unit)
      if (file%has_error) then
         write (error,'(a,a,i0,a,a)') trim(path),', line ',file%iline,': ',trim(file%error_message)
      elseif (.not.file%eof) then
         if (associated(root)) then
            select type (root)
               class is (type_dictionary)
                  write (error,'(a,a,i0,a)') trim(path),', line ',file%iline,': unexpected decrease in indentation.'
               class is (type_scalar)
                  write (error,'(a,a,i0,a)') trim(path),', line ',file%iline,': expected end of file after reading &
                                             &one scalar value.'
               class default
                  write (error,'(a,a,i0,a)') trim(path),', line ',file%iline,': expected end of file.'
            end select
         else
            write (error,'(a,a,i0,a)') trim(path),', line ',file%iline,': expected end of file.'
         end if
      end if

      if (associated(root)) call root%set_path('')

      return

90    error = 'Unable to open '//trim(path)//' for reading.'
   end function

   subroutine next_line(file)
      class (type_file),intent(inout) :: file
      integer                         :: i
      logical                         :: done

      done = .false.
      do while (.not.done)
         ! Read entire line
         read (file%unit,'(A)',end=91) file%line
         file%iline = file%iline + 1

         ! Determine indentation and strip this.
         file%indent = len(file%line)
         do i=1,len(file%line)
            if (file%line(i:i)==achar(9)) then
               ! Found tabs in indentation: not allowed.
               call file%set_error('tab in indentation is not allowed.')
               return
            elseif (file%line(i:i)/=' ') then
               ! Found non-space: indentation ends here.
               file%indent = i-1
               exit
            end if
         end do
         file%line = file%line(file%indent+1:)

         ! If the line starts with comment character; move to next.
         if (file%line(1:1)=='#') cycle

         ! Search for whitespace delimited comment within the string; remove if found.
         do i=1,len_trim(file%line)-1
            if (is_whitespace(file%line(i:i)).and.file%line(i+1:i+1)=='#') then
               file%line = file%line(:i-1)
               exit
            end if
         end do

         ! Strip trailing whitespace
         do i=len(file%line),1,-1
            if (.not.is_whitespace(file%line(i:i))) then
               ! We found a non-whitespace character. Strip trailing whitespace and report we have a valid line.
               file%line = file%line(:i)
               done = .true.
               exit
            end if
         end do
      end do

      ! Check for unsupported YAML features.
      do i=1,len_trim(file%line)
         if (file%line(i:i)=='['.or.file%line(i:i)==']'.or.file%line(i:i)=='{'.or.file%line(i:i)=='}') then
            call file%set_error('flow mappings and sequences using []{} are not supported.')
            return
         end if
         if (file%line(i:i)=='"'.or.file%line(i:i)=='''') then
            call file%set_error('single- and double-quoted strings are not supported.')
            return
         end if
      end do

      return

91    file%indent = 0
      file%eof = .true.
   end subroutine

   recursive function read_value(file) result(node)
      class (type_file),intent(inout) :: file
      class (type_node),pointer       :: node

      integer                    :: icolon,icolon_stop,firstindent
      type (type_key_value_pair) :: pair
      class (type_node), pointer :: list_item

      nullify(node)
      if (file%eof) return

      if (file%line(1:2)=='- ') then
         allocate(type_list::node)
         firstindent = file%indent
         do
            file%line = file%line(3:)
            file%indent = file%indent + 2
            list_item => read_value(file)
            if (file%has_error) return
            select type (node)
               class is (type_list)
                  call node%append(list_item)
            end select

            ! Check indentation of next line.
            if (file%indent>firstindent) then
               call file%set_error('unexpected increase in indentation following list item.')
               return
            elseif (file%eof .or. file%indent<firstindent) then
               ! End-of-file or decrease in indentation signifies that the list has ended.
               return
            end if
         end do
      end if

      ! Find the first colon (if any)
      call find_mapping_character(file%line,icolon,icolon_stop)

      if (icolon==-1) then
         ! No colon found: item is a value
         allocate(type_scalar::node)
         select type (node)
            class is (type_scalar)
               node%string = trim(file%line)
         end select
         call file%next_line()
      else
         ! Colon found: item starts a mapping
         allocate(type_dictionary::node)
         firstindent = file%indent
         do
            pair = read_key_value_pair(file,icolon,icolon_stop)
            if (file%has_error) return
            select type (node)
               class is (type_dictionary)
                  call node%set(pair%key,pair%value)
            end select

            ! Check indentation of next line.
            if (file%indent>firstindent) then
               call file%set_error('unexpected increase in indentation following key-value pair "'//trim(pair%key)//'".')
               return
            elseif (file%eof .or. file%indent<firstindent) then
               ! End-of-file or decrease in indentation signifies that the mapping has ended.
               exit
            end if

            ! We are expecting a new key-value pair, since indentation has not changed. Find position of colon.
            call find_mapping_character(file%line,icolon,icolon_stop)
            if (icolon==-1) then
               call file%set_error('expected a key indicated by inline ": " or trailing :')
               return
            end if
         end do
      end if
   end function

   recursive function read_key_value_pair(file,icolon,icolon_stop) result(pair)
      class (type_file),intent(inout) :: file
      integer,          intent(in)    :: icolon,icolon_stop
      type (type_key_value_pair)      :: pair

      integer :: istop,baseindent

      istop = len_trim(file%line)

      pair%key = file%line(:icolon-1)
      if (icolon_stop==istop) then
         ! Colon ends the line; we need to read the value from the next line.
         baseindent = file%indent
         call file%next_line()
         if (file%has_error) return
         if (file%eof .or. file%indent<baseindent .or. (file%indent==baseindent .and. file%line(1:2)/='- ')) then
            ! Indentation equal to, or below, that of label (or file ends after label).
            ! That implies the value of the key-value pair is null.
            ! See YAML specification, section 7.2. Empty Nodes.
            allocate(type_null::pair%value)
         else
            ! Value on next line with higher indentation - read it.
            pair%value => read_value(file)
         end if
      else
         ! Value follows colon-space. Skip the label and read the value.
         file%line = file%line(icolon_stop+1:)
         file%indent = file%indent + icolon_stop
         pair%value => read_value(file)
      end if
   end function

   subroutine find_mapping_character(string,istart,istop)
      character(len=*),intent(in)  :: string
      integer,         intent(out) :: istart,istop
      integer                      :: i,length

      ! Default: mapping indicator not found.
      istart = -1
      istop = -1

      ! Search for mapping indicator
      length = len_trim(string)
      do i=1,length-1
         if (string(i:i+1)==': ') then
            ! Found "colon space" mapping indicator
            istart = i
            exit
         end if
      end do

      ! No mapping indicator found yet; check whether string ends with colon.
      if (istart==-1 .and. string(length:length)==':') istart = length

      ! If we have not found a mapping indicator by now, there isn't one: return.
      if (istart==-1) return

      ! Eliminate all trailing whitespace
      istop = istart
      do i=istart+1,length
         if (.not.is_whitespace(string(i:i))) then
            istop = i-1
            exit
         end if
      end do

      ! Eliminate all preceding whitespace
      do i=istart-1,1,-1
         if (.not.is_whitespace(string(i:i))) then
            istart = i+1
            exit
         end if
      end do
   end subroutine

   logical function is_whitespace(string)
      character(len=*),intent(in) :: string
      ! White space in YAML includes spaces and tabs only (NB tabs are not allowed in indentation!)
      is_whitespace = (string(1:1)==' '.or.string(1:1)==achar(9))
   end function

   subroutine set_error(file,error)
      class (type_file),intent(inout) :: file
      character(len=*), intent(in)    :: error
      file%error_message = error
      file%has_error = .true.
   end subroutine


   subroutine node_finalize(self)
      class (type_node),intent(inout) :: self
   end subroutine

   subroutine dictionary_reset_accessed(self)
      class (type_dictionary),intent(in) :: self
      type (type_key_value_pair),pointer :: pair
      pair => self%first
      do while (associated(pair))
         pair%accessed = .false.
         pair => pair%next
      end do
   end subroutine

   function dictionary_get(self,key) result(value)
      class (type_dictionary),intent(in) :: self
      character(len=*),       intent(in) :: key
      class(type_node),pointer           :: value

      type (type_key_value_pair),pointer :: pair

      nullify(value)
      pair => self%first
      do while (associated(pair))
         if (pair%key==key) exit
         pair => pair%next
      end do
      if (associated(pair)) then
         value => pair%value
         pair%accessed = .true.
      end if
   end function

   subroutine dictionary_set(self,key,value)
      class (type_dictionary),intent(inout) :: self
      character(len=*),       intent(in)    :: key
      class(type_node),pointer              :: value

      type (type_key_value_pair),pointer :: pair

      if (.not.associated(self%first)) then
         ! This will be the first pair.
         allocate(self%first)
         pair => self%first
      else
         ! Try to find a pair with the same key, or failing that, the last pair.
         pair => self%first
         do while (associated(pair%next))
            if (pair%key==key) exit
            pair => pair%next
         end do
         if (.not.pair%key==key) then
            ! Key did not exist yet, which must mean we are operating on the last existing pair.
            ! Append a new pair.
            allocate(pair%next)
            pair => pair%next
         else
            deallocate(pair%value)
         end if
      end if

      ! Store key and value.
      pair%key = key
      pair%value => value
   end subroutine

   subroutine dictionary_set_string(self,key,value)
      class (type_dictionary),intent(inout) :: self
      character(len=*),       intent(in)    :: key,value

      class (type_scalar),pointer :: scalar_node
      class (type_node),  pointer :: node

      allocate(scalar_node)
      scalar_node%string = value
      node => scalar_node
      call self%set(key,node)
   end subroutine

   subroutine value_dump(self,unit,indent)
      class (type_scalar),intent(in) :: self
      integer,            intent(in) :: unit,indent
      write (unit,'(a)') trim(self%string)
   end subroutine

   subroutine null_dump(self,unit,indent)
      class (type_null),intent(in) :: self
      integer,          intent(in) :: unit,indent
      write (unit,'(a)') 'null'
   end subroutine

   recursive subroutine dictionary_dump(self,unit,indent)
      class (type_dictionary),intent(in) :: self
      integer,                intent(in) :: unit,indent
      type (type_key_value_pair),pointer :: pair

      logical :: first

      first = .true.
      pair => self%first
      do while (associated(pair))
         if (first) then
            first = .false.
         else
            write (unit,'(a)',advance='NO') repeat(' ',indent)
         end if

         select type (value=>pair%value)
            class is (type_dictionary)
               write (unit,'(a)') trim(pair%key)//':'
               write (unit,'(a)',advance='NO') repeat(' ',indent+2)
               call value%dump(unit,indent+2)
            class is (type_list)
               write (unit,'(a)') trim(pair%key)//':'
               write (unit,'(a)',advance='NO') repeat(' ',indent+2)
               call value%dump(unit,indent+2)
            class default
               write (unit,'(a)',advance='NO') trim(pair%key)//': '
               call value%dump(unit,indent+len_trim(pair%key)+2)
         end select
         pair => pair%next
      end do
   end subroutine

   recursive subroutine dictionary_flatten(self,target,prefix)
      class (type_dictionary),intent(in)    :: self
      type (type_dictionary), intent(inout) :: target
      character(len=*),       intent(in)    :: prefix

      type (type_key_value_pair),pointer :: pair

      pair => self%first
      do while (associated(pair))
         select type (value=>pair%value)
            class is (type_scalar)
               call target%set_string(prefix//trim(pair%key),value%string)
            class is (type_dictionary)
               call value%flatten(target,prefix=prefix//trim(pair%key)//'/')
         end select
         pair => pair%next
      end do
   end subroutine

   function scalar_to_logical(self,default,success) result(value)
      class (type_scalar),intent(in)  :: self
      logical,            intent(in)  :: default
      logical,optional,   intent(out) :: success
      logical                         :: value

      integer :: ios

      value = default
      read(self%string,*,iostat=ios) value
      if (present(success)) success = (ios == 0)
   end function

   function scalar_to_integer(self,default,success) result(value)
      class (type_scalar),intent(in)  :: self
      integer,            intent(in)  :: default
      logical,optional,   intent(out) :: success
      integer                         :: value

      integer :: ios

      value = default
      read(self%string,*,iostat=ios) value
      if (present(success)) success = (ios == 0)
   end function

   function scalar_to_real(self,default,success) result(value)
      class (type_scalar),intent(in)  :: self
      real(real_kind),    intent(in)  :: default
      logical,optional,   intent(out) :: success
      real(real_kind)                 :: value

      integer :: ios

      value = default
      read(self%string,*,iostat=ios) value
      if (present(success)) success = (ios == 0)
   end function

   recursive subroutine node_set_path(self,path)
      class (type_node),intent(inout) :: self
      character(len=*), intent(in)    :: path
      self%path = path
   end subroutine

   recursive subroutine dictionary_set_path(self,path)
      class (type_dictionary),intent(inout) :: self
      character(len=*),       intent(in)    :: path

      type (type_key_value_pair),pointer :: pair

      self%path = path
      pair => self%first
      do while (associated(pair))
         call pair%value%set_path(trim(self%path)//'/'//trim(pair%key))
         pair => pair%next
      end do
   end subroutine

   function dictionary_get_scalar(self,key,required,error) result(scalar)
      class (type_dictionary),  intent(in) :: self
      character(len=*),         intent(in) :: key
      logical,                  intent(in) :: required
      type(type_error),pointer             :: error
      class (type_scalar),pointer          :: scalar

      class (type_node),pointer          :: node

      nullify(error)
      nullify(scalar)
      node => self%get(key)
      if (required.and..not.associated(node)) then
         allocate(error)
         error%message = trim(self%path)//' does not contain key "'//trim(key)//'".'
      end if
      if (associated(node)) then
         select type (node)
            class is (type_scalar)
               scalar => node
            class is (type_null)
               allocate(error)
               error%message = trim(node%path)//' must be set to a scalar value, not to null.'
            class is (type_dictionary)
               allocate(error)
               error%message = trim(node%path)//' must be set to a scalar value, not to a dictionary.'
            class is (type_list)
               allocate(error)
               error%message = trim(node%path)//' must be set to a scalar value, not to a list.'
         end select
      end if
   end function

   function dictionary_get_dictionary(self,key,required,error) result(dictionary)
      class (type_dictionary),  intent(in) :: self
      character(len=*),         intent(in) :: key
      logical,                  intent(in) :: required
      type(type_error),pointer             :: error
      class (type_dictionary),pointer      :: dictionary

      class (type_node),pointer :: node

      nullify(error)
      nullify(dictionary)
      node => self%get(key)
      if (required.and..not.associated(node)) then
         allocate(error)
         error%message = trim(self%path)//' does not contain key "'//trim(key)//'".'
      end if
      if (associated(node)) then
         select type (typed_node=>node)
            class is (type_null)
               allocate(dictionary)
               dictionary%path = node%path
            class is (type_dictionary)
               dictionary => typed_node
            class default
               allocate(error)
               error%message = trim(node%path)//' must be a dictionary.'
         end select
      end if
   end function

   function dictionary_get_list(self,key,required,error) result(list)
      class (type_dictionary),  intent(in) :: self
      character(len=*),         intent(in) :: key
      logical,                  intent(in) :: required
      type(type_error),pointer             :: error
      class (type_list),pointer            :: list

      class (type_node),pointer :: node

      nullify(error)
      nullify(list)
      node => self%get(key)
      if (required.and..not.associated(node)) then
         allocate(error)
         error%message = trim(self%path)//' does not contain key "'//trim(key)//'".'
      end if
      if (associated(node)) then
         select type (typed_node=>node)
            class is (type_null)
               allocate(list)
            class is (type_list)
               list => typed_node
            class default
               allocate(error)
               error%message = trim(node%path)//' must be a list.'
         end select
      end if
   end function

   function dictionary_get_string(self,key,default,error) result(value)
      class (type_dictionary),  intent(in) :: self
      character(len=*),         intent(in) :: key
      character(len=*),optional,intent(in) :: default
      type(type_error),pointer             :: error
      character(len=string_length)         :: value

      class(type_scalar),pointer           :: node

      if (present(default)) value = default
      node => self%get_scalar(key,.not.present(default),error)
      if (associated(node)) value = node%string
   end function

   function dictionary_get_logical(self,key,default,error) result(value)
      class (type_dictionary),  intent(in) :: self
      character(len=*),         intent(in) :: key
      logical,         optional,intent(in) :: default
      type(type_error),pointer             :: error
      logical                              :: value

      class (type_scalar),pointer          :: node
      logical                              :: success

      if (present(default)) value = default
      node => self%get_scalar(key,.not.present(default),error)
      if (associated(node)) then
         value = node%to_logical(value,success)
         if (.not.success) then
            allocate(error)
            error%message = trim(node%path)//' is set to "'//trim(node%string) &
                          //'", which cannot be interpreted as a Boolean value.'
         end if
      end if
   end function

   function dictionary_get_integer(self,key,default,error) result(value)
      class (type_dictionary),  intent(in) :: self
      character(len=*),         intent(in) :: key
      integer,         optional,intent(in) :: default
      type(type_error),pointer             :: error
      integer                              :: value

      class (type_scalar),pointer          :: node
      logical                              :: success

      if (present(default)) value = default
      node => self%get_scalar(key,.not.present(default),error)
      if (associated(node)) then
         value = node%to_integer(value,success)
         if (.not.success) then
            allocate(error)
            error%message = trim(node%path)//' is set to "'//trim(node%string)//'", which cannot be interpreted as an integer.'
         end if
      end if
   end function

   function dictionary_get_real(self,key,default,error) result(value)
      class (type_dictionary),  intent(in) :: self
      character(len=*),         intent(in) :: key
      real(real_kind), optional,intent(in) :: default
      type(type_error),pointer             :: error
      real(real_kind)                      :: value

      class (type_scalar),pointer          :: node
      logical                              :: success

      if (present(default)) value = default
      node => self%get_scalar(key,.not.present(default),error)
      if (associated(node)) then
         value = node%to_real(value,success)
         if (.not.success) then
            allocate(error)
            error%message = trim(node%path)//' is set to "'//trim(node%string)//'", which cannot be interpreted as a real number.'
         end if
      end if
   end function

   subroutine dictionary_finalize(self)
      class (type_dictionary),intent(inout) :: self

      type (type_key_value_pair),pointer :: pair, next

      pair => self%first
      do while (associated(pair))
         next => pair%next
         call pair%value%finalize()
         deallocate(pair%value)
         deallocate(pair)
         pair => next
      end do
      nullify(self%first)
   end subroutine dictionary_finalize

   subroutine list_append(self,node)
      class (type_list),intent(inout) :: self
      class(type_node),target         :: node

      type (type_list_item),pointer :: item

      if (.not.associated(self%first)) then
         ! This will be the first pair.
         allocate(self%first)
         self%first%node => node
      else
         ! Try to find a pair with the same key, or failing that, the last pair.
         item => self%first
         do while (associated(item%next))
            item => item%next
         end do
         allocate(item%next)
         item%next%node => node
      end if
   end subroutine list_append

   recursive subroutine list_dump(self,unit,indent)
      class (type_list),intent(in) :: self
      integer,          intent(in) :: unit,indent

      type (type_list_item),pointer :: item
      logical :: first

      first = .true.
      item => self%first
      do while (associated(item))
         if (first) then
            first = .false.
         else
            write (unit,'(a)',advance='NO') repeat(' ',indent)
         end if
         write (unit,'(a)',advance='NO') '- '
         call item%node%dump(unit,indent+2)
         item => item%next
      end do
   end subroutine list_dump

   recursive subroutine list_set_path(self,path)
      class (type_list),intent(inout) :: self
      character(len=*), intent(in)    :: path

      type (type_list_item),pointer :: item
      integer :: inode
      character(len=6) :: strindex

      self%path = path
      inode = 0
      item => self%first
      do while (associated(item))
         write (strindex,'(i0)') inode
         call item%node%set_path(trim(self%path)//'['//trim(strindex)//']')
         inode = inode + 1
         item => item%next
      end do
   end subroutine list_set_path

end module yaml
