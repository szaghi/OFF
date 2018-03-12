!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_XML_TagDerivedType Data_Type_XML_Tag
!> Module definition of Type_XML_Tag
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_XML_TagPrivateProcedure Data_Type_XML_Tag
!> Module definition of Type_XML_Tag
!> @}

!> @brief Module Data_Type_XML_Tag contains the definition of Type_XML_Tag and useful procedures for its handling.
module Data_Type_XML_Tag
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                        ! Integers and reals precision definition.
USE Data_Type_Varying_String, only: Type_Varying_String ! Definition of Type_Varying_String.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: xml_tags_count
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type defining Type_XML_Tag, a useful type for parsing XML file.
!> @note A valid XML tag must have the following syntax for a tag without a value (with only attributes):
!> @code
!>   <Tag_Name att#1_Name="att#1_val" att#2_Name="att#2_val"... att#Nt_Name="att#Nt_val"/>
!> @endcode
!> while a tag with a value must have the following syntax:
!> @code
!>   <Tag_Name att#1_Name="att#1_val" att#2_Name="att#2_val"... att#Nt_Name="att#Nt_val">Tag_value</Tag_Name>
!> @endcode
!> It is worth noting that the syntax is case sensitive and that the attributes are optional. Each attribute name must be followed
!> by '="' without any additional white spaces and its value must be termined by '"'. Each attribute is separated by a white
!> space. If the string member does not contain the tag_name no attributes are parsed.
!> @ingroup Data_Type_XML_TagDerivedType
type, public:: Type_XML_Tag
  type(Type_Varying_String)::              string              !< String containing the complete XML tag.
  type(Type_Varying_String)::              tag_name            !< Tag name.
  type(Type_Varying_String)::              tag_val             !< Tag value.
  type(Type_Varying_String), allocatable:: att_name(:)         !< Attributes names.
  type(Type_Varying_String), allocatable:: att_val(:)          !< Attributes values.
  type(Type_XML_Tag), pointer::            nested(:) => null() !< Nested tags into tag_val.
  contains
    procedure:: free                   ! Procedure for freeing dynamic memory.
    procedure:: alloc                  ! Procedure for allocating dynamic memory.
    procedure:: set                    ! Procedure for setting tag data.
    procedure:: get                    ! Procedure for getting the tag value and attributes from tag%string.
    procedure:: get_value              ! Procedure for getting the tag value from tag%string.
    procedure:: get_attributes         ! Procedure for getting the attributes values from tag%string.
    procedure:: parse_tag_name         ! Procedure for parsing the tag name contained into a string.
    procedure:: parse_attributes_names ! Procedure for parsing the tag attributes names contained into a string.
    procedure:: parse                  ! Procedure for parsing the tag contained into a string.
    procedure:: search_into_string     ! Procedure for searching the tag into a string.
    procedure:: search_into_other_tag  ! Procedure for searching the tag into the value of another tag.
    procedure:: get_nested             ! Procedure for getting XML tags nested to the current tag value.
    procedure:: stringify              ! Procedure for converting the whole tag insto a single string.
    final::     finalize               ! Procedure for freeing dynamic memory when finalizing.
    ! generic interfaces
    generic:: search => search_into_string, search_into_other_tag
    ! operators overloading
    generic:: assignment(=) => assign_tag
    ! private procedures
    procedure, pass(tag1), private:: assign_tag
endtype Type_XML_Tag
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_XML_TagPublicProcedure
  !> @{
  !> @brief Procedure for counting the number of xml tags (with respect the actual nesting level) into a string.
  function xml_tags_count(string) result(ntags)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: string      !< String containing the input.
  type(Type_XML_Tag)::       tag         !< Dummy XML tag.
  integer(I4P)::             ntags       !< Number of tags.
  integer(I4P)::             tstart,tend !< Counters for tracking string parsing.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  tstart = 1
  ntags = 0
  do while(tstart<len(string))
    call tag%parse(string=string(tstart:),tend=tend)
    tstart = tstart + tend
    if (allocated(tag%tag_name%vs)) then
      ntags = ntags + 1
      call tag%free
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction xml_tags_count
  !> @}

  !> @ingroup Data_Type_XML_TagPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  recursive subroutine free(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag !< XML tag.
  integer(I4P)::                       n   !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tag%string%free
  call tag%tag_name%free
  call tag%tag_val%free
  if (allocated(tag%att_name)) then
    call tag%att_name%free ; deallocate(tag%att_name)
  endif
  if (allocated(tag%att_val)) then
    call tag%att_val%free ; deallocate(tag%att_val)
  endif
  if (associated(tag%nested)) then
    do n=1,size(tag%nested)
      call tag%nested(n)%free
    enddo
    deallocate(tag%nested)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Procedure for freeing dynamic memory when finalizing.
  subroutine finalize(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_XML_tag), intent(INOUT):: tag !< XML tag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tag%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc(tag,att_name,att_val,Na)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag      !< XML tag.
  logical, optional,   intent(IN)::    att_name !< Flag for freeing attributes names array.
  logical, optional,   intent(IN)::    att_val  !< Flag for freeing attributes values array.
  integer(I4P),        intent(IN)::    Na       !< Number of attributes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(att_name)) then
    if (att_name) then
      if (allocated(tag%att_name)) then
        call tag%att_name%free ; deallocate(tag%att_name)
      endif
      allocate(tag%att_name(1:Na))
      ! call tag%att_name%free
    endif
  endif
  if (present(att_val)) then
    if (att_val) then
      if (allocated(tag%att_val)) then
        call tag%att_val%free ; deallocate(tag%att_val)
      endif
      allocate(tag%att_val(1:Na))
      ! call tag%att_val%free
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc

  !> @brief Procedure for setting tag data.
  !> @note If the att_name argument is passed all the attributes names must have the same number of characters.
  pure subroutine set(tag,string,tag_name,att_name,tstart,tend,search_into_string,search_into_other_tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag),          intent(INOUT):: tag                   !< XML tag.
  character(*),       optional, intent(IN)::    string                !< String containing the complete XML tag.
  character(*),       optional, intent(IN)::    tag_name              !< Tag name.
  character(*),       optional, intent(IN)::    att_name(:)           !< Attributes names.
  integer(I4P),       optional, intent(OUT)::   tstart                !< Starting index of tag inside the search string.
  integer(I4P),       optional, intent(OUT)::   tend                  !< Ending index of tag inside the search string.
  character(*),       optional, intent(IN)::    search_into_string    !< String into which searching for tag data.
  type(Type_XML_Tag), optional, intent(IN)::    search_into_other_tag !< Other tag into which searching for tag data.
  integer(I4P)::                                tstartd               !< Starting index of tag inside the search string.
  integer(I4P)::                                tendd                 !< Ending index of tag inside the search string.
  integer(I4P)::                                a                     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(string)) tag%string%vs = string
  if (present(tag_name)) tag%tag_name%vs = tag_name
  if (present(att_name)) then
    call tag%alloc(att_name=.true.,Na=size(att_name,dim=1))
    do a=1,size(att_name,dim=1)
      tag%att_name(a)%vs = att_name(a)
    enddo
  endif
  if (allocated(tag%tag_name%vs)) then
    if (present(search_into_string)) then
      call tag%search(string=search_into_string,tstart=tstartd,tend=tendd)
    elseif (present(search_into_other_tag)) then
      call tag%search(otag=search_into_other_tag,tstart=tstartd,tend=tendd)
    endif
  endif
  if (present(tstart)) tstart = tstartd
  if (present(tend  )) tend   = tendd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  !> @brief Procedure for getting the tag value and attributes from tag%string after tag_name and att_name have been set.
  !> @note It is worth noting that the leading and trailing white spaces of tag value and attributes are removed.
  elemental subroutine get(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag   !< XML tag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tag%get_value()
  call tag%get_attributes()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get

  !> @brief Procedure for getting the tag value from tag%string after tag_name and att_name have been set.
  !> @note It is worth noting that the leading and trailing white spaces of tag value are removed.
  elemental subroutine get_value(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag   !< XML tag.
  integer::                            c1,c2 !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(tag%string%vs)) then
    if (index(string=tag%string%vs,substring='<'//tag%tag_name%vs)>0) then
      c2 = index(string=tag%string%vs,substring='</'//tag%tag_name%vs//'>')
      if (c2>0) then ! parsing tag value
        c1 = index(string=tag%string%vs,substring='>')
        tag%tag_val%vs = trim(adjustl(tag%string%vs(c1+1:c2-1)))
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_value

  !> @brief Procedure for getting the attributes values from tag%string after tag_name and att_name have been set.
  !> @note It is worth noting that the leading and trailing white spaces of attributes values are removed.
  elemental subroutine get_attributes(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag     !< XML tag.
  integer::                            a,c1,c2 !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(tag%string%vs)) then
    if (index(string=tag%string%vs,substring='<'//tag%tag_name%vs)>0) then
      if (allocated(tag%att_name)) then ! parsing attributes
        call tag%alloc(att_val=.true.,Na=size(tag%att_name,dim=1))
        do a=1,size(tag%att_name,dim=1)
          c1 = index(string=tag%string%vs,substring=tag%att_name(a)%vs//'="')+len(tag%att_name(a)%vs)+2
          if (c1>len(tag%att_name(a)%vs)+2) then
            c2 = index(string=tag%string%vs(c1:),substring='"')
            if (c2>0) then
              tag%att_val(a)%vs = trim(adjustl(tag%string%vs(c1:c1+c2-2)))
            else
              call tag%att_val(a)%free
            endif
          else
            call tag%att_val(a)%free
          endif
        enddo
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_attributes

  !> @brief Procedure for parsing the tag name contained into a string.
  !> @note It is assumed that the first tag contained into the string is parsed, the others eventually present are omitted.
  !> Valid syntax are:
  !> + <tag_name att1="att1 val" att2="att2 val"...>...</tag_name>
  !> + <tag_name att1="att1 val" att2="att2 val".../>
  !> Inside the attributes value the symbols '<' and '>' are not allowed.
  elemental subroutine parse_tag_name(tag,tstart,tend,string)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag),    intent(INOUT):: tag     !< XML tag.
  integer(I4P), optional, intent(OUT)::   tstart  !< Starting index of tag inside the string.
  integer(I4P), optional, intent(OUT)::   tend    !< Ending index of tag inside the string.
  character(*),           intent(IN)::    string  !< String containing the input.
  integer(I4P)::                          tstartd !< Starting index of tag inside the string.
  integer(I4P)::                          tendd   !< Ending index of tag inside the string.
  character(len=1)::                      c1      !< Dummy string for parsing file.
  character(len=:), allocatable::         c2      !< Dummy string for parsing file.
  integer(I4P)::                          c,s     !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  tstartd = 0
  tendd   = 0
  c = 1
  Tag_Search: do while(c<=len(string))
    c1 = string(c:c)
    if (c1=='<') then
      tstartd = c
      c2 = c1
      Tag_Name: do while(c<len(string))
        c = c + 1 ; c1 = string(c:c)
        c2 = c2//c1
        if (c1=='>') then
          tendd = c
          exit Tag_Name
        endif
      enddo Tag_Name
      s = index(string=c2,substring=' ')
      if (s>0) then ! there are attributes
        tag%tag_name%vs = c2(2:s-1)
      else
        if (index(string=c2,substring='/>')>0) then ! self closing tag
          tag%tag_name%vs = c2(2:len(c2)-2)
        else
          tag%tag_name%vs = c2(2:len(c2)-1)
        endif
      endif
      exit Tag_Search
    endif
    c = c + 1
  enddo Tag_Search
  if (present(tstart)) tstart = tstartd
  if (present(tend  )) tend   = tendd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine parse_tag_name

  !> @brief Procedure for parsing the tag attributes names contained into a string.
  !> @note Valid syntax is:
  !> + att1="att1 val" att2="att2 val"...
  !> Inside the attributes value the symbols '<' and '>' are not allowed.
  elemental subroutine parse_attributes_names(tag,string)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag    !< XML tag.
  character(*),        intent(IN)::    string !< String containing the input.
  character(len=:), allocatable::      att    !< Dummy string for parsing file.
  integer(I4P)::                       c,a,Na !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Na = 0
  c = 1
  Att_Count: do while(c<=len(string))
    if (string(c:c)=='=') Na = Na + 1
    c = c + 1
  enddo Att_Count
  if (Na>0) then
    if (allocated(tag%att_name)) deallocate(tag%att_name) ; allocate(tag%att_name(1:Na))
    att = string
    c = 1
    a = 1
    Att_Search: do while(c<=len(att))
      if (att(c:c)=='=') then
        tag%att_name(a)%vs = trim(adjustl(att(:c-1)))
        att = att(c:)
        c = 1
        a = a + 1
      endif
      c = c + 1
    enddo Att_Search
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine parse_attributes_names

  !> @brief Procedure for parsing the tag contained into a string.
  !> @note It is assumed that the first tag contained into the string is parsed, the others eventually present are omitted.
  !> Valid syntax are:
  !> + <tag_name att1="att1 val" att2="att2 val"...>...</tag_name>
  !> + <tag_name att1="att1 val" att2="att2 val".../>
  !> Inside the attributes value the symbols '<' and '>' are not allowed.
  elemental subroutine parse(tag,tstart,tend,string)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag),    intent(INOUT):: tag     !< XML tag.
  integer(I4P), optional, intent(OUT)::   tstart  !< Starting index of tag inside the string.
  integer(I4P), optional, intent(OUT)::   tend    !< Ending index of tag inside the string.
  character(*),           intent(IN)::    string  !< String containing the input.
  integer(I4P)::                          tstartd !< Starting index of tag inside the string.
  integer(I4P)::                          tendd   !< Ending index of tag inside the string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  tstartd = 0
  tendd   = 0
  call tag%parse_tag_name(string=string,tstart=tstartd,tend=tendd)
  if (allocated(tag%tag_name%vs)) then
    if (index(string=string(tstartd:tendd),substring='=')>0) call tag%parse_attributes_names(string=string(tstartd:tendd))
    if (index(string=string,substring='</'//tag%tag_name%vs//'>')>0) &
      tendd = index(string=string,substring='</'//tag%tag_name%vs//'>') + len('</'//tag%tag_name%vs//'>') - 1
    tag%string%vs = string(tstartd:tendd)
    call tag%get
  endif
  if (present(tstart)) tstart = tstartd
  if (present(tend  )) tend   = tendd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine parse

  !> @brief Procedure for searching the tag (which name as been previously set) into a string.
  elemental subroutine search_into_string(tag,tstart,tend,string)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag),    intent(INOUT):: tag     !< XML tag.
  integer(I4P), optional, intent(OUT)::   tstart  !< Starting index of tag inside the string.
  integer(I4P), optional, intent(OUT)::   tend    !< Ending index of tag inside the string.
  character(*),           intent(IN)::    string  !< String containing the input.
  integer(I4P)::                          tstartd !< Starting index of tag inside the string.
  integer(I4P)::                          tendd   !< Ending index of tag inside the string.
  character(len=:), allocatable::         dstring !< Dummy string.
  type(Type_XML_Tag)::                    dtag    !< Dummy XML tag.
  logical::                               found   !< Flag for inquiring search result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  tstartd = 0
  tendd   = 0
  if (allocated(tag%tag_name%vs)) then
    dstring = string
    found = .false.
    Tag_Search: do while ((.not.found).or.(len(dstring)<len(tag%tag_name%vs)))
      call dtag%parse(string=dstring,tstart=tstartd,tend=tendd)
      if (tstartd==0.and.tendd==0) then
        exit Tag_Search ! no tag found
      else
        if (allocated(dtag%tag_name%vs)) then
          if (dtag%tag_name%vs==tag%tag_name%vs) then
            found = .true.
          else
            dstring = dstring(tendd+1:)
          endif
        endif
      endif
    enddo Tag_Search
  endif
  if (found) tag = dtag
  if (present(tstart)) tstart = tstartd
  if (present(tend  )) tend   = tendd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine search_into_string

  !> @brief Procedure for searching the tag (which name as been previously set) into the value of another tag.
  elemental subroutine search_into_other_tag(tag,tstart,tend,otag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag),    intent(INOUT):: tag     !< XML tag.
  integer(I4P), optional, intent(OUT)::   tstart  !< Starting index of tag inside the string.
  integer(I4P), optional, intent(OUT)::   tend    !< Ending index of tag inside the string.
  type(Type_XML_Tag),     intent(IN)::    otag    !< Other tag.
  integer(I4P)::                          tstartd !< Starting index of tag inside the string.
  integer(I4P)::                          tendd   !< Ending index of tag inside the string.
  character(len=:), allocatable::         dstring !< Dummy string.
  type(Type_XML_Tag)::                    dtag    !< Dummy XML tag.
  logical::                               found   !< Flag for inquiring search result.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  tstartd = 0
  tendd   = 0
  if (allocated(tag%tag_name%vs).and.(allocated(otag%tag_val%vs))) then
    dstring = otag%tag_val%vs
    found = .false.
    Tag_Search: do while ((.not.found).or.(len(dstring)<len(tag%tag_name%vs)))
      call dtag%parse(string=dstring,tstart=tstartd,tend=tendd)
      if (tstartd==0.and.tendd==0) then
        exit Tag_Search ! no tag found
      else
        if (allocated(dtag%tag_name%vs)) then
          if (dtag%tag_name%vs==tag%tag_name%vs) then
            found = .true.
          else
            dstring = dstring(tendd+1:)
          endif
        endif
      endif
    enddo Tag_Search
  endif
  if (found) tag = dtag
  if (present(tstart)) tstart = tstartd
  if (present(tend  )) tend   = tendd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine search_into_other_tag

  !> @brief Procedure for getting XML tags nested to the current tag value.
  subroutine get_nested(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag         !< XML tag.
  integer(I4P)::                       Nn          !< Number of nested tags.
  integer(I4P)::                       tstart,tend !< Counters for tracking string parsing.
  integer(I4P)::                       n           !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(tag%tag_val%vs)) then
    Nn = xml_tags_count(string=tag%tag_val%vs)
    if (Nn>0) then
      if (associated(tag%nested)) deallocate(tag%nested) ; allocate(tag%nested(1:Nn))
      tstart = 1
      n = 1
      do while(tstart<len(tag%tag_val%vs))
        call tag%nested(n)%parse(string=tag%tag_val%vs(tstart:),tend=tend)
        tstart = tstart + tend
        n = n + 1
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_nested

  !> @brief Procedure for converting the whole tag insto a single string.
  pure subroutine stringify(tag,pref,string)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag),           intent(IN)::  tag    !< XML tag.
  character(*),     optional,    intent(IN)::  pref   !< Prefixing string.
  character(len=:), allocatable, intent(OUT):: string !< Output string containing the whole tag.
  integer(I4P)::                               a      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = '' ; if (present(pref)) string = pref
  if (allocated(tag%string%vs)) then
    string = string//tag%string%vs
  else
    if (allocated(tag%tag_name%vs)) then
      string = string//'<'//tag%tag_name%vs
      if (allocated(tag%att_name).and.allocated(tag%att_val)) then
        if (size(tag%att_name)==size(tag%att_val)) then ! consistency check
          do a=1,size(tag%att_name)
            if (allocated(tag%att_name(a)%vs).and.allocated(tag%att_val(a)%vs)) &
              string = string//' '//tag%att_name(a)%vs//'="'//tag%att_val(a)%vs//'"'
          enddo
        endif
      endif
      if (allocated(tag%tag_val%vs)) then
        string = string//'>'//tag%tag_val%vs//'</'//tag%tag_name%vs//'>'
      else
        string = string//'/>'
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine stringify

  ! Assignment (=)
  !> @brief Procedure for assignment between two boundary conditions variables.
  elemental subroutine assign_tag(tag1,tag2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag1 !< LHS.
  type(Type_XML_Tag),  intent(IN)::    tag2 !< RHS.
  integer(I4P)::                       a    !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(tag2%string%vs  )) tag1%string%vs   = tag2%string%vs
  if (allocated(tag2%tag_name%vs)) tag1%tag_name%vs = tag2%tag_name%vs
  if (allocated(tag2%tag_val%vs )) tag1%tag_val%vs  = tag2%tag_val%vs
  if (allocated(tag2%att_name)) then
    if (allocated(tag1%att_name)) deallocate(tag1%att_name) ; allocate(tag1%att_name(1:size(tag2%att_name)))
    do a=1,size(tag2%att_name)
      tag1%att_name(a) = tag2%att_name(a)
    enddo
  endif
  if (allocated(tag2%att_val)) then
    if (allocated(tag1%att_val)) deallocate(tag1%att_val) ; allocate(tag1%att_val(1:size(tag2%att_val)))
    do a=1,size(tag2%att_val)
      tag1%att_val(a) = tag2%att_val(a)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_tag
  !> @}
endmodule Data_Type_XML_Tag
