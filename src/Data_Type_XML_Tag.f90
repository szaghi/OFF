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
  type(Type_Varying_String)::              string      !< String containing the complete XML tag.
  type(Type_Varying_String)::              tag_name    !< Tag name.
  type(Type_Varying_String)::              tag_val     !< Tag value.
  type(Type_Varying_String), allocatable:: att_name(:) !< Attributes names.
  type(Type_Varying_String), allocatable:: att_val(:)  !< Attributes values.
  contains
    procedure:: free           ! Procedure for freeing dynamic memory.
    procedure:: alloc          ! Procedure for allocating dynamic memory.
    procedure:: set            ! Procedure for setting tag data.
    procedure:: get_value      ! Procedure for getting the tag value.
    procedure:: get_attributes ! Procedure for getting the attributes values.
    final::     finalize       ! Procedure for freeing dynamic memory when finalizing.
endtype Type_XML_Tag
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_XML_TagPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag !< XML tag.
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

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
      call tag%att_name%free
    endif
  endif
  if (present(att_val)) then
    if (att_val) then
      if (allocated(tag%att_val)) then
        call tag%att_val%free ; deallocate(tag%att_val)
      endif
      allocate(tag%att_val(1:Na))
      call tag%att_val%free
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc

  !> @brief Procedure for setting tag data.
  !> @note If the att_name argument is passed all the attributes names must have the same number of charcters.
  pure subroutine set(tag,string,tag_name,att_name)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag),    intent(INOUT):: tag         !< XML tag.
  character(*), optional, intent(IN)::    string      !< String containing the complete XML tag.
  character(*), optional, intent(IN)::    tag_name    !< Tag name.
  character(*), optional, intent(IN)::    att_name(:) !< Attributes names.
  integer(I4P)::                          a           !< Counter.
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  !> @brief Procedure for getting the tag value.
  !> @note It is worth noting that the leading and trailing white spaces of tag value are removed.
  elemental subroutine get_value(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag   !< XML tag.
  integer::                            c1,c2 !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (index(string=tag%string%vs,substring='<'//tag%tag_name%vs)>0) then
    c2 = index(string=tag%string%vs,substring='</'//tag%tag_name%vs//'>')
    if (c2>0) then ! parsing tag value
      c1 = index(string=tag%string%vs,substring='>')
      tag%tag_val%vs = trim(adjustl(tag%string%vs(c1+1:c2-1)))
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_value

  !> @brief Procedure for getting the attributes values.
  !> @note It is worth noting that the leading and trailing white spaces of attributes values are removed.
  elemental subroutine get_attributes(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_XML_Tag), intent(INOUT):: tag     !< XML tag.
  integer::                            a,c1,c2 !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_attributes

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(tag)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_XML_tag), intent(INOUT):: tag !< XML tag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tag%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
  !> @}
endmodule Data_Type_XML_Tag
