!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_SpecieDerivedType Data_Type_Specie
!> Module definition of fluid specie, Type_Specie
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_SpeciePrivateProcedure Data_Type_Specie
!> Module definition of fluid specie, Type_Specie
!> @}

!> @brief Module Data_Type_Specie contains the definition of Type_Specie, that defines the fluid specie properties.
module Data_Type_Specie
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                          ! Integers and reals precision definition.
USE Data_Type_XML_Tag, only: Type_XML_Tag ! Definition of Type_XML_Tag.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_Specie.
!> @ingroup Data_Type_SpecieDerivedType
type, public:: Type_Specie
  real(R8P):: cp = 0._R8P !< Specific heat at constant pressure.
  real(R8P):: cv = 0._R8P !< Specific heat at constant volume.
  contains
    procedure:: load_str_xml               ! Procedure for loading specie from a string in XML format.
    procedure:: save_str_xml               ! Procedure for saving specie into a string in XML format.
    procedure:: print => print_specie_self ! Procedure for printing specie with a pretty format.
    ! operators overloading
    generic:: assignment(=) => assign_specie
    ! private procedures
    procedure, pass(spec1), private:: assign_specie
endtype Type_Specie
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_SpeciePrivateProcedure
  !> @{
  !> @brief Procedure for loading specie from a string formatted in XML syntax.
  !> @note The XML string of specie data must have the following syntax:
  !> @code
  !>   <specie Cp="..." Cv="..."/>
  !> @endcode
  !> The tag name is 'specie' and it contains the specie definition. This definition is performed by means of a combination of two
  !> independent parameters (defined into tag attributes) that can be arbitrarly chosen in the following list:
  !> @code
  !>   Cp="..."
  !>   Cv="..."
  !>   Gamma="..."
  !>   R="..."
  !> @endcode
  !> corresponding to the specific heat a constant pressure (Cp), the specific heat a constant volume (Cv), the specific heats
  !> ratio (Gamma=Cp/Cv) and the specific heats difference (the fluid constant, R=cp-cv).
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_str_xml(specie,str_xml)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Specie),     intent(INOUT):: specie  !< Specie data.
  character(*),           intent(IN)::    str_xml !< String containing primitives data in XML syntax.
  type(Type_XML_Tag)::                    tag     !< Tag parsing data.
  real(R8P)::                             cp,cv   !< Specific heats.
  real(R8P)::                             g       !< Specific heats ratio (cp/cv).
  real(R8P)::                             R       !< Specific heats difference (cp-cv, fluid constant).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tag%set(string=str_xml,tag_name='specie')
  call tag%alloc(att_name=.true.,Na=4)
  tag%att_name(1)%vs='Cp'
  tag%att_name(2)%vs='Cv'
  tag%att_name(3)%vs='Gamma'
  tag%att_name(4)%vs='R'
  call tag%get_attributes
  if     (allocated(tag%att_val(1)%vs).and.allocated(tag%att_val(2)%vs)) then
      ! Cp and Cv
      cp = cton(str=trim(adjustl(tag%att_val(1)%vs)),knd=1._R8P)
      cv = cton(str=trim(adjustl(tag%att_val(2)%vs)),knd=1._R8P)
  elseif (allocated(tag%att_val(3)%vs).and.allocated(tag%att_val(4)%vs)) then
      ! Gamma and R
      g  = cton(str=trim(adjustl(tag%att_val(3)%vs)),knd=1._R8P)
      R  = cton(str=trim(adjustl(tag%att_val(4)%vs)),knd=1._R8P)
      cv = R/(g-1._R_P)
      cp = g*cv
  elseif (allocated(tag%att_val(3)%vs).and.allocated(tag%att_val(1)%vs)) then
      ! Gamma and Cp
      g  = cton(str=trim(adjustl(tag%att_val(3)%vs)),knd=1._R8P)
      cp = cton(str=trim(adjustl(tag%att_val(1)%vs)),knd=1._R8P)
      cv = cp/g
  elseif (allocated(tag%att_val(3)%vs).and.allocated(tag%att_val(2)%vs)) then
      ! Gamma and Cv
      g  = cton(str=trim(adjustl(tag%att_val(3)%vs)),knd=1._R8P)
      cv = cton(str=trim(adjustl(tag%att_val(2)%vs)),knd=1._R8P)
      cp = cv*g
  elseif (allocated(tag%att_val(4)%vs).and.allocated(tag%att_val(1)%vs)) then
      ! R and Cp
      R  = cton(str=trim(adjustl(tag%att_val(4)%vs)),knd=1._R8P)
      cp = cton(str=trim(adjustl(tag%att_val(1)%vs)),knd=1._R8P)
      cv = cp-R
  elseif (allocated(tag%att_val(4)%vs).and.allocated(tag%att_val(2)%vs)) then
      ! R and Cv
      R  = cton(str=trim(adjustl(tag%att_val(4)%vs)),knd=1._R8P)
      cv = cton(str=trim(adjustl(tag%att_val(2)%vs)),knd=1._R8P)
      cp = cv+R
  endif
  specie%cp = cp
  specie%cv = cv
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_str_xml

  !> @brief Procedure for saving specie into a string formatted in XML syntax.
  !> @note The definition of the XML string syntax is described in the 'load_str_xml' documentation.
  elemental function save_str_xml(specie) result(str_xml)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Specie), intent(IN):: specie  !< Specie data.
  character(len=:), allocatable::  str_xml !< String containing primitives data in XML syntax.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  str_xml = '<specie Cp="'//trim(str(n=specie%cp))//'" Cv="'//trim(str(n=specie%cv))//'"/>'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_str_xml

  !> @brief Procedure for printing specie with a pretty format.
  subroutine print_specie_self(spec,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Specie),     intent(IN)::  spec    !< Specie.
  character(*), optional, intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::       prefd   !< Prefixing string.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' cp='//str(n=spec%cp)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' cv='//str(n=spec%cv)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_specie_self

  ! Assignment (=)
  !> @brief Procedure for assignment between two species variables.
  elemental subroutine assign_specie(spec1,spec2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Specie), intent(INOUT):: spec1
  type(Type_Specie),  intent(IN)::    spec2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  spec1%cp = spec2%cp
  spec1%cv = spec2%cv
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_specie
  !> @}
endmodule Data_Type_Specie
