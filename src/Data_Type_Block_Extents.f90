!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Block_ExtentsDerivedType Data_Type_Block_Extents
!> Module definition of Type_Block_Extents
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_Block_ExtentsInterface Data_Type_Block_Extents
!> Module definition of Type_Block_Extents
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Block_ExtentsPrivateProcedure Data_Type_Block_Extents
!> Module definition of Type_Block_Extents
!> @}

!> @brief Module Data_Type_Block_Extents contains the definition of Type_Block_Extents and useful procedures for its handling.
!> Type_Block_Extents contains all the data defining the main extents of block.
module Data_Type_Block_Extents
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision              ! Integers and reals precision definition.
USE Data_Type_Vector,  only: Type_Vector ! Definition of Type_Vector.
USE Data_Type_XML_Tag, only: Type_XML_Tag ! Definition of Type_XML_Tag.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the main extents of a block.
!> @ingroup Data_Type_Block_ExtentsDerivedType
type, public:: Type_Block_Extents
  type(Type_Vector):: emin !< Coordinates of minimum abscissa.
  type(Type_Vector):: emax !< Coordinates of maximum abscissa.
  contains
    procedure:: load_str_xml               ! Procedure for loading extents from a string in XML format.
    procedure:: save_str_xml               ! Procedure for saving extents into a string in XML format.
    procedure:: print  => print_block_exts ! Procedure for printing extents with a pretty format.
    ! operators overloading
    generic:: assignment(=) => assign_blk_exts
    ! private procedures
    procedure, pass(blk1), private:: assign_blk_exts
  endtype Type_Block_Extents
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Block_ExtentsPrivateProcedure
  !> @{
  !> @brief Procedure for loading block extents from a string formatted in XML syntax.
  !> @note The XML string of block extents data must have the following syntax:
  !> @code
  !>   <extents xmin="#xmin" ymin="#ymin" zmin="#zmin" xmax="#xmax" ymax="#ymax" zmax="#zmax"/>
  !> @endcode
  !> The attributes (x/y/z)min and (x/y/z)max define block extentes, i.e. the minimum and maximum abscissa coordinates of block.
  !> The order of attributes in not influent. It is worth nothing that the syntax of tag names and attributes is case
  !> sensitive, whereas the syntax of their values is case insensitive.
  subroutine load_str_xml(block_exts,str_xml)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Extents), intent(INOUT):: block_exts !< Block extents.
  character(*),              intent(IN)::    str_xml    !< String containing block extents data in XML syntax.
  type(Type_XML_Tag)                         tag        !< XML tag data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tag%set(string=str_xml,tag_name='extents',att_name=['xmin','ymin','zmin',&
                                                           'xmax','ymax','zmax']) ! setting xml tag data
  call tag%get_attributes                                                         ! parsing xml tag attributes
  ! loading extents data from xml tag attributes values
  block_exts%emin%x = cton(str=tag%att_val(1)%vs,knd=1._R8P)
  block_exts%emin%y = cton(str=tag%att_val(2)%vs,knd=1._R8P)
  block_exts%emin%z = cton(str=tag%att_val(3)%vs,knd=1._R8P)
  block_exts%emax%x = cton(str=tag%att_val(4)%vs,knd=1._R8P)
  block_exts%emax%y = cton(str=tag%att_val(5)%vs,knd=1._R8P)
  block_exts%emax%z = cton(str=tag%att_val(6)%vs,knd=1._R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_str_xml

  !> @brief Procedure for saving block extents into a string formatted in XML syntax.
  !> @note The definition of the XML string syntax is described in the 'load_str_xml' documentation.
  pure function save_str_xml(block_exts) result(str_xml)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Extents), intent(IN):: block_exts !< Block extents.
  character(len=:), allocatable::         str_xml    !< String containing block extents data in XML syntax.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  str_xml = '<extents xmin="'//trim(str(n=block_exts%emin%x))//&
                   '" ymin="'//trim(str(n=block_exts%emin%y))//&
                   '" zmin="'//trim(str(n=block_exts%emin%z))//&
                   '" xmax="'//trim(str(n=block_exts%emax%x))//&
                   '" ymax="'//trim(str(n=block_exts%emax%y))//&
                   '" zmax="'//trim(str(n=block_exts%emax%z))//'"/>'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_str_xml

  !> @brief Procedure for printing in a pretty ascii format the components of type Type_Block_Extents.
  subroutine print_block_exts(block_exts,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Extents), intent(IN)::  block_exts !< Block extents.
  character(*), optional,    intent(IN)::  pref       !< Prefixing string.
  integer(I4P), optional,    intent(OUT):: iostat     !< IO error.
  character(*), optional,    intent(OUT):: iomsg      !< IO error message.
  integer(I4P),              intent(IN)::  unit       !< Logic unit.
  character(len=:), allocatable::          prefd      !< Prefixing string for outputs.
  integer(I4P)::                           iostatd    !< IO error.
  character(500)::                         iomsgd     !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' xmin='//trim(str(n=block_exts%emin%x))//&
                                                               ' ymin='//trim(str(n=block_exts%emin%y))//&
                                                               ' zmin='//trim(str(n=block_exts%emin%z))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' xmax='//trim(str(n=block_exts%emax%x))//&
                                                               ' ymax='//trim(str(n=block_exts%emax%y))//&
                                                               ' zmax='//trim(str(n=block_exts%emax%z))
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_block_exts

  ! Assignment (=)
  !> @brief Procedure for assignment between two blocks extents variables.
  elemental subroutine assign_blk_exts(blk1,blk2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Extents), intent(INOUT):: blk1
  type(Type_Block_Extents),  intent(IN)::    blk2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  blk1%emin = blk2%emin
  blk1%emax = blk2%emax
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_blk_exts
  !> @}
endmodule Data_Type_Block_Extents
