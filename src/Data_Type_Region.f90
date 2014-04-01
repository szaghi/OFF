!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_RegionDerivedType Data_Type_Region
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_RegionInterface Data_Type_Region
!> Module definition of Type_Region
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_RegionPrivateProcedure Data_Type_Region
!> Module definition of Type_Region
!> @}

!> @brief This module contains the definition of Type_Region and its procedures.
!> Type_Region is a derived type containing the description of a particular region useful for initial conditions construction.
module Data_Type_Region
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                              ! Integers and reals precision definition.
USE Data_Type_Primitive, only: Type_Primitive ! Definition of Type_Primitive.
USE Data_Type_Vector,    only: Type_Vector    ! Definition of Type_Vector.
USE Data_Type_XML_Tag,   only: Type_XML_Tag   ! Definition of Type_XML_Tag.
USE Lib_IO_Misc,         only: Upper_Case     ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing the description of a particular region useful for initial conditions construction.
!> @ingroup Data_Type_RegionDerivedType
type, public:: Type_Region
  character(len=:), allocatable:: shtype          !< Type of the region shape: 'CYLINDER-X', 'CYLINDER-Y', 'CYLINDER-Z','SPHERE',...
  type(Type_Vector)::             center          !< Center of the region.
  real(R8P)::                     radius = 0._R8P !< Radius of the region.
  real(R8P)::                     height = 0._R8P !< Height of the region.
  type(Type_Primitive)::          prim            !< Primitive variables.
  contains
    procedure:: load_shape_str_xml    ! Procedure for loading region shape from a string in XML format.
    procedure:: save_shape_str_xml    ! Procedure for saving region shape into a string in XML format.
    procedure:: print => print_region ! Procedure for printing region with a pretty format.
    ! operators overloading
    generic:: assignment(=) => assign_self
    ! private procedures
    procedure, pass(self1), private:: assign_self
endtype Type_Region
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_RegionPrivateProcedure
  !> @{
  !> @brief Procedure for loading region shape from a string formatted in XML syntax.
  !> @note The XML string of region shape data must have the following syntax:
  !> @code
  !>   <shape tp="type_of_region" center="#x #y #z" radius="#radius" height="#height"/>
  !> @endcode
  !> The order of attributes in not influent. It is worth nothing that the syntax of tag names and attributes is case
  !> sensitive, whereas the syntax of their values is case insensitive.
  subroutine load_shape_str_xml(region,str_xml)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Region), intent(INOUT):: region  !< Region data.
  character(*),       intent(IN)::    str_xml !< String containing primitives data in XML syntax.
  type(Type_XML_Tag)                  tag     !< XML tag data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tag%set(string=str_xml,tag_name='shape',att_name=['shtype','center','radius','height']) ! setting xml tag data
  call tag%get_attributes                                                                      ! parsing xml tag attributes
  ! loading region data from xml tag attributes values
  region%shtype = trim(adjustl(Upper_Case(tag%att_val(1)%vs)))
  read(tag%att_val(2)%vs,*)region%center%x,region%center%y,region%center%z
  region%radius = cton(str=tag%att_val(3)%vs,knd=1._R8P)
  region%height = cton(str=tag%att_val(4)%vs,knd=1._R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_shape_str_xml

  !> @brief Procedure for saving region shape into a string formatted in XML syntax.
  !> @note The definition of the XML string syntax is described in the 'load_shape_str_xml' documentation.
  elemental function save_shape_str_xml(region) result(str_xml)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Region), intent(IN):: region  !< Region data.
  character(len=:), allocatable::  str_xml !< String containing region data in XML syntax.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  str_xml='<shape  shtype="'//region%shtype//                                                                                     &
                '" center="'//trim(str(n=region%center%x))//' '//trim(str(n=region%center%y))//' '//trim(str(n=region%center%z))//&
                '" radius="'//trim(str(n=region%radius))//'" height="'//trim(str(n=region%height))//'"/>'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_shape_str_xml

  !> @brief Procedure for printing region with a pretty format.
  subroutine print_region(region,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Region),     intent(IN)::  region  !< Region data.
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
  if (allocated(region%shtype)) write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' shape type   ='//region%shtype
                                write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' center(x,y,z)='//&
                                  str(n=region%center%x)//' '//str(n=region%center%y)//' '//str(n=region%center%z)
                                write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' radius       ='//str(n=region%radius)
                                write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' height       ='//str(n=region%height)
               call region%prim%print(unit=unit,          iostat=iostatd,iomsg=iomsgd,pref=prefd)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_region

  ! Assignment (=)
  !> @brief Procedure for assignment between two selfs.
  elemental subroutine assign_self(self1,self2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Region), intent(INOUT):: self1
  type(Type_Region),  intent(IN)::    self2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self2%shtype)) self1%shtype = self2%shtype
                               self1%center = self2%center
                               self1%radius = self2%radius
                               self1%prim   = self2%prim
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_self
  !> @}
endmodule Data_Type_Region
