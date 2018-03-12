!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_SpeciesDerivedType Data_Type_File_Species
!> Module definition of Type_File_Species
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_SpeciesPrivateProcedure Data_Type_File_Species
!> Module definition of Type_File_Species
!> @}

!> @brief Module Data_Type_File_Species contains the definition of Type_File_Species.
module Data_Type_File_Species
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                              ! Integers and reals precision definition.
USE Data_Type_File_Base, only: Type_File_Base ! Definition of Type_File_Base.
USE Data_Type_Species,   only: Type_Species   ! Definition of Type_Species.
USE Data_Type_XML_Tag,   only: Type_XML_Tag   ! Definition of Type_XML_Tag.
USE Lib_IO_Misc,         only: iostat_eor     ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Species.
  !> @note This options file is formatted by means of XML syntax. A validated file should be formatted as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <Species>
  !>     <Specie s="1"> <specie .../></Specie>
  !>     <Specie s="2"> <specie .../></Specie>
  !>     ...
  !>     <Specie s="Ns"><specie .../></Specie>
  !>   </Species>
  !> @endcode
  !> The main tag is 'Species' that contains the definition of each specie. Each specie if defined by a 'Specie' tag that is an
  !> 'overloading' tag of 'specie' one actually defining the specie components. The overloading tag 'Specie' introduces the
  !> attribute 's' defining the numeration of each specie into the global numeration. The total number of species defined is
  !> automatically computed parsing the whole file, thus it is not necessary to explicitly define it.
  !> The tag 'specie' is defined into 'Data_Type_Specie' module.
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
!> @ingroup Data_Type_File_SpeciesDerivedType
type, public, extends(Type_File_Base):: Type_File_Species
  contains
    procedure:: load  => load_species ! Procedure for loading species file.
    procedure:: save  => save_species ! Procedure for saving species file.
endtype Type_File_Species
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_SpeciesPrivateProcedure
  !> @{
  !> @brief Procedure for loading species file.
  subroutine load_species(file_d,species)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Species), intent(INOUT):: file_d  !< File data.
  type(Type_Species),       intent(INOUT):: species !< Species.
  character(len=:), allocatable::           stream  !< String containing the file data as a single stream.
  type(Type_XML_Tag)::                      tag     !< Main XML tag.
  integer(I4P)::                            s,ss    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg)
    call file_d%get_stream(stream=stream)!; if (iostat/=0) return
    call tag%set(tag_name='Species',search_into_string=stream)
    call tag%get_nested
    if (associated(tag%nested)) then
      species%Ns = size(tag%nested) ; call species%alloc
      do s=1,size(tag%nested)
        call tag%nested(s)%get_nested
        ss = cton(str=tag%nested(s)%att_val(1)%vs,knd=1_I4P)
        call species%heats(s)%load_str_xml(str_xml=tag%nested(s)%nested(1)%string%vs)
      enddo
    endif
    call species%compute_Npc
  endassociate
  call file_d%fallback
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_species

  !> @brief Procedure for saving species file.
  subroutine save_species(file_d,species)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Species), intent(INOUT):: file_d  !< File data.
  type(Type_Species),       intent(IN)::    species !< Species.
  integer(I4P)::                            s       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg)
    call file_d%open(ascii=.true.,action='WRITE') ; if (iostat/=0) return
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<?xml version="1.0"?>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<Species>'
    do s=1,species%Ns
      write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <Specie s="'//trim(str(n=s))//'">'//&
                                                          species%heats(s)%save_str_xml()//'</Specie>'
    enddo
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'</Species>'
    call file_d%close
  endassociate
  call file_d%fallback
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_species
  !> @}
endmodule Data_Type_File_Species
