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
  subroutine load_species(file_d,species)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Species), intent(INOUT):: file_d  !< File data.
  type(Type_Species),       intent(INOUT):: species !< Species.
  type(Type_XML_Tag)::                      tag     !< Tag parsing data.
  integer(I4P)::                            s       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg)
    allocate(character(len=1000):: tag%string%vs)
    call file_d%open(ascii=.true.,action='READ') ; if (iostat/=0) return
    species%Ns = 0
    Species_Count: do
      read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=10)tag%string%vs ; if (iostat/=0.and.iostat/=iostat_eor) return
      if (index(string=tag%string%vs,substring='<Specie ')>0) then
        species%Ns = species%Ns + 1
      endif
    enddo Species_Count
    10 rewind(unit=unit) ; call species%alloc
    call tag%set(tag_name='Specie',att_name=['s'])
    Read_Loop: do
      read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=20)tag%string%vs ; if (iostat/=0.and.iostat/=iostat_eor) return
      if (index(string=tag%string%vs,substring='<Specie ')>0) then
        call tag%get_attributes ; call tag%get_value
        s = cton(str=trim(adjustl(tag%att_val(1)%vs)),knd=1_I4P)
        call species%heats(s)%load_str_xml(str_xml=tag%tag_val%vs)
      endif
    enddo Read_Loop
    20 call file_d%close
    call species%compute_Npc
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_species

  !> @brief Procedure for saving species file.
  !> @note The definition of the species file syntax is described in the 'load_species' documentation.
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_species
  !> @}
endmodule Data_Type_File_Species
