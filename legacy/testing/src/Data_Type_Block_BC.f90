!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Block_BCDerivedType Data_Type_Block_BC
!> Module definition of Type_Block_BC
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_Block_BCInterface Data_Type_Block_BC
!> Module definition of Type_Block_BC
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Block_BCPrivateProcedure Data_Type_Block_BC
!> Module definition of Type_Block_BC
!> @}

!> @brief Module Data_Type_Block_BC contains the definition of Type_Block_BC and useful procedures for its handling.
!> Type_Block_BC contains all the data defining the block-face boundary conditions of a block.
module Data_Type_Block_BC
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                         ! Integers and reals precision definition.
USE Data_Type_BC                         ! Definition of Type_BC.
USE Data_Type_XML_Tag                    ! Definition of Type_XML_Tag.
USE Lib_Strings, only: Upper_Case,unique ! Library for strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the block-faces boundary conditions.
!> @ingroup Data_Type_Block_BCDerivedType
type, public:: Type_Block_BC
  type(Type_BC)::                 F(1:6)        !< Block face boundary conditions.
  integer(I4P)::                  Nijk_adj(1:6) !< N(ijk) of adjacent blocks.
  character(len=:), allocatable:: orientation   !< IJK axis orientation.
  contains
    procedure:: free  => free_block_bc  ! Procedure for freeing dynamic memory.
    procedure:: alloc => alloc_block_bc ! Procedure for allocating dynamic memory.
    procedure:: load  => load_block_bc  ! Procedure for loading boundary conditions.
    procedure:: save  => save_block_bc  ! Procedure for saving boundary conditions.
    procedure:: print => print_block_bc ! Procedure for printing boundary conditions with a pretty format.
    procedure:: load_str_xml            ! Procedure for loading boundary conditions from a string in XML format.
    procedure:: save_str_xml            ! Procedure for saving boundary conditions into a string in XML format.
    final::     finalize                ! Procedure for freeing dynamic memory when finalizing.
    ! operators overloading
    generic:: assignment(=) => assign_blk_bc
    ! private procedures
    procedure, pass(blk1), private:: assign_blk_bc
  endtype Type_Block_BC
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Block_BCPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_block_bc(block_bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_BC), intent(INOUT):: block_bc !< Block boundary conditions.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call block_bc%F%free
  if (allocated(block_bc%orientation)) deallocate(block_bc%orientation)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_block_bc

  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine finalize(block_bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Block_BC), intent(INOUT):: block_bc !< Block boundary conditions.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call block_bc%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc_block_bc(block_bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_BC), intent(INOUT):: block_bc !< Block boundary conditions.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call block_bc%F%alloc
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_block_bc

  !> @brief Procedure for loading boundary conditions.
  subroutine load_block_bc(block_bc,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Bc),   intent(INOUT):: block_bc !< Block boundary conditions.
  integer(I8P), optional, intent(IN)::    pos      !< Position specifier.
  integer(I4P), optional, intent(OUT)::   iostat   !< IO error.
  character(*), optional, intent(OUT)::   iomsg    !< IO error message.
  integer(I4P),           intent(IN)::    unit     !< Logic unit.
  integer(I4P)::                          iostatd  !< IO error.
  character(500)::                        iomsgd   !< Temporary variable for IO error message.
  integer(I4P)::                          f        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    call block_bc%F(1)%load(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)
  else
    call block_bc%F(1)%load(unit=unit,        iostat=iostatd,iomsg=iomsgd)
  endif
  do f=2,size(block_bc%F,dim=1)
    call block_bc%F(f)%load(unit=unit,        iostat=iostatd,iomsg=iomsgd)
  enddo
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_block_bc

  !> @brief Procedure for saving boundary conditions.
  subroutine save_block_bc(block_bc,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_block_bc),   intent(IN)::  block_bc !< Block boundary conditions.
  integer(I8P), optional, intent(IN)::  pos      !< Position specifier.
  integer(I4P), optional, intent(OUT):: iostat   !< IO error.
  character(*), optional, intent(OUT):: iomsg    !< IO error message.
  integer(I4P),           intent(IN)::  unit     !< Logic unit.
  integer(I4P)::                        iostatd  !< IO error.
  character(500)::                      iomsgd   !< Temporary variable for IO error message.
  integer(I4P)::                        f        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    call block_bc%F(1)%save(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)
  else
    call block_bc%F(1)%save(unit=unit,        iostat=iostatd,iomsg=iomsgd)
  endif
  do f=2,size(block_bc%F,dim=1)
    call block_bc%F(f)%save(unit=unit,        iostat=iostatd,iomsg=iomsgd)
  enddo
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_block_bc

  !> @brief Procedure for loading block boundary conditions from a string formatted in XML syntax.
  !> @note The XML string of block boundary conditions data must have the following syntax:
  !> @code
  !>   <boundaryconditions imin="#imin [d]" imax="#imax [d]" jmin="#jmin [d]" jmax="#jmax [d]" kmin="#kmin [d]" kmax="#kmax [d]"/>
  !> @endcode
  !> The attributes (i/j/k)min and (i/j/k)max define block-faces boundary conditions. Their values must contain a valid BC id ascii
  !> descriptor (defined into Data_Type_BC module). For some types of boundary conditions the values must also contain an additional
  !> data ([d]). For example a valid xml string is the following:
  !> @code
  !>   <boundaryconditions imin="EXT" imax="ADJ 2" jmin="EXT" jmax="EXT" kmin="EXT" kmax="EXT"/>
  !> @endcode
  !> The order of attributes in not influent. It is worth nothing that the syntax of tag names and attributes is case
  !> sensitive, whereas the syntax of their values is case insensitive.
  !> sensitive, whereas the syntax of their values is case insensitive.
  subroutine load_str_xml(block_bc,str_xml)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_BC), intent(INOUT):: block_bc !< Block boundary conditions.
  character(*),         intent(IN)::    str_xml  !< String containing block boundary conditions data in XML syntax.
  type(Type_XML_Tag)                    tag      !< XML tag data.
  integer(I4P)::                        f        !< Faces counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tag%set(string=str_xml,tag_name='boundaryconditions',&
                              att_name=['imin','imax','jmin','jmax','kmin','kmax']) ! setting xml tag data
  call tag%get_attributes                                                           ! parsing xml tag attributes
  do f=1,6
    tag%att_val(f)%vs = trim(unique(string=tag%att_val(f)%vs,substring=' '))        ! removing the eventual multiples white spaces
    ! loading boundary conditions data from xml tag attributes values
    call block_bc%F(f)%str2id(bc_str=Upper_Case(tag%att_val(f)%vs(1:3)))
    call block_bc%F(f)%alloc
    if     (block_bc%F(f)%is_adj()) then
      block_bc%F(f)%adj%ID = cton(str=tag%att_val(f)%vs(4:),knd=1_I8P)
    elseif (block_bc%F(f)%is_inf()) then
      block_bc%F(f)%inf   = cton(str=tag%att_val(f)%vs(4:),knd=1_I4P)
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_str_xml

  !> @brief Procedure for saving block boundary conditions into a string formatted in XML syntax.
  !> @note The definition of the XML string syntax is described in the 'load_str_xml' documentation.
  pure function save_str_xml(block_bc) result(str_xml)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_BC), intent(IN):: block_bc      !< Block boundary conditions.
  character(len=:), allocatable::    str_xml       !< String containing block boundary conditions data in XML syntax.
  integer(I4P)::                     f             !< Faces counters.
  character(4)::                     att_name(1:6) !< List of attributes names.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  att_name = ['imin','imax','jmin','jmax','kmin','kmax']
  str_xml = '<boundaryconditions'
  do f=1,6
    if     (block_bc%F(f)%is_adj()) then
      str_xml = trim(str_xml)//' '//att_name(f)//'="'//block_bc%F(f)%id2str()//' '//trim(str(n=block_bc%F(f)%adj%ID))//'"'
    elseif (block_bc%F(f)%is_inf()) then
      str_xml = trim(str_xml)//' '//att_name(f)//'="'//block_bc%F(f)%id2str()//' '//trim(str(n=block_bc%F(f)%inf  ))//'"'
    else
      str_xml = trim(str_xml)//' '//att_name(f)//'="'//block_bc%F(f)%id2str()//'"'
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_str_xml

  !> @brief Procedure for printing in a pretty ascii format the components of type Type_Block_BC.
  subroutine print_block_bc(block_bc,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_BC),   intent(IN)::  block_bc !< Block boundary conditions.
  character(*), optional, intent(IN)::  pref     !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat   !< IO error.
  character(*), optional, intent(OUT):: iomsg    !< IO error message.
  integer(I4P),           intent(IN)::  unit     !< Logic unit.
  character(len=:), allocatable::       prefd    !< Prefixing string for outputs.
  integer(I4P)::                        iostatd  !< IO error.
  character(500)::                      iomsgd   !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' imin='//block_bc%F(1)%id2str()//' imax='//block_bc%F(2)%id2str()
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' jmin='//block_bc%F(3)%id2str()//' jmax='//block_bc%F(4)%id2str()
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' kmin='//block_bc%F(5)%id2str()//' kmax='//block_bc%F(6)%id2str()
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_block_bc

  ! Assignment (=)
  !> @brief Procedure for assignment between two blocks extents variables.
  elemental subroutine assign_blk_bc(blk1,blk2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_BC), intent(INOUT):: blk1
  type(Type_Block_BC),  intent(IN)::    blk2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  blk1%F = blk2%F
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_blk_bc
  !> @}
endmodule Data_Type_Block_BC
