!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Mesh_DimensionsDerivedType Data_Type_Mesh_Dimensions
!> Module definition of Type_Mesh_Dimensions
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Mesh_DimensionsPrivateProcedure Data_Type_Mesh_Dimensions
!> Module definition of Type_Mesh_Dimensions
!> @}

!> @brief Module Data_Type_Mesh_Dimensions contains the definition of Type_Mesh_Dimensions and useful procedures for its handling.
!> Type_Mesh_Dimensions contains all the data defining the main dimensions of a mesh.
module Data_Type_Mesh_Dimensions
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                            ! Integers and reals precision definition.
USE Data_Type_Block_Dimensions, only: Type_Block_Dimensions ! Definition of Type_Block_Dimensions.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the main dimensions of a mesh.
!> @ingroup Data_Type_Mesh_DimensionsDerivedType
type, public:: Type_Mesh_Dimensions
  integer(I4P)::                             Nl     = 1_I4P   !< Number of grid levels.
  integer(I4P)::                             Nb     = 0_I4P   !< Number of blocks of myrank.
  integer(I4P)::                             Nb_tot = 0_I4P   !< Number of total blocks (sum over each process).
  type(Type_Block_Dimensions), allocatable:: block_dims(:,:)  !< Blocks dimensions: a global view of all blocks [1:Nb_tot,1:Nl]
  logical::                                  set    = .false. !< Flag for inquiring if the dimensions are set up.
  contains
    procedure:: free  => free_mesh_dims           ! Procedure for freeing dynamic memory.
    procedure:: alloc => alloc_mesh_dims          ! Procedure for allocating dynamic memory.
    procedure:: load  => load_mesh_dims           ! Procedure for loading dimensions from mesh file.
    procedure:: save  => save_mesh_dims           ! Procedure for saving dimensions into mesh file.
    procedure:: print => print_mesh_dims          ! Procedure for printing dimensions with a pretty format.
    generic::   assignment(=) => assign_mesh_dims ! Assignment operator overloading.
    final::     finalize                          ! Procedure for allocating dynamic memory when finalizing.
    ! private procedures
    procedure:: assign_mesh_dims                  ! Procedure for assignment.
endtype Type_Mesh_Dimensions
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Mesh_DimensionsPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_mesh_dims(mesh_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Mesh_Dimensions), intent(INOUT):: mesh_dims !< Mesh dimensions data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(mesh_dims%block_dims)) deallocate(mesh_dims%block_dims)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_mesh_dims

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc_mesh_dims(mesh_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Mesh_Dimensions), intent(INOUT):: mesh_dims !< Mesh dimensions data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call mesh_dims%free ; allocate(mesh_dims%block_dims(1:mesh_dims%Nb_tot,1:mesh_dims%Nl))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_mesh_dims

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(mesh_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Mesh_Dimensions), intent(INOUT):: mesh_dims !< Mesh dimensions data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call mesh_dims%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for loading mesh dimensions from file.
  subroutine load_mesh_dims(mesh_dims,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Mesh_Dimensions), intent(INOUT):: mesh_dims !< Mesh dimensions.
  integer(I8P), optional,      intent(IN)::    pos       !< Position specifier.
  integer(I4P), optional,      intent(OUT)::   iostat    !< IO error.
  character(*), optional,      intent(OUT)::   iomsg     !< IO error message.
  integer(I4P),                intent(IN)::    unit      !< Logic unit.
  integer(I4P)::                               iostatd   !< IO error.
  character(500)::                             iomsgd    !< Temporary variable for IO error message.
  integer(I4P)::                               l,b       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    read(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)mesh_dims%Nl,mesh_dims%Nb,mesh_dims%Nb_tot
  else
    read(unit=unit,        iostat=iostatd,iomsg=iomsgd)mesh_dims%Nl,mesh_dims%Nb,mesh_dims%Nb_tot
  endif
  call mesh_dims%alloc
  do l=1,mesh_dims%Nl ; do b=1,mesh_dims%Nb_tot
    if (present(pos)) then
      call mesh_dims%block_dims(b,l)%load(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)
    else
      call mesh_dims%block_dims(b,l)%load(unit=unit,        iostat=iostatd,iomsg=iomsgd)
    endif
  enddo ; enddo
  mesh_dims%set = .true.
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = adjustl(trim(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_mesh_dims

  !> @brief Procedure for saving mesh dimensions into file.
  subroutine save_mesh_dims(mesh_dims,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Mesh_Dimensions), intent(IN)::  mesh_dims !< Mesh dimensions.
  integer(I8P), optional,      intent(IN)::  pos       !< Position specifier.
  integer(I4P), optional,      intent(OUT):: iostat    !< IO error.
  character(*), optional,      intent(OUT):: iomsg     !< IO error message.
  integer(I4P),                intent(IN)::  unit      !< Logic unit.
  integer(I4P)::                             iostatd   !< IO error.
  character(500)::                           iomsgd    !< Temporary variable for IO error message.
  integer(I4P)::                             l,b       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    write(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)mesh_dims%Nl,mesh_dims%Nb,mesh_dims%Nb_tot
  else
    write(unit=unit,        iostat=iostatd,iomsg=iomsgd)mesh_dims%Nl,mesh_dims%Nb,mesh_dims%Nb_tot
  endif
  do l=1,mesh_dims%Nl ; do b=1,mesh_dims%Nb_tot
    call mesh_dims%block_dims(b,l)%save(unit=unit,iostat=iostatd,iomsg=iomsgd)
  enddo ; enddo
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = adjustl(trim(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_mesh_dims

  !> @brief Procedure for printing dimensions to stdout.
  subroutine print_mesh_dims(mesh_dims,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Mesh_Dimensions), intent(IN)::  mesh_dims  !< Global-level data.
  character(*), optional,      intent(IN)::  pref       !< Prefixing string for outputs.
  integer(I4P), optional,      intent(OUT):: iostat     !< IO error.
  character(*), optional,      intent(OUT):: iomsg      !< IO error message.
  integer(I4P),                intent(IN)::  unit       !< Logic unit.
  character(len=:), allocatable::            prefd      !< Prefixing string for outputs.
  integer(I4P)::                             iostatd    !< IO error.
  character(500)::                           iomsgd     !< Temporary variable for IO error message.
  integer(I4P)::                             l,b        !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Number of grid levels '//trim(str(.true.,mesh_dims%Nl))
  write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Number of local  blocks (Nb) '//trim(str(.true.,mesh_dims%Nb))
  write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Number of global blocks (Nb_tot) '//trim(str(.true.,mesh_dims%Nb_tot))
  do l=1,mesh_dims%Nl
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Level '//trim(str(.true.,l))
    do b=1,mesh_dims%Nb_tot
      write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'     Block '//trim(str(.true.,b))
      call mesh_dims%block_dims(b,l)%print(pref=prefd//repeat(' ',7),iostat=iostatd,unit=unit,iomsg=iomsgd)
    enddo
  enddo
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_mesh_dims

  ! Assignment (=)
  !> @brief Procedure for assignment between two mesh dimensions variables.
  elemental subroutine assign_mesh_dims(mesh1,mesh2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Mesh_Dimensions), intent(INOUT):: mesh1
  type(Type_Mesh_Dimensions),  intent(IN)::    mesh2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mesh1%Nl         = mesh2%Nl
  mesh1%Nb         = mesh2%Nb
  mesh1%Nb_tot     = mesh2%Nb_tot
  mesh1%block_dims = mesh2%block_dims
  mesh1%set        = mesh2%set
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_mesh_dims
  !> @}
endmodule Data_Type_Mesh_Dimensions
