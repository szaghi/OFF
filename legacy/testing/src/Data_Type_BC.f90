!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_BCDerivedType Data_Type_BC
!> Module definition of Type_BC
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_BCGlobalVarPar Data_Type_BC
!> Module definition of Type_BC
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_BCInterface Data_Type_BC
!> Module definition of Type_BC
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_BCPrivateProcedure Data_Type_BC
!> Module definition of Type_BC
!> @}

!> @brief This module contains the definition of Type_BC and its procedures.
!> Type_BC is a derived type containing all boundary conditions informations.
module Data_Type_BC
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                    ! Integers and reals precision definition.
USE Data_Type_Cell_Indexes, only: Type_Cell_Indexes ! Definition of Type_Cell_Indexes.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Data_Type_BCGlobalVarPar
!> @{
character(3), parameter, public:: bc_nan_str = 'NAN' !< Definition of non-assigned boundary condition parameter string.
integer(I1P), parameter, public:: bc_nan     =-1_I1P !< Definition of non-assigned boundary condition parameter id.
character(3), parameter, public:: bc_ref_str = 'REF' !< Definition of reflective boundary condition parameter string.
integer(I1P), parameter, public:: bc_ref     = 1_I1P !< Definition of reflective boundary condition parameter id.
character(3), parameter, public:: bc_ext_str = 'EXT' !< Definition of extrapolation boundary condition parameter string.
integer(I1P), parameter, public:: bc_ext     = 2_I1P !< Definition of extrapolation boundary condition parameter id.
character(3), parameter, public:: bc_per_str = 'PER' !< Definition of periodic boundary condition parameter string.
integer(I1P), parameter, public:: bc_per     = 3_I1P !< Definition of periodic boundary condition parameter id.
character(3), parameter, public:: bc_adj_str = 'ADJ' !< Definition of adjacent boundary condition parameter string.
integer(I1P), parameter, public:: bc_adj     = 4_I1P !< Definition of adjacent boundary condition parameter id.
character(3), parameter, public:: bc_in1_str = 'IN1' !< Definition of inflow 1 boundary condition parameter string.
integer(I1P), parameter, public:: bc_in1     = 5_I1P !< Definition of inflow 1 boundary condition parameter id.
character(3), parameter, public:: bc_in2_str = 'IN2' !< Definition of inflow 2 boundary condition parameter string.
integer(I1P), parameter, public:: bc_in2     = 6_I1P !< Definition of inflow 2 boundary condition parameter id.
integer(I1P), parameter, public:: Nbc = 7_I1P        !< Number of possible boundary conditions.
character(3), parameter, public:: bc_list_str(1:Nbc) = &
                                                     (/ bc_nan_str, &
                                                        bc_ref_str, &
                                                        bc_ext_str, &
                                                        bc_per_str, &
                                                        bc_adj_str, &
                                                        bc_in1_str, &
                                                        bc_in2_str  &
                                                      /) !< Boundary conditions string list.
integer(I1P), parameter, public:: bc_list    (1:Nbc) = &
                                                     (/ bc_nan, &
                                                        bc_ref, &
                                                        bc_ext, &
                                                        bc_per, &
                                                        bc_adj, &
                                                        bc_in1, &
                                                        bc_in2  &
                                                      /) !< Boundary conditions list.
!> @}
!> Derived type containing boundary conditions informations.
!> @note
!>   meaning of \b inf component \n
!>   - \b bc_in1 type condition: supersonic inflow steady conditions \n
!>     the component \b inf is the index array of inflow boundary conditions; inflow 1 conditions (primitive variables) are stored
!>     in array "in1" so the boundary conditions can be accessed by: @code in1(1:Np,inf) @endcode
!>   - \b bc_in2 type condition: supersonic inflow unsteady conditions \n
!>     the component \b inf is the index array of inflow boundary conditions; inflow 2 conditions (primitive variables) are stored
!>     in array "in2" so the boundary conditions can be accessed by: @code in2(1:Np,n,inf) @endcode
!>     where n is the time step counter.
!> @ingroup Data_Type_BCDerivedType
type, public:: Type_BC
  integer(I1P)::                         tp = bc_ext !< Type of boundary condition (bc_nan,bc_ref,bc_ext...).
  integer(I4P),            allocatable:: inf         !< Auxiliary informations for inflow-type boundary condition.
  type(Type_Cell_Indexes), allocatable:: adj         !< Connection indexes for adjacent boundary condition.
  contains
    procedure:: free                    ! Procedure for freeing dynamic memory.
    procedure:: alloc                   ! Procedure for allocating dynamic memory.
    procedure:: set                     ! Procedure for setting bc members.
    procedure:: load   => load_bc_self  ! Procedure for loading boundary conditions.
    procedure:: save   => save_bc_self  ! Procedure for saving boundary conditions.
    procedure:: print  => print_bc_self ! Procedure for printing boundary conditions with a pretty format.
    procedure:: str2id                  ! Procedure for setting integer id from string id.
    procedure:: id2str                  ! Procedure for converting integer id to string id.
    procedure:: is_ref                  ! Procedure for checking if bc is reflective.
    procedure:: is_ext                  ! Procedure for checking if bc is extrapolation.
    procedure:: is_per                  ! Procedure for checking if bc is periodic.
    procedure:: is_adj                  ! Procedure for checking if bc is adjacent.
    procedure:: is_in1                  ! Procedure for checking if bc is inflow 1.
    procedure:: is_in2                  ! Procedure for checking if bc is inflow 2.
    procedure:: is_inf                  ! Procedure for checking if bc is inflow.
    final::     finalize                ! Procedure for freeing dynamic memory when finalizing.
    ! operators overloading
    generic:: assignment(=) => assign_bc
    ! private procedures
    procedure, pass(bc1), private:: assign_bc
endtype Type_BC
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_BCPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free(bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(INOUT):: bc !< Boundary conditions data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(bc%inf)) deallocate(bc%inf)
  if (allocated(bc%adj)) deallocate(bc%adj)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine finalize(bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(INOUT):: bc !< Boundary conditions data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call bc%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc(bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(INOUT):: bc !< Boundary conditions data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(bc%tp)
  case(bc_in1,bc_in2)
    if (.not.allocated(bc%inf))   allocate(bc%inf) ; bc%inf = 0_I4P
    if (     allocated(bc%adj)) deallocate(bc%adj)
  case(bc_adj)
    if (.not.allocated(bc%adj))   allocate(bc%adj)
    if (     allocated(bc%inf)) deallocate(bc%inf)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc

  !> @brief Procedure for setting members of Type_BC variable.
  elemental subroutine set(bc,tp,inf,adj)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC),          intent(INOUT)::        bc  !< Vector.
  integer(I1P),            intent(IN), optional:: tp  !< Type of boundary condition (bc_nan,bc_ref,bc_ext...).
  integer(I4P),            intent(IN), optional:: inf !< Auxiliary informations for inflow-type boundary condition.
  type(Type_Cell_Indexes), intent(IN), optional:: adj !< Connection indexes for adjacent boundary condition.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(tp))  bc%tp  = tp
  if (present(inf)) bc%inf = inf
  if (present(adj)) bc%adj = adj
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  !> @brief Procedure for loading boundary conditions.
  subroutine load_bc_self(bc,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC),         intent(INOUT):: bc      !< Boundary conditions data.
  integer(I8P), optional, intent(IN)::    pos     !< Position specifier.
  integer(I4P), optional, intent(OUT)::   iostat  !< IO error.
  character(*), optional, intent(OUT)::   iomsg   !< IO error message.
  integer(I4P),           intent(IN)::    unit    !< Logic unit.
  integer(I4P)::                          iostatd !< IO error.
  character(500)::                        iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    read(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)bc%tp
  else
    read(unit=unit,        iostat=iostatd,iomsg=iomsgd)bc%tp
  endif
  call bc%alloc
  select case(bc%tp)
  case(bc_in1,bc_in2)
    read(unit=unit,        iostat=iostatd,iomsg=iomsgd)bc%inf
  case(bc_adj)
    read(unit=unit,        iostat=iostatd,iomsg=iomsgd)bc%adj
  endselect
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = adjustl(trim(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_bc_self

  !> @brief Procedure for saving boundary conditions.
  subroutine save_bc_self(bc,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC),         intent(IN)::  bc      !< Boundary conditions data.
  integer(I8P), optional, intent(IN)::  pos     !< Position specifier.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    write(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)bc%tp
  else
    write(unit=unit,        iostat=iostatd,iomsg=iomsgd)bc%tp
  endif
  select case(bc%tp)
  case(bc_in1,bc_in2)
    write(unit=unit,        iostat=iostatd,iomsg=iomsgd)bc%inf
  case(bc_adj)
    write(unit=unit,        iostat=iostatd,iomsg=iomsgd)bc%adj
  endselect
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = adjustl(trim(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_bc_self

  !> @brief Procedure for printing boundary conditions with a pretty format.
  subroutine print_bc_self(bc,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC),         intent(IN)::  bc      !< Boundary conditions data.
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
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' type='//bc%id2str()
  select case(bc%tp)
  case(bc_adj)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   adjacent cell(b,i,j,k)='&
      //trim(str(.true.,bc%adj%ID))//','//&
      trim(str(.true.,bc%adj%i))//','//trim(str(.true.,bc%adj%j))//','//trim(str(.true.,bc%adj%k))
  case(bc_in1,bc_in2)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   inflow='//trim(str(.true.,bc%inf))
  endselect
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_bc_self

  !> @brief Procedure for setting integer id from string id.
  elemental subroutine str2id(bc,bc_str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(INOUT):: bc     !< BC data.
  character(3),   intent(IN)::    bc_str !< String of boundary condition.
  integer(I1P)::                  b      !< Boundary conditions counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do b=1_I1P,Nbc
    if (trim(adjustl(bc_str))==trim(adjustl(bc_list_str(b)))) then
      bc%tp = bc_list(b)
      exit
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine str2id

  !> @brief Procedure for converting integer id to string id.
  elemental function id2str(bc) result(bc_str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(IN):: bc     !< BC data.
  character(3)::               bc_str !< String of boundary condition.
  integer(I1P)::               b      !< Boundary conditions counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do b=1_I1P,Nbc
    if (bc%tp==bc_list(b)) then
      bc_str = trim(adjustl(bc_list_str(b)))
      exit
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction id2str

  !> @brief Procedure for checking if bc is reflective.
  elemental function is_ref(bc) result(it_is)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(IN):: bc     !< BC data.
  logical::                    it_is  !< Logical flag for checking if bc is of type reflective.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  it_is = (bc%tp==bc_ref)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_ref

  !> @brief Procedure for checking if bc is extrapolation.
  elemental function is_ext(bc) result(it_is)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(IN):: bc     !< BC data.
  logical::                    it_is  !< Logical flag for checking if bc is of type extrapolation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  it_is = (bc%tp==bc_ext)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_ext

  !> @brief Procedure for checking if bc is periodic.
  elemental function is_per(bc) result(it_is)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(IN):: bc     !< BC data.
  logical::                    it_is  !< Logical flag for checking if bc is of type periodic.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  it_is = (bc%tp==bc_per)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_per

  !> @brief Procedure for checking if bc is adjacent.
  elemental function is_adj(bc) result(it_is)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(IN):: bc     !< BC data.
  logical::                    it_is  !< Logical flag for checking if bc is of type adjacent.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  it_is = (bc%tp==bc_adj)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_adj

  !> @brief Procedure for checking if bc is inflow 1.
  elemental function is_in1(bc) result(it_is)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(IN):: bc     !< BC data.
  logical::                    it_is  !< Logical flag for checking if bc is of type inflow.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  it_is = (bc%tp==bc_in1)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_in1

  !> @brief Procedure for checking if bc is inflow 2.
  elemental function is_in2(bc) result(it_is)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(IN):: bc     !< BC data.
  logical::                    it_is  !< Logical flag for checking if bc is of type inflow.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  it_is = (bc%tp==bc_in2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_in2

  !> @brief Procedure for checking if bc is inflow.
  elemental function is_inf(bc) result(it_is)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(IN):: bc     !< BC data.
  logical::                    it_is  !< Logical flag for checking if bc is of type inflow.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  it_is = ((bc%is_in1()).or.(bc%is_in2()))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_inf

  ! Assignment (=)
  !> @brief Procedure for assignment between two boundary conditions variables.
  elemental subroutine assign_bc(bc1,bc2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(INOUT):: bc1 !< LHS.
  type(Type_BC),  intent(IN)::    bc2 !< RHS.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bc1%tp = bc2%tp
  select case(bc1%tp)
  case(bc_in1,bc_in2)
    if (.not.allocated(bc1%inf))   allocate(bc1%inf) ; bc1%inf = bc2%inf
    if (     allocated(bc1%adj)) deallocate(bc1%adj)
  case(bc_adj)
    if (.not.allocated(bc1%adj))   allocate(bc1%adj) ; bc1%adj = bc2%adj
    if (     allocated(bc1%inf)) deallocate(bc1%inf)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_bc
  !> @}
endmodule Data_Type_BC
