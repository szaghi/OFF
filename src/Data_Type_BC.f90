!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_BCDerivedType Data_Type_BC
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_BCGlobalVarPar Data_Type_BC
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_BCPublicProcedure Data_Type_BC
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_BCPrivateProcedure Data_Type_BC
!> @}

!> @brief This module contains the definition of Type_BC and its procedures.
!> Type_BC is a derived type containing all boundary conditions informations.
module Data_Type_BC
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: bc_nan_str,bc_nan
public:: bc_ref_str,bc_ref
public:: bc_ext_str,bc_ext
public:: bc_per_str,bc_per
public:: bc_adj_str,bc_adj
public:: bc_in1_str,bc_in1
public:: bc_in2_str,bc_in2
public:: Nbc
public:: bc_list,bc_list_str
public:: write_bc,read_bc
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Data_Type_BCGlobalVarPar
!> @{
character(3), parameter:: bc_nan_str = 'NAN' !< Definition of non-assigned boundary condition parameter string.
integer(I1P), parameter:: bc_nan     =-1_I1P !< Definition of non-assigned boundary condition parameter id.
character(3), parameter:: bc_ref_str = 'REF' !< Definition of reflective boundary condition parameter string.
integer(I1P), parameter:: bc_ref     = 1_I1P !< Definition of reflective boundary condition parameter id.
character(3), parameter:: bc_ext_str = 'EXT' !< Definition of extrapolation boundary condition parameter string.
integer(I1P), parameter:: bc_ext     = 2_I1P !< Definition of extrapolation boundary condition parameter id.
character(3), parameter:: bc_per_str = 'PER' !< Definition of periodic boundary condition parameter string.
integer(I1P), parameter:: bc_per     = 3_I1P !< Definition of periodic boundary condition parameter id.
character(3), parameter:: bc_adj_str = 'ADJ' !< Definition of adjacent boundary condition parameter string.
integer(I1P), parameter:: bc_adj     = 4_I1P !< Definition of adjacent boundary condition parameter id.
character(3), parameter:: bc_in1_str = 'IN1' !< Definition of inflow 1 boundary condition parameter string.
integer(I1P), parameter:: bc_in1     = 5_I1P !< Definition of inflow 1 boundary condition parameter id.
character(3), parameter:: bc_in2_str = 'IN2' !< Definition of inflow 2 boundary condition parameter string.
integer(I1P), parameter:: bc_in2     = 6_I1P !< Definition of inflow 2 boundary condition parameter id.
integer(I1P), parameter:: Nbc = 7_I1P        !< Number of possible boundary conditions.
character(3), parameter:: bc_list_str(1:Nbc) = &
                                             (/ bc_nan_str, &
                                                bc_ref_str, &
                                                bc_ext_str, &
                                                bc_per_str, &
                                                bc_adj_str, &
                                                bc_in1_str, &
                                                bc_in2_str  &
                                              /) !< Boundary conditions string list.
integer(I1P), parameter:: bc_list    (1:Nbc) = &
                                             (/ bc_nan, &
                                                bc_ref, &
                                                bc_ext, &
                                                bc_per, &
                                                bc_adj, &
                                                bc_in1, &
                                                bc_in2  &
                                              /) !< Boundary conditions list.
!> @}
!> Derived type containing adjacent boundary condition.
!> @ingroup Data_Type_BCDerivedType
type, public:: Type_Adj
  integer(I4P):: b = 0_I4P !< b index of adjacent block.
  integer(I4P):: i = 0_I4P !< i index of adjacent cell in the b block.
  integer(I4P):: j = 0_I4P !< j index of adjacent cell in the b block.
  integer(I4P):: k = 0_I4P !< k index of adjacent cell in the b block.
endtype Type_Adj
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
  integer(I1P)::                tp = bc_ext !< Type of boundary condition (bc_nan,bc_ref,bc_ext...).
  integer(I4P),   allocatable:: inf         !< Auxiliary informations for inflow-type boundary condition.
  type(Type_Adj), allocatable:: adj         !< Connection indexes for adjacent boundary condition.
  contains
    procedure:: init            ! Procedure for initilizing allocatable variables.
    procedure:: free => free_bc ! Procedure for freeing the memory of allocatable variables.
    procedure:: set             ! Procedure for setting bc members.
    procedure:: str2id          ! Procedure for setting integer id from string id.
    procedure:: id2str          ! Procedure for converting integer id to string id.
endtype Type_BC
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_BCPublicProcedure
  !> @{
  !> @brief Function for writing Type_BC data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function write_bc(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(IN), optional:: scalar         !< Scalar bc data.
  type(Type_BC), intent(IN), optional:: array1D(:)     !< One dimensional array bc data.
  type(Type_BC), intent(IN), optional:: array2D(:,:)   !< Two dimensional array bc data.
  type(Type_BC), intent(IN), optional:: array3D(:,:,:) !< Three dimensional array bc data.
  character(*),  intent(IN), optional:: format         !< Format specifier.
  integer(I4P),  intent(IN)::           unit           !< Logic unit.
  integer(I_P)::                        err            !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I4P)::                        i1,i2,i3       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(format)) then
    select case(adjustl(trim(format)))
    case('*')
      if (present(scalar)) then
        write(unit,*,iostat=err)scalar%tp
        select case(scalar%tp)
        case(bc_in1,bc_in2)
          write(unit,*,iostat=err)scalar%inf
        case(bc_adj)
          write(unit,*,iostat=err)scalar%adj
        endselect
      elseif (present(array1D)) then
        write(unit,*,iostat=err)array1D%tp
        do i1=lbound(array1D,dim=1),ubound(array1D,dim=1)
          select case(array1D(i1)%tp)
          case(bc_in1,bc_in2)
            write(unit,*,iostat=err)array1D(i1)%inf
          case(bc_adj)
            write(unit,*,iostat=err)array1D(i1)%adj
          endselect
        enddo
      elseif (present(array2D)) then
        write(unit,*,iostat=err)array2D%tp
        do i2=lbound(array2D,dim=2),ubound(array2D,dim=2)
          do i1=lbound(array2D,dim=1),ubound(array2D,dim=1)
            select case(array2D(i1,i2)%tp)
            case(bc_in1,bc_in2)
              write(unit,*,iostat=err)array2D(i1,i2)%inf
            case(bc_adj)
              write(unit,*,iostat=err)array2D(i1,i2)%adj
            endselect
          enddo
        enddo
      elseif (present(array3D)) then
        write(unit,*,iostat=err)array3D%tp
        do i3=lbound(array3D,dim=3),ubound(array3D,dim=3)
          do i2=lbound(array3D,dim=2),ubound(array3D,dim=2)
            do i1=lbound(array3D,dim=1),ubound(array3D,dim=1)
              select case(array3D(i1,i2,i3)%tp)
              case(bc_in1,bc_in2)
                write(unit,*,iostat=err)array3D(i1,i2,i3)%inf
              case(bc_adj)
                write(unit,*,iostat=err)array3D(i1,i2,i3)%adj
              endselect
            enddo
          enddo
        enddo
      endif
    case default
      if (present(scalar)) then
        write(unit,adjustl(trim(format)),iostat=err)scalar%tp
        select case(scalar%tp)
        case(bc_in1,bc_in2)
          write(unit,adjustl(trim(format)),iostat=err)scalar%inf
        case(bc_adj)
          write(unit,adjustl(trim(format)),iostat=err)scalar%adj
        endselect
      elseif (present(array1D)) then
        write(unit,adjustl(trim(format)),iostat=err)array1D%tp
        do i1=lbound(array1D,dim=1),ubound(array1D,dim=1)
          select case(array1D(i1)%tp)
          case(bc_in1,bc_in2)
            write(unit,adjustl(trim(format)),iostat=err)array1D(i1)%inf
          case(bc_adj)
            write(unit,adjustl(trim(format)),iostat=err)array1D(i1)%adj
          endselect
        enddo
      elseif (present(array2D)) then
        write(unit,adjustl(trim(format)),iostat=err)array2D%tp
        do i2=lbound(array2D,dim=2),ubound(array2D,dim=2)
          do i1=lbound(array2D,dim=1),ubound(array2D,dim=1)
            select case(array2D(i1,i2)%tp)
            case(bc_in1,bc_in2)
              write(unit,adjustl(trim(format)),iostat=err)array2D(i1,i2)%inf
            case(bc_adj)
              write(unit,adjustl(trim(format)),iostat=err)array2D(i1,i2)%adj
            endselect
          enddo
        enddo
      elseif (present(array3D)) then
        write(unit,adjustl(trim(format)),iostat=err)array3D%tp
        do i3=lbound(array3D,dim=3),ubound(array3D,dim=3)
          do i2=lbound(array3D,dim=2),ubound(array3D,dim=2)
            do i1=lbound(array3D,dim=1),ubound(array3D,dim=1)
              select case(array3D(i1,i2,i3)%tp)
              case(bc_in1,bc_in2)
                write(unit,adjustl(trim(format)),iostat=err)array3D(i1,i2,i3)%inf
              case(bc_adj)
                write(unit,adjustl(trim(format)),iostat=err)array3D(i1,i2,i3)%adj
              endselect
            enddo
          enddo
        enddo
      endif
    endselect
  else
    if (present(scalar)) then
      write(unit,iostat=err)scalar%tp
      select case(scalar%tp)
      case(bc_in1,bc_in2)
        write(unit,iostat=err)scalar%inf
      case(bc_adj)
        write(unit,iostat=err)scalar%adj
      endselect
    elseif (present(array1D)) then
      write(unit,iostat=err)array1D%tp
      do i1=lbound(array1D,dim=1),ubound(array1D,dim=1)
        select case(array1D(i1)%tp)
        case(bc_in1,bc_in2)
          write(unit,iostat=err)array1D(i1)%inf
        case(bc_adj)
          write(unit,iostat=err)array1D(i1)%adj
        endselect
      enddo
    elseif (present(array2D)) then
      write(unit,iostat=err)array2D%tp
      do i2=lbound(array2D,dim=2),ubound(array2D,dim=2)
        do i1=lbound(array2D,dim=1),ubound(array2D,dim=1)
          select case(array2D(i1,i2)%tp)
          case(bc_in1,bc_in2)
            write(unit,iostat=err)array2D(i1,i2)%inf
          case(bc_adj)
            write(unit,iostat=err)array2D(i1,i2)%adj
          endselect
        enddo
      enddo
    elseif (present(array3D)) then
      write(unit,iostat=err)array3D%tp
      do i3=lbound(array3D,dim=3),ubound(array3D,dim=3)
        do i2=lbound(array3D,dim=2),ubound(array3D,dim=2)
          do i1=lbound(array3D,dim=1),ubound(array3D,dim=1)
            select case(array3D(i1,i2,i3)%tp)
            case(bc_in1,bc_in2)
              write(unit,iostat=err)array3D(i1,i2,i3)%inf
            case(bc_adj)
              write(unit,iostat=err)array3D(i1,i2,i3)%adj
            endselect
          enddo
        enddo
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction write_bc

  !> @brief Function for reading Type_BC data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function read_bc(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(INOUT), optional:: scalar         !< Scalar bc data.
  type(Type_BC), intent(INOUT), optional:: array1D(:)     !< One dimensional array bc data.
  type(Type_BC), intent(INOUT), optional:: array2D(:,:)   !< Two dimensional array bc data.
  type(Type_BC), intent(INOUT), optional:: array3D(:,:,:) !< Three dimensional array bc data.
  character(*),  intent(IN),    optional:: format         !< Format specifier.
  integer(I4P),  intent(IN)::              unit           !< Logic unit.
  integer(I_P)::                           err            !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I4P)::                           i1,i2,i3       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(format)) then
    select case(adjustl(trim(format)))
    case('*')
      if (present(scalar)) then
        read(unit,*,iostat=err)scalar%tp ; call scalar%init
        select case(scalar%tp)
        case(bc_in1,bc_in2)
          read(unit,*,iostat=err)scalar%inf
        case(bc_adj)
          read(unit,*,iostat=err)scalar%adj
        endselect
      elseif (present(array1D)) then
        read(unit,*,iostat=err)array1D%tp ; call array1D%init
        do i1=lbound(array1D,dim=1),ubound(array1D,dim=1)
          select case(array1D(i1)%tp)
          case(bc_in1,bc_in2)
            read(unit,*,iostat=err)array1D(i1)%inf
          case(bc_adj)
            read(unit,*,iostat=err)array1D(i1)%adj
          endselect
        enddo
      elseif (present(array2D)) then
        read(unit,*,iostat=err)array2D%tp ; call array2D%init
        do i2=lbound(array2D,dim=2),ubound(array2D,dim=2)
          do i1=lbound(array2D,dim=1),ubound(array2D,dim=1)
            select case(array2D(i1,i2)%tp)
            case(bc_in1,bc_in2)
              read(unit,*,iostat=err)array2D(i1,i2)%inf
            case(bc_adj)
              read(unit,*,iostat=err)array2D(i1,i2)%adj
            endselect
          enddo
        enddo
      elseif (present(array3D)) then
        read(unit,*,iostat=err)array3D%tp ; call array3D%init
        do i3=lbound(array3D,dim=3),ubound(array3D,dim=3)
          do i2=lbound(array3D,dim=2),ubound(array3D,dim=2)
            do i1=lbound(array3D,dim=1),ubound(array3D,dim=1)
              select case(array3D(i1,i2,i3)%tp)
              case(bc_in1,bc_in2)
                read(unit,*,iostat=err)array3D(i1,i2,i3)%inf
              case(bc_adj)
                read(unit,*,iostat=err)array3D(i1,i2,i3)%adj
              endselect
            enddo
          enddo
        enddo
      endif
    case default
      if (present(scalar)) then
        read(unit,adjustl(trim(format)),iostat=err)scalar%tp ; call scalar%init
        select case(scalar%tp)
        case(bc_in1,bc_in2)
          read(unit,adjustl(trim(format)),iostat=err)scalar%inf
        case(bc_adj)
          read(unit,adjustl(trim(format)),iostat=err)scalar%adj
        endselect
      elseif (present(array1D)) then
        read(unit,adjustl(trim(format)),iostat=err)array1D%tp ; call array1D%init
        do i1=lbound(array1D,dim=1),ubound(array1D,dim=1)
          select case(array1D(i1)%tp)
          case(bc_in1,bc_in2)
            read(unit,adjustl(trim(format)),iostat=err)array1D(i1)%inf
          case(bc_adj)
            read(unit,adjustl(trim(format)),iostat=err)array1D(i1)%adj
          endselect
        enddo
      elseif (present(array2D)) then
        read(unit,adjustl(trim(format)),iostat=err)array2D%tp ; call array2D%init
        do i2=lbound(array2D,dim=2),ubound(array2D,dim=2)
          do i1=lbound(array2D,dim=1),ubound(array2D,dim=1)
            select case(array2D(i1,i2)%tp)
            case(bc_in1,bc_in2)
              read(unit,adjustl(trim(format)),iostat=err)array2D(i1,i2)%inf
            case(bc_adj)
              read(unit,adjustl(trim(format)),iostat=err)array2D(i1,i2)%adj
            endselect
          enddo
        enddo
      elseif (present(array3D)) then
        read(unit,adjustl(trim(format)),iostat=err)array3D%tp ; call array3D%init
        do i3=lbound(array3D,dim=3),ubound(array3D,dim=3)
          do i2=lbound(array3D,dim=2),ubound(array3D,dim=2)
            do i1=lbound(array3D,dim=1),ubound(array3D,dim=1)
              select case(array3D(i1,i2,i3)%tp)
              case(bc_in1,bc_in2)
                read(unit,adjustl(trim(format)),iostat=err)array3D(i1,i2,i3)%inf
              case(bc_adj)
                read(unit,adjustl(trim(format)),iostat=err)array3D(i1,i2,i3)%adj
              endselect
            enddo
          enddo
        enddo
      endif
    endselect
  else
    if (present(scalar)) then
      read(unit,iostat=err)scalar%tp ; call scalar%init
      select case(scalar%tp)
      case(bc_in1,bc_in2)
        read(unit,iostat=err)scalar%inf
      case(bc_adj)
        read(unit,iostat=err)scalar%adj
      endselect
    elseif (present(array1D)) then
      read(unit,iostat=err)array1D%tp ; call array1D%init
      do i1=lbound(array1D,dim=1),ubound(array1D,dim=1)
        select case(array1D(i1)%tp)
        case(bc_in1,bc_in2)
          read(unit,iostat=err)array1D(i1)%inf
        case(bc_adj)
          read(unit,iostat=err)array1D(i1)%adj
        endselect
      enddo
    elseif (present(array2D)) then
      read(unit,iostat=err)array2D%tp ; call array2D%init
      do i2=lbound(array2D,dim=2),ubound(array2D,dim=2)
        do i1=lbound(array2D,dim=1),ubound(array2D,dim=1)
          select case(array2D(i1,i2)%tp)
          case(bc_in1,bc_in2)
            read(unit,iostat=err)array2D(i1,i2)%inf
          case(bc_adj)
            read(unit,iostat=err)array2D(i1,i2)%adj
          endselect
        enddo
      enddo
    elseif (present(array3D)) then
      read(unit,iostat=err)array3D%tp ; call array3D%init
      do i3=lbound(array3D,dim=3),ubound(array3D,dim=3)
        do i2=lbound(array3D,dim=2),ubound(array3D,dim=2)
          do i1=lbound(array3D,dim=1),ubound(array3D,dim=1)
            select case(array3D(i1,i2,i3)%tp)
            case(bc_in1,bc_in2)
              read(unit,iostat=err)array3D(i1,i2,i3)%inf
            case(bc_adj)
              read(unit,iostat=err)array3D(i1,i2,i3)%adj
            endselect
          enddo
        enddo
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction read_bc
  !> @}

  !> @ingroup Data_Type_BCPrivateProcedure
  !> @{
  !> @brief Subroutine for initializing Type_BC allocatable variables.
  elemental subroutine init(bc,bc0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(INOUT)::        bc  !< Boundary conditions data.
  type(Type_BC),  intent(IN), optional:: bc0 !< Optional initialization data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(bc0)) then
    bc%tp = bc0%tp
    select case(bc%tp)
    case(bc_in1,bc_in2)
      if (.not.allocated(bc%inf)) allocate(bc%inf) ; bc%inf = bc0%inf
      if (     allocated(bc%adj)) deallocate(bc%adj)
    case(bc_adj)
      if (.not.allocated(bc%adj)) allocate(bc%adj) ; bc%adj = bc0%adj
      if (     allocated(bc%inf)) deallocate(bc%inf)
    endselect
  else
    select case(bc%tp)
    case(bc_in1,bc_in2)
      if (.not.allocated(bc%inf)) allocate(bc%inf) ; bc%inf = 0_I4P
      if (     allocated(bc%adj)) deallocate(bc%adj)
    case(bc_adj)
      if (.not.allocated(bc%adj)) allocate(bc%adj)
      if (     allocated(bc%inf)) deallocate(bc%inf)
    endselect
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  !> @brief Subroutine for freeing the memory of Type_BC allocatable variables.
  elemental subroutine free_bc(bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(INOUT):: bc  !< Boundary conditions data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
   if (allocated(bc%inf)) deallocate(bc%inf)
   if (allocated(bc%adj)) deallocate(bc%adj)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_bc

  !> @brief Subroutine for setting members of Type_BC variable.
  elemental subroutine set(bc,tp,inf,adj)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(INOUT)::        bc  !< Vector.
  integer(I1P),   intent(IN), optional:: tp  !< Type of boundary condition (bc_nan,bc_ref,bc_ext...).
  integer(I4P),   intent(IN), optional:: inf !< Auxiliary informations for inflow-type boundary condition.
  type(Type_Adj), intent(IN), optional:: adj !< Connection indexes for adjacent boundary condition.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(tp))  bc%tp  = tp
  if (present(inf)) bc%inf = inf
  if (present(adj)) bc%adj = adj
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  !> @brief Subroutine for setting integer id from string id.
  elemental subroutine str2id(bc,bc_str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC), intent(INOUT):: bc     !< BC data.
  character(3),   intent(IN)::    bc_str !< String of boundary condition.
  integer(I1P)::                  b      !< Boundary conditions counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do b=1_I1P,Nbc
    if (adjustl(trim(bc_str))==adjustl(trim(bc_list_str(b)))) then
      bc%tp = bc_list(b)
      exit
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine str2id

  !> @brief Function for converting integer id to string id.
  !> @return \b bc_str character(3) variable.
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
      bc_str = adjustl(trim(bc_list_str(b)))
      exit
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction id2str
  !> @}
endmodule Data_Type_BC
