!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_BC Data_Type_BC
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
USE IR_Precision                                             ! Integers and reals precision definition.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stderr => ERROR_UNIT ! Standard output/error logical units.
#ifdef MPI2
USE MPI                                                      ! MPI runtime library.
#endif
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
public:: init,set,get,free
public:: write,read
public:: get_bc_id,get_bc_str
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Data_Type_BC
!> @{
character(3), parameter:: bc_nan_str = 'NAN' !< Definition of non-assigned boundary condition parameter string.
integer(I_P), parameter:: bc_nan     =-1_I_P !< Definition of non-assigned boundary condition parameter id.
character(3), parameter:: bc_ref_str = 'REF' !< Definition of reflective boundary condition parameter string.
integer(I_P), parameter:: bc_ref     = 1_I_P !< Definition of reflective boundary condition parameter id.
character(3), parameter:: bc_ext_str = 'EXT' !< Definition of extrapolation boundary condition parameter string.
integer(I_P), parameter:: bc_ext     = 2_I_P !< Definition of extrapolation boundary condition parameter id.
character(3), parameter:: bc_per_str = 'PER' !< Definition of periodic boundary condition parameter string.
integer(I_P), parameter:: bc_per     = 3_I_P !< Definition of periodic boundary condition parameter id.
character(3), parameter:: bc_adj_str = 'ADJ' !< Definition of adjacent boundary condition parameter string.
integer(I_P), parameter:: bc_adj     = 4_I_P !< Definition of adjacent boundary condition parameter id.
character(3), parameter:: bc_in1_str = 'IN1' !< Definition of inflow 1 boundary condition parameter string.
integer(I_P), parameter:: bc_in1     = 5_I_P !< Definition of inflow 1 boundary condition parameter id.
character(3), parameter:: bc_in2_str = 'IN2' !< Definition of inflow 2 boundary condition parameter string.
integer(I_P), parameter:: bc_in2     = 6_I_P !< Definition of inflow 2 boundary condition parameter id.
integer(I_P), parameter:: Nbc = 7            !< Number of possible boundary conditions.
character(3), parameter:: bc_list_str(1:Nbc) = &
                                             (/ bc_nan_str, &
                                                bc_ref_str, &
                                                bc_ext_str, &
                                                bc_per_str, &
                                                bc_adj_str, &
                                                bc_in1_str, &
                                                bc_in2_str  &
                                              /) !< Boundary conditions string list.
integer(I_P), parameter:: bc_list    (1:Nbc) = &
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
!> @ingroup DerivedType
type, public:: Type_Adj
  sequence
  integer(I_P):: b = 0_I_P !< b index of adjacent block.
  integer(I_P):: i = 0_I_P !< i index of adjacent block.
  integer(I_P):: j = 0_I_P !< j index of adjacent block.
  integer(I_P):: k = 0_I_P !< k index of adjacent block.
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
!> @ingroup DerivedType
type, public:: Type_BC
  integer(I_P)::                tp = bc_ext !< Type of boundary condition (bc_nan,bc_ref,bc_ext...).
  integer(I_P),   allocatable:: inf         !< Auxiliary informations for inflow-type  boundary condition.
  type(Type_Adj), allocatable:: adj         !< Connection indexes for adjacent boundary condition.
endtype Type_BC
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Function for freeing the memory of Type_BC \em dynamic components.
!> This is a generic interface to 4 functions as it can be used for scalar variables, 1D/2D or 3D arrays. The calling signatures
!> are:
!> @code ...
!> integer(I4P):: err
!> type(Type_BC):: bc_scal,bc_1D(10),bc_2D(10,2),bc_3D(10,2,3)
!> ...
!> ! freeing dynamic components memory of bc_scal, bc_1D, bc_2D and bc_3D
!> err = free(bc_scal)
!> err = free(bc_1D)
!> err = free(bc_2D)
!> err = free(bc_3D)
!> ... @endcode
!> @ingroup Interface,Data_Type_BCPublicProcedure
interface free
  module procedure Free_Scalar,Free_Array1D,Free_Array2D,Free_Array3D
endinterface
!> @brief Write overloading of Type_BC variable.
!> This is a generic interface to 8 functions: there are 2 functions (one binary and another ascii) for writing scalar variables,
!> 1D/2D or 3D arrays. The functions return an error integer code. The calling signatures are:
!> @code ...
!> integer(I4P):: err,unit
!> character(1):: format="*"
!> type(Type_BC):: bc_scal,bc_1D(10),bc_2D(10,2),bc_3D(10,2,3)
!> ...
!> ! formatted writing of bc_scal, bc_1D, bc_2D and bc_3D
!> err = write(unit,format,bc_scal)
!> err = write(unit,format,bc_1D)
!> err = write(unit,format,bc_2D)
!> err = write(unit,format,bc_3D)
!> ! binary writing of bc_scal, bc_1D, bc_2D and bc_3D
!> err = write(unit,bc_scal)
!> err = write(unit,bc_1D)
!> err = write(unit,bc_2D)
!> err = write(unit,bc_3D)
!> ... @endcode
!> @ingroup Interface,Data_Type_BCPublicProcedure
interface write
  module procedure Write_Bin_Scalar,     Write_Ascii_Scalar
  module procedure Write_Bin_Array1D,Write_Ascii_Array1D
  module procedure Write_Bin_Array2D,Write_Ascii_Array2D
  module procedure Write_Bin_Array3D,Write_Ascii_Array3D
endinterface
!> @brief Read overloading of Type_BC variable.
!> This is a generic interface to 8 functions: there are 2 functions (one binary and another ascii) for reading scalar variables,
!> 1D/2D or 3D arrays. The functions return an error integer code. The calling signatures are:
!> @code ...
!> integer(I4P):: err,unit
!> character(1):: format="*"
!> type(Type_BC):: bc_scal,bc_1D(10),bc_2D(10,2),bc_3D(10,2,3)
!> ...
!> ! formatted reading of bc_scal, bc_1D, bc_2D and bc_3D
!> err = read(unit,format,bc_scal)
!> err = read(unit,format,bc_1D)
!> err = read(unit,format,bc_2D)
!> err = read(unit,format,bc_3D)
!> ! binary reading of bc_scal, bc_1D, bc_2D and bc_3D
!> err = read(unit,bc_scal)
!> err = read(unit,bc_1D)
!> err = read(unit,bc_2D)
!> err = read(unit,bc_3D)
!> ... @endcode
!> @ingroup Interface,Data_Type_BCPublicProcedure
interface read
  module procedure Read_Bin_Scalar,     Read_Ascii_Scalar
  module procedure Read_Bin_Array1D,Read_Ascii_Array1D
  module procedure Read_Bin_Array2D,Read_Ascii_Array2D
  module procedure Read_Bin_Array3D,Read_Ascii_Array3D
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_BCPublicProcedure
  !> @{
  !> Function for initializing components of Type_BC variable.
  !> @return \b bc Type_BC variable.
  elemental function init(tp,inf,adj) result(bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),   intent(IN), optional:: tp  !< Type of boundary condition.
  integer(I_P),   intent(IN), optional:: inf !< Auxiliary informations for inflow-type  boundary condition.
  type(Type_Adj), intent(IN), optional:: adj !< Connection indexes for adjacent boundary condition.
  type(Type_BC)::                        bc  !< Boundary conditions data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(tp))  call set(tp  = tp,  bc = bc)
  if (present(inf)) call set(inf = inf, bc = bc)
  if (present(adj)) call set(adj = adj, bc = bc)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction init

  !> Subroutine for setting components of Type_BC variable.
  !> @return \b bc Type_BC variable.
  elemental subroutine set(tp,inf,adj,bc) !result(bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),   intent(IN), optional:: tp  !< Type of boundary condition.
  integer(I_P),   intent(IN), optional:: inf !< Auxiliary informations for inflow-type  boundary condition.
  type(Type_Adj), intent(IN), optional:: adj !< Connection indexes for adjacent boundary condition.
  type(Type_BC),  intent(INOUT)::        bc  !< Boundary conditions data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(tp))  bc%tp  = tp
  if (present(inf)) then
   if (allocated(bc%inf)) deallocate(bc%inf) ; allocate(bc%inf) ; bc%inf = inf
  endif
  if (present(adj)) then
   if (allocated(bc%adj)) deallocate(bc%adj) ; allocate(bc%adj) ; bc%adj = adj
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  !> Subroutine for extracting Type_BC variable components.
  elemental subroutine get(tp,inf,adj,bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),   intent(OUT), optional:: tp  !< Type of boundary condition.
  integer(I_P),   intent(OUT), optional:: inf !< Auxiliary informations for inflow-type  boundary condition.
  type(Type_Adj), intent(OUT), optional:: adj !< Connection indexes for adjacent boundary condition.
  type(Type_BC),  intent(IN)::            bc  !< Boundary conditions data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(tp )) tp  = bc%tp
  if (present(inf)) inf = bc%inf
  if (present(adj)) adj = bc%adj
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get
  !> @}

  !> @ingroup Data_Type_BCPrivateProcedure
  !> @{
  ! free dynamic memory
  !>Function for freeing the memory of Type_BC \em dynamic components (scalar).
  !> @return \b err integer(I4P) variable.
  function Free_Scalar(bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(INOUT):: bc  !< Boundary conditions data.
  integer(I4P)::                 err !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  if (allocated(bc%inf)) deallocate(bc%inf,stat=err)
  if (allocated(bc%adj)) deallocate(bc%adj,stat=err)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Scalar

  !>Function for freeing the memory of Type_BC \em dynamic components (array 1D).
  !> @return \b err integer(I4P) variable.
  function Free_Array1D(bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(INOUT):: bc(:) !< Boundary conditions data.
  integer(I4P)::                 err   !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I4P)::                 i     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do i=lbound(bc,dim=1),ubound(bc,dim=1)
    if (allocated(bc(i)%inf)) deallocate(bc(i)%inf,stat=err)
    if (allocated(bc(i)%adj)) deallocate(bc(i)%adj,stat=err)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array1D

  !>Function for freeing the memory of Type_BC \em dynamic components (array 2D).
  !> @return \b err integer(I4P) variable.
  function Free_Array2D(bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(INOUT):: bc(:,:) !< Boundary conditions data.
  integer(I4P)::                 err     !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I4P)::                 i,j     !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do j=lbound(bc,dim=2),ubound(bc,dim=2)
    do i=lbound(bc,dim=1),ubound(bc,dim=1)
      if (allocated(bc(i,j)%inf)) deallocate(bc(i,j)%inf,stat=err)
      if (allocated(bc(i,j)%adj)) deallocate(bc(i,j)%adj,stat=err)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array2D

  !>Function for freeing the memory of Type_BC \em dynamic components (array 3D).
  !> @return \b err integer(I4P) variable.
  function Free_Array3D(bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(INOUT):: bc(:,:,:) !< Boundary conditions data.
  integer(I4P)::                 err       !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I4P)::                 i,j,k     !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do k=lbound(bc,dim=3),ubound(bc,dim=3)
    do j=lbound(bc,dim=2),ubound(bc,dim=2)
      do i=lbound(bc,dim=1),ubound(bc,dim=1)
        if (allocated(bc(i,j,k)%inf)) deallocate(bc(i,j,k)%inf,stat=err)
        if (allocated(bc(i,j,k)%adj)) deallocate(bc(i,j,k)%adj,stat=err)
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array3D

  ! write
  !>Function for writing (binary, scalar) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Write_Bin_Scalar(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit  !< Logic unit.
  type(Type_bc), intent(IN):: bc    !< Boundary conditions data.
  integer(I_P)::              err   !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err) bc%tp
  select case(bc%tp)
  case(bc_in1,bc_in2)
    write(unit,iostat=err)bc%inf
  case(bc_adj)
    write(unit,iostat=err)bc%adj
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Scalar

  !>Function for writing (ascii, scalar) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Write_Ascii_Scalar(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit   !< Logic unit.
  character(*),  intent(IN):: format !< Format specifier.
  type(Type_bc), intent(IN):: bc     !< Boundary conditions data.
  integer(I_P)::              err    !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err) bc%tp
    select case(bc%tp)
    case(bc_in1,bc_in2)
      write(unit,*,iostat=err)bc%inf
    case(bc_adj)
      write(unit,*,iostat=err)bc%adj
    endselect
  case default
    write(unit,adjustl(trim(format)),iostat=err) bc%tp
    select case(bc%tp)
    case(bc_in1,bc_in2)
      write(unit,adjustl(trim(format)),iostat=err)bc%inf
    case(bc_adj)
      write(unit,adjustl(trim(format)),iostat=err)bc%adj
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Scalar

  !>Function for writing (binary, array 1D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Write_Bin_Array1D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit  !< Logic unit.
  type(Type_bc), intent(IN):: bc(:) !< Boundary conditions data.
  integer(I_P)::              err   !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::              n     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)bc%tp
  do n=lbound(bc,dim=1),ubound(bc,dim=1)
    select case(bc(n)%tp)
    case(bc_in1,bc_in2)
      write(unit,iostat=err)bc(n)%inf
    case(bc_adj)
      write(unit,iostat=err)bc(n)%adj
    endselect
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array1D

  !>Function for writing (ascii, array 1D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Write_Ascii_Array1D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit   !< Logic unit.
  character(*),  intent(IN):: format !< Format specifier.
  type(Type_bc), intent(IN):: bc(:)  !< Boundary conditions data.
  integer(I_P)::              err    !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::              n      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)bc%tp
    do n=lbound(bc,dim=1),ubound(bc,dim=1)
      select case(bc(n)%tp)
      case(bc_in1,bc_in2)
        write(unit,*,iostat=err)bc(n)%inf
      case(bc_adj)
        write(unit,*,iostat=err)bc(n)%adj
      endselect
    enddo
  case default
    write(unit,adjustl(trim(format)),iostat=err)bc%tp
    do n=lbound(bc,dim=1),ubound(bc,dim=1)
      select case(bc(n)%tp)
      case(bc_in1,bc_in2)
        write(unit,adjustl(trim(format)),iostat=err)bc(n)%inf
      case(bc_adj)
        write(unit,adjustl(trim(format)),iostat=err)bc(n)%adj
      endselect
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array1D

  !>Function for writing (binary, array 2D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Write_Bin_Array2D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit    !< Logic unit.
  type(Type_bc), intent(IN):: bc(:,:) !< Boundary conditions data.
  integer(I_P)::              err     !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::              n1,n2   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)bc%tp
  do n2=lbound(bc,dim=2),ubound(bc,dim=2)
    do n1=lbound(bc,dim=1),ubound(bc,dim=1)
      select case(bc(n1,n2)%tp)
      case(bc_in1,bc_in2)
        write(unit,iostat=err)bc(n1,n2)%inf
      case(bc_adj)
        write(unit,iostat=err)bc(n1,n2)%adj
      endselect
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array2D

  !>Function for writing (ascii, array 2D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Write_Ascii_Array2D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit    !< Logic unit.
  character(*),  intent(IN):: format  !< Format specifier.
  type(Type_bc), intent(IN):: bc(:,:) !< Boundary conditions data.
  integer(I_P)::              err     !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::              n1,n2   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)bc%tp
    do n2=lbound(bc,dim=2),ubound(bc,dim=2)
      do n1=lbound(bc,dim=1),ubound(bc,dim=1)
        select case(bc(n1,n2)%tp)
        case(bc_in1,bc_in2)
          write(unit,*,iostat=err)bc(n1,n2)%inf
        case(bc_adj)
          write(unit,*,iostat=err)bc(n1,n2)%adj
        endselect
      enddo
    enddo
  case default
    write(unit,adjustl(trim(format)),iostat=err)bc%tp
    do n2=lbound(bc,dim=2),ubound(bc,dim=2)
      do n1=lbound(bc,dim=1),ubound(bc,dim=1)
        select case(bc(n1,n2)%tp)
        case(bc_in1,bc_in2)
          write(unit,adjustl(trim(format)),iostat=err)bc(n1,n2)%inf
        case(bc_adj)
          write(unit,adjustl(trim(format)),iostat=err)bc(n1,n2)%adj
        endselect
      enddo
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array2D

  !>Function for writing (binary, array 3D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Write_Bin_Array3D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit      !< Logic unit.
  type(Type_bc), intent(IN):: bc(:,:,:) !< Boundary conditions data.
  integer(I_P)::              err       !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::              n1,n2,n3  !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)bc%tp
  do n3=lbound(bc,dim=3),ubound(bc,dim=3)
    do n2=lbound(bc,dim=2),ubound(bc,dim=2)
      do n1=lbound(bc,dim=1),ubound(bc,dim=1)
        select case(bc(n1,n2,n3)%tp)
        case(bc_in1,bc_in2)
          write(unit,iostat=err)bc(n1,n2,n3)%inf
        case(bc_adj)
          write(unit,iostat=err)bc(n1,n2,n3)%adj
        endselect
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array3D

  !>Function for writing (ascii, array 3D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Write_Ascii_Array3D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit      !< Logic unit.
  character(*),  intent(IN):: format    !< Format specifier.
  type(Type_bc), intent(IN):: bc(:,:,:) !< Boundary conditions data.
  integer(I_P)::              err       !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::              n1,n2,n3  !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)bc%tp
    do n3=lbound(bc,dim=3),ubound(bc,dim=3)
      do n2=lbound(bc,dim=2),ubound(bc,dim=2)
        do n1=lbound(bc,dim=1),ubound(bc,dim=1)
          select case(bc(n1,n2,n3)%tp)
          case(bc_in1,bc_in2)
            write(unit,*,iostat=err)bc(n1,n2,n3)%inf
          case(bc_adj)
            write(unit,*,iostat=err)bc(n1,n2,n3)%adj
          endselect
        enddo
      enddo
    enddo
  case default
    write(unit,adjustl(trim(format)),iostat=err)bc%tp
    do n3=lbound(bc,dim=3),ubound(bc,dim=3)
      do n2=lbound(bc,dim=2),ubound(bc,dim=2)
        do n1=lbound(bc,dim=1),ubound(bc,dim=1)
          select case(bc(n1,n2,n3)%tp)
          case(bc_in1,bc_in2)
            write(unit,adjustl(trim(format)),iostat=err)bc(n1,n2,n3)%inf
          case(bc_adj)
            write(unit,adjustl(trim(format)),iostat=err)bc(n1,n2,n3)%adj
          endselect
        enddo
      enddo
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array3D

  ! read
  !>Function for reading (binary, scalar) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Read_Bin_Scalar(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit  !< Logic unit.
  type(Type_bc), intent(INOUT):: bc    !< Boundary conditions data.
  integer(I_P)::                 err   !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)bc%tp
  select case(bc%tp)
  case(bc_in1,bc_in2)
    if (allocated(bc%inf)) deallocate(bc%inf) ; allocate(bc%inf) ; read(unit,iostat=err)bc%inf
  case(bc_adj)
    if (allocated(bc%adj)) deallocate(bc%adj) ; allocate(bc%adj) ; read(unit,iostat=err)bc%adj
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Scalar

  !>Function for reading (ascii, scalar) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Read_Ascii_Scalar(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit   !< Logic unit.
  character(*),  intent(IN)::    format !< Format specifier.
  type(Type_bc), intent(INOUT):: bc     !< Boundary conditions data.
  integer(I_P)::                 err    !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)bc%tp
    select case(bc%tp)
    case(bc_in1,bc_in2)
      if (allocated(bc%inf)) deallocate(bc%inf) ; allocate(bc%inf) ; read(unit,*,iostat=err)bc%inf
    case(bc_adj)
      if (allocated(bc%adj)) deallocate(bc%adj) ; allocate(bc%adj) ; read(unit,*,iostat=err)bc%adj
    endselect
  case default
    read(unit,adjustl(trim(format)),iostat=err) bc%tp
    select case(bc%tp)
    case(bc_in1,bc_in2)
      if (allocated(bc%inf)) deallocate(bc%inf) ; allocate(bc%inf) ; read(unit,adjustl(trim(format)),iostat=err)bc%inf
    case(bc_adj)
      if (allocated(bc%adj)) deallocate(bc%adj) ; allocate(bc%adj) ; read(unit,adjustl(trim(format)),iostat=err)bc%adj
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Scalar

  !>Function for reading (binary, array 1D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Read_Bin_Array1D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit  !< Logic unit.
  type(Type_bc), intent(INOUT):: bc(:) !< Boundary conditions data.
  integer(I_P)::                 err   !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                 n     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)bc%tp
  do n=lbound(bc,dim=1),ubound(bc,dim=1)
    select case(bc(n)%tp)
    case(bc_in1,bc_in2)
      if (allocated(bc(n)%inf)) deallocate(bc(n)%inf) ; allocate(bc(n)%inf) ; read(unit,iostat=err)bc(n)%inf
    case(bc_adj)
      if (allocated(bc(n)%adj)) deallocate(bc(n)%adj) ; allocate(bc(n)%adj) ; read(unit,iostat=err)bc(n)%adj
    endselect
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array1D

  !>Function for reading (ascii, array 1D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Read_Ascii_Array1D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit   !< Logic unit.
  character(*),  intent(IN)::    format !< Format specifier.
  type(Type_bc), intent(INOUT):: bc(:)  !< Boundary conditions data.
  integer(I_P)::                 err    !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                 n      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)bc%tp
    do n=lbound(bc,dim=1),ubound(bc,dim=1)
      select case(bc(n)%tp)
      case(bc_in1,bc_in2)
        if (allocated(bc(n)%inf)) deallocate(bc(n)%inf) ; allocate(bc(n)%inf) ; read(unit,*,iostat=err)bc(n)%inf
      case(bc_adj)
        if (allocated(bc(n)%adj)) deallocate(bc(n)%adj) ; allocate(bc(n)%adj) ; read(unit,*,iostat=err)bc(n)%adj
      endselect
    enddo
  case default
    read(unit,adjustl(trim(format)),iostat=err)bc%tp
    do n=lbound(bc,dim=1),ubound(bc,dim=1)
      select case(bc(n)%tp)
      case(bc_in1,bc_in2)
        if (allocated(bc(n)%inf)) deallocate(bc(n)%inf) ; allocate(bc(n)%inf)
        read(unit,adjustl(trim(format)),iostat=err)bc(n)%inf
      case(bc_adj)
        if (allocated(bc(n)%adj)) deallocate(bc(n)%adj) ; allocate(bc(n)%adj)
        read(unit,adjustl(trim(format)),iostat=err)bc(n)%adj
      endselect
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array1D

  !>Function for reading (binary, array 2D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Read_Bin_Array2D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit    !< logic unit.
  type(Type_bc), intent(INOUT):: bc(:,:) !< boundary conditions data.
  integer(I_P)::                 err     !< error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                 n1,n2   !< counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)bc%tp
  do n2=lbound(bc,dim=2),ubound(bc,dim=2)
    do n1=lbound(bc,dim=1),ubound(bc,dim=1)
      select case(bc(n1,n2)%tp)
      case(bc_in1,bc_in2)
        if (allocated(bc(n1,n2)%inf)) deallocate(bc(n1,n2)%inf) ; allocate(bc(n1,n2)%inf)
        read(unit,iostat=err)bc(n1,n2)%inf
      case(bc_adj)
        if (allocated(bc(n1,n2)%adj)) deallocate(bc(n1,n2)%adj) ; allocate(bc(n1,n2)%adj)
        read(unit,iostat=err)bc(n1,n2)%adj
      endselect
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array2D

  !>Function for reading (ascii, array 2D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Read_Ascii_Array2D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit    !< Logic unit.
  character(*),  intent(IN)::    format  !< Format specifier.
  type(Type_bc), intent(INOUT):: bc(:,:) !< Boundary conditions data.
  integer(I_P)::                 err     !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                 n1,n2   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)bc%tp
    do n2=lbound(bc,dim=2),ubound(bc,dim=2)
      do n1=lbound(bc,dim=1),ubound(bc,dim=1)
        select case(bc(n1,n2)%tp)
        case(bc_in1,bc_in2)
          if (allocated(bc(n1,n2)%inf)) deallocate(bc(n1,n2)%inf) ; allocate(bc(n1,n2)%inf)
          read(unit,*,iostat=err)bc(n1,n2)%inf
        case(bc_adj)
          if (allocated(bc(n1,n2)%adj)) deallocate(bc(n1,n2)%adj) ; allocate(bc(n1,n2)%adj)
          read(unit,*,iostat=err)bc(n1,n2)%adj
        endselect
      enddo
    enddo
  case default
    read(unit,adjustl(trim(format)),iostat=err)bc%tp
    do n2=lbound(bc,dim=2),ubound(bc,dim=2)
      do n1=lbound(bc,dim=1),ubound(bc,dim=1)
        select case(bc(n1,n2)%tp)
        case(bc_in1,bc_in2)
          if (allocated(bc(n1,n2)%inf)) deallocate(bc(n1,n2)%inf) ; allocate(bc(n1,n2)%inf)
          read(unit,adjustl(trim(format)),iostat=err)bc(n1,n2)%inf
        case(bc_adj)
          if (allocated(bc(n1,n2)%adj)) deallocate(bc(n1,n2)%adj) ; allocate(bc(n1,n2)%adj)
          read(unit,adjustl(trim(format)),iostat=err)bc(n1,n2)%adj
        endselect
      enddo
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array2D

  !>Function for reading (binary, array 3D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Read_Bin_Array3D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit      !< logic unit.
  type(Type_bc), intent(INOUT):: bc(:,:,:) !< boundary conditions data.
  integer(I_P)::                 err       !< error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                 n1,n2,n3  !< counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)bc%tp
  do n3=lbound(bc,dim=3),ubound(bc,dim=3)
    do n2=lbound(bc,dim=2),ubound(bc,dim=2)
      do n1=lbound(bc,dim=1),ubound(bc,dim=1)
        select case(bc(n1,n2,n3)%tp)
        case(bc_in1,bc_in2)
          if (allocated(bc(n1,n2,n3)%inf)) deallocate(bc(n1,n2,n3)%inf) ; allocate(bc(n1,n2,n3)%inf)
          read(unit,iostat=err)bc(n1,n2,n3)%inf
        case(bc_adj)
          if (allocated(bc(n1,n2,n3)%adj)) deallocate(bc(n1,n2,n3)%adj) ; allocate(bc(n1,n2,n3)%adj)
          read(unit,iostat=err)bc(n1,n2,n3)%adj
        endselect
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array3D

  !>Function for reading (ascii, array 3D) Type_BC variable.
  !> @return \b err integer(I4P) variable.
  function Read_Ascii_Array3D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit      !< Logic unit.
  character(*),  intent(IN)::    format    !< Format specifier.
  type(Type_bc), intent(INOUT):: bc(:,:,:) !< Boundary conditions data.
  integer(I_P)::                 err       !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                 n1,n2,n3  !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)bc%tp
    do n3=lbound(bc,dim=3),ubound(bc,dim=3)
      do n2=lbound(bc,dim=2),ubound(bc,dim=2)
        do n1=lbound(bc,dim=1),ubound(bc,dim=1)
          select case(bc(n1,n2,n3)%tp)
          case(bc_in1,bc_in2)
            if (allocated(bc(n1,n2,n3)%inf)) deallocate(bc(n1,n2,n3)%inf) ; allocate(bc(n1,n2,n3)%inf)
            read(unit,*,iostat=err)bc(n1,n2,n3)%inf
          case(bc_adj)
            if (allocated(bc(n1,n2,n3)%adj)) deallocate(bc(n1,n2,n3)%adj) ; allocate(bc(n1,n2,n3)%adj)
            read(unit,*,iostat=err)bc(n1,n2,n3)%adj
          endselect
        enddo
      enddo
    enddo
  case default
    read(unit,adjustl(trim(format)),iostat=err)bc%tp
    do n3=lbound(bc,dim=3),ubound(bc,dim=3)
      do n2=lbound(bc,dim=2),ubound(bc,dim=2)
        do n1=lbound(bc,dim=1),ubound(bc,dim=1)
          select case(bc(n1,n2,n3)%tp)
          case(bc_in1,bc_in2)
            if (allocated(bc(n1,n2,n3)%inf)) deallocate(bc(n1,n2,n3)%inf) ; allocate(bc(n1,n2,n3)%inf)
            read(unit,adjustl(trim(format)),iostat=err)bc(n1,n2,n3)%inf
          case(bc_adj)
            if (allocated(bc(n1,n2,n3)%adj)) deallocate(bc(n1,n2,n3)%adj) ; allocate(bc(n1,n2,n3)%adj)
            read(unit,adjustl(trim(format)),iostat=err)bc(n1,n2,n3)%adj
          endselect
        enddo
      enddo
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array3D
  !> @}

  !> @ingroup Data_Type_BCPublicProcedure
  !> @{
  !> Function for getting integer id of a boundary condition from the corresponding boundary condition string.
  !> @return \b bc_id integer(I_P) variable.
  function get_bc_id(myrank,bc_str) result(bc_id)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: myrank !< Actual rank process.
  character(*), intent(IN):: bc_str !< String id of boundary condition.
  integer(I_P)::             bc_id  !< Integer id of boundary condition.
  integer(I_P)::             b      !< Boundary conditions counter.
#ifdef MPI2
  integer(I_P)::             err    !< Error for MPI communications.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do b=1,Nbc
    if (adjustl(trim(bc_str))==adjustl(trim(bc_list_str(b)))) then
      bc_id = bc_list(b)
      exit
    elseif (b==Nbc) then
      write(stderr,'(A,I3)')' My RANK is: ',myrank
      write(stderr,'(A)')   ' Attention!'
      write(stderr,'(A)')   ' The boundary condition:'
      write(stderr,'(A)')   ' '//adjustl(trim(bc_str))
      write(stderr,'(A)')   ' is unknown!'
#ifdef MPI2
      call MPI_FINALIZE(err)
#endif
      stop
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction get_bc_id

  !> Function for getting string id of a boundary condition from the corresponding boundary condition integer id.
  !> @return \b bc_str character(3) variable.
  function get_bc_str(myrank,bc_id) result(bc_str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: myrank !< Actual rank process.
  integer(I_P), intent(IN):: bc_id  !< Integer id of boundary condition.
  character(3)::             bc_str !< String of boundary condition.
  integer(I_P)::             b      !< Boundary conditions counter.
#ifdef MPI2
  integer(I_P)::             err    !< Error for MPI communications.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do b=1,Nbc
    if (bc_id==bc_list(b)) then
      bc_str = adjustl(trim(bc_list_str(b)))
      exit
    elseif (b==Nbc) then
      write(stderr,'(A,I3)')' My RANK is: ',myrank
      write(stderr,'(A)')   ' Attention!'
      write(stderr,'(A)')   ' The boundary condition:'
      write(stderr,FI_P)    bc_id
      write(stderr,'(A)')   ' is unknown!'
#ifdef MPI2
      call MPI_FINALIZE(err)
#endif
      stop
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction get_bc_str
  !> @}
endmodule Data_Type_BC
