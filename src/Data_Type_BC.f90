module Data_Type_BC
!-----------------------------------------------------------------------------------------------------------------------------------
! The module Data_Type_BC contains the definition of Type_BC and its functions and subroutines.
!-----------------------------------------------------------------------------------------------------------------------------------

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
public:: Type_Adj
public:: Type_BC
public:: init,set,get,free
public:: write,read
public:: get_bc_id,get_bc_str
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Type id paramters
character(3), parameter:: bc_nan_str = 'NAN' ! Definition of non-assigned boundary condition parameter string.
integer(I_P), parameter:: bc_nan     =-1_I_P ! Definition of non-assigned boundary condition parameter id.
character(3), parameter:: bc_ref_str = 'REF' ! Definition of reflective boundary condition parameter string.
integer(I_P), parameter:: bc_ref     = 1_I_P ! Definition of reflective boundary condition parameter id.
character(3), parameter:: bc_ext_str = 'EXT' ! Definition of extrapolation boundary condition parameter string.
integer(I_P), parameter:: bc_ext     = 2_I_P ! Definition of extrapolation boundary condition parameter id.
character(3), parameter:: bc_per_str = 'PER' ! Definition of periodic boundary condition parameter string.
integer(I_P), parameter:: bc_per     = 3_I_P ! Definition of periodic boundary condition parameter id.
character(3), parameter:: bc_adj_str = 'ADJ' ! Definition of adjacent boundary condition parameter string.
integer(I_P), parameter:: bc_adj     = 4_I_P ! Definition of adjacent boundary condition parameter id.
character(3), parameter:: bc_in1_str = 'IN1' ! Definition of inflow 1 boundary condition parameter string.
integer(I_P), parameter:: bc_in1     = 5_I_P ! Definition of inflow 1 boundary condition parameter id.
character(3), parameter:: bc_in2_str = 'IN2' ! Definition of inflow 2 boundary condition parameter string.
integer(I_P), parameter:: bc_in2     = 6_I_P ! Definition of inflow 2 boundary condition parameter id.
integer(I_P), parameter:: Nbc = 7            ! Number of possible boundary conditions.
character(3), parameter:: bc_list_str(1:Nbc) = (/ bc_nan_str, &
                                                  bc_ref_str, &
                                                  bc_ext_str, &
                                                  bc_per_str, &
                                                  bc_adj_str, &
                                                  bc_in1_str, &
                                                  bc_in2_str /) ! Boundary conditions string list.
integer(I_P), parameter:: bc_list    (1:Nbc) = (/ bc_nan, &
                                                  bc_ref, &
                                                  bc_ext, &
                                                  bc_per, &
                                                  bc_adj, &
                                                  bc_in1, &
                                                  bc_in2 /) ! Boundary conditions list.
! Definition of adjacent boundary condition
type:: Type_Adj
  sequence
  integer(I_P):: b = 0_I_P ! b index of connected block.
  integer(I_P):: i = 0_I_P ! i index of connected block.
  integer(I_P):: j = 0_I_P ! j index of connected block.
  integer(I_P):: k = 0_I_P ! k index of connected block.
endtype Type_Adj
! Definition of Type_BC
type:: Type_BC
  integer(I_P)::                tp = bc_ext ! Type of boundary condition (bc_nan,bc_ref,bc_ext...).
  integer(I_P),   allocatable:: inf         ! Auxiliary infos for inlfow-type  boundary condition:
                                            ! 1)"tp = bc_in1": supersonic inflow steady conditions
                                            !   it is the index array of inflow boundary conditions; inflow 1 conditions (primitive
                                            !   variables) are stored in array "in1" so the boundary conditions can be accessed by:
                                            !   in1(1:Np,inf);
                                            ! 2)"tp = bc_in2": supersonic inflow unsteady conditions
                                            !   it is the index array of inflow boundary conditions; inflow 2 conditions (primitive
                                            !   variables) are stored in array "in2" so the boundary conditions can be accessed by:
                                            !   in2(1:Np,n,inf)
  type(Type_Adj), allocatable:: adj         ! Connection indexes for adjacent boundary condition.
endtype Type_BC
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!!Free overloading.
interface free
  module procedure Free_Scalar,  & ! scalar
                   Free_Array1D, & ! array1D
                   Free_Array2D, & ! array2D
                   Free_Array3D    ! array3D
endinterface
!!Write overloading.
interface write
  module procedure Write_Bin_Scalar,        & ! binary scalar
                   Write_Ascii_Scalar,      & ! ascii scalar
                   Write_Bin_Vectorial1D,   & ! binary vectorial 1D
                   Write_Ascii_Vectorial1D, & ! ascii vectorial 1D
                   Write_Bin_Vectorial2D,   & ! binary vectorial 2D
                   Write_Ascii_Vectorial2D, & ! ascii vectorial 2D
                   Write_Bin_Vectorial3D,   & ! binary vectorial 3D
                   Write_Ascii_Vectorial3D    ! ascii vectorial 3D
endinterface
!!Read overloading.
interface read
  module procedure Read_Bin_Scalar,        & ! binary scalar
                   Read_Ascii_Scalar,      & ! ascii scalar
                   Read_Bin_Vectorial1D,   & ! binary vectorial 1D
                   Read_Ascii_Vectorial1D, & ! ascii vectorial 1D
                   Read_Bin_Vectorial2D,   & ! binary vectorial 2D
                   Read_Ascii_Vectorial2D, & ! ascii vectorial 2D
                   Read_Bin_Vectorial3D,   & ! binary vectorial 3D
                   Read_Ascii_Vectorial3D    ! ascii vectorial 3D
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  elemental function init(tp,inf,adj) result(bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for initializing Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),   intent(IN), optional:: tp
  integer(I_P),   intent(IN), optional:: inf
  type(Type_Adj), intent(IN), optional:: adj
  type(Type_BC)::                        bc
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(tp))  call set(tp  = tp,  bc = bc)
  if (present(inf)) call set(inf = inf, bc = bc)
  if (present(adj)) call set(adj = adj, bc = bc)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction init

  elemental subroutine set(tp,inf,adj,bc) !result(bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for setting Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),   intent(IN), optional:: tp
  integer(I_P),   intent(IN), optional:: inf
  type(Type_Adj), intent(IN), optional:: adj
  type(Type_BC),  intent(INOUT)::        bc
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

  elemental subroutine get(inf,adj,bc)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for extraction Type_BC attributes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),   intent(OUT), optional:: inf
  type(Type_Adj), intent(OUT), optional:: adj
  type(Type_BC),  intent(IN)::            bc
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(inf)) inf = bc%inf
  if (present(adj)) adj = bc%adj
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get

  ! free dynamic memory
  function Free_Scalar(bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_BC (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(INOUT):: bc
  integer(I4P)::                 err  ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  if (allocated(bc%inf)) deallocate(bc%inf,stat=err)
  if (allocated(bc%adj)) deallocate(bc%adj,stat=err)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Scalar

  function Free_Array1D(bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_BC (array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(INOUT):: bc(:)
  integer(I4P)::                 err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                 i
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

  function Free_Array2D(bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_BC (array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(INOUT):: bc(:,:)
  integer(I4P)::                 err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                 i,j
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

  function Free_Array3D(bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_BC (array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC), intent(INOUT):: bc(:,:,:)
  integer(I4P)::                 err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                 i,j,k
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
  function Write_Bin_Scalar(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, scalar) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit  ! logic unit
  type(Type_bc), intent(IN):: bc
  integer(I_P)::              err   ! Error traping flag: 0 no errors, >0 error occours.
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

  function Write_Ascii_Scalar(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, scalar) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit   ! logic unit
  character(*),  intent(IN):: format ! format specifier
  type(Type_bc), intent(IN):: bc
  integer(I_P)::              err    ! Error traping flag: 0 no errors, >0 error occours.
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

  function Write_Bin_Vectorial1D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, vectorial 1D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit  ! logic unit
  type(Type_bc), intent(IN):: bc(:)
  integer(I_P)::              err   ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::              n
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
  endfunction Write_Bin_Vectorial1D

  function Write_Ascii_Vectorial1D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, vectorial 1D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit   ! logic unit
  character(*),  intent(IN):: format ! format specifier
  type(Type_bc), intent(IN):: bc(:)
  integer(I_P)::              err    ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::              n
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
  endfunction Write_Ascii_Vectorial1D

  function Write_Bin_Vectorial2D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, vectorial 2D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit    ! logic unit
  type(Type_bc), intent(IN):: bc(:,:)
  integer(I_P)::              err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::              n1,n2
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
  endfunction Write_Bin_Vectorial2D

  function Write_Ascii_Vectorial2D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, vectorial 2D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit    ! logic unit
  character(*),  intent(IN):: format  ! format specifier
  type(Type_bc), intent(IN):: bc(:,:)
  integer(I_P)::              err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::              n1,n2
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
  endfunction Write_Ascii_Vectorial2D

  function Write_Bin_Vectorial3D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, vectorial 3D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit      ! logic unit
  type(Type_bc), intent(IN):: bc(:,:,:)
  integer(I_P)::              err       ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::              n1,n2,n3
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
  endfunction Write_Bin_Vectorial3D

  function Write_Ascii_Vectorial3D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, vectorial 3D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN):: unit      ! logic unit
  character(*),  intent(IN):: format    ! format specifier
  type(Type_bc), intent(IN):: bc(:,:,:)
  integer(I_P)::              err       ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::              n1,n2,n3
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
  endfunction Write_Ascii_Vectorial3D

  ! read
  function Read_Bin_Scalar(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, scalar) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit  ! logic unit
  type(Type_bc), intent(INOUT):: bc
  integer(I_P)::                 err   ! Error traping flag: 0 no errors, >0 error occours.
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

  function Read_Ascii_Scalar(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, scalar) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit   ! logic unit
  character(*),  intent(IN)::    format ! format specifier
  type(Type_bc), intent(INOUT):: bc
  integer(I_P)::                 err    ! Error traping flag: 0 no errors, >0 error occours.
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

  function Read_Bin_Vectorial1D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, vectorial 1D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit  ! logic unit
  type(Type_bc), intent(INOUT):: bc(:)
  integer(I_P)::                 err   ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::                 n
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
  endfunction Read_Bin_Vectorial1D

  function Read_Ascii_Vectorial1D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, vectorial 1D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit   ! logic unit
  character(*),  intent(IN)::    format ! format specifier
  type(Type_bc), intent(INOUT):: bc(:)
  integer(I_P)::                 err    ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::                 n
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
  endfunction Read_Ascii_Vectorial1D

  function Read_Bin_Vectorial2D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, vectorial 2D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit    ! logic unit
  type(Type_bc), intent(INOUT):: bc(:,:)
  integer(I_P)::                 err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::                 n1,n2
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
  endfunction Read_Bin_Vectorial2D

  function Read_Ascii_Vectorial2D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, vectorial 2D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit    ! logic unit
  character(*),  intent(IN)::    format  ! format specifier
  type(Type_bc), intent(INOUT):: bc(:,:)
  integer(I_P)::                 err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::                 n1,n2
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
  endfunction Read_Ascii_Vectorial2D

  function Read_Bin_Vectorial3D(unit,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, vectorial 3D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit      ! logic unit
  type(Type_bc), intent(INOUT):: bc(:,:,:)
  integer(I_P)::                 err       ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::                 n1,n2,n3
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
  endfunction Read_Bin_Vectorial3D

  function Read_Ascii_Vectorial3D(unit,format,bc) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, vectorial 3D) Type_BC.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN)::    unit      ! logic unit
  character(*),  intent(IN)::    format    ! format specifier
  type(Type_bc), intent(INOUT):: bc(:,:,:)
  integer(I_P)::                 err       ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::                 n1,n2,n3
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
  endfunction Read_Ascii_Vectorial3D

  ! get_bc
  function get_bc_id(myrank,bc_str) result(bc_id)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for getting integer id of a bc type from the corresponding bc string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: myrank ! Actual rank process.
  character(*), intent(IN):: bc_str ! String of bc.
  integer(I_P)::             bc_id  ! Id of bc.
  integer(I_P)::             b      ! Boundary conditions counter.
#ifdef MPI2
  integer(I_P)::             err    ! Error for MPI comunications.
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

  function get_bc_str(myrank,bc_id) result(bc_str)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for getting string id of a bc type from the corresponding bc integer id.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: myrank ! Actual rank process.
  integer(I_P), intent(IN):: bc_id  ! Id of bc.
  character(3)::             bc_str ! String of bc.
  integer(I_P)::             b      ! Boundary conditions counter.
#ifdef MPI2
  integer(I_P)::             err    ! Error for MPI comunications.
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
endmodule Data_Type_BC
