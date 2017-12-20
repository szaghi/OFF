!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_CellDerivedType Data_Type_Cell
!> Module definition of Type_Cell
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_CellInterface Data_Type_Cell
!> Module definition of Type_Cell
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_CellPrivateProcedure Data_Type_Cell
!> Module definition of Type_Cell
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_CellPublicProcedure Data_Type_Cell
!> Module definition of Type_Cell
!> @}

!> This module contains the definition of Type_Cell and its procedures.
!> @todo \b DocWriteRead: Complete the documentation of write and read functions
module Data_Type_Cell
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                         ! Integers and reals precision definition.
USE Data_Type_Conservative,    only: Type_Conservative   ! Definition of Type_Conservative.
USE Data_Type_Primitive,       only: Type_Primitive      ! Definition of Type_Primitive.
USE Data_Type_Species,         only: Type_Species        ! Definition of Type_Species.
USE Data_Type_Vector,          only: Type_Vector         ! Definition of Type_Vector.
USE Lib_Variables_Conversions, only: prim2cons,cons2prim ! Pocedures for varibles set conversions.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: cells_bc_set_ext,cells_bc_set_ref,cells_bc_set_per
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing cell-level data.
!> Cell-level type contains data (mesh, boundary conditions and fluid dynamic data) of each cell of the numerical grid.
!> @ingroup Data_Type_CellDerivedType
type, public:: Type_Cell
  real(R8P)::                            V = 0._R8P  !< Cell volume.
  type(Type_Vector)::                    cent        !< Cell center coordinates.
  real(R8P)::                            Dt = 0._R8P !< Local time step.
  type(Type_Primitive)::                 P           !< Primitive variables.
  type(Type_Conservative)::              U           !< Conservative variables.
  type(Type_Conservative), allocatable:: KS(:)       !< Runge-Kutta stages of conservative variables [1:rk_ord].
  contains
    procedure:: free      => free_cell       ! Procedure for freeing dynamic memory.
    procedure:: alloc     => alloc_cell      ! Procedure for allocating dynamic memory.
    procedure:: prim2cons => prim2cons_cell  ! Procedure for converting primitive to conservative variables.
    procedure:: cons2prim => cons2prim_cell  ! Procedure for converting conservative to primitive variables.
    procedure:: load      => load_cell_self  ! Procedure for loading cell data.
    procedure:: save      => save_cell_self  ! Procedure for saving cell data.
    procedure:: print     => print_cell_self ! Procedure for printing cell data with a pretty format.
    final::     finalize                     ! Procedure for freeing dynamic memory when finalizing.
    ! operators overloading
    generic:: assignment(=) => assign_cell
    ! private procedures
    procedure, pass(cell1), private:: assign_cell
endtype Type_Cell
!> @brief Pointer of Type_SCell for creating array of pointers of Type_SCell.
!> @ingroup Data_Type_CellDerivedType
type, public:: Type_Cell_Ptr
  type(Type_Cell),  pointer:: p => null()
endtype Type_Cell_Ptr
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_CellPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_cell(cell)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell), intent(INOUT):: cell !< Cell data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call cell%P%free
  call cell%U%free
  if (allocated(cell%KS)) then
    call cell%KS%free ; deallocate(cell%KS)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_cell

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(cell)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Cell), intent(INOUT):: cell !< Cell data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call cell%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc_cell(cell,Ns,Nrk)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell), intent(INOUT):: cell !< Cell data.
  integer(I4P),     intent(IN)::    Ns   !< Number of species.
  integer(I1P),     intent(IN)::    Nrk  !< Number of Runge-Kutta stages.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call cell%free ; allocate(cell%KS(1:Nrk))
  call cell%P%alloc(Ns=Ns)
  call cell%U%alloc(Ns=Ns)
  call cell%KS%alloc(Ns=Ns)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_cell

  !> @brief Procedure for converting primitive variables to conservative variables.
  elemental subroutine prim2cons_cell(cell)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell), intent(INOUT):: cell !< Cell data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call prim2cons(prim=cell%P,cons=cell%U)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prim2cons_cell

  !> @brief Procedure for converting conservative variables to primitive variables.
  elemental subroutine cons2prim_cell(cell,species0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell),   intent(INOUT):: cell     !< Cell data.
  type(Type_Species), intent(IN)::    species0 !< Initial species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call cons2prim(cons=cell%U,prim=cell%P,species0=species0)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine cons2prim_cell

  !> @brief Procedure for loading of Type_Cell data.
  subroutine load_cell_self(cell,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell),       intent(INOUT):: cell    !< Cell data.
  integer(I8P), optional, intent(IN)::    pos     !< Position specifier.
  integer(I4P), optional, intent(OUT)::   iostat  !< IO error.
  character(*), optional, intent(OUT)::   iomsg   !< IO error message.
  integer(I4P),           intent(IN)::    unit    !< Logic unit.
  integer(I4P)::                          iostatd !< IO error.
  character(500)::                        iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    read(            unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)cell%Dt
  else
    read(            unit=unit,        iostat=iostatd,iomsg=iomsgd)cell%Dt
  endif
    call cell%P%load(unit=unit,        iostat=iostatd,iomsg=iomsgd)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = adjustl(trim(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_cell_self

  !> @brief Procedure for saving Type_Cell data.
  subroutine save_cell_self(cell,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell),       intent(IN)::  cell    !< Cell data.
  integer(I8P), optional, intent(IN)::  pos     !< Position specifier.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    write(           unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)cell%Dt
  else
    write(           unit=unit,        iostat=iostatd,iomsg=iomsgd)cell%Dt
  endif
    call cell%P%save(unit=unit,        iostat=iostatd,iomsg=iomsgd)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = adjustl(trim(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_cell_self

  !> @brief Procedure for printing primitives with a pretty format.
  subroutine print_cell_self(cell,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell),       intent(IN)::  cell    !< Cell data.
  character(*), optional, intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::       prefd   !< Prefixing string.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  integer(I4P)::                        s       !< RK stages counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(               unit=unit,fmt='(A)',       iostat=iostatd,iomsg=iomsgd)prefd//' Volume='//str(n=cell%V)
  write(               unit=unit,fmt='(A)',       iostat=iostatd,iomsg=iomsgd)prefd//' Center Coordinates:'
  call cell%cent%print(unit=unit,pref=prefd//'  ',iostat=iostatd,iomsg=iomsgd)
  write(               unit=unit,fmt='(A)',       iostat=iostatd,iomsg=iomsgd)prefd//' Dt='//str(n=cell%Dt)
  write(               unit=unit,fmt='(A)',       iostat=iostatd,iomsg=iomsgd)prefd//' Primitives:'
  call cell%P%print(   unit=unit,pref=prefd//'  ',iostat=iostatd,iomsg=iomsgd)
  write(               unit=unit,fmt='(A)',       iostat=iostatd,iomsg=iomsgd)prefd//' Conservatives:'
  call cell%U%print(   unit=unit,pref=prefd//'  ',iostat=iostatd,iomsg=iomsgd)
  if (allocated(cell%KS)) then
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' RK Conservatives Stages:'
    do s=lbound(cell%KS,dim=1),ubound(cell%KS,dim=1)
      write(                unit=unit,fmt='(A)',         iostat=iostatd,iomsg=iomsgd)prefd//'   stage s('//trim(str(.true.,s))//')'
      call cell%KS(s)%print(unit=unit,pref=prefd//'    ',iostat=iostatd,iomsg=iomsgd)
    enddo
  endif
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_cell_self

  ! Assignment (=)
  !> @brief Procedure for assignment between two cells variables.
  elemental subroutine assign_cell(cell1,cell2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell), intent(INOUT):: cell1
  type(Type_Cell),  intent(IN)::    cell2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
                           cell1%V    = cell2%V
                           cell1%cent = cell2%cent
                           cell1%Dt   = cell2%Dt
                           cell1%P    = cell2%P
                           cell1%U    = cell2%U
  if (allocated(cell2%KS)) cell1%KS   = cell2%KS
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_cell
  !> @}

  !> @ingroup Data_Type_CellPublicProcedure
  !> @{
  !> @brief Procedure for imposing extrapolation of ghost cells from internal ones along a generic direction.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) the array \b C is declared as assumed-shape
  !> with only the lower bound defined. Its extentions is: C [1-gc(1):N+gc(2)].
  pure subroutine cells_bc_set_ext(gc,ic,N,boundary,C)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),    intent(IN)::    gc(1:2)     !< Number of ghost cells.
  integer(I4P),    intent(IN)::    ic          !< Number of internal cells used for extrapolation (1 or gc).
  integer(I4P),    intent(IN)::    N           !< Number of internal cells.
  character(1),    intent(IN)::    boundary    !< Boundary left ('l') or right ('r').
  type(Type_Cell), intent(INOUT):: C(1-gc(1):) !< Cells data [1-gc(1):N+gc(2)].
  integer(I4P)::                   i           !< Cell counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (boundary=='l') then
    if (ic==1.OR.N<gc(1)) then
      ! extrapolation using only the cell 1
      do i=1-gc(1),0
        C(i)%P = C(1)%P
      enddo
    else
      ! extrapolation using the cells 1,2,...,gc
      do i=1-gc(1),0
        C(i)%P = C(-i+1)%P
      enddo
    endif
  endif
  if (boundary=='r') then
    if (ic==1.OR.N<gc(2)) then
      ! extrapolation using only the cell N
      do i=N+1,N+gc(2)
        C(i)%P = C(N)%P
      enddo
    else
      ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
      do i=N+1,N+gc(2)
        C(i)%P = C(N+1-(i-N))%P
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine cells_bc_set_ext

  !> @brief Procedure for imposing reflective boundary conditions along a generic direction.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) the array \b C is declared as assumed-shape
  !> with only the lower bound defined. Its extentions is: C [1-gc(1):N+gc(2)].
  pure subroutine cells_bc_set_ref(gc,ic,N,NF,boundary,C)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN)::    gc(1:2)     !< Number of ghost cells.
  integer(I4P),      intent(IN)::    ic          !< Number of internal cells used for extrapolation (1 or gc).
  integer(I4P),      intent(IN)::    N           !< Number of internal cells.
  type(Type_Vector), intent(IN)::    NF          !< Face normal.
  character(1),      intent(IN)::    boundary    !< Boundary left ('l') or right ('r').
  type(Type_Cell),   intent(INOUT):: C(1-gc(1):) !< Cells data [1-gc(1):N+gc(2)].
  integer(I4P)::                     i           !< Cell counter.
  type(Type_Vector)::                vr          !< Reflected velocity vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (boundary=='l') then
    if (ic==1.OR.N<gc(1)) then
      ! reflection using only the cell 1
      vr = C(1)%P%v - (2._R_P*(C(1)%P%v.paral.NF)) ! reflected velocity
      do i=1-gc(1),0
        C(i)%P%r = C(1)%P%r
        C(i)%P%v = vr
        C(i)%P%p = C(1)%P%p
        C(i)%P%d = C(1)%P%d
        C(i)%P%g = C(1)%P%g
      enddo
    else
      ! reflection using the cells 1,2,...,gc
      do i=1-gc(1),0
        vr = C(-i+1)%P%v - (2._R_P*(C(-i+1)%P%v.paral.NF)) ! reflected velocity
        C(i)%P%r = C(-i+1)%P%r
        C(i)%P%v = vr
        C(i)%P%p = C(-i+1)%P%p
        C(i)%P%d = C(-i+1)%P%d
        C(i)%P%g = C(-i+1)%P%g
      enddo
    endif
  endif
  if (boundary=='r') then
    if (ic==1.OR.N<gc(2)) then
      ! reflection using only the cell N
      vr = C(N)%P%v - (2._R_P*(C(N)%P%v.paral.NF)) ! reflected velocity
      do i=N+1,N+gc(2)
        C(i)%P%r = C(N)%P%r
        C(i)%P%v = vr
        C(i)%P%p = C(N)%P%p
        C(i)%P%d = C(N)%P%d
        C(i)%P%g = C(N)%P%g
      enddo
    else
      ! reflection using the cells N-gc,N-gc+1,N-gc+2,...,N
      do i=N+1,N+gc(2)
        vr = C(N+1-(i-N))%P%v - (2._R_P*(C(N+1-(i-N))%P%v.paral.NF)) ! reflected velocity
        C(i)%P%r = C(N+1-(i-N))%P%r
        C(i)%P%v = vr
        C(i)%P%p = C(N+1-(i-N))%P%p
        C(i)%P%d = C(N+1-(i-N))%P%d
        C(i)%P%g = C(N+1-(i-N))%P%g
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine cells_bc_set_ref

  !> @brief Procedure for imposing periodic boundary conditions along a generic direction.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) the array \b C is declared as assumed-shape
  !> with only the lower bound defined. Its extentions is: C [1-gc(1):N+gc(2)].
  pure subroutine cells_bc_set_per(gc,ic,N,boundary,C)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),    intent(IN)::    gc(1:2)     !< Number of ghost cells.
  integer(I_P),    intent(IN)::    ic          !< Number of internal cells used for extrapolation (1 or gc).
  integer(I_P),    intent(IN)::    N           !< Number of internal cells.
  character(1),    intent(IN)::    boundary    !< Boundary left ('l') or right ('r').
  type(Type_Cell), intent(INOUT):: C(1-gc(1):) !< Cells data [1-gc(1):N+gc(2)].
  integer(I_P)::                   i           !< Cell counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (boundary=='l') then
    if (ic==1.OR.N<gc(1)) then
      ! extrapolation using only the cell N
      do i=1-gc(1),0
        C(i)%P = C(N)%P
      enddo
    else
      ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
      do i=1-gc(1),0
        C(i)%P = C(i+N)%P
      enddo
    endif
  endif
  if (boundary=='r') then
    if (ic==1.OR.N<gc(2)) then
      ! extrapolation using only the cell 1
      do i=N+1,N+gc(2)
        C(i)%P = C(1)%P
      enddo
    else
      ! extrapolation using the cells 1,2,...,gc
      do i=N+1,N+gc(2)
        C(i)%P = C(i-N)%P
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine cells_bc_set_per
  !> @}
endmodule Data_Type_Cell
