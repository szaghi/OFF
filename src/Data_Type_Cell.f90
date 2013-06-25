!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_CellDerivedType Data_Type_Cell
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_CellPublicProcedure Data_Type_Cell
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_CellPrivateProcedure Data_Type_Cell
!> @}

!> This module contains the definition of Type_Cell and its procedures.
!> @todo \b DocWriteRead: Complete the documentation of write and read functions
module Data_Type_Cell
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision           ! Integers and reals precision definition.
USE Data_Type_Conservative !< Definition of Type_Conservative.
USE Data_Type_Primitive    !< Definition of Type_Primitive.
USE Data_Type_Vector       !< Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
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
    procedure, non_overridable:: init => init_cell ! Procedure for initilizing allocatable variables.
    procedure, non_overridable:: free => free_cell ! Procedure for freeing the memory of allocatable variables.
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
  !> Subroutine for freeing dynamic data of Type_Cell variables.
  elemental subroutine free_cell(cell)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell), intent(INOUT):: cell !< Cell data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call cell%P%free
  call cell%U%free
  if (allocated(cell%KS)) deallocate(cell%KS)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_cell

  !> Subroutine for initializing dynamic data of Type_Cell variables.
  elemental subroutine init_cell(cell,Ns,prim0,cons0,rk_ord)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell),        intent(INOUT)::        cell   !< Cell data.
  integer(I_P),            intent(IN), optional:: Ns     !< Number of species.
  type(Type_Primitive),    intent(IN), optional:: prim0  !< Primitive inizialization data.
  type(Type_Conservative), intent(IN), optional:: cons0  !< Conservative inizialization data.
  integer(I1P),            intent(IN)::           rk_ord !< Number of Runge-Kutta stages.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cell%KS)) deallocate(cell%KS) ; allocate(cell%KS(1:rk_ord))
  if (present(cons0)) then
    call cell%U%init(cons0=cons0)
    call cell%KS%init(cons0=cons0)
  elseif (present(Ns)) then
    call cell%U%init(Ns=Ns)
    call cell%KS%init(Ns=Ns)
  endif
  if (present(prim0)) then
    call cell%P%init(prim0=prim0)
  elseif (present(Ns)) then
    call cell%P%init(Ns=Ns)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_cell
  !> @}
endmodule Data_Type_Cell
