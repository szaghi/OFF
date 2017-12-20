!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_BC_in1DerivedType Data_Type_BC_in1
!> Module definition of Type_BC_in1
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_BC_in1PrivateProcedure Data_Type_BC_in1
!> Module definition of Type_BC_in1
!> @}

!> @brief This module contains the definition of Type_BC_in1 and its procedures.
!> Type_BC_in1 is a derived type containing inflow 1 boundary conditions.
module Data_Type_BC_in1
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                              ! Integers and reals precision definition.
USE Data_Type_Primitive, only: Type_Primitive ! Definition of Type_Primitive.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing inflow 1 boundary conditions.
!> @ingroup Data_Type_BC_in1DerivedType
type, public:: Type_BC_in1
  integer(I4P)::                      Nin1 = 0_I4P !< Number of inflow 1 boundary conditions.
  type(Type_Primitive), allocatable:: P(:)         !< Inflow 1 boundary conditions primitive variables [1:Nin1].
  contains
    procedure:: free  => free_bc  ! Procedure for freeing dynamic memory.
    procedure:: alloc => alloc_bc ! Procedure for allocating dynamic memory.
    final::     finalize          ! Procedure for freeing dynamic memory when finalizing.
    ! operators overloading
    generic:: assignment(=) => assign_bc
    ! private procedures
    procedure, pass(bc1), private:: assign_bc
endtype Type_BC_in1
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_BC_in1PrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_bc(bc_in1)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC_in1), intent(INOUT):: bc_in1 !< Inflow 1 boundary conditions.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(bc_in1%P)) then
    call bc_in1%P%free ; deallocate(bc_in1%P)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_bc

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc_bc(bc_in1,Ns)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC_in1), intent(INOUT):: bc_in1 !< Inflow 1 boundary conditions.
  integer(I4P),       intent(IN)::    Ns     !< Number of species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call bc_in1%free
  if (bc_in1%Nin1>0) then
    allocate(bc_in1%P(1:bc_in1%Nin1))
    call bc_in1%P%alloc(Ns=Ns)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_bc

  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine finalize(bc_in1)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_BC_in1), intent(INOUT):: bc_in1 !< Inflow 1 boundary conditions.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call bc_in1%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  ! Assignment (=)
  !> @brief Procedure for assignment between two boundary conditions variables.
  elemental subroutine assign_bc(bc1,bc2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_BC_in1), intent(INOUT):: bc1
  type(Type_BC_in1),  intent(IN)::    bc2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bc1%Nin1 = bc2%Nin1
  bc1%P    = bc2%P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_bc
  !> @}
endmodule Data_Type_BC_in1
