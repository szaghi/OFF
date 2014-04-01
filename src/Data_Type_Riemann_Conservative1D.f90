!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Riemann_Conservative1DDerivedType Data_Type_Riemann_Conservative1D
!> Module definition of Type_Riemann_Conservative1D
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Riemann_Conservative1DPrivateProcedure Data_Type_Riemann_Conservative1D
!> Module definition of Type_Riemann_Conservative1D
!> @}

!> @brief Module Data_Type_Riemann_Conservative1D contains the definition of Type_Riemann_Conservative1D, that defines the Riemann
!> conservative fluxes for 1D Riemann problem.
module Data_Type_Riemann_Conservative1D
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                  ! Integers and reals precision definition.
USE Data_Type_Riemann_Primitive1D, only: Type_Riemann_Primitive1D ! Definition of Type_Riemann_Primitive1D.
USE Lib_Thermodynamic_Laws_Ideal,  only: H                        ! Library of thermodynamic laws for ideal gas.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_Riemann_Conservative1D.
!> @ingroup Data_Type_Riemann_Conservative1DDerivedType
type, public:: Type_Riemann_Conservative1D
  real(R8P):: r = 0._R8P !< Flux of mass conservation.
  real(R8P):: u = 0._R8P !< Flux of momentum conservation.
  real(R8P):: E = 0._R8P !< Flux of energy conservation.
  contains
    procedure:: compute => compute_self ! Procedure for computing Riemann flux from Riemann state.
    procedure:: print   => print_self   ! Procedure for printing self with a pretty format.
    ! operators overloading
#include 'Data_Type_Bounds_Proc_OpOverloading.inc'
endtype Type_Riemann_Conservative1D
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Riemann_Conservative1DPrivateProcedure
  !> @{
  !> @brief Procedure for computing Riemann flux from Riemann state.
  elemental subroutine compute_self(self,state)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Conservative1D), intent(INOUT):: self    !< Riemann flux.
  type(Type_Riemann_Primitive1D),     intent(IN)::    state   !< Riemann state.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%r = state%r*state%u
  self%u = self%r*state%u + state%p
  self%E = self%r*H(p=state%p,r=state%r,u=state%u,g=state%g)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_self

  !> @brief Procedure for printing self with a pretty format.
  subroutine print_self(self,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Conservative1D), intent(IN)::  self    !< Riemann flux.
  character(*), optional,             intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional,             intent(OUT):: iostat  !< IO error.
  character(*), optional,             intent(OUT):: iomsg   !< IO error message.
  integer(I4P),                       intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::                   prefd   !< Prefixing string.
  integer(I4P)::                                    iostatd !< IO error.
  character(500)::                                  iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' rs='//str(n=self%r)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' ru='//str(n=self%u)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' rE='//str(n=self%E)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_self

  ! Operators overloading.
  ! Operator (=)
#define self_type_ class(Type_Riemann_Conservative1D), intent(INOUT):: self
#define ass_scal_ self%r=real(scal,R8P);self%u=real(scal,R8P);self%E=real(scal,R8P)
  !> @brief Procedure for assignment between two selfs.
  elemental subroutine assign_self(self1,self2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Conservative1D), intent(INOUT):: self1
  type(Type_Riemann_Conservative1D),  intent(IN)::    self2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self1%r = self2%r
  self1%u = self2%u
  self1%E = self2%E
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_self
#include 'Data_Type_Bounds_Proc_AssDefinitions.inc'
#undef self_type_
#undef ass_scal_

  ! Operator (*)
#define self_type_ class(Type_Riemann_Conservative1D), intent(IN):: self
#define mul_type_ type(Type_Riemann_Conservative1D):: mul
#define mul_scal_ mul%r=real(scal,R8P)*self%r;mul%u=real(scal,R8P)*self%u;mul%E=real(scal,R8P)*self%E
  !> @brief Procedure for multiply (by components) two selfs.
  elemental function self_mul_self(self1,self2) result(mul)

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Conservative1D), intent(IN):: self1
  type(Type_Riemann_Conservative1D),  intent(IN):: self2
  type(Type_Riemann_Conservative1D)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%r = self1%r * self2%r
  mul%u = self1%u * self2%u
  mul%E = self1%E * self2%E
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_mul_self
#include 'Data_Type_Bounds_Proc_MulDefinitions.inc'
#undef self_type_
#undef mul_type_
#undef mul_scal_

  ! Operator (/)
#define self_type_ class(Type_Riemann_Conservative1D), intent(IN):: self
#define div_type_ type(Type_Riemann_Conservative1D):: div
#define div_scal_ div%r=self%r/real(scal,R8P);div%u=self%u/real(scal,R8P);div%E=self%E/real(scal,R8P)
  !> @brief Procedure for divide self for self.
  elemental function self_div_self(self1,self2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Conservative1D), intent(IN):: self1
  type(Type_Riemann_Conservative1D),  intent(IN):: self2
  type(Type_Riemann_Conservative1D)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%r = self1%r / self2%r
  div%u = self1%u / self2%u
  div%E = self1%E / self2%E
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_div_self
#include 'Data_Type_Bounds_Proc_DivDefinitions.inc'
#undef self_type_
#undef div_type_
#undef div_scal_

  ! Operator (+)
#define self_type_ class(Type_Riemann_Conservative1D), intent(IN):: self
#define summ_type_ type(Type_Riemann_Conservative1D):: summ
#define sum_scal_ summ%r=self%r+real(scal,R8P);summ%u=self%u+real(scal,R8P);summ%E=self%E+real(scal,R8P)
  !> @brief Procedure for applay unary + to a self.
  elemental function positive_self(self) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Conservative1D), intent(IN):: self
  type(Type_Riemann_Conservative1D)::              pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pos%r = + self%r
  pos%u = + self%u
  pos%E = + self%E
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_self

  !> @brief Procedure for sum self and self.
  elemental function self_sum_self(self1,self2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Conservative1D), intent(IN):: self1
  type(Type_Riemann_Conservative1D),  intent(IN):: self2
  type(Type_Riemann_Conservative1D)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%r = self1%r + self2%r
  summ%u = self1%u + self2%u
  summ%E = self1%E + self2%E
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sum_self
#include 'Data_Type_Bounds_Proc_SumDefinitions.inc'
#undef self_type_
#undef summ_type_
#undef sum_scal_

  ! Operator (-)
#define self_type_ class(Type_Riemann_Conservative1D), intent(IN):: self
#define sub_type_ type(Type_Riemann_Conservative1D):: sub
#define self_sub_scal_ sub%r=self%r-real(scal,R8P);sub%u=self%u-real(scal,R8P);sub%E=self%E-real(scal,R8P)
#define scal_sub_self_ sub%r=real(scal,R8P)-self%r;sub%u=real(scal,R8P)-self%u;sub%E=real(scal,R8P)-self%E
  !> @brief Procedure for applay unary - to a self.
  elemental function negative_self(self) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Conservative1D), intent(IN):: self
  type(Type_Riemann_Conservative1D)::              neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neg%r = - self%r
  neg%u = - self%u
  neg%E = - self%E
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_self

  !> @brief Procedure for subtract self and self.
  elemental function self_sub_self(self1,self2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Conservative1D), intent(IN):: self1
  type(Type_Riemann_Conservative1D),  intent(IN):: self2
  type(Type_Riemann_Conservative1D)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%r = self1%r - self2%r
  sub%u = self1%u - self2%u
  sub%E = self1%E - self2%E
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sub_self
#include 'Data_Type_Bounds_Proc_SubDefinitions.inc'
#undef self_type_
#undef sub_type_
#undef self_sub_scal_
#undef scal_sub_self_
  !> @}
endmodule Data_Type_Riemann_Conservative1D
