!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Riemann_Primitive1DDerivedType Data_Type_Riemann_Primitive1D
!> Module definition of Riemann state in 1D, Type_Riemann_Primitive1D
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Riemann_Primitive1DPrivateProcedure Data_Type_Riemann_Primitive1D
!> Module definition of Riemann state in 1D, Type_Riemann_Primitive1D
!> @}

!> @brief Module Data_Type_Riemann_Primitive1D contains the definition of Type_Riemann_Primitive1D, that defines the Riemann
!> primitive variables of state (1 and 4, left and right) for 1D Riemann problem.
module Data_Type_Riemann_Primitive1D
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                          ! Integers and reals precision definition.
USE Lib_Thermodynamic_Laws_Ideal, only: a ! Library of thermodynamic laws for ideal calorically perfect gas.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_Riemann_Primitive1D.
!> @ingroup Data_Type_Riemann_Primitive1DDerivedType
type, public:: Type_Riemann_Primitive1D
  real(R8P):: r = 0._R8P !< Density.
  real(R8P):: u = 0._R8P !< Velocity.
  real(R8P):: p = 0._R8P !< Pressure.
  real(R8P):: g = 0._R8P !< Specific heats ratio.
  ! dependent members
  real(R8P):: a     = 0._R8P !< Speed of sound.
  real(R8P):: delta = 0._R8P !< (g-1)/2.
  real(R8P):: eta   = 0._R8P !< 2g/(g-1).
  contains
    procedure:: set               => set_self               ! Procedure for setting self data.
    procedure:: compute_dependent => compute_dependent_self ! Procedure for computing dependent members.
    procedure:: print             => print_self             ! Procedure for printing self with a pretty format.
    ! operators overloading
#include "Data_Type_Bounds_Proc_OpOverloading.inc"
endtype Type_Riemann_Primitive1D
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Riemann_Primitive1DPrivateProcedure
  !> @{
  !> @brief Procedure for setting self data.
  elemental subroutine set_self(self,r,u,p,g,a,delta,eta)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Primitive1D), intent(INOUT):: self  !< Riemann state.
  real(R8P), optional,             intent(IN)::    r     !< Density.
  real(R8P), optional,             intent(IN)::    u     !< Velocity.
  real(R8P), optional,             intent(IN)::    p     !< Pressure.
  real(R8P), optional,             intent(IN)::    g     !< Specific heats ratio.
  real(R8P), optional,             intent(IN)::    a     !< Speed of sound.
  real(R8P), optional,             intent(IN)::    delta !< (g-1)/2.
  real(R8P), optional,             intent(IN)::    eta   !< 2g/(g-1).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(r    )) self%r     = r
  if (present(u    )) self%u     = u
  if (present(p    )) self%p     = p
  if (present(g    )) self%g     = g
  if (present(a    )) self%a     = a
  if (present(delta)) self%delta = delta
  if (present(eta  )) self%eta   = eta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_self

  !> @brief Procedure for computing dependent members.
  elemental subroutine compute_dependent_self(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Primitive1D), intent(INOUT):: self !< Riemann state.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%a     = a(p=self%p,r=self%r,g=self%g)
  self%delta = 0.5_R8P*(self%g-1.0_R8P)
  self%eta   = 2._R8P*self%g/(self%g-1.0_R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_dependent_self

  !> @brief Procedure for printing self with a pretty format.
  subroutine print_self(self,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Primitive1D), intent(IN)::  self    !< Riemann state.
  character(*), optional,          intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional,          intent(OUT):: iostat  !< IO error.
  character(*), optional,          intent(OUT):: iomsg   !< IO error message.
  integer(I4P),                    intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::                prefd   !< Prefixing string.
  integer(I4P)::                                 iostatd !< IO error.
  character(500)::                               iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' r='//str(n=self%r)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' u='//str(n=self%u)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' p='//str(n=self%p)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' g='//str(n=self%g)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_self

  ! Operators overloading.
  ! Operator (=)
#define self_type_ class(Type_Riemann_Primitive1D), intent(INOUT):: self
#define ass_scal_1_ self%r=real(scal,R8P);self%u=real(scal,R8P);self%p=real(scal,R8P);self%g=real(scal,R8P);
#define ass_scal_2_ self%a=real(scal,R8P);self%delta=real(scal,R8P);self%eta=real(scal,R8P)
#define ass_scal_ ass_scal_1_ ass_scal_2_
  !> @brief Procedure for assignment between two selfs.
  elemental subroutine assign_self(self1,self2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Primitive1D), intent(INOUT):: self1
  type(Type_Riemann_Primitive1D),  intent(IN)::    self2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self1%r     = self2%r
  self1%u     = self2%u
  self1%p     = self2%p
  self1%g     = self2%g
  self1%a     = self2%a
  self1%delta = self2%delta
  self1%eta   = self2%eta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_self
#include "Data_Type_Bounds_Proc_AssDefinitions.inc"
#undef self_type_
#undef ass_scal_1_
#undef ass_scal_2_
#undef ass_scal_

  ! Operator (*)
#define self_type_ class(Type_Riemann_Primitive1D), intent(IN):: self
#define mul_type_ type(Type_Riemann_Primitive1D):: mul
#define mul_scal_1_ mul%r=real(scal,R8P)*self%r;mul%u=real(scal,R8P)*self%u;mul%p=real(scal,R8P)*self%p;mul%g=real(scal,R8P)*self%g;
#define mul_scal_2_ mul%a=real(scal,R8P)*self%a;mul%delta=real(scal,R8P)*self%delta;mul%eta=real(scal,R8P)*self%eta
#define mul_scal_ mul_scal_1_ mul_scal_2_
  !> @brief Procedure for multiply (by components) two selfs.
  elemental function self_mul_self(self1,self2) result(mul)

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Primitive1D), intent(IN):: self1
  type(Type_Riemann_Primitive1D),  intent(IN):: self2
  type(Type_Riemann_Primitive1D)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%r     = self1%r     * self2%r
  mul%u     = self1%u     * self2%u
  mul%p     = self1%p     * self2%p
  mul%g     = self1%g     * self2%g
  mul%a     = self1%a     * self2%a
  mul%delta = self1%delta * self2%delta
  mul%eta   = self1%eta   * self2%eta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_mul_self
#include "Data_Type_Bounds_Proc_MulDefinitions.inc"
#undef self_type_
#undef mul_type_
#undef mul_scal_1_
#undef mul_scal_2_
#undef mul_scal_

  ! Operator (/)
#define self_type_ class(Type_Riemann_Primitive1D), intent(IN):: self
#define div_type_ type(Type_Riemann_Primitive1D):: div
#define div_scal_1_ div%r=self%r/real(scal,R8P);div%u=self%u/real(scal,R8P);div%p=self%p/real(scal,R8P);div%g=self%g/real(scal,R8P);
#define div_scal_2_ div%a=self%a/real(scal,R8P);div%delta=self%delta/real(scal,R8P);div%eta=self%eta/real(scal,R8P)
#define div_scal_ div_scal_1_ div_scal_2_
  !> @brief Procedure for divide self for self.
  elemental function self_div_self(self1,self2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Primitive1D), intent(IN):: self1
  type(Type_Riemann_Primitive1D),  intent(IN):: self2
  type(Type_Riemann_Primitive1D)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%r     = self1%r     / self2%r
  div%u     = self1%u     / self2%u
  div%p     = self1%p     / self2%p
  div%g     = self1%g     / self2%g
  div%a     = self1%a     / self2%a
  div%delta = self1%delta / self2%delta
  div%eta   = self1%eta   / self2%eta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_div_self
#include "Data_Type_Bounds_Proc_DivDefinitions.inc"
#undef self_type_
#undef div_type_
#undef div_scal_1_
#undef div_scal_2_
#undef div_scal_

  ! Operator (+)
#define self_type_ class(Type_Riemann_Primitive1D), intent(IN):: self
#define summ_type_ type(Type_Riemann_Primitive1D):: summ
#define sum_scal_1_ summ%r=self%r+real(scal,R8P);summ%u=self%u+real(scal,R8P);summ%p=self%p+real(scal,R8P);
#define sum_scal_2_ summ%g=self%g+real(scal,R8P);summ%a=self%a+real(scal,R8P);summ%delta=self%delta+real(scal,R8P);
#define sum_scal_3_ summ%eta=self%eta+real(scal,R8P)
#define sum_scal_ sum_scal_1_ sum_scal_2_ sum_scal_3_
  !> @brief Procedure for applay unary + to a self.
  elemental function positive_self(self) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Primitive1D), intent(IN):: self
  type(Type_Riemann_Primitive1D)::              pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pos%r     = + self%r
  pos%u     = + self%u
  pos%p     = + self%p
  pos%g     = + self%g
  pos%a     = + self%a
  pos%delta = + self%delta
  pos%eta   = + self%eta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_self

  !> @brief Procedure for sum self and self.
  elemental function self_sum_self(self1,self2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Primitive1D), intent(IN):: self1
  type(Type_Riemann_Primitive1D),  intent(IN):: self2
  type(Type_Riemann_Primitive1D)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%r     = self1%r     + self2%r
  summ%u     = self1%u     + self2%u
  summ%p     = self1%p     + self2%p
  summ%g     = self1%g     + self2%g
  summ%a     = self1%a     + self2%a
  summ%delta = self1%delta + self2%delta
  summ%eta   = self1%eta   + self2%eta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sum_self
#include "Data_Type_Bounds_Proc_SumDefinitions.inc"
#undef self_type_
#undef summ_type_
#undef sum_scal_1_
#undef sum_scal_2_
#undef sum_scal_3_
#undef sum_scal_

  ! Operator (-)
#define self_type_ class(Type_Riemann_Primitive1D), intent(IN):: self
#define sub_type_ type(Type_Riemann_Primitive1D):: sub
#define self_sub_scal_1_ sub%r=self%r-real(scal,R8P);sub%u=self%u-real(scal,R8P);sub%p=self%p-real(scal,R8P);
#define self_sub_scal_2_ sub%g=self%g-real(scal,R8P);sub%a=self%a-real(scal,R8P);sub%delta=self%delta-real(scal,R8P);
#define self_sub_scal_3_ sub%eta=self%eta-real(scal,R8P)
#define scal_sub_self_1_ sub%r=real(scal,R8P)-self%r;sub%u=real(scal,R8P)-self%u;sub%p=real(scal,R8P)-self%p;
#define scal_sub_self_2_ sub%g=real(scal,R8P)-self%g;sub%a=real(scal,R8P)-self%a;sub%delta=real(scal,R8P)-self%delta;
#define scal_sub_self_3_ sub%eta=real(scal,R8P)-self%eta
#define self_sub_scal_ self_sub_scal_1_ self_sub_scal_2_ self_sub_scal_3_
#define scal_sub_self_ scal_sub_self_1_ scal_sub_self_2_ scal_sub_self_3_
  !> @brief Procedure for applay unary - to a self.
  elemental function negative_self(self) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Primitive1D), intent(IN):: self
  type(Type_Riemann_Primitive1D)::              neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neg%r     = - self%r
  neg%u     = - self%u
  neg%p     = - self%p
  neg%g     = - self%g
  neg%a     = - self%a
  neg%delta = - self%delta
  neg%eta   = - self%eta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_self

  !> @brief Procedure for subtract self and self.
  elemental function self_sub_self(self1,self2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_Primitive1D), intent(IN):: self1
  type(Type_Riemann_Primitive1D),  intent(IN):: self2
  type(Type_Riemann_Primitive1D)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%r     = self1%r     - self2%r
  sub%u     = self1%u     - self2%u
  sub%p     = self1%p     - self2%p
  sub%g     = self1%g     - self2%g
  sub%a     = self1%a     - self2%a
  sub%delta = self1%delta - self2%delta
  sub%eta   = self1%eta   - self2%eta
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sub_self
#include "Data_Type_Bounds_Proc_SubDefinitions.inc"
#undef self_type_
#undef sub_type_
#undef self_sub_scal_1_
#undef self_sub_scal_2_
#undef self_sub_scal_3_
#undef self_sub_scal_
#undef scal_sub_self_1_
#undef scal_sub_self_2_
#undef scal_sub_self_3_
#undef scal_sub_self_
  !> @}
endmodule Data_Type_Riemann_Primitive1D
