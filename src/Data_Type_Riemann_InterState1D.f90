!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Riemann_InterState1DDerivedType Data_Type_Riemann_InterState1D
!> Module definition of Riemann intermediate state in 1D, Type_Riemann_InterState1D
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Riemann_InterState1DPrivateProcedure Data_Type_Riemann_InterState1D
!> Module definition of Riemann intermediate state in 1D, Type_Riemann_InterState1D
!> @}

!> @brief Module Data_Type_Riemann_InterState1D contains the definition of Type_Riemann_InterState1D, that defines the Riemann
!> intermediate state (state 2 and 3) for 1D Riemann problem.
module Data_Type_Riemann_InterState1D
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                  ! Integers and reals precision definition.
USE Data_Type_Riemann_Primitive1D, only: Type_Riemann_Primitive1D ! Definition of Type_Riemann_Primitive1D.
USE Lib_Thermodynamic_Laws_Ideal,  only: E,H,a,r                  ! Library of thermodynamic laws for ideal calorically perfect gas.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_Riemann_InterState1D.
!> @ingroup Data_Type_Riemann_InterState1DDerivedType
type, public:: Type_Riemann_InterState1D
  real(R8P):: r2 = 0._R8P !< Density of state 2.
  real(R8P):: r3 = 0._R8P !< Density of state 3.
  real(R8P):: u  = 0._R8P !< Velocity (u23=S): contact discontinuity velocity.
  real(R8P):: p  = 0._R8P !< Pressure (p23).
  real(R8P):: S1 = 0._R8P !< Left-front of left fan.
  real(R8P):: S2 = 0._R8P !< Right-front of left fan.
  real(R8P):: S3 = 0._R8P !< Left-front of right fan.
  real(R8P):: S4 = 0._R8P !< Right-front of right fan.
  contains
    procedure:: set                   => set_self                   ! Procedure for setting self data.
    procedure:: compute_from_up       => compute_from_up_self       ! Procedure for computing intermediates from u and p.
    procedure:: compute_w14_from_u    => compute_w14_from_u_self    ! Procedure for computing wawes 1,4 from u.
    procedure:: compute_w14_from_up   => compute_w14_from_up_self   ! Procedure for computing wawes 1,4 from u and p.
    procedure:: compute_w1234_from_u  => compute_w1234_from_u_self  ! Procedure for computing wawes 1,2,3,4 from u.
    procedure:: compute_w1234_from_up => compute_w1234_from_up_self ! Procedure for computing wawes 1,2,3,4 from u and p.
    procedure:: print                 => print_self                 ! Procedure for printing self with a pretty format.
    ! operators overloading
#include "Data_Type_Bounds_Proc_OpOverloading.inc"
endtype Type_Riemann_InterState1D
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Riemann_InterState1DPrivateProcedure
  !> @{
  !> @brief Procedure for setting self data.
  elemental subroutine set_self(self,r2,r3,u,p,S1,S2,S3,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(INOUT):: self !< Riemann state.
  real(R8P), optional,              intent(IN)::    r2   !< Density of state 2.
  real(R8P), optional,              intent(IN)::    r3   !< Density of state 3.
  real(R8P), optional,              intent(IN)::    u    !< Velocity (u=S).
  real(R8P), optional,              intent(IN)::    p    !< Pressure (p).
  real(R8P), optional,              intent(IN)::    S1   !< Left-front of left fan.
  real(R8P), optional,              intent(IN)::    S2   !< Right-front of left fan.
  real(R8P), optional,              intent(IN)::    S3   !< Left-front of right fan.
  real(R8P), optional,              intent(IN)::    S4   !< Right-front of right fan.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(r2)) self%r2 = r2
  if (present(r3)) self%r3 = r3
  if (present(u )) self%u  = u
  if (present(p )) self%p  = p
  if (present(S1)) self%S1 = S1
  if (present(S2)) self%S2 = S2
  if (present(S3)) self%S3 = S3
  if (present(S4)) self%S4 = S4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_self

  !> @brief Procedure for computing intermediates states knowing the value of speed and pressure (u, p) of intermediates states.
  elemental subroutine compute_from_up_self(state23,state1,state4)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state (2,3).
  type(Type_Riemann_Primitive1D),   intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),   intent(IN)::    state4  !< State 4 (right).
  real(R8P)::                                       a2,a3   !< Speed of sound of state 2 and 3.
  real(R8P)::                                       gm1,gp1 !< g - 1, g + 1.
  real(R8P)::                                       p_p     !< Pressure ratio.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left state
  if (state23%u<state1%u) then
    ! shock
    gm1        = state1%g - 1._R8P
    gp1        = state1%g + 1._R8P
    p_p        = state23%p/state1%p
    a2         = state1%a*sqrt((gp1 + gm1*p_p)/(gp1 + gm1/p_p))
    state23%r2 = r(p=state23%p,a=a2,g=state1%g)
    state23%S1 = state1%u - state1%a*sqrt(1._R8P + 0.5_R8P*gp1/state1%g*(p_p-1._R8P))
    state23%S2 = state23%S1
  else
    ! rarefaction
    a2         = state1%a - 0.5_R8P*(state1%g-1._R8P)*(state23%u-state1%u)
    state23%r2 = r(p=state23%p,a=a2,g=state1%g)
    state23%S1 = state1%u  - state1%a
    state23%S2 = state23%u - a2
  endif
  ! computing right state
  if (state23%u>state4%u) then
    ! shock
    gm1        = state4%g - 1._R8P
    gp1        = state4%g + 1._R8P
    p_p        = state23%p/state4%p
    a3         = state4%a*sqrt((gp1 + gm1*p_p)/(gp1 + gm1/p_p))
    state23%r3 = r(p=state23%p,a=a3,g=state4%g)
    state23%S4 = state4%u + state4%a*sqrt(1._R8P + 0.5_R8P*gp1/state4%g*(p_p-1._R8P))
    state23%S3 = state23%S4
  else
    ! rarefaction
    a3         = state4%a + 0.5_R8P*(state4%g-1._R8P)*(state23%u-state4%u)
    state23%r3 = r(p=state23%p,a=a3,g=state4%g)
    state23%S4 = state4%u  + state4%a
    state23%S3 = state23%u + a3
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_from_up_self

  !> @brief Procedure for computing waves speed 1 and 4 knowing the value of speed (u) of intermediates states.
  elemental subroutine compute_w14_from_u_self(state23,state1,state4)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state (2,3).
  type(Type_Riemann_Primitive1D),   intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),   intent(IN)::    state4  !< State 4 (right).
  real(R8P)::                                       x       !< Dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left state
  if (state23%u<state1%u) then
    ! shock
    x  = 0.25_R8P*(state1%g + 1._R8P)*(state23%u-state1%u)/state1%a
    state23%S1 = state1%u + state1%a*(x - sqrt(1.0_R8P+x*x))
  else
    ! rarefaction
    state23%S1 = state1%u - state1%a
  endif
  ! computing right state
  if (state23%u>state4%u) then
    ! shock
    x  = 0.25_R8P*(state4%g + 1._R8P)*(state23%u-state4%u)/state4%a
    state23%S4 = state4%u + state4%a*(x + sqrt(1.0_R8P+x*x))
  else
    ! rarefaction
    state23%S4 = state4%u  + state4%a
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_w14_from_u_self

  !> @brief Subroutine for computing waves speed 1 and 4 knowing the value of speed and pressure (u, p) of intermediates states.
  elemental subroutine compute_w14_from_up_self(state23,state1,state4)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state (2,3).
  type(Type_Riemann_Primitive1D),   intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),   intent(IN)::    state4  !< State 4 (right).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left state
  if (state23%u<state1%u) then
    ! shock
    state23%S1 = state1%u - state1%a*sqrt(1._R8P + 0.5_R8P*(state1%g + 1._R8P)/state1%g*(state23%p/state1%p-1._R8P))
  else
    ! rarefaction
    state23%S1 = state1%u - state1%a
  endif
  ! computing right state
  if (state23%u>state4%u) then
    ! shock
    state23%S4 = state4%u + state4%a*sqrt(1._R8P + 0.5_R8P*(state4%g + 1._R8P)/state4%g*(state23%p/state4%p-1._R8P))
  else
    ! rarefaction
    state23%S4 = state4%u  + state4%a
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_w14_from_up_self

  !> @brief Procedure for computing waves speed knowing the value of speed (u) of intermediates states.
  elemental subroutine compute_w1234_from_u_self(state23,state1,state4)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state (2,3).
  type(Type_Riemann_Primitive1D),   intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),   intent(IN)::    state4  !< State 4 (right).
  real(R8P)::                                       x       !< Dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left state
  if (state23%u<state1%u) then
    ! shock
    x  = 0.25_R8P*(state1%g + 1._R8P)*(state23%u-state1%u)/state1%a
    state23%S1 = state1%u + state1%a*(x - sqrt(1.0_R8P+x*x))
    state23%S2 = state23%S1
  else
    ! rarefaction
    state23%S1 = state1%u  - state1%a
    state23%S2 = state23%u - state1%a - 0.5_R8P*(state1%g - 1._R8P)*(state23%u - state1%u)
  endif
  ! computing right state
  if (state23%u>state4%u) then
    ! shock
    x  = 0.25_R8P*(state4%g + 1._R8P)*(state23%u-state4%u)/state4%a
    state23%S4 = state4%u + state4%a*(x + sqrt(1.0_R8P+x*x))
    state23%S3 = state23%S4
  else
    ! rarefaction
    state23%S4 = state4%u  + state4%a
    state23%S3 = state23%u + state4%a + 0.5_R8P*(state4%g - 1._R8P)*(state23%u - state4%u)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_w1234_from_u_self

  !> @brief Procedure for computing waves speed knowing the value of speed and pressure (u, p) of intermediates states.
  elemental subroutine compute_w1234_from_up_self(state23,state1,state4)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state (2,3).
  type(Type_Riemann_Primitive1D),   intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),   intent(IN)::    state4  !< State 4 (right).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left state
  if (state23%u<state1%u) then
    ! shock
    state23%S1 = state1%u - state1%a*sqrt(1._R8P + 0.5_R8P*(state1%g + 1._R8P)/state1%g*(state23%p/state1%p-1._R8P))
    state23%S2 = state23%S1
  else
    ! rarefaction
    state23%S1 = state1%u  - state1%a
    state23%S2 = state23%u - state1%a - 0.5_R8P*(state1%g - 1._R8P)*(state23%u - state1%u)
  endif
  ! computing right state
  if (state23%u>state4%u) then
    ! shock
    state23%S4 = state4%u + state4%a*sqrt(1._R8P + 0.5_R8P*(state4%g + 1._R8P)/state4%g*(state23%p/state4%p-1._R8P))
    state23%S3 = state23%S4
  else
    ! rarefaction
    state23%S4 = state4%u  + state4%a
    state23%S3 = state23%u + state4%a + 0.5_R8P*(state4%g - 1._R8P)*(state23%u - state4%u)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_w1234_from_up_self

  !> @brief Procedure for printing self with a pretty format.
  subroutine print_self(self,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(IN)::  self    !< Riemann state.
  character(*), optional,           intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional,           intent(OUT):: iostat  !< IO error.
  character(*), optional,           intent(OUT):: iomsg   !< IO error message.
  integer(I4P),                     intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::                 prefd   !< Prefixing string.
  integer(I4P)::                                  iostatd !< IO error.
  character(500)::                                iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' r2='//str(n=self%r2)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' r3='//str(n=self%r3)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' u ='//str(n=self%u )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' p ='//str(n=self%p )
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' S1='//str(n=self%S1)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' S2='//str(n=self%S2)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' S3='//str(n=self%S3)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' S4='//str(n=self%S4)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_self

  ! Operators overloading.
  ! Operator (=)
#define self_type_ class(Type_Riemann_InterState1D), intent(INOUT):: self
#define ass_scal_1_ self%r2=real(scal,R8P);self%r3=real(scal,R8P);self%u=real(scal,R8P);self%p=real(scal,R8P);
#define ass_scal_2_ self%S1=real(scal,R8P);self%s2=real(scal,R8P);self%S3=real(scal,R8P);self%S4=real(scal,R8P)
#define ass_scal_ ass_scal_1_ ass_scal_2_
  !> @brief Procedure for assignment between two selfs.
  elemental subroutine assign_self(self1,self2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(INOUT):: self1
  type(Type_Riemann_InterState1D),  intent(IN)::    self2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self1%r2 = self2%r2
  self1%r3 = self2%r3
  self1%u  = self2%u
  self1%p  = self2%p
  self1%S1 = self2%S1
  self1%S2 = self2%S2
  self1%S3 = self2%S3
  self1%S4 = self2%S4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_self
#include "Data_Type_Bounds_Proc_AssDefinitions.inc"
#undef self_type_
#undef ass_scal_1_
#undef ass_scal_2_
#undef ass_scal_

  ! Operator (*)
#define self_type_ class(Type_Riemann_InterState1D), intent(IN):: self
#define mul_type_ type(Type_Riemann_InterState1D):: mul
#define mul_scal_1_ mul%r2=real(scal,R8P)*self%r2;mul%r3=real(scal,R8P)*self%r3;
#define mul_scal_2_ mul%u=real(scal,R8P)*self%u;mul%p=real(scal,R8P)*self%p;
#define mul_scal_3_ mul%S1=real(scal,R8P)*self%S1;mul%S2=real(scal,R8P)*self%S2;
#define mul_scal_4_ mul%S3=real(scal,R8P)*self%S3;mul%S4=real(scal,R8P)*self%S4
#define mul_scal_ mul_scal_1_ mul_scal_2_ mul_scal_3_ mul_scal_4_
  !> @brief Procedure for multiply (by components) two selfs.
  elemental function self_mul_self(self1,self2) result(mul)

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(IN):: self1
  type(Type_Riemann_InterState1D),  intent(IN):: self2
  type(Type_Riemann_InterState1D)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%r2 = self1%r2 * self2%r2
  mul%r3 = self1%r3 * self2%r3
  mul%u  = self1%u  * self2%u
  mul%p  = self1%p  * self2%p
  mul%S1 = self1%S1 * self2%S1
  mul%S2 = self1%S2 * self2%S2
  mul%S3 = self1%S3 * self2%S3
  mul%S4 = self1%S4 * self2%S4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_mul_self
#include "Data_Type_Bounds_Proc_MulDefinitions.inc"
#undef self_type_
#undef mul_type_
#undef mul_scal_1_
#undef mul_scal_2_
#undef mul_scal_3_
#undef mul_scal_4_
#undef mul_scal_

  ! Operator (/)
#define self_type_ class(Type_Riemann_InterState1D), intent(IN):: self
#define div_type_ type(Type_Riemann_InterState1D):: div
#define div_scal_1_ div%r2=self%r2/real(scal,R8P);div%r3=self%r3/real(scal,R8P);
#define div_scal_2_ div%u=self%u/real(scal,R8P);div%p=self%p/real(scal,R8P);
#define div_scal_3_ div%S1=self%S1/real(scal,R8P);div%S2=self%S2/real(scal,R8P);
#define div_scal_4_ div%S3=self%S3/real(scal,R8P);div%S4=self%S4/real(scal,R8P)
#define div_scal_ div_scal_1_ div_scal_2_ div_scal_3_ div_scal_4_
  !> @brief Procedure for divide self for self.
  elemental function self_div_self(self1,self2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(IN):: self1
  type(Type_Riemann_InterState1D),  intent(IN):: self2
  type(Type_Riemann_InterState1D)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%r2 = self1%r2 / self2%r2
  div%r3 = self1%r3 / self2%r3
  div%u  = self1%u  / self2%u
  div%p  = self1%p  / self2%p
  div%S1 = self1%S1 / self2%S1
  div%S2 = self1%S2 / self2%S2
  div%S3 = self1%S3 / self2%S3
  div%S4 = self1%S4 / self2%S4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_div_self
#include "Data_Type_Bounds_Proc_DivDefinitions.inc"
#undef self_type_
#undef div_type_
#undef div_scal_1_
#undef div_scal_2_
#undef div_scal_3_
#undef div_scal_4_
#undef div_scal_

  ! Operator (+)
#define self_type_ class(Type_Riemann_InterState1D), intent(IN):: self
#define summ_type_ type(Type_Riemann_InterState1D):: summ
#define sum_scal_1_ summ%r2=self%r2+real(scal,R8P);summ%r3=self%r3+real(scal,R8P);
#define sum_scal_2_ summ%u=self%u+real(scal,R8P);summ%p=self%p+real(scal,R8P);
#define sum_scal_3_ summ%S1=self%S1+real(scal,R8P);summ%S2=self%S2+real(scal,R8P);
#define sum_scal_4_ summ%S3=self%S3+real(scal,R8P);summ%S4=self%S4+real(scal,R8P)
#define sum_scal_ sum_scal_1_ sum_scal_2_ sum_scal_3_ sum_scal_4_
  !> @brief Procedure for applay unary + to a self.
  elemental function positive_self(self) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(IN):: self
  type(Type_Riemann_InterState1D)::              pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pos%r2 = + self%r2
  pos%r3 = + self%r3
  pos%u  = + self%u
  pos%p  = + self%p
  pos%S1 = + self%S1
  pos%S2 = + self%S2
  pos%S3 = + self%S3
  pos%S4 = + self%S4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_self

  !> @brief Procedure for sum self and self.
  elemental function self_sum_self(self1,self2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(IN):: self1
  type(Type_Riemann_InterState1D),  intent(IN):: self2
  type(Type_Riemann_InterState1D)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%r2 = self1%r2 + self2%r2
  summ%r3 = self1%r3 + self2%r3
  summ%u  = self1%u  + self2%u
  summ%p  = self1%p  + self2%p
  summ%S1 = self1%S1 + self2%S1
  summ%S2 = self1%S2 + self2%S2
  summ%S3 = self1%S3 + self2%S3
  summ%S4 = self1%S4 + self2%S4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sum_self
#include "Data_Type_Bounds_Proc_SumDefinitions.inc"
#undef self_type_
#undef summ_type_
#undef sum_scal_1_
#undef sum_scal_2_
#undef sum_scal_3_
#undef sum_scal_4_
#undef sum_scal_

  ! Operator (-)
#define self_type_ class(Type_Riemann_InterState1D), intent(IN):: self
#define sub_type_ type(Type_Riemann_InterState1D):: sub
#define self_sub_scal_1_ sub%r2=self%r2-real(scal,R8P);sub%r3=self%r3-real(scal,R8P);
#define self_sub_scal_2_ sub%u=self%u-real(scal,R8P);sub%p=self%p-real(scal,R8P);
#define self_sub_scal_3_ sub%S1=self%S1-real(scal,R8P);sub%S2=self%S2-real(scal,R8P);
#define self_sub_scal_4_ sub%S3=self%S3-real(scal,R8P);sub%S4=self%S4-real(scal,R8P)
#define scal_sub_self_1_ sub%r2=real(scal,R8P)-self%r2;sub%r3=real(scal,R8P)-self%r3;
#define scal_sub_self_2_ sub%u=real(scal,R8P)-self%u;sub%p=real(scal,R8P)-self%p ;
#define scal_sub_self_3_ sub%S1=real(scal,R8P)-self%S1;sub%S2=real(scal,R8P)-self%S2;
#define scal_sub_self_4_ sub%S3=real(scal,R8P)-self%S3;sub%S4=real(scal,R8P)-self%S4
#define self_sub_scal_ self_sub_scal_1_ self_sub_scal_2_ self_sub_scal_3_ self_sub_scal_4_
#define scal_sub_self_ scal_sub_self_1_ scal_sub_self_2_ scal_sub_self_3_ scal_sub_self_4_
  !> @brief Procedure for applay unary - to a self.
  elemental function negative_self(self) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(IN):: self
  type(Type_Riemann_InterState1D)::              neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neg%r2 = - self%r2
  neg%r3 = - self%r3
  neg%u  = - self%u
  neg%p  = - self%p
  neg%S1 = - self%S1
  neg%S2 = - self%S2
  neg%S3 = - self%S3
  neg%S4 = - self%S4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_self

  !> @brief Procedure for subtract self and self.
  elemental function self_sub_self(self1,self2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Riemann_InterState1D), intent(IN):: self1
  type(Type_Riemann_InterState1D),  intent(IN):: self2
  type(Type_Riemann_InterState1D)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%r2 = self1%r2 - self2%r2
  sub%r3 = self1%r3 - self2%r3
  sub%u  = self1%u  - self2%u
  sub%p  = self1%p  - self2%p
  sub%S1 = self1%S1 - self2%S1
  sub%S2 = self1%S2 - self2%S2
  sub%S3 = self1%S3 - self2%S3
  sub%S4 = self1%S4 - self2%S4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sub_self
#include "Data_Type_Bounds_Proc_SubDefinitions.inc"
#undef self_type_
#undef sub_type_
#undef self_sub_scal_1_
#undef self_sub_scal_2_
#undef self_sub_scal_3_
#undef self_sub_scal_4_
#undef self_sub_scal_
#undef scal_sub_self_1_
#undef scal_sub_self_2_
#undef scal_sub_self_3_
#undef scal_sub_self_4_
#undef scal_sub_self_
  !> @}
endmodule Data_Type_Riemann_InterState1D
