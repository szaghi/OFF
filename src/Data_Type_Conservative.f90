!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_ConservativeDerivedType Data_Type_Conservative
!> Module definition of Type_Conservative
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_ConservativePublicProcedure Data_Type_Conservative
!> Module definition of Type_Conservative
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_ConservativePrivateProcedure Data_Type_Conservative
!> Module definition of Type_Conservative
!> @}

!> This module contains the definition of Type_Conservative and its procedures.
!> Type_Conservative is a derived type that handles conservative fluid dynamic variables.
!> @note The operators of assignment (=), multiplication (*), division (/), sum (+) and subtraction (-) have been overloaded.
!> Therefore this module provides a far-complete algebra based on Type_Conservative derived type. This algebra simplifies the
!> numerical integration of Partial Differential Equations (PDE) systems based on conservative formulation.
!> @todo \b DocComplete: Complete the documentation of internal procedures
module Data_Type_Conservative
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                        ! Integers and reals precision definition.
USE Data_Type_Vector, only: Type_Vector ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: array2consf
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing conservative variables.
!> @note This derived type can represent multi species fluids. The density component, \b rs, is a dynamic memory component defined
!> as an allocatable 1D array. \b rs is allocated at runtime with the number of initial species that constitute the initial fluid
!> mixture. Due to the presence of a dynamic component a freeing memory "method" for this component is necessary. Before deallocate
!> a variable defined as Type_Conservative the free function must be invoked to free the memory of the dynamic component.
!> @ingroup Data_Type_ConservativeDerivedType
type, public:: Type_Conservative
  real(R8P), allocatable:: rs(:)       !< Density of single species [1:Ns].
  type(Type_Vector)::      rv          !< Momentum vector.
  real(R8P)::              re = 0._R8P !< Product of density for specific total internal energy (sum(r)*E).
  contains
    procedure:: free                             ! Procedure for freeing dynamic memory.
    procedure:: alloc                            ! Procedure for allocating dynamic memory.
    procedure:: cons2array                       ! Procedure for converting derived type Type_Conservative to array.
    procedure:: array2cons                       ! Procedure for converting array to derived type Type_Conservative.
    procedure:: print => print_conservative_self ! Procedure for printing conservatives with a pretty format.
    final::     finalize                         ! Procedure for freeing dynamic memory when finalizing.
    ! operators overloading
#include 'Data_Type_Bounds_Proc_OpOverloading.inc'
endtype Type_Conservative
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_ConservativePublicProcedure
  !> @{
  !> Function for converting array to derived type Type_Conservative.
  pure function array2consf(array) result(cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN)::   array(:) !< Conservative data in the form of array.
  type(Type_Conservative):: cons     !< Derived type conservative data.
  integer(I_P)::            Ns       !< Number of species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(array)-4
  if (allocated(cons%rs)) deallocate(cons%rs) ; allocate(cons%rs(1:Ns))
  cons%rs   = array(1:Ns)
  cons%rv%x = array(Ns+1)
  cons%rv%y = array(Ns+2)
  cons%rv%z = array(Ns+3)
  cons%re   = array(Ns+4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction array2consf
  !> @}

  !> @ingroup Data_Type_ConservativePrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free(cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(INOUT):: cons !< Conservative data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons%rs)) deallocate(cons%rs)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc(cons,Ns)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(INOUT):: cons  !< Conservative initialized data.
  integer(I4P),             intent(IN)::    Ns    !< Number of species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call cons%free ; allocate(cons%rs(1:Ns)) ; cons%rs = 0._R8P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons !< Conservative data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call cons%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for converting derived type Type_Conservative to array.
  !> @return \b array real(R8P), dimension(1:size(cons\%rs)+4) variable.
  pure function cons2array(cons) result(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(IN):: cons                     !< Derived type conservative data.
  real(R8P)::                            array(1:size(cons%rs)+4) !< Conservative data in the form of array.
  integer(I_P)::                         Ns                       !< Number of species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(cons%rs)
  array(1:Ns) = cons%rs
  array(Ns+1) = cons%rv%x
  array(Ns+2) = cons%rv%y
  array(Ns+3) = cons%rv%z
  array(Ns+4) = cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons2array

  !> Subroutine for converting array to derived type Type_Conservative.
  pure subroutine array2cons(cons,array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(INOUT):: cons     !< Derived type conservative data.
  real(R8P),                intent(IN)::    array(:) !< Conservative data in the form of array.
  integer(I_P)::                            Ns       !< Number of species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(array)-4
  if (allocated(cons%rs)) deallocate(cons%rs) ; allocate(cons%rs(1:Ns))
  cons%rs   = array(1:Ns)
  cons%rv%x = array(Ns+1)
  cons%rv%y = array(Ns+2)
  cons%rv%z = array(Ns+3)
  cons%re   = array(Ns+4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array2cons

  !> @brief Procedure for printing in a pretty ascii format the components of type Type_Conservative.
  !> @return \b err integer(I_P) variable for error trapping.
  subroutine print_conservative_self(cons,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(IN)::  cons    !< Conservatives.
  character(*), optional,   intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional,   intent(OUT):: iostat  !< IO error.
  character(*), optional,   intent(OUT):: iomsg   !< IO error message.
  integer(I4P),             intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::         prefd   !< Prefixing string.
  integer(I4P)::                          iostatd !< IO error.
  character(500)::                        iomsgd  !< Temporary variable for IO error message.
  integer(I4P)::                          s       !< Species counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  do s= lbound(cons%rs,dim=1),ubound(cons%rs,dim=1)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' rs('//trim(str(.true.,s))//')='//str(n=cons%rs(s))
  enddo
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' rv(x)='//str(n=cons%rv%x)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' rv(y)='//str(n=cons%rv%y)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' rv(z)='//str(n=cons%rv%z)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' re   ='//str(n=cons%re  )
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_conservative_self

  ! Operators overloading.
  ! Operator (=)
#define self_type_ class(Type_Conservative), intent(INOUT):: self
#define ass_scal_ self%rs = real(scal,R8P) ; self%rv = real(scal,R8P) ; self%re = real(scal,R8P)
  !> @brief Procedure for assignment between two selfs.
  elemental subroutine assign_self(self1,self2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(INOUT):: self1
  type(Type_Conservative),  intent(IN)::    self2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self2%rs)) self1%rs = self2%rs
                           self1%rv = self2%rv
                           self1%re = self2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_self
#include 'Data_Type_Bounds_Proc_AssDefinitions.inc'
#undef self_type_
#undef ass_scal_

  ! Operator (*)
#define self_type_ class(Type_Conservative), intent(IN):: self
#define mul_type_ type(Type_Conservative):: mul
#define mul_alloc_ allocate(mul%rs(1:size(self%rs))) ;
#define mul_scal_ mul_alloc_ mul%rs = real(scal,R8P)*self%rs ; mul%rv = real(scal,R8P)*self%rv ; mul%re = real(scal,R8P)*self%re
  !> @brief Procedure for multiply (by components) two selfs.
  elemental function self_mul_self(self1,self2) result(mul)

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(IN):: self1
  type(Type_Conservative),  intent(IN):: self2
  type(Type_Conservative)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(self1%rs)))
  mul%rs = self1%rs * self2%rs
  mul%rv = self1%rv * self2%rv
  mul%re = self1%re * self2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_mul_self
#include 'Data_Type_Bounds_Proc_MulDefinitions.inc'
#undef self_type_
#undef mul_type_
#undef mul_alloc_
#undef mul_scal_

  ! Operator (/)
#define self_type_ class(Type_Conservative), intent(IN):: self
#define div_type_ type(Type_Conservative):: div
#define div_alloc_ allocate(div%rs(1:size(self%rs))) ;
#define div_scal_ div_alloc_ div%rs = self%rs/real(scal,R8P) ; div%rv = self%rv/real(scal,R8P) ; div%re = self%re/real(scal,R8P)
  !> @brief Procedure for divide self for self.
  elemental function self_div_self(self1,self2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(IN):: self1
  type(Type_Conservative),  intent(IN):: self2
  type(Type_Conservative)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%rs(1:size(self1%rs)))
  div%rs = self1%rs / self2%rs
  div%rv = self1%rv / self2%rv
  div%re = self1%re / self2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_div_self
#include 'Data_Type_Bounds_Proc_DivDefinitions.inc'
#undef self_type_
#undef div_type_
#undef div_alloc_
#undef div_scal_

  ! Operator (+)
#define self_type_ class(Type_Conservative), intent(IN):: self
#define summ_type_ type(Type_Conservative):: summ
#define summ_alloc_ allocate(summ%rs(1:size(self%rs))) ;
#define sum_scal_ summ_alloc_ summ%rs = real(scal,R8P)+self%rs ; summ%rv = real(scal,R8P)+self%rv ; summ%re = real(scal,R8P)+self%re
  !> @brief Procedure for applay unary + to a self.
  elemental function positive_self(self) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(IN):: self
  type(Type_Conservative)::              pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(pos%rs(1:size(self%rs)))
  pos%rs =  + self%rs
  pos%rv =  + self%rv
  pos%re =  + self%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_self

  !> @brief Procedure for sum self and self.
  elemental function self_sum_self(self1,self2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(IN):: self1
  type(Type_Conservative),  intent(IN):: self2
  type(Type_Conservative)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(self1%rs)))
  summ%rs = self1%rs + self2%rs
  summ%rv = self1%rv + self2%rv
  summ%re = self1%re + self2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sum_self
#include 'Data_Type_Bounds_Proc_SumDefinitions.inc'
#undef self_type_
#undef summ_type_
#undef summ_alloc_
#undef sum_scal_

  ! Operator (-)
#define self_type_ class(Type_Conservative), intent(IN):: self
#define sub_type_ type(Type_Conservative):: sub
#define sub_alloc_ allocate(sub%rs(1:size(self%rs))) ;
#define self_sub_scal_ sub_alloc_ sub%rs=self%rs-real(scal,R8P) ; sub%rv=self%rv-real(scal,R8P) ; sub%re=self%re-real(scal,R8P)
#define scal_sub_self_ sub_alloc_ sub%rs=real(scal,R8P)-self%rs ; sub%rv=real(scal,R8P)-self%rv ; sub%re=real(scal,R8P)-self%re
  !> @brief Procedure for applay unary - to a self.
  elemental function negative_self(self) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(IN):: self
  type(Type_Conservative)::              neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(neg%rs(1:size(self%rs)))
  neg%rs =  - self%rs
  neg%rv =  - self%rv
  neg%re =  - self%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_self

  !> @brief Procedure for subtract self and self.
  elemental function self_sub_self(self1,self2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(IN):: self1
  type(Type_Conservative),  intent(IN):: self2
  type(Type_Conservative)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(self1%rs)))
  sub%rs = self1%rs - self2%rs
  sub%rv = self1%rv - self2%rv
  sub%re = self1%re - self2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sub_self
#include 'Data_Type_Bounds_Proc_SubDefinitions.inc'
#undef self_type_
#undef sub_type_
#undef sub_alloc_
#undef self_sub_scal_
#undef scal_sub_self_
  !> @}
endmodule Data_Type_Conservative
