!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Primitive1DDerivedType Data_Type_Primitive1D
!> Module definition of Type_Primitive1D
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_Primitive1DInterface Data_Type_Primitive1D
!> Module definition of Type_Primitive1D
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Primitive1DPrivateProcedure Data_Type_Primitive1D
!> Module definition of Type_Primitive1D
!> @}

!> This module contains the definition of Type_Primitive1D and its procedures.
!> Type_Primitive1D is a derived type that handles primitive fluid dynamic variables in 1D dimension.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Data_Type_Primitive1DDerivedType
module Data_Type_Primitive1D
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                     ! Integers and reals precision definition.
USE Data_Type_Species,            only: Type_Species ! Definition of Type_Species.
USE Lib_Thermodynamic_Laws_Ideal, only: a,E,H        ! Library for thermodynamic laws for ideal calorically perfect gas.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing primitive variables.
!> @note This derived type can represent multi species fluids. The density component, \b r, is a dynamic memory component defined
!> as an allocatable 1D array. \b r is allocated at runtime with the number of initial species that constitute the initial fluid
!> mixture. Due to the presence of a dynamic component a freeing memory "method" for this component is necessary. Before deallocate
!> a variable defined as Type_Primitive1D the free function must be invoked to free the memory of the dynamic component.
!> @note It is worth noting that the component r, u and p are independent whereas d and g are dependent from the previous ones.
!> Consequently, there are two procedures, namely compute_d and compute_g, for their computation from the independent components.
!> @ingroup Data_Type_Primitive1DDerivedType
type, public:: Type_Primitive1D
  real(R8P), allocatable:: r(:)       !< Density of single species [1:Ns].
  real(R8P)::              u = 0._R8P !< Velocity.
  real(R8P)::              p = 0._R8P !< Pressure.
  real(R8P)::              d = 0._R8P !< Density = sum(r(1:Ns)).
  real(R8P)::              g = 0._R8P !< Specific heats ratio \f$ \gamma = \frac{c_p}{c_v} \f$; it depends on partial densities
                                      !< r and on the initial species specific heats.
  contains
    procedure:: free                          ! Procedure for freeing dynamic memory.
    procedure:: alloc                         ! Procedure for allocating dynamic memory.
    procedure:: compute_d                     ! Procedure for computing the density from the densities of single species.
    procedure:: compute_g                     ! Procedure for computing the specific heats ratio.
    procedure:: prim2array                    ! Procedure for converting derived type Type_Primitive to array.
    procedure:: array2prim                    ! Procedure for converting array to derived type Type_Primitive.
    procedure:: a => a_self                   ! Procedure for computing the speed of sound.
    procedure:: E => E_self                   ! Procedure for computing the total specific energy.
    procedure:: H => H_self                   ! Procedure for computing the total specific entalpy.
    procedure:: eigenvectL                    ! Procedure for computing the left eigenvectors matrix.
    procedure:: eigenvectR                    ! Procedure for computing the right eigenvectors matrix.
    procedure:: print => print_primitive_self ! Procedure for printing primitives with a pretty format.
    final::     finalize                      ! Procedure for freeing dynamic memory when finalizing.
    ! operators overloading
#include "Data_Type_Bounds_Proc_OpOverloading.inc"
endtype Type_Primitive1D
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Primitive1DPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free(prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(INOUT):: prim !< Primitive data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim%r)) deallocate(prim%r)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive1D), intent(INOUT):: prim !< Primitive data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call prim%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc(prim,Ns)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(INOUT):: prim  !< Primitive data.
  integer(I4P),            intent(IN)::    Ns    !< Number of species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call prim%free ; allocate(prim%r(1:Ns)) ; prim%r = 0._R8P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc

  !> @brief Procedure for computing the density from the densities of single species.
  elemental subroutine compute_d(prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(INOUT):: prim !< Primitive data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%d = sum(prim%r(:))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_d

  !> @brief Procedure for computing the specific heats ratio.
  elemental subroutine compute_g(prim,species0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(INOUT):: prim     !< Primitive data.
  type(Type_Species),      intent(IN)::    species0 !< Initial species.
  real(R8P), allocatable::                 c(:)     !< Species concentration.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  c = prim%r/prim%d ; prim%g = dot_product(c,species0%heats(:)%cp)/dot_product(c,species0%heats(:)%cv)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_g

  !> @brief Procedure for converting derived type Type_Primitive1D to array.
  pure function prim2array(prim) result(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN):: prim                    !< Derived type primitive data.
  real(R8P)::                           array(1:size(prim%r)+4) !< Primitive data in the form of array.
  integer(I4P)::                        Ns                      !< Number of species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(prim%r)
  array(1:Ns) = prim%r
  array(Ns+1) = prim%u
  array(Ns+2) = prim%p
  array(Ns+3) = prim%d
  array(Ns+4) = prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim2array

  !> @brief Procedure for converting array to derived type Type_Primitive1D.
  pure subroutine array2prim(prim,array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(INOUT):: prim     !< Derived type primitive data.
  real(R8P),               intent(IN)::    array(:) !< Primitive data in the form of array.
  integer(I4P)::                           Ns       !< Number of species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(array)-4
  if (allocated(prim%r)) deallocate(prim%r)  ; allocate(prim%r(1:Ns))
  prim%r = array(1:Ns)
  prim%u = array(Ns+1)
  prim%p = array(Ns+2)
  prim%d = array(Ns+3)
  prim%g = array(Ns+4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array2prim

  !> @brief Procedure for computing the speed of sound.
  elemental function a_self(prim) result(ss)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN):: prim !< Primitive data.
  real(R8P)::                           ss   !< Speed of sound.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ss = a(p=prim%p,r=prim%d,g=prim%g)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction a_self

  !> @brief Procedure for computing the total specific energy.
  elemental function E_self(prim) result(energy)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN):: prim   !< Primitive data.
  real(R8P)::                           energy !< Total specific energy.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  energy = E(p=prim%p,r=prim%d,u=prim%u,g=prim%g)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction E_self

  !> @brief Procedure for computing the total specific entalpy.
  elemental function H_self(prim) result(entalpy)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN):: prim    !< Primitive data.
  real(R8P)::                           entalpy !< Total specific entalpy.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  entalpy = H(p=prim%p,r=prim%d,u=prim%u,g=prim%g)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction H_self

  !> @brief Procedure for computing the left eigenvectors matrix.
  pure function eigenvectL(prim) result(L)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN)::  prim                                 !< Primitive data.
  real(R8P)::                            L(1:size(prim%r)+4,1:size(prim%r)+4) !< Left eigenvectors matrix.
  real(R8P)::                            gp                                   !< prim%g*prim%p.
  real(R8P)::                            a                                    !< Speed of sound.
  integer(I4P)::                         s,Ns                                 !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing
  Ns = size(prim%r,dim=1)
  gp = prim%g*prim%p
  a  = sqrt(gp/prim%d)
  L  = 0._R8P
  ! assigning non-zero values of L
                       L(1,   Ns+1) = -gp/a           ; L(1,   Ns+2) =  1._R8P
  do s=2,Ns+1
    if (prim%r(s-1)>0) L(s,   s-1 ) =  gp/prim%r(s-1) ; L(s,   Ns+2) = -1._R8P
  enddo
                       L(Ns+2,Ns+1) =  gp/a           ; L(Ns+2,Ns+2) =  1._R8P
    !L(Ns+3:Ns+4,:) = 1._R8P
    !L(:,Ns+3:Ns+4) = 1._R8P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction eigenvectL

  !> @brief Procedure for computing the right eigenvectors matrix.
  pure function eigenvectR(prim) result(R)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN)::  prim                                 !< Primitive data.
  real(R8P)::                            R(1:size(prim%r)+4,1:size(prim%r)+4) !< Right eigenvectors matrix.
  real(R8P)::                            gp                                   !< prim%g*prim%p.
  real(R8P)::                            gp_inv                               !< 1/(gp).
  real(R8P)::                            a                                    !< Speed of sound.
  integer(I4P)::                         s,Ns                                 !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing
  Ns     = size(prim%r,dim=1)
  gp     = prim%g*prim%p
  a      = sqrt(gp/prim%d)
  gp_inv = 1._R8P/gp
  R      = 0._R8P
  ! assigning non-zero values of R
  do s=1,Ns
    R(s,   1) =  0.5_R8P*prim%r(s)*gp_inv ; R(s,s+1) = prim%r(s)*gp_inv ; R(s,   Ns+2) = R(s,1)
  enddo
    R(Ns+1,1) = -0.5_R8P*a*gp_inv ;                                       R(Ns+1,Ns+2) = 0.5_R8P*a*gp_inv
    R(Ns+2,1) =  0.5_R8P          ;                                       R(Ns+2,Ns+2) = 0.5_R8P
    !R(Ns+3:Ns+4,:) = 1._R8P
    !R(:,Ns+3:Ns+4) = 1._R8P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction eigenvectR

  !> @brief Procedure for printing primitives with a pretty format.
  subroutine print_primitive_self(prim,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D),  intent(IN)::  prim    !< Primitives.
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
  do s=lbound(prim%r,dim=1),ubound(prim%r,dim=1)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' r('//trim(str(.true.,s))//')='//str(n=prim%r(s))
  enddo
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' u   ='//str(n=prim%u)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' p   ='//str(n=prim%p)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' d   ='//str(n=prim%d)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' g   ='//str(n=prim%g)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_primitive_self

  ! Operators overloading.
  ! Operator (=)
#define self_type_ class(Type_Primitive1D), intent(INOUT):: self
#define ass_scal_ self%r=real(scal,R8P);self%u=real(scal,R8P);self%p=real(scal,R8P);self%d=real(scal,R8P);self%g=real(scal,R8P)
  !> @brief Procedure for assignment between two selfs.
  elemental subroutine assign_self(self1,self2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(INOUT):: self1
  type(Type_Primitive1D),  intent(IN)::    self2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self1%r = self2%r
  self1%u = self2%u
  self1%p = self2%p
  self1%d = self2%d
  self1%g = self2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_self
#include "Data_Type_Bounds_Proc_AssDefinitions.inc"
#undef self_type_
#undef ass_scal_

  ! Operator (*)
#define self_type_ class(Type_Primitive1D), intent(IN):: self
#define mul_type_ type(Type_Primitive1D):: mul
#define mul_alloc_ allocate(mul%r(1:size(self%r))) ;
#define mul_scal_1_ mul%r=real(scal,R8P)*self%r;mul%u=real(scal,R8P)*self%u;mul%p=real(scal,R8P)*self%p;
#define mul_scal_2_ mul%d=real(scal,R8P)*self%d;mul%g=real(scal,R8P)*self%g
#define mul_scal_ mul_alloc_ mul_scal_1_ mul_scal_2_
  !> @brief Procedure for multiply (by components) two selfs.
  elemental function self_mul_self(self1,self2) result(mul)

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN):: self1
  type(Type_Primitive1D),  intent(IN):: self2
  type(Type_Primitive1D)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(self1%r)))
  mul%r = self1%r * self2%r
  mul%u = self1%u * self2%u
  mul%p = self1%p * self2%p
  mul%d = self1%d * self2%d
  mul%g = self1%g * self2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_mul_self
#include "Data_Type_Bounds_Proc_MulDefinitions.inc"
#undef self_type_
#undef mul_type_
#undef mul_alloc_
#undef mul_scal_1_
#undef mul_scal_2_
#undef mul_scal_

  ! Operator (/)
#define self_type_ class(Type_Primitive1D), intent(IN):: self
#define div_type_ type(Type_Primitive1D):: div
#define div_alloc_ allocate(div%r(1:size(self%r))) ;
#define div_scal_1_ div%r=self%r/real(scal,R8P);div%u=self%u/real(scal,R8P);div%p=self%p/real(scal,R8P);
#define div_scal_2_ div%d=self%d/real(scal,R8P);div%g=self%g/real(scal,R8P)
#define div_scal_ div_alloc_ div_scal_1_ div_scal_2_
  !> @brief Procedure for divide self for self.
  elemental function self_div_self(self1,self2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN):: self1
  type(Type_Primitive1D),  intent(IN):: self2
  type(Type_Primitive1D)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%r(1:size(self1%r)))
  div%r = self1%r / self2%r
  div%u = self1%u / self2%u
  div%p = self1%p / self2%p
  div%d = self1%d / self2%d
  div%g = self1%g / self2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_div_self
#include "Data_Type_Bounds_Proc_DivDefinitions.inc"
#undef self_type_
#undef div_type_
#undef div_alloc_
#undef div_scal_1_
#undef div_scal_2_
#undef div_scal_

  ! Operator (+)
#define self_type_ class(Type_Primitive1D), intent(IN):: self
#define summ_type_ type(Type_Primitive1D):: summ
#define summ_alloc_ allocate(summ%r(1:size(self%r))) ;
#define sum_scal_1_ summ%r=real(scal,R8P)+self%r;summ%u=real(scal,R8P)+self%u;summ%p=real(scal,R8P)+self%p;
#define sum_scal_2_ summ%d=real(scal,R8P)+self%d;summ%g=real(scal,R8P)+self%g
#define sum_scal_ summ_alloc_ sum_scal_1_ sum_scal_2_
  !> @brief Procedure for applay unary + to a self.
  elemental function positive_self(self) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN):: self
  type(Type_Primitive1D)::              pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(pos%r(1:size(self%r)))
  pos%r =  + self%r
  pos%u =  + self%u
  pos%p =  + self%p
  pos%d =  + self%d
  pos%g =  + self%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_self

  !> @brief Procedure for sum self and self.
  elemental function self_sum_self(self1,self2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN):: self1
  type(Type_Primitive1D),  intent(IN):: self2
  type(Type_Primitive1D)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(self1%r)))
  summ%r = self1%r + self2%r
  summ%u = self1%u + self2%u
  summ%p = self1%p + self2%p
  summ%d = self1%d + self2%d
  summ%g = self1%g + self2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sum_self
#include "Data_Type_Bounds_Proc_SumDefinitions.inc"
#undef self_type_
#undef summ_type_
#undef summ_alloc_
#undef sum_scal_1_
#undef sum_scal_2_
#undef sum_scal_

  ! Operator (-)
#define self_type_ class(Type_Primitive1D), intent(IN):: self
#define sub_type_ type(Type_Primitive1D):: sub
#define sub_alloc_ allocate(sub%r(1:size(self%r))) ;
#define self_sub_scal_1_ sub%r=self%r-real(scal,R8P);sub%u=self%u-real(scal,R8P);sub%p=self%p-real(scal,R8P);
#define self_sub_scal_2_ sub%d=self%d-real(scal,R8P);sub%g=self%g-real(scal,R8P)
#define scal_sub_self_1_ sub%r=real(scal,R8P)-self%r;sub%u=real(scal,R8P)-self%u;sub%p=real(scal,R8P)-self%p;
#define scal_sub_self_2_ sub%d=real(scal,R8P)-self%d;sub%g=real(scal,R8P)-self%g
#define self_sub_scal_ sub_alloc_ self_sub_scal_1_ self_sub_scal_2_
#define scal_sub_self_ sub_alloc_ scal_sub_self_1_ scal_sub_self_2_
  !> @brief Procedure for applay unary - to a self.
  elemental function negative_self(self) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN):: self
  type(Type_Primitive1D)::              neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(neg%r(1:size(self%r)))
  neg%r =  - self%r
  neg%u =  - self%u
  neg%p =  - self%p
  neg%d =  - self%d
  neg%g =  - self%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_self

  !> @brief Procedure for subtract self and self.
  elemental function self_sub_self(self1,self2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive1D), intent(IN):: self1
  type(Type_Primitive1D),  intent(IN):: self2
  type(Type_Primitive1D)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(self1%r)))
  sub%r = self1%r - self2%r
  sub%u = self1%u - self2%u
  sub%p = self1%p - self2%p
  sub%d = self1%d - self2%d
  sub%g = self1%g - self2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sub_self
#include "Data_Type_Bounds_Proc_SubDefinitions.inc"
#undef self_type_
#undef sub_type_
#undef sub_alloc_
#undef self_sub_scal_1_
#undef self_sub_scal_2_
#undef self_sub_scal_
#undef scal_sub_self_1_
#undef scal_sub_self_2_
#undef scal_sub_self_
  !> @}
endmodule Data_Type_Primitive1D
