!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_PrimitiveDerivedType Data_Type_Primitive
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_PrimitiveInterface Data_Type_Primitive
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_PrimitivePublicProcedure Data_Type_Primitive
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_PrimitivePrivateProcedure Data_Type_Primitive
!> @}

!> This module contains the definition of Type_Primitive and its procedures.
!> Type_Primitive is a derived type that handles primitive fluid dynamic variables.
!> @note The operators of assignment (=), multiplication (*), division (/), sum (+) and subtraction (-) have been overloaded.
!> Therefore this module provides a far-complete algebra based on Type_Primitive derived type. This algebra simplifies the
!> numerical integration of Partial Differential Equations (PDE) systems based on primitive formulation.
!> @todo \b DocComplete: Complete the documentation of internal procedures
module Data_Type_Primitive
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision     ! Integers and reals precision definition.
USE Data_Type_Vector ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: write_primitive,read_primitive
public:: assignment (=)
public:: operator (*)
public:: operator (/)
public:: operator (+)
public:: operator (-)
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing primitive variables.
!> @note This derived type can represent multi species fluids. The density component, \b r, is a dynamic memory component defined
!> as an allocatable 1D array. \b r is allocated at runtime with the number of initial species that constitute the initial fluid
!> mixture. Due to the presence of a dynamic component a freeing memory "method" for this component is necessary. Before deallocate
!> a variable defined as Type_Primitive the free function must be invoked to free the memory of the dynamic component.
!> @ingroup Data_Type_PrimitiveDerivedType
type, public:: Type_Primitive
  real(R8P), allocatable:: r(:)       !< Density of single species [1:Ns].
  type(Type_Vector)::      v          !< Velocity vector.
  real(R8P)::              p = 0._R8P !< Pressure.
  real(R8P)::              d = 0._R8P !< Density = sum(r(1:Ns)).
  real(R8P)::              g = 0._R8P !< Specific heats ratio \f$ \gamma = \frac{c_p}{c_v} \f$.
  contains
    procedure, non_overridable:: init       ! Procedure for initilizing allocatable variables.
    procedure, non_overridable:: free       ! Procedure for freeing the memory of allocatable variables.
    procedure, non_overridable:: prim2array ! Procedure for converting derived type Type_Primitive to array.
    procedure, non_overridable:: array2prim ! Procedure for converting array to derived type Type_Primitive.
    procedure, non_overridable:: pprint     ! Procedure for printing Type_Primitive components with a "pretty" format.
endtype Type_Primitive
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Assignment operator (=) overloading.
!> @ingroup Data_Type_PrimitiveInterface
interface assignment (=)
  module procedure assign_prim
#ifdef r16p
  module procedure assign_ScalR16P
#endif
  module procedure assign_ScalR8P
  module procedure assign_ScalR4P
  module procedure assign_ScalI8P
  module procedure assign_ScalI4P
  module procedure assign_ScalI2P
  module procedure assign_ScalI1P
endinterface
!> @brief Multiplication operator (*) overloading.
!> @note The admissible multiplications are:
!>       - Type_Primitive * Type_Primitive: each component of first primitive variable (prim1) is multiplied for the
!>         corresponding component of the second one (prim2), i.e. \n
!>         \f$ {\rm result\%r = prim1\%r*prim2\%r} \f$ \n
!>         \f$ {\rm result\%v = prim1\%v*prim2\%v} \f$ \n
!>         \f$ {\rm result\%p = prim1\%p*prim2\%p} \f$ \n
!>         \f$ {\rm result\%d = prim1\%d*prim2\%d} \f$ \n
!>         \f$ {\rm result\%g = prim1\%g*prim2\%g} \f$ \n
!>       - scalar number (real or integer of any kinds defined in IR_Precision module) * Type_Primitive: each component of
!>         Type_Primitive is multiplied for the scalar, i.e. \n
!>         \f$ {\rm result\%r = prim\%r*scalar} \f$ \n
!>         \f$ {\rm result\%v = prim\%v*scalar} \f$ \n
!>         \f$ {\rm result\%p = prim\%p*scalar} \f$ \n
!>         \f$ {\rm result\%d = prim\%d*scalar} \f$ \n
!>         \f$ {\rm result\%g = prim\%g*scalar} \f$ \n
!>       - Type_Primitive * scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Primitive is multiplied for the scalar, i.e. \n
!>         \f$ {\rm result\%r = prim\%r*scalar} \f$ \n
!>         \f$ {\rm result\%v = prim\%v*scalar} \f$ \n
!>         \f$ {\rm result\%p = prim\%p*scalar} \f$ \n
!>         \f$ {\rm result\%d = prim\%d*scalar} \f$ \n
!>         \f$ {\rm result\%g = prim\%g*scalar} \f$ \n
!> @ingroup Data_Type_PrimitiveInterface
interface operator (*)
  module procedure prim_mul_prim
#ifdef r16p
  module procedure ScalR16P_mul_prim
#endif
  module procedure ScalR8P_mul_prim
  module procedure ScalR4P_mul_prim
  module procedure ScalI8P_mul_prim
  module procedure ScalI4P_mul_prim
  module procedure ScalI2P_mul_prim
  module procedure ScalI1P_mul_prim
#ifdef r16p
  module procedure prim_mul_ScalR16P
#endif
  module procedure prim_mul_ScalR8P
  module procedure prim_mul_ScalR4P
  module procedure prim_mul_ScalI8P
  module procedure prim_mul_ScalI4P
  module procedure prim_mul_ScalI2P
  module procedure prim_mul_ScalI1P
endinterface
!> @brief Division operator (/) overloading.
!> @note The admissible divisions are:
!>       - Type_Primitive / Type_Primitive: each component of first primitive variable (prim1) is divided for the
!>         corresponding component of the second one (prim2), i.e. \n
!>         \f$ {\rm result\%r = \frac{prim1\%r}{prim2\%r}} \f$ \n
!>         \f$ {\rm result\%v = \frac{prim1\%v}{prim2\%v}} \f$ \n
!>         \f$ {\rm result\%p = \frac{prim1\%p}{prim2\%p}} \f$ \n
!>         \f$ {\rm result\%d = \frac{prim1\%d}{prim2\%d}} \f$ \n
!>         \f$ {\rm result\%g = \frac{prim1\%g}{prim2\%g}} \f$ \n
!>       - Type_Primitive / scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Primitive is divided for the scalar, i.e. \n
!>         \f$ {\rm result\%r = \frac{prim\%r}{scalar}} \f$ \n
!>         \f$ {\rm result\%v = \frac{prim\%v}{scalar}} \f$ \n
!>         \f$ {\rm result\%p = \frac{prim\%p}{scalar}} \f$ \n
!>         \f$ {\rm result\%d = \frac{prim\%d}{scalar}} \f$ \n
!>         \f$ {\rm result\%g = \frac{prim\%g}{scalar}} \f$ \n
!> @ingroup Data_Type_PrimitiveInterface
interface operator (/)
  module procedure prim_div_prim
#ifdef r16p
  module procedure prim_div_ScalR16P
#endif
  module procedure prim_div_ScalR8P
  module procedure prim_div_ScalR4P
  module procedure prim_div_ScalI8P
  module procedure prim_div_ScalI4P
  module procedure prim_div_ScalI2P
  module procedure prim_div_ScalI1P
endinterface
!> @brief Sum operator (+) overloading.
!> @note The admissible summations are:
!>       - Type_Primitive + Type_Primitive: each component of first primitive variable (prim1) is summed with the
!>         corresponding component of the second one (prim2), i.e. \n
!>         \f$ {\rm result\%r = prim1\%r+prim2\%r} \f$ \n
!>         \f$ {\rm result\%v = prim1\%v+prim2\%v} \f$ \n
!>         \f$ {\rm result\%p = prim1\%p+prim2\%p} \f$ \n
!>         \f$ {\rm result\%d = prim1\%d+prim2\%d} \f$ \n
!>         \f$ {\rm result\%g = prim1\%g+prim2\%g} \f$ \n
!>       - scalar number (real or integer of any kinds defined in IR_Precision module) + Type_Primitive: each component of
!>         Type_Primitive is summed with the scalar, i.e. \n
!>         \f$ {\rm result\%r = prim\%r+scalar} \f$ \n
!>         \f$ {\rm result\%v = prim\%v+scalar} \f$ \n
!>         \f$ {\rm result\%p = prim\%p+scalar} \f$ \n
!>         \f$ {\rm result\%d = prim\%d+scalar} \f$ \n
!>         \f$ {\rm result\%g = prim\%g+scalar} \f$ \n
!>       - Type_Primitive + scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Primitive is summed with the scalar, i.e. \n
!>         \f$ {\rm result\%r = prim\%r+scalar} \f$ \n
!>         \f$ {\rm result\%v = prim\%v+scalar} \f$ \n
!>         \f$ {\rm result\%p = prim\%p+scalar} \f$ \n
!>         \f$ {\rm result\%d = prim\%d+scalar} \f$ \n
!>         \f$ {\rm result\%g = prim\%g+scalar} \f$ \n
!> @ingroup Data_Type_PrimitiveInterface
interface operator (+)
  module procedure positive_prim
  module procedure prim_sum_prim
#ifdef r16p
  module procedure ScalR16P_sum_prim
#endif
  module procedure ScalR8P_sum_prim
  module procedure ScalR4P_sum_prim
  module procedure ScalI8P_sum_prim
  module procedure ScalI4P_sum_prim
  module procedure ScalI2P_sum_prim
  module procedure ScalI1P_sum_prim
#ifdef r16p
  module procedure prim_sum_ScalR16P
#endif
  module procedure prim_sum_ScalR8P
  module procedure prim_sum_ScalR4P
  module procedure prim_sum_ScalI8P
  module procedure prim_sum_ScalI4P
  module procedure prim_sum_ScalI2P
  module procedure prim_sum_ScalI1P
endinterface
!> @brief Subtraction operator (-) overloading.
!> @note The admissible subtractions are:
!>       - Type_Primitive - Type_Primitive: each component of first primitive variable (prim1) is subtracted with the
!>         corresponding component of the second one (prim2), i.e. \n
!>         \f$ {\rm result\%r = prim1\%r-prim2\%r} \f$ \n
!>         \f$ {\rm result\%v = prim1\%v-prim2\%v} \f$ \n
!>         \f$ {\rm result\%p = prim1\%p-prim2\%p} \f$ \n
!>         \f$ {\rm result\%d = prim1\%d-prim2\%d} \f$ \n
!>         \f$ {\rm result\%g = prim1\%g-prim2\%g} \f$ \n
!>       - scalar number (real or integer of any kinds defined in IR_Precision module) - Type_Primitive: each component of
!>         Type_Primitive is subtracted with the scalar, i.e. \n
!>         \f$ {\rm result\%r = scalar-prim\%r} \f$ \n
!>         \f$ {\rm result\%v = scalar-prim\%v} \f$ \n
!>         \f$ {\rm result\%p = scalar-prim\%p} \f$ \n
!>         \f$ {\rm result\%d = scalar-prim\%d} \f$ \n
!>         \f$ {\rm result\%g = scalar-prim\%g} \f$ \n
!>       - Type_Primitive - scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Primitive is subtracted with the scalar, i.e. \n
!>         \f$ {\rm result\%r = prim\%r-scalar} \f$ \n
!>         \f$ {\rm result\%v = prim\%v-scalar} \f$ \n
!>         \f$ {\rm result\%p = prim\%p-scalar} \f$ \n
!>         \f$ {\rm result\%d = prim\%d-scalar} \f$ \n
!>         \f$ {\rm result\%g = prim\%g-scalar} \f$ \n
!> @ingroup Data_Type_PrimitiveInterface
interface operator (-)
  module procedure negative_prim
  module procedure prim_sub_prim
#ifdef r16p
  module procedure ScalR16P_sub_prim
#endif
  module procedure ScalR8P_sub_prim
  module procedure ScalR4P_sub_prim
  module procedure ScalI8P_sub_prim
  module procedure ScalI4P_sub_prim
  module procedure ScalI2P_sub_prim
  module procedure ScalI1P_sub_prim
#ifdef r16p
  module procedure prim_sub_ScalR16P
#endif
  module procedure prim_sub_ScalR8P
  module procedure prim_sub_ScalR4P
  module procedure prim_sub_ScalI8P
  module procedure prim_sub_ScalI4P
  module procedure prim_sub_ScalI2P
  module procedure prim_sub_ScalI1P
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_PrimitivePublicProcedure
  !> @{
  !> @brief Function for writing Type_Primitive data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function write_primitive(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN), optional:: scalar                !< Scalar primitive data.
  type(Type_Primitive), intent(IN), optional:: array1D(:)            !< One dimensional array primitive data.
  type(Type_Primitive), intent(IN), optional:: array2D(:,:)          !< Two dimensional array primitive data.
  type(Type_Primitive), intent(IN), optional:: array3D(:,:,:)        !< Three dimensional array primitive data.
  character(*),         intent(IN), optional:: format                !< Format specifier.
  integer(I4P),         intent(IN)::           unit                  !< Logic unit.
  integer(I_P)::                               err                   !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                               i1,i2,i3,N(1:2,1:3),j !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(format)) then
    select case(adjustl(trim(format)))
    case('*')
      if (present(scalar)) then
        write(unit,*,iostat=err)scalar%r(:)
        err = write_vector(scalar=scalar%v,format=format,unit=unit)
        write(unit,*,iostat=err)scalar%p
        write(unit,*,iostat=err)scalar%d
        write(unit,*,iostat=err)scalar%g
      elseif (present(array1D)) then
        do j=1,1
          N(1,j) = lbound(array1D,dim=j)
          N(2,j) = ubound(array1D,dim=j)
        enddo
        write(unit,*,iostat=err)(array1D(i1)%r(:),i1=N(1,1),N(2,1))
        err = write_vector(array1D=array1D%v,format=format,unit=unit)
        write(unit,*,iostat=err)array1D%p
        write(unit,*,iostat=err)array1D%d
        write(unit,*,iostat=err)array1D%g
      elseif (present(array2D)) then
        do j=1,2
          N(1,j) = lbound(array2D,dim=j)
          N(2,j) = ubound(array2D,dim=j)
        enddo
        write(unit,*,iostat=err)((array2D(i1,i2)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
        err = write_vector(array2D=array2D%v,format=format,unit=unit)
        write(unit,*,iostat=err)array2D%p
        write(unit,*,iostat=err)array2D%d
        write(unit,*,iostat=err)array2D%g
      elseif (present(array3D)) then
        do j=1,3
          N(1,j) = lbound(array3D,dim=j)
          N(2,j) = ubound(array3D,dim=j)
        enddo
        write(unit,*,iostat=err)(((array3D(i1,i2,i3)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
        err = write_vector(array3D=array3D%v,format=format,unit=unit)
        write(unit,*,iostat=err)array3D%p
        write(unit,*,iostat=err)array3D%d
        write(unit,*,iostat=err)array3D%g
      endif
    case default
      if (present(scalar)) then
        write(unit,adjustl(trim(format)),iostat=err)scalar%r(:)
        err = write_vector(scalar=scalar%v,format=format,unit=unit)
        write(unit,adjustl(trim(format)),iostat=err)scalar%p
        write(unit,adjustl(trim(format)),iostat=err)scalar%d
        write(unit,adjustl(trim(format)),iostat=err)scalar%g
      elseif (present(array1D)) then
        do j=1,1
          N(1,j) = lbound(array1D,dim=j)
          N(2,j) = ubound(array1D,dim=j)
        enddo
        write(unit,adjustl(trim(format)),iostat=err)(array1D(i1)%r(:),i1=N(1,1),N(2,1))
        err = write_vector(array1D=array1D%v,format=format,unit=unit)
        write(unit,adjustl(trim(format)),iostat=err)array1D%p
        write(unit,adjustl(trim(format)),iostat=err)array1D%d
        write(unit,adjustl(trim(format)),iostat=err)array1D%g
      elseif (present(array2D)) then
        do j=1,2
          N(1,j) = lbound(array2D,dim=j)
          N(2,j) = ubound(array2D,dim=j)
        enddo
        write(unit,adjustl(trim(format)),iostat=err)((array2D(i1,i2)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
        err = write_vector(array2D=array2D%v,format=format,unit=unit)
        write(unit,adjustl(trim(format)),iostat=err)array2D%p
        write(unit,adjustl(trim(format)),iostat=err)array2D%d
        write(unit,adjustl(trim(format)),iostat=err)array2D%g
      elseif (present(array3D)) then
        do j=1,3
          N(1,j) = lbound(array3D,dim=j)
          N(2,j) = ubound(array3D,dim=j)
        enddo
        write(unit,adjustl(trim(format)),iostat=err)(((array3D(i1,i2,i3)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
        err = write_vector(array3D=array3D%v,format=format,unit=unit)
        write(unit,adjustl(trim(format)),iostat=err)array3D%p
        write(unit,adjustl(trim(format)),iostat=err)array3D%d
        write(unit,adjustl(trim(format)),iostat=err)array3D%g
      endif
    endselect
  else
    if (present(scalar)) then
      write(unit,iostat=err)scalar%r(:)
      err = write_vector(scalar=scalar%v,unit=unit)
      write(unit,iostat=err)scalar%p
      write(unit,iostat=err)scalar%d
      write(unit,iostat=err)scalar%g
    elseif (present(array1D)) then
      do j=1,1
        N(1,j) = lbound(array1D,dim=j)
        N(2,j) = ubound(array1D,dim=j)
      enddo
      write(unit,iostat=err)(array1D(i1)%r(:),i1=N(1,1),N(2,1))
      err = write_vector(array1D=array1D%v,unit=unit)
      write(unit,iostat=err)array1D%p
      write(unit,iostat=err)array1D%d
      write(unit,iostat=err)array1D%g
    elseif (present(array2D)) then
      do j=1,2
        N(1,j) = lbound(array2D,dim=j)
        N(2,j) = ubound(array2D,dim=j)
      enddo
      write(unit,iostat=err)((array2D(i1,i2)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
      err = write_vector(array2D=array2D%v,unit=unit)
      write(unit,iostat=err)array2D%p
      write(unit,iostat=err)array2D%d
      write(unit,iostat=err)array2D%g
    elseif (present(array3D)) then
      do j=1,3
        N(1,j) = lbound(array3D,dim=j)
        N(2,j) = ubound(array3D,dim=j)
      enddo
      write(unit,iostat=err)(((array3D(i1,i2,i3)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
      err = write_vector(array3D=array3D%v,unit=unit)
      write(unit,iostat=err)array3D%p
      write(unit,iostat=err)array3D%d
      write(unit,iostat=err)array3D%g
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction write_primitive

  !> @brief Function for reading Type_Primitive data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function read_primitive(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT), optional:: scalar                !< Scalar primitive data.
  type(Type_Primitive), intent(INOUT), optional:: array1D(:)            !< One dimensional array primitive data.
  type(Type_Primitive), intent(INOUT), optional:: array2D(:,:)          !< Two dimensional array primitive data.
  type(Type_Primitive), intent(INOUT), optional:: array3D(:,:,:)        !< Three dimensional array primitive data.
  character(*),         intent(IN),    optional:: format                !< Format specifier.
  integer(I4P),         intent(IN)::              unit                  !< Logic unit.
  integer(I_P)::                                  err                   !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                                  i1,i2,i3,N(1:2,1:3),j !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(format)) then
    select case(adjustl(trim(format)))
    case('*')
      if (present(scalar)) then
        read(unit,*,iostat=err)scalar%r(:)
        err = read_vector(scalar=scalar%v,format=format,unit=unit)
        read(unit,*,iostat=err)scalar%p
        read(unit,*,iostat=err)scalar%d
        read(unit,*,iostat=err)scalar%g
      elseif (present(array1D)) then
        do j=1,1
          N(1,j) = lbound(array1D,dim=j)
          N(2,j) = ubound(array1D,dim=j)
        enddo
        read(unit,*,iostat=err)(array1D(i1)%r(:),i1=N(1,1),N(2,1))
        err = read_vector(array1D=array1D%v,format=format,unit=unit)
        read(unit,*,iostat=err)array1D%p
        read(unit,*,iostat=err)array1D%d
        read(unit,*,iostat=err)array1D%g
      elseif (present(array2D)) then
        do j=1,2
          N(1,j) = lbound(array2D,dim=j)
          N(2,j) = ubound(array2D,dim=j)
        enddo
        read(unit,*,iostat=err)((array2D(i1,i2)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
        err = read_vector(array2D=array2D%v,format=format,unit=unit)
        read(unit,*,iostat=err)array2D%p
        read(unit,*,iostat=err)array2D%d
        read(unit,*,iostat=err)array2D%g
      elseif (present(array3D)) then
        do j=1,3
          N(1,j) = lbound(array3D,dim=j)
          N(2,j) = ubound(array3D,dim=j)
        enddo
        read(unit,*,iostat=err)(((array3D(i1,i2,i3)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
        err = read_vector(array3D=array3D%v,format=format,unit=unit)
        read(unit,*,iostat=err)array3D%p
        read(unit,*,iostat=err)array3D%d
        read(unit,*,iostat=err)array3D%g
      endif
    case default
      if (present(scalar)) then
        read(unit,adjustl(trim(format)),iostat=err)scalar%r(:)
        err = read_vector(scalar=scalar%v,format=format,unit=unit)
        read(unit,adjustl(trim(format)),iostat=err)scalar%p
        read(unit,adjustl(trim(format)),iostat=err)scalar%d
        read(unit,adjustl(trim(format)),iostat=err)scalar%g
      elseif (present(array1D)) then
        do j=1,1
          N(1,j) = lbound(array1D,dim=j)
          N(2,j) = ubound(array1D,dim=j)
        enddo
        read(unit,adjustl(trim(format)),iostat=err)(array1D(i1)%r(:),i1=N(1,1),N(2,1))
        err = read_vector(array1D=array1D%v,format=format,unit=unit)
        read(unit,adjustl(trim(format)),iostat=err)array1D%p
        read(unit,adjustl(trim(format)),iostat=err)array1D%d
        read(unit,adjustl(trim(format)),iostat=err)array1D%g
      elseif (present(array2D)) then
        do j=1,2
          N(1,j) = lbound(array2D,dim=j)
          N(2,j) = ubound(array2D,dim=j)
        enddo
        read(unit,adjustl(trim(format)),iostat=err)((array2D(i1,i2)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
        err = read_vector(array2D=array2D%v,format=format,unit=unit)
        read(unit,adjustl(trim(format)),iostat=err)array2D%p
        read(unit,adjustl(trim(format)),iostat=err)array2D%d
        read(unit,adjustl(trim(format)),iostat=err)array2D%g
      elseif (present(array3D)) then
        do j=1,3
          N(1,j) = lbound(array3D,dim=j)
          N(2,j) = ubound(array3D,dim=j)
        enddo
        read(unit,adjustl(trim(format)),iostat=err)(((array3D(i1,i2,i3)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
        err = read_vector(array3D=array3D%v,format=format,unit=unit)
        read(unit,adjustl(trim(format)),iostat=err)array3D%p
        read(unit,adjustl(trim(format)),iostat=err)array3D%d
        read(unit,adjustl(trim(format)),iostat=err)array3D%g
      endif
    endselect
  else
    if (present(scalar)) then
      read(unit,iostat=err)scalar%r(:)
      err = read_vector(scalar=scalar%v,unit=unit)
      read(unit,iostat=err)scalar%p
      read(unit,iostat=err)scalar%d
      read(unit,iostat=err)scalar%g
    elseif (present(array1D)) then
      do j=1,1
        N(1,j) = lbound(array1D,dim=j)
        N(2,j) = ubound(array1D,dim=j)
      enddo
      read(unit,iostat=err)(array1D(i1)%r(:),i1=N(1,1),N(2,1))
      err = read_vector(array1D=array1D%v,unit=unit)
      read(unit,iostat=err)array1D%p
      read(unit,iostat=err)array1D%d
      read(unit,iostat=err)array1D%g
    elseif (present(array2D)) then
      do j=1,2
        N(1,j) = lbound(array2D,dim=j)
        N(2,j) = ubound(array2D,dim=j)
      enddo
      read(unit,iostat=err)((array2D(i1,i2)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
      err = read_vector(array2D=array2D%v,unit=unit)
      read(unit,iostat=err)array2D%p
      read(unit,iostat=err)array2D%d
      read(unit,iostat=err)array2D%g
    elseif (present(array3D)) then
      do j=1,3
        N(1,j) = lbound(array3D,dim=j)
        N(2,j) = ubound(array3D,dim=j)
      enddo
      read(unit,iostat=err)(((array3D(i1,i2,i3)%r(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
      err = read_vector(array3D=array3D%v,unit=unit)
      read(unit,iostat=err)array3D%p
      read(unit,iostat=err)array3D%d
      read(unit,iostat=err)array3D%g
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction read_primitive
  !> @}

  !> @ingroup Data_Type_PrimitivePrivateProcedure
  !> @{
  !> @brief Subroutine for initializing Type_Primitive allocatable variables.
  elemental subroutine init(prim,Ns,prim0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive), intent(INOUT)::        prim  !< Primitive initialized data.
  integer(I_P),          intent(IN), optional:: Ns    !< Number of species.
  type(Type_Primitive),  intent(IN), optional:: prim0 !< Primitive initialization data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(Ns)) then
    if (allocated(prim%r)) deallocate(prim%r) ; allocate(prim%r(1:Ns)) ; prim%r = 0._R8P
  elseif (present(prim0)) then
    if (allocated(prim%r)) deallocate(prim%r) ; allocate(prim%r(1:size(prim0%r))) ; prim = prim0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  !> @brief Subroutine for freeing the memory of Type_Primitive allocatable variables.
  elemental subroutine free(prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive), intent(INOUT):: prim !< Primitive data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim%r)) deallocate(prim%r)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Function for converting derived type Type_Primitive to array.
  !> @return \b array real(R8P), dimension(1:size(prim\%r)+6) variable.
  pure function prim2array(prim) result(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive), intent(IN):: prim                    !< Derived type primitive data.
  real(R8P)::                         array(1:size(prim%r)+6) !< Primitive data in the form of array.
  integer(I_P)::                      Ns                      !< Number of species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(prim%r)
  array(1:Ns) = prim%r
  array(Ns+1) = prim%v%x
  array(Ns+2) = prim%v%y
  array(Ns+3) = prim%v%z
  array(Ns+4) = prim%p
  array(Ns+5) = prim%d
  array(Ns+6) = prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim2array

  !> @brief Function for converting array to derived type Type_Primitive.
  !> @return \b prim type(Type_Primitive) variable.
  pure subroutine array2prim(prim,array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive), intent(INOUT):: prim     !< Derived type primitive data.
  real(R8P),             intent(IN)::    array(:) !< Primitive data in the form of array.
  integer(I_P)::                         Ns       !< Number of species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(array)-6
  if (allocated(prim%r)) deallocate(prim%r)  ; allocate(prim%r(1:Ns))
  prim%r   = array(1:Ns)
  prim%v%x = array(Ns+1)
  prim%v%y = array(Ns+2)
  prim%v%z = array(Ns+3)
  prim%p   = array(Ns+4)
  prim%d   = array(Ns+5)
  prim%g   = array(Ns+6)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array2prim

  !> @brief Function for printing in a pretty ascii format the components of type Type_Primitive.
  !> @return \b err integer(I_P) variable for error trapping.
  function pprint(prim,myrank,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Primitive), intent(IN)::           prim   !< Primitives.
  integer(I_P),          intent(IN), optional:: myrank !< Actual rank process.
  integer(I4P),          intent(IN)::           unit   !< Logic unit.
  integer(I_P)::                                err    !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                                Ns     !< Number of species.
  integer(I_P)::                                s      !< Species counter.
  character(DI_P)::                             rks    !< String containing myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,0_I_P)) ; if (present(myrank)) rks = 'rank'//trim(str(.true.,myrank))
  Ns  = size(prim%r)
  do s=1,Ns
    write(unit,'(A)',iostat=err)trim(rks)//' r('//trim(str(.true.,s))//')='//str(n=prim%r(s))
  enddo
    write(unit,'(A)',iostat=err)trim(rks)//' v(x)='//str(n=prim%v%x)
    write(unit,'(A)',iostat=err)trim(rks)//' v(y)='//str(n=prim%v%y)
    write(unit,'(A)',iostat=err)trim(rks)//' v(z)='//str(n=prim%v%z)
    write(unit,'(A)',iostat=err)trim(rks)//' p   ='//str(n=prim%p  )
    write(unit,'(A)',iostat=err)trim(rks)//' d   ='//str(n=prim%d  )
    write(unit,'(A)',iostat=err)trim(rks)//' g   ='//str(n=prim%g  )
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction pprint

  ! Assignment (=)
  elemental subroutine assign_prim(prim1,prim2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between two prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim1
  type(Type_Primitive), intent(IN)::    prim2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim1%r = prim2%r
  prim1%v = prim2%v
  prim1%p = prim2%p
  prim1%d = prim2%d
  prim1%g = prim2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_prim

#ifdef r16p
  elemental subroutine assign_ScalR16P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R16P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  real(R16P),           intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%r = real(scal,R8P)
  prim%v = real(scal,R8P)
  prim%p = real(scal,R8P)
  prim%d = real(scal,R8P)
  prim%g = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR16P
#endif

  elemental subroutine assign_ScalR8P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R8P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  real(R8P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%r = real(scal,R8P)
  prim%v = real(scal,R8P)
  prim%p = real(scal,R8P)
  prim%d = real(scal,R8P)
  prim%g = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR8P

  elemental subroutine assign_ScalR4P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R4P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  real(R4P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%r = real(scal,R8P)
  prim%v = real(scal,R8P)
  prim%p = real(scal,R8P)
  prim%d = real(scal,R8P)
  prim%g = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR4P

  elemental subroutine assign_ScalI8P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I8P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  integer(I8P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%r = real(scal,R8P)
  prim%v = real(scal,R8P)
  prim%p = real(scal,R8P)
  prim%d = real(scal,R8P)
  prim%g = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI8P

  elemental subroutine assign_ScalI4P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I4P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  integer(I4P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%r = real(scal,R8P)
  prim%v = real(scal,R8P)
  prim%p = real(scal,R8P)
  prim%d = real(scal,R8P)
  prim%g = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI4P

  elemental subroutine assign_ScalI2P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I2P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  integer(I2P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%r = real(scal,R8P)
  prim%v = real(scal,R8P)
  prim%p = real(scal,R8P)
  prim%d = real(scal,R8P)
  prim%g = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI2P

  elemental subroutine assign_ScalI1P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I1P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  integer(I1P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%r = real(scal,R8P)
  prim%v = real(scal,R8P)
  prim%p = real(scal,R8P)
  prim%d = real(scal,R8P)
  prim%g = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI1P

  ! Multiplication (*)
  elemental function prim_mul_prim(prim1,prim2) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply (by components) Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN)::  prim1 ! First prim obj.
  type(Type_Primitive), intent(IN)::  prim2 ! Second prim obj.
  type(Type_Primitive)::              mul   ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim1%r)))
  mul%r = prim1%r * prim2%r
  mul%v = prim1%v * prim2%v
  mul%p = prim1%p * prim2%p
  mul%d = prim1%d * prim2%d
  mul%g = prim1%g * prim2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_prim

#ifdef r16p
  elemental function ScalR16P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R16P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),           intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_mul_prim

  elemental function prim_mul_ScalR16P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R16P),           intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalR16P
#endif

  elemental function ScalR8P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_mul_prim

  elemental function prim_mul_ScalR8P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R8P),            intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalR8P

  elemental function ScalR4P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_mul_prim

  elemental function prim_mul_ScalR4P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R4P),            intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalR4P

  elemental function ScalI8P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_mul_prim

  elemental function prim_mul_ScalI8P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalI8P

  elemental function ScalI4P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_mul_prim

  elemental function prim_mul_ScalI4P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalI4P

  elemental function ScalI2P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I2P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_mul_prim

  elemental function prim_mul_ScalI2P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalI2P

  elemental function ScalI1P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I1P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_mul_prim

  elemental function prim_mul_ScalI1P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R8P) * prim%r
  mul%v = real(scal,R8P) * prim%v
  mul%p = real(scal,R8P) * prim%p
  mul%d = real(scal,R8P) * prim%d
  mul%g = real(scal,R8P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalI1P

  ! Division (/)
  elemental function prim_div_prim(prim1,prim2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide (by components) Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN)::  prim1 ! First prim obj.
  type(Type_Primitive), intent(IN)::  prim2 ! Second prim obj.
  type(Type_Primitive)::              div   ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%r(1:size(prim1%r)))
  div%r = prim1%r / prim2%r
  div%v = prim1%v / prim2%v
  div%p = prim1%p / prim2%p
  div%d = prim1%d / prim2%d
  div%g = prim1%g / prim2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_prim

#ifdef r16p
  elemental function prim_div_ScalR16P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R16P),           intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%r(1:size(prim1%r)))
  div%r = prim%r / real(scal,R8P)
  div%v = prim%v / real(scal,R8P)
  div%p = prim%p / real(scal,R8P)
  div%d = prim%d / real(scal,R8P)
  div%g = prim%g / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalR16P
#endif

  elemental function prim_div_ScalR8P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R8P),            intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R8P)
  div%v = prim%v / real(scal,R8P)
  div%p = prim%p / real(scal,R8P)
  div%d = prim%d / real(scal,R8P)
  div%g = prim%g / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalR8P

  elemental function prim_div_ScalR4P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R4P),            intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R8P)
  div%v = prim%v / real(scal,R8P)
  div%p = prim%p / real(scal,R8P)
  div%d = prim%d / real(scal,R8P)
  div%g = prim%g / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalR4P

  elemental function prim_div_ScalI8P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R8P)
  div%v = prim%v / real(scal,R8P)
  div%p = prim%p / real(scal,R8P)
  div%d = prim%d / real(scal,R8P)
  div%g = prim%g / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalI8P

  elemental function prim_div_ScalI4P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R8P)
  div%v = prim%v / real(scal,R8P)
  div%p = prim%p / real(scal,R8P)
  div%d = prim%d / real(scal,R8P)
  div%g = prim%g / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalI4P

  elemental function prim_div_ScalI2P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R8P)
  div%v = prim%v / real(scal,R8P)
  div%p = prim%p / real(scal,R8P)
  div%d = prim%d / real(scal,R8P)
  div%g = prim%g / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalI2P

  elemental function prim_div_ScalI1P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R8P)
  div%v = prim%v / real(scal,R8P)
  div%p = prim%p / real(scal,R8P)
  div%d = prim%d / real(scal,R8P)
  div%g = prim%g / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalI1P

  ! Sum (+)
  elemental function positive_prim(prim) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary + to a Primitive objecy.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(pos%r(1:size(prim%r)))
  pos%r =  + prim%r
  pos%v =  + prim%v
  pos%p =  + prim%p
  pos%d =  + prim%d
  pos%g =  + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_prim

  elemental function prim_sum_prim(prim1,prim2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum (by components) Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN)::  prim1 ! First prim obj.
  type(Type_Primitive), intent(IN)::  prim2 ! Second prim obj.
  type(Type_Primitive)::              summ  ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim1%r)))
  summ%r = prim1%r + prim2%r
  summ%v = prim1%v + prim2%v
  summ%p = prim1%p + prim2%p
  summ%d = prim1%d + prim2%d
  summ%g = prim1%g + prim2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_prim

#ifdef r16p
  elemental function ScalR16P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R16P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),           intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sum_prim

  elemental function prim_sum_ScalR16P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R16P),           intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalR16P
#endif

  elemental function ScalR8P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sum_prim

  elemental function prim_sum_ScalR8P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R8P),            intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalR8P

  elemental function ScalR4P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sum_prim

  elemental function prim_sum_ScalR4P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R4P),            intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalR4P

  elemental function ScalI8P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sum_prim

  elemental function prim_sum_ScalI8P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalI8P

  elemental function ScalI4P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sum_prim

  elemental function prim_sum_ScalI4P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalI4P

  elemental function ScalI2P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I2P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sum_prim

  elemental function prim_sum_ScalI2P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalI2P

  elemental function ScalI1P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I1P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sum_prim

  elemental function prim_sum_ScalI1P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R8P) + prim%r
  summ%v = real(scal,R8P) + prim%v
  summ%p = real(scal,R8P) + prim%p
  summ%d = real(scal,R8P) + prim%d
  summ%g = real(scal,R8P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalI1P

  ! Subtraction (-)
  elemental function negative_prim(prim) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary - to a Primitive objecy.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(neg%r(1:size(prim%r)))
  neg%r =  - prim%r
  neg%v =  - prim%v
  neg%p =  - prim%p
  neg%d =  - prim%d
  neg%g =  - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_prim

  elemental function prim_sub_prim(prim1,prim2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract (by components) Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN)::  prim1 ! First prim obj.
  type(Type_Primitive), intent(IN)::  prim2 ! Second prim obj.
  type(Type_Primitive)::              sub  ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim1%r)))
  sub%r = prim1%r - prim2%r
  sub%v = prim1%v - prim2%v
  sub%p = prim1%p - prim2%p
  sub%d = prim1%d - prim2%d
  sub%g = prim1%g - prim2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_prim

#ifdef r16p
  elemental function ScalR16P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R16P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),           intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R8P) - prim%r
  sub%v = real(scal,R8P) - prim%v
  sub%p = real(scal,R8P) - prim%p
  sub%d = real(scal,R8P) - prim%d
  sub%g = real(scal,R8P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sub_prim

  elemental function prim_sub_ScalR16P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R16P),           intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R8P)
  sub%v = prim%v - real(scal,R8P)
  sub%p = prim%p - real(scal,R8P)
  sub%d = prim%d - real(scal,R8P)
  sub%g = prim%g - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalR16P
#endif

  elemental function ScalR8P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R8P) - prim%r
  sub%v = real(scal,R8P) - prim%v
  sub%p = real(scal,R8P) - prim%p
  sub%d = real(scal,R8P) - prim%d
  sub%g = real(scal,R8P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sub_prim

  elemental function prim_sub_ScalR8P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R8P),            intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R8P)
  sub%v = prim%v - real(scal,R8P)
  sub%p = prim%p - real(scal,R8P)
  sub%d = prim%d - real(scal,R8P)
  sub%g = prim%g - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalR8P

  elemental function ScalR4P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R8P) - prim%r
  sub%v = real(scal,R8P) - prim%v
  sub%p = real(scal,R8P) - prim%p
  sub%d = real(scal,R8P) - prim%d
  sub%g = real(scal,R8P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sub_prim

  elemental function prim_sub_ScalR4P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R4P),            intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R8P)
  sub%v = prim%v - real(scal,R8P)
  sub%p = prim%p - real(scal,R8P)
  sub%d = prim%d - real(scal,R8P)
  sub%g = prim%g - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalR4P

  elemental function ScalI8P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R8P) - prim%r
  sub%v = real(scal,R8P) - prim%v
  sub%p = real(scal,R8P) - prim%p
  sub%d = real(scal,R8P) - prim%d
  sub%g = real(scal,R8P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sub_prim

  elemental function prim_sub_ScalI8P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R8P)
  sub%v = prim%v - real(scal,R8P)
  sub%p = prim%p - real(scal,R8P)
  sub%d = prim%d - real(scal,R8P)
  sub%g = prim%g - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalI8P

  elemental function ScalI4P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R8P) - prim%r
  sub%v = real(scal,R8P) - prim%v
  sub%p = real(scal,R8P) - prim%p
  sub%d = real(scal,R8P) - prim%d
  sub%g = real(scal,R8P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sub_prim

  elemental function prim_sub_ScalI4P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R8P)
  sub%v = prim%v - real(scal,R8P)
  sub%p = prim%p - real(scal,R8P)
  sub%d = prim%d - real(scal,R8P)
  sub%g = prim%g - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalI4P

  elemental function ScalI2P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I2P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R8P) - prim%r
  sub%v = real(scal,R8P) - prim%v
  sub%p = real(scal,R8P) - prim%p
  sub%d = real(scal,R8P) - prim%d
  sub%g = real(scal,R8P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sub_prim

  elemental function prim_sub_ScalI2P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R8P)
  sub%v = prim%v - real(scal,R8P)
  sub%p = prim%p - real(scal,R8P)
  sub%d = prim%d - real(scal,R8P)
  sub%g = prim%g - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalI2P

  elemental function ScalI1P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I1P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R8P) - prim%r
  sub%v = real(scal,R8P) - prim%v
  sub%p = real(scal,R8P) - prim%p
  sub%d = real(scal,R8P) - prim%d
  sub%g = real(scal,R8P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sub_prim

  elemental function prim_sub_ScalI1P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R8P)
  sub%v = prim%v - real(scal,R8P)
  sub%p = prim%p - real(scal,R8P)
  sub%d = prim%d - real(scal,R8P)
  sub%g = prim%g - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalI1P
  !> @}
endmodule Data_Type_Primitive
