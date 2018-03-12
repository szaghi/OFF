!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_ConservativeDerivedType Data_Type_Conservative
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_ConservativeInterface Data_Type_Conservative
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_ConservativePublicProcedure Data_Type_Conservative
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_ConservativePrivateProcedure Data_Type_Conservative
!> @}

!> This module contains the definition of Type_Conservative and its procedures.
!> Type_Conservative is a derived type that handles conservative fluid dynamic variables.
!> @note The operators of assignment (=), multiplication (*), division (/), sum (+) and subtraction (-) have been overloaded.
!> Therefore this module provides a far-complete algebra based on Type_Conservative derived type. This algebra simplifies the
!> numerical integration of Partial Differential Equations (PDE) systems based on conservative formulation.
!> @todo \b DocComplete: Complete the documentation of internal procedures
module Data_Type_Conservative
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision     ! Integers and reals precision definition.
USE Data_Type_Vector ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: write_conservative,read_conservative
public:: array2consf
public:: assignment (=)
public:: operator (*)
public:: operator (/)
public:: operator (+)
public:: operator (-)
public:: operator (.dot.)
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
    procedure:: init              ! Procedure for initializing allocatable variables.
    procedure:: free => free_cons ! Procedure for freeing the memory of allocatable variables.
    procedure:: cons2array        ! Procedure for converting derived type Type_Conservative to array.
    procedure:: array2cons        ! Procedure for converting array to derived type Type_Conservative.
    procedure:: pprint            ! Procedure for printing Type_Conservative components with a "pretty" format.
endtype Type_Conservative
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Assignment operator (=) overloading.
!> @ingroup Data_Type_ConservativeInterface
interface assignment (=)
  module procedure assign_cons
#ifdef r16p
  module procedure assign_ScalR16P
#endif
  module procedure assign_ScalR8P
  module procedure assign_ScalR4P
  module procedure assign_ScalI8P
  module procedure assign_ScalI4P
  module procedure assign_ScalI2P
  module procedure assign_ScalI1P
end interface
!> @brief Multiplication operator (*) overloading.
!> @note The admissible multiplications are:
!>       - Type_Conservative * Type_Conservative: each component of first conservative variable (cons1) is multiplied for the
!>         corresponding component of the second one (cons2), i.e. \n
!>         \f$ {\rm result\%rs = cons1\%rs*cons2\%rs} \f$ \n
!>         \f$ {\rm result\%rv = cons1\%rv*cons2\%rv} \f$ \n
!>         \f$ {\rm result\%re = cons1\%re*cons2\%re} \f$ \n
!>       - scalar number (real or integer of any kinds defined in IR_Precision module) * Type_Conservative: each component of
!>         Type_Conservative is multiplied for the scalar, i.e. \n
!>         \f$ {\rm result\%rs = cons\%rs*scalar} \f$ \n
!>         \f$ {\rm result\%rv = cons\%rv*scalar} \f$ \n
!>         \f$ {\rm result\%re = cons\%re*scalar} \f$ \n
!>       - Type_Conservative * scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Conservative is multiplied for the scalar, i.e. \n
!>         \f$ {\rm result\%rs = cons\%rs*scalar} \f$ \n
!>         \f$ {\rm result\%rv = cons\%rv*scalar} \f$ \n
!>         \f$ {\rm result\%re = cons\%re*scalar} \f$ \n
!> @ingroup Data_Type_ConservativeInterface
interface operator (*)
  module procedure cons_mul_cons
#ifdef r16p
  module procedure ScalR16P_mul_cons
#endif
  module procedure ScalR8P_mul_cons
  module procedure ScalR4P_mul_cons
  module procedure ScalI8P_mul_cons
  module procedure ScalI4P_mul_cons
  module procedure ScalI2P_mul_cons
  module procedure ScalI1P_mul_cons
#ifdef r16p
  module procedure cons_mul_ScalR16P
#endif
  module procedure cons_mul_ScalR8P
  module procedure cons_mul_ScalR4P
  module procedure cons_mul_ScalI8P
  module procedure cons_mul_ScalI4P
  module procedure cons_mul_ScalI2P
  module procedure cons_mul_ScalI1P
endinterface
!> @brief Division operator (/) overloading.
!> @note The admissible divisions are:
!>       - Type_Conservative / Type_Conservative: each component of first conservative variable (cons1) is divided for the
!>         corresponding component of the second one (cons2): i.e. \n
!>         \f$ {\rm result\%rs = \frac{cons1\%rs}{cons2\%rs}} \f$ \n
!>         \f$ {\rm result\%rv = \frac{cons1\%rv}{cons2\%rv}} \f$ \n
!>         \f$ {\rm result\%re = \frac{cons1\%re}{cons2\%re}} \f$ \n
!>       - Type_Conservative / scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Conservative is divided for the scalar, i.e. \n
!>         \f$ {\rm result\%rs = \frac{cons\%rs}{scalar}} \f$ \n
!>         \f$ {\rm result\%rv = \frac{cons\%rv}{scalar}} \f$ \n
!>         \f$ {\rm result\%re = \frac{cons\%re}{scalar}} \f$ \n
!> @ingroup Data_Type_ConservativeInterface
interface operator (/)
  module procedure cons_div_cons
#ifdef r16p
  module procedure cons_div_ScalR16P
#endif
  module procedure cons_div_ScalR8P
  module procedure cons_div_ScalR4P
  module procedure cons_div_ScalI8P
  module procedure cons_div_ScalI4P
  module procedure cons_div_ScalI2P
  module procedure cons_div_ScalI1P
endinterface
!> @brief Sum operator (+) overloading.
!> @note The admissible summations are:
!>       - Type_Conservative + Type_Conservative: each component of first conservative variable (cons1) is summed with the
!>         corresponding component of the second one (cons2): i.e. \n
!>         \f$ {\rm result\%rs = cons1\%rs+cons2\%rs} \f$ \n
!>         \f$ {\rm result\%rv = cons1\%rv+cons2\%rv} \f$ \n
!>         \f$ {\rm result\%re = cons1\%re+cons2\%re} \f$ \n
!>       - scalar number (real or integer of any kinds defined in IR_Precision module) + Type_Conservative: each component of
!>         Type_Conservative is summed for the scalar, i.e. \n
!>         \f$ {\rm result\%rs = cons\%rs+scalar} \f$ \n
!>         \f$ {\rm result\%rv = cons\%rv+scalar} \f$ \n
!>         \f$ {\rm result\%re = cons\%re+scalar} \f$ \n
!>       - Type_Conservative + scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Conservative is summed for the scalar, i.e. \n
!>         \f$ {\rm result\%rs = cons\%rs+scalar} \f$ \n
!>         \f$ {\rm result\%rv = cons\%rv+scalar} \f$ \n
!>         \f$ {\rm result\%re = cons\%re+scalar} \f$ \n
!> @ingroup Data_Type_ConservativeInterface
interface operator (+)
  module procedure positive_cons
  module procedure cons_sum_cons
#ifdef r16p
  module procedure ScalR16P_sum_cons
#endif
  module procedure ScalR8P_sum_cons
  module procedure ScalR4P_sum_cons
  module procedure ScalI8P_sum_cons
  module procedure ScalI4P_sum_cons
  module procedure ScalI2P_sum_cons
  module procedure ScalI1P_sum_cons
#ifdef r16p
  module procedure cons_sum_ScalR16P
#endif
  module procedure cons_sum_ScalR8P
  module procedure cons_sum_ScalR4P
  module procedure cons_sum_ScalI8P
  module procedure cons_sum_ScalI4P
  module procedure cons_sum_ScalI2P
  module procedure cons_sum_ScalI1P
endinterface
!> @brief Subtraction operator (-) overloading.
!> @note The admissible subtractions are:
!>       - Type_Conservative - Type_Conservative: each component of first conservative variable (cons1) is subtracted with the
!>         corresponding component of the second one (cons2): i.e. \n
!>         \f$ {\rm result\%rs = cons1\%rs-cons2\%rs} \f$ \n
!>         \f$ {\rm result\%rv = cons1\%rv-cons2\%rv} \f$ \n
!>         \f$ {\rm result\%re = cons1\%re-cons2\%re} \f$ \n
!>       - scalar number (real or integer of any kinds defined in IR_Precision module) - Type_Conservative: each component of
!>         Type_Conservative is subtracted with the scalar, i.e. \n
!>         \f$ {\rm result\%rs = scalar-cons\%rs} \f$ \n
!>         \f$ {\rm result\%rv = scalar-cons\%rv} \f$ \n
!>         \f$ {\rm result\%re = scalar-cons\%re} \f$ \n
!>       - Type_Conservative - scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Conservative is subtracted with the scalar, i.e. \n
!>         \f$ {\rm result\%rs = cons\%rs-scalar} \f$ \n
!>         \f$ {\rm result\%rv = cons\%rv-scalar} \f$ \n
!>         \f$ {\rm result\%re = cons\%re-scalar} \f$ \n
!> @ingroup Data_Type_ConservativeInterface
interface operator (-)
  module procedure negative_cons
  module procedure cons_sub_cons
#ifdef r16p
  module procedure ScalR16P_sub_cons
#endif
  module procedure ScalR8P_sub_cons
  module procedure ScalR4P_sub_cons
  module procedure ScalI8P_sub_cons
  module procedure ScalI4P_sub_cons
  module procedure ScalI2P_sub_cons
  module procedure ScalI1P_sub_cons
#ifdef r16p
  module procedure cons_sub_ScalR16P
#endif
  module procedure cons_sub_ScalR8P
  module procedure cons_sub_ScalR4P
  module procedure cons_sub_ScalI8P
  module procedure cons_sub_ScalI4P
  module procedure cons_sub_ScalI2P
  module procedure cons_sub_ScalI1P
endinterface
!> @brief Dot product operator (.dot.) definition.
!> This operator is defined as:
!> \f$ {\rm D}= x_1 \cdot x_2 + y_1 \cdot y_2 + z_1 \cdot z_2 \f$
!> @ingroup Data_Type_ConservativeInterface
interface operator (.dot.)
  module procedure cons_dot_cons
  module procedure scalR8_dot_cons,scalR4_dot_cons,scalI8_dot_cons,scalI4_dot_cons,scalI2_dot_cons,scalI1_dot_cons
  module procedure cons_dot_scalR8,cons_dot_scalR4,cons_dot_scalI8,cons_dot_scalI4,cons_dot_scalI2,cons_dot_scalI1
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_ConservativePublicProcedure
  !> @{
  !> @brief Function for writing Type_Conservative data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function write_conservative(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN), optional:: scalar                !< Scalar conservative data.
  type(Type_Conservative), intent(IN), optional:: array1D(:)            !< One dimensional array conservative data.
  type(Type_Conservative), intent(IN), optional:: array2D(:,:)          !< Two dimensional array conservative data.
  type(Type_Conservative), intent(IN), optional:: array3D(:,:,:)        !< Three dimensional array conservative data.
  character(*),            intent(IN), optional:: format                !< Format specifier.
  integer(I4P),            intent(IN)::           unit                  !< Logic unit.
  integer(I_P)::                                  err                   !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                                  i1,i2,i3,N(1:2,1:3),j !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(format)) then
    select case(adjustl(trim(format)))
    case('*')
      if (present(scalar)) then
        write(unit,*,iostat=err)scalar%rs(:)
        err = write_vector(scalar=scalar%rv,format=format,unit=unit)
        write(unit,*,iostat=err)scalar%re
      elseif (present(array1D)) then
        do j=1,1
          N(1,j) = lbound(array1D,dim=j)
          N(2,j) = ubound(array1D,dim=j)
        enddo
        write(unit,*,iostat=err)(array1D(i1)%rs(:),i1=N(1,1),N(2,1))
        err = write_vector(array1D=array1D%rv,format=format,unit=unit)
        write(unit,*,iostat=err)array1D%re
      elseif (present(array2D)) then
        do j=1,2
          N(1,j) = lbound(array2D,dim=j)
          N(2,j) = ubound(array2D,dim=j)
        enddo
        write(unit,*,iostat=err)((array2D(i1,i2)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
        err = write_vector(array2D=array2D%rv,format=format,unit=unit)
        write(unit,*,iostat=err)array2D%re
      elseif (present(array3D)) then
        do j=1,3
          N(1,j) = lbound(array3D,dim=j)
          N(2,j) = ubound(array3D,dim=j)
        enddo
        write(unit,*,iostat=err)(((array3D(i1,i2,i3)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
        err = write_vector(array3D=array3D%rv,format=format,unit=unit)
        write(unit,*,iostat=err)array3D%re
      endif
    case default
      if (present(scalar)) then
        write(unit,adjustl(trim(format)),iostat=err)scalar%rs(:)
        err = write_vector(scalar=scalar%rv,format=format,unit=unit)
        write(unit,adjustl(trim(format)),iostat=err)scalar%re
      elseif (present(array1D)) then
        do j=1,1
          N(1,j) = lbound(array1D,dim=j)
          N(2,j) = ubound(array1D,dim=j)
        enddo
        write(unit,adjustl(trim(format)),iostat=err)(array1D(i1)%rs(:),i1=N(1,1),N(2,1))
        err = write_vector(array1D=array1D%rv,format=format,unit=unit)
        write(unit,adjustl(trim(format)),iostat=err)array1D%re
      elseif (present(array2D)) then
        do j=1,2
          N(1,j) = lbound(array2D,dim=j)
          N(2,j) = ubound(array2D,dim=j)
        enddo
        write(unit,adjustl(trim(format)),iostat=err)((array2D(i1,i2)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
        err = write_vector(array2D=array2D%rv,format=format,unit=unit)
        write(unit,adjustl(trim(format)),iostat=err)array2D%re
      elseif (present(array3D)) then
        do j=1,3
          N(1,j) = lbound(array3D,dim=j)
          N(2,j) = ubound(array3D,dim=j)
        enddo
        write(unit,adjustl(trim(format)),iostat=err)(((array3D(i1,i2,i3)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
        err = write_vector(array3D=array3D%rv,format=format,unit=unit)
        write(unit,adjustl(trim(format)),iostat=err)array3D%re
      endif
    endselect
  else
    if (present(scalar)) then
      write(unit,iostat=err)scalar%rs(:)
      err = write_vector(scalar=scalar%rv,unit=unit)
      write(unit,iostat=err)scalar%re
    elseif (present(array1D)) then
      do j=1,1
        N(1,j) = lbound(array1D,dim=j)
        N(2,j) = ubound(array1D,dim=j)
      enddo
      write(unit,iostat=err)(array1D(i1)%rs(:),i1=N(1,1),N(2,1))
      err = write_vector(array1D=array1D%rv,unit=unit)
      write(unit,iostat=err)array1D%re
    elseif (present(array2D)) then
      do j=1,2
        N(1,j) = lbound(array2D,dim=j)
        N(2,j) = ubound(array2D,dim=j)
      enddo
      write(unit,iostat=err)((array2D(i1,i2)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
      err = write_vector(array2D=array2D%rv,unit=unit)
      write(unit,iostat=err)array2D%re
    elseif (present(array3D)) then
      do j=1,3
        N(1,j) = lbound(array3D,dim=j)
        N(2,j) = ubound(array3D,dim=j)
      enddo
      write(unit,iostat=err)(((array3D(i1,i2,i3)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
      err = write_vector(array3D=array3D%rv,unit=unit)
      write(unit,iostat=err)array3D%re
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction write_conservative

  !> @brief Function for reading Type_Conservative data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function read_conservative(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT), optional:: scalar                !< Scalar conservative data.
  type(Type_Conservative), intent(INOUT), optional:: array1D(:)            !< One dimensional array conservative data.
  type(Type_Conservative), intent(INOUT), optional:: array2D(:,:)          !< Two dimensional array conservative data.
  type(Type_Conservative), intent(INOUT), optional:: array3D(:,:,:)        !< Three dimensional array conservative data.
  character(*),            intent(IN),    optional:: format                !< Format specifier.
  integer(I4P),            intent(IN)::              unit                  !< Logic unit.
  integer(I_P)::                                     err                   !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                                     i1,i2,i3,N(1:2,1:3),j !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(format)) then
    select case(adjustl(trim(format)))
    case('*')
      if (present(scalar)) then
        read(unit,*,iostat=err)scalar%rs(:)
        err = read_vector(scalar=scalar%rv,format=format,unit=unit)
        read(unit,*,iostat=err)scalar%re
      elseif (present(array1D)) then
        do j=1,1
          N(1,j) = lbound(array1D,dim=j)
          N(2,j) = ubound(array1D,dim=j)
        enddo
        read(unit,*,iostat=err)(array1D(i1)%rs(:),i1=N(1,1),N(2,1))
        err = read_vector(array1D=array1D%rv,format=format,unit=unit)
        read(unit,*,iostat=err)array1D%re
      elseif (present(array2D)) then
        do j=1,2
          N(1,j) = lbound(array2D,dim=j)
          N(2,j) = ubound(array2D,dim=j)
        enddo
        read(unit,*,iostat=err)((array2D(i1,i2)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
        err = read_vector(array2D=array2D%rv,format=format,unit=unit)
        read(unit,*,iostat=err)array2D%re
      elseif (present(array3D)) then
        do j=1,3
          N(1,j) = lbound(array3D,dim=j)
          N(2,j) = ubound(array3D,dim=j)
        enddo
        read(unit,*,iostat=err)(((array3D(i1,i2,i3)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
        err = read_vector(array3D=array3D%rv,format=format,unit=unit)
        read(unit,*,iostat=err)array3D%re
      endif
    case default
      if (present(scalar)) then
        read(unit,adjustl(trim(format)),iostat=err)scalar%rs(:)
        err = read_vector(scalar=scalar%rv,format=format,unit=unit)
        read(unit,adjustl(trim(format)),iostat=err)scalar%re
      elseif (present(array1D)) then
        do j=1,1
          N(1,j) = lbound(array1D,dim=j)
          N(2,j) = ubound(array1D,dim=j)
        enddo
        read(unit,adjustl(trim(format)),iostat=err)(array1D(i1)%rs(:),i1=N(1,1),N(2,1))
        err = read_vector(array1D=array1D%rv,format=format,unit=unit)
        read(unit,adjustl(trim(format)),iostat=err)array1D%re
      elseif (present(array2D)) then
        do j=1,2
          N(1,j) = lbound(array2D,dim=j)
          N(2,j) = ubound(array2D,dim=j)
        enddo
        read(unit,adjustl(trim(format)),iostat=err)((array2D(i1,i2)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
        err = read_vector(array2D=array2D%rv,format=format,unit=unit)
        read(unit,adjustl(trim(format)),iostat=err)array2D%re
      elseif (present(array3D)) then
        do j=1,3
          N(1,j) = lbound(array3D,dim=j)
          N(2,j) = ubound(array3D,dim=j)
        enddo
        read(unit,adjustl(trim(format)),iostat=err)(((array3D(i1,i2,i3)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
        err = read_vector(array3D=array3D%rv,format=format,unit=unit)
        read(unit,adjustl(trim(format)),iostat=err)array3D%re
      endif
    endselect
  else
    if (present(scalar)) then
      read(unit,iostat=err)scalar%rs(:)
      err = read_vector(scalar=scalar%rv,unit=unit)
      read(unit,iostat=err)scalar%re
    elseif (present(array1D)) then
      do j=1,1
        N(1,j) = lbound(array1D,dim=j)
        N(2,j) = ubound(array1D,dim=j)
      enddo
      read(unit,iostat=err)(array1D(i1)%rs(:),i1=N(1,1),N(2,1))
      err = read_vector(array1D=array1D%rv,unit=unit)
      read(unit,iostat=err)array1D%re
    elseif (present(array2D)) then
      do j=1,2
        N(1,j) = lbound(array2D,dim=j)
        N(2,j) = ubound(array2D,dim=j)
      enddo
      read(unit,iostat=err)((array2D(i1,i2)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2))
      err = read_vector(array2D=array2D%rv,unit=unit)
      read(unit,iostat=err)array2D%re
    elseif (present(array3D)) then
      do j=1,3
        N(1,j) = lbound(array3D,dim=j)
        N(2,j) = ubound(array3D,dim=j)
      enddo
      read(unit,iostat=err)(((array3D(i1,i2,i3)%rs(:),i1=N(1,1),N(2,1)),i2=N(1,2),N(2,2)),i3=N(1,3),N(2,3))
      err = read_vector(array3D=array3D%rv,unit=unit)
      read(unit,iostat=err)array3D%re
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction read_conservative

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
  !> @brief Subroutine for initializing Type_Conservative allocatable variables.
  elemental subroutine init(cons,Ns,cons0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(INOUT)::        cons  !< Conservative initialized data.
  integer(I_P),             intent(IN), optional:: Ns    !< Number of species.
  type(Type_Conservative),  intent(IN), optional:: cons0 !< Conservative inizialization data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(Ns)) then
    if (allocated(cons%rs)) deallocate(cons%rs) ; allocate(cons%rs(1:Ns)) ; cons%rs = 0._R8P
  elseif (present(cons0)) then
    if (allocated(cons%rs)) deallocate(cons%rs) ; allocate(cons%rs(1:size(cons0%rs))) ; cons = cons0
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  !> @brief Subroutine for freeing the memory of Type_Conservative allocatable variables.
  elemental subroutine free_cons(cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(INOUT):: cons !< Conservative data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons%rs)) deallocate(cons%rs)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_cons

  !> @brief Function for converting derived type Type_Conservative to array.
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

  !> @brief Function for printing in a pretty ascii format the components of type Type_Conservative.
  !> @return \b err integer(I_P) variable for error trapping.
  function pprint(cons,myrank,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Conservative), intent(IN)::           cons   !< Conservatives.
  integer(I_P),             intent(IN), optional:: myrank !< Actual rank process.
  integer(I4P),             intent(IN)::           unit   !< Logic unit.
  integer(I_P)::                                   err    !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                                   Ns     !< Number of species.
  integer(I_P)::                                   s      !< Species counter.
  character(DI_P)::                                rks    !< String containing myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,0_I_P)) ; if (present(myrank)) rks = 'rank'//trim(str(.true.,myrank))
  Ns  = size(cons%rs)
  do s=1,Ns
    write(unit,'(A)',iostat=err)trim(rks)//' rs('//trim(str(.true.,s))//')='//str(n=cons%rs(s))
  enddo
    write(unit,'(A)',iostat=err)trim(rks)//' rv(x)='//str(n=cons%rv%x)
    write(unit,'(A)',iostat=err)trim(rks)//' rv(y)='//str(n=cons%rv%y)
    write(unit,'(A)',iostat=err)trim(rks)//' rv(z)='//str(n=cons%rv%z)
    write(unit,'(A)',iostat=err)trim(rks)//' re   ='//str(n=cons%re  )
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction pprint

  ! Assignment (=)
  elemental subroutine assign_cons(cons1,cons2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between two cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons1
  type(Type_Conservative), intent(IN)::    cons2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons1%rs = cons2%rs
  cons1%rv = cons2%rv
  cons1%re = cons2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_cons

#ifdef r16p
  elemental subroutine assign_ScalR16P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R16P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  real(R16P),              intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons%rs = real(scal,R8P)
  cons%rv = real(scal,R8P)
  cons%re = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR16P
#endif

  elemental subroutine assign_ScalR8P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R8P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  real(R8P),               intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons%rs = real(scal,R8P)
  cons%rv = real(scal,R8P)
  cons%re = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR8P

  elemental subroutine assign_ScalR4P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R4P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  real(R4P),               intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons%rs = real(scal,R8P)
  cons%rv = real(scal,R8P)
  cons%re = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR4P

  elemental subroutine assign_ScalI8P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I8P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  integer(I8P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons%rs = real(scal,R8P)
  cons%rv = real(scal,R8P)
  cons%re = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI8P

  elemental subroutine assign_ScalI4P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I4P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  integer(I4P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons%rs = real(scal,R8P)
  cons%rv = real(scal,R8P)
  cons%re = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI4P

  elemental subroutine assign_ScalI2P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I2P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  integer(I2P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons%rs = real(scal,R8P)
  cons%rv = real(scal,R8P)
  cons%re = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI2P

  elemental subroutine assign_ScalI1P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I1P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  integer(I1P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons%rs = real(scal,R8P)
  cons%rv = real(scal,R8P)
  cons%re = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI1P

  ! Multiplication (*)
  elemental function cons_mul_cons(cons1,cons2) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply (by components) conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN)::  cons1 ! First cons obj.
  type(Type_Conservative), intent(IN)::  cons2 ! Second cons obj.
  type(Type_Conservative)::              mul   ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons1%rs)))
  mul%rs = cons1%rs * cons2%rs
  mul%rv = cons1%rv * cons2%rv
  mul%re = cons1%re * cons2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_cons

#ifdef r16p
  elemental function ScalR16P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R16P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),              intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_mul_cons

  elemental function cons_mul_ScalR16P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R16P),              intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalR16P
#endif

  elemental function ScalR8P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_mul_cons

  elemental function cons_mul_ScalR8P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R8P),               intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalR8P

  elemental function ScalR4P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_mul_cons

  elemental function cons_mul_ScalR4P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R4P),               intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalR4P

  elemental function ScalI8P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_mul_cons

  elemental function cons_mul_ScalI8P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalI8P

  elemental function ScalI4P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_mul_cons

  elemental function cons_mul_ScalI4P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalI4P

  elemental function ScalI2P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I2P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_mul_cons

  elemental function cons_mul_ScalI2P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalI2P

  elemental function ScalI1P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I1P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_mul_cons

  elemental function cons_mul_ScalI1P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R8P) * cons%rs
  mul%rv = real(scal,R8P) * cons%rv
  mul%re = real(scal,R8P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalI1P

  ! Division (/)
  elemental function cons_div_cons(cons1,cons2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide (by components) conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN)::  cons1 ! First cons obj.
  type(Type_Conservative), intent(IN)::  cons2 ! Second cons obj.
  type(Type_Conservative)::              div   ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%rs(1:size(cons1%rs)))
  div%rs = cons1%rs / cons2%rs
  div%rv = cons1%rv / cons2%rv
  div%re = cons1%re / cons2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_cons

#ifdef r16p
  elemental function cons_div_ScalR16P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R16P),              intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R8P)
  div%rv = cons%rv / real(scal,R8P)
  div%re = cons%re / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalR16P
#endif

  elemental function cons_div_ScalR8P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R8P),               intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R8P)
  div%rv = cons%rv / real(scal,R8P)
  div%re = cons%re / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalR8P

  elemental function cons_div_ScalR4P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R4P),               intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R8P)
  div%rv = cons%rv / real(scal,R8P)
  div%re = cons%re / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalR4P

  elemental function cons_div_ScalI8P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R8P)
  div%rv = cons%rv / real(scal,R8P)
  div%re = cons%re / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalI8P

  elemental function cons_div_ScalI4P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R8P)
  div%rv = cons%rv / real(scal,R8P)
  div%re = cons%re / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalI4P

  elemental function cons_div_ScalI2P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R8P)
  div%rv = cons%rv / real(scal,R8P)
  div%re = cons%re / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalI2P

  elemental function cons_div_ScalI1P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R8P)
  div%rv = cons%rv / real(scal,R8P)
  div%re = cons%re / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalI1P

  ! Sum (+)
  elemental function positive_cons(cons) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary + to a conservative objecy.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(pos%rs(1:size(cons%rs)))
  pos%rs =  + cons%rs
  pos%rv =  + cons%rv
  pos%re =  + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_cons

  elemental function cons_sum_cons(cons1,cons2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum (by components) conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN)::  cons1 ! First cons obj.
  type(Type_Conservative), intent(IN)::  cons2 ! Second cons obj.
  type(Type_Conservative)::              summ  ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons1%rs)))
  summ%rs = cons1%rs + cons2%rs
  summ%rv = cons1%rv + cons2%rv
  summ%re = cons1%re + cons2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_cons

#ifdef r16p
  elemental function ScalR16P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R16P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),              intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sum_cons

  elemental function cons_sum_ScalR16P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R16P),              intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalR16P
#endif

  elemental function ScalR8P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sum_cons

  elemental function cons_sum_ScalR8P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R8P),               intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalR8P

  elemental function ScalR4P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sum_cons

  elemental function cons_sum_ScalR4P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R4P),               intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalR4P

  elemental function ScalI8P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sum_cons

  elemental function cons_sum_ScalI8P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalI8P

  elemental function ScalI4P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sum_cons

  elemental function cons_sum_ScalI4P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalI4P

  elemental function ScalI2P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I2P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sum_cons

  elemental function cons_sum_ScalI2P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalI2P

  elemental function ScalI1P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I1P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sum_cons

  elemental function cons_sum_ScalI1P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R8P) + cons%rs
  summ%rv = real(scal,R8P) + cons%rv
  summ%re = real(scal,R8P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalI1P

  ! Subtraction (-)
    elemental function negative_cons(cons) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary - to a conservative objecy.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(neg%rs(1:size(cons%rs)))
  neg%rs =  - cons%rs
  neg%rv =  - cons%rv
  neg%re =  - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_cons

  elemental function cons_sub_cons(cons1,cons2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract (by components) conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN)::  cons1 ! First cons obj.
  type(Type_Conservative), intent(IN)::  cons2 ! Second cons obj.
  type(Type_Conservative)::              sub  ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons1%rs)))
  sub%rs = cons1%rs - cons2%rs
  sub%rv = cons1%rv - cons2%rv
  sub%re = cons1%re - cons2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_cons

#ifdef r16p
  elemental function ScalR16P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R16P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),              intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R8P) - cons%rs
  sub%rv = real(scal,R8P) - cons%rv
  sub%re = real(scal,R8P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sub_cons

  elemental function cons_sub_ScalR16P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R16P),              intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R8P)
  sub%rv = cons%rv - real(scal,R8P)
  sub%re = cons%re - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalR16P
#endif

  elemental function ScalR8P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R8P) - cons%rs
  sub%rv = real(scal,R8P) - cons%rv
  sub%re = real(scal,R8P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sub_cons

  elemental function cons_sub_ScalR8P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R8P),               intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R8P)
  sub%rv = cons%rv - real(scal,R8P)
  sub%re = cons%re - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalR8P

  elemental function ScalR4P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R8P) - cons%rs
  sub%rv = real(scal,R8P) - cons%rv
  sub%re = real(scal,R8P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sub_cons

  elemental function cons_sub_ScalR4P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R4P),               intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R8P)
  sub%rv = cons%rv - real(scal,R8P)
  sub%re = cons%re - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalR4P

  elemental function ScalI8P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R8P) - cons%rs
  sub%rv = real(scal,R8P) - cons%rv
  sub%re = real(scal,R8P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sub_cons

  elemental function cons_sub_ScalI8P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R8P)
  sub%rv = cons%rv - real(scal,R8P)
  sub%re = cons%re - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalI8P

  elemental function ScalI4P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R8P) - cons%rs
  sub%rv = real(scal,R8P) - cons%rv
  sub%re = real(scal,R8P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sub_cons

  elemental function cons_sub_ScalI4P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R8P)
  sub%rv = cons%rv - real(scal,R8P)
  sub%re = cons%re - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalI4P

  elemental function ScalI2P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I2P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R8P) - cons%rs
  sub%rv = real(scal,R8P) - cons%rv
  sub%re = real(scal,R8P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sub_cons

  elemental function cons_sub_ScalI2P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R8P)
  sub%rv = cons%rv - real(scal,R8P)
  sub%re = cons%re - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalI2P

  elemental function ScalI1P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I1P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R8P) - cons%rs
  sub%rv = real(scal,R8P) - cons%rv
  sub%re = real(scal,R8P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sub_cons

  elemental function cons_sub_ScalI1P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R8P)
  sub%rv = cons%rv - real(scal,R8P)
  sub%re = cons%re - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalI1P

  !> Function for computing the dot product of 2 conservative variables arrays.
  pure function cons_dot_cons(cons1,cons2) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons1(:) ! First conservative variables.
  type(Type_Conservative), intent(IN):: cons2(:) ! Second conservative variables.
  type(Type_Conservative)::             dot      ! Dot product.
  integer(I4P)::                        i        ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons1(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(cons1),size(cons2))
    dot = dot + cons1(i)*cons2(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_dot_cons

  !> Function for computing the dot product of a scalar (R8P) array for a conservative variables one.
  pure function scalR8_dot_cons(scal,cons) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),               intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction scalR8_dot_cons

  !> Function for computing the dot product of a scalar (R4P) array for a conservative variables one.
  pure function scalR4_dot_cons(scal,cons) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),               intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction scalR4_dot_cons

  !> Function for computing the dot product of a scalar (I8P) array for a conservative variables one.
  pure function scalI8_dot_cons(scal,cons) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),            intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction scalI8_dot_cons

  !> Function for computing the dot product of a scalar (I4P) array for a conservative variables one.
  pure function scalI4_dot_cons(scal,cons) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),            intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction scalI4_dot_cons

  !> Function for computing the dot product of a scalar (I2P) array for a conservative variables one.
  pure function scalI2_dot_cons(scal,cons) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),            intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction scalI2_dot_cons

  !> Function for computing the dot product of a scalar (I1P) array for a conservative variables one.
  pure function scalI1_dot_cons(scal,cons) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),            intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction scalI1_dot_cons

  !> Function for computing the dot product of a conservative variables one for a scalar (R8P) array.
  pure function cons_dot_scalR8(cons,scal) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  real(R8P),               intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_dot_scalR8

  !> Function for computing the dot product of a conservative variables one for a scalar (R4P) array.
  pure function cons_dot_scalR4(cons,scal) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  real(R4P),               intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_dot_scalR4

  !> Function for computing the dot product of a conservative variables one for a scalar (I8P) array.
  pure function cons_dot_scalI8(cons,scal) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  integer(I8P),            intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_dot_scalI8

  !> Function for computing the dot product of a conservative variables one for a scalar (I4P) array.
  pure function cons_dot_scalI4(cons,scal) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  integer(I4P),            intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_dot_scalI4

  !> Function for computing the dot product of a conservative variables one for a scalar (I2P) array.
  pure function cons_dot_scalI2(cons,scal) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  integer(I2P),            intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_dot_scalI2

  !> Function for computing the dot product of a conservative variables one for a scalar (I1P) array.
  pure function cons_dot_scalI1(cons,scal) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons(:) ! Conservative variables.
  integer(I1P),            intent(IN):: scal(:) ! Scalar array.
  type(Type_Conservative)::             dot     ! Dot product.
  integer(I4P)::                        i       ! Array elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(dot%rs(1:size(cons(1)%rs)))
  dot = 0._R8P
  do i=1,min(size(scal),size(cons))
    dot = dot + scal(i)*cons(i)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_dot_scalI1
  !> @}
endmodule Data_Type_Conservative
