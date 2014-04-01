!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_CompiledCodeDerivedType Data_Type_CompiledCode
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_CompiledCodePrivateProcedure Data_Type_CompiledCode
!> @}

!> This module contains the definition of Type_CompiledCode and its procedures.
!> This derived type has useful checking the compilation options used for compiling the codes.
module Data_Type_CompiledCode
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type contains the compilation options used for compiling the codes.
!> @ingroup Data_Type_CompiledCodeDerivedType
type, public:: Type_CompiledCode
  character(len=:), allocatable:: compiler !< Compiler used.
  character(len=:), allocatable:: compflag !< Compilation flags.
  character(len=:), allocatable:: compproc !< Pre-processing compilation flags.
  character(len=:), allocatable:: complibs !< External libraries used.
  contains
    procedure:: init              ! Procedure for initializing code compilation options.
    procedure:: print => print_cc ! Procedure for printing code compilation options.
endtype Type_CompiledCode
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_CompiledCodePrivateProcedure
  !> @{
  !> @brief Procedure for initializing code compilation options.
#define _STR_EXPAND(tok) #tok
#define _STR(tok) _STR_EXPAND(tok)
  subroutine init(cc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_CompiledCode), intent(INOUT):: cc !< Code compilation options.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef _COMPILER
  cc%compiler = _STR(_COMPILER)
#endif
#ifdef _COMPFLAG
  cc%compflag = _STR(_COMPFLAG)
#endif
#ifdef _COMPPROC
  cc%compproc = _STR(_COMPPROC)
#endif
#ifdef _COMPLIBS
  cc%complibs = _STR(_COMPLIBS)
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  !> @brief Procedure for printing code compilation options.
  subroutine print_cc(cc,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_CompiledCode), intent(INOUT):: cc      !< Code compilation options.
  character(*), optional,   intent(IN)::    pref    !< Prefixing string.
  integer(I4P), optional,   intent(OUT)::   iostat  !< IO error.
  character(*), optional,   intent(OUT)::   iomsg   !< IO error message.
  integer(I4P),             intent(IN)::    unit    !< Logic unit.
  character(len=:), allocatable::           prefd   !< Prefixing string.
  integer(I4P)::                            iostatd !< IO error.
  character(500)::                          iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  if (allocated(cc%compiler)) write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Compiler       : '//trim(cc%compiler)
  if (allocated(cc%compflag)) write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Comp. flags    : '//trim(cc%compflag)
  if (allocated(cc%compproc)) write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Pre-proc flags : '//trim(cc%compproc)
  if (allocated(cc%complibs)) write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Libraries used : '//trim(cc%complibs)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_cc
  !> @}
endmodule Data_Type_CompiledCode
