!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_ResidualsDerivedType Data_Type_File_Residuals
!> Module definition of Type_File_Residuals
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_File_ResidualsInterface Data_Type_File_Residuals
!> Module definition of Type_File_Residuals
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_File_ResidualsPublicProcedure Data_Type_File_Residuals
!> Module definition of Type_File_Residuals
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_ResidualsPrivateProcedure Data_Type_File_Residuals
!> Module definition of Type_File_Residuals
!> @}

!> @brief Module Data_Type_File_Residuals contains the definition of Type_File_Residuals, that is the main IBM options file.
module Data_Type_File_Residuals
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision        ! Integers and reals precision definition.
USE Data_Type_File_Base ! Definition of Type_File_Base.
USE Data_Type_OS        ! Definition of Type_OS.
USE Lib_IO_Misc         ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Residuals.
!> @ingroup Data_Type_File_ResidualsDerivedType
type, public, extends(Type_File_Base):: Type_File_Residuals
  contains
    procedure:: init => init_residuals ! Procedure for initializing residuals file.
endtype Type_File_Residuals
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_ResidualsPrivateProcedure
  !> @{
  !> @brief Procedure for initializing residuals file.
  function init_residuals(file_res,stderrpref,iomsg,OS,Nc,n,date) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Residuals), intent(INOUT):: file_res    !< File data.
  character(*), optional,     intent(IN)::    stderrpref  !< Prefixing string for stderr outputs.
  character(*), optional,     intent(OUT)::   iomsg       !< IO error message.
  type(Type_OS),              intent(IN)::    OS          !< Running architecture.
  integer(I4P),               intent(IN)::    Nc          !< Number of conservative variables (Nc = Ns + 4).
  integer(I8P),               intent(IN)::    n           !< Time steps counter.
  character(*),               intent(IN)::    date        !< Current date.
  integer(I4P)::                              err         !< Error trapping flag: 0 no errors, >0 error occurs.
  character(500)::                            iomsgd = '' !< Temporary variable for IO error message.
  character(500)::                            varname_res !< Variables name for the gnuplot residuals file.
  integer(I4P)::                              c           !< Counter.
  !--------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! creating the gnuplot script file for visualizing the residuals log file
  open(unit=Get_Unit(file_res%fileunit),file=trim(file_res%path_out)//'gplot_res',iostat=err,iomsg=iomsgd)
  if (err/=0) then
    if (present(iomsg)) iomsg=trim(iomsgd)
    return
  endif
  write(file_res%fileunit,'(A)',iostat=err,iomsg=iomsgd)'set xlabel "Iteration"'
  write(file_res%fileunit,'(A)',iostat=err,iomsg=iomsgd)'set ylabel "Residuals"'
  write(file_res%fileunit,'(A)',iostat=err,iomsg=iomsgd)'set log y'
  write(file_res%fileunit,'(A)',iostat=err,iomsg=iomsgd)"p 'residuals.log' u 1:3 w l title 'R1', "//char(92)
  do c=2,Nc-1
    write(file_res%fileunit,'(A)',iostat=err,iomsg=iomsgd)"  'residuals.log' u 1:"//trim(str(.true.,2+c))//" w l title 'R"//&
                                                           trim(str(.true.,c))//"', "//char(92)
  enddo
    write(file_res%fileunit,'(A)',iostat=err,iomsg=iomsgd)"  'residuals.log' u 1:"//trim(str(.true.,2+Nc))//" w l title 'R"//&
                                                          trim(str(.true.,Nc))//"'"
  close(file_res%fileunit,iostat=err,iomsg=iomsgd)
  ! initialize gnuplot log file of residuals
  file_res%filename='residuals.log'
  if (n>0) then
    if (present(stderrpref)) then
      err=file_res%open(stderrpref=stderrpref,iomsg=iomsgd,ascii=.true.,append=.true.,action='WRITE')
    else
      err=file_res%open(                      iomsg=iomsgd,ascii=.true.,append=.true.,action='WRITE')
    endif
  else
    if (present(stderrpref)) then
      err=file_res%open(stderrpref=stderrpref,iomsg=iomsgd,ascii=.true.,action='WRITE')
    else
      err=file_res%open(                      iomsg=iomsgd,ascii=.true.,action='WRITE')
    endif
    if (err/=0) then
      if (present(iomsg)) iomsg=trim(iomsgd)
      return
    endif
    ! initialize header
    do c=1,Nc-1
      varname_res = trim(varname_res)//'"R'//trim(str(.true.,c))//'",'
    enddo
    varname_res = trim(varname_res)//'"R'//trim(str(.true.,Nc))//'"'
    write(file_res%fileunit,'(A)',iostat=err,iomsg=iomsgd)'# L2 norm of residuals'
    write(file_res%fileunit,'(A)',iostat=err,iomsg=iomsgd)'# Simulation started on '//date
    write(file_res%fileunit,'(A)',iostat=err,iomsg=iomsgd)'# "n","t",'//adjustl(trim(varname_res))
  endif
  if (present(iomsg)) iomsg=trim(iomsgd)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction init_residuals
  !> @}
endmodule Data_Type_File_Residuals
