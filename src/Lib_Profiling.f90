!> @ingroup Library
!> @{
!> @defgroup Lib_ProfilingLibrary Lib_Profiling
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup Lib_ProfilingPrivateVarPar Lib_Profiling
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_ProfilingPublicProcedure Lib_Profiling
!> @}

!> This module contains procedures and data for profiling parts of the code.
!> This is a library module.
!> @ingroup Lib_ProfilingLibrary
module Lib_Profiling
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision   ! Integers and reals precision definition. module Lib_Profiling
USE Data_Type_Time ! Definition of Type_Time.
USE Lib_IO_Misc    ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: profile
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Lib_ProfilingPrivateVarPar
!> @{
real(R8P)::                 instant0 = 0._R8P !< The Crono starting instant used for profing the code.
integer(I_P)::              Np                !< Number of parts of the code profiled.
real(R8P),    allocatable:: partial(:,:)      !< Partial times used for profing the code [1:2,1:Np].
integer(I_P), allocatable:: unitprofile(:)    !< Logic units for profiling files [1:Np].
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Subroutine for profiling the code.
!> @ingroup Lib_ProfilingPublicProcedure
  subroutine profile(Np,header,p,pstart,pstop,finalize,myrank)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN), optional:: Np        !< Number of parts of the code profiled. Initializing the profiling if passed.
  character(*), intent(IN), optional:: header(:) !< Headers of profiling files: can be passed only with Np.
  integer(I_P), intent(IN), optional:: p         !< Number of part of the code currently profiled.
  logical,      intent(IN), optional:: pstart    !< Flag for start the profiling of the current profiled part of the code.
  logical,      intent(IN), optional:: pstop     !< Flag for stop the profiling of the current profiled part of the code.
  logical,      intent(IN), optional:: finalize  !< Flag for finilizing the profiling of the code.
  integer(I_P), intent(IN)::           myrank    !< Current rank process.
  integer(I_P)::                       u         !< Unit files counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(Np)) then ! the profiling is starting
    if (allocated(partial    )) deallocate(partial    ) ; allocate(partial(1:2,1:Np))
    if (allocated(unitprofile)) deallocate(unitprofile) ; allocate(unitprofile(1:Np))
    instant0 = Crono(start=.true.)
    do u=1,Np
      open(unit=Get_Unit(unitprofile(u)),file='profiling-'//trim(strz(3,u))//'.p'//trim(strz(3,myrank))//'.dat')
      write(unitprofile(u),'(A)')'VARIABLES="Time Consumed"'
      if (present(header)) write(unitprofile(u),'(A)')trim(header(u))
    enddo
  else
    if (present(pstart)) then
      partial(1,p) = Crono()
    elseif (present(pstop)) then
      partial(2,p) = Crono(instant1=partial(1,p))
      write(unitprofile(p),'(A)')trim(str(n=partial(2,p)))
    elseif (present(finalize)) then
      if (allocated(partial)) deallocate(partial)
      do u=1,Np
        close(unitprofile(u))
      enddo
      if (allocated(unitprofile)) deallocate(unitprofile)
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine profile
endmodule Lib_Profiling
