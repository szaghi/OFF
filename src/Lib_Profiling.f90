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
USE IR_Precision                       ! Integers and reals precision definition.
USE Data_Type_Time                     ! Definition of Type_Time.
USE Lib_IO_Misc                        ! Procedures for IO and strings operations.
USE Lib_Parallel, only: Nthreads,Nproc ! Number of threads and procs.
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
!real(R8P)::                 instant0 = 0._R8P !< The Crono starting instant used for profing the code.
integer(I8P)::              instant0 = 0_I8P  !< The Crono starting instant used for profing the code.
integer(I_P)::              Npp = 0_I_P       !< Number of parts of the code profiled.
!real(R8P),    allocatable:: partial(:,:)      !< Partial times used for profing the code [1:2,1:Npp].
integer(I8P), allocatable:: partial(:,:)      !< Partial times used for profing the code [1:2,1:Npp].
integer(I_P), allocatable:: unitprofile(:)    !< Logic units for profiling files [1:Npp].
integer(I_P), allocatable:: tictoc(:)         !< Counter of files access [1:Npp].
integer(I_P)::              pp = 0_I_P        !< Index of current profiled part of the code.
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Subroutine for profiling the code.
  !> @ingroup Lib_ProfilingPublicProcedure
  subroutine profile(Np,fnamep,header,p,pstart,pstop,finalize,myrank)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN), optional:: Np        !< Number of parts of the code profiled. Initializing the profiling if passed.
  character(*), intent(IN), optional:: fnamep    !< Prefix name of profiling files: can be passed only with Np.
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
    Npp = Np
    !if (allocated(partial    )) deallocate(partial    ) ; allocate(partial(1:2,1:Np)) ; partial     = 0._R8P
    if (allocated(partial    )) deallocate(partial    ) ; allocate(partial(1:2,1:Np)) ; partial     = 0_I8P
    if (allocated(unitprofile)) deallocate(unitprofile) ; allocate(unitprofile(1:Np)) ; unitprofile = 0_I_P
    if (allocated(tictoc     )) deallocate(tictoc     ) ; allocate(tictoc(     1:Np)) ; tictoc      = 0_I_P
    !instant0 = Crono(start=.true.)
    call system_clock(instant0)
    do u=1,Np
      if (present(fnamep)) then
        open(unit=Get_Unit(unitprofile(u)),file=trim(fnamep)//'-'//trim(strz(3,u))//'.p'//trim(strz(3,myrank))//'.dat')
      else
        open(unit=Get_Unit(unitprofile(u)),file='profiling-'//trim(strz(3,u))//'.p'//trim(strz(3,myrank))//'.dat')
      endif
      write(unitprofile(u),'(A)')'VARIABLES="Nthreads" "Nprocs" "TicToc" "Time Consumed"'
      if (present(header)) write(unitprofile(u),'(A)')trim(header(u))
    enddo
  else
    if (present(pstart)) then
      if (present(p)) pp = p
      !partial(1,pp) = Crono()
      call system_clock(partial(1,pp))
      tictoc(pp) = tictoc(pp) + 1_I_P
    elseif (present(pstop)) then
      if (present(p)) pp = p
      !partial(2,pp) = Crono(instant1=partial(1,pp))
      call system_clock(partial(2,pp)) ; partial(2,pp) = partial(2,pp) - partial(1,pp)
      write(unitprofile(pp),'(A)')trim(str(n=Nthreads))//' '//trim(str(n=Nproc))//' '//&
                                  trim(str(n=tictoc(pp)))//' '//trim(str(n=partial(2,pp)))
      pp = pp + 1 ; if (pp > Npp) pp = 1_I_P
    elseif (present(finalize)) then
      do u=1,Npp
        close(unitprofile(u))
      enddo
      if (allocated(partial)) deallocate(partial)
      if (allocated(unitprofile)) deallocate(unitprofile)
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine profile
endmodule Lib_Profiling
