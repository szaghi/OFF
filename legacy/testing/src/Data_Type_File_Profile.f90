!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_ProfileDerivedType Data_Type_File_Profile
!> Module definition of Type_File_Profile
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_File_ProfileInterface Data_Type_File_Profile
!> Module definition of Type_File_Profile
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_File_ProfilePublicProcedure Data_Type_File_Profile
!> Module definition of Type_File_Profile
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_ProfilePrivateProcedure Data_Type_File_Profile
!> Module definition of Type_File_Profile
!> @}

!> @ingroup Data_Type_File_ProfileDerivedType
module Data_Type_File_Profile
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                              ! Integers and reals precision definition.
USE Data_Type_File_Base, only: Type_File_Base ! Definition of Type_File_Base.
USE Lib_IO_Misc,         only: Get_Unit       ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: profile
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Profile.
!> @ingroup Data_Type_File_ProfileDerivedType
type, public, extends(Type_File_Base):: Type_File_Profile
  integer(I4P)::              Npp = 0_I4P    !< Number of parts of the code profiled.
  integer(I8P), allocatable:: partial(:,:)   !< Partial times used for profing the code [1:2,1:Npp].
  integer(I4P), allocatable:: unitprofile(:) !< Logic units for profiling files [1:Npp].
  integer(I4P), allocatable:: tictoc(:)      !< Counter of files access [1:Npp].
  integer(I4P)::              pp = 0_I4P     !< Index of current profiled part of the code.
  contains
    procedure:: free  => free_profile  ! Procedure for freeing memory.
    procedure:: alloc => alloc_profile ! Procedure for allocating memory.
    procedure:: profile                ! Procedure for profiling the code.
endtype Type_File_Profile
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_ProfilePrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic data of Type_File_Profile variables.
  elemental subroutine free_profile(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Profile), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free_base
  if (allocated(file_d%partial    )) deallocate(file_d%partial    )
  if (allocated(file_d%unitprofile)) deallocate(file_d%unitprofile)
  if (allocated(file_d%tictoc     )) deallocate(file_d%tictoc     )
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_profile

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize_profile(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_File_Profile), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize_profile

  !> Subroutine for allocating dynamic data of Type_File_Profile variables.
  elemental subroutine alloc_profile(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Profile), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(file_d%partial    )) deallocate(file_d%partial    )
  if (allocated(file_d%unitprofile)) deallocate(file_d%unitprofile)
  if (allocated(file_d%tictoc     )) deallocate(file_d%tictoc     )
  allocate(file_d%partial(1:2,1:file_d%Npp)) ; file_d%partial     = 0_I8P
  allocate(file_d%unitprofile(1:file_d%Npp)) ; file_d%unitprofile = 0_I4P
  allocate(file_d%tictoc(     1:file_d%Npp)) ; file_d%tictoc      = 0_I4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_profile

  !> @brief Subroutine for profiling the code.
  subroutine profile(file_d,Np,header,p,pstart,pstop,finalize,myrank,Nthreads,Nproc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Profile), intent(INOUT):: file_d    !< File data.
  integer(I4P), optional,   intent(IN)::    Np        !< Number of parts of the code profiled. Initializing the profiling if passed.
  character(*), optional,   intent(IN)::    header(:) !< Headers of profiling files: can be passed only with Np.
  integer(I4P), optional,   intent(IN)::    p         !< Number of part of the code currently profiled.
  logical,      optional,   intent(IN)::    pstart    !< Flag for start the profiling of the current profiled part of the code.
  logical,      optional,   intent(IN)::    pstop     !< Flag for stop the profiling of the current profiled part of the code.
  logical,      optional,   intent(IN)::    finalize  !< Flag for finilizing the profiling of the code.
  integer(I4P),             intent(IN)::    myrank    !< Current rank process.
  integer(I4P),             intent(IN)::    Nthreads  !< Number of OpenMP threads.
  integer(I4P),             intent(IN)::    Nproc     !< Number of MPI processes.
  integer(I4P)::                            u         !< Unit files counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(Np)) then ! the profiling is starting
    file_d%Npp = Np
    call file_d%alloc
    do u=1,Np
        open(unit=Get_Unit(file_d%unitprofile(u)),&
             file=file_d%path_out//file_d%name//'-'//trim(strz(3,u))//'.p'//trim(strz(3,myrank))//'.prf')
      write(file_d%unitprofile(u),'(A)')'VARIABLES="Nthreads" "Nprocs" "TicToc" "Time Consumed"'
      if (present(header)) write(file_d%unitprofile(u),'(A)')trim(header(u))
    enddo
  else
    if (present(pstart)) then
      if (present(p)) file_d%pp = p
      call system_clock(file_d%partial(1,file_d%pp))
      file_d%tictoc(file_d%pp) = file_d%tictoc(file_d%pp) + 1_I4P
    elseif (present(pstop)) then
      if (present(p)) file_d%pp = p
      call system_clock(file_d%partial(2,file_d%pp))
      file_d%partial(2,file_d%pp) = file_d%partial(2,file_d%pp) - file_d%partial(1,file_d%pp)
      write(file_d%unitprofile(file_d%pp),'(A)')trim(str(n=Nthreads))//' '//&
                                                trim(str(n=Nproc))//' '//&
                                                trim(str(n=file_d%tictoc(file_d%pp)))//' '//&
                                                trim(str(n=file_d%partial(2,file_d%pp)))
      file_d%pp = file_d%pp + 1 ; if (file_d%pp > file_d%Npp) file_d%pp = 1_I4P
    elseif (present(finalize)) then
      do u=1,file_d%Npp
        close(file_d%unitprofile(u))
      enddo
      call file_d%free
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine profile
  !> @}
endmodule Data_Type_File_Profile
