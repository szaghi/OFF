!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_ParallelDerivedType Data_Type_Parallel
!> Module definition of Type_Parallel
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_ParallelInterface Data_Type_Parallel
!> Module definition of Type_Parallel
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_ParallelPublicProcedure Data_Type_Parallel
!> Module definition of Type_Parallel
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_ParallelPrivateProcedure Data_Type_Parallel
!> Module definition of Type_Parallel
!> @}

!> @brief Module Data_Type_Parallel contains the definition of Type_Parallel, that defines the main parallel data informations.
module Data_Type_Parallel
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_Parallel.
!> @ingroup Data_Type_ParallelDerivedType
type, public:: Type_Parallel
  integer(I4P)::                  myrank   = 0_I4P !< Rank of the actual process.
  integer(I4P)::                  Nthreads = 1_I4P !< Number of OpenMP threads.
  integer(I4P)::                  Nproc    = 1_I4P !< Number of MPI processes.
  integer(I4P),     allocatable:: procmap(:)       !< Processes/blocks map    [1:Nb_tot].
  integer(I4P),     allocatable:: blockmap(:)      !< Local/global blocks map [1:Nb].
  character(len=:), allocatable:: rks              !< String containing myrank with a pretty format.
  contains
    procedure:: free       => free_parallel       ! Procedure for freeing dynamic memory.
    procedure:: alloc      => alloc_parallel      ! Procedure for allocating dynamic memory.
    procedure:: set_serial => set_serial_parallel ! Procedure for setting data in serial framework [Nproc=1].
    procedure:: print      => print_procmap       ! Procedure for printing the processes/blocks map and local/global blocks map.
    procedure:: set_rks                           ! Procedure for setting rks string accordingly to parallel data.
    final::     finalize                          ! Procedure for freeing memory when finalizing.
endtype Type_Parallel
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_ParallelPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_parallel(parallel,procmap,blockmap,rks)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Parallel), intent(INOUT):: parallel !< Parallel data.
  logical, optional,    intent(IN)::    procmap  !< Inquiring flag.
  logical, optional,    intent(IN)::    blockmap !< Inquiring flag.
  logical, optional,    intent(IN)::    rks      !< Inquiring flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(procmap ).AND.allocated(parallel%procmap )) deallocate(parallel%procmap )
  if (present(blockmap).AND.allocated(parallel%blockmap)) deallocate(parallel%blockmap)
  if (present(rks     ).AND.allocated(parallel%rks     )) deallocate(parallel%rks     )
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_parallel

  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine finalize(parallel)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Parallel), intent(INOUT):: parallel !< Parallel data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call parallel%free(procmap=.true.,blockmap=.true.,rks=.true.)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc_parallel(parallel,Nb_tot,Nb)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Parallel),   intent(INOUT):: parallel !< Parallel data.
  integer(I4P), optional, intent(IN)::    Nb_tot   !< Number of total blocks (sum over each process).
  integer(I4P), optional, intent(IN)::    Nb       !< Number of blocks of myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(Nb_tot)) then
    call parallel%free(procmap=.true. ) ; allocate(parallel%procmap (1:Nb_tot)) ; parallel%procmap  = 0_I4P
  endif
  if (present(Nb)) then
    call parallel%free(blockmap=.true.) ; allocate(parallel%blockmap(1:Nb    )) ; parallel%blockmap = 0_I4P
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_parallel

  !> @brief Procedure for setting data in serial framework [Nproc=1].
  !> @note In a serial framework (number of MPI processes equals to 1) the procmap array has all elements equal to 0, namely all
  !> blocks are assigned to processes of rank 0, while the blockmap array is the identity, namely blockmap(b)=b for all blocks.
  elemental subroutine set_serial_parallel(parallel,Nb_tot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Parallel), intent(INOUT):: parallel !< Parallel data.
  integer(I4P),         intent(IN)::    Nb_tot   !< Number of total blocks (sum over each process).
  integer(I4P)::                        b        !< Blocks counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call parallel%alloc(Nb_tot=Nb_tot,Nb=Nb_tot)
  do b=1,Nb_tot
    parallel%blockmap(b) = b ! the blocks map is identity
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_serial_parallel

  !> @brief Procedure for printing to stdout the processes/blocks map.
  subroutine print_procmap(parallel,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Parallel),   intent(IN)::  parallel !< parallel-level data.
  character(*), optional, intent(IN)::  pref     !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat   !< IO error.
  character(*), optional, intent(OUT):: iomsg    !< IO error message.
  integer(I4P),           intent(IN)::  unit     !< Logic unit.
  character(len=:), allocatable::       prefd    !< Prefixing string.
  integer(I4P)::                        iostatd  !< IO error.
  character(500)::                      iomsgd   !< Temporary variable for IO error message.
  integer(I4P)::                        Nb_tot   !< Totabl blocks number.
  integer(I4P)::                        Nb       !< Local blocks number.
  integer(I4P)::                        b        !< Blocks counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  associate(Nproc=>parallel%Nproc,procmap=>parallel%procmap,blockmap=>parallel%blockmap)
    Nb_tot = size(procmap,dim=1) ; Nb = size(blockmap,dim=1)
    write(  unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Processes/blocks map'
    do b=1,Nb_tot
      write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Block '//trim(strz(Nb_tot,b))//&
                                                                   ' process '//trim(strz(Nproc,procmap(b)))
    enddo
    write(  unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Local/global blocks numeration'
    do b=1,Nb
      write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Local block '//trim(strz(Nb,b))//&
                                                                   ' global numeration '//trim(strz(Nb_tot,blockmap(b)))
    enddo
  endassociate
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_procmap

  !> @brief Procedure for setting rks string accordingly to parallel data.
  subroutine set_rks(parallel)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Parallel), intent(INOUT):: parallel !< Parallel data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  parallel%rks = 'rank'//trim(strz(parallel%Nproc,parallel%myrank))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_rks
  !> @}
endmodule Data_Type_Parallel
