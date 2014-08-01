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
USE IR_Precision                                            ! Integers and reals precision definition.
USE Data_Type_Block_Dimensions, only: Type_Block_Dimensions ! Definition of Type_Block_Dimensions.
USE Data_Type_Tree,             only: Type_Tree             ! Definition of Type_Tree.
USE Lib_Math,                   only: digit                 ! Procedure for computing the significant digits of a number.
#ifdef OPENMP
USE OMP_LIB                                                 ! OpenMP runtime library.
#endif
#ifdef _MPI
USE MPI                                                     ! MPI runtime library.
#endif
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
  character(len=:), allocatable:: rks              !< String containing myrank with a pretty format.
  type(Type_Tree)::               BPmap            !< Tree of block/process map, ID-proc map.
  contains
    procedure:: free          => free_parallel          ! Procedure for freeing dynamic memory.
    procedure:: init          => init_parallel          ! Procedure for initializing dynamic memory.
    procedure:: is_master     => is_master_parallel     ! Procedure for checking if current process is the master one (0).
    procedure:: is_process    => is_process_parallel    ! Procedure for checking if current process is the p-th one.
    procedure:: compute_BPmap => compute_BPmap_parallel ! Procedure for computing BPmap.
    !procedure:: print         => print_parallel         ! Procedure for printing parallel data.
    final::     finalize                                ! Procedure for freeing memory when finalizing.
endtype Type_Parallel
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_ParallelPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_parallel(parallel)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Parallel), intent(INOUT):: parallel !< Parallel data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  parallel%myrank   = 0_I4P
  parallel%Nthreads = 1_I4P
  parallel%Nproc    = 1_I4P
  if (allocated(parallel%rks)) deallocate(parallel%rks)
  call parallel%BPmap%free
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
  call parallel%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for allocating dynamic memory, except BPmap one.
  elemental subroutine init_parallel(parallel)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Parallel), intent(INOUT):: parallel !< Parallel data.
  integer(I4P)::                        err      !< Error traping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call parallel%free
  ! initializing parallel environments
  associate(myrank=>parallel%myrank,Nthreads=>parallel%Nthreads,Nproc=>parallel%Nproc)
#ifdef _MPI
  call MPI_INIT(err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,err)
#endif
#ifdef OPENMP
  !$OMP PARALLEL      &
  !$OMP DEFAULT(none) &
  !$OMP SHARED(Nthreads)
  Nthreads = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
#endif
  endassociate
  parallel%rks = 'rank'//trim(strz(digit(parallel%Nproc),parallel%myrank))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_parallel

  !> @brief Procedure for checking if current process is the master one (0).
  elemental function is_master_parallel(parallel) result(I_am_master)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Parallel), intent(IN):: parallel    !< Parallel data.
  logical::                          I_am_master !< Flag for inquiring if current process is the master one.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  I_am_master = (parallel%myrank==0)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_master_parallel

  !> @brief Procedure for checking if current process is the p-th one.
  elemental function is_process_parallel(parallel,p) result(I_am_pth)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Parallel), intent(IN):: parallel !< Parallel data.
  integer(I4P),         intent(IN):: p        !< Rank of p-th process.
  logical::                          I_am_pth !< Flag for inquiring if current process is the p-th one.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  I_am_pth = (parallel%myrank==p)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_process_parallel

  !> @brief Procedure for computing BPmap.
  subroutine compute_BPmap_parallel(parallel,block_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Parallel), intent(INOUT)::  parallel         !< Parallel data.
  type(Type_Tree),      intent(IN)::     block_dims       !< Mesh dimensions.
  type(Type_Block_Dimensions), pointer:: blkdims          !< Block dimensions pointer for scanning block_dims tree.
  integer(I8P)::                         ID               !< Counter.
  integer(I8P)::                         mesh_weight      !< Weight of mesh, sum over all blocks.
  integer(I8P)::                         ideal_prc_weight !< Ideal Weight of each processes, optimal balance.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(Nproc=>parallel%Nproc,BPmap=>parallel%BPmap)
    call BPmap%init(source=block_dims)
    mesh_weight = 0
    do while(block_dims%loopID(ID=ID))
      blkdims => block_dims%dat(ID=ID)
      mesh_weight = mesh_weight + blkdims%Ncells()
    enddo
    ideal_prc_weight = mesh_weight / Nproc
    ! cazzo inserire load balance
    do while(block_dims%loopID(ID=ID))
      call BPmap%put(ID=ID,d=0_I4P) ! cazzo tutti i blocchi sul proc 0... seriale
    enddo
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_BPmap_parallel

 !!> @brief Procedure for pretty printing parallel data.
 !subroutine print_procmap(parallel,pref,iostat,iomsg,unit)
 !!---------------------------------------------------------------------------------------------------------------------------------
 !implicit none
 !class(Type_Parallel),   intent(IN)::  parallel !< parallel-level data.
 !character(*), optional, intent(IN)::  pref     !< Prefixing string.
 !integer(I4P), optional, intent(OUT):: iostat   !< IO error.
 !character(*), optional, intent(OUT):: iomsg    !< IO error message.
 !integer(I4P),           intent(IN)::  unit     !< Logic unit.
 !character(len=:), allocatable::       prefd    !< Prefixing string.
 !integer(I4P)::                        iostatd  !< IO error.
 !character(500)::                      iomsgd   !< Temporary variable for IO error message.
 !integer(I4P)::                        Nb_tot   !< Totabl blocks number.
 !integer(I4P)::                        Nb       !< Local blocks number.
 !integer(I4P)::                        b        !< Blocks counter.
 !!---------------------------------------------------------------------------------------------------------------------------------
 !
 !!---------------------------------------------------------------------------------------------------------------------------------
 !prefd = '' ; if (present(pref)) prefd = pref
 !associate(Nproc=>parallel%Nproc,procmap=>parallel%procmap,blockmap=>parallel%blockmap)
 !  Nb_tot = size(procmap,dim=1) ; Nb = size(blockmap,dim=1)
 !  write(  unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Processes/blocks map'
 !  do b=1,Nb_tot
 !    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Block '//trim(strz(Nb_tot,b))//&
 !                                                                 ' process '//trim(strz(Nproc,procmap(b)))
 !  enddo
 !  write(  unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Local/global blocks numeration'
 !  do b=1,Nb
 !    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Local block '//trim(strz(Nb,b))//&
 !                                                                 ' global numeration '//trim(strz(Nb_tot,blockmap(b)))
 !  enddo
 !endassociate
 !if (present(iostat)) iostat = iostatd
 !if (present(iomsg))  iomsg  = iomsgd
 !return
 !!---------------------------------------------------------------------------------------------------------------------------------
 !endsubroutine print_procmap
  !> @}
endmodule Data_Type_Parallel
