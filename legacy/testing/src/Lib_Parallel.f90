!> @ingroup Library
!> @{
!> @defgroup Lib_ParallelLibrary Lib_Parallel
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Lib_ParallelGlobalVarPar Lib_Parallel
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup Lib_ParallelPrivateVarPar Lib_Parallel
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_ParallelPublicProcedure Lib_Parallel
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_ParallelPrivateProcedure Lib_Parallel
!> @}

!> This module contains the definition of procedures for send/receive data among processes for parallel (MPI) operations.
!> It is based on MPI library.
!> @note The communications have a tag-shift (for make them unique) that assumes a maximum number of processes of 10000.
!> Increment this parameter if using more processes than 10000.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_ParallelLibrary
module Lib_Parallel
!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef _MPI
USE IR_Precision               ! Integers and reals precision definition.
USE Data_Type_BC, only: bc_adj ! Definition of Type_BC.
USE Data_Type_Parallel         ! Definition of Type_Parallel.
USE Data_Type_SBlock           ! Definition of Type_SBlock.
USE Lib_IO_Misc                ! Procedures for IO and strings operations.
USE MPI                        ! MPI runtime library.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
#ifdef _MPI
public:: init_MPI_maps
public:: print_MPI_maps
public:: prim_sendrecv
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef _MPI
!> @ingroup Lib_ParallelPrivateVarPar
!> @{
integer(I4P), parameter::   maxproc = 10000 !< Maximum number of processes used for communications tag shift.
integer(I4P)             :: gNcR            !< Global number of receive cells (sum(NcR)).
integer(I4P)             :: gNcS            !< Global number of send   cells (sum(NcS)).
integer(I4P), allocatable:: NcR(:,:)        !< Number of receive cells from each process [    0:Nproc-1,1:Nl].
integer(I4P), allocatable:: NcS(:,:)        !< Number of send   cells for  each process  [    0:Nproc-1,1:Nl].
integer(I4P), allocatable:: bbR(:,:,:)      !< Processes bounds of receive cells         [1:2,0:Nproc-1,1:Nl].
integer(I4P), allocatable:: bbS(:,:,:)      !< Processes bounds of send   cells          [1:2,0:Nproc-1,1:Nl].
integer(I4P), allocatable:: recvmap(:,:)    !< Receiving cells map of   myrank from other processes [1:4,1:gNcR].
integer(I4P), allocatable:: reqsmap(:,:)    !< Querying  cells map of   myrank for  other processes [1:4,1:gNcR].
integer(I4P), allocatable:: sendmap(:,:)    !< Sending  cells map from myrank to   other processes [1:4,1:gNcS].
real(R8P),    allocatable:: Precv(:,:)      !< Receiving buffer of primitive variable of myrank from other processes [1:Np,1:gNcR].
real(R8P),    allocatable:: Psend(:,:)      !< Sending  buffer of primitive variable of myrank for other  processes [1:Np,1:gNcS].
!> @}
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef _MPI
contains
  !> @ingroup Lib_ParallelPrivateProcedure
  !> @{
  !> @brief Procedure for safety allocation of NcR, NcS, bbR and bbS.
  subroutine alloc_SR(parallel,mesh_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Parallel),        intent(IN):: parallel  !< Parallel data.
  type(Type_Mesh_Dimensions), intent(IN):: mesh_dims !< Mesh dimensions.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(Nproc=>parallel%Nproc,Nl=>mesh_dims%Nl)
    if (allocated(NcR)) deallocate(NcR) ; allocate(NcR(    0:Nproc,1:Nl)) ; NcR = 0_I_P
    if (allocated(NcS)) deallocate(NcS) ; allocate(NcS(    0:Nproc,1:Nl)) ; NcS = 0_I_P
    if (allocated(bbR)) deallocate(bbR) ; allocate(bbR(1:2,0:Nproc,1:Nl)) ; bbR = 0_I_P
    if (allocated(bbS)) deallocate(bbS) ; allocate(bbS(1:2,0:Nproc,1:Nl)) ; bbS = 0_I_P
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_SR

  !> @brief Procedure for computing the number of cells that must be received from other processes than myrank.
  subroutine compute_NcR(parallel,mesh_dims,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Parallel),        intent(IN):: parallel     !< Parallel data.
  type(Type_Mesh_Dimensions), intent(IN):: mesh_dims    !< Mesh dimensions.
  type(Type_SBlock),          intent(IN):: block(1:,1:) !< Block-level data.
  integer(I_P)::                           l            !< Grid levels counter.
  integer(I_P)::                           proc         !< Processes counter.
  integer(I_P)::                           c            !< Cells counter.
  integer(I_P)::                           b            !< Blocks counter.
  integer(I_P)::                           i,j,k        !< Space counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing the number of receive cells of myrank for other processes; checking the equivalence of proc and myrank:
  !  if proc=myrank => NcR=0 => myrank doesn't communicate with itself
  do l=1,mesh_dims%Nl
    do b=1,mesh_dims%Nb
      associate(gc=>block(b,l)%dims%gc,Ni=>block(b,l)%dims%Ni,Nj=>block(b,l)%dims%Nj,Nk=>block(b,l)%dims%Nk,&
                Fi=>block(b,l)%Fi,Fj=>block(b,l)%Fj,Fk=>block(b,l)%Fk)
      ! i interfaces
        do k=1,Nk
          do j=1,Nj
          ! left i
            do i=0-gc(1),0
              if (Fi(i,j,k)%BC%tp==bc_adj) then
                if (parallel%procmap(Fi(i,j,k)%BC%adj%b)/=parallel%myrank) &
                  NcR(parallel%procmap(Fi(i,j,k)%BC%adj%b),l) = NcR(parallel%procmap(Fi(i,j,k)%BC%adj%b),l) + 1
            endif
          enddo
          ! right i
            do i=Ni,Ni+gc(2)
              if (Fi(i,j,k)%BC%tp==bc_adj) then
                if (parallel%procmap(Fi(i,j,k)%BC%adj%b)/=parallel%myrank) &
                  NcR(parallel%procmap(Fi(i,j,k)%BC%adj%b),l) = NcR(parallel%procmap(Fi(i,j,k)%BC%adj%b),l) + 1
              endif
          enddo
        enddo
      enddo
      ! j interfaces
        do k=1,Nk
        ! left j
          do j=0-gc(3),0
            do i=1,Ni
              if (Fj(i,j,k)%BC%tp==bc_adj) then
                if (parallel%procmap(Fj(i,j,k)%BC%adj%b)/=parallel%myrank) &
                  NcR(parallel%procmap(Fj(i,j,k)%BC%adj%b),l) = NcR(parallel%procmap(Fj(i,j,k)%BC%adj%b),l) + 1
              endif
          enddo
        enddo
        ! right j
          do j=Nj,Nj+gc(4)
            do i=1,Ni
              if (Fj(i,j,k)%BC%tp==bc_adj) then
                if (parallel%procmap(Fj(i,j,k)%BC%adj%b)/=parallel%myrank) &
                  NcR(parallel%procmap(Fj(i,j,k)%BC%adj%b),l) = NcR(parallel%procmap(Fj(i,j,k)%BC%adj%b),l) + 1
              endif
          enddo
        enddo
      enddo
      ! k interfaces
      ! left k
        do k=0-gc(5),0
          do j=1,Nj
            do i=1,Ni
              if (Fk(i,j,k)%BC%tp==bc_adj) then
                if (parallel%procmap(Fk(i,j,k)%BC%adj%b)/=parallel%myrank) &
                  NcR(parallel%procmap(Fk(i,j,k)%BC%adj%b),l) = NcR(parallel%procmap(Fk(i,j,k)%BC%adj%b),l) + 1
              endif
          enddo
        enddo
      enddo
      ! right k
        do k=Nk,Nk+gc(6)
          do j=1,Nj
            do i=1,Ni
              if (Fk(i,j,k)%BC%tp==bc_adj) then
                if (parallel%procmap(Fk(i,j,k)%BC%adj%b)/=parallel%myrank) &
                  NcR(parallel%procmap(Fk(i,j,k)%BC%adj%b),l) = NcR(parallel%procmap(Fk(i,j,k)%BC%adj%b),l) + 1
              endif
          enddo
        enddo
      enddo
      endassociate
    enddo
  enddo
  gNcR = sum(NcR)
  c = 0
  do l=1,mesh_dims%Nl
    do proc=0,parallel%Nproc-1
      if (NcR(proc,l)>0) then
        c             = c + NcR(proc,l)
        bbR(1,proc,l) = c - NcR(proc,l) + 1
        bbR(2,proc,l) = c
      endif
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_NcR

  !> @brief Procedure for computing the bounding boxes of sending cells of myrank to other processes.
  subroutine compute_bbS(parallel,mesh_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Parallel),        intent(IN):: parallel  !< Parallel data.
  type(Type_Mesh_Dimensions), intent(IN):: mesh_dims !< Mesh dimensions.
  integer(I_P)::                           l         !< Grid levels counter.
  integer(I_P)::                           proc      !< Processes counters.
  integer(I_P)::                           c         !< Cells counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  gNcS = sum(NcS)
  c = 0
  do l=1,mesh_dims%Nl
    do proc=0,parallel%Nproc-1
      if (NcS(proc,l)>0) then
        c             = c + NcS(proc,l)
        bbS(1,proc,l) = c - NcS(proc,l) + 1
        bbS(2,proc,l) = c
      endif
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_bbS

  !> @brief Procedure for safety allocation of send/receive variables.
  subroutine alloc_sendrecv(Np)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: Np !< Number of primitive variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(recvmap)) deallocate(recvmap) ; allocate(recvmap(1:4, 1:gNcR))        ; recvmap = 0_I_P
  if (allocated(reqsmap)) deallocate(reqsmap) ; allocate(reqsmap(1:4, 1:gNcR))        ; reqsmap = 0_I_P
  if (allocated(sendmap)) deallocate(sendmap) ; allocate(sendmap(1:4, 1:gNcS))        ; sendmap = 0_I_P
  if (allocated(Precv  )) deallocate(Precv  ) ; allocate(Precv  (1:Np,1:gNcR)) ; Precv   = 0._R_P
  if (allocated(Psend  )) deallocate(Psend  ) ; allocate(Psend  (1:Np,1:gNcS)) ; Psend   = 0._R_P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_sendrecv

  !> @brief Procedure for computing querying and receiving maps of actual process.
  subroutine compute_recv_maps(parallel,mesh_dims,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Parallel),        intent(IN):: parallel     !< Parallel data.
  type(Type_Mesh_Dimensions), intent(IN):: mesh_dims    !< Mesh dimensions.
  type(Type_SBlock),          intent(IN):: block(1:,1:) !< Block-level data.
  integer(I_P)::                           l            !< Grid levels counter.
  integer(I_P)::                           proc         !< Processes counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do l=1,mesh_dims%Nl ! loop pver grid levels
    do proc=0,parallel%Nproc-1 ! the communications is organized by processes sequence
      if (proc==parallel%myrank) cycle ! myrank doesn't communicate with itself
      if (NcR(proc,l)==0) cycle ! there are no data to communicate to process proc
      call scan_proc(NcR     = NcR(proc,l),                              &
                     proc    = proc,                                     &
                     block   = block(:,l),                               &
                     reqsmap = reqsmap(1:4,bbR(1,proc,l):bbR(2,proc,l)), &
                     recvmap = recvmap(1:4,bbR(1,proc,l):bbR(2,proc,l)))
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Procedure for searching receiving cells of myrank from process proc.
    subroutine scan_proc(NcR,proc,block,reqsmap,recvmap)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),      intent(IN)::  NcR                !< Number of send/receive cells.
    integer(I_P),      intent(IN)::  proc               !< Other process than myrank to send/receive cells.
    type(Type_SBlock), intent(IN)::  block(1:)          !< Block-level data.
    integer(I_P),      intent(OUT):: reqsmap(1:4,1:NcR) !< Querying cells map.
    integer(I_P),      intent(OUT):: recvmap(1:4,1:NcR) !< Receiving cells map.
    integer(I_P)::                   c                  !< Cells counter.
    integer(I_P)::                   b                  !< Blocks counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    c = 0_I_P ! initialize cells counter
    do b=1,mesh_dims%Nb ! loop over blocks
      call scan_block(NcR     = NcR,                &
                      block   = block(b),           &
                      proc    = proc,               &
                      b       = b,                  &
                      c       = c,                  &
                      reqsmap = reqsmap(1:4,1:NcR), &
                      recvmap = recvmap(1:4,1:NcR))
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine scan_proc

    !> @brief Procedure for searching receiving cells of myrank from process proc into the actual block "b".
    subroutine scan_block(NcR,block,proc,b,c,reqsmap,recvmap)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),      intent(IN)::    NcR                ! Number of send/receive cells.
    type(Type_SBlock), intent(IN)::    block              ! Block-level data.
    integer(I_P),      intent(IN)::    proc               ! Other process than myrank to s/r cells.
    integer(I_P),      intent(IN)::    b                  ! Actual block number.
    integer(I_P),      intent(INOUT):: c                  ! Actual cell counter.
    integer(I_P),      intent(OUT)::   reqsmap(1:4,1:NcR) ! Querying cells map.
    integer(I_P),      intent(OUT)::   recvmap(1:4,1:NcR) ! Receiving cells map.
    integer(I_P)::                     Ni,Nj,Nk,gc(1:6)   ! Temp var for storing block dimensions.
    integer(I_P)::                     i,j,k              ! Spaces counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    gc = block%dims%gc
    Ni = block%dims%Ni
    Nj = block%dims%Nj
    Nk = block%dims%Nk
    do k=1,Nk
      do j=1,Nj
        ! left i
        if (block%Fi(0,j,k)%BC%tp==bc_adj) then
          if (parallel%procmap(block%Fi(0,j,k)%BC%adj%b)==proc) then
            do i=1-gc(1),0
              c = c + 1
              reqsmap(1,c) = block%Fi(i,j,k)%BC%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%Fi(i,j,k)%BC%adj%i ; recvmap(2,c) = i
              reqsmap(3,c) = block%Fi(i,j,k)%BC%adj%j ; recvmap(3,c) = j
              reqsmap(4,c) = block%Fi(i,j,k)%BC%adj%k ; recvmap(4,c) = k
            enddo
          endif
        endif
        ! right i
        if (block%Fi(Ni,j,k)%BC%tp==bc_adj) then
          if (parallel%procmap(block%Fi(Ni,j,k)%BC%adj%b)==proc) then
            do i=Ni,Ni+gc(2)-1
              c = c + 1
              reqsmap(1,c) = block%Fi(i,j,k)%BC%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%Fi(i,j,k)%BC%adj%i ; recvmap(2,c) = i+1
              reqsmap(3,c) = block%Fi(i,j,k)%BC%adj%j ; recvmap(3,c) = j
              reqsmap(4,c) = block%Fi(i,j,k)%BC%adj%k ; recvmap(4,c) = k
            enddo
          endif
        endif
      enddo
    enddo
    do k=1,Nk
      do i=1,Ni
        ! left j
        if (block%Fj(i,0,k)%BC%tp==bc_adj) then
          if (parallel%procmap(block%Fj(i,0,k)%BC%adj%b)==proc) then
            do j=1-gc(3),0
              c = c + 1
              reqsmap(1,c) = block%Fj(i,j,k)%BC%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%Fj(i,j,k)%BC%adj%i ; recvmap(2,c) = i
              reqsmap(3,c) = block%Fj(i,j,k)%BC%adj%j ; recvmap(3,c) = j
              reqsmap(4,c) = block%Fj(i,j,k)%BC%adj%k ; recvmap(4,c) = k
            enddo
          endif
        endif
        ! right j
        if (block%Fj(i,Nj,k)%BC%tp==bc_adj) then
          if (parallel%procmap(block%Fj(i,Nj,k)%BC%adj%b)==proc) then
            do j=Nj,Nj+gc(4)-1
              c = c + 1
              reqsmap(1,c) = block%Fj(i,j,k)%BC%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%Fj(i,j,k)%BC%adj%i ; recvmap(2,c) = i
              reqsmap(3,c) = block%Fj(i,j,k)%BC%adj%j ; recvmap(3,c) = j+1
              reqsmap(4,c) = block%Fj(i,j,k)%BC%adj%k ; recvmap(4,c) = k
            enddo
          endif
        endif
      enddo
    enddo
    do j=1,Nj
      do i=1,Ni
        ! left k
        if (block%Fk(i,j,0)%BC%tp==bc_adj) then
          if (parallel%procmap(block%Fk(i,j,0)%BC%adj%b)==proc) then
            do k=1-gc(5),0
              c = c + 1
              reqsmap(1,c) = block%Fk(i,j,k)%BC%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%Fk(i,j,k)%BC%adj%i ; recvmap(2,c) = i
              reqsmap(3,c) = block%Fk(i,j,k)%BC%adj%j ; recvmap(3,c) = j
              reqsmap(4,c) = block%Fk(i,j,k)%BC%adj%k ; recvmap(4,c) = k
            enddo
          endif
        endif
        ! right k
        if (block%Fk(i,j,Nk)%BC%tp==bc_adj) then
          if (parallel%procmap(block%Fk(i,j,Nk)%BC%adj%b)==proc) then
            do k=Nk,Nk+gc(6)-1
              c = c + 1
              reqsmap(1,c) = block%Fk(i,j,k)%BC%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%Fk(i,j,k)%BC%adj%i ; recvmap(2,c) = i
              reqsmap(3,c) = block%Fk(i,j,k)%BC%adj%j ; recvmap(3,c) = j
              reqsmap(4,c) = block%Fk(i,j,k)%BC%adj%k ; recvmap(4,c) = k+1
            enddo
          endif
        endif
      enddo
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine scan_block
  endsubroutine compute_recv_maps

  !> @brief Procedure for communicate the number of cells that myrank must receive from other processes.
  subroutine NcRsendrecv(parallel,mesh_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Parallel),        intent(IN):: parallel                   !< Parallel data.
  type(Type_Mesh_Dimensions), intent(IN):: mesh_dims                  !< Mesh dimensions.
  integer(I_P)::                           l                          !< Grid levels counter.
  integer(I_P)::                           proc                       !< Processes counters.
  integer(I_P)::                           ierr,stat(MPI_STATUS_SIZE) !< MPI error flags.
  integer(I_P), parameter::                tagshift=0*maxproc         !< Shift for tags (to isolate these kind of communications).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do l=1,mesh_dims%Nl ! loop over grid levels
    do proc=0,parallel%Nproc-1 ! the communications is organized by processes sequence
      ! sending querying (receiving) cells number of myrank to proc and
      ! using the querying number cells of proc for building the sending number of cells of myrank
      call MPI_SENDRECV(NcR(proc,l),1,MPI_INTEGER,proc,tagshift+parallel%Nproc*(parallel%myrank+1), &
                        NcS(proc,l),1,MPI_INTEGER,proc,tagshift+parallel%Nproc*(proc           +1), &
                        MPI_COMM_WORLD,stat,ierr)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine NcRsendrecv

  !> @brief Procedure for communicate the querying cells map of myrank to other processes and for building the sending cells maps
  !> of myrank for other processes (the building is done by receiving the querying cells map of other processes).
  subroutine mapsendrecv(parallel,mesh_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Parallel),        intent(IN):: parallel                   !< Parallel data.
  type(Type_Mesh_Dimensions), intent(IN):: mesh_dims                  !< Mesh dimensions.
  integer(I_P)::                           l                          !< Grid levels counter.
  integer(I_P)::                           proc                       !< Processes counter.
  integer(I_P)::                           ierr,stat(MPI_STATUS_SIZE) !< MPI error flags.
  integer(I_P), parameter::                tagshift=1*maxproc         !< Shift for tags (to isolate these kind of communications).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do l=1,mesh_dims%Nl ! loop over grid levels
    do proc=0,parallel%Nproc-1 ! the communications is organized by processes sequence
      if ((NcR(proc,l)==0).AND.(NcS(proc,l)==0)) cycle ! there are no data to communicate to process proc
      ! sending querying (reqsmap) cells map of myrank to proc and
      ! using the querying cells map of proc for building the sending (sendmap) cells map of myrank
      call MPI_SENDRECV(reqsmap(1,bbR(1,proc,l)),4*NcR(proc,l),MPI_INTEGER,proc,tagshift+parallel%Nproc*(parallel%myrank+1), &
                        sendmap(1,bbS(1,proc,l)),4*NcS(proc,l),MPI_INTEGER,proc,tagshift+parallel%Nproc*(proc           +1), &
                        MPI_COMM_WORLD,stat,ierr)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine mapsendrecv
  !> @}

  !> @ingroup Lib_ParallelPublicProcedure
  !> @{
  !> @brief Procedure for initializing MPI send/receive maps.
  subroutine init_MPI_maps(parallel,mesh_dims,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Parallel),        intent(IN):: parallel     !< Parallel data.
  type(Type_Mesh_Dimensions), intent(IN):: mesh_dims    !< Mesh dimensions.
  type(Type_SBlock),          intent(IN):: block(1:,1:) !< Block-level data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call alloc_SR(parallel=parallel,mesh_dims=mesh_dims)
  call compute_NcR(parallel=parallel,mesh_dims=mesh_dims,block=block)
  call NcRsendrecv(parallel=parallel,mesh_dims=mesh_dims)
  call compute_bbS(parallel=parallel,mesh_dims=mesh_dims)
  call alloc_sendrecv(Np=block(1,1)%dims%Np)
  call compute_recv_maps(parallel=parallel,mesh_dims=mesh_dims,block=block)
  call mapsendrecv(parallel=parallel,mesh_dims=mesh_dims)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_MPI_maps

  !> @brief Procedure for pretty printing MPI send/receive maps.
  subroutine print_MPI_maps(parallel,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Parallel),    intent(IN)::  parallel !< Parallel data.
  character(*), optional, intent(IN)::  pref     !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat   !< IO error.
  character(*), optional, intent(OUT):: iomsg    !< IO error message.
  integer(I4P),           intent(IN)::  unit     !< Logic unit.
  character(len=:), allocatable::       prefd    !< Prefixing string.
  integer(I4P)::                        iostatd  !< IO error.
  character(500)::                      iomsgd   !< Temporary variable for IO error message.
  integer(I4P)::                        proc     !< Processes counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  do proc=0,parallel%Nproc-1 ! the communications is organized by processes sequence
    if (proc==parallel%myrank) cycle ! myrank doesn't communicate with itself
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//&
      ' Process '//trim(str(.true.,parallel%myrank))//           &
      ' must send '//trim(str(.true.,NcS(proc,1)))//             &
      ' and receive '//trim(str(.true.,NcR(proc,1)))//           &
      ' finite voluems to/from process '//trim(str(.true.,proc))
  enddo
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_MPI_maps

  !> Subroutine for performing send/receive of primitive variables of myrank to other processes.
  subroutine prim_sendrecv(l,parallel,mesh_dims,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),               intent(IN)::    l                          !< Grid level.
  type(Type_Parallel),        intent(IN)::    parallel                   !< Parallel data.
  type(Type_Mesh_Dimensions), intent(IN)::    mesh_dims                  !< Mesh dimensions.
  type(Type_SBlock),          intent(INOUT):: block(1:)                  !< Block-level data.
  integer(I_P)::                              proc                       !< Processes counter.
  integer(I_P)::                              c                          !< Cells counter.
  integer(I_P)::                              b                          !< Blocks counter.
  integer(I_P)::                              ierr,stat(MPI_STATUS_SIZE) !< MPI error flags.
  integer(I_P), parameter::                   tagshift=2*maxproc         !< Shift for tags (to isolate these communications).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do proc=0,parallel%Nproc-1 ! the communications is organized by processes sequence
    if ((NcR(proc,l)==0).AND.(NcR(proc,l)==0)) cycle ! there are no data to communicate to process proc
    ! building the send buffer (Psend) of myrank for proc using local var P
    do c=bbS(1,proc,l),bbS(2,proc,l)
      b = minloc(array=parallel%blockmap,dim=1,mask=parallel%blockmap==sendmap(1,c))
      Psend(:,c) = block(b)%C(sendmap(2,c),sendmap(3,c),sendmap(4,c))%P%prim2array()
    enddo
    ! sending Psend of myrank to proc and storing in Precv of proc
    associate(Np=>block(1)%dims%Np,Nproc=>parallel%Nproc,myrank=>parallel%myrank)
      call MPI_SENDRECV(Psend(1,bbS(1,proc,l)),Np*NcS(proc,l),MPI_REAL8,proc,tagshift+Nproc*(myrank+1), &
                        Precv(1,bbR(1,proc,l)),Np*NcR(proc,l),MPI_REAL8,proc,tagshift+Nproc*(proc  +1), &
                        MPI_COMM_WORLD,stat,ierr)
    endassociate
    ! coping the receive buffer (Precv) of myrank from proc in the local var P
    do c=bbR(1,proc,l),bbR(2,proc,l)
       call block(recvmap(1,c))%C(recvmap(2,c),recvmap(3,c),recvmap(4,c))%P%array2prim(Precv(:,c))
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prim_sendrecv
  !> @}
#endif
endmodule Lib_Parallel
