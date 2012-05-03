!> This module contains the definition of procedures for send/receive data among processes for parallel (MPI) operations.
!> It is based on MPI library.
!> @note The communications have a tag-shift (for make them unique) that assumes a maximum number of processes of 10000.
!> Increment this parameter if using more processes than 10000.
!> @todo \b DocComplete: Complete the documentation of internal procedures
module Lib_Parallel
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision        ! Integers and reals precision definition.
USE Data_Type_BC        ! Definition of Type_BC.
USE Data_Type_Globals   ! Definition of Type_Global and Type_Block.
USE Data_Type_Primitive ! Definition of Type_Primitive.
USE Lib_IO_Misc         ! Procedures for IO and strings operations.
#ifdef MPI2
USE MPI                 ! MPI runtime library.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: Nthreads
public:: Nproc
public:: procmap
public:: blockmap
public:: Init_sendrecv
public:: Psendrecv
public:: procmap_load
public:: procmap_save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
integer(I_P)::              Nthreads = 1_I_P !< Number of OpenMP threads.
integer(I_P)::              Nproc    = 1_I_P !< Number of processes (for MPI parallelization).
integer(I_P), allocatable:: procmap(:)       !< Processes/blocks map    [1:Nb_tot].
integer(I_P), allocatable:: blockmap(:)      !< Local/global blocks map [1:Nb].
integer(I_P), parameter::   maxproc = 10000  !< Maximum number of processes used for communications tag shift.
integer(I_P)             :: gNcR             !< Global number of receive cells (sum(NcR)).
integer(I_P)             :: gNcS             !< Global number of send   cells (sum(NcS)).
integer(I_P), allocatable:: NcR(:,:)         !< Number of receive cells from each process [    0:Nproc-1,1:Nl].
integer(I_P), allocatable:: NcS(:,:)         !< Number of send   cells for  each process  [    0:Nproc-1,1:Nl].
integer(I_P), allocatable:: bbR(:,:,:)       !< Processes bounds of receive cells         [1:2,0:Nproc-1,1:Nl].
integer(I_P), allocatable:: bbS(:,:,:)       !< Processes bounds of send   cells          [1:2,0:Nproc-1,1:Nl].
integer(I_P), allocatable:: recvmap(:,:)     !< Receiving cells map of   myrank from other processes [1:4,1:gNcR].
integer(I_P), allocatable:: reqsmap(:,:)     !< Querying  cells map of   myrank for  other processes [1:4,1:gNcR].
integer(I_P), allocatable:: sendmap(:,:)     !< Sending  cells map from myrank to   other processes [1:4,1:gNcS].
real(R_P),    allocatable:: Precv(:,:)       !< Receiving buffer of primitive variable of myrank from other processes [1:Np,1:gNcR].
real(R_P),    allocatable:: Psend(:,:)       !< Sending  buffer of primitive variable of myrank for other  processes [1:Np,1:gNcS].
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine alloc_SR(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for safety allocation of NcR, NcS, bbR and bbS.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(IN):: global ! Global-level data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(NcR)) deallocate(NcR) ; allocate(NcR(    0:Nproc,1:global%mesh%Nl)) ; NcR = 0_I_P
  if (allocated(NcS)) deallocate(NcS) ; allocate(NcS(    0:Nproc,1:global%mesh%Nl)) ; NcS = 0_I_P
  if (allocated(bbR)) deallocate(bbR) ; allocate(bbR(1:2,0:Nproc,1:global%mesh%Nl)) ; bbR = 0_I_P
  if (allocated(bbS)) deallocate(bbS) ; allocate(bbS(1:2,0:Nproc,1:global%mesh%Nl)) ; bbS = 0_I_P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_SR

  subroutine compute_NcR(myrank,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for computing the number of cells that must be received from other processes than myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN):: myrank                                   ! Actual rank process.
  type(Type_Global), intent(IN):: global                                   ! Global-level data.
  type(Type_Block),  intent(IN):: block(1:global%mesh%Nb,1:global%mesh%Nl) ! Block-level data.
  integer(I_P)::                  l                                        ! Grid levels counter.
  integer(I_P)::                  proc                                     ! Processes counter.
  integer(I_P)::                  c                                        ! Cells counter.
  integer(I_P)::                  b                                        ! Blocks counter.
  integer(I_P)::                  i,j,k                                    ! Space counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing the number of receive cells of myrank for other processes; checking the equivalence of proc and myrank:
  !  if proc=myrank => NcR=0 => myrank doesn't communicate with itself
  do l=1,global%mesh%Nl
    do b=1,global%mesh%Nb
      ! i interfaces
      do k=1,block(b,l)%mesh%Nk
        do j=1,block(b,l)%mesh%Nj
          ! left i
          do i=0-block(b,l)%mesh%gc(1),0
            if (block(b,l)%bc%BCi(i,j,k)%tp==bc_adj) then
              if (procmap(block(b,l)%bc%BCi(i,j,k)%adj%b)/=myrank) &
              NcR(procmap(block(b,l)%bc%BCi(i,j,k)%adj%b),l) = NcR(procmap(block(b,l)%bc%BCi(i,j,k)%adj%b),l) + 1
            endif
          enddo
          ! right i
          do i=block(b,l)%mesh%Ni,block(b,l)%mesh%Ni+block(b,l)%mesh%gc(2)
            if (block(b,l)%bc%BCi(i,j,k)%tp==bc_adj) then
              if (procmap(block(b,l)%bc%BCi(i,j,k)%adj%b)/=myrank) &
              NcR(procmap(block(b,l)%bc%BCi(i,j,k)%adj%b),l) = NcR(procmap(block(b,l)%bc%BCi(i,j,k)%adj%b),l) + 1
            endif
          enddo
        enddo
      enddo
      ! j interfaces
      do k=1,block(b,l)%mesh%Nk
        ! left j
        do j=0-block(b,l)%mesh%gc(3),0
          do i=1,block(b,l)%mesh%Ni
            if (block(b,l)%bc%BCj(i,j,k)%tp==bc_adj) then
              if (procmap(block(b,l)%bc%BCj(i,j,k)%adj%b)/=myrank) &
              NcR(procmap(block(b,l)%bc%BCj(i,j,k)%adj%b),l) = NcR(procmap(block(b,l)%bc%BCj(i,j,k)%adj%b),l) + 1
            endif
          enddo
        enddo
        ! right j
        do j=block(b,l)%mesh%Nj,block(b,l)%mesh%Nj+block(b,l)%mesh%gc(4)
          do i=1,block(b,l)%mesh%Ni
            if (block(b,l)%bc%BCj(i,j,k)%tp==bc_adj) then
              if (procmap(block(b,l)%bc%BCj(i,j,k)%adj%b)/=myrank) &
              NcR(procmap(block(b,l)%bc%BCj(i,j,k)%adj%b),l) = NcR(procmap(block(b,l)%bc%BCj(i,j,k)%adj%b),l) + 1
            endif
          enddo
        enddo
      enddo
      ! k interfaces
      ! left k
      do k=0-block(b,l)%mesh%gc(5),0
        do j=1,block(b,l)%mesh%Nj
          do i=1,block(b,l)%mesh%Ni
            if (block(b,l)%bc%BCk(i,j,k)%tp==bc_adj) then
              if (procmap(block(b,l)%bc%BCk(i,j,k)%adj%b)/=myrank) &
              NcR(procmap(block(b,l)%bc%BCk(i,j,k)%adj%b),l) = NcR(procmap(block(b,l)%bc%BCk(i,j,k)%adj%b),l) + 1
            endif
          enddo
        enddo
      enddo
      ! right k
      do k=block(b,l)%mesh%Nk,block(b,l)%mesh%Nk+block(b,l)%mesh%gc(6)
        do j=1,block(b,l)%mesh%Nj
          do i=1,block(b,l)%mesh%Ni
            if (block(b,l)%bc%BCk(i,j,k)%tp==bc_adj) then
              if (procmap(block(b,l)%bc%BCk(i,j,k)%adj%b)/=myrank) &
              NcR(procmap(block(b,l)%bc%BCk(i,j,k)%adj%b),l) = NcR(procmap(block(b,l)%bc%BCk(i,j,k)%adj%b),l) + 1
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  gNcR = sum(NcR)
  c = 0
  do l=1,global%mesh%Nl
    do proc=0,Nproc-1
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

  subroutine NcRsendrecv(myrank,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for communicate the number of cells that myrank must receive from other processes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN):: myrank                     ! Actual rank process.
  type(Type_Global), intent(IN):: global                     ! Global-level data.
  integer(I_P)::                  l                          ! Grid levels counter.
  integer(I_P)::                  proc                       ! Processes counters.
#ifdef MPI2
  integer(I_P)::                  ierr,stat(MPI_STATUS_SIZE) ! MPI error flags.
  integer(I_P), parameter::       tagshift=0*maxproc         ! Shift for tags (to isolate these kind of communications).
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do l=1,global%mesh%Nl ! loop over grid levels
    do proc=0,Nproc-1 ! the communications is organized by processes sequence
#ifdef MPI2
      ! sending querying (receiving) cells number of myrank to proc and
      ! using the querying number cells of proc for building the sending number of cells of myrank
      call MPI_SENDRECV(NcR(proc,l),1,MPI_INTEGER,proc,tagshift+Nproc*(myrank+1), &
                        NcS(proc,l),1,MPI_INTEGER,proc,tagshift+Nproc*(proc  +1), &
                        MPI_COMM_WORLD,stat,ierr)
#endif
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine NcRsendrecv

  subroutine compute_bbS(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for computing the bounding boxes of sending cells of myrank to other processes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(IN):: global ! Global-level data.
  integer(I_P)::                  l      ! Grid levels counter.
  integer(I_P)::                  proc   ! Processes counters.
  integer(I_P)::                  c      ! Cells counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  gNcS = sum(NcS)
  c = 0
  do l=1,global%mesh%Nl
    do proc=0,Nproc-1
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

  subroutine alloc_sendrecv(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for safety allocation of send/receive variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(IN):: global ! Global-level data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(recvmap)) deallocate(recvmap) ; allocate(recvmap(1:4, 1:gNcR))              ; recvmap = 0_I_P
  if (allocated(reqsmap)) deallocate(reqsmap) ; allocate(reqsmap(1:4, 1:gNcR))              ; reqsmap = 0_I_P
  if (allocated(sendmap)) deallocate(sendmap) ; allocate(sendmap(1:4, 1:gNcS))              ; sendmap = 0_I_P
  if (allocated(Precv  )) deallocate(Precv  ) ; allocate(Precv  (1:global%fluid%Np,1:gNcR)) ; Precv   = 0._R_P
  if (allocated(Psend  )) deallocate(Psend  ) ; allocate(Psend  (1:global%fluid%Np,1:gNcS)) ; Psend   = 0._R_P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_sendrecv

  subroutine compute_recv_maps(myrank,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing querying and receiving maps of actual process.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN):: myrank                                   ! Actual rank process.
  type(Type_Global), intent(IN):: global                                   ! Global-level data.
  type(Type_Block),  intent(IN):: block(1:global%mesh%Nb,1:global%mesh%Nl) ! Block-level data.
  integer(I_P)::                  l                                        ! Grid levels counter.
  integer(I_P)::                  proc                                     ! Processes counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do l=1,global%mesh%Nl ! loop pver grid levels
    do proc=0,Nproc-1 ! the communications is organized by processes sequence
      if (proc==myrank) cycle ! myrank doesn't communicate with itself
      if (NcR(proc,l)==0) cycle ! there are no data to communicate to process proc
      call scan_proc(NcR     = NcR(proc,l),                              &
                     proc    = proc,                                     &
                     global  = global,                                   &
                     block   = block(:,l),                               &
                     reqsmap = reqsmap(1:4,bbR(1,proc,l):bbR(2,proc,l)), &
                     recvmap = recvmap(1:4,bbR(1,proc,l):bbR(2,proc,l)))
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    subroutine scan_proc(NcR,proc,global,block,reqsmap,recvmap)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for searching receiving cells of myrank from process proc.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),      intent(IN)::  NcR                     ! Number of send/receive cells.
    integer(I_P),      intent(IN)::  proc                    ! Other process than myrank to send/receive cells.
    type(Type_Global), intent(IN)::  global                  ! Global-level data.
    type(Type_Block),  intent(IN)::  block(1:global%mesh%Nb) ! Block-level data.
    integer(I_P),      intent(OUT):: reqsmap(1:4,1:NcR)      ! Querying cells map.
    integer(I_P),      intent(OUT):: recvmap(1:4,1:NcR)      ! Receiving cells map.
    integer(I_P)::                   c                       ! Cells counter.
    integer(I_P)::                   b                       ! Blocks counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    c = 0_I_P ! initialize cells counter
    do b=1,global%mesh%Nb ! loop over blocks
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

    subroutine scan_block(NcR,block,proc,b,c,reqsmap,recvmap)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for searching receiving cells of myrank from process proc into the actual block "b".
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),      intent(IN)::    NcR                ! Number of send/receive cells.
    type(Type_Block),  intent(IN)::    block              ! Block-level data.
    integer(I_P),      intent(IN)::    proc               ! Other process than myrank to s/r cells.
    integer(I_P),      intent(IN)::    b                  ! Actual block number.
    integer(I_P),      intent(INOUT):: c                  ! Actual cell counter.
    integer(I_P),      intent(OUT)::   reqsmap(1:4,1:NcR) ! Querying cells map.
    integer(I_P),      intent(OUT)::   recvmap(1:4,1:NcR) ! Receiving cells map.
    integer(I_P)::                     Ni,Nj,Nk,gc(1:6)   ! Temp var for storing block dimensions.
    integer(I_P)::                     i,j,k              ! Spaces counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    Ni = block%mesh%Ni
    Nj = block%mesh%Nj
    Nk = block%mesh%Nk
    gc = block%mesh%gc
    do k=1,Nk
      do j=1,Nj
        ! left i
        if (block%bc%BCi(0,j,k)%tp==bc_adj) then
          if (procmap(block%bc%BCi(0,j,k)%adj%b)==proc) then
            do i=1-gc(1),0
              c = c + 1
              reqsmap(1,c) = block%bc%BCi(i,j,k)%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%bc%BCi(i,j,k)%adj%i ; recvmap(2,c) = i
              reqsmap(3,c) = block%bc%BCi(i,j,k)%adj%j ; recvmap(3,c) = j
              reqsmap(4,c) = block%bc%BCi(i,j,k)%adj%k ; recvmap(4,c) = k
            enddo
          endif
        endif
        ! right i
        if (block%bc%BCi(Ni,j,k)%tp==bc_adj) then
          if (procmap(block%bc%BCi(Ni,j,k)%adj%b)==proc) then
            do i=Ni,Ni+gc(2)-1
              c = c + 1
              reqsmap(1,c) = block%bc%BCi(i,j,k)%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%bc%BCi(i,j,k)%adj%i ; recvmap(2,c) = i+1
              reqsmap(3,c) = block%bc%BCi(i,j,k)%adj%j ; recvmap(3,c) = j
              reqsmap(4,c) = block%bc%BCi(i,j,k)%adj%k ; recvmap(4,c) = k
            enddo
          endif
        endif
      enddo
    enddo
    do k=1,Nk
      do i=1,Ni
        ! left j
        if (block%bc%BCj(i,0,k)%tp==bc_adj) then
          if (procmap(block%bc%BCj(i,0,k)%adj%b)==proc) then
            do j=1-gc(3),0
              c = c + 1
              reqsmap(1,c) = block%bc%BCj(i,j,k)%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%bc%BCj(i,j,k)%adj%i ; recvmap(2,c) = i
              reqsmap(3,c) = block%bc%BCj(i,j,k)%adj%j ; recvmap(3,c) = j
              reqsmap(4,c) = block%bc%BCj(i,j,k)%adj%k ; recvmap(4,c) = k
            enddo
          endif
        endif
        ! right j
        if (block%bc%BCj(i,Nj,k)%tp==bc_adj) then
          if (procmap(block%bc%BCj(i,Nj,k)%adj%b)==proc) then
            do j=Nj,Nj+gc(4)-1
              c = c + 1
              reqsmap(1,c) = block%bc%BCj(i,j,k)%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%bc%BCj(i,j,k)%adj%i ; recvmap(2,c) = i
              reqsmap(3,c) = block%bc%BCj(i,j,k)%adj%j ; recvmap(3,c) = j+1
              reqsmap(4,c) = block%bc%BCj(i,j,k)%adj%k ; recvmap(4,c) = k
            enddo
          endif
        endif
      enddo
    enddo
    do j=1,Nj
      do i=1,Ni
        ! left k
        if (block%bc%BCk(i,j,0)%tp==bc_adj) then
          if (procmap(block%bc%BCk(i,j,0)%adj%b)==proc) then
            do k=1-gc(5),0
              c = c + 1
              reqsmap(1,c) = block%bc%BCk(i,j,k)%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%bc%BCk(i,j,k)%adj%i ; recvmap(2,c) = i
              reqsmap(3,c) = block%bc%BCk(i,j,k)%adj%j ; recvmap(3,c) = j
              reqsmap(4,c) = block%bc%BCk(i,j,k)%adj%k ; recvmap(4,c) = k
            enddo
          endif
        endif
        ! right k
        if (block%bc%BCk(i,j,Nk)%tp==bc_adj) then
          if (procmap(block%bc%BCk(i,j,Nk)%adj%b)==proc) then
            do k=Nk,Nk+gc(6)-1
              c = c + 1
              reqsmap(1,c) = block%bc%BCk(i,j,k)%adj%b ; recvmap(1,c) = b
              reqsmap(2,c) = block%bc%BCk(i,j,k)%adj%i ; recvmap(2,c) = i
              reqsmap(3,c) = block%bc%BCk(i,j,k)%adj%j ; recvmap(3,c) = j
              reqsmap(4,c) = block%bc%BCk(i,j,k)%adj%k ; recvmap(4,c) = k+1
            enddo
          endif
        endif
      enddo
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine scan_block
  endsubroutine compute_recv_maps

  subroutine mapsendrecv(myrank,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for communicate the querying cells map of myrank to other processes and for building the sending cells maps of myrank
  ! for other processes (the building is done by receiving the querying cells map of other processes).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN):: myrank                     ! Actual rank process.
  type(Type_Global), intent(IN):: global                     ! Global-level data.
  integer(I_P)::                  l                          ! Grid levels counter.
  integer(I_P)::                  proc                       ! Processes counter.
#ifdef MPI2
  integer(I_P)::                  ierr,stat(MPI_STATUS_SIZE) ! MPI error flags.
  integer(I_P), parameter::       tagshift=1*maxproc         ! Shift for tags (to isolate these kind of communications).
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do l=1,global%mesh%Nl ! loop over grid levels
    do proc=0,Nproc-1 ! the communications is organized by processes sequence
      if ((NcR(proc,l)==0).AND.(NcS(proc,l)==0)) cycle ! there are no data to communicate to process proc
#ifdef MPI2
      ! sending querying (reqsmap) cells map of myrank to proc and
      ! using the querying cells map of proc for building the sending (sendmap) cells map of myrank
      call MPI_SENDRECV(reqsmap(1,bbR(1,proc,l)),4*NcR(proc,l),MPI_INTEGER,proc,tagshift+Nproc*(myrank+1), &
                        sendmap(1,bbS(1,proc,l)),4*NcS(proc,l),MPI_INTEGER,proc,tagshift+Nproc*(proc  +1), &
                        MPI_COMM_WORLD,stat,ierr)
#endif
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine mapsendrecv

  !> Subroutine for initializing the send/receive communications data.
  subroutine Init_sendrecv(myrank,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN):: myrank                                   !< Actual rank process.
  type(Type_Global), intent(IN):: global                                   !< Global-level data.
  type(Type_Block),  intent(IN):: block(1:global%mesh%Nb,1:global%mesh%Nl) !< Block-level data.
  integer(I_P)::                  err                                      !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call alloc_SR(global=global)
  call compute_NcR(myrank=myrank,global=global,block=block)
  call NcRsendrecv(myrank=myrank,global=global)
  call compute_bbS(global=global)
  call alloc_sendrecv(global=global)
  call compute_recv_maps(myrank=myrank,global=global,block=block)
  call mapsendrecv(myrank=myrank,global=global)
  err = Printsendrecvmaps(myrank)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_sendrecv

  !> Subroutine for performing send/receive of primitive variables of myrank to other processes.
  subroutine Psendrecv(myrank,l,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    myrank                     !< Actual rank process.
  integer(I_P),      intent(IN)::    l                          !< Grid level.
  type(Type_Global), intent(IN)::    global                     !< Global-level data.
  type(Type_Block),  intent(INOUT):: block(1:global%mesh%Nb)    !< Block-level data.
  integer(I_P)::                     proc                       !< Processes counter.
  integer(I_P)::                     c                          !< Cells counter.
  integer(I_P)::                     b                          !< Blocks counter.
#ifdef MPI2
  integer(I_P)::                     ierr,stat(MPI_STATUS_SIZE) !< MPI error flags.
  integer(I_P), parameter::          tagshift=2*maxproc         !< Shift for tags (to isolate these communications).
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do proc=0,Nproc-1 ! the communications is organized by processes sequence
    if ((NcR(proc,l)==0).AND.(NcR(proc,l)==0)) cycle ! there are no data to communicate to process proc
    ! building the send buffer (Psend) of myrank for proc using local var P
    do c=bbS(1,proc,l),bbS(2,proc,l)
      b = minloc(array=blockmap,dim=1,mask=blockmap==sendmap(1,c))
      Psend(:,c) = prim2array(block(b)%fluid%P(sendmap(2,c),sendmap(3,c),sendmap(4,c)))
    enddo
#ifdef MPI2
    ! sending Psend of myrank to proc and storing in Precv of proc
    call MPI_SENDRECV(Psend(1,bbS(1,proc,l)),global%fluid%Np*NcS(proc,l),MPI_REAL8,proc,tagshift+Nproc*(myrank+1), &
                      Precv(1,bbR(1,proc,l)),global%fluid%Np*NcR(proc,l),MPI_REAL8,proc,tagshift+Nproc*(proc  +1), &
                      MPI_COMM_WORLD,stat,ierr)
#endif
    ! coping the receive buffer (Precv) of myrank from proc in the local var P
    do c=bbR(1,proc,l),bbR(2,proc,l)
       block(recvmap(1,c))%fluid%P(recvmap(2,c),recvmap(3,c),recvmap(4,c)) = array2prim(Precv(:,c))
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Psendrecv

  function Printsendrecvmaps(myrank) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for send/receive primitive variables of myrank to other processes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: myrank ! Actual rank process.
  integer(I_P)::             err    ! Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::             proc   ! Processes counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(stdout,'(A)',iostat=err)'rank'//trim(str(.true.,myrank))//&
                                '----------------------------------------------------------------------'
  do proc=0,Nproc-1 ! the communications is organized by processes sequence
    if (proc==myrank) cycle ! myrank doesn't communicate with itself
    write(stdout,'(A)',IOSTAT=err) 'rank'//trim(str(.true.,myrank))//' Process '//trim(str(.true.,myrank))// &
                                   ' must send '//trim(str(.true.,NcS(proc,1)))//  &
                                   ' and receive '//trim(str(.true.,NcR(proc,1)))// &
                                   ' finite voluems with process '//trim(str(.true.,proc))
  enddo
  write(stdout,'(A)',iostat=err)'rank'//trim(str(.true.,myrank))//&
                                '----------------------------------------------------------------------'
  write(stdout,*)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Printsendrecvmaps

  !> Function for loading the processes/blocks map and local/global blocks map.
  !> @return \b err integer(I4P) variable.
  function procmap_load(myrank,filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< File name processes/blocks map.
  type(Type_Global), intent(INOUT):: global   !< Global-level data.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of procmap file.
  integer(I_P)::                     b,bb     !< Blocks counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=trim(filename),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'procmap_load')
  UnitFree = Get_Unit() ; open(unit = UnitFree, file = trim(filename), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*,iostat=err) global%mesh%Nb_tot
  if (allocated(procmap)) deallocate(procmap) ; allocate(procmap(1:global%mesh%Nb_tot)) ; procmap  = 0_I_P
  read(UnitFree,*,iostat=err)
  do b=1,global%mesh%Nb_tot
    read(UnitFree,*,iostat=err) procmap(b) ! reading the process number of bth block
  enddo
  close(UnitFree)
  ! computing the local/global blocks map
  if (Nproc==1_I_P) then
    ! there is no MPI environment thus all blocks are loaded by process 0
    global%mesh%Nb = global%mesh%Nb_tot
    if (allocated(blockmap)) deallocate(blockmap) ; allocate(blockmap(1:global%mesh%Nb)) ; blockmap = 0_I_P
    do b=1,global%mesh%Nb_tot
      blockmap(b) = b  ! the blocks map is identity
    enddo
  else
    ! computing the local (of myrank) number of blocks
    global%mesh%Nb = 0
    do b=1,global%mesh%Nb_tot
      if (procmap(b)==myrank) global%mesh%Nb = global%mesh%Nb + 1
    enddo
    if (allocated(blockmap)) deallocate(blockmap) ; allocate(blockmap(1:global%mesh%Nb )) ; blockmap = 0_I_P
    bb = 0
    do b=1,global%mesh%Nb_tot
      if (procmap(b)==myrank) then
        bb = bb + 1
        blockmap(bb) = b
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction procmap_load

  !> Function for saving the processes/blocks map.
  !> @return \b err integer(I4P) variable.
  function procmap_save(filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),      intent(IN):: filename !< File name processes/blocks map.
  type(Type_Global), intent(IN):: global   !< Global-level data.
  integer(I_P)::                  err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                  UnitFree !< Free logic unit.
  integer(I_P)::                  b        !< Block counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! saving the map
  UnitFree = Get_Unit()
  open(unit = UnitFree, file = trim(filename), form = 'FORMATTED')
  write(UnitFree,'(A)',iostat=err) str(.true.,global%mesh%Nb_tot)//' Nb_tot = number of total blocks'
  write(UnitFree,*,iostat=err)
  do b=1,global%mesh%Nb_tot
    write(UnitFree,'(A)',iostat=err) str(.true.,procmap(b))//' proc(bth) = rank of process where bth block is loaded'
  enddo
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction procmap_save
endmodule Lib_Parallel
