!< OFF, parallel object definition.

module off_parallel_object
!< OFF, parallel object definition.

implicit none
public :: parallel_object

type :: parallel_object
   !< OFF, parallel object class.
   integer(I4P)              :: myrank=0                     !< Myrank.
   integer(I4P)              :: threads_number=1             !< Number of OpenMP threads.
   integer(I4P)              :: procs_number=1               !< Number of processes (for MPI parallelization).
   integer(I4P), allocatable :: block_to_process_map(:)      !< Maps of processes ID for all blocks, [1:blocks_global_number].
   integer(I4P), allocatable :: block_local_to_global_map(:) !< Maps of local block ID to global,    [1:blocks_number].
   integer(I4P)              :: recv_cells_number_global=0   !< Global number of receive cells (sum(recv_cells_number)).
   integer(I4P)              :: send_cells_number_global=0   !< Global number of send    cells (sum(send_cells_number)).

   integer(I4P), allocatable :: recv_cells_number(:)         !< Number of receive cells from each process [0:procs_number-1].
   integer(I4P), allocatable :: send_cells_number(:)         !< Number of send    cells from each process [0:procs_number-1].
   integer(I4P), allocatable :: recv_bb(:,:)                 !< Processes bounds of receive cells [1:2,0:procs_number-1].
   integer(I4P), allocatable :: send_bb(:,:)                 !< Processes bounds of send   cells  [1:2,0:procs_number-1].

   ! integer(I4P), allocatable :: recvmap(:,:)     !< Receiving cells map of  myrank from other processes [1:4,1:gNcR].
   ! integer(I4P), allocatable :: reqsmap(:,:)     !< Querying  cells map of  myrank for  other processes [1:4,1:gNcR].
   ! integer(I4P), allocatable :: sendmap(:,:)     !< Sending  cells map from myrank to   other processes [1:4,1:gNcS].

   real(R8P), allocatable :: field_recv(:) !< Receiving buffer field of myrank from all other processes[1:recv_cells_number_global].
   real(R8P), allocatable :: field_send(:) !< Sending   buffer field of myrank for  all other processes[1:send_cells_number_global].

endtype parallel_object

contains
   pure subroutine initialize(self, procs_number, myrank)
   !< Initialize parallel.
   class(parallel_object), intent(inout) :: self         !< Parallel.
   integer(I4P),           intent(in)    :: procs_number !< Number of processes.
   integer(I4P),           intent(in)    :: myrank       !< Myrank.

   self%procs_number = procs_number
   self%myrank = myrank
   endsubroutine initialize

   pure subroutine alloc(self)
   !< Allocate maps.
   class(parallel_object), intent(inout) :: self !< Parallel.

   associate(procs_number=>self%procs_number,                                                       &
             recv_cells_number=>self%recv_cells_number, send_cells_number=>send%recv_cells_number,  &
             recv_bb=>self%recv_bb, send_bb=>send%recv_bb)
      if (allocated(recv_cells_number)) deallocate(recv_cells_number)
      allocate(recv_cells_number(0:procs_number-1)) ; recv_cells_number = 0
      if (allocated(send_cells_number)) deallocate(send_cells_number)
      allocate(send_cells_number(0:procs_number-1)) ; send_cells_number = 0
      if (allocated(recv_bb)) deallocate(recv_bb)
      allocate(recv_bb(1:2, 0:procs_number-1)) ; recv_bb = 0
      if (allocated(send_bb)) deallocate(send_bb)
      allocate(send_bb(1:2, 0:procs_number-1)) ; send_bb = 0
   endassociate
   endsubroutine alloc

   pure subroutine compute_recv_cells_number(blocks)
   !< Compute receive cells number from all other processes.
   class(parallel_object), intent(inout) :: self       !< Parallel.
   type(block_object),     intent(in)    :: blocks(1:) !< Block data.
   integer(I4P)                          :: proc       !< Processes counter.
   integer(I4P)                          :: c          !< Cells counter.
   integer(I4P)                          :: b, Nb      !< Blocks counter.
   integer(I4P)                          :: i, j, k    !< Space counters.

   Nb = size(blocks, dim=1)
   ! if proc=myrank => recv_cells_number=0 => myrank doesn't communicate with itself
   do b=1, Nb
      ! i interfaces
      do k=1, blocks(b)%signature%Nk
          do j=1,blocks(b)%signature%Nj
             ! left i
             do i=0-blocks(b)%gc(1),0
                if (blocks(b)%Fi(i,j,k)%BC%tp==bc_adj) then
                   if (block_to_process_map(block(b,l)%Fi(i,j,k)%BC%adj%b)/=global%myrank) &
                      recv_cells_number(block_to_process_map(block(b,l)%Fi(i,j,k)%BC%adj%b),l) = recv_cells_number(block_to_process_map(block(b,l)%Fi(i,j,k)%BC%adj%b),l) + 1
                endif
             enddo
             ! right i
             do i=block(b,l)%Ni,block(b,l)%Ni+block(b,l)%gc(2)
               if (block(b,l)%Fi(i,j,k)%BC%tp==bc_adj) then
                  if (block_to_process_map(block(b,l)%Fi(i,j,k)%BC%adj%b)/=global%myrank) &
                     recv_cells_number(block_to_process_map(block(b,l)%Fi(i,j,k)%BC%adj%b),l) = recv_cells_number(block_to_process_map(block(b,l)%Fi(i,j,k)%BC%adj%b),l) + 1
                 endif
             enddo
          enddo
      enddo

      ! j interfaces
      do k=1,block(b,l)%Nk
        ! left j
        do j=0-block(b,l)%gc(3),0
          do i=1,block(b,l)%Ni
            if (block(b,l)%Fj(i,j,k)%BC%tp==bc_adj) then
               if (block_to_process_map(block(b,l)%Fj(i,j,k)%BC%adj%b)/=global%myrank) &
                  recv_cells_number(block_to_process_map(block(b,l)%Fj(i,j,k)%BC%adj%b),l) = recv_cells_number(block_to_process_map(block(b,l)%Fj(i,j,k)%BC%adj%b),l) + 1
              endif
          enddo
        enddo
        ! right j
        do j=block(b,l)%Nj,block(b,l)%Nj+block(b,l)%gc(4)
          do i=1,block(b,l)%Ni
            if (block(b,l)%Fj(i,j,k)%BC%tp==bc_adj) then
               if (block_to_process_map(block(b,l)%Fj(i,j,k)%BC%adj%b)/=global%myrank) &
                  recv_cells_number(block_to_process_map(block(b,l)%Fj(i,j,k)%BC%adj%b),l) = recv_cells_number(block_to_process_map(block(b,l)%Fj(i,j,k)%BC%adj%b),l) + 1
              endif
          enddo
        enddo
      enddo
      ! k interfaces
      ! left k
      do k=0-block(b,l)%gc(5),0
        do j=1,block(b,l)%Nj
          do i=1,block(b,l)%Ni
            if (block(b,l)%Fk(i,j,k)%BC%tp==bc_adj) then
               if (block_to_process_map(block(b,l)%Fk(i,j,k)%BC%adj%b)/=global%myrank) &
                  recv_cells_number(block_to_process_map(block(b,l)%Fk(i,j,k)%BC%adj%b),l) = recv_cells_number(block_to_process_map(block(b,l)%Fk(i,j,k)%BC%adj%b),l) + 1
              endif
          enddo
        enddo
      enddo
      ! right k
      do k=block(b,l)%Nk,block(b,l)%Nk+block(b,l)%gc(6)
        do j=1,block(b,l)%Nj
          do i=1,block(b,l)%Ni
            if (block(b,l)%Fk(i,j,k)%BC%tp==bc_adj) then
               if (block_to_process_map(block(b,l)%Fk(i,j,k)%BC%adj%b)/=global%myrank) &
                  recv_cells_number(block_to_process_map(block(b,l)%Fk(i,j,k)%BC%adj%b),l) = recv_cells_number(block_to_process_map(block(b,l)%Fk(i,j,k)%BC%adj%b),l) + 1
              endif
          enddo
        enddo
      enddo
   enddo

   recv_cells_number_global = sum(recv_cells_number)
   c = 0
   do proc=0, procs_number - 1
      if (recv_cells_number(proc) > 0) then
         c                = c + recv_cells_number(proc)
         recv_bb(1, proc) = c - recv_cells_number(proc) + 1
         recv_bb(2, proc) = c
      endif
   enddo
   endsubroutine compute_recv_cells_number

   subroutine compute_send_bb(self)
   !< Compute the bounding boxes of sending cells of myrank to all other processes.
   class(parallel_object), intent(inout) :: self !< Parallel.
   integer(I4P)                          :: proc !< Processes counters.
   integer(I4P)                          :: c    !< Cells counters.

   self%send_cells_number_global = sum(self%send_cells_number)
   c = 0
   do proc=0, procs_number - 1
      if (self%send_cells_number(proc)>0) then
         c               = c + send_cells_number(proc)
         send_bb(1,proc) = c - send_cells_number(proc) + 1
         send_bb(2,proc) = c
      endif
   enddo
   endsubroutine compute_send_bb

   subroutine send_recv_cells_number(self)
   !< Send receive cells number to all other processes for building send cells number locally.
   !< @note It could be better to exploit Alltoall!
   !< @note Result is to transpose matrix of receive blocks to the send one.
   class(parallel_object), intent(in) :: self !< Parallel.
   integer(I4P)                       :: l                           !< Grid levels counter.
   integer(I4P)                       :: proc                        !< Processes counters.
   integer(I4P)                       :: ierr, stat(MPI_STATUS_SIZE) !< MPI error flags.
   integer(I4P), parameter            :: tagshift=0*MAXPROC          !< Shift for tags (to isolate these kind of communications).

   do proc=0, procs_number - 1 ! the communications is organized by processes sequence
      ! sending receiving cells number of myrank to proc and
      ! using the receiving number cells of proc for building the sending number of cells of myrank
      call MPI_SENDRECV(recv_cells_number(proc), 1, MPI_INTEGER, proc, tagshift+Nproc*(myrank + 1), &
                        send_cells_number(proc), 1, MPI_INTEGER, proc, tagshift+Nproc*(proc   + 1), &
                        MPI_COMM_WORLD, stat, ierr)
   enddo
   endsubroutine send_recv_cells_number

   subroutine alloc_sendrecv_nuffers(global)
   !< Allocate send/receive buffers.
   class(parallel_object), intent(inout) :: self !< Parallel.

   ! if (allocated(recvmap)) deallocate(recvmap) ; allocate(recvmap(1:4, 1:gNcR))        ; recvmap = 0_I_P
   ! if (allocated(reqsmap)) deallocate(reqsmap) ; allocate(reqsmap(1:4, 1:gNcR))        ; reqsmap = 0_I_P
   ! if (allocated(sendmap)) deallocate(sendmap) ; allocate(sendmap(1:4, 1:gNcS))        ; sendmap = 0_I_P
   if (allocated(field_recv)) deallocate(field_recv) ; allocate(Precv  (1:recv_cells_number_blobal)) ; field_recv = 0._R_P
   if (allocated(field_send)) deallocate(field_send) ; allocate(Psend  (1:send_cells_number_blobal)) ; field_send = 0._R_P
   endsubroutine alloc_sendrecv_buffers

  subroutine compute_recv_maps(block)
  ! Subroutine for computing querying and receiving maps of actual process.
  type(Type_SBlock), intent(IN):: block(1:,1:) ! Block-level data.
  type(Type_Global), pointer::    global       ! Global-level data.
  integer(I_P)::                  l            ! Grid levels counter.
  integer(I_P)::                  proc         ! Processes counter.

  do proc=0,Nproc-1 ! the communications is organized by processes sequence
    if (proc==global%myrank) cycle ! myrank doesn't communicate with itself
    if (NcR(proc)==0) cycle ! there are no data to communicate to process proc
    call scan_proc(NcR     = NcR(proc),                              &
                   proc    = proc,                                     &
                   block   = block(:,l),                               &
                   reqsmap = reqsmap(1:4,bbR(1,proc,l):bbR(2,proc,l)), &
                   recvmap = recvmap(1:4,bbR(1,proc,l):bbR(2,proc,l)))
  enddo
  contains
    subroutine scan_proc(NcR,proc,global,block,reqsmap,recvmap)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for searching receiving cells of myrank from process proc.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),      intent(IN)::  NcR                ! Number of send/receive cells.
    integer(I_P),      intent(IN)::  proc               ! Other process than myrank to send/receive cells.
    type(Type_Global), intent(IN)::  global             ! Global-level data.
    type(Type_SBlock), intent(IN)::  block(1:global%Nb) ! Block-level data.
    integer(I_P),      intent(OUT):: reqsmap(1:4,1:NcR) ! Querying cells map.
    integer(I_P),      intent(OUT):: recvmap(1:4,1:NcR) ! Receiving cells map.
    integer(I_P)::                   c                  ! Cells counter.
    integer(I_P)::                   b                  ! Blocks counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    c = 0_I_P ! initialize cells counter
    do b=1,global%Nb ! loop over blocks
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
    gc = block%gc
    Ni = block%Ni
    Nj = block%Nj
    Nk = block%Nk
    do k=1,Nk
      do j=1,Nj
        ! left i
        if (block%Fi(0,j,k)%BC%tp==bc_adj) then
          if (procmap(block%Fi(0,j,k)%BC%adj%b)==proc) then
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
          if (procmap(block%Fi(Ni,j,k)%BC%adj%b)==proc) then
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
          if (procmap(block%Fj(i,0,k)%BC%adj%b)==proc) then
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
          if (procmap(block%Fj(i,Nj,k)%BC%adj%b)==proc) then
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
          if (procmap(block%Fk(i,j,0)%BC%adj%b)==proc) then
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
          if (procmap(block%Fk(i,j,Nk)%BC%adj%b)==proc) then
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

  subroutine mapsendrecv(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for communicate the querying cells map of myrank to other processes and for building the sending cells maps of myrank
  ! for other processes (the building is done by receiving the querying cells map of other processes).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(IN):: global                     ! Global-level data.
  integer(I_P)::                  l                          ! Grid levels counter.
  integer(I_P)::                  proc                       ! Processes counter.
  integer(I_P)::                  ierr,stat(MPI_STATUS_SIZE) ! MPI error flags.
  integer(I_P), parameter::       tagshift=1*maxproc         ! Shift for tags (to isolate these kind of communications).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do l=1,global%Nl ! loop over grid levels
    do proc=0,Nproc-1 ! the communications is organized by processes sequence
      if ((NcR(proc,l)==0).AND.(NcS(proc,l)==0)) cycle ! there are no data to communicate to process proc
#ifdef MPI2
      ! sending querying (reqsmap) cells map of myrank to proc and
      ! using the querying cells map of proc for building the sending (sendmap) cells map of myrank
      call MPI_SENDRECV(reqsmap(1,bbR(1,proc,l)),4*NcR(proc,l),MPI_INTEGER,proc,tagshift+Nproc*(global%myrank+1), &
                        sendmap(1,bbS(1,proc,l)),4*NcS(proc,l),MPI_INTEGER,proc,tagshift+Nproc*(proc         +1), &
                        MPI_COMM_WORLD,stat,ierr)
#endif
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine mapsendrecv

  !> Subroutine for initializing the send/receive communications data.
  subroutine Init_sendrecv(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_SBlock), intent(IN):: block(1:,1:) !< Block-level data.
  type(Type_Global), pointer::    global       ! Global-level data.
  integer(I_P)::                  err          !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  global => block(1,1)%global
  call alloc_SR(global=global)
  call compute_NcR(block=block)
  call NcRsendrecv(global=global)
  call compute_bbS(global=global)
  call alloc_sendrecv(global=global)
  call compute_recv_maps(block=block)
  call mapsendrecv(global=global)
  err = Printsendrecvmaps(global%myrank)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_sendrecv

  !> Subroutine for performing send/receive of primitive variables of myrank to other processes.
  subroutine Psendrecv(l,block)
     ! usare pack invece di minloc
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    l                          !< Grid level.
  type(Type_SBlock), intent(INOUT):: block(1:)                  !< Block-level data.
  type(Type_Global), pointer::       global                     ! Global-level data.
  integer(I_P)::                     proc                       !< Processes counter.
  integer(I_P)::                     c                          !< Cells counter.
  integer(I_P)::                     b                          !< Blocks counter.
  integer(I_P)::                     ierr,stat(MPI_STATUS_SIZE) !< MPI error flags.
  integer(I_P), parameter::          tagshift=2*maxproc         !< Shift for tags (to isolate these communications).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  global => block(1)%global
  do proc=0,Nproc-1 ! the communications is organized by processes sequence
    if ((NcR(proc,l)==0).AND.(NcR(proc,l)==0)) cycle ! there are no data to communicate to process proc
    ! building the send buffer (Psend) of myrank for proc using local var P
    do c=bbS(1,proc,l),bbS(2,proc,l)
      b = minloc(array=blockmap,dim=1,mask=blockmap==sendmap(1,c))
      Psend(:,c) = block(b)%C(sendmap(2,c),sendmap(3,c),sendmap(4,c))%P%prim2array()
    enddo
    ! sending Psend of myrank to proc and storing in Precv of proc
    call MPI_SENDRECV(Psend(1,bbS(1,proc,l)),global%Np*NcS(proc,l),MPI_REAL8,proc,tagshift+Nproc*(global%myrank+1), &
                      Precv(1,bbR(1,proc,l)),global%Np*NcR(proc,l),MPI_REAL8,proc,tagshift+Nproc*(proc         +1), &
                      MPI_COMM_WORLD,stat,ierr)
    ! coping the receive buffer (Precv) of myrank from proc in the local var P
    do c=bbR(1,proc,l),bbR(2,proc,l)
       call block(recvmap(1,c))%C(recvmap(2,c),recvmap(3,c),recvmap(4,c))%P%array2prim(Precv(:,c))
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Psendrecv
  !> @}

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
#endif

  !> @ingroup Lib_ParallelPublicProcedure
  !> @{
  !> Function for loading the processes/blocks map and local/global blocks map.
  !> @return \b err integer(I4P) variable.
  function procmap_load(filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),      intent(IN)::    filename !< File name processes/blocks map.
  type(Type_Global), intent(INOUT):: global   !< Global-level data.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of procmap file.
  integer(I_P)::                     b,bb     !< Blocks counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=trim(filename),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(global%myrank,filename,'procmap_load')
  open(unit = Get_Unit(UnitFree), file = trim(filename), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*,iostat=err) global%Nb_tot
  if (allocated(procmap)) deallocate(procmap) ; allocate(procmap(1:global%Nb_tot)) ; procmap  = 0_I_P
  read(UnitFree,*,iostat=err)
  do b=1,global%Nb_tot
    read(UnitFree,*,iostat=err) procmap(b) ! reading the process number of bth block
  enddo
  close(UnitFree)
  ! computing the local/global blocks map
  if (Nproc==1_I_P) then
    ! there is no MPI environment thus all blocks are loaded by process 0
    global%Nb = global%Nb_tot
    if (allocated(blockmap)) deallocate(blockmap) ; allocate(blockmap(1:global%Nb)) ; blockmap = 0_I_P
    do b=1,global%Nb_tot
      blockmap(b) = b  ! the blocks map is identity
    enddo
  else
    ! computing the local (of myrank) number of blocks
    global%Nb = 0
    do b=1,global%Nb_tot
      if (procmap(b)==global%myrank) global%Nb = global%Nb + 1
    enddo
    if (allocated(blockmap)) deallocate(blockmap) ; allocate(blockmap(1:global%Nb )) ; blockmap = 0_I_P
    bb = 0
    do b=1,global%Nb_tot
      if (procmap(b)==global%myrank) then
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
  open(unit = Get_Unit(UnitFree), file = trim(filename), form = 'FORMATTED')
  write(UnitFree,'(A)',iostat=err) str(.true.,global%Nb_tot)//' Nb_tot = number of total blocks'
  write(UnitFree,*,iostat=err)
  do b=1,global%Nb_tot
    write(UnitFree,'(A)',iostat=err) str(.true.,procmap(b))//' proc(bth) = rank of process where bth block is loaded'
  enddo
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction procmap_save
  !> @}

endmodule off_parallel_object
