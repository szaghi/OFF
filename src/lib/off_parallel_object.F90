#include "preprocessor_macros.h"
!< OFF, parallel object definition.

module off_parallel_object
!< OFF, parallel object definition.

#ifdef _MPI_
use mpi
#endif
#ifdef _OPENMP_
use omp_lib
#endif
use off_block_object, only : block_object
use penf, only : I4P, R8P

implicit none
public :: parallel_object

type :: parallel_object
   !< OFF, parallel object class.
   ! dimensions
   integer(I4P)              :: me=0                       !< ID of current process.
   integer(I4P)              :: threads_number=1           !< Number of threads (OpenMP).
   integer(I4P)              :: processes_number=1         !< Number of processes (MPI).
   integer(I4P)              :: fields_number=0            !< Number of fields sent/received.
   integer(I4P)              :: recv_cells_number_global=0 !< Global number of receive cells (sum(recv_cells_number)).
   integer(I4P)              :: send_cells_number_global=0 !< Global number of send    cells (sum(send_cells_number)).
   integer(I4P), allocatable :: recv_cells_number(:)       !< Number of receive cells from each process [0:processes_number-1].
   integer(I4P), allocatable :: send_cells_number(:)       !< Number of send    cells from each process [0:processes_number-1].
   integer(I4P), allocatable :: recv_bb(:,:)               !< Processes bounds of receive cells [1:2,0:processes_number-1].
   integer(I4P), allocatable :: send_bb(:,:)               !< Processes bounds of send   cells  [1:2,0:processes_number-1].
   ! maps
   integer(I4P), allocatable :: block_to_process_map(:)      !< Maps of processes ID for all blocks, [1:blocks_global_number].
   integer(I4P), allocatable :: block_local_to_global_map(:) !< Maps of local block ID to global,    [1:blocks_number].
   integer(I4P), allocatable :: recv_map(:,:)                !< Receive map from other processes [1:4,1:recv_cells_number_global].
   integer(I4P), allocatable :: reqs_map(:,:)                !< Query   map to   other processes [1:4,1:recv_cells_number_global].
   integer(I4P), allocatable :: send_map(:,:)                !< Send    map to   other processes [1:4,1:send_cells_number_global].
   ! send/receive buffers of me from all other processes
   real(R8P), allocatable :: field_recv(:,:) !< Receiving buffer field [1:fields_number,1:recv_cells_number_global].
   real(R8P), allocatable :: field_send(:,:) !< Sending   buffer field [1:fields_number,1:send_cells_number_global].
   contains
      ! public methods
      procedure, pass(self) :: allocate_dimensions            !< Allocate (main) dimensions.
      procedure, pass(self) :: allocate_sendrecv_maps_buffers !< Allocate send/receive maps and buffers.
      procedure, pass(self) :: compute_recv_cells_number      !< Compute receive cells number from all other processes.
      procedure, pass(self) :: compute_recv_maps              !< Compute querying and receiving maps.
      procedure, pass(self) :: compute_send_bb                !< Compute bounding boxes of sending cells of me to all other procs.
      procedure, pass(self) :: destroy                        !< Destroy parallel.
      procedure, pass(self) :: initialize                     !< Initialize parallel.
      procedure, pass(self) :: send_recv_cells_number         !< Send recve cells number to all other processes.
      procedure, pass(self) :: send_recv_field                !< Send/receive of field variables to/from all other processes.
      procedure, pass(self) :: send_reqs_map                  !< Send query cells map to all other processes.
      ! operators
      generic :: assignment(=) => parallel_assign_parallel !< Overload `=`.
      ! private methods
      procedure, pass(lhs), private :: parallel_assign_parallel !< Operator `=`.
endtype parallel_object

contains
   pure subroutine allocate_dimensions(self)
   !< Allocate (main) dimensions.
   class(parallel_object), intent(inout) :: self !< Parallel.

   associate(pn=>self%processes_number)
      if (allocated(self%recv_cells_number)) deallocate(self%recv_cells_number) ; allocate(self%recv_cells_number(0:pn-1))
      if (allocated(self%send_cells_number)) deallocate(self%send_cells_number) ; allocate(self%send_cells_number(0:pn-1))
      if (allocated(self%recv_bb          )) deallocate(self%recv_bb          ) ; allocate(self%recv_bb(1:2,      0:pn-1))
      if (allocated(self%send_bb          )) deallocate(self%send_bb          ) ; allocate(self%send_bb(1:2,      0:pn-1))
      self%recv_cells_number = 0
      self%send_cells_number = 0
      self%recv_bb = 0
      self%send_bb = 0
   endassociate
   endsubroutine allocate_dimensions

   subroutine allocate_sendrecv_maps_buffers(self)
   !< Allocate send/receive maps and buffers.
   class(parallel_object), intent(inout) :: self !< Parallel.

   associate(rcng=>self%recv_cells_number_global, scng=>self%send_cells_number_global, fn=>self%fields_number)
      if (allocated(self%recv_map  )) deallocate(self%recv_map  ) ; allocate(self%recv_map(  1:4, 1:rcng))
      if (allocated(self%reqs_map  )) deallocate(self%reqs_map  ) ; allocate(self%reqs_map(  1:4, 1:rcng))
      if (allocated(self%send_map  )) deallocate(self%send_map  ) ; allocate(self%send_map(  1:4, 1:scng))
      if (allocated(self%field_recv)) deallocate(self%field_recv) ; allocate(self%field_recv(1:fn,1:rcng))
      if (allocated(self%field_send)) deallocate(self%field_send) ; allocate(self%field_send(1:fn,1:scng))
      self%recv_map   = 0_I4P
      self%reqs_map   = 0_I4P
      self%send_map   = 0_I4P
      self%field_recv = 0._R8P
      self%field_send = 0._R8P
   endassociate
   endsubroutine allocate_sendrecv_maps_buffers

   _PURE_ subroutine destroy(self)
   !< Destroy parallel.
   class(parallel_object), intent(inout) :: self !< Parallel.
   type(parallel_object)                 :: fresh !< Fresh instance of parallel object.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, blocks)
   !< Initialize parallel.
   class(parallel_object), intent(inout) :: self       !< Parallel.
   type(block_object),     intent(in)    :: blocks(1:) !< Block data.
   integer(I4P)                          :: mpi_error  !< MPI error flag.

   call self%destroy
#ifdef _MPI_
   call MPI_INIT(error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%me, error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%processes_number, mpi_error)
#endif

#ifdef _OPENMP_
   !$OMP PARALLEL      &
   !$OMP DEFAULT(none) &
   !$OMP SHARED(Nthreads)
   self%threads_number = OMP_GET_NUM_THREADS()
   !$OMP END PARALLEL
#endif

   self%fields_number = size(blocks(1)%cell(1,1,1)%U%array(), dim=1)
   call self%allocate_dimensions
   call self%compute_recv_cells_number(blocks=blocks)
   call self%send_recv_cells_number
   call self%compute_send_bb
   call self%allocate_sendrecv_maps_buffers
   call self%compute_recv_maps(blocks=blocks)
   call self%send_reqs_map
   !err = Printsendrecvmaps(global%myrank)
   endsubroutine initialize

   _PURE_ subroutine compute_recv_cells_number(self, blocks)
   !< Compute receive cells number from all other processes.
   class(parallel_object), intent(inout) :: self       !< Parallel.
   type(block_object),     intent(in)    :: blocks(1:) !< Block data.
   integer(I4P)                          :: proc       !< Processes counter.
   integer(I4P)                          :: c          !< Cells counter.
   integer(I4P)                          :: b          !< Blocks counter.

   do b=1, size(blocks)
      call blocks(b)%update_recv_cells_number(me=self%me,                                     &
                                              block_to_process_map=self%block_to_process_map, &
                                              recv_cells_number=self%recv_cells_number)
   enddo
   self%recv_cells_number_global = sum(self%recv_cells_number)
   c = 0
   do proc=0, self%processes_number - 1
      if (self%recv_cells_number(proc) > 0) then
         c                     = c + self%recv_cells_number(proc)
         self%recv_bb(1, proc) = c - self%recv_cells_number(proc) + 1
         self%recv_bb(2, proc) = c
      endif
   enddo
   endsubroutine compute_recv_cells_number

   _PURE_ subroutine compute_send_bb(self)
   !< Compute the bounding boxes of sending cells of me to all other processes.
   class(parallel_object), intent(inout) :: self !< Parallel.
   integer(I4P)                          :: proc !< Processes counters.
   integer(I4P)                          :: c    !< Cells counters.

   self%send_cells_number_global = sum(self%send_cells_number)
   c = 0
   do proc=0, self%processes_number - 1
      if (self%send_cells_number(proc)>0) then
         c                    = c + self%send_cells_number(proc)
         self%send_bb(1,proc) = c - self%send_cells_number(proc) + 1
         self%send_bb(2,proc) = c
      endif
   enddo
   endsubroutine compute_send_bb

   subroutine send_recv_cells_number(self)
   !< Send receive cells number to all other processes for building send cells number locally.
   !<
   !< @note It could be better to exploit `Alltoall`.
   !< @note Result is to transpose matrix of receive blocks to the send one.
   class(parallel_object), intent(inout) :: self                        !< Parallel.
#ifdef _MPI_
   integer(I4P)                          :: proc                        !< Processes counters.
   integer(I4P)                          :: mpi_error                   !< MPI error flag.
   integer(I4P)                          :: mpi_status(MPI_STATUS_SIZE) !< MPI status flags.
   integer(I4P), parameter               :: tagshift=0*MAXPROC          !< Shift for tags (to isolate these kind of communications).

   do proc=0, self%processes_number - 1 ! the communications is organized by processes sequence
      ! sending receiving cells number of me to proc and
      ! using the receiving number cells of proc for building the sending number of cells of me
      call MPI_SENDRECV(self%recv_cells_number(proc), 1, MPI_INTEGER, proc, tagshift + self%processe_number * (self%me + 1), &
                        self%send_cells_number(proc), 1, MPI_INTEGER, proc, tagshift + self%processe_number * (proc    + 1), &
                        MPI_COMM_WORLD, mpi_status, mpi_error)
   enddo
#endif
   endsubroutine send_recv_cells_number

   subroutine compute_recv_maps(self, blocks)
   !< Compute querying and receiving maps.
   class(parallel_object), intent(inout) :: self       !< Parallel.
   type(block_object),     intent(in)    :: blocks(1:) !< Block data.
   integer(I4P)                          :: proc       !< Processes counter.
   integer(I4P)                          :: c          !< Cells counter.
   integer(I4P)                          :: b          !< Blocks counter.

   do proc=0, self%processes_number - 1 ! the communications is organized by processes sequence
      if (proc==self%me) cycle ! me doesn't communicate with itself
      if (self%recv_cells_number(proc)>0) then
         c = 0_I4P ! initialize cells counter
         do b=1, size(blocks, dim=1)
            call blocks(b)%compute_recv_maps(b=b,                                                                   &
                                             process=proc,                                                          &
                                             block_to_process_map=self%block_to_process_map,                        &
                                             c=c,                                                                   &
                                             reqs_map=self%reqs_map(1:4,self%recv_bb(1,proc):self%recv_bb(2,proc)), &
                                             recv_map=self%recv_map(1:4,self%recv_bb(1,proc):self%recv_bb(2,proc)))

         enddo

      endif
   enddo
   endsubroutine compute_recv_maps

   subroutine send_reqs_map(self)
   !< Send query cells map to all other processes for building the sending cells maps of me.
   class(parallel_object), intent(inout) :: self                        !< Parallel.
#ifdef _MPI_
   integer(I4P)                          :: proc                        !< Processes counters.
   integer(I4P)                          :: mpi_error                   !< MPI error flag.
   integer(I4P)                          :: mpi_status(MPI_STATUS_SIZE) !< MPI status flags.
   integer(I_P), parameter               :: tagshift=1*maxproc          !< Shift for tags (to isolate these kind of communications).

   associate(processes_number=>self%processes_number, recv_cells_number=>self%recv_cells_number,          &
             send_cells_number=>self%send_cells_number, reqs_map=>self%reqs_map, send_map=>self%send_map, &
             recv_bb=>self%recv_bb, send_bb=>self%send_bb)
      do proc=0, processes_number - 1 ! the communications is organized by processes sequence
         if ((recv_cells_number(proc)==0).and.(send_cells_number(proc)==0)) cycle ! there are no data to communicate to process proc
         ! sending query cells map of myrank to proc and
         ! using the querying cells map of proc for building the sending (sendmap) cells map of myrank
         call MPI_SENDRECV(reqs_map(1,recv_bb(1,proc)), 4 * recv_cells_number(proc), MPI_INTEGER, &
                           proc, tagshift + processes_number * (me   + 1), &
                           send_map(1,send_bb(1,proc)), 4 * send_cells_number(proc), MPI_INTEGER, &
                           proc, tagshift + processes_number * (proc + 1), &
                           MPI_COMM_WORLD, mpi_status, mpi_error)
      enddo
   endassociate
#endif
   endsubroutine send_reqs_map

   subroutine send_recv_field(self, blocks)
   !< Send/receive of field variables to/from all other processes.
   class(parallel_object), intent(in)    :: self                        !< Parallel.
   type(block_object),     intent(inout) :: blocks(1:)                  !< Block data.
#ifdef _MPI_
   integer(I4P)                          :: proc                        !< Processes counter.
   integer(I4P)                          :: c                           !< Cells counter.
   integer(I4P)                          :: b                           !< Blocks counter.
   integer(I4P)                          :: mpi_error                   !< MPI error flag.
   integer(I4P)                          :: mpi_status(MPI_STATUS_SIZE) !< MPI status flags.
   integer(I4P), parameter               :: tagshift=2*maxproc          !< Shift for tags (to isolate these communications).

   associate(processes_number=>self%processes_number, recv_cells_number=>self%recv_cells_number,          &
             send_cells_number=>self%send_cells_number, reqs_map=>self%reqs_map, send_map=>self%send_map, &
             recv_bb=>self%recv_bb, send_bb=>self%send_bb)
      do proc=0, processes_number - 1 ! the communications is organized by processes sequence
         if ((recv_cells_number(proc)==0).and.(send_cells_number(proc)==0)) cycle ! there are no data to communicate to process proc
         ! build the send buffer for proc
         do c=send_bb(1,proc),send_bb(2,proc)
            b = minloc(array=block_local_to_global_map, dim=1, mask=block_local_to_global_map==send_map(1,c))
            field_send(:,c) = blocks(b)%cell(send_map(2,c),send_map(3,c),send_map(4,c))%U%array()
         enddo
         ! send field_send of me to proc and store in field_recv of proc
         call MPI_SENDRECV(field_send(1,send_bb(1,proc)), fields_number * send_cells_number(proc), MPI_REAL8, proc, &
                           tagshift + processes_number * (me   + 1),                                                &
                           field_recv(1,recv_bb(1,proc)), fields_number * recv_cells_number(proc), MPI_REAL8, proc, &
                           tagshift + processes_number * (proc + 1),                                                &
                           MPI_COMM_WORLD, mpi_status, mpi_error)
         ! copy the receive buffer from proc
         do c=recv_bb(1,proc),recv_bb(2,proc)
            blocks(recv_map(1,c))%cell(recv_map(2,c),recv_map(3,c),recv_map(4,c))%U = field_recv(:,c)
         enddo
      enddo
   endassociate
#endif
   endsubroutine send_recv_field

! function Printsendrecvmaps(myrank) result(err)
! ! Subroutine for send/receive primitive variables of myrank to other processes.
!
! integer(I_P), intent(IN):: myrank ! Actual rank process.
! integer(I_P)::             err    ! Error trapping flag: 0 no errors, >0 error occurs.
! integer(I_P)::             proc   ! Processes counters.
!
! write(stdout,'(A)',iostat=err)'rank'//trim(str(.true.,myrank))//&
!                               '----------------------------------------------------------------------'
! do proc=0,Nproc-1 ! the communications is organized by processes sequence
!   if (proc==myrank) cycle ! myrank doesn't communicate with itself
!   write(stdout,'(A)',IOSTAT=err) 'rank'//trim(str(.true.,myrank))//' Process '//trim(str(.true.,myrank))// &
!                                  ' must send '//trim(str(.true.,NcS(proc,1)))//  &
!                                  ' and receive '//trim(str(.true.,NcR(proc,1)))// &
!                                  ' finite voluems with process '//trim(str(.true.,proc))
! enddo
! write(stdout,'(A)',iostat=err)'rank'//trim(str(.true.,myrank))//&
!                               '----------------------------------------------------------------------'
! write(stdout,*)
! endfunction Printsendrecvmaps

! function procmap_load(filename,global) result(err)
! !> Function for loading the processes/blocks map and local/global blocks map.
! character(*),      intent(IN)::    filename !< File name processes/blocks map.
! type(Type_Global), intent(INOUT):: global   !< Global-level data.
! integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
! integer(I_P)::                     UnitFree !< Free logic unit.
! logical::                          is_file  !< Flag for inquiring the presence of procmap file.
! integer(I_P)::                     b,bb     !< Blocks counter.
!
! inquire(file=trim(filename),exist=is_file,iostat=err)
! if (.NOT.is_file) call File_Not_Found(global%myrank,filename,'procmap_load')
! open(unit = Get_Unit(UnitFree), file = trim(filename), status = 'OLD', action = 'READ', form = 'FORMATTED')
! read(UnitFree,*,iostat=err) global%Nb_tot
! if (allocated(procmap)) deallocate(procmap) ; allocate(procmap(1:global%Nb_tot)) ; procmap  = 0_I_P
! read(UnitFree,*,iostat=err)
! do b=1,global%Nb_tot
!   read(UnitFree,*,iostat=err) procmap(b) ! reading the process number of bth block
! enddo
! close(UnitFree)
! ! computing the local/global blocks map
! if (Nproc==1_I_P) then
!   ! there is no MPI environment thus all blocks are loaded by process 0
!   global%Nb = global%Nb_tot
!   if (allocated(blockmap)) deallocate(blockmap) ; allocate(blockmap(1:global%Nb)) ; blockmap = 0_I_P
!   do b=1,global%Nb_tot
!     blockmap(b) = b  ! the blocks map is identity
!   enddo
! else
!   ! computing the local (of myrank) number of blocks
!   global%Nb = 0
!   do b=1,global%Nb_tot
!     if (procmap(b)==global%myrank) global%Nb = global%Nb + 1
!   enddo
!   if (allocated(blockmap)) deallocate(blockmap) ; allocate(blockmap(1:global%Nb )) ; blockmap = 0_I_P
!   bb = 0
!   do b=1,global%Nb_tot
!     if (procmap(b)==global%myrank) then
!       bb = bb + 1
!       blockmap(bb) = b
!     endif
!   enddo
! endif
! endfunction procmap_load

! function procmap_save(filename,global) result(err)
! !> Function for saving the processes/blocks map.
! character(*),      intent(IN):: filename !< File name processes/blocks map.
! type(Type_Global), intent(IN):: global   !< Global-level data.
! integer(I_P)::                  err      !< Error trapping flag: 0 no errors, >0 error occurs.
! integer(I_P)::                  UnitFree !< Free logic unit.
! integer(I_P)::                  b        !< Block counter.
!
! ! saving the map
! open(unit = Get_Unit(UnitFree), file = trim(filename), form = 'FORMATTED')
! write(UnitFree,'(A)',iostat=err) str(.true.,global%Nb_tot)//' Nb_tot = number of total blocks'
! write(UnitFree,*,iostat=err)
! do b=1,global%Nb_tot
!   write(UnitFree,'(A)',iostat=err) str(.true.,procmap(b))//' proc(bth) = rank of process where bth block is loaded'
! enddo
! close(UnitFree)
! endfunction procmap_save

   ! private methods
   pure subroutine parallel_assign_parallel(lhs, rhs)
   !< Operator `=`.
   class(parallel_object), intent(inout) :: lhs !< Left hand side.
   type(parallel_object),  intent(in)    :: rhs !< Right hand side.

   lhs%me                       = rhs%me
   lhs%threads_number           = rhs%threads_number
   lhs%processes_number         = rhs%processes_number
   lhs%fields_number            = rhs%fields_number
   lhs%recv_cells_number_global = rhs%recv_cells_number_global
   lhs%send_cells_number_global = rhs%send_cells_number_global
   if (allocated(rhs%recv_cells_number)) then
      lhs%recv_cells_number = rhs%recv_cells_number
   else
      if (allocated(lhs%recv_cells_number)) deallocate(lhs%recv_cells_number)
   endif
   if (allocated(rhs%send_cells_number)) then
      lhs%send_cells_number = rhs%send_cells_number
   else
      if (allocated(lhs%send_cells_number)) deallocate(lhs%send_cells_number)
   endif
   if (allocated(rhs%recv_bb)) then
      lhs%recv_bb = rhs%recv_bb
   else
      if (allocated(lhs%recv_bb)) deallocate(lhs%recv_bb)
   endif
   if (allocated(rhs%send_bb)) then
      lhs%send_bb = rhs%send_bb
   else
      if (allocated(lhs%send_bb)) deallocate(lhs%send_bb)
   endif
   if (allocated(rhs%block_to_process_map)) then
      lhs%block_to_process_map = rhs%block_to_process_map
   else
      if (allocated(lhs%block_to_process_map)) deallocate(lhs%block_to_process_map)
   endif
   if (allocated(rhs%block_local_to_global_map)) then
      lhs%block_local_to_global_map = rhs%block_local_to_global_map
   else
      if (allocated(lhs%block_local_to_global_map)) deallocate(lhs%block_local_to_global_map)
   endif
   if (allocated(rhs%recv_map)) then
      lhs%recv_map = rhs%recv_map
   else
      if (allocated(lhs%recv_map)) deallocate(lhs%recv_map)
   endif
   if (allocated(rhs%reqs_map)) then
      lhs%reqs_map = rhs%reqs_map
   else
      if (allocated(lhs%reqs_map)) deallocate(lhs%reqs_map)
   endif
   if (allocated(rhs%send_map)) then
      lhs%send_map = rhs%send_map
   else
      if (allocated(lhs%send_map)) deallocate(lhs%send_map)
   endif
   if (allocated(rhs%field_recv)) then
      lhs%field_recv = rhs%field_recv
   else
      if (allocated(lhs%field_recv)) deallocate(lhs%field_recv)
   endif
   if (allocated(rhs%field_send)) then
      lhs%field_send = rhs%field_send
   else
      if (allocated(lhs%field_send)) deallocate(lhs%field_send)
   endif
   endsubroutine parallel_assign_parallel
endmodule off_parallel_object
