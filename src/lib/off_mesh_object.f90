!< OFF mesh object definition and implementation.

module off_mesh_object
!< OFF mesh object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_bc_object, only : BC_WALL, BC_PERIODIC, BC_EXTRAPOLATED, BC_ADJACENT
use off_block_object, only : block_object
use off_cell_object, only : cell_object
use off_file_grid_object, only : file_grid_object
use off_file_parametric_grid_object, only : file_parametric_grid_object
use off_file_parametric_ic_object, only : file_parametric_ic_object
use off_grid_dimensions_object, only : grid_dimensions_object
use flow, only : eos_compressible
use penf, only : I4P, R8P, str
use stringifor, only : string
use vecfor, only : vector

implicit none
private
public :: mesh_object

type :: mesh_object
   !< mesh object class.
   !<
   !< [[mesh_object]] is a container for all mesh data.
   type(file_grid_object)            :: file_grid            !< Grid file handler.
   type(file_parametric_grid_object) :: file_parametric_grid !< Parametric grid file handler.
   type(file_parametric_ic_object)   :: file_parametric_ic   !< Parametric initial conditions file handler.
   type(grid_dimensions_object)      :: grid_dimensions      !< Grid dimensions.
   type(block_object), allocatable   :: blocks(:)            !< Blocks list.
   contains
      ! public methods
      procedure, pass(self) :: allocate_blocks            !< Allocate blocks accordingly to grid dimensions.
      procedure, pass(self) :: compute_residuals          !< Compute residuals.
      procedure, pass(self) :: conservative_to_primitive  !< Convert conservative variables to primitive ones.
      procedure, pass(self) :: description                !< Return a pretty-formatted description of the mesh.
      procedure, pass(self) :: destroy                    !< Destroy mesh.
      procedure, pass(self) :: impose_boundary_conditions !< Impose boundary conditions.
      procedure, pass(self) :: initialize                 !< Initialize mesh.
      procedure, pass(self) :: load_grid_from_file        !< Load grid from file.
      procedure, pass(self) :: load_ic_from_file          !< Load initial conditions from file.
      procedure, pass(self) :: save_grid_into_file        !< Save grid into file.
endtype mesh_object

contains
   ! public methods
   subroutine allocate_blocks(self)
   !< Allocate blocks accordingly to grid dimensions.
   class(mesh_object), intent(inout) :: self !< Mesh.
   integer(I4P)                      :: b    !< Counter.

   if (self%grid_dimensions%blocks_number > 0) then
      if (allocated(self%blocks)) then
         call self%blocks%destroy
         deallocate(self%blocks)
      endif
      allocate(self%blocks(1:self%grid_dimensions%blocks_number))
      do b=1, self%grid_dimensions%blocks_number
         call self%blocks(b)%initialize(signature=self%grid_dimensions%block_signature(b))
      enddo
   endif
   endsubroutine allocate_blocks

   subroutine compute_residuals(self, gcu)
   !< Compute residuals.
   class(mesh_object), intent(inout) :: self !< Mesh.
   integer(I4P),       intent(in)    :: gcu  !< Number of ghost cells used (depend on space accuracy).
   integer(I4P)                      :: b    !< Counter.

   do b=1, self%grid_dimensions%blocks_number
      call self%blocks(b)%compute_residuals(gcu=gcu)
   enddo
   endsubroutine compute_residuals

   elemental subroutine conservative_to_primitive(self)
   !< Convert conservative variables to primitive one.
   class(mesh_object), intent(inout) :: self !< Mesh.
   integer(I4P)                      :: b

   if (self%grid_dimensions%blocks_number>0) then
      do b=1,self%grid_dimensions%blocks_number
         call self%blocks(b)%conservative_to_primitive
      enddo
   endif
   endsubroutine conservative_to_primitive

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the mesh.
   class(mesh_object), intent(in)           :: self    !< Mesh.
   character(*),       intent(in), optional :: prefix  !< Prefixing string.
   character(len=:), allocatable            :: desc    !< Description.
   character(len=:), allocatable            :: prefix_ !< Prefixing string, local variable.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//self%grid_dimensions%description(prefix=prefix_//'  ')
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy mesh.
   class(mesh_object), intent(inout) :: self !< Mesh.

   call self%file_grid%destroy
   call self%file_parametric_grid%destroy
   call self%file_parametric_ic%destroy
   call self%grid_dimensions%destroy
   if (allocated(self%blocks)) then
      call self%blocks%destroy
      deallocate(self%blocks)
   endif
   endsubroutine destroy

   subroutine impose_boundary_conditions(self)
   !< Impose boundary conditions on all blocks of the mesh.
   class(mesh_object), intent(inout) :: self       !< Mesh.
   integer(I4P)                      :: b, i, j, k !< Counter.

   do b=1, self%grid_dimensions%blocks_number
      associate(block_b=>self%blocks(b), gc=>self%blocks(b)%signature%gc, &
                Ni=>self%blocks(b)%signature%Ni, Nj=>self%blocks(b)%signature%Nj, Nk=>self%blocks(b)%signature%Nk)
         ! i direction
         do k=1, Nk
            do j=1, Nj
               ! left
               if     (self%blocks(b)%cell(0,j,k)%bc%is(BC_WALL)) then
                  call impose_boundary_conditions_wall(gc=gc(1:2), ic=gc(1), N=Ni, normal=block_b%face_i(0,j,k)%normal, &
                                                       boundary='l', stride=block_b%cell(:,j,k))
               elseif (self%blocks(b)%cell(0,j,k)%bc%is(BC_PERIODIC)) then
                  call impose_boundary_conditions_periodic(gc=gc(1:2), ic=gc(1), N=Ni, boundary='l', stride=block_b%cell(:,j,k))
               elseif (self%blocks(b)%cell(0,j,k)%bc%is(BC_EXTRAPOLATED)) then
                  call impose_boundary_conditions_extrapolation(gc=gc(1:2), ic=gc(1), N=Ni, boundary='l', &
                                                                stride=block_b%cell(:,j,k))
               elseif (self%blocks(b)%cell(0,j,k)%bc%is(BC_ADJACENT)) then
                  call impose_boundary_conditions_adjacent(gc=gc(1), frame=block_b%cell(1-gc(1):0,j,k))
               endif
               ! right
               if     (self%blocks(b)%cell(Ni+1,j,k)%bc%is(BC_WALL)) then
                  call impose_boundary_conditions_wall(gc=gc(1:2), ic=gc(2), N=Ni, normal=block_b%face_i(Ni,j,k)%normal, &
                                                       boundary='r', stride=block_b%cell(:,j,k))
               elseif (self%blocks(b)%cell(Ni+1,j,k)%bc%is(BC_PERIODIC)) then
                  call impose_boundary_conditions_periodic(gc=gc(1:2), ic=gc(2), N=Ni, boundary='r', stride=block_b%cell(:,j,k))
               elseif (self%blocks(b)%cell(Ni+1,j,k)%bc%is(BC_EXTRAPOLATED)) then
                  call impose_boundary_conditions_extrapolation(gc=gc(1:2), ic=gc(2), N=Ni, boundary='r', &
                                                                stride=block_b%cell(:,j,k))
               elseif (self%blocks(b)%cell(Ni+1,j,k)%bc%is(BC_ADJACENT)) then
                  call impose_boundary_conditions_adjacent(gc=gc(2), frame=block_b%cell(Ni+1:Ni+gc(2),j,k))
               endif
            enddo
         enddo
         ! j direction
         do k=1, Nk
            do i=1, Ni
               ! left
               if     (self%blocks(b)%cell(i,0,k)%bc%is(BC_WALL)) then
                  call impose_boundary_conditions_wall(gc=gc(3:4), ic=gc(3), N=Nj, normal=block_b%face_i(i,0,k)%normal, &
                                                       boundary='l', stride=block_b%cell(i,:,k))
               elseif (self%blocks(b)%cell(i,0,k)%bc%is(BC_PERIODIC)) then
                  call impose_boundary_conditions_periodic(gc=gc(3:4), ic=gc(3), N=Nj, boundary='l', stride=block_b%cell(i,:,k))
               elseif (self%blocks(b)%cell(i,0,k)%bc%is(BC_EXTRAPOLATED)) then
                  call impose_boundary_conditions_extrapolation(gc=gc(3:4), ic=gc(3), N=Nj, boundary='l', &
                                                                stride=block_b%cell(i,:,k))
               elseif (self%blocks(b)%cell(i,0,k)%bc%is(BC_ADJACENT)) then
                  call impose_boundary_conditions_adjacent(gc=gc(3), frame=block_b%cell(i,1-gc(3):0,k))
               endif
               ! right
               if     (self%blocks(b)%cell(i,Nj+1,k)%bc%is(BC_WALL)) then
                  call impose_boundary_conditions_wall(gc=gc(3:4), ic=gc(4), N=Nj, normal=block_b%face_i(i,Nj,k)%normal, &
                                                       boundary='r', stride=block_b%cell(i,:,k))
               elseif (self%blocks(b)%cell(i,Nj+1,k)%bc%is(BC_PERIODIC)) then
                  call impose_boundary_conditions_periodic(gc=gc(3:4), ic=gc(4), N=Nj, boundary='r', stride=block_b%cell(i,:,k))
               elseif (self%blocks(b)%cell(i,Nj+1,k)%bc%is(BC_EXTRAPOLATED)) then
                  call impose_boundary_conditions_extrapolation(gc=gc(3:4), ic=gc(4), N=Nj, boundary='r', &
                                                                stride=block_b%cell(i,:,k))
               elseif (self%blocks(b)%cell(i,Nj+1,k)%bc%is(BC_ADJACENT)) then
                  call impose_boundary_conditions_adjacent(gc=gc(4), frame=block_b%cell(i,Nj+1:Nj+gc(4),k))
               endif
            enddo
         enddo
         ! k direction
         do j=1, Nj
            do i=1, Ni
               ! left
               if     (self%blocks(b)%cell(i,j,0)%bc%is(BC_WALL)) then
                  call impose_boundary_conditions_wall(gc=gc(5:6), ic=gc(5), N=Nk, normal=block_b%face_i(i,j,0)%normal, &
                                                       boundary='l', stride=block_b%cell(i,j,:))
               elseif (self%blocks(b)%cell(i,j,0)%bc%is(BC_PERIODIC)) then
                  call impose_boundary_conditions_periodic(gc=gc(5:6), ic=gc(5), N=Nk, boundary='l', stride=block_b%cell(i,j,:))
               elseif (self%blocks(b)%cell(i,j,0)%bc%is(BC_EXTRAPOLATED)) then
                  call impose_boundary_conditions_extrapolation(gc=gc(5:6), ic=gc(5), N=Nk, boundary='l', &
                                                                stride=block_b%cell(i,j,:))
               elseif (self%blocks(b)%cell(i,j,0)%bc%is(BC_ADJACENT)) then
                  call impose_boundary_conditions_adjacent(gc=gc(5), frame=block_b%cell(i,j,1-gc(5):0))
               endif
               ! right
               if     (self%blocks(b)%cell(i,j,Nk+1)%bc%is(BC_WALL)) then
                  call impose_boundary_conditions_wall(gc=gc(5:6), ic=gc(6), N=Nk, normal=block_b%face_i(i,Nj,k)%normal, &
                                                       boundary='r', stride=block_b%cell(i,j,:))
               elseif (self%blocks(b)%cell(i,j,Nk+1)%bc%is(BC_PERIODIC)) then
                  call impose_boundary_conditions_periodic(gc=gc(5:6), ic=gc(6), N=Nk, boundary='r', stride=block_b%cell(i,j,:))
               elseif (self%blocks(b)%cell(i,j,Nk+1)%bc%is(BC_EXTRAPOLATED)) then
                  call impose_boundary_conditions_extrapolation(gc=gc(5:6), ic=gc(6), N=Nk, boundary='r', &
                                                                stride=block_b%cell(i,j,:))
               elseif (self%blocks(b)%cell(i,j,Nk+1)%bc%is(BC_ADJACENT)) then
                  call impose_boundary_conditions_adjacent(gc=gc(6), frame=block_b%cell(i,j,Nk+1:Nk+gc(6)))
               endif
            enddo
         enddo
      endassociate
   enddo
   contains
      pure subroutine impose_boundary_conditions_wall(gc, ic, N, normal, boundary, stride)
      !< Impose wall boundary conditions on a stride of cells along a direction.
      integer(I4P),      intent(in)    :: gc(1:2)          !< Number of ghost cells.
      integer(I4P),      intent(in)    :: ic               !< Number of internal cells used for extrapolation (1 or gc).
      integer(I4P),      intent(in)    :: N                !< Number of internal cells.
      type(vector),      intent(in)    :: normal           !< Face normal.
      character(1),      intent(in)    :: boundary         !< Boundary left ('l') or right ('r').
      type(cell_object), intent(inout) :: stride(1-gc(1):) !< Cells stride [1-gc(1):N+gc(2)].
      integer(I4P)                     :: i                !< Counter.
      type(vector)                     :: vr               !< Reflected velocity vector.

      if (boundary=='l') then
         if (ic==1.or.N<gc(1)) then ! reflection using only the cell 1
            vr = stride(1)%P%velocity - (2._R8P*(stride(1)%P%velocity.paral.normal)) ! reflected velocity
            do i=1-gc(1), 0
               stride(i)%P%velocity = vr
               stride(i)%P%pressure = stride(1)%P%pressure
               stride(i)%P%density = stride(1)%P%density
            enddo
         else ! reflection using the cells 1,2,...,gc
            do i=1-gc(1), 0
               vr = stride(-i+1)%P%velocity - (2._R8P*(stride(-i+1)%P%velocity.paral.normal)) ! reflected velocity
               stride(i)%P%velocity = vr
               stride(i)%P%pressure = stride(-i+1)%P%pressure
               stride(i)%P%density = stride(-i+1)%P%density
            enddo
         endif
      endif
      if (boundary=='r') then
         if (ic==1.or.N<gc(2)) then ! reflection using only the cell N
            vr = stride(N)%P%velocity - (2._R8P*(stride(N)%P%velocity.paral.normal)) ! reflected velocity
            do i=N+1, N+gc(2)
               stride(i)%P%velocity = vr
               stride(i)%P%pressure = stride(N)%P%pressure
               stride(i)%P%density = stride(N)%P%density
            enddo
         else ! reflection using the cells N-gc,N-gc+1,N-gc+2,...,N
            do i=N+1, N+gc(2)
               vr = stride(N+1-(i-N))%P%velocity - (2._R8P*(stride(N+1-(i-N))%P%velocity.paral.normal)) ! reflected velocity
               stride(i)%P%velocity = vr
               stride(i)%P%pressure = stride(N+1-(i-N))%P%pressure
               stride(i)%P%density = stride(N+1-(i-N))%P%density
            enddo
         endif
      endif
      endsubroutine impose_boundary_conditions_wall

      pure subroutine impose_boundary_conditions_periodic(gc, ic, N, boundary, stride)
      !< Impose periodic boundary conditions on a stride of cells along a direction.
      integer(I4P),      intent(in)    :: gc(1:2)          !< Number of ghost cells.
      integer(I4P),      intent(in)    :: ic               !< Number of internal cells used for extrapolation (1 or gc).
      integer(I4P),      intent(in)    :: N                !< Number of internal cells.
      character(1),      intent(in)    :: boundary         !< Boundary left ('l') or right ('r').
      type(cell_object), intent(inout) :: stride(1-gc(1):) !< Cells stride [1-gc(1):N+gc(2)].
      integer(I4P)                     :: i                !< Cell counter.

      if (boundary=='l') then
         if (ic==1.or.N<gc(1)) then ! extrapolation using only the cell N
            do i=1-gc(1), 0
               stride(i)%P = stride(N)%P
            enddo
         else ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
            do i=1-gc(1), 0
               stride(i)%P = stride(i+N)%P
            enddo
         endif
      endif
      if (boundary=='r') then
         if (ic==1.or.N<gc(2)) then ! extrapolation using only the cell 1
            do i=N+1, N+gc(2)
               stride(i)%P = stride(1)%P
            enddo
         else ! extrapolation using the cells 1,2,...,gc
            do i=N+1, N+gc(2)
               stride(i)%P = stride(i-N)%P
            enddo
         endif
      endif
      endsubroutine impose_boundary_conditions_periodic

      pure subroutine impose_boundary_conditions_extrapolation(gc, ic, N, boundary, stride)
      !< Impose boundary conditions of extrapolation on a stride of cells along a direction.
      integer(I4P),      intent(in)    :: gc(1:2)          !< Number of ghost cells, 1 => left, 2 => right.
      integer(I4P),      intent(in)    :: ic               !< Number of internal cells used for extrapolation (1 or gc).
      integer(I4P),      intent(in)    :: N                !< Number of internal cells.
      character(1),      intent(in)    :: boundary         !< Boundary left ('l') or right ('r').
      type(cell_object), intent(inout) :: stride(1-gc(1):) !< Cells stride [1-gc(1):N+gc(2)].
      integer(I4P)                     :: i                !< Counter.

      if (boundary=='l') then
         if (ic==1.or.N<gc(1)) then ! extrapolation using only the cell 1
            do i=1-gc(1), 0
               stride(i)%P = stride(1)%P
            enddo
         else ! extrapolation using the cells 1,2,...,gc
            do i=1-gc(1), 0
               stride(i)%P = stride(-i+1)%P
            enddo
         endif
      endif
      if (boundary=='r') then
         if (ic==1.or.N<gc(2)) then ! extrapolation using only the cell N
            do i=N+1, N+gc(2)
              stride(i)%P = stride(N)%P
            enddo
         else ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
            do i=N+1, N+gc(2)
              stride(i)%P = stride(N+1-(i-N))%P
            enddo
         endif
      endif
      endsubroutine impose_boundary_conditions_extrapolation

      pure subroutine impose_boundary_conditions_adjacent(gc, frame)
      !< Impose adjacent boundary conditions on a frmae of cells along a direction.
      integer(I4P),      intent(in)    :: gc           !< Number of ghost cells.
      type(cell_object), intent(inout) :: frame(1-gc:) !< Cells frame [1-gc:0].
      integer(I4P)                     :: i            !< Counter.

      do i=1-gc, 0
         associate(adj_b=>frame(i)%bc%adj(1), adj_i=>frame(i)%bc%adj(2), adj_j=>frame(i)%bc%adj(3), adj_k=>frame(i)%bc%adj(4))
            frame(i)%P = self%blocks(adj_b)%cell(adj_i, adj_j, adj_k)%P
         endassociate
      enddo
      endsubroutine impose_boundary_conditions_adjacent
   endsubroutine impose_boundary_conditions

   subroutine initialize(self, grid_file_name, is_grid_file_parametric, ic_file_name, is_ic_file_parametric, mesh)
   !< Initialize mesh.
   !<
   !< The initialization can be done also loading from files (standard or paramatric ones).
   class(mesh_object), intent(inout)        :: self                    !< Mesh.
   character(*),       intent(in), optional :: grid_file_name          !< Grid file name.
   logical,            intent(in), optional :: is_grid_file_parametric !< Sentinel for parametric grid file.
   character(*),       intent(in), optional :: ic_file_name            !< Initial conditions file name.
   logical,            intent(in), optional :: is_ic_file_parametric   !< Sentinel for initial condition file.
   type(mesh_object),  intent(in), optional :: mesh                    !< mesh data.

   call self%destroy
   call self%grid_dimensions%initialize
   if (present(mesh)) then
      call self%grid_dimensions%initialize(block_signature=mesh%blocks%signature)
      allocate(self%blocks(1:size(mesh%blocks, dim=1)), source=mesh%blocks)
   else
      if (present(grid_file_name)) call self%load_grid_from_file(is_parametric=is_grid_file_parametric, file_name=grid_file_name)
      ! if (present(ic_file_name)) call self%load_ic_from_file(is_parametric=is_ic_file_parametric, file_name=ic_file_name)
   endif
   endsubroutine initialize

   subroutine load_grid_from_file(self, is_parametric, file_name)
   !< Load grid from file.
   class(mesh_object), intent(inout)        :: self           !< Mesh.
   logical,            intent(in), optional :: is_parametric  !< Sentinel to load grid parametric grid file.
   character(*),       intent(in), optional :: file_name      !< File name.
   logical                                  :: is_parametric_ !< Sentinel to load grid parametric grid file, local variable.
   integer(I4P)                             :: b              !< Counter.

   is_parametric_ = .false. ;  if (present(is_parametric)) is_parametric_ = is_parametric
   if (is_parametric_) then
      call self%file_parametric_grid%get_grid_dimensions(file_name=file_name, grid_dimensions=self%grid_dimensions)
      if (self%grid_dimensions%blocks_number>0) then
         call self%allocate_blocks
         do b=1, self%grid_dimensions%blocks_number
            call self%blocks(b)%create_linspace
         enddo
      endif
   else
      call self%file_grid%initialize(file_name=file_name)
      call self%file_grid%load_grid_dimensions_from_file(grid_dimensions=self%grid_dimensions)
      call self%allocate_blocks
      call self%file_grid%load_nodes_from_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
   endif
   endsubroutine load_grid_from_file

   subroutine load_ic_from_file(self, is_parametric, eos, file_name)
   !< Load initial conditions from file.
   class(mesh_object),     intent(inout)        :: self           !< Mesh.
   logical,                intent(in), optional :: is_parametric  !< Sentinel to load grid parametric grid file.
   type(eos_compressible), intent(in), optional :: eos            !< Equation of state.
   character(*),           intent(in), optional :: file_name      !< File name.
   logical                                      :: is_parametric_ !< Sentinel to load grid parametric grid file, local variable.

   is_parametric_ = .false. ;  if (present(is_parametric)) is_parametric_ = is_parametric
   if (is_parametric_.and.present(eos)) then
      call self%file_parametric_ic%load_file(file_name=file_name)
      if (self%grid_dimensions%blocks_number>0) then
         call self%file_parametric_ic%get_initial_conditions(eos=eos, blocks=self%blocks, file_name=file_name)
      else
         error stop 'error: mesh blocks have not been allocated!'
      endif
   else
      error stop 'error: loading initial conditions from restart files to be implemented!'
   endif
   endsubroutine load_ic_from_file

   subroutine save_grid_into_file(self, is_parametric, file_name, ascii, metrics, off, tecplot, vtk)
   !< Save grid into file.
   class(mesh_object), intent(inout)        :: self           !< Mesh.
   logical,            intent(in), optional :: is_parametric  !< Sentinel to load grid parametric grid file.
   character(*),       intent(in), optional :: file_name      !< File name.
   logical,            intent(in), optional :: ascii          !< Ascii/binary output.
   logical,            intent(in), optional :: metrics        !< Save also metrics data.
   logical,            intent(in), optional :: off            !< Save in OFF format sentinel.
   logical,            intent(in), optional :: tecplot        !< Tecplot output format sentinel.
   logical,            intent(in), optional :: vtk            !< VTK output format sentinel.
   logical                                  :: is_parametric_ !< Sentinel to load grid parametric grid file, local variable.
   logical                                  :: off_           !< OFF format sentinel, local variable.
   logical                                  :: vtk_           !< VTK format sentinel, local variable.
   type(string)                             :: file_name_     !< File name buffer.
   integer(I4P)                             :: b              !< Counter.

   is_parametric_ = .false. ;  if (present(is_parametric)) is_parametric_ = is_parametric
   if (is_parametric_) then
      error stop 'error: mesh_object%save_grid_into_file(is_parametric=.true., ...) to be implemented'
   else
      off_ = .true.  ; if (present(off)) off_ = off
      vtk_ = .false. ; if (present(vtk)) vtk_ = vtk

      if (off_) then
         if (present(file_name)) call self%file_grid%initialize(file_name=file_name)
         call self%file_grid%save_grid_dimensions_into_file(grid_dimensions=self%grid_dimensions)
         call self%file_grid%save_nodes_into_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
      endif

      if (vtk_) then
         do b=1, self%grid_dimensions%blocks_number
            file_name_ = trim(adjustl(file_name))
            if (file_name_%basename(strip_last_extension=.true.)/='') file_name_ = file_name_%basename(strip_last_extension=.true.)
            file_name_ = file_name_//'-block'//                                             &
                         '-id_'//trim(str(n=self%blocks(b)%signature%id, no_sign=.true.))// &
                         '-lv_'//trim(str(n=self%blocks(b)%signature%level, no_sign=.true.))//'.vts'
            call self%blocks(b)%save_file_grid(file_name=file_name_%chars(), ascii=ascii, metrics=metrics, vtk=vtk)
         enddo
      endif
      endif
   endsubroutine save_grid_into_file
endmodule off_mesh_object
