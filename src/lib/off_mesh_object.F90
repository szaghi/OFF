#include "preprocessor_macros.h"
!< OFF mesh object definition and implementation.

module off_mesh_object
!< OFF mesh object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, stdout=>output_unit
use off_bc_object, only : BC_WALL, BC_PERIODIC, BC_EXTRAPOLATED, BC_ADJACENT
use off_block_object, only : block_object
use off_cell_object, only : cell_object
use off_file_grid_object, only : file_grid_object
use off_file_solution_object, only : file_solution_object
use off_grid_dimensions_object, only : grid_dimensions_object
use flow, only : eos_compressible
use off_solver_object, only : solver_object
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
   type(grid_dimensions_object)    :: grid_dimensions !< Grid dimensions.
   type(block_object), allocatable :: blocks(:)       !< Blocks list.
   contains
      ! public methods
      procedure, pass(self) :: allocate_blocks                    !< Allocate blocks accordingly to grid dimensions.
      procedure, pass(self) :: compute_residuals                  !< Compute residuals.
      procedure, pass(self) :: conservative_to_primitive          !< Convert conservative variables to primitive ones.
      procedure, pass(self) :: description                        !< Return a pretty-formatted description of the mesh.
      procedure, pass(self) :: destroy                            !< Destroy mesh.
      procedure, pass(self) :: impose_boundary_conditions         !< Impose boundary conditions.
      procedure, pass(self) :: initialize                         !< Initialize mesh.
      procedure, pass(self) :: load_grid_from_file                !< Load grid from file.
      procedure, pass(self) :: load_ic_from_file                  !< Load initial conditions from file.
      procedure, pass(self) :: primitive_to_conservative          !< Convert primitive variables to conservative ones.
      procedure, pass(self) :: save_grid_into_file                !< Save grid into file.
      procedure, pass(self) :: save_solution_into_file            !< Save solution into file.
      procedure, pass(self) :: set_parametric_boundary_conditions !< Set boundary conditions from parametric input file.
      ! fast operators
      procedure, pass(opr) :: conservative_add_conservatives_fast      !< `+` fast operator.
      procedure, pass(opr) :: conservative_multiply_conservatives_fast !< `*` fast operator.
      procedure, pass(opr) :: conservative_multiply_real_scalar_fast   !< `* real_scalar` fast operator.
      procedure, pass(opr) :: conservative_subtract_conservatives_fast !< `-` fast operator.
      ! operators
      generic :: assignment(=) => mesh_assign_mesh !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: mesh_assign_mesh !< Operator `=`.
endtype mesh_object

contains
   ! public methods
   _PURE_ subroutine allocate_blocks(self)
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

   subroutine compute_residuals(self, solver, gcu)
   !< Compute residuals.
   class(mesh_object),   intent(inout) :: self   !< Mesh.
   class(solver_object), intent(in)    :: solver !< Riemann solver.
   integer(I4P),         intent(in)    :: gcu    !< Number of ghost cells used (depend on space accuracy).
   integer(I4P)                        :: b      !< Counter.

   do b=1, self%grid_dimensions%blocks_number
      call self%blocks(b)%compute_residuals(solver=solver, gcu=gcu)
   enddo
   endsubroutine compute_residuals

   _ELEMENTAL_ subroutine conservative_to_primitive(self)
   !< Convert conservative variables to primitive one.
   class(mesh_object), intent(inout) :: self !< Mesh.
   integer(I4P)                      :: b    !< Counter.

   if (self%grid_dimensions%blocks_number>0) then
      do b=1,self%grid_dimensions%blocks_number
         call self%blocks(b)%conservative_to_primitive
      enddo
   endif
   endsubroutine conservative_to_primitive

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the mesh.
   class(mesh_object), intent(in)           :: self             !< Mesh.
   character(*),       intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable            :: desc             !< Description.
   character(len=:), allocatable            :: prefix_          !< Prefixing string, local variable.
   integer(I4P)                             :: b                !< Counter.
   character(len=1), parameter              :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//prefix_//'Grid dimensions'//NL
   desc = desc//self%grid_dimensions%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Blocks data'//NL
   if (self%grid_dimensions%blocks_number>0) then
      do b=1, self%grid_dimensions%blocks_number - 1
         desc = desc//prefix_//'  block id ['//trim(str(no_sign=.true., n=self%blocks(b)%signature%id))//']'//NL
         desc = desc//self%blocks(b)%description(prefix=prefix_//'    ')//NL
      enddo
      desc = desc//prefix_//'  block id ['//trim(str(no_sign=.true., n=self%blocks(b)%signature%id))//']'//NL
      desc = desc//self%blocks(b)%description(prefix=prefix_//'    ')
   endif
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy mesh.
   class(mesh_object), intent(inout) :: self !< Mesh.

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

   subroutine initialize(self, mesh, eos, file_grid, file_ic)
   !< Initialize mesh.
   !<
   !< The initialization can be done also loading from files (standard or paramatric ones).
   class(mesh_object),         intent(inout)           :: self      !< Mesh.
   type(mesh_object),          intent(in),    optional :: mesh      !< mesh data.
   type(eos_compressible),     intent(in),    optional :: eos       !< EOS data.
   type(file_grid_object),     intent(inout), optional :: file_grid !< Grid file.
   type(file_solution_object), intent(inout), optional :: file_ic   !< Initial conditions file.

   call self%destroy
   call self%grid_dimensions%initialize
   if (present(mesh)) then
      call self%grid_dimensions%initialize(block_signature=mesh%blocks%signature)
      allocate(self%blocks(1:size(mesh%blocks, dim=1)), source=mesh%blocks)
   else
      if (present(file_grid)) call self%load_grid_from_file(file_grid=file_grid)
      if (present(eos).and.allocated(self%blocks)) call self%blocks%set_eos(eos=eos)
      if (present(file_ic)) call self%load_ic_from_file(file_ic=file_ic)
   endif
   endsubroutine initialize

   subroutine load_grid_from_file(self, file_grid)
   !< Load grid from file.
   class(mesh_object),     intent(inout) :: self      !< Mesh.
   type(file_grid_object), intent(inout) :: file_grid !< Grid file.
   integer(I4P)                          :: b         !< Counter.

   call file_grid%load_grid_dimensions_from_file(grid_dimensions=self%grid_dimensions)
   if (self%grid_dimensions%blocks_number>0) then
      call self%allocate_blocks
      if (file_grid%is_parametric) then
         do b=1, self%grid_dimensions%blocks_number
            if (self%blocks(b)%signature%is_cartesian) then
               write(stdout, '(A)') 'create linear-space-grid for [block-'//trim(str(self%blocks(b)%signature%id, no_sign=.true.)) &
                  //']'
               call self%blocks(b)%create_linspace
            endif
         enddo
         call self%set_parametric_boundary_conditions
      else
         call file_grid%load_nodes_from_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
      endif
   endif
   endsubroutine load_grid_from_file

   subroutine load_ic_from_file(self, file_ic)
   !< Load initial conditions from file.
   class(mesh_object),         intent(inout) :: self    !< Mesh.
   type(file_solution_object), intent(inout) :: file_ic !< Initial conditions file.

   call file_ic%load_conservatives_from_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
   endsubroutine load_ic_from_file

   _ELEMENTAL_ subroutine primitive_to_conservative(self)
   !< Convert primitive variables to conservative one.
   class(mesh_object), intent(inout) :: self !< Mesh.
   integer(I4P)                      :: b    !< Counter.

   if (self%grid_dimensions%blocks_number>0) then
      do b=1,self%grid_dimensions%blocks_number
         call self%blocks(b)%primitive_to_conservative
      enddo
   endif
   endsubroutine primitive_to_conservative

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
         ! if (present(file_name)) call self%file_grid%initialize(file_name=file_name)
         ! call self%file_grid%save_grid_dimensions_into_file(grid_dimensions=self%grid_dimensions)
         ! call self%file_grid%save_nodes_into_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
      endif

      if (vtk_) then
         do b=1, self%grid_dimensions%blocks_number
            file_name_ = trim(adjustl(file_name))
            if (file_name_%basename(strip_last_extension=.true.)/='') file_name_ = file_name_%basename(strip_last_extension=.true.)
            file_name_ = file_name_//'-grid-block'//                                        &
                         '-id_'//trim(str(n=self%blocks(b)%signature%id, no_sign=.true.))// &
                         '-lv_'//trim(str(n=self%blocks(b)%signature%level, no_sign=.true.))//'.vts'
            call self%blocks(b)%save_file_grid(file_name=file_name_%chars(), ascii=ascii, metrics=metrics, vtk=vtk_)
         enddo
      endif
      endif
   endsubroutine save_grid_into_file

   subroutine save_solution_into_file(self, is_parametric, file_name, ascii, off, tecplot, vtk)
   !< Save solution into file.
   class(mesh_object), intent(inout)        :: self           !< Mesh.
   logical,            intent(in), optional :: is_parametric  !< Sentinel to load grid parametric grid file.
   character(*),       intent(in), optional :: file_name      !< File name.
   logical,            intent(in), optional :: ascii          !< Ascii/binary output.
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
      error stop 'error: mesh_object%save_solution_into_file(is_parametric=.true., ...) to be implemented'
   else
      off_ = .true.  ; if (present(off)) off_ = off
      vtk_ = .false. ; if (present(vtk)) vtk_ = vtk

      if (off_) then
         ! if (present(file_name)) call self%file_grid%initialize(file_name=file_name)
         ! call self%file_grid%save_grid_dimensions_into_file(grid_dimensions=self%grid_dimensions)
         ! call self%file_grid%save_nodes_into_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
      endif

      if (vtk_) then
         do b=1, self%grid_dimensions%blocks_number
            file_name_ = trim(adjustl(file_name))
            if (file_name_%basename(strip_last_extension=.true.)/='') file_name_ = file_name_%basename(strip_last_extension=.true.)
            file_name_ = file_name_//'-solution-block'//                                    &
                         '-id_'//trim(str(n=self%blocks(b)%signature%id, no_sign=.true.))// &
                         '-lv_'//trim(str(n=self%blocks(b)%signature%level, no_sign=.true.))//'.vts'
            call self%blocks(b)%save_file_solution(file_name=file_name_%chars(), ascii=ascii, vtk=vtk_)
         enddo
      endif
      endif
   endsubroutine save_solution_into_file

   subroutine set_parametric_boundary_conditions(self)
   !< Set boundary conditions from parametric input file.
   class(mesh_object), intent(inout) :: self           !< Mesh.
   integer(I4P)                      :: b, ba, i, j, k !< Counter.
   character(len=:), allocatable     :: buffer_c       !< Buffer character.
   type(string)                      :: buffer_s       !< Buffer string.

   do b=1, self%grid_dimensions%blocks_number
      associate(block_c=>self%blocks(b)) ! current block
         ! i left frame
         buffer_s = trim(block_c%signature%faces_bc(1))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         if (buffer_c(1:3) == 'ADJ') then
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=1, block_c%signature%nk
                  do j=1, block_c%signature%nj
                     do i=1 - block_c%signature%gc(1), 0
                        if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = block_a%signature%ni + i
                        block_c%cell(i,j,k)%bc%adj(3) = j
                        block_c%cell(i,j,k)%bc%adj(4) = k
                     enddo
                  enddo
               enddo
            endassociate
         else
            do k=1, block_c%signature%nk
               do j=1, block_c%signature%nj
                  do i=1 - block_c%signature%gc(1), 0
                     if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         endif

         ! i right frame
         buffer_s = trim(block_c%signature%faces_bc(2))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         if (buffer_c(1:3) == 'ADJ') then
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=1, block_c%signature%nk
                  do j=1, block_c%signature%nj
                     do i=block_c%signature%ni + 1, block_c%signature%ni + block_c%signature%gc(2)
                        if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = i - block_c%signature%ni
                        block_c%cell(i,j,k)%bc%adj(3) = j
                        block_c%cell(i,j,k)%bc%adj(4) = k
                     enddo
                  enddo
               enddo
            endassociate
         else
            do k=1, block_c%signature%nk
               do j=1, block_c%signature%nj
                  do i=block_c%signature%ni + 1, block_c%signature%ni + block_c%signature%gc(2)
                     if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         endif

         ! j left frame
         buffer_s = trim(block_c%signature%faces_bc(3))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         if (buffer_c(1:3) == 'ADJ') then
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=1, block_c%signature%nk
                  do j=1 - block_c%signature%gc(3), 0
                     do i=1, block_c%signature%ni
                        if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = i
                        block_c%cell(i,j,k)%bc%adj(3) = block_a%signature%nj + j
                        block_c%cell(i,j,k)%bc%adj(4) = k
                     enddo
                  enddo
               enddo
            endassociate
         else
            do k=1, block_c%signature%nk
               do j=1 - block_c%signature%gc(3), 0
                  do i=1, block_c%signature%ni
                     if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         endif

         ! j right frame
         buffer_s = trim(block_c%signature%faces_bc(4))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         if (buffer_c(1:3) == 'ADJ') then
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=1, block_c%signature%nk
                  do j=block_c%signature%nj + 1, block_c%signature%nj + block_c%signature%gc(4)
                     do i=1, block_c%signature%ni
                        if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = i
                        block_c%cell(i,j,k)%bc%adj(3) = j - block_c%signature%nj
                        block_c%cell(i,j,k)%bc%adj(4) = k
                     enddo
                  enddo
               enddo
            endassociate
         else
            do k=1, block_c%signature%nk
               do j=block_c%signature%nj + 1, block_c%signature%nj + block_c%signature%gc(4)
                  do i=1, block_c%signature%ni
                     if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         endif

         ! k left frame
         buffer_s = trim(block_c%signature%faces_bc(5))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         if (buffer_c(1:3) == 'ADJ') then
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=1 - block_c%signature%gc(5), 0
                  do j=1, block_c%signature%nj
                     do i=1, block_c%signature%ni
                        if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = i
                        block_c%cell(i,j,k)%bc%adj(3) = j
                        block_c%cell(i,j,k)%bc%adj(4) = block_a%signature%nk + k
                     enddo
                  enddo
               enddo
            endassociate
         else
            do k=1 - block_c%signature%gc(5), 0
               do j=1, block_c%signature%nj
                  do i=1, block_c%signature%ni
                     if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         endif

         ! k right frame
         buffer_s = trim(block_c%signature%faces_bc(6))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         if (buffer_c(1:3) == 'ADJ') then
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=block_c%signature%nk + 1, block_c%signature%nk + block_c%signature%gc(6)
                  do j=1, block_c%signature%nj
                     do i=1, block_c%signature%ni
                        if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = i
                        block_c%cell(i,j,k)%bc%adj(3) = j
                        block_c%cell(i,j,k)%bc%adj(4) = k - block_c%signature%nk
                     enddo
                  enddo
               enddo
            endassociate
         else
            do k=block_c%signature%nk + 1, block_c%signature%nk + block_c%signature%gc(6)
               do j=1, block_c%signature%nj
                  do i=1, block_c%signature%ni
                     if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         endif
      endassociate
   enddo
   contains
      pure function adjacent_block_index(bc_code)
      !< Return the index of adjacent block.
      character(*), intent(in) :: bc_code              !< Boundary conditions code.
      integer(I4P)             :: adjacent_block_index !< Index of adjacent block.

      ! getting adjacent block index
      if     (bc_code(1:4) == 'ADJ_') then      ! abbreviated syntax
         read(bc_code(5:), *) adjacent_block_index
      elseif (bc_code(1:9) == 'ADJACENT_') then ! verobose syntax
         read(bc_code(10:), *) adjacent_block_index
      endif
      endfunction adjacent_block_index
   endsubroutine set_parametric_boundary_conditions

   ! fast operators
   ! +
   pure subroutine conservative_add_conservatives_fast(opr, lhs, rhs)
   !< `+` fast operator.
   class(mesh_object), intent(inout) :: opr !< Operator result.
   type(mesh_object),  intent(in)    :: lhs !< Left hand side.
   type(mesh_object),  intent(in)    :: rhs !< Right hand side.
   integer(I4P)                      :: b   !< Counter.

   do b=1, opr%grid_dimensions%blocks_number
      call opr%blocks(b)%conservative_add_conservatives_fast(lhs=lhs%blocks(b), rhs=rhs%blocks(b))
   enddo
   endsubroutine conservative_add_conservatives_fast

   ! *
   pure subroutine conservative_multiply_conservatives_fast(opr, lhs, rhs)
   !< `*` fast operator.
   class(mesh_object), intent(inout) :: opr !< Operator result.
   type(mesh_object),  intent(in)    :: lhs !< Left hand side.
   type(mesh_object),  intent(in)    :: rhs !< Right hand side.
   integer(I4P)                      :: b   !< Counter.

   do b=1, opr%grid_dimensions%blocks_number
      call opr%blocks(b)%conservative_multiply_conservatives_fast(lhs=lhs%blocks(b), rhs=rhs%blocks(b))
   enddo
   endsubroutine conservative_multiply_conservatives_fast

   pure subroutine conservative_multiply_real_scalar_fast(opr, lhs, rhs)
   !< `* real_scalar` fast operator.
   class(mesh_object), intent(inout) :: opr !< Operator result.
   type(mesh_object),  intent(in)    :: lhs !< Left hand side.
   real(R8P),          intent(in)    :: rhs !< Right hand side.
   integer(I4P)                      :: b   !< Counter.

   do b=1, opr%grid_dimensions%blocks_number
      call opr%blocks(b)%conservative_multiply_real_scalar_fast(lhs=lhs%blocks(b), rhs=rhs)
   enddo
   endsubroutine conservative_multiply_real_scalar_fast

   ! -
   pure subroutine conservative_subtract_conservatives_fast(opr, lhs, rhs)
   !< `-` fast operator.
   class(mesh_object), intent(inout) :: opr !< Operator result.
   type(mesh_object),  intent(in)    :: lhs !< Left hand side.
   type(mesh_object),  intent(in)    :: rhs !< Right hand side.
   integer(I4P)                      :: b   !< Counter.

   do b=1, opr%grid_dimensions%blocks_number
      call opr%blocks(b)%conservative_subtract_conservatives_fast(lhs=lhs%blocks(b), rhs=rhs%blocks(b))
   enddo
   endsubroutine conservative_subtract_conservatives_fast

   ! private methods
   _PURE_ subroutine mesh_assign_mesh(lhs, rhs)
   !< Operator `=`.
   class(mesh_object), intent(inout) :: lhs !< Left hand side.
   type(mesh_object),  intent(in)    :: rhs !< Right hand side.
   integer(I4P)                      :: b   !< Counter.

   if (lhs%grid_dimensions%blocks_number/=rhs%grid_dimensions%blocks_number) then
      call lhs%destroy
      lhs%grid_dimensions = rhs%grid_dimensions
      call lhs%allocate_blocks
   endif
   do b=1, lhs%grid_dimensions%blocks_number
      lhs%blocks(b) = rhs%blocks(b)
   enddo
   endsubroutine mesh_assign_mesh
endmodule off_mesh_object
