#include "preprocessor_macros.h"
!< OFF mesh object definition and implementation.

module off_mesh_object
!< OFF mesh object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, stdout=>output_unit
use off_bc_object, only : BC_WALL, BC_PERIODIC, BC_EXTRAPOLATED, BC_ADJACENT, BC_INLET_SUPERSONIC, &
                          BC_OUTLET_SUBSONIC, BC_INLET_SUBSONIC
use off_block_object, only : block_object
use off_cell_object, only : cell_object
use off_error_object, only : error_object
use off_file_grid_object, only : file_grid_object
use off_file_solution_object, only : file_solution_object
use off_grid_dimensions_object, only : grid_dimensions_object
use finer, only : file_ini
use flow, only : eos_compressible, primitive_compressible, primitive_to_conservative_compressible, &
                 conservative_to_primitive_compressible
use fossil, only : file_stl_object, surface_stl_object
use off_solver_object, only : solver_object
use penf, only : I4P, I8P, R8P, str, strz
use stringifor, only : string
use vecfor, only : vector, ex, ey, ez
use vtk_fortran, only : vtm_file

implicit none
private
public :: mesh_object

type :: mesh_object
   !< mesh object class.
   !<
   !< [[mesh_object]] is a container for all mesh data.
   type(error_object)              :: error           !< Errors handler.
   type(grid_dimensions_object)    :: grid_dimensions !< Grid dimensions.
   type(block_object), allocatable :: blocks(:)       !< Blocks list.
   contains
      ! public methods
      procedure, pass(self) :: allocate_blocks                    !< Allocate blocks accordingly to grid dimensions.
      procedure, pass(self) :: compute_residuals                  !< Compute residuals.
      procedure, pass(self) :: description                        !< Return a pretty-formatted description of the mesh.
      procedure, pass(self) :: destroy                            !< Destroy mesh.
      procedure, pass(self) :: impose_boundary_conditions         !< Impose boundary conditions.
      procedure, pass(self) :: initialize                         !< Initialize mesh.
      procedure, pass(self) :: load_file_solution                 !< Load file solution.
      procedure, pass(self) :: load_grid_from_file                !< Load grid from file.
      procedure, pass(self) :: load_ic_from_file                  !< Load initial conditions from file.
      procedure, pass(self) :: load_stl_geometries                !< Load STL geometries and *immerge* them into the block grid.
      procedure, pass(self) :: save_file_grid                     !< Save file grid.
      procedure, pass(self) :: save_file_solution                 !< Save file solution.
      procedure, pass(self) :: set_parametric_boundary_conditions !< Set boundary conditions from parametric input file.
      ! fast operators
      procedure, pass(opr) :: conservative_add_conservatives_fast      !< `+` fast operator.
      procedure, pass(opr) :: conservative_multiply_conservatives_fast !< `*` fast operator.
      procedure, pass(opr) :: conservative_multiply_real_scalar_fast   !< `* real_scalar` fast operator.
      procedure, pass(opr) :: conservative_subtract_conservatives_fast !< `-` fast operator.
      ! operators
      generic :: assignment(=) => mesh_assign_mesh !< Overload `=`.
      ! private methods
      procedure, pass(lhs)  :: mesh_assign_mesh                          !< Operator `=`.
      procedure, pass(self) :: set_parametric_boundary_conditions_custom !< Set boundary conditions from parametric custom file.
endtype mesh_object

contains
   ! public methods
   _PURE_ subroutine allocate_blocks(self, interfaces_number)
   !< Allocate blocks accordingly to grid dimensions.
   class(mesh_object), intent(inout)        :: self              !< Mesh.
   integer(I4P),       intent(in), optional :: interfaces_number !< Number of different interfaces.
   integer(I4P)                             :: b                 !< Counter.

   if (self%grid_dimensions%blocks_number > 0) then
      if (allocated(self%blocks)) then
         call self%blocks%destroy
         deallocate(self%blocks)
      endif
      allocate(self%blocks(1:self%grid_dimensions%blocks_number))
      do b=1, self%grid_dimensions%blocks_number
         call self%blocks(b)%initialize(signature=self%grid_dimensions%block_signature(b), interfaces_number=interfaces_number)
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

   call self%error%destroy
   call self%grid_dimensions%destroy
   if (allocated(self%blocks)) then
      call self%blocks%destroy
      deallocate(self%blocks)
   endif
   endsubroutine destroy

   _PURE_ subroutine impose_boundary_conditions(self)
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
               select case(self%blocks(b)%cell(0,j,k)%bc%id)
               case(BC_WALL)
                  call impose_boundary_conditions_wall(gc=gc(1:2), ic=gc(1), N=Ni, normal=block_b%face_i(0,j,k)%normal, &
                                                       boundary='l', stride=block_b%cell(:,j,k))
               case(BC_PERIODIC)
                  call impose_boundary_conditions_periodic(gc=gc(1:2), ic=gc(1), N=Ni, boundary='l', stride=block_b%cell(:,j,k))
               case(BC_EXTRAPOLATED)
                  call impose_boundary_conditions_extrapolation(gc=gc(1:2), ic=gc(1), N=Ni, boundary='l', &
                                                                stride=block_b%cell(:,j,k))
               case(BC_ADJACENT)
                  call impose_boundary_conditions_adjacent(gc=gc(1), frame=block_b%cell(1-gc(1):0,j,k))
               case(BC_INLET_SUPERSONIC)
                  call impose_boundary_conditions_inlet_supersonic(gc=gc(1), frame=block_b%cell(1-gc(1):0,j,k))
               case(BC_OUTLET_SUBSONIC)
                  call impose_boundary_conditions_outlet_subsonic(gc=gc(1:2), ic=gc(1), N=Ni, normal=block_b%face_i(0,j,k)%normal, &
                                                                  boundary='l', stride=block_b%cell(:,j,k))
               case(BC_INLET_SUBSONIC)
                  call impose_boundary_conditions_inlet_subsonic(gc=gc(1:2), ic=gc(1), N=Ni, normal=block_b%face_i(0,j,k)%normal, &
                                                                 boundary='l', stride=block_b%cell(:,j,k))
               endselect
               ! right
               select case(self%blocks(b)%cell(Ni+1,j,k)%bc%id)
               case(BC_WALL)
                  call impose_boundary_conditions_wall(gc=gc(1:2), ic=gc(2), N=Ni, normal=block_b%face_i(Ni,j,k)%normal, &
                                                       boundary='r', stride=block_b%cell(:,j,k))
               case(BC_PERIODIC)
                  call impose_boundary_conditions_periodic(gc=gc(1:2), ic=gc(2), N=Ni, boundary='r', stride=block_b%cell(:,j,k))
               case(BC_EXTRAPOLATED)
                  call impose_boundary_conditions_extrapolation(gc=gc(1:2), ic=gc(2), N=Ni, boundary='r', &
                                                                stride=block_b%cell(:,j,k))
               case(BC_ADJACENT)
                  call impose_boundary_conditions_adjacent(gc=gc(2), frame=block_b%cell(Ni+1:Ni+gc(2),j,k))
               case(BC_INLET_SUPERSONIC)
                  call impose_boundary_conditions_inlet_supersonic(gc=gc(2), frame=block_b%cell(Ni+1:Ni+gc(2),j,k))
               case(BC_OUTLET_SUBSONIC)
                  call impose_boundary_conditions_outlet_subsonic(gc=gc(1:2), ic=gc(2), N=Ni, normal=block_b%face_i(Ni,j,k)%normal,&
                                                                  boundary='r', stride=block_b%cell(:,j,k))
               case(BC_INLET_SUBSONIC)
                  call impose_boundary_conditions_inlet_subsonic(gc=gc(1:2), ic=gc(2), N=Ni, normal=block_b%face_i(Ni,j,k)%normal,&
                                                                 boundary='r', stride=block_b%cell(:,j,k))
               endselect
            enddo
         enddo
         ! j direction
         do k=1, Nk
            do i=1, Ni
               ! left
               select case(self%blocks(b)%cell(i,0,k)%bc%id)
               case(BC_WALL)
                  call impose_boundary_conditions_wall(gc=gc(3:4), ic=gc(3), N=Nj, normal=block_b%face_j(i,0,k)%normal, &
                                                       boundary='l', stride=block_b%cell(i,:,k))
               case(BC_PERIODIC)
                  call impose_boundary_conditions_periodic(gc=gc(3:4), ic=gc(3), N=Nj, boundary='l', stride=block_b%cell(i,:,k))
               case(BC_EXTRAPOLATED)
                  call impose_boundary_conditions_extrapolation(gc=gc(3:4), ic=gc(3), N=Nj, boundary='l', &
                                                                stride=block_b%cell(i,:,k))
               case(BC_ADJACENT)
                  call impose_boundary_conditions_adjacent(gc=gc(3), frame=block_b%cell(i,1-gc(3):0,k))
               case(BC_INLET_SUPERSONIC)
                  call impose_boundary_conditions_inlet_supersonic(gc=gc(3), frame=block_b%cell(i,1-gc(3):0,k))
               case(BC_OUTLET_SUBSONIC)
                  call impose_boundary_conditions_outlet_subsonic(gc=gc(3:4), ic=gc(3), N=Nj, normal=block_b%face_j(i,0,k)%normal, &
                                                                  boundary='l', stride=block_b%cell(i,:,k))
               case(BC_INLET_SUBSONIC)
                  call impose_boundary_conditions_inlet_subsonic(gc=gc(3:4), ic=gc(3), N=Nj, normal=block_b%face_j(i,0,k)%normal, &
                                                                 boundary='l', stride=block_b%cell(i,:,k))
               endselect
               ! right
               select case(self%blocks(b)%cell(i,Nj+1,k)%bc%id)
               case(BC_WALL)
                  call impose_boundary_conditions_wall(gc=gc(3:4), ic=gc(4), N=Nj, normal=block_b%face_j(i,Nj,k)%normal, &
                                                       boundary='r', stride=block_b%cell(i,:,k))
               case(BC_PERIODIC)
                  call impose_boundary_conditions_periodic(gc=gc(3:4), ic=gc(4), N=Nj, boundary='r', stride=block_b%cell(i,:,k))
               case(BC_EXTRAPOLATED)
                  call impose_boundary_conditions_extrapolation(gc=gc(3:4), ic=gc(4), N=Nj, boundary='r', &
                                                                stride=block_b%cell(i,:,k))
               case(BC_ADJACENT)
                  call impose_boundary_conditions_adjacent(gc=gc(4), frame=block_b%cell(i,Nj+1:Nj+gc(4),k))
               case(BC_INLET_SUPERSONIC)
                  call impose_boundary_conditions_inlet_supersonic(gc=gc(4), frame=block_b%cell(i,Nj+1:Nj+gc(4),k))
               case(BC_OUTLET_SUBSONIC)
                  call impose_boundary_conditions_outlet_subsonic(gc=gc(3:4), ic=gc(4), N=Nj, normal=block_b%face_j(i,Nj,k)%normal,&
                                                                  boundary='r', stride=block_b%cell(i,:,k))
               case(BC_INLET_SUBSONIC)
                  call impose_boundary_conditions_inlet_subsonic(gc=gc(3:4), ic=gc(4), N=Nj, normal=block_b%face_j(i,Nj,k)%normal,&
                                                                 boundary='r', stride=block_b%cell(i,:,k))
               endselect
            enddo
         enddo
         ! k direction
         do j=1, Nj
            do i=1, Ni
               ! left
               select case(self%blocks(b)%cell(i,j,0)%bc%id)
               case(BC_WALL)
                  call impose_boundary_conditions_wall(gc=gc(5:6), ic=gc(5), N=Nk, normal=block_b%face_k(i,j,0)%normal, &
                                                       boundary='l', stride=block_b%cell(i,j,:))
               case(BC_PERIODIC)
                  call impose_boundary_conditions_periodic(gc=gc(5:6), ic=gc(5), N=Nk, boundary='l', stride=block_b%cell(i,j,:))
               case(BC_EXTRAPOLATED)
                  call impose_boundary_conditions_extrapolation(gc=gc(5:6), ic=gc(5), N=Nk, boundary='l', &
                                                                stride=block_b%cell(i,j,:))
               case(BC_ADJACENT)
                  call impose_boundary_conditions_adjacent(gc=gc(5), frame=block_b%cell(i,j,1-gc(5):0))
               case(BC_INLET_SUPERSONIC)
                  call impose_boundary_conditions_inlet_supersonic(gc=gc(5), frame=block_b%cell(i,j,1-gc(5):0))
               case(BC_OUTLET_SUBSONIC)
                  call impose_boundary_conditions_outlet_subsonic(gc=gc(5:6), ic=gc(5), N=Nk, normal=block_b%face_k(i,j,0)%normal, &
                                                                  boundary='l', stride=block_b%cell(i,j,:))
               case(BC_INLET_SUBSONIC)
                  call impose_boundary_conditions_inlet_subsonic(gc=gc(5:6), ic=gc(5), N=Nk, normal=block_b%face_k(i,j,0)%normal, &
                                                                 boundary='l', stride=block_b%cell(i,j,:))
               endselect
               ! right
               select case(self%blocks(b)%cell(i,j,Nk+1)%bc%id)
               case(BC_WALL)
                  call impose_boundary_conditions_wall(gc=gc(5:6), ic=gc(6), N=Nk, normal=block_b%face_k(i,j,Nk)%normal, &
                                                       boundary='r', stride=block_b%cell(i,j,:))
               case(BC_PERIODIC)
                  call impose_boundary_conditions_periodic(gc=gc(5:6), ic=gc(6), N=Nk, boundary='r', stride=block_b%cell(i,j,:))
               case(BC_EXTRAPOLATED)
                  call impose_boundary_conditions_extrapolation(gc=gc(5:6), ic=gc(6), N=Nk, boundary='r', &
                                                                stride=block_b%cell(i,j,:))
               case(BC_ADJACENT)
                  call impose_boundary_conditions_adjacent(gc=gc(6), frame=block_b%cell(i,j,Nk+1:Nk+gc(6)))
               case(BC_INLET_SUPERSONIC)
                  call impose_boundary_conditions_inlet_supersonic(gc=gc(6), frame=block_b%cell(i,j,Nk+1:Nk+gc(6)))
               case(BC_OUTLET_SUBSONIC)
                  call impose_boundary_conditions_outlet_subsonic(gc=gc(5:6), ic=gc(6), N=Nk, normal=block_b%face_k(i,j,Nk)%normal,&
                                                                  boundary='r', stride=block_b%cell(i,j,:))
               case(BC_INLET_SUBSONIC)
                  call impose_boundary_conditions_inlet_subsonic(gc=gc(5:6), ic=gc(6), N=Nk, normal=block_b%face_k(i,j,Nk)%normal,&
                                                                 boundary='r', stride=block_b%cell(i,j,:))
               endselect
            enddo
         enddo

         ! evolve immersed boundaries, if any
         do i=1, 3!max(3, maxval(gc))
            call block_b%evolve_immersed_boundary
         enddo
      endassociate
   enddo
   contains
      _PURE_ subroutine impose_boundary_conditions_wall(gc, ic, N, normal, boundary, stride)
      !< Impose wall boundary conditions on a stride of cells along a direction.
      integer(I4P),      intent(in)    :: gc(1:2)          !< Number of ghost cells.
      integer(I4P),      intent(in)    :: ic               !< Number of internal cells used for extrapolation (1 or gc).
      integer(I4P),      intent(in)    :: N                !< Number of internal cells.
      type(vector),      intent(in)    :: normal           !< Face normal.
      character(1),      intent(in)    :: boundary         !< Boundary left ('l') or right ('r').
      type(cell_object), intent(inout) :: stride(1-gc(1):) !< Cells stride [1-gc(1):N+gc(2)].
      integer(I4P)                     :: i                !< Counter.

      if (boundary=='l') then
         if (ic==1.or.N<gc(1)) then ! reflection using only the cell 1
            do i=1-gc(1), 0
               stride(i)%U = stride(1)%U
               stride(i)%U%momentum = stride(i)%U%momentum - (2._R8P*(stride(i)%U%momentum.paral.normal))
            enddo
         else ! reflection using the cells 1,2,...,gc
            do i=1-gc(1), 0
               stride(i)%U = stride(-i+1)%U
               stride(i)%U%momentum = stride(i)%U%momentum - (2._R8P*(stride(i)%U%momentum.paral.normal))
            enddo
         endif
      endif
      if (boundary=='r') then
         if (ic==1.or.N<gc(2)) then ! reflection using only the cell N
            do i=N+1, N+gc(2)
               stride(i)%U = stride(N)%U
               stride(i)%U%momentum = stride(i)%U%momentum - (2._R8P*(stride(i)%U%momentum.paral.normal))
            enddo
         else ! reflection using the cells N-gc,N-gc+1,N-gc+2,...,N
            do i=N+1, N+gc(2)
               stride(i)%U = stride(N+1-(i-N))%U
               stride(i)%U%momentum = stride(i)%U%momentum - (2._R8P*(stride(i)%U%momentum.paral.normal))
            enddo
         endif
      endif
      endsubroutine impose_boundary_conditions_wall

      _PURE_ subroutine impose_boundary_conditions_periodic(gc, ic, N, boundary, stride)
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
               stride(i)%U = stride(N)%U
            enddo
         else ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
            do i=1-gc(1), 0
               stride(i)%U = stride(i+N)%U
            enddo
         endif
      endif
      if (boundary=='r') then
         if (ic==1.or.N<gc(2)) then ! extrapolation using only the cell 1
            do i=N+1, N+gc(2)
               stride(i)%U = stride(1)%U
            enddo
         else ! extrapolation using the cells 1,2,...,gc
            do i=N+1, N+gc(2)
               stride(i)%U = stride(i-N)%U
            enddo
         endif
      endif
      endsubroutine impose_boundary_conditions_periodic

      _PURE_ subroutine impose_boundary_conditions_extrapolation(gc, ic, N, boundary, stride)
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
               stride(i)%U = stride(1)%U
            enddo
         else ! extrapolation using the cells 1,2,...,gc
            do i=1-gc(1), 0
               stride(i)%U = stride(-i+1)%U
            enddo
         endif
      elseif (boundary=='r') then
         if (ic==1.or.N<gc(2)) then ! extrapolation using only the cell N
            do i=N+1, N+gc(2)
              stride(i)%U = stride(N)%U
            enddo
         else ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
            do i=N+1, N+gc(2)
              stride(i)%U = stride(N+1-(i-N))%U
            enddo
         endif
      endif
      endsubroutine impose_boundary_conditions_extrapolation

      _PURE_ subroutine impose_boundary_conditions_adjacent(gc, frame)
      !< Impose adjacent boundary conditions on a frmae of cells along a direction.
      integer(I4P),      intent(in)    :: gc           !< Number of ghost cells.
      type(cell_object), intent(inout) :: frame(1-gc:) !< Cells frame [1-gc:0].
      integer(I4P)                     :: i            !< Counter.

      do i=1-gc, 0
         associate(adj_b=>frame(i)%bc%adj(1), adj_i=>frame(i)%bc%adj(2), adj_j=>frame(i)%bc%adj(3), adj_k=>frame(i)%bc%adj(4))
            frame(i)%U = self%blocks(adj_b)%cell(adj_i, adj_j, adj_k)%U
         endassociate
      enddo
      endsubroutine impose_boundary_conditions_adjacent

      _PURE_ subroutine impose_boundary_conditions_inlet_supersonic(gc, frame)
      !< Impose boundary conditions of extrapolation on a stride of cells along a direction.
      integer(I4P),      intent(in)    :: gc           !< Number of ghost cells.
      type(cell_object), intent(inout) :: frame(1-gc:) !< Cells frame [1-gc:0].
      integer(I4P)                     :: i            !< Counter.

      do i=1-gc, 0
         frame(i)%U = frame(i)%bc%U
      enddo
      endsubroutine impose_boundary_conditions_inlet_supersonic

      _PURE_ subroutine impose_boundary_conditions_outlet_subsonic(gc, ic, N, normal, boundary, stride)
      !< Impose boundary conditions of subsonice outlet (fixed pressure) on a stride of cells along a direction.
      integer(I4P),      intent(in)    :: gc(1:2)          !< Number of ghost cells.
      integer(I4P),      intent(in)    :: ic               !< Number of internal cells used for extrapolation (1 or gc).
      integer(I4P),      intent(in)    :: N                !< Number of internal cells.
      type(vector),      intent(in)    :: normal           !< Face normal.
      character(1),      intent(in)    :: boundary         !< Boundary left ('l') or right ('r').
      type(cell_object), intent(inout) :: stride(1-gc(1):) !< Cells stride [1-gc(1):N+gc(2)].
      integer(I4P)                     :: i                !< Counter.
      type(primitive_compressible)     :: p_0              !< Primitive variables of free stream.
      type(primitive_compressible)     :: p_in             !< Primitive variables inside domain.
      type(primitive_compressible)     :: p_out            !< Primitive variables outside domain.

      if (boundary=='l') then
         if (ic==1.or.N<gc(1)) then ! extrapolation using only the cell 1
            do i=1-gc(1), 0
               p_in = conservative_to_primitive_compressible(conservative=stride(1)%U,    eos=self%blocks(1)%eos)
               p_0  = conservative_to_primitive_compressible(conservative=stride(i)%bc%U, eos=self%blocks(1)%eos)
               call compute_p_out_outlet_subsonic(eos=self%blocks(1)%eos, normal=normal, p_0=p_0, p_in=p_in, p_out=p_out)
               stride(i)%U = primitive_to_conservative_compressible(primitive=p_out, eos=self%blocks(1)%eos)
            enddo
         else ! extrapolation using the cells 1,2,...,gc
            do i=1-gc(1), 0
               p_in = conservative_to_primitive_compressible(conservative=stride(-i+1)%U, eos=self%blocks(1)%eos)
               p_0  = conservative_to_primitive_compressible(conservative=stride(i)%bc%U, eos=self%blocks(1)%eos)
               call compute_p_out_outlet_subsonic(eos=self%blocks(1)%eos, normal=normal, p_0=p_0, p_in=p_in, p_out=p_out)
               stride(i)%U = primitive_to_conservative_compressible(primitive=p_out, eos=self%blocks(1)%eos)
            enddo
         endif
      elseif (boundary=='r') then
         if (ic==1.or.N<gc(2)) then ! extrapolation using only the cell N
            do i=N+1, N+gc(2)
               p_in = conservative_to_primitive_compressible(conservative=stride(N)%U,    eos=self%blocks(1)%eos)
               p_0  = conservative_to_primitive_compressible(conservative=stride(i)%bc%U, eos=self%blocks(1)%eos)
               call compute_p_out_outlet_subsonic(eos=self%blocks(1)%eos, normal=normal, p_0=p_0, p_in=p_in, p_out=p_out)
               stride(i)%U = primitive_to_conservative_compressible(primitive=P_out, eos=self%blocks(1)%eos)
            enddo
         else ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
            do i=N+1, N+gc(2)
               p_in = conservative_to_primitive_compressible(conservative=stride(N+1-(i-N))%U, eos=self%blocks(1)%eos)
               p_0  = conservative_to_primitive_compressible(conservative=stride(i)%bc%U,      eos=self%blocks(1)%eos)
               call compute_p_out_outlet_subsonic(eos=self%blocks(1)%eos, normal=normal, p_0=p_0, p_in=p_in, p_out=p_out)
               stride(i)%U = primitive_to_conservative_compressible(primitive=P_out, eos=self%blocks(1)%eos)
            enddo
         endif
      endif
      endsubroutine impose_boundary_conditions_outlet_subsonic

      _PURE_ subroutine compute_p_out_outlet_subsonic(eos, normal, p_0, p_in, p_out)
      !< Compute primitive variables in the outside of domain for subsonic outlet.
      type(eos_compressible),       intent(in)  :: eos      !< Equation of state.
      type(vector),                 intent(in)  :: normal   !< Face normal.
      type(primitive_compressible), intent(in)  :: p_0      !< Primitive variables of free stream.
      type(primitive_compressible), intent(in)  :: p_in     !< Primitive variables inside domain.
      type(primitive_compressible), intent(out) :: p_out    !< Primitive variables outside domain.
      real(R8P)                                 :: a_in     !< Speed of sound inside domain.
      real(R8P)                                 :: a_0      !< Speed of sound of free stream.
      real(R8P)                                 :: a_mean   !< Mean speed of sound.
      real(R8P)                                 :: rho_mean !< Mean density.
      real(R8P)                                 :: delta    !< `(p_e-p_i)/a_mean`.

      a_in = self%blocks(1)%eos%speed_of_sound(density=p_in%density, pressure=p_in%pressure )
      a_0  = self%blocks(1)%eos%speed_of_sound(density=p_0%density,  pressure=p_0%pressure)
      a_mean   = 0.5_R8P * (a_in         + a_0        )
      rho_mean = 0.5_R8P * (p_in%density + p_0%density)
      delta = (p_0%pressure - p_in%pressure) / a_mean
      p_out%pressure = p_0%pressure
      p_out%density  = p_in%density  +          delta / a_mean
      p_out%velocity = p_in%velocity - normal * delta / rho_mean
      endsubroutine compute_p_out_outlet_subsonic

      _PURE_ subroutine impose_boundary_conditions_inlet_subsonic(gc, ic, N, normal, boundary, stride)
      !< Impose boundary conditions of subsonic inlet (fixed velocity and density) on a stride of cells along a direction.
      integer(I4P),      intent(in)    :: gc(1:2)          !< Number of ghost cells.
      integer(I4P),      intent(in)    :: ic               !< Number of internal cells used for extrapolation (1 or gc).
      integer(I4P),      intent(in)    :: N                !< Number of internal cells.
      type(vector),      intent(in)    :: normal           !< Face normal.
      character(1),      intent(in)    :: boundary         !< Boundary left ('l') or right ('r').
      type(cell_object), intent(inout) :: stride(1-gc(1):) !< Cells stride [1-gc(1):N+gc(2)].
      integer(I4P)                     :: i                !< Counter.
      type(primitive_compressible)     :: p_0              !< Primitive variables of free stream.
      type(primitive_compressible)     :: p_in             !< Primitive variables inside domain.
      type(primitive_compressible)     :: p_out            !< Primitive variables outside domain.

      if (boundary=='l') then
         if (ic==1.or.N<gc(1)) then ! extrapolation using only the cell 1
            do i=1-gc(1), 0
               p_in = conservative_to_primitive_compressible(conservative=stride(1)%U,    eos=self%blocks(1)%eos)
               p_0  = conservative_to_primitive_compressible(conservative=stride(i)%bc%U, eos=self%blocks(1)%eos)
               call compute_p_out_inlet_subsonic(eos=self%blocks(1)%eos, normal=normal, p_0=p_0, p_in=p_in, p_out=p_out)
               stride(i)%U = primitive_to_conservative_compressible(primitive=P_out, eos=self%blocks(1)%eos)
            enddo
         else ! extrapolation using the cells 1,2,...,gc
            do i=1-gc(1), 0
               p_in = conservative_to_primitive_compressible(conservative=stride(-i+1)%U, eos=self%blocks(1)%eos)
               p_0  = conservative_to_primitive_compressible(conservative=stride(i)%bc%U, eos=self%blocks(1)%eos)
               call compute_p_out_inlet_subsonic(eos=self%blocks(1)%eos, normal=normal, p_0=p_0, p_in=p_in, p_out=p_out)
               stride(i)%U = primitive_to_conservative_compressible(primitive=P_out, eos=self%blocks(1)%eos)
            enddo
         endif
      elseif (boundary=='r') then
         if (ic==1.or.N<gc(2)) then ! extrapolation using only the cell N
            do i=N+1, N+gc(2)
               p_in  = conservative_to_primitive_compressible(conservative=stride(N)%U,    eos=self%blocks(1)%eos)
               p_0   = conservative_to_primitive_compressible(conservative=stride(i)%bc%U, eos=self%blocks(1)%eos)
               call compute_p_out_inlet_subsonic(eos=self%blocks(1)%eos, normal=normal, p_0=p_0, p_in=p_in, p_out=p_out)
               stride(i)%U = primitive_to_conservative_compressible(primitive=P_out, eos=self%blocks(1)%eos)
            enddo
         else ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
            do i=N+1, N+gc(2)
               p_in  = conservative_to_primitive_compressible(conservative=stride(N+1-(i-N))%U, eos=self%blocks(1)%eos)
               p_0   = conservative_to_primitive_compressible(conservative=stride(i)%bc%U,      eos=self%blocks(1)%eos)
               call compute_p_out_inlet_subsonic(eos=self%blocks(1)%eos, normal=normal, p_0=p_0, p_in=p_in, p_out=p_out)
               stride(i)%U = primitive_to_conservative_compressible(primitive=P_out, eos=self%blocks(1)%eos)
            enddo
         endif
      endif
      endsubroutine impose_boundary_conditions_inlet_subsonic

      _PURE_ subroutine compute_p_out_inlet_subsonic(eos, normal, p_0, p_in, p_out)
      !< Compute primitive variables in the outside of domain for subsonic inlet.
      type(eos_compressible),       intent(in)  :: eos      !< Equation of state.
      type(vector),                 intent(in)  :: normal   !< Face normal.
      type(primitive_compressible), intent(in)  :: p_0      !< Primitive variables of free stream.
      type(primitive_compressible), intent(in)  :: p_in     !< Primitive variables inside domain.
      type(primitive_compressible), intent(out) :: p_out    !< Primitive variables outside domain.
      real(R8P)                                 :: a_in     !< Speed of sound inside domain.
      real(R8P)                                 :: a_0      !< Speed of sound of free stream.
      real(R8P)                                 :: a_mean   !< Mean speed of sound.
      real(R8P)                                 :: rho_mean !< Mean density.
      real(R8P)                                 :: delta    !< `(p_e-p_i)/a_mean`.

      a_in  = eos%speed_of_sound(density=P_in%density,  pressure=P_in%pressure )
      a_0 = eos%speed_of_sound(density=P_0%density, pressure=P_0%pressure)
      a_mean   = 0.5_R8P * (a_in         + a_0        )
      rho_mean = 0.5_R8P * (p_in%density + p_0%density)
      p_out%pressure = 0.5_R8P * (p_0%pressure + p_in%pressure - rho_mean * a_mean * ((p_0%velocity - p_in%velocity).dot.normal))
      delta = (p_out%pressure - p_0%pressure) / a_mean
      p_out%density  = p_0%density  +          delta / a_mean
      p_out%velocity = p_0%velocity + normal * delta / rho_mean
      endsubroutine compute_p_out_inlet_subsonic
   endsubroutine impose_boundary_conditions

   subroutine initialize(self, mesh, eos, interfaces_number, file_grid, file_ic)
   !< Initialize mesh.
   !<
   !< The initialization can be done also loading from files (standard or paramatric ones).
   class(mesh_object),         intent(inout)           :: self              !< Mesh.
   type(mesh_object),          intent(in),    optional :: mesh              !< mesh data.
   type(eos_compressible),     intent(in),    optional :: eos               !< EOS data.
   integer(I4P),               intent(in),    optional :: interfaces_number !< Number of different interfaces.
   type(file_grid_object),     intent(inout), optional :: file_grid         !< Grid file.
   type(file_solution_object), intent(inout), optional :: file_ic           !< Initial conditions file.

   call self%destroy
   call self%grid_dimensions%initialize
   if (present(mesh)) then
      call self%grid_dimensions%initialize(block_signature=mesh%blocks%signature)
      allocate(self%blocks(1:size(mesh%blocks, dim=1)), source=mesh%blocks)
   else
      if (present(file_grid)) call self%load_grid_from_file(file_grid=file_grid, eos=eos, interfaces_number=interfaces_number)
      if (present(file_ic)) call self%load_ic_from_file(file_ic=file_ic)
   endif
   endsubroutine initialize

   subroutine load_file_solution(self, file_solution, file_name, n)
   !< Load file solution.
   class(mesh_object),         intent(inout)        :: self           !< Mesh.
   type(file_solution_object), intent(in), optional :: file_solution  !< File solution handler.
   character(*),               intent(in), optional :: file_name      !< File name.
   integer(I8P),               intent(in), optional :: n              !< Time step.
   type(file_solution_object)                       :: file_solution_ !< File solution handler, local variable.
   type(string)                                     :: base_name      !< Base file name.
   type(string)                                     :: file_name_     !< File name, local variable.
   type(string)                                     :: n_             !< Time step string.

   ! building file base name
   if (present(file_name)) then
      file_name_ = trim(adjustl(file_name))
   else
      if (present(file_solution)) then
         if (allocated(file_solution%file_name)) then
            file_name_ = trim(adjustl(file_solution%file_name))
         else
            error stop 'error: mesh_object%load_file_solution needs a file name to be passed (explicitely or via file_solution)'
         endif
      else
         error stop 'error: mesh_object%load_file_solution needs a file name to be passed (explicitely or via file_solution)'
      endif
   endif
   base_name = file_name_
   if (base_name%basename(strip_last_extension=.true.)/='') base_name = base_name%basename(strip_last_extension=.true.)
   n_ = '' ; if (present(n)) n_ = '-n_'//trim(strz(n=n, nz_pad=15))

   if (present(file_solution)) then
      file_solution_ = file_solution
   else
      call file_solution_%initialize(file_name=file_name_%chars())
   endif

   call file_solution_%load_grid_dimensions_from_file(grid_dimensions=self%grid_dimensions)
   call file_solution_%load_conservatives_from_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
   endsubroutine load_file_solution

   subroutine load_grid_from_file(self, eos, interfaces_number, file_grid)
   !< Load grid from file.
   class(mesh_object),     intent(inout)           :: self              !< Mesh.
   type(eos_compressible), intent(in),    optional :: eos               !< EOS data.
   integer(I4P),           intent(in),    optional :: interfaces_number !< Number of different interfaces.
   type(file_grid_object), intent(inout)           :: file_grid         !< Grid file.
   integer(I4P)                                    :: b                 !< Counter.

   call file_grid%load_grid_dimensions_from_file(grid_dimensions=self%grid_dimensions)
   if (self%grid_dimensions%blocks_number>0) then
      call self%allocate_blocks(interfaces_number=interfaces_number)
      if (file_grid%is_parametric) then
         if (present(eos).and.allocated(self%blocks)) call self%blocks%set_eos(eos=eos)
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

   subroutine load_stl_geometries(self, fini, go_on_fail)
   !< Load geometries described by STL triangulated surface and *immerge* them into the block grid.
   class(mesh_object), intent(inout)        :: self                    !< Mesh.
   type(file_ini),     intent(in)           :: fini                    !< Simulation parameters ini file handler.
   logical,            intent(in), optional :: go_on_fail              !< Go on if load fails.
   logical                                  :: go_on_fail_             !< Go on if load fails, local variable.
   type(file_stl_object)                    :: file_stl                !< STL file handler.
   type(surface_stl_object)                 :: surface_stl             !< STL surface handler.
   character(len=:), allocatable            :: section                 !< Section of INI file containing STL files data.
   integer(I4P)                             :: files_number            !< Number of IB body files.
   character(999)                           :: file_name               !< File name of IB body files.
   integer(I4P)                             :: excluded_blocks_number  !< Excluded blocks number from geometry immerse.
   integer(I4P), allocatable                :: excluded_blocks_list(:) !< Excluded blocks list.
   real(R8P)                                :: resize_factor(3)        !< Resize factor.
   logical                                  :: distance_sign_inverse   !< Flag to invert sign of geometries distances.
   integer(I4P)                             :: aabb_ref_levels         !< AABB refinement levels.
   integer(I4P)                             :: b, f                    !< Counter.

   if (self%grid_dimensions%blocks_number>0) then
      go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
      section = 'stl_geometries'
      if (fini%has_section(section)) then
         call fini%get(section_name=section, option_name='files_number', val=files_number, error=self%error%status)
         call self%error%check(message='failed to load ['//section//'].(files_number)', is_severe=.not.go_on_fail_)
         if (files_number>0) then
            do f=1, files_number
               call fini%get(section_name=section, option_name='file_stl_'//trim(str(n=f, no_sign=.true.)), &
                             val=file_name, error=self%error%status)
               call self%error%check(message='failed to load ['//section//'].('//'file_stl_'//trim(str(n=f, no_sign=.true.))//')', &
                                     is_severe=.not.go_on_fail_)

               excluded_blocks_number = 0
               if (allocated(excluded_blocks_list)) deallocate(excluded_blocks_list)
               call fini%get(section_name=section, option_name='excluded_blocks_number_stl_'//trim(str(n=f, no_sign=.true.)), &
                             val=excluded_blocks_number, error=self%error%status)
               if (excluded_blocks_number>0) then
                  allocate(excluded_blocks_list(1:excluded_blocks_number))
                  call fini%get(section_name=section, option_name='excluded_blocks_list_stl_'//trim(str(n=f, no_sign=.true.)), &
                                val=excluded_blocks_list, error=self%error%status)
               endif

               call fini%get(section_name=section, option_name='distance_sign_inverse_stl_'//trim(str(n=f, no_sign=.true.)), &
                             val=distance_sign_inverse, error=self%error%status)
               if (self%error%status /= 0) distance_sign_inverse = .false.

               call fini%get(section_name=section, option_name='resize_factor_stl_'//trim(str(n=f, no_sign=.true.)), &
                             val=resize_factor, error=self%error%status)
               if (self%error%status /= 0) resize_factor = 1._R8P

               call fini%get(section_name=section, option_name='aabb_ref_levels_stl_'//trim(str(n=f, no_sign=.true.)), &
                             val=aabb_ref_levels, error=self%error%status)
               if (self%error%status /= 0) aabb_ref_levels = 3

               call file_stl%load_from_file(facet=surface_stl%facet, file_name=trim(adjustl(file_name)), guess_format=.true.)
               print '(A)', file_stl%statistics()
               call surface_stl%resize(factor=resize_factor(1) * ex + resize_factor(2) * ey + resize_factor(3) * ez)
               ! call surface_stl%analize(aabb_refinement_levels=aabb_ref_levels)
               ! call surface_stl%sanitize
               call surface_stl%analize(aabb_refinement_levels=aabb_ref_levels)
               print '(A)', surface_stl%statistics()

               do b=1, self%grid_dimensions%blocks_number
                  if (allocated(excluded_blocks_list)) then
                     if (any(excluded_blocks_list==b)) cycle
                  endif
                  call self%blocks(b)%immerge_stl_geometry(surface_stl=surface_stl, n=f,                &
                                                           distance_sign_inverse=distance_sign_inverse)
               enddo
            enddo
            do b=1, self%grid_dimensions%blocks_number
               call self%blocks(b)%update_level_set_distance
            enddo
         endif
      endif
   endif
   endsubroutine load_stl_geometries

   subroutine save_file_grid(self, file_grid, file_name, is_parametric, metrics, ascii, off, tecplot, vtk, n, force)
   !< Save file grid.
   class(mesh_object),     intent(inout)        :: self           !< Mesh.
   type(file_grid_object), intent(in), optional :: file_grid      !< File grid handler.
   character(*),           intent(in), optional :: file_name      !< File name.
   logical,                intent(in), optional :: is_parametric  !< Sentinel to load grid parametric grid file.
   logical,                intent(in), optional :: metrics        !< Save also metrics data.
   logical,                intent(in), optional :: ascii          !< Ascii/binary output.
   logical,                intent(in), optional :: off            !< Save in OFF format sentinel.
   logical,                intent(in), optional :: tecplot        !< Tecplot output format sentinel.
   logical,                intent(in), optional :: vtk            !< VTK output format sentinel.
   integer(I8P),           intent(in), optional :: n              !< Time step.
   logical,                intent(in), optional :: force          !< Sentinel to force saving.
   type(file_grid_object)                       :: file_grid_     !< File grid handler, local variable.
   logical                                      :: is_parametric_ !< Sentinel to load grid parametric grid file, local variable.
   type(string)                                 :: base_name      !< Base file name.
   type(string)                                 :: file_name_     !< File name, local variable.
   logical                                      :: metrics_       !< Save metrics sentinel, local variable.
   logical                                      :: ascii_         !< Ascii/binary output, local variable.
   logical                                      :: off_           !< OFF format sentinel, local variable.
   logical                                      :: vtk_           !< VTK format sentinel, local variable.
   type(string)                                 :: n_             !< Time step string.
   logical                                      :: force_         !< Sentinel to force saving.
   integer(I4P)                                 :: b              !< Counter.

   ! initialize sentinel for forcing saving
   force_ = .false. ; if (present(force)) force_ = force

   ! check save frequency if it has a meaning
   if (present(file_grid).and.present(n).and.(.not.force_)) then
      if (file_grid%save_frequency <= 0) return
      if (mod(n, file_grid%save_frequency) /= 0) return
   endif

   ! building file base name
   if (present(file_name)) then
      file_name_ = trim(adjustl(file_name))
   else
      if (present(file_grid)) then
         if (allocated(file_grid%file_name)) then
            file_name_ = trim(adjustl(file_grid%file_name))
         else
            error stop 'error: mesh_object%save_file_grid needs a file name to be passed (explicitely or via file_grid)'
         endif
      else
         error stop 'error: mesh_object%save_file_grid needs a file name to be passed (explicitely or via file_grid)'
      endif
   endif
   base_name = file_name_
   if (base_name%basename(strip_last_extension=.true.)/='') base_name = base_name%basename(strip_last_extension=.true.)
   n_ = '' ; if (present(n)) n_ = '-n_'//trim(strz(n=n, nz_pad=15))

   if (present(file_grid)) then
      file_grid_ = file_grid
      is_parametric_ = file_grid%is_parametric
      metrics_ = file_grid%save_metrics
      ascii_ = file_grid%ascii_format
      off_ = file_grid%off_format
      vtk_ = file_grid%vtk_format
   else
      is_parametric_ = .false. ; if (present(is_parametric)) is_parametric_ = is_parametric
      metrics_ = .false. ; if (present(metrics)) metrics_ = metrics
      ascii_ = .true. ; if (present(ascii)) ascii_ = ascii
      off_ = .true. ; if (present(off)) off_ = off
      vtk_ = .false. ; if (present(vtk)) vtk_ = vtk
      call file_grid_%initialize(file_name=file_name_%chars(), is_parametric=is_parametric_)
      file_grid_%save_metrics = metrics_
      file_grid_%ascii_format = ascii_
      file_grid_%off_format   = off_
      file_grid_%vtk_format   = vtk_
   endif

   if (is_parametric_) then
      error stop 'error: mesh_object%save_file_grid(is_parametric=.true., ...) to be implemented'
   else
      if (off_) then
         call file_grid_%save_grid_dimensions_into_file(grid_dimensions=self%grid_dimensions)
         call file_grid_%save_nodes_into_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
      endif

      if (vtk_) then
         do b=1, self%grid_dimensions%blocks_number
            file_name_ = base_name//'-grid-block'//                                         &
                         '-id_'//trim(str(n=self%blocks(b)%signature%id, no_sign=.true.))// &
                         '-lv_'//trim(str(n=self%blocks(b)%signature%level, no_sign=.true.))//n_//'.vts'
            call self%blocks(b)%save_file_grid(file_name=file_name_%chars(), metrics=metrics, ascii=ascii, vtk=vtk_)
         enddo
      endif
   endif
   endsubroutine save_file_grid

   subroutine save_file_solution(self, file_solution, file_name, is_parametric, metrics, ascii, gc, off, tecplot, vtk, n, force)
   !< Save file solution.
   class(mesh_object),         intent(inout)        :: self           !< Mesh.
   type(file_solution_object), intent(in), optional :: file_solution  !< File solution handler.
   character(*),               intent(in), optional :: file_name      !< File name.
   logical,                    intent(in), optional :: is_parametric  !< Sentinel to load grid parametric grid file.
   logical,                    intent(in), optional :: metrics        !< Save metrics sentinel.
   logical,                    intent(in), optional :: gc             !< Save ghost cells sentinel.
   logical,                    intent(in), optional :: ascii          !< Ascii/binary output.
   logical,                    intent(in), optional :: off            !< Save in OFF format sentinel.
   logical,                    intent(in), optional :: tecplot        !< Tecplot output format sentinel.
   logical,                    intent(in), optional :: vtk            !< VTK output format sentinel.
   integer(I8P),               intent(in), optional :: n              !< Time step.
   logical,                    intent(in), optional :: force          !< Sentinel to force saving.
   type(file_solution_object)                       :: file_solution_ !< File solution handler, local variable.
   logical                                          :: is_parametric_ !< Sentinel to load grid parametric grid file, local variable.
   type(string)                                     :: base_name      !< Base file name.
   type(string)                                     :: file_name_     !< File name, local variable.
   type(string)                                     :: file_names_    !< List of file names
   logical                                          :: metrics_       !< Save metrics sentinel, local variable.
   logical                                          :: gc_            !< Save ghost cells sentinel, local variable.
   logical                                          :: ascii_         !< Ascii/binary output, local variable.
   logical                                          :: off_           !< OFF format sentinel, local variable.
   logical                                          :: vtk_           !< VTK format sentinel, local variable.
   type(string)                                     :: n_             !< Time step string.
   logical                                          :: force_         !< Sentinel to force saving, local variable.
   type(vtm_file)                                   :: vtm_file_      !< VTM file handler.
   integer(I4P)                                     :: b              !< Counter.

   ! initialize sentinel for forcing saving
   force_ = .false. ; if (present(force)) force_ = force

   ! check save frequency if it has a meaning
   if (present(file_solution).and.present(n).and.(.not.force_)) then
      if (file_solution%save_frequency <= 0) return
      if (mod(n, file_solution%save_frequency) /= 0) return
   endif

   ! building file base name
   if (present(file_name)) then
      file_name_ = trim(adjustl(file_name))
   else
      if (present(file_solution)) then
         if (allocated(file_solution%file_name)) then
            file_name_ = trim(adjustl(file_solution%file_name))
         else
            error stop 'error: mesh_object%save_file_solution needs a file name to be passed (explicitely or via file_solution)'
         endif
      else
         error stop 'error: mesh_object%save_file_solution needs a file name to be passed (explicitely or via file_solution)'
      endif
   endif
   base_name = file_name_
   if (base_name%basename(strip_last_extension=.true.)/='') base_name = base_name%basename(strip_last_extension=.true.)
   n_ = '' ; if (present(n)) n_ = '-n_'//trim(strz(n=n, nz_pad=15))

   if (present(file_solution)) then
      file_solution_ = file_solution
      is_parametric_ = file_solution%is_parametric
      metrics_ = file_solution%save_metrics
      gc_ = file_solution%save_ghost_cells
      ascii_ = file_solution%ascii_format
      off_ = file_solution%off_format
      vtk_ = file_solution%vtk_format
   else
      is_parametric_ = .false. ; if (present(is_parametric)) is_parametric_ = is_parametric
      metrics_ = .false. ; if (present(metrics)) metrics_ = metrics
      gc_ = .false. ; if (present(gc)) gc_ = gc
      ascii_ = .false. ; if (present(ascii)) ascii_ = ascii
      off_ = .false. ; if (present(off)) off_ = off
      vtk_ = .false. ; if (present(vtk)) vtk_ = vtk
      call file_solution_%initialize(file_name=file_name_%chars(), is_parametric=is_parametric_)
      file_solution_%save_metrics = metrics_
      file_solution_%ascii_format = ascii_
      file_solution_%off_format   = off_
      file_solution_%vtk_format   = vtk_
   endif

   if (is_parametric_) then
      error stop 'error: mesh_object%save_file_solution(is_parametric=.true., ...) to be implemented'
   else
      if (off_) then
         call file_solution_%save_grid_dimensions_into_file(grid_dimensions=self%grid_dimensions)
         call file_solution_%save_conservatives_into_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
      endif

      if (vtk_) then
         file_names_ = ''
         do b=1, self%grid_dimensions%blocks_number
            file_name_ = base_name//'-solution-block'//                                     &
                         '-id_'//trim(str(n=self%blocks(b)%signature%id, no_sign=.true.))// &
                         '-lv_'//trim(str(n=self%blocks(b)%signature%level, no_sign=.true.))//n_//'.vts'
            call self%blocks(b)%save_file_solution(file_name=file_name_%chars(), metrics=metrics_, gc=gc_, ascii=ascii_, vtk=vtk_)
            file_names_ = file_names_//' '//file_name_
         enddo
         self%error%status = vtm_file_%initialize(filename=base_name//'-solution'//n_//'.vtm')
         self%error%status = vtm_file_%write_block(filenames=trim(file_names_%chars()), name='blocks')
         self%error%status = vtm_file_%finalize()
      endif
   endif
   endsubroutine save_file_solution

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
         select case(buffer_c(1:3))
         case('ADJ')
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=1, block_c%signature%nk
                  do j=1, block_c%signature%nj
                     do i=1 - block_c%signature%gc(1), 0
                        ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = block_a%signature%ni + i
                        block_c%cell(i,j,k)%bc%adj(3) = j
                        block_c%cell(i,j,k)%bc%adj(4) = k
                     enddo
                  enddo
               enddo
            endassociate
         case('FRE', 'WAL', 'PER', 'EXT')
            do k=1, block_c%signature%nk
               do j=1, block_c%signature%nj
                  do i=1 - block_c%signature%gc(1), 0
                     ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         case default
            ! bc is defined into a custom INI file
            call self%set_parametric_boundary_conditions_custom(file_name=trim(block_c%signature%faces_bc(1)))
         endselect

         ! i right frame
         buffer_s = trim(block_c%signature%faces_bc(2))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         select case(buffer_c(1:3))
         case('ADJ')
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=1, block_c%signature%nk
                  do j=1, block_c%signature%nj
                     do i=block_c%signature%ni + 1, block_c%signature%ni + block_c%signature%gc(2)
                        ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = i - block_c%signature%ni
                        block_c%cell(i,j,k)%bc%adj(3) = j
                        block_c%cell(i,j,k)%bc%adj(4) = k
                     enddo
                  enddo
               enddo
            endassociate
         case('FRE', 'WAL', 'PER', 'EXT')
            do k=1, block_c%signature%nk
               do j=1, block_c%signature%nj
                  do i=block_c%signature%ni + 1, block_c%signature%ni + block_c%signature%gc(2)
                     ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         case default
            ! bc is defined into a custom INI file
            call self%set_parametric_boundary_conditions_custom(file_name=trim(block_c%signature%faces_bc(2)))
         endselect

         ! j left frame
         buffer_s = trim(block_c%signature%faces_bc(3))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         select case(buffer_c(1:3))
         case('ADJ')
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=1, block_c%signature%nk
                  do j=1 - block_c%signature%gc(3), 0
                     do i=1, block_c%signature%ni
                        ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = i
                        block_c%cell(i,j,k)%bc%adj(3) = block_a%signature%nj + j
                        block_c%cell(i,j,k)%bc%adj(4) = k
                     enddo
                  enddo
               enddo
            endassociate
         case('FRE', 'WAL', 'PER', 'EXT')
            do k=1 - block_c%signature%gc(5), block_c%signature%nk + block_c%signature%gc(6)
               do j=1 - block_c%signature%gc(3), 0
                  do i=1 - block_c%signature%gc(1), block_c%signature%ni + block_c%signature%gc(2)
                     ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         case default
            ! bc is defined into a custom INI file
            call self%set_parametric_boundary_conditions_custom(file_name=trim(block_c%signature%faces_bc(3)))
         endselect

         ! j right frame
         buffer_s = trim(block_c%signature%faces_bc(4))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         select case(buffer_c(1:3))
         case('ADJ')
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=1, block_c%signature%nk
                  do j=block_c%signature%nj + 1, block_c%signature%nj + block_c%signature%gc(4)
                     do i=1, block_c%signature%ni
                        ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = i
                        block_c%cell(i,j,k)%bc%adj(3) = j - block_c%signature%nj
                        block_c%cell(i,j,k)%bc%adj(4) = k
                     enddo
                  enddo
               enddo
            endassociate
         case('FRE', 'WAL', 'PER', 'EXT')
            do k=1 - block_c%signature%gc(5), block_c%signature%nk + block_c%signature%gc(6)
               do j=block_c%signature%nj + 1, block_c%signature%nj + block_c%signature%gc(4)
                  do i=1 - block_c%signature%gc(1), block_c%signature%ni + block_c%signature%gc(2)
                     ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         case default
            ! bc is defined into a custom INI file
            call self%set_parametric_boundary_conditions_custom(file_name=trim(block_c%signature%faces_bc(4)))
         endselect

         ! k left frame
         buffer_s = trim(block_c%signature%faces_bc(5))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         select case(buffer_c(1:3))
         case('ADJ')
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=1 - block_c%signature%gc(5), 0
                  do j=1, block_c%signature%nj
                     do i=1, block_c%signature%ni
                        ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = i
                        block_c%cell(i,j,k)%bc%adj(3) = j
                        block_c%cell(i,j,k)%bc%adj(4) = block_a%signature%nk + k
                     enddo
                  enddo
               enddo
            endassociate
         case('FRE', 'WAL', 'PER', 'EXT')
            do k=1 - block_c%signature%gc(5), 0
               do j=1, block_c%signature%nj
                  do i=1, block_c%signature%ni
                     ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         case default
            ! bc is defined into a custom INI file
            call self%set_parametric_boundary_conditions_custom(file_name=trim(block_c%signature%faces_bc(5)))
         endselect

         ! k right frame
         buffer_s = trim(block_c%signature%faces_bc(6))
         buffer_s = buffer_s%upper()
         buffer_c = buffer_s%chars()
         select case(buffer_c(1:3))
         case('ADJ')
            ba = adjacent_block_index(bc_code=buffer_c)
            associate(block_a=>self%blocks(ba)) ! adjacent block
               do k=block_c%signature%nk + 1, block_c%signature%nk + block_c%signature%gc(6)
                  do j=1, block_c%signature%nj
                     do i=1, block_c%signature%ni
                        ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                        call block_c%cell(i,j,k)%bc%initialize(code='ADJACENT')
                        block_c%cell(i,j,k)%bc%adj(1) = ba
                        block_c%cell(i,j,k)%bc%adj(2) = i
                        block_c%cell(i,j,k)%bc%adj(3) = j
                        block_c%cell(i,j,k)%bc%adj(4) = k - block_c%signature%nk
                     enddo
                  enddo
               enddo
            endassociate
         case('FRE', 'WAL', 'PER', 'EXT')
            do k=block_c%signature%nk + 1, block_c%signature%nk + block_c%signature%gc(6)
               do j=1, block_c%signature%nj
                  do i=1, block_c%signature%ni
                     ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                     call block_c%cell(i,j,k)%bc%initialize(code=buffer_c)
                  enddo
               enddo
            enddo
         case default
            ! bc is defined into a custom INI file
            call self%set_parametric_boundary_conditions_custom(file_name=trim(block_c%signature%faces_bc(6)))
         endselect
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

   lhs%error = rhs%error
   if (lhs%grid_dimensions%blocks_number/=rhs%grid_dimensions%blocks_number) then
      call lhs%destroy
      lhs%grid_dimensions = rhs%grid_dimensions
      call lhs%allocate_blocks
   endif
   do b=1, lhs%grid_dimensions%blocks_number
      lhs%blocks(b) = rhs%blocks(b)
   enddo
   endsubroutine mesh_assign_mesh

   subroutine set_parametric_boundary_conditions_custom(self, file_name)
   !< Set boundary conditions from parametric custom file.
   !<
   !< Blocks EOS must be already set entering this procedure.
   class(mesh_object), intent(inout) :: self                   !< Mesh.
   character(*),       intent(in)    :: file_name              !< File name of custom parametric bc.
   character(99)                     :: code_c                 !< BC code character buffer.
   type(string)                      :: code_s                 !< BC code string buffer.
   character(99)                     :: frame                  !< Block frame to be set.
   type(file_ini)                    :: fini                   !< Parametric file handler.
   type(primitive_compressible)      :: P                      !< Primitive variables.
   real(R8P)                         :: velocity(3)            !< Velocity temporary array.
   integer(I4P)                      :: b, i, j, k             !< Counter.
   integer(I4P)                      :: i1, i2, j1, j2, k1, k2 !< Counter.

   call fini%initialize(filename=file_name)
   call fini%load(error=self%error%status)
   if (self%error%status==0) then
      call fini%get(section_name='boundary_conditions', option_name='type', val=code_c, error=self%error%status)
      call self%error%check(message='failed to load [boundary_conditions].(type)', is_severe=.true.)
      code_s = trim(code_c)
      code_s = code_s%upper()
      call fini%get(section_name='boundary_conditions', option_name='block', val=b, error=self%error%status)
      call self%error%check(message='failed to load [boundary_conditions].(block)', is_severe=.true.)
      call fini%get(section_name='boundary_conditions', option_name='frame', val=frame, error=self%error%status)
      call self%error%check(message='failed to load [boundary_conditions].(frame)', is_severe=.true.)
      call fini%get(section_name='boundary_conditions', option_name='pressure', val=P%pressure, error=self%error%status)
      call self%error%check(message='failed to load [boundary_conditions].(pressure)', is_severe=.false.)
      call fini%get(section_name='boundary_conditions', option_name='density', val=P%density, error=self%error%status)
      call self%error%check(message='failed to load [boundary_conditions].(density)', is_severe=.false.)
      call fini%get(section_name='boundary_conditions', option_name='velocity', val=velocity, error=self%error%status)
      call self%error%check(message='failed to load [boundary_conditions].(velocity)', is_severe=.false.)
      P%velocity%x = velocity(1)
      P%velocity%y = velocity(2)
      P%velocity%z = velocity(3)
      associate(block_c=>self%blocks(b)) ! current block
         select case(trim(adjustl(frame)))
         case('i_left')
            i1 = 1 - block_c%signature%gc(1) ; i2 = 0
            j1 = 1                           ; j2 = block_c%signature%nj
            k1 = 1                           ; k2 = block_c%signature%nk
         case('i_right')
            i1 = block_c%signature%ni + 1    ; i2 = block_c%signature%ni + block_c%signature%gc(2)
            j1 = 1                           ; j2 = block_c%signature%nj
            k1 = 1                           ; k2 = block_c%signature%nk
         case('j_left')
            i1 = 1                           ; i2 = block_c%signature%ni
            j1 = 1 - block_c%signature%gc(3) ; j2 = 0
            k1 = 1                           ; k2 = block_c%signature%nk
         case('j_right')
            i1 = 1                           ; i2 = block_c%signature%ni
            j1 = block_c%signature%nj + 1    ; j2 = block_c%signature%nj + block_c%signature%gc(4)
            k1 = 1                           ; k2 = block_c%signature%nk
         case('k_left')
            i1 = 1                           ; i2 = block_c%signature%ni
            j1 = 1                           ; j2 = block_c%signature%nj
            k1 = 1 - block_c%signature%gc(5) ; k2 = 0
         case('k_right')
            i1 = 1                           ; i2 = block_c%signature%ni
            j1 = 1                           ; j2 = block_c%signature%nj
            k1 = block_c%signature%nk + 1    ; k2 = block_c%signature%nk + block_c%signature%gc(6)
         endselect
         do k=k1, k2
            do j=j1, j2
               do i=i1, i2
                  ! if (.not.allocated(block_c%cell(i,j,k)%bc)) allocate(block_c%cell(i,j,k)%bc)
                  call block_c%cell(i,j,k)%bc%initialize(code=code_s%chars())
                  block_c%cell(i, j, k)%bc%U = primitive_to_conservative_compressible(primitive=P, eos=block_c%eos)
               enddo
            enddo
         enddo
      endassociate
   else
      error stop 'boundary conditions file "'//trim(adjustl(file_name))//'" not found!'
   endif
   endsubroutine set_parametric_boundary_conditions_custom
endmodule off_mesh_object
