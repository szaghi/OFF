#include "preprocessor_macros.h"
!< OFF block object definition and implementation.

module off_block_object
!< OFF block object definition and implementation.
!<
!< [[block_object]] is a Finite Volume block-structured class.
!<
!< It allows the easy handling of metrics data for the robust and efficient computation of numerical spatial
!< operators in the framework of Finite Volume Methods (FVM).
!<
!< Let us assume that the fluid domain \(D\) is decomposed in \(N_b\) structured blocks \(D^b\), each subdivided in
!< \(N_i \times N_j \times N_k\) disjoint hexahedrons \(D_{ijk}^b\) such that \(\bigcup D_{ijk}^b = D^b\).
!< The block class is designed to aid the computations of the spatial operators into each block:
!< $$
!<\frac{\partial}{{\partial t}}\int\limits_{V_{ijk}} {\overrightarrow U dV}  =
!<-\sum\limits_{s = 1}^6 {\int\limits_{S_s} {\left(\overline{\overline {F}}\right) \cdot \overrightarrow n dS}} +
!< \int\limits_{V_{ijk}} {\overrightarrow {{Q}} dV}\label{eq:rans-cons-num}
!< $$
!< where \(S_s\) is the \(s^{th}\) face of the finite volume \(D_{ijk}\) whose measure is \(V_{ijk}\).
!<
!< A structured block is composed of hexahedron finite volumes with quadrilateral faces using the
!< following internal numeration for nodes and faces:
!<```
!< /|\Z
!<  |                            F(4)         _ F(6)
!<  |                            /|\          /!
!<  |                        7    |          /    8
!<  |                         *------------------*
!<  |                        /|   |        /    /|
!<  |                       / |   |       /    / |
!<  |                      /  |   |      /    /  |
!<  |                     /   |   |     /    /   |
!<  |                    /    |   |    +    /    |
!<  |                   /     |   |        /     |
!<  |                  /      |   +       /      |
!<  |                 /      3|          /       |4
!<  |                /        * --------/--------*
!<  |      F(1)<----/----+   /         /        /
!<  |              *------------------*    +-------->F(2)
!<  |             5|       /          |6      /
!<  |              |      /           |      /
!<  |              |     /        +   |     /
!<  |              |    /         |   |    /
!<  |              |   /      +   |   |   /
!<  |              |  /      /    |   |  /
!<  |              | /      /     |   | /
!<  |              |/      /      |   |/
!<  |              *------------------*
!<  |             1      /        |    2
!<  |                   /        \|/
!<  |   _ Y           |/_       F(3)
!<  |   /|         F(5)
!<  |  /
!<  | /
!<  |/                                                    X
!<  O----------------------------------------------------->
!<```
!< Each hexadron cells is faces-connected to its neighboring, thus the cells build a structured block with implicit
!< connectivity, e.g. in 2D space a block could be as the following:
!<```
!<                 _ J
!<                 /|                          _____
!<               5+ ...*----*----*----*----*...     |
!<               /    /    /    /    /    /         |
!<              /    /    /    /    /    /          |
!<            4+ ...*----*----*----*----*...        |
!<            /    /    /    /    /    /            |
!<           /    /    /    /    /    /             |
!<         3+ ...*----*----*----*----*...           |  Structured block of 4x4 Finite Volumes
!<         /    /    / FV /    /    /               |
!<        /    /    /    /    /    /                |
!<      2+ ...*----*----*----*----*...              |
!<      /    /    /    /    /    /                  |
!<     /    /    /    /    /    /                   |
!<   1+ ...*----*----*----*----*...                 |
!<   /     .    .    .    .    .                    |
!<  /      .    .    .    .    .               _____
!< O-------+----+----+----+----+-------------> I
!<         1    2    3    4    5
!<```
!< The nodes of cells are not required to be on the Cartesian coordinates, thus allowing a general
!< curvilinear mesh: the are 3 implicit coordinate lines, *i*, *j* and *k* that are not required to be orthogonal.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_bc_object, only : BC_ADJACENT
use off_block_signature_object, only : block_signature_object
use off_cell_object, only : cell_object
use off_error_object, only : error_object
use off_face_object, only : face_object
use off_level_set_object, only : level_set_object, level_set_normal
use off_node_object, only : node_object
use off_solver_object, only : solver_object
use flow, only : conservative_compressible, primitive_compressible,                              &
                 conservative_to_primitive_compressible, primitive_to_conservative_compressible, &
                 eos_compressible
use foreseer, only : riemann_solver_object
use fossil, only : surface_stl_object
use penf, only : FR8P, FI4P, I1P, I4P, I8P, MaxR8P, MinR8P, R8P, str
use vecfor, only : vector, ex, ey, ez, normL2, sq_norm
use vtk_fortran, only : vtk_file

implicit none
private
public :: block_object

integer(I4P), parameter :: NO_ERROR                           = 0 !< No errors occurred.
integer(I4P), parameter :: ERROR_BLOCK_COMPUTE_EXTENTS_FAILED = 1 !< Failed to compute block extents.
integer(I4P), parameter :: ERROR_BLOCK_CREATE_FAILED          = 2 !< Failed to create block.
integer(I4P), parameter :: ERROR_BLOCK_DESTROY_FAILED         = 3 !< Failed to destroy block.
integer(I4P), parameter :: ERROR_BLOCK_CREATE_LINSPACE_FAILED = 4 !< Failed to create a uniform-spaced linear block.

type :: block_object
   !< Block object class.
   type(error_object)             :: error         !< Errors handler.
   type(block_signature_object)   :: signature     !< Signature, namely id, level, dimensions, etc...
   type(cell_object), allocatable :: cell(:,:,:)   !< Cell.
   type(face_object), allocatable :: face_i(:,:,:) !< Faces along I direction.
   type(face_object), allocatable :: face_j(:,:,:) !< Faces along I direction.
   type(face_object), allocatable :: face_k(:,:,:) !< Faces along I direction.
   type(node_object), allocatable :: node(:,:,:)   !< Cell.
   ! fluid dynamic data
   type(eos_compressible) :: eos !< Equation of state.
   contains
      ! public methods
      procedure, pass(self) :: cells_number                 !< Return the number of cells.
      procedure, pass(self) :: compute_space_operator       !< Compute space operator.
      procedure, pass(self) :: compute_dt                   !< Compute the current time step by means of CFL condition.
      procedure, pass(self) :: compute_recv_maps            !< Compute querying and receiving maps.
      procedure, pass(self) :: compute_residuals            !< Compute residuals.
      procedure, pass(self) :: create_linspace              !< Create a Cartesian block with linearly spaced nodes.
      procedure, pass(self) :: description                  !< Return a pretty-formatted description of the block.
      procedure, pass(self) :: destroy                      !< Destroy block.
      procedure, pass(self) :: dt_min                       !< Return the minimum Dt into internal cells.
      procedure, pass(self) :: evolve_immersed_boundary     !< Evolve immersed boundary from domain into the immersed boundary.
      procedure, pass(self) :: immerge_stl_geometry         !< *Immerge* geometry (described by STL file) into the block grid.
      procedure, pass(self) :: interpolate_at_nodes         !< Interpolate cell-centered variable at nodes.
      procedure, pass(self) :: initialize                   !< Initialize block.
      procedure, pass(self) :: load_conservatives_from_file !< Load nodes from file.
      procedure, pass(self) :: load_nodes_from_file         !< Load nodes from file.
      procedure, pass(self) :: nodes_number                 !< Return the number of nodes.
      procedure, pass(self) :: set_eos                      !< Set EOS.
      procedure, pass(self) :: save_conservatives_into_file !< Save conservatives into file.
      procedure, pass(self) :: save_file_grid               !< Save grid file.
      procedure, pass(self) :: save_file_solution           !< Save solution file.
      procedure, pass(self) :: save_nodes_into_file         !< Save nodes into file.
      procedure, pass(self) :: update_recv_cells_number     !< Update receive cells number for multi-processes communication.
      procedure, pass(self) :: update_level_set_distance    !< Update level set distance.
      ! operators
      generic :: assignment(=) => block_assign_block !< Overload `=`.
      ! fast operators
      procedure, pass(opr) :: conservative_add_conservatives_fast      !< `+` fast operator.
      procedure, pass(opr) :: conservative_multiply_conservatives_fast !< `*` fast operator.
      procedure, pass(opr) :: conservative_multiply_real_scalar_fast   !< `* real_scalar` fast operator.
      procedure, pass(opr) :: conservative_subtract_conservatives_fast !< `-` fast operator.
      ! private methods
      procedure, pass(lhs),  private :: block_assign_block     !< Operator `=`.
      procedure, pass(self), private :: compute_cells_center   !< Compute cells center.
      procedure, pass(self), private :: compute_extents        !< Compute block extents.
      procedure, pass(self), private :: compute_faces_metrics  !< Compute block faces metrics.
      procedure, pass(self), private :: compute_metrics        !< Compute block metrics.
      procedure, pass(self), private :: compute_volumes        !< Compute block volumes.
      procedure, pass(self), private :: correct_metrics        !< Correct block metrics.
      procedure, pass(self), private :: node_to_center         !< Compute cell centers coordinates from cell nodes.
      procedure, pass(self), private :: nullify_normals        !< Nullify normals for 2D or 1D domains.
      procedure, pass(self), private :: save_file_grid_tec     !< Save mesh grid into Tecplot file.
      procedure, pass(self), private :: save_file_grid_vtk     !< Save mesh grid into VTK file.
      procedure, pass(self), private :: save_file_solution_vtk !< Save grid file in VTK format.
endtype block_object

contains
   ! public methods
   elemental function cells_number(self, with_ghosts) result(cells_number_)
   !< Return the number of cells.
   class(block_object), intent(in)           :: self          !< Block.
   logical,             intent(in), optional :: with_ghosts   !< Take into account ghost cells.
   integer(I4P)                              :: cells_number_ !< Number of cells.

   cells_number_ = self%signature%cells_number(with_ghosts=with_ghosts)
   endfunction cells_number

   pure subroutine compute_dt(self, CFL)
   !< Compute the current time step by means of CFL condition.
   class(block_object), intent(inout) :: self      !< Block.
   real(R8P),           intent(in)    :: CFL       !< CFL value.
   type(vector)                       :: u(0:6)    !< Velocity vectors.
   real(R8P)                          :: a         !< Speed of sound.
   real(R8P)                          :: umax(0:6) !< Maximum propagation speed of signals on faces and cell.
   integer(I4P)                       :: i, j, k   !< Counter.
   real(R8P)                          :: dt_max    !< Maximum block time step.

   dt_max = MinR8P
   associate(Ni=>self%signature%Ni, Nj=>self%signature%Nj, Nk=>self%signature%Nk)
      do k=1, Nk
         do j=1, Nj
            do i=1, Ni
               if (self%cell(i, j, k)%level_set%distance<0._R8P) then
                  ! cell into an immersed boundary
                  self%cell(i, j, k)%Dt = 0._R8P
               else
                  a = self%eos%speed_of_sound(density=self%cell(i, j, k)%U%density, &
                                              pressure=self%cell(i, j, k)%U%pressure(eos=self%eos))

                  ! velocity vector surrounding the cell
                  u(0) = self%cell(i  , j, k)%U%velocity()
                  u(1) = self%cell(i-1, j, k)%U%velocity()
                  u(2) = self%cell(i+1, j, k)%U%velocity()
                  u(3) = self%cell(i, j-1, k)%U%velocity()
                  u(4) = self%cell(i, j+1, k)%U%velocity()
                  u(5) = self%cell(i, j, k-1)%U%velocity()
                  u(6) = self%cell(i, j, k+1)%U%velocity()

                  ! maximum of normal velocities means
                  umax(1) = (abs((0.5_R8P * (u(1) + u(0))).dot.self%face_i(i-1, j,  k)%normal) + a) * self%face_i(i-1, j,  k)%area
                  umax(2) = (abs((0.5_R8P * (u(2) + u(0))).dot.self%face_i(i  , j,  k)%normal) + a) * self%face_i(i  , j,  k)%area
                  umax(3) = (abs((0.5_R8P * (u(3) + u(0))).dot.self%face_j(i ,j-1,  k)%normal) + a) * self%face_j(i ,j-1,  k)%area
                  umax(4) = (abs((0.5_R8P * (u(4) + u(0))).dot.self%face_j(i ,j  ,  k)%normal) + a) * self%face_j(i ,j  ,  k)%area
                  umax(5) = (abs((0.5_R8P * (u(5) + u(0))).dot.self%face_k(i ,j  ,k-1)%normal) + a) * self%face_k(i ,j  ,k-1)%area
                  umax(6) = (abs((0.5_R8P * (u(6) + u(0))).dot.self%face_k(i ,j  ,k  )%normal) + a) * self%face_k(i ,j  ,k  )%area
                  umax(0) = maxval(umax(1:6))

                  ! time step
                  self%cell(i, j, k)%Dt = self%cell(i, j, k)%volume * CFL / umax(0)
                  dt_max = max(dt_max, self%cell(i, j, k)%Dt)
               endif
            enddo
         enddo
      enddo
      ! setting immersed boundary cell to maximum time step
      do k=1, Nk
         do j=1, Nj
            do i=1, Ni
               if (self%cell(i, j, k)%level_set%distance<0._R8P) self%cell(i, j, k)%Dt = dt_max
            enddo
         enddo
      enddo
   endassociate
   endsubroutine compute_dt

   subroutine compute_recv_maps(self, b, process, block_to_process_map, c, reqs_map, recv_map)
   !< Compute querying and receiving maps.
   class(block_object), intent(in)    :: self                     !< Block.
   integer(I4P),        intent(in)    :: b                        !< Block ID in local numeration.
   integer(I4P),        intent(in)    :: process                  !< Other process than me to send/receive cells.
   integer(I4P),        intent(in)    :: block_to_process_map(1:) !< Maps of processes ID for all blocks[1:blocks_global_number].
   integer(I4P),        intent(inout) :: c                        !< Cells counter.
   integer(I4P),        intent(out)   :: reqs_map(1:,1:)          !< Querying cells map  [1:4,1:recv_cells_number].
   integer(I4P),        intent(out)   :: recv_map(1:,1:)          !< Receiving cells map [1:4,1:recv_cells_number].
   integer(I4P)                       :: i, j, k                  !< Counter.

   associate(gc=>self%signature%gc, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk, cell=>self%cell)
      do k=1, nk
         do j=1, nj
            if (cell(0   ,j,k)%bc%is(BC_ADJACENT)) call update_frame(gc=gc(1), frame=cell(1-gc(1):0    ,j,k), c=c, &
                                                                     reqs_map=reqs_map, recv_map=recv_map)
            if (cell(ni+1,j,k)%bc%is(BC_ADJACENT)) call update_frame(gc=gc(2), frame=cell(ni+1:ni+gc(2),j,k), c=c, &
                                                                     reqs_map=reqs_map, recv_map=recv_map)
         enddo
      enddo
      do k=1, nk
         do i=1, ni
            if (cell(i,0   ,k)%bc%is(BC_ADJACENT)) call update_frame(gc=gc(3), frame=cell(i,1-gc(3):0    ,k), c=c, &
                                                                     reqs_map=reqs_map, recv_map=recv_map)
            if (cell(i,nj+1,k)%bc%is(BC_ADJACENT)) call update_frame(gc=gc(4), frame=cell(i,nj+1:nj+gc(4),k), c=c, &
                                                                     reqs_map=reqs_map, recv_map=recv_map)
         enddo
      enddo
      do j=1, nj
         do i=1, ni
            if (cell(i,j,0   )%bc%is(BC_ADJACENT)) call update_frame(gc=gc(5), frame=cell(i,j,1-gc(5):0    ), c=c, &
                                                                     reqs_map=reqs_map, recv_map=recv_map)
            if (cell(i,j,nk+1)%bc%is(BC_ADJACENT)) call update_frame(gc=gc(6), frame=cell(i,j,nk+1:nk+gc(6)), c=c, &
                                                                     reqs_map=reqs_map, recv_map=recv_map)
         enddo
      enddo
   endassociate
   contains
      _PURE_ subroutine update_frame(gc, frame, c, reqs_map, recv_map)
      !< Update maps for current frame.
      integer(I4P),      intent(in)    :: gc              !< Number of ghost cells.
      type(cell_object), intent(in)    :: frame(1-gc:)    !< Cells frame [1-gc:0].
      integer(I4P),      intent(inout) :: c               !< Cells counter.
      integer(I4P),      intent(out)   :: reqs_map(1:,1:) !< Querying cells map  [1:4,1:recv_cells_number].
      integer(I4P),      intent(out)   :: recv_map(1:,1:) !< Receiving cells map [1:4,1:recv_cells_number].
      integer(I4P)                     :: ii              !< Counter.
      integer(I4P)                     :: proc            !< Processor ID.

      do ii=1-gc, 0
         associate(adj_b=>frame(ii)%bc%adj(1), adj_i=>frame(ii)%bc%adj(2), adj_j=>frame(ii)%bc%adj(3), adj_k=>frame(ii)%bc%adj(4))
         proc = block_to_process_map(adj_b)
         if ( proc == process) then
            c = c + 1
            reqs_map(1, c) = adj_b ; recv_map(1, c) = b
            reqs_map(2, c) = adj_i ; recv_map(2, c) = i
            reqs_map(3, c) = adj_j ; recv_map(3, c) = j
            reqs_map(4, c) = adj_k ; recv_map(4, c) = k
         endif
         endassociate
      enddo
      endsubroutine update_frame
  endsubroutine compute_recv_maps

   subroutine compute_residuals(self, solver, gcu)
   !< Compute residuals.
   class(block_object),  intent(inout)          :: self                !< Block.
   class(solver_object), intent(in)             :: solver              !< Solver.
   integer(I4P),         intent(in)             :: gcu                 !< Number of ghost cells used (depend on space accuracy).
   type(conservative_compressible), allocatable :: fluxes_con_i(:,:,:) !< Convective fluxes along i-direction.
   type(conservative_compressible), allocatable :: fluxes_con_j(:,:,:) !< Convective fluxes along j-direction.
   type(conservative_compressible), allocatable :: fluxes_con_k(:,:,:) !< Convective fluxes along k-direction.
   integer(I4P)                                 :: i, j, k, gcm        !< Counter.

   associate(gc=>self%signature%gc, Ni=>self%signature%Ni, Nj=>self%signature%Nj, Nk=>self%signature%Nk, &
             cell=>self%cell, face_i=>self%face_i, face_j=>self%face_j, face_k=>self%face_k)
      allocate(fluxes_con_i(0:Ni,1:Nj,1:Nk))
      allocate(fluxes_con_j(1:Ni,0:Nj,1:Nk))
      allocate(fluxes_con_k(1:Ni,1:Nj,0:Nk))
      gcm = min(gcu, gc(1), gc(2))
      do k=1, Nk
         do j=1, Nj
            call compute_fluxes_convective(solver = solver,                           &
                                           eos    = self%eos,                         &
                                           gc     = gcm,                              &
                                           N      = Ni,                               &
                                           faces  = self%face_i (0-gcm:Ni+gcm, j, k), &
                                           cells  = self%cell   (1-gcm:Ni+gcm, j, k), &
                                           fluxes = fluxes_con_i(    0:Ni    , j, k))
         enddo
      enddo
      gcm = min(gcu, gc(3), gc(4))
      do k=1, Nk
         do i=1, Ni
            call compute_fluxes_convective(solver = solver,                           &
                                           eos    = self%eos,                         &
                                           gc     = gcm,                              &
                                           N      = Nj,                               &
                                           faces  = self%face_j (i, 0-gcm:Nj+gcm, k), &
                                           cells  = self%cell   (i, 1-gcm:Nj+gcm, k), &
                                           fluxes = fluxes_con_j(i,     0:Nj    , k))
         enddo
      enddo
      gcm = min(gcu, gc(5), gc(6))
      do j=1, Nj
         do i=1, Ni
            call compute_fluxes_convective(solver = solver,                           &
                                           eos    = self%eos,                         &
                                           gc     = gcm,                              &
                                           N      = Nk,                               &
                                           faces  = self%face_k (i, j, 0-gcm:Nk+gcm), &
                                           cells  = self%cell   (i, j, 1-gcm:Nk+gcm), &
                                           fluxes = fluxes_con_k(i, j,     0:Nk    ))
         enddo
      enddo
      do k=1, Nk
         do j=1, Nj
            do i=1, Ni

               if (cell(i, j, k)%level_set%distance<0._R8P) then
                  cell(i,j,k)%U%density = 0._R8P
                  cell(i,j,k)%U%momentum = 0._R8P
                  cell(i,j,k)%U%energy = 0._R8P
                  cycle
               endif

               cell(i,j,k)%U%density =                                          &
               (face_i(i-1,j,  k  )%area * fluxes_con_i(i-1,j,  k  )%density -  &
                face_i(i,  j,  k  )%area * fluxes_con_i(i,  j,  k  )%density +  &
                face_j(i,  j-1,k  )%area * fluxes_con_j(i,  j-1,k  )%density -  &
                face_j(i,  j,  k  )%area * fluxes_con_j(i,  j,  k  )%density +  &
                face_k(i,  j,  k-1)%area * fluxes_con_k(i,  j,  k-1)%density -  &
                face_k(i,  j,  k  )%area * fluxes_con_k(i,  j,  k  )%density) / &
               cell(i,j,k)%volume

               cell(i,j,k)%U%momentum =                                          &
               (face_i(i-1,j,  k  )%area * fluxes_con_i(i-1,j,  k  )%momentum -  &
                face_i(i,  j,  k  )%area * fluxes_con_i(i,  j,  k  )%momentum +  &
                face_j(i,  j-1,k  )%area * fluxes_con_j(i,  j-1,k  )%momentum -  &
                face_j(i,  j,  k  )%area * fluxes_con_j(i,  j,  k  )%momentum +  &
                face_k(i,  j,  k-1)%area * fluxes_con_k(i,  j,  k-1)%momentum -  &
                face_k(i,  j,  k  )%area * fluxes_con_k(i,  j,  k  )%momentum) / &
               cell(i,j,k)%volume

               cell(i,j,k)%U%energy =                                          &
               (face_i(i-1,j,  k  )%area * fluxes_con_i(i-1,j,  k  )%energy -  &
                face_i(i,  j,  k  )%area * fluxes_con_i(i,  j,  k  )%energy +  &
                face_j(i,  j-1,k  )%area * fluxes_con_j(i,  j-1,k  )%energy -  &
                face_j(i,  j,  k  )%area * fluxes_con_j(i,  j,  k  )%energy +  &
                face_k(i,  j,  k-1)%area * fluxes_con_k(i,  j,  k-1)%energy -  &
                face_k(i,  j,  k  )%area * fluxes_con_k(i,  j,  k  )%energy) / &
               cell(i,j,k)%volume
            enddo
         enddo
      enddo
   endassociate
   endsubroutine compute_residuals

   subroutine compute_space_operator(self)!, space_operator)
   !< Compute space operator.
   !<
   !< @TODO implement space operator.
   !<
   !< @TODO re-add elemental attribute.
   class(block_object),        intent(in)    :: self                     !< Block.
   ! class(conservative_object), intent(inout) :: space_operator(1:,1:,1:) !< Space operator.

   error stop 'error: block space operator to be implemented'
   endsubroutine compute_space_operator

   subroutine create_linspace(self, emin, emax)
   !< Create a Cartesian block with linearly spaced nodes.
   !<
   !< @note If the extents (emin, emax) of the block are not passed, the values already saved into the block are used.
   !<
   !< @TODO re-add elemental attribute.
   class(block_object), intent(inout)        :: self    !< Block.
   type(vector),        intent(in), optional :: emin    !< Coordinates of minimum abscissa of the block.
   type(vector),        intent(in), optional :: emax    !< Coordinates of maximum abscissa of the block.
   type(vector)                              :: delta   !< Diagonal of block bounding-box.
   real(R8P)                                 :: delta_x !< X component of diagonal of block bounding-box.
   real(R8P)                                 :: delta_y !< Y component of diagonal of block bounding-box.
   real(R8P)                                 :: delta_z !< Z component of diagonal of block bounding-box.
   integer(I4P)                              :: i       !< Counter.
   integer(I4P)                              :: j       !< Counter.
   integer(I4P)                              :: k       !< Counter.

   self%error%status = ERROR_BLOCK_CREATE_LINSPACE_FAILED

   self%signature%is_cartesian = .true.
   if (present(emin)) self%signature%emin = emin
   if (present(emax)) self%signature%emax = emax
   associate(gc=>self%signature%gc, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk)
      if (self%signature%emin/=self%signature%emax) then
         delta = (self%signature%emax - self%signature%emin) / (ni * ex + nj * ey + nk * ez)
         delta_x = delta.dot.ex
         delta_y = delta.dot.ey
         delta_z = delta.dot.ez
         do k=0 - gc(5), nk + gc(6)
            do j=0 - gc(3), nj + gc(4)
               do i=0 - gc(1), ni + gc(2)
                  self%node(i, j, k)%vertex = self%signature%emin + (i * delta_x) * ex + (j * delta_y) * ey + (k * delta_z) * ez
               enddo
            enddo
         enddo
         call self%compute_metrics
         self%error%status = NO_ERROR
      endif
   endassociate
   endsubroutine create_linspace

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the block.
   class(block_object), intent(in)           :: self             !< Block.
   character(*),        intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable             :: desc             !< Description.
   character(len=:), allocatable             :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter               :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
         desc = desc//prefix_//'error status      : '//trim(str(n=self%error%status))//NL
   if (allocated(self%cell)) then
         desc = desc//prefix_//'max density       : '//trim(str(n=maxval(self%cell%U%density)))//NL
         desc = desc//prefix_//'max energy        : '//trim(str(n=maxval(self%cell%U%energy)))//NL
         desc = desc//prefix_//'max momentum      : '//trim(str(n=maxval(self%cell%U%momentum%x)))//' '// &
                                                       trim(str(n=maxval(self%cell%U%momentum%y)))//' '// &
                                                       trim(str(n=maxval(self%cell%U%momentum%z)))//NL
      ! if (allocated(self%cell(0,1,1)%bc).and.allocated(self%cell(self%signature%Ni+1,1,1)%bc).and. &
      !     allocated(self%cell(1,0,1)%bc).and.allocated(self%cell(1,self%signature%Nj+1,1)%bc).and. &
      !     allocated(self%cell(1,1,0)%bc).and.allocated(self%cell(1,1,self%signature%Nk+1)%bc)) then
         desc = desc//prefix_//'faces-i bc (L,R)  : '//trim(str(self%cell(0,1,1)%bc%id))//', '//&
                                                       trim(str(self%cell(self%signature%Ni+1,1,1)%bc%id))//NL
         desc = desc//prefix_//'faces-j bc (L,R)  : '//trim(str(self%cell(1,0,1)%bc%id))//', '//&
                                                       trim(str(self%cell(1,self%signature%Nj+1,1)%bc%id))//NL
         desc = desc//prefix_//'faces-k bc (L,R)  : '//trim(str(self%cell(1,1,0)%bc%id))//', '//&
                                                       trim(str(self%cell(1,1,self%signature%Nk+1)%bc%id))//NL
      ! endif
   endif
         desc = desc//prefix_//'equations of state: '//NL//self%eos%description(prefix=prefix_//'  ')
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy block.
   class(block_object), intent(inout) :: self      !< Block.
   type(eos_compressible)             :: eos_fresh !< Fresh EOS state.

   call self%error%destroy
   call self%signature%destroy
   if (allocated(self%cell)) then
      call self%cell%destroy
      deallocate(self%cell)
   endif
   if (allocated(self%face_i)) then
      call self%face_i%destroy
      deallocate(self%face_i)
   endif
   if (allocated(self%face_j)) then
      call self%face_j%destroy
      deallocate(self%face_j)
   endif
   if (allocated(self%face_k)) then
      call self%face_k%destroy
      deallocate(self%face_k)
   endif
   if (allocated(self%node)) then
      call self%node%destroy
      deallocate(self%node)
   endif
   self%eos = eos_fresh
   endsubroutine destroy

   elemental function dt_min(self)
   !< Return the minimum Dt into internal cells.
   class(block_object), intent(in) :: self    !< Block.
   real(R8P)                       :: dt_min  !< Minimum Dt into internal cells.
   integer(I4P)                    :: i, j, k !< Counter.

   dt_min = MaxR8P
   associate(Ni=>self%signature%Ni, Nj=>self%signature%Nj, Nk=>self%signature%Nk)
      do k=1, Nk
         do j=1, Nj
            do i=1, Ni
               dt_min = min(dt_min, self%cell(i,j,k)%Dt)
            enddo
         enddo
      enddo
   endassociate
   endfunction dt_min

   _PURE_ subroutine evolve_immersed_boundary(self)
   !< Evolve immersed boundary from domain into the immersed boundary.
   class(block_object), intent(inout)           :: self       !< Block.
   type(conservative_compressible), allocatable :: U(:,:,:)   !< Conservative variables temporary buffer.
   type(conservative_compressible)              :: Ux, Uy, Uz !< Conservative variables variations.
   real(R8P)                                    :: dt         !< Pseudo-time step.
   integer(I4P)                                 :: i, j, k    !< Counter.

   associate(ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk, cell=>self%cell, &
             is_null_x=>self%signature%is_null_x, is_null_y=>self%signature%is_null_y, is_null_z=>self%signature%is_null_z)
      allocate(U(1:ni, 1:nj, 1:nk))
      do k=1, nk
         do j=1, nj
            do i=1, ni
               if (cell(i,j,k)%level_set%distance<0._R8P) then
                  if (is_null_x) then
                     Ux = 0._R8P * Ux
                  else
                     if (cell(i,j,k)%level_set%normal%x<0._R8P) then
                        Ux = cell(i,j,k)%U - cell(i-1,j,k)%U
                     else
                        Ux = cell(i+1,j,k)%U - cell(i,j,k)%U
                     endif
                  endif
                  if (is_null_y) then
                     Uy = 0._R8P * Uy
                  else
                     if (cell(i,j,k)%level_set%normal%y<0._R8P) then
                        Uy = cell(i,j,k)%U - cell(i-1,j,k)%U
                     else
                        Uy = cell(i,j+1,k)%U - cell(i,j,k)%U
                     endif
                  endif
                  if (is_null_z) then
                     Uz = 0._R8P * Uz
                  else
                     if (cell(i,j,k)%level_set%normal%z<0._R8P) then
                        Uz = cell(i,j,k)%U - cell(i-1,j,k)%U
                     else
                        Uz = cell(i,j,k+1)%U - cell(i,j,k)%U
                     endif
                  endif
                  dt = (1._R8P / (abs(cell(i,j,k)%level_set%normal%x) + &
                                  abs(cell(i,j,k)%level_set%normal%y) + &
                                  abs(cell(i,j,k)%level_set%normal%z))) * 0.3_R8P
                  U(i,j,k) = cell(i,j,k)%U + dt * (Ux * cell(i,j,k)%level_set%normal%x + &
                                                   Uy * cell(i,j,k)%level_set%normal%y + &
                                                   Uz * cell(i,j,k)%level_set%normal%z)
               endif
            enddo
         enddo
      enddo
      do k=1, nk
         do j=1, nj
            do i=1, ni
               if (cell(i,j,k)%level_set%distance<0._R8P) cell(i,j,k)%U = U(i,j,k)
            enddo
         enddo
      enddo
   endassociate
   endsubroutine evolve_immersed_boundary

   subroutine immerge_stl_geometry(self, surface_stl, n, distance_sign_inverse, distance_sign_algorithm)
   !< *Immerge* geometry (described into a STL file) into the block grid.
   class(block_object),      intent(inout)        :: self                     !< Block.
   type(surface_stl_object), intent(in)           :: surface_stl              !< STL surface handler.
   integer(I4P),             intent(in)           :: n                        !< Number of geometry in the global numbering.
   logical,                  intent(in), optional :: distance_sign_inverse    !< Invert sign distance.
   character(*),             intent(in), optional :: distance_sign_algorithm  !< Algorithm used for *sign* of distance computation.
   integer(I4P)                                   :: distance_sign            !< Distance sing.
   character(len=:), allocatable                  :: distance_sign_algorithm_ !< Algorithm used for *sign* of distance computation.
   integer(I4P)                                   :: i, j, k                  !< Counter.

   if (allocated(self%cell)) then
      distance_sign = 1
      if (present(distance_sign_inverse)) then
         if (distance_sign_inverse) distance_sign = -1
      endif
      distance_sign_algorithm_ = 'ray_intersections'
      if (present(distance_sign_algorithm)) distance_sign_algorithm_ = distance_sign_algorithm

      if ((self%signature%emin%x >= surface_stl%bmax%x).or.(self%signature%emax%x <= surface_stl%bmin%x).or.&
          (self%signature%emin%y >= surface_stl%bmax%y).or.(self%signature%emax%y <= surface_stl%bmin%y).or.&
          (self%signature%emin%z >= surface_stl%bmax%z).or.(self%signature%emax%z <= surface_stl%bmin%z)) then
         return
      else
         do k=1, self%signature%nk
            do j=1, self%signature%nj
               do i=1, self%signature%ni
                  self%cell(i,j,k)%level_set%distances(n) = distance_sign * &
                                                            surface_stl%distance(point=self%cell(i,j,k)%center, &
                                                                                 is_signed=.true.,              &
                                                                                 sign_algorithm=trim(distance_sign_algorithm_))
               enddo
            enddo
         enddo
      endif
   endif
   endsubroutine immerge_stl_geometry

   _PURE_ subroutine initialize(self, signature,                                           &
                                id, level, gc, ni, nj, nk,                                 &
                                emin, emax, is_cartesian, is_null_x, is_null_y, is_null_z, &
                                interfaces_number, distances, eos, U0)
   !< Initialize block.
   !<
   !< Assign block signature, allocate dynamic memory and set block features.
   class(block_object),             intent(inout)        :: self                    !< Block.
   type(block_signature_object),    intent(in), optional :: signature               !< Signature, namely id, level, dimensions, etc.
   integer(I8P),                    intent(in), optional :: id                      !< Unique (Morton) identification code.
   integer(I4P),                    intent(in), optional :: level                   !< Grid refinement level.
   integer(I4P),                    intent(in), optional :: gc(1:)                  !< Number of ghost cells along each frame.
   integer(I4P),                    intent(in), optional :: ni                      !< Number of cells in I direction.
   integer(I4P),                    intent(in), optional :: nj                      !< Number of cells in J direction.
   integer(I4P),                    intent(in), optional :: nk                      !< Number of cells in K direction.
   type(vector),                    intent(in), optional :: emin                    !< Coordinates of minimum abscissa of the block.
   type(vector),                    intent(in), optional :: emax                    !< Coordinates of maximum abscissa of the block.
   logical,                         intent(in), optional :: is_cartesian            !< Flag for checking if the block is Cartesian.
   logical,                         intent(in), optional :: is_null_x               !< Nullify X direction (2D yz, 1D y/z domain).
   logical,                         intent(in), optional :: is_null_y               !< Nullify Y direction (2D xy, 1D x/y domain).
   logical,                         intent(in), optional :: is_null_z               !< Nullify Z direction (2D xy, 1D x/y domain).
   integer(I4P),                    intent(in), optional :: interfaces_number       !< Number of different interfaces.
   real(R8P),                       intent(in), optional :: distances(:)            !< Distance from all interfaces.
   type(eos_compressible),          intent(in), optional :: eos                     !< EOS data.
   type(conservative_compressible), intent(in), optional :: U0                      !< Initial state of conservative variables.
   integer(I4P)                                          :: i                       !< Counter.
   integer(I4P)                                          :: j                       !< Counter.
   integer(I4P)                                          :: k                       !< Counter.

   self%error%status = ERROR_BLOCK_CREATE_FAILED

   call self%destroy

   call self%signature%initialize(signature=signature,                             &
                                  id=id, level=level, gc=gc, ni=ni, nj=nj, nk=nk,  &
                                  interfaces_number=interfaces_number,             &
                                  U0=U0,                                           &
                                  emin=emin, emax=emax, is_cartesian=is_cartesian, &
                                  is_null_x=is_null_x, is_null_y=is_null_y, is_null_z=is_null_z)

   associate(gc_=>self%signature%gc, ni_=>self%signature%ni, nj_=>self%signature%nj, nk_=>self%signature%nk)
      allocate(self%cell(  1 - gc_(1) : ni_ + gc_(2), 1 - gc_(3) : nj_ + gc_(4), 1 - gc_(5) : nk_ + gc_(6)))
      allocate(self%face_i(0 - gc_(1) : ni_ + gc_(2), 1 - gc_(3) : nj_ + gc_(4), 1 - gc_(5) : nk_ + gc_(6)))
      allocate(self%face_j(1 - gc_(1) : ni_ + gc_(2), 0 - gc_(3) : nj_ + gc_(4), 1 - gc_(5) : nk_ + gc_(6)))
      allocate(self%face_k(1 - gc_(1) : ni_ + gc_(2), 1 - gc_(3) : nj_ + gc_(4), 0 - gc_(5) : nk_ + gc_(6)))
      allocate(self%node(  0 - gc_(1) : ni_ + gc_(2), 0 - gc_(3) : nj_ + gc_(4), 0 - gc_(5) : nk_ + gc_(6)))
   endassociate

   associate(gc=>self%signature%gc, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk)
      do k=1 - gc(5), nk + gc(6)
         do j=1 - gc(3), nj + gc(4)
            do i=1 - gc(1), ni + gc(2)
               call self%cell(i,j,k)%initialize(interfaces_number=self%signature%interfaces_number, distances=distances, U=U0)
            enddo
         enddo
      enddo
   endassociate

   call self%face_i%initialize
   call self%face_j%initialize
   call self%face_k%initialize

   if (present(eos)) self%eos = eos

   self%error%status = NO_ERROR
   endsubroutine initialize

   pure subroutine interpolate_at_nodes(self, var_cell, var_node)
   !< Interpolate cell-centered variable at nodes.
   !<
   !< @note The interpolation is linear and based on the volume-weights.
   !<
   !< @note Only internal cells are considered, ghost ones are trimmed.
   class(block_object), intent(in)  :: self                                !< Block.
   real(R8P),           intent(in)  :: var_cell(1-self%signature%gc(1):, &
                                                1-self%signature%gc(3):, &
                                                1-self%signature%gc(5):)   !< Cell-centered variable.
   real(R8P),           intent(out) :: var_node(0-self%signature%gc(1):, &
                                                0-self%signature%gc(3):, &
                                                0-self%signature%gc(5):)   !< Node-centered variable.
   real(R8P), allocatable           :: var_cell_framed(:,:,:)              !< Cell-centered var framed.
   real(R8P), allocatable           :: volume_framed(:,:,:)                !< Volume framed.
   integer(I4P)                     :: i                                   !< Counter.
   integer(I4P)                     :: j                                   !< Counter.
   integer(I4P)                     :: k                                   !< Counter.

   associate(gc=>self%signature%gc, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk)
      ! building framed variable and volume
      allocate(var_cell_framed(0:ni+1, 0:nj+1, 0:nk+1)) ; var_cell_framed = 0._R8P
      var_cell_framed(1:ni, 1:nj, 1:nk) = var_cell(1:ni, 1:nj, 1:nk)
      allocate(volume_framed(0:ni+1, 0:nj+1, 0:nk+1)) ; volume_framed = 0._R8P
      volume_framed(1:ni, 1:nj, 1:nk) = self%cell(1:ni, 1:nj, 1:nk)%volume
      ! check frames
      if (gc(1)>0) then
         var_cell_framed(0, 1:nj, 1:nk) = var_cell( 0, 1:nj, 1:nk)
         volume_framed(  0, 1:nj, 1:nk) = self%cell(0, 1:nj, 1:nk)%volume
      endif
      if (gc(2)>0) then
         var_cell_framed(ni+1, 1:nj, 1:nk) = var_cell( ni+1, 1:nj, 1:nk)
         volume_framed(  ni+1, 1:nj, 1:nk) = self%cell(ni+1, 1:nj, 1:nk)%volume
      endif
      if (gc(3)>0) then
         var_cell_framed(1:ni, 0, 1:nk) = var_cell( 1:ni, 0, 1:nk)
         volume_framed(  1:ni, 0, 1:nk) = self%cell(1:ni, 0, 1:nk)%volume
      endif
      if (gc(4)>0) then
         var_cell_framed(ni, 1:nj+1, 1:nk) = var_cell( ni, 1:nj+1, 1:nk)
         volume_framed(  ni, 1:nj+1, 1:nk) = self%cell(ni, 1:nj+1, 1:nk)%volume
      endif
      if (gc(5)>0) then
         var_cell_framed(1:ni, 1:nj, 0) = var_cell( 1:ni, 1:nj, 0)
         volume_framed(  1:ni, 1:nj, 0) = self%cell(1:ni, 1:nj, 0)%volume
      endif
      if (gc(6)>0) then
         var_cell_framed(ni, 1:nj, 1:nk+1) = var_cell( ni, 1:nj, 1:nk+1)
         volume_framed(  ni, 1:nj, 1:nk+1) = self%cell(ni, 1:nj, 1:nk+1)%volume
      endif
      ! interpolate on nodes
      do k=0, nk
         do j=0, nj
            do i=0, ni
               var_node(i, j, k) = (var_cell_framed(i+1, j+1, k+1) * volume_framed(i+1, j+1, k+1) &
                                 +  var_cell_framed(i  , j+1, k+1) * volume_framed(i  , j+1, k+1) &
                                 +  var_cell_framed(i+1, j  , k+1) * volume_framed(i+1, j  , k+1) &
                                 +  var_cell_framed(i  , j  , k+1) * volume_framed(i  , j  , k+1) &
                                 +  var_cell_framed(i+1, j+1, k  ) * volume_framed(i+1, j+1, k  ) &
                                 +  var_cell_framed(i  , j+1, k  ) * volume_framed(i  , j+1, k  ) &
                                 +  var_cell_framed(i+1, j  , k  ) * volume_framed(i+1, j  , k  ) &
                                 +  var_cell_framed(i  , j  , k  ) * volume_framed(i  , j  , k  ))&
                                 / sum(volume_framed(i:i+1, j:j+1, k:k+1))
            enddo
         enddo
      enddo
   endassociate
   endsubroutine interpolate_at_nodes

   subroutine load_conservatives_from_file(self, file_unit, pos)
   !< Load conservative variables from file.
   class(block_object), intent(inout) :: self      !< Block.
   integer(I4P),        intent(in)    :: file_unit !< File unit.
   integer(I4P),        intent(in)    :: pos       !< Position to start the loading.

   read(file_unit, pos=pos, iostat=self%error%status) self%cell%U%density, &
                                                      self%cell%U%momentum%x, self%cell%U%momentum%y, self%cell%U%momentum%z, &
                                                      self%cell%U%energy
   endsubroutine load_conservatives_from_file

   subroutine load_nodes_from_file(self, file_unit, pos)
   !< Load nodes from file.
   class(block_object), intent(inout) :: self      !< Block.
   integer(I4P),        intent(in)    :: file_unit !< File unit.
   integer(I4P),        intent(in)    :: pos       !< Position to start the loading.

   read(file_unit, pos=pos, iostat=self%error%status) self%node%vertex%x, self%node%vertex%y, self%node%vertex%z
   call self%compute_extents
   endsubroutine load_nodes_from_file

   elemental function nodes_number(self, with_ghosts) result(nodes_number_)
   !< Return the number of nodes.
   class(block_object), intent(in)           :: self          !< Block.
   logical,             intent(in), optional :: with_ghosts   !< Take into account ghost cells.
   integer(I4P)                              :: nodes_number_ !< Number of cells.

   nodes_number_ = self%signature%nodes_number(with_ghosts=with_ghosts)
   endfunction nodes_number

   subroutine save_conservatives_into_file(self, file_unit, pos)
   !< Save conservative variables into file.
   class(block_object), intent(inout) :: self      !< Block.
   integer(I4P),        intent(in)    :: file_unit !< File unit.
   integer(I4P),        intent(in)    :: pos       !< Position to start the loading.

   write(file_unit, pos=pos, iostat=self%error%status) self%cell%U%density, &
                                                       self%cell%U%momentum%x, self%cell%U%momentum%y, self%cell%U%momentum%z, &
                                                       self%cell%U%energy
   endsubroutine save_conservatives_into_file

   subroutine save_file_grid(self, file_name, ascii, metrics, tecplot, vtk)
   !< Save grid file.
   class(block_object), intent(inout)        :: self         !< Block.
   character(*),        intent(in)           :: file_name    !< File name.
   logical,             intent(in), optional :: ascii        !< Ascii/binary output.
   logical,             intent(in), optional :: metrics      !< Save metrics sentinel.
   logical,             intent(in), optional :: tecplot      !< Tecplot output format sentinel.
   logical,             intent(in), optional :: vtk          !< VTK output format sentinel.
   logical                                   :: tecplot_     !< Tecplot format sentinel, local variable.
   logical                                   :: vtk_         !< VTK format sentinel, local variable.

   tecplot_ = .false. ; if (present(tecplot)) tecplot_ = tecplot
   vtk_     = .false. ; if (present(vtk    )) vtk_     = vtk

   if (vtk_) call self%save_file_grid_vtk(file_name=file_name, ascii=ascii, metrics=metrics)
   endsubroutine save_file_grid

   subroutine save_file_solution(self, file_name, metrics, gc, ascii, tecplot, vtk)
   !< Save solution file.
   class(block_object), intent(inout)        :: self         !< Block.
   character(*),        intent(in)           :: file_name    !< File name.
   logical,             intent(in), optional :: metrics      !< Save metrics sentinel.
   logical,             intent(in), optional :: gc           !< Save ghost cells sentinel.
   logical,             intent(in), optional :: ascii        !< Ascii/binary output.
   logical,             intent(in), optional :: tecplot      !< Tecplot output format sentinel.
   logical,             intent(in), optional :: vtk          !< VTK output format sentinel.
   logical                                   :: tecplot_     !< Tecplot format sentinel, local variable.
   logical                                   :: vtk_         !< VTK format sentinel, local variable.

   tecplot_ = .false. ; if (present(tecplot)) tecplot_ = tecplot
   vtk_     = .false. ; if (present(vtk    )) vtk_     = vtk

   if (vtk_) call self%save_file_solution_vtk(file_name=file_name, metrics=metrics, gc=gc, ascii=ascii)
   endsubroutine save_file_solution

   subroutine save_bc_into_file(self, file_unit, pos)
   !< Save bc into file.
   class(block_object), intent(inout) :: self      !< Block.
   integer(I4P),        intent(in)    :: file_unit !< File unit.
   integer(I4P),        intent(in)    :: pos       !< Position to start the loading.
   integer(I4P)                       :: i, j, k   !< Counter.

   associate(gc=>self%signature%gc, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk)
      do k=1-gc(5), nk+gc(6)
         do j=1-gc(3), nj+gc(4)
            do i=1-gc(1), ni+gc(2)
            error stop 'to be implemented'
            enddo
         enddo
      enddo
   endassociate
   endsubroutine save_bc_into_file

   subroutine save_nodes_into_file(self, file_unit, pos)
   !< Save nodes into file.
   class(block_object), intent(inout) :: self      !< Block.
   integer(I4P),        intent(in)    :: file_unit !< File unit.
   integer(I4P),        intent(in)    :: pos       !< Position to start the loading.

   write(file_unit, pos=pos, iostat=self%error%status) self%node%vertex%x, self%node%vertex%y, self%node%vertex%z
   endsubroutine save_nodes_into_file

   elemental subroutine set_eos(self, eos)
   !< Set EOS.
   class(block_object),    intent(inout) :: self !< Block.
   type(eos_compressible), intent(in)    :: eos  !< EOS data.

   self%eos = eos
   endsubroutine set_eos

   _PURE_ subroutine update_recv_cells_number(self, me, block_to_process_map, recv_cells_number)
   !< Update receive cells number for multi-processes communication.
   class(block_object), intent(in)    :: self                     !< Block.
   integer(I4P),        intent(in)    :: me                       !< ID of current process.
   integer(I4P),        intent(in)    :: block_to_process_map(1:) !< Maps of processes ID for all blocks, [1:blocks_global_number].
   integer(I4P),        intent(inout) :: recv_cells_number(0:)    !< Number of receive cells from each process [0:procs_number-1].
   integer(I4P)                       :: i, j, k                  !< Counter.

   associate(gc=>self%signature%gc, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk, cell=>self%cell)
      do k=1, nk
         do j=1, nj
            if (cell(0   ,j,k)%bc%is(BC_ADJACENT)) call update_frame(gc=gc(1), frame=cell(1-gc(1):0    ,j,k), rcn=recv_cells_number)
            if (cell(ni+1,j,k)%bc%is(BC_ADJACENT)) call update_frame(gc=gc(2), frame=cell(ni+1:ni+gc(2),j,k), rcn=recv_cells_number)
         enddo
      enddo
      do k=1, nk
         do i=1, ni
            if (cell(i,0   ,k)%bc%is(BC_ADJACENT)) call update_frame(gc=gc(3), frame=cell(i,1-gc(3):0    ,k), rcn=recv_cells_number)
            if (cell(i,nj+1,k)%bc%is(BC_ADJACENT)) call update_frame(gc=gc(4), frame=cell(i,nj+1:nj+gc(4),k), rcn=recv_cells_number)
         enddo
      enddo
      do j=1, nj
         do i=1, ni
            if (cell(i,j,0   )%bc%is(BC_ADJACENT)) call update_frame(gc=gc(5), frame=cell(i,j,1-gc(5):0    ), rcn=recv_cells_number)
            if (cell(i,j,nk+1)%bc%is(BC_ADJACENT)) call update_frame(gc=gc(6), frame=cell(i,j,nk+1:nk+gc(6)), rcn=recv_cells_number)
         enddo
      enddo
   endassociate
   contains
      _PURE_ subroutine update_frame(gc, frame, rcn)
      !< Update counter for current frame.
      integer(I4P),      intent(in)    :: gc           !< Number of ghost cells.
      type(cell_object), intent(in)    :: frame(1-gc:) !< Cells frame [1-gc:0].
      integer(I4P),      intent(inout) :: rcn(0:)      !< Number of receive cells from each process [0:procs_number-1].
      integer(I4P)                     :: ii           !< Counter.
      integer(I4P)                     :: b            !< Block ID in the global numeration.
      integer(I4P)                     :: proc         !< Processor ID.

      do ii=1-gc, 0
         b = frame(ii)%bc%adj(1)
         proc = block_to_process_map(b)
         if ( proc /= me) then ! check if proc/=me, me doesn't communicate with itself
            rcn(proc) = rcn(proc) + 1
         endif
      enddo
      endsubroutine update_frame
   endsubroutine update_recv_cells_number

   elemental subroutine update_level_set_distance(self)
   !< Update level set distance.
   !<
   !< @note Also recompute level set normal.
   class(block_object), intent(inout) :: self    !< Block.
   integer(I4P)                       :: i, j, k !< Counter.

   associate(ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk, cell=>self%cell)
      ! update level set distance
      do k=1, nk
         do j=1, nj
            do i=1, ni
               call cell(i,j,k)%level_set%update_distance
            enddo
         enddo
      enddo
      ! recompute level set normal
      do k=1, nk
         do j=1, nj
            do i=1, ni
               ! if (cell(i)%level_set%distance<0._R8P.and.cell(i+1)%level_set%distance>0._R8P) then
                  cell(i,j,k)%level_set%normal = level_set_normal(cell=[cell(i-1,j  ,k  )%level_set, cell(i+1,j  ,k  )%level_set,  &
                                                                        cell(i  ,j-1,k  )%level_set, cell(i  ,j+1,k  )%level_set,  &
                                                                        cell(i  ,j  ,k-1)%level_set, cell(i  ,j  ,k+1)%level_set], &
                                                                  dx=cell(i+1,j,k)%center%x-cell(i-1,j,k)%center%x,                &
                                                                  dy=cell(i,j+1,k)%center%y-cell(i,j-1,k)%center%y,                &
                                                                  dz=cell(i,j,k+1)%center%z-cell(i,j,k-1)%center%z)
               ! endif
            enddo
         enddo
      enddo
   endassociate
   endsubroutine update_level_set_distance

   ! fast operators
   ! +
   pure subroutine conservative_add_conservatives_fast(opr, lhs, rhs)
   !< `+` fast operator.
   class(block_object), intent(inout) :: opr   !< Operator result.
   type(block_object),  intent(in)    :: lhs   !< Left hand side.
   type(block_object),  intent(in)    :: rhs   !< Right hand side.
   integer(I4P)                       :: i,j,k !< Counter.

   do k=1, opr%signature%nk
      do j=1, opr%signature%nj
         do i=1, opr%signature%ni
            if (opr%cell(i,j,k)%level_set%distance<0._R8P) then
               cycle
            else
               call opr%cell(i,j,k)%U%field_add_field_fast(lhs=lhs%cell(i,j,k)%U, rhs=rhs%cell(i,j,k)%U)
            endif
         enddo
      enddo
   enddo
   endsubroutine conservative_add_conservatives_fast

   ! *
   pure subroutine conservative_multiply_conservatives_fast(opr, lhs, rhs)
   !< `*` fast operator.
   class(block_object), intent(inout) :: opr   !< Operator result.
   type(block_object),  intent(in)    :: lhs   !< Left hand side.
   type(block_object),  intent(in)    :: rhs   !< Right hand side.
   integer(I4P)                       :: i,j,k !< Counter.

   do k=1, opr%signature%nk
      do j=1, opr%signature%nj
         do i=1, opr%signature%ni
            if (opr%cell(i,j,k)%level_set%distance<0._R8P) then
               cycle
            else
               opr%cell(i,j,k)%U = lhs%cell(i,j,k)%U * rhs%cell(i,j,k)%U
            endif
         enddo
      enddo
   enddo
   endsubroutine conservative_multiply_conservatives_fast

   pure subroutine conservative_multiply_real_scalar_fast(opr, lhs, rhs)
   !< `* real_scalar` fast operator.
   class(block_object), intent(inout) :: opr   !< Operator result.
   type(block_object),  intent(in)    :: lhs   !< Left hand side.
   real(R8P),           intent(in)    :: rhs   !< Right hand side.
   integer(I4P)                       :: i,j,k !< Counter.

   do k=1, opr%signature%nk
      do j=1, opr%signature%nj
         do i=1, opr%signature%ni
            if (opr%cell(i,j,k)%level_set%distance<0._R8P) then
               cycle
            else
               opr%cell(i,j,k)%U = lhs%cell(i,j,k)%U * rhs
            endif
         enddo
      enddo
   enddo
   endsubroutine conservative_multiply_real_scalar_fast

   ! -
   pure subroutine conservative_subtract_conservatives_fast(opr, lhs, rhs)
   !< `-` fast operator.
   class(block_object), intent(inout) :: opr   !< Operator result.
   type(block_object),  intent(in)    :: lhs   !< Left hand side.
   type(block_object),  intent(in)    :: rhs   !< Right hand side.
   integer(I4P)                       :: i,j,k !< Counter.

   do k=1, opr%signature%nk
      do j=1, opr%signature%nj
         do i=1, opr%signature%ni
            if (opr%cell(i,j,k)%level_set%distance<0._R8P) then
               cycle
            else
               call opr%cell(i,j,k)%U%field_subtract_field_fast(lhs=lhs%cell(i,j,k)%U, rhs=rhs%cell(i,j,k)%U)
            endif
         enddo
      enddo
   enddo
   endsubroutine conservative_subtract_conservatives_fast

   ! private methods
   pure subroutine block_assign_block(lhs, rhs)
   !< Operator `=`.
   class(block_object), intent(inout) :: lhs !< Left hand side.
   type(block_object),  intent(in)    :: rhs !< Right hand side.

   lhs%error = rhs%error
   lhs%signature = rhs%signature
   if (lhs%signature%cells_number()==rhs%signature%cells_number()) then
      lhs%cell = rhs%cell
      lhs%face_i = rhs%face_i
      lhs%face_j = rhs%face_j
      lhs%face_k = rhs%face_k
      lhs%node = rhs%node
   endif
   lhs%eos = rhs%eos
   endsubroutine block_assign_block

   elemental subroutine compute_cells_center(self)
   !< Compute cells center.
   class(block_object), intent(inout) :: self !< Block.
   integer(I4P)                       :: i    !< Counter.
   integer(I4P)                       :: j    !< Counter.
   integer(I4P)                       :: k    !< Counter.

   associate(gc=>self%signature%gc, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk)
      do k=1 - gc(5), nk + gc(6)
         do j=1 - gc(3), nj + gc(4)
            do i=1 - gc(1), ni + gc(2)
               self%cell(i, j, k)%center = (self%node(i,   j,   k  )%vertex + &
                                            self%node(i-1, j,   k  )%vertex + &
                                            self%node(i  , j-1, k  )%vertex + &
                                            self%node(i  , j  , k-1)%vertex + &
                                            self%node(i-1, j-1, k-1)%vertex + &
                                            self%node(i  , j-1, k-1)%vertex + &
                                            self%node(i-1, j  , k-1)%vertex + &
                                            self%node(i-1, j-1, k  )%vertex) * 0.125_R8P
            enddo
         enddo
      enddo
   endassociate
   endsubroutine compute_cells_center

   elemental subroutine compute_extents(self)
   !< Compute block extents.
   class(block_object), intent(inout) :: self !< Block.

   self%error%status = ERROR_BLOCK_COMPUTE_EXTENTS_FAILED

   associate(ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk)
      if (allocated(self%node)) then
         self%signature%emin%x = minval(self%node(1:ni, 1:nj, 1:nk)%vertex%x)
         self%signature%emin%y = minval(self%node(1:ni, 1:nj, 1:nk)%vertex%y)
         self%signature%emin%z = minval(self%node(1:ni, 1:nj, 1:nk)%vertex%z)

         self%signature%emax%x = maxval(self%node(1:ni, 1:nj, 1:nk)%vertex%x)
         self%signature%emax%y = maxval(self%node(1:ni, 1:nj, 1:nk)%vertex%y)
         self%signature%emax%z = maxval(self%node(1:ni, 1:nj, 1:nk)%vertex%z)

         self%error%status = NO_ERROR
      endif
   endassociate
   endsubroutine compute_extents

   elemental subroutine compute_faces_metrics(self)
   !< Compute block faces metrics.
   class(block_object), intent(inout) :: self         !< Block.
   type(vector)                       :: triplet(1:4) !< Dummy vectors.
   real(R8P)                          :: signi        !< Sign of direction of normals along I coordinate.
   real(R8P)                          :: signj        !< Sign of direction of normals along J coordinate.
   real(R8P)                          :: signk        !< Sign of direction of normals along K coordinate.
   integer(I4P)                       :: i            !< Counter.
   integer(I4P)                       :: j            !< Counter.
   integer(I4P)                       :: k            !< Counter.

   associate(gc=>self%signature%gc, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk, &
             node=>self%node)
      i = max(1, ni)
      j = max(1, nj)
      k = max(1, nk)
      triplet(1) = node(i, j,   k)%vertex - node(i, j-1, k-1)%vertex
      triplet(2) = node(i, j-1, k)%vertex - node(i, j,   k-1)%vertex
      triplet(3) = triplet(1).cross.triplet(2)
      triplet(4) = (0.25_R8P * (node(i,   j, k  )%vertex + node(i,   j-1, k  )%vertex   + &
                                node(i,   j, k-1)%vertex + node(i,   j-1, k-1)%vertex)) - &
                   (0.25_R8P * (node(i-1, j, k  )%vertex + node(i-1, j-1, k  )%vertex   + &
                                node(i-1, j, k-1)%vertex + node(i-1, j-1, k-1)%vertex))
      signi = sign(1._R8P, (triplet(3).dot.triplet(4)))

      triplet(1) = node(i, j, k  )%vertex - node(i-1, j, k-1)%vertex
      triplet(2) = node(i, j, k-1)%vertex - node(i-1, j, k  )%vertex
      triplet(3) = triplet(1).cross.triplet(2)
      triplet(4) = (0.25_R8P * (node(i, j,   k  )%vertex + node(i-1, j,   k  )%vertex   + &
                                node(i, j,   k-1)%vertex + node(i-1, j,   k-1)%vertex)) - &
                   (0.25_R8P * (node(i, j-1, k  )%vertex + node(i-1, j-1, k  )%vertex   + &
                                node(i, j-1, k-1)%vertex + node(i-1, j-1, k-1)%vertex))
      signj = sign(1._R8P, (triplet(3).dot.triplet(4)))

      triplet(1) = node(i,   j, k)%vertex - node(i-1, j-1, k)%vertex
      triplet(2) = node(i-1, j, k)%vertex - node(i,   j-1, k)%vertex
      triplet(3) = triplet(1).cross.triplet(2)
      triplet(4) = (0.25_R8P * (node(i, j, k    )%vertex + node(i-1, j, k    )%vertex   + &
                                node(i, j-1, k  )%vertex + node(i-1, j-1, k  )%vertex)) - &
                   (0.25_R8P * (node(i, j, k-1  )%vertex + node(i-1, j, k-1  )%vertex   + &
                                node(i, j-1, k-1)%vertex + node(i-1, j-1, k-1)%vertex))
      signk = sign(1._R8P, (triplet(3).dot.triplet(4)))

      do k=1-gc(5), nk+gc(6)
         do j=1-gc(3), nj+gc(4)
            do i=0-gc(1)+1, ni+gc(2)-1
               call self%face_i(i, j, k)%compute_metrics(pt1 = node(i, j-1, k-1)%vertex, &
                                                         pt2 = node(i, j  , k-1)%vertex, &
                                                         pt3 = node(i, j  , k  )%vertex, &
                                                         pt4 = node(i, j-1, k  )%vertex, signd=signi)
            enddo
         enddo
      enddo

      do k=1-gc(5), nk+gc(6)
         do j=0-gc(3)+1, nj+gc(4)-1
            do i=1-gc(1), ni+gc(2)
               call self%face_j(i, j, k)%compute_metrics(pt1 = node(i-1, j, k-1)%vertex, &
                                                         pt2 = node(i-1, j, k  )%vertex, &
                                                         pt3 = node(i  , j, k  )%vertex, &
                                                         pt4 = node(i  , j, k-1)%vertex, signd=signj)
            enddo
         enddo
      enddo

      do k=0-gc(5)+1, nk+gc(6)-1
         do j=1-gc(3), nj+gc(4)
            do i=1-gc(1), ni+gc(2)
               call self%face_k(i, j, k)%compute_metrics(pt1 = node(i-1, j-1, k)%vertex, &
                                                         pt2 = node(i  , j-1, k)%vertex, &
                                                         pt3 = node(i  , j  , k)%vertex, &
                                                         pt4 = node(i-1, j  , k)%vertex, signd=signk)
            enddo
         enddo
      enddo
   endassociate
   endsubroutine compute_faces_metrics

   subroutine compute_metrics(self)
   !< Compute block metrics.
   !<
   !< @TODO re-add elemental attribute.
   !<
   !< @TODO re-add metrics correction call.
   class(block_object), intent(inout) :: self !< Block.

   call self%compute_cells_center
   call self%compute_faces_metrics
   call self%compute_volumes
   ! call self%correct_metrics
   call self%nullify_normals
   endsubroutine compute_metrics

   elemental subroutine compute_volumes(self)
   !< Compute block volumes.
   !<
   !< The volume of each cell is computed using the formula:
   !< $$
   !< v = [(\vec n_7 - \vec n_1) + (\vec n_6 - \vec n_0), (\vec n_7 - \vec n_2), (\vec n_3 - \vec n_0)] +
   !<     [(\vec n_6 - \vec n_0), (\vec n_7 - \vec n_2) + (\vec n_5 - \vec n_0), (\vec n_7 - \vec n_4)] +
   !<     [(\vec n_7 - \vec n_1), (\vec n_5 - \vec n_0), (\vec n_7 - \vec n_4) + (\vec n_3 - \vec n_0)]
   !< $$
   !< where \([\vec A, \vec B, \vec C]=\vec A \cdot (\vec B \times \vec C)\) is the triple product.
   !<
   !<### References
   !<
   !< [1] *Efficient computation of volume of hexahedral cells*, Grandy J., 1997.
   class(block_object), intent(inout) :: self         !< Block.
   type(vector)                       :: triplet(1:9) !< Dummy vectors.
   integer(I4P)                       :: i            !< Counter.
   integer(I4P)                       :: j            !< Counter.
   integer(I4P)                       :: k            !< Counter.

   associate(gc=>self%signature%gc, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk, &
             node=>self%node)
      do k=1-gc(5)+1, nk+gc(6)-1
         do j=1-gc(3)+1, nj+gc(4)-1
            do i=1-gc(1)+1, ni+gc(2)-1
               triplet(1) = node(i  , j  , k  )%vertex - node(i  , j-1, k-1)%vertex &
                          + node(i-1, j  , k  )%vertex - node(i-1, j-1, k-1)%vertex
               triplet(2) = node(i  , j  , k  )%vertex - node(i-1, j  , k-1)%vertex
               triplet(3) = node(i  , j  , k-1)%vertex - node(i-1, j-1, k-1)%vertex
               triplet(4) = node(i-1, j  , k  )%vertex - node(i-1, j-1, k-1)%vertex
               triplet(5) = node(i  , j  , k  )%vertex - node(i-1, j  , k-1)%vertex &
                          + node(i  , j-1, k  )%vertex - node(i-1, j-1, k-1)%vertex
               triplet(6) = node(i  , j  , k  )%vertex - node(i-1, j-1, k  )%vertex
               triplet(7) = node(i  , j  , k  )%vertex - node(i  , j-1, k-1)%vertex
               triplet(8) = node(i  , j-1, k  )%vertex - node(i-1, j-1, k-1)%vertex
               triplet(9) = node(i  , j  , k  )%vertex - node(i-1, j-1, k  )%vertex &
                          + node(i  , j  , k-1)%vertex - node(i-1, j-1, k-1)%vertex
               self%cell(i, j, k)%volume = ((triplet(1).dot.(triplet(2).cross.triplet(3))) + &
                                            (triplet(4).dot.(triplet(5).cross.triplet(6))) + &
                                            (triplet(7).dot.(triplet(8).cross.triplet(9)))) / 12.0_R8P
            enddo
         enddo
      enddo
   endassociate
   endsubroutine compute_volumes

   subroutine correct_metrics(self)
   !> Correct the metrics.
   !<
   !< Check for boundary conditions and volumes issues.
   !<
   !< @TODO Implement correction.
   class(block_object), intent(inout) :: self !< Block.

   error stop 'error: block metrics correction to be implemented'
   endsubroutine correct_metrics

   pure function node_to_center(self) result(center)
   !< Compute cell centers coordinates from cell nodes.
   class(block_object), intent(in) :: self          !< Block.
   type(vector), allocatable       :: center(:,:,:) !< Cell centers coordinates.
   integer(I4P)                    :: i             !< Counter.
   integer(I4P)                    :: j             !< Counter.
   integer(I4P)                    :: k             !< Counter.

   associate(gc=>self%signature%gc, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk)
      allocate(center(1 - gc(1) : ni + gc(2), 1 - gc(3) : nj + gc(4), 1 - gc(5) : nk + gc(6)))
      do k=1 - gc(5), nk + gc(6)
         do j=1 - gc(3), nj + gc(4)
            do i=1 - gc(1), ni + gc(2)
               center(i, j, k) = (self%node(i,   j,   k  )%vertex + &
                                  self%node(i-1, j,   k  )%vertex + &
                                  self%node(i  , j-1, k  )%vertex + &
                                  self%node(i  , j  , k-1)%vertex + &
                                  self%node(i-1, j-1, k-1)%vertex + &
                                  self%node(i  , j-1, k-1)%vertex + &
                                  self%node(i-1, j  , k-1)%vertex + &
                                  self%node(i-1, j-1, k  )%vertex) * 0.125_R8P
            enddo
         enddo
      enddo
   endassociate
   endfunction node_to_center

   elemental subroutine nullify_normals(self)
   !< Nullify normals for 2D or 1D domains.
   class(block_object), intent(inout) :: self !< Block.

   if (self%signature%is_null_x) then
      self%face_i%normal = (self%face_i%normal.paral.ex) + (0._R8P * ey                ) + (0._R8P * ez                )
      self%face_j%normal = (0._R8P * ex                ) + (self%face_j%normal.paral.ey) + (self%face_j%normal.paral.ez)
      self%face_k%normal = (0._R8P * ex                ) + (self%face_k%normal.paral.ey) + (self%face_k%normal.paral.ez)
   endif
   if (self%signature%is_null_y) then
      self%face_i%normal = (self%face_i%normal.paral.ex) + (0._R8P * ey                ) + (self%face_i%normal.paral.ez)
      self%face_j%normal = (0._R8P * ex                ) + (self%face_j%normal.paral.ey) + (0._R8P * ez                )
      self%face_k%normal = (self%face_k%normal.paral.ex) + (0._R8P * ey                ) + (self%face_k%normal.paral.ez)
   endif
   if (self%signature%is_null_z) then
      self%face_i%normal = (self%face_i%normal.paral.ex) + (self%face_i%normal.paral.ey) + (0._R8P * ez                )
      self%face_j%normal = (self%face_j%normal.paral.ex) + (self%face_j%normal.paral.ey) + (0._R8P * ez                )
      self%face_k%normal = (0._R8P * ex                ) + (0._R8P * ey                ) + (self%face_k%normal.paral.ez)
   endif
   endsubroutine nullify_normals

   subroutine save_file_grid_tec(self, file_name, ascii, metrics)
   !< Save mesh grid into Tecplot file.
   !<
   !< @TODO implement Tecplot output.
   class(block_object), intent(inout)        :: self      !< Block.
   character(*),        intent(in)           :: file_name !< Output file name.
   logical,             intent(in), optional :: ascii     !< Ascii/binary output.
   logical,             intent(in), optional :: metrics   !< Save also metrics data.

   error stop 'error: block mesh Tecplot output to be implemented'
   endsubroutine save_file_grid_tec

   subroutine save_file_grid_vtk(self, file_name, ascii, metrics)
   !< Save mesh grid into VTK file.
   class(block_object), intent(inout)        :: self      !< Block.
   character(*),        intent(in)           :: file_name !< Output file name.
   logical,             intent(in), optional :: ascii     !< Ascii/binary output.
   logical,             intent(in), optional :: metrics   !< Save also metrics data.
   logical                                   :: ascii_    !< Ascii/binary output.
   logical                                   :: metrics_  !< Save also metrics data.
   type(vtk_file)                            :: vtk       !< VTK file.

   ascii_   = .false. ; if (present(ascii  )) ascii_   = ascii
   metrics_ = .false. ; if (present(metrics)) metrics_ = metrics

   associate(node=>self%node, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk, &
             nn=>self%nodes_number(with_ghosts=.false.))
      if (ascii_) then
         self%error%status = vtk%initialize(format='ascii',                                                    &
                                            filename=trim(adjustl(file_name)), mesh_topology='StructuredGrid', &
                                            nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
      else
         self%error%status = vtk%initialize(format='raw',                                                      &
                                            filename=trim(adjustl(file_name)), mesh_topology='StructuredGrid', &
                                            nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
      endif

      self%error%status = vtk%xml_writer%write_piece(nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
      self%error%status = vtk%xml_writer%write_geo(n=nn, x=node(0:ni, 0:nj, 0:nk)%vertex%x, &
                                                         y=node(0:ni, 0:nj, 0:nk)%vertex%y, &
                                                         z=node(0:ni, 0:nj, 0:nk)%vertex%z)
      if (metrics_) then
         self%error%status = vtk%xml_writer%write_dataarray(location='cell', action='open')
         self%error%status = vtk%xml_writer%write_dataarray(data_name='volume', x=self%cell(1:ni, 1:nj, 1:nk)%volume, &
                                                            one_component=.true.)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='area_i', x=self%face_i(1:ni, 1:nj, 1:nk)%area, &
                                                            one_component=.true.)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='area_j', x=self%face_j(1:ni, 1:nj, 1:nk)%area, &
                                                            one_component=.true.)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='area_k', x=self%face_k(1:ni, 1:nj, 1:nk)%area, &
                                                            one_component=.true.)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='normals_i', x=self%face_i(1:ni, 1:nj, 1:nk)%normal%x, &
                                                                                   y=self%face_i(1:ni, 1:nj, 1:nk)%normal%y, &
                                                                                   z=self%face_i(1:ni, 1:nj, 1:nk)%normal%z)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='normals_j', x=self%face_j(1:ni, 1:nj, 1:nk)%normal%x, &
                                                                                   y=self%face_j(1:ni, 1:nj, 1:nk)%normal%y, &
                                                                                   z=self%face_j(1:ni, 1:nj, 1:nk)%normal%z)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='normals_k', x=self%face_k(1:ni, 1:nj, 1:nk)%normal%x, &
                                                                                   y=self%face_k(1:ni, 1:nj, 1:nk)%normal%y, &
                                                                                   z=self%face_k(1:ni, 1:nj, 1:nk)%normal%z)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='distance',                             &
                                                            x=self%cell(1:ni, 1:nj, 1:nk)%level_set%distance, &
                                                            one_component=.true.)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='ls_normal', &
                                                            x=self%cell(1:ni, 1:nj, 1:nk)%level_set%normal%x, &
                                                            y=self%cell(1:ni, 1:nj, 1:nk)%level_set%normal%y, &
                                                            z=self%cell(1:ni, 1:nj, 1:nk)%level_set%normal%z)
         self%error%status = vtk%xml_writer%write_dataarray(location='cell', action='close')
      endif
      self%error%status = vtk%xml_writer%write_piece()
      self%error%status = vtk%finalize()
   endassociate
   endsubroutine save_file_grid_vtk

   subroutine save_file_solution_vtk(self, file_name, metrics, gc, ascii)
   !< Save mesh solution into VTK file.
   class(block_object), intent(inout)        :: self                      !< Block.
   character(*),        intent(in)           :: file_name                 !< Output file name.
   logical,             intent(in), optional :: metrics                   !< Save metrics sentinel.
   logical,             intent(in), optional :: gc                        !< Save ghost cells sentinel.
   logical,             intent(in), optional :: ascii                     !< Ascii/binary output.
   logical                                   :: metrics_                  !< Save metrics sentinel, local variable.
   logical                                   :: gc_                       !< Save ghost cells sentinel, local variable.
   logical                                   :: ascii_                    !< Ascii/binary output.
   type(vtk_file)                            :: vtk                       !< VTK file.
   real(R8P), allocatable                    :: distances(:,:,:)          !< Interfaces distances.
   integer(I4P)                              :: i1, i2, j1, j2, k1, k2, d !< Counter.
   integer(I4P)                              :: i, j, k                   !< Counter.

   metrics_ = .false. ; if (present(metrics)) metrics_ = metrics
   gc_      = .false. ; if (present(gc     )) gc_      = gc
   ascii_   = .false. ; if (present(ascii  )) ascii_   = ascii

   if (gc_) then
      ! save also ghost cells frames
      i1 = 1 - self%signature%gc(1) ; i2 = self%signature%ni + self%signature%gc(2)
      j1 = 1 - self%signature%gc(3) ; j2 = self%signature%nj + self%signature%gc(4)
      k1 = 1 - self%signature%gc(5) ; k2 = self%signature%nk + self%signature%gc(6)
   else
      ! save only internal cells
      i1 = 1                        ; i2 = self%signature%ni
      j1 = 1                        ; j2 = self%signature%nj
      k1 = 1                        ; k2 = self%signature%nk
   endif

   associate(node=>self%node, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk, &
             nn=>self%nodes_number(with_ghosts=gc_))
      if (ascii_) then
         self%error%status = vtk%initialize(format='ascii',                                                    &
                                            filename=trim(adjustl(file_name)), mesh_topology='StructuredGrid', &
                                            nx1=i1-1, nx2=i2, ny1=j1-1, ny2=j2, nz1=k1-1, nz2=k2)
      else
         self%error%status = vtk%initialize(format='raw',                                                      &
                                            filename=trim(adjustl(file_name)), mesh_topology='StructuredGrid', &
                                            nx1=i1-1, nx2=i2, ny1=j1-1, ny2=j2, nz1=k1-1, nz2=k2)
      endif

      self%error%status = vtk%xml_writer%write_piece(nx1=i1-1, nx2=i2, ny1=j1-1, ny2=j2, nz1=k1-1, nz2=k2)

      self%error%status = vtk%xml_writer%write_geo(n=nn, x=node(i1-1:i2, j1-1:j2, k1-1:k2)%vertex%x, &
                                                         y=node(i1-1:i2, j1-1:j2, k1-1:k2)%vertex%y, &
                                                         z=node(i1-1:i2, j1-1:j2, k1-1:k2)%vertex%z)

      self%error%status = vtk%xml_writer%write_dataarray(location='cell', action='open')

      self%error%status = vtk%xml_writer%write_dataarray(data_name='density', x=self%cell(i1:i2,j1:j2,k1:k2)%U%density, &
                                                         one_component=.true.)
      self%error%status = vtk%xml_writer%write_dataarray(data_name='energy', x=self%cell(i1:i2,j1:j2,k1:k2)%U%energy, &
                                                         one_component=.true.)
      self%error%status = vtk%xml_writer%write_dataarray(data_name='momentum', x=self%cell(i1:i2,j1:j2,k1:k2)%U%momentum%x, &
                                                                               y=self%cell(i1:i2,j1:j2,k1:k2)%U%momentum%y, &
                                                                               z=self%cell(i1:i2,j1:j2,k1:k2)%U%momentum%z)

      if (metrics_) then
         self%error%status = vtk%xml_writer%write_dataarray(data_name='volume', x=self%cell(i1:i2,j1:j2,k1:k2)%volume, &
                                                            one_component=.true.)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='area_i', x=self%face_i(i1:i2,j1:j2,k1:k2)%area, &
                                                            one_component=.true.)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='area_j', x=self%face_j(i1:i2,j1:j2,k1:k2)%area, &
                                                            one_component=.true.)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='area_k', x=self%face_k(i1:i2,j1:j2,k1:k2)%area, &
                                                            one_component=.true.)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='normals_i', x=self%face_i(i1:i2,j1:j2,k1:k2)%normal%x, &
                                                                                   y=self%face_i(i1:i2,j1:j2,k1:k2)%normal%y, &
                                                                                   z=self%face_i(i1:i2,j1:j2,k1:k2)%normal%z)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='normals_j', x=self%face_j(i1:i2,j1:j2,k1:k2)%normal%x, &
                                                                                   y=self%face_j(i1:i2,j1:j2,k1:k2)%normal%y, &
                                                                                   z=self%face_j(i1:i2,j1:j2,k1:k2)%normal%z)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='normals_k', x=self%face_k(i1:i2,j1:j2,k1:k2)%normal%x, &
                                                                                   y=self%face_k(i1:i2,j1:j2,k1:k2)%normal%y, &
                                                                                   z=self%face_k(i1:i2,j1:j2,k1:k2)%normal%z)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='distance',                              &
                                                            x=self%cell(i1:i2,j1:j2,k1:k2)%level_set%distance, &
                                                            one_component=.true.)
         self%error%status = vtk%xml_writer%write_dataarray(data_name='ls_normal', &
                                                            x=self%cell(1:ni, 1:nj, 1:nk)%level_set%normal%x, &
                                                            y=self%cell(1:ni, 1:nj, 1:nk)%level_set%normal%y, &
                                                            z=self%cell(1:ni, 1:nj, 1:nk)%level_set%normal%z)
         do d=1, self%signature%interfaces_number
            allocate(distances(i1:i2,j1:j2,k1:k2))
            do K=k1, k2
               do j=j1, j2
                  do i=i1, i2
                     distances(i, j, k) = self%cell(i, j, k)%level_set%distances(d)
                  enddo
               enddo
            enddo
            self%error%status = vtk%xml_writer%write_dataarray(data_name='distance_'//trim(str(d,.true.)), &
                                                               x=distances, one_component=.true.)
            deallocate(distances)
         enddo
      endif

      self%error%status = vtk%xml_writer%write_dataarray(location='cell', action='close')
      self%error%status = vtk%xml_writer%write_piece()
      self%error%status = vtk%finalize()
   endassociate
   endsubroutine save_file_solution_vtk

   ! private non TBP
   subroutine compute_fluxes_convective(solver, eos, gc, N, faces, cells, fluxes)
   !< Compute the conservative fluxes.
   class(solver_object),            intent(in)    :: solver          !< Solver.
   type(eos_compressible),          intent(in)    :: eos             !< EOS data.
   integer(I4P),                    intent(in)    :: gc              !< Number of ghost cells used.
   integer(I4P),                    intent(in)    :: N               !< Number of cells.
   type(face_object),               intent(in)    :: faces(0-gc:)    !< Faces data [0-gc:N+gc].
   type(cell_object),               intent(in)    :: cells(1-gc:)    !< Cells data [1-gc:N+gc].
   type(conservative_compressible), intent(inout) :: fluxes(0:)      !< Convective fluxes of 3D conservative variables [0:N].
   type(conservative_compressible), allocatable   :: UR(:,:)         !< Reconstructed conservative variables.
   type(vector), allocatable                      :: tangential(:,:) !< Interface tangential component of velocity.
   integer(I4P)                                   :: i               !< Counter.

   allocate(UR(1:2,0:N+1))
   allocate(tangential(1:2,0:N+1))
   call reconstruct_interfaces_characteristic(solver=solver, eos=eos, gc=gc, N=N, U=cells%U, level_set=cells%level_set, &
                                              normal=faces%normal, UR=UR, tangential=tangential)
   do i=0, N
      ! computing normal fluxes
      call solver%riemann_solver%solve(eos_left=eos,  state_left=UR( 2, i  ), &
                                       eos_right=eos, state_right=UR(1, i+1), normal=faces(i)%normal, fluxes=fluxes(i))
      ! uptdating fluxes with tangential components
      if (fluxes(i)%density > 0._R8P) then
         fluxes(i)%momentum = fluxes(i)%momentum + fluxes(i)%density *         tangential(2,i)
         fluxes(i)%energy   = fluxes(i)%energy   + fluxes(i)%density * sq_norm(tangential(2,i)) * 0.5_R8P
      else
         fluxes(i)%momentum = fluxes(i)%momentum + fluxes(i)%density *         tangential(1,i+1)
         fluxes(i)%energy   = fluxes(i)%energy   + fluxes(i)%density * sq_norm(tangential(1,i+1)) * 0.5_R8P
      endif
   enddo
   endsubroutine compute_fluxes_convective

   subroutine reconstruct_interfaces_characteristic(solver, eos, gc, N, U, level_set, normal, UR, tangential)
   !< Reconstruct interfaces states.
   !<
   !< The reconstruction is done in pseudo characteristic variables.
   class(solver_object),            intent(in)    :: solver                  !< Solver.
   type(eos_compressible),          intent(in)    :: eos                     !< EOS data.
   integer(I4P),                    intent(in)    :: gc                      !< Number of ghost cells used.
   integer(I4P),                    intent(in)    :: N                       !< Number of cells.
   type(conservative_compressible), intent(in)    :: U(1-gc:)                !< Conservative variables.
   type(level_set_object),          intent(in)    :: level_set(1-gc:)        !< Level set cells data.
   type(vector),                    intent(in)    :: normal(0:)              !< Face normals.
   type(conservative_compressible), intent(inout) :: UR(1:, 0:)              !< Reconstructed conservative vars.
   type(vector),                    intent(inout) :: tangential(1:, 0:)      !< Interface tangential component of velocity.
   type(conservative_compressible)                :: U_(1-gc:N+gc)           !< Conservative variables, local.
   type(primitive_compressible)                   :: P(1-gc:N+gc)            !< Primitive variables.
   type(primitive_compressible)                   :: PR(1:2, 0:N+1)          !< Reconstructed primitive variables.
   type(primitive_compressible)                   :: Pm(1:2)                 !< Mean of primitive variables.
   real(R8P)                                      :: LPm(1:3, 1:3, 1:2)      !< Mean left eigenvectors matrix.
   real(R8P)                                      :: RPm(1:3, 1:3, 1:2)      !< Mean right eigenvectors matrix.
   real(R8P)                                      :: C(1:2, 1-gc:-1+gc, 1:3) !< Pseudo characteristic variables.
   real(R8P)                                      :: CR(1:2, 1:3)            !< Pseudo characteristic reconst.
   real(R8P)                                      :: buffer(1:3)             !< Dummy buffer.
   integer(I4P)                                   :: i, j, f, v              !< Counter.

   ! store temporary conservative variables
   U_ = U

   ! impose immersed boundaries
   !< @TODO generalize this
   do i=0, N
      if (level_set(i)%distance<0._R8P.and.level_set(i+1)%distance>0._R8P) then
         ! left interface is an immersed boundary: note that there are necessary gc cells immersed to the left...
         do j=i+1-gc, i
            if (level_set(j)%distance<0._R8P) then
               ! U_(j) = U_(i+(i-j)+1)
               U_(j)%momentum = U_(j)%momentum - (2._R8P*(U_(j)%momentum.paral.normal(i)))
               ! U_(j)%momentum = U_(j)%momentum - (2._R8P*(U_(j)%momentum.paral.level_set(i+1)%normal))
            endif
         enddo
      elseif (level_set(i)%distance>0._R8P.and.level_set(i+1)%distance<0._R8P) then
         ! right interface is an immersed boundary: note that there are necessary gc cells immersed to the right...
         do j=i+1, i+gc
            if (level_set(j)%distance<0._R8P) then
               ! U_(j) = U_(i+(i-j)+1)
               U_(j)%momentum = U_(j)%momentum - (2._R8P*(U_(j)%momentum.paral.normal(i)))
               ! U_(j)%momentum = U_(j)%momentum - (2._R8P*(U_(j)%momentum.paral.level_set(i)%normal))
            endif
         enddo
      endif
   enddo

   ! compute tangential component of velocity
   do i=0, N+1
      do f=1, 2
         if (i==0  .and.f==1) cycle
         if (i==N+1.and.f==2) cycle
         tangential(f,i) = U_(i)%velocity() - (U_(i)%velocity().paral.normal(i+f-1))
      enddo
   enddo

   ! compute reconstruction of conservative variables
   select case(gc)
   case(1) ! 1st order piecewise constant reconstruction
      do i=0, N+1
         UR(1, i) = U_(i)
         UR(2, i) = UR(1, i)
      enddo
   case default ! 3rd-17th order WENO reconstruction
      ! compute primitive variables
      do i=1-gc, N+gc
         P(i) = conservative_to_primitive_compressible(conservative=U_(i), eos=eos)
      enddo

      ! compute WENO reconstruction
      do i=0, N+1
         ! compute pseudo characteristic variables
         do f=1, 2
            if (i==0  .and.f==1) cycle
            if (i==N+1.and.f==2) cycle
            ! Pm(f) = 0.5_R8P * (P(i+f-2) + P(i+f-1))
            call Pm(f)%field_add_field_fast(lhs=P(i+f-2), rhs=P(i+f-1))
            call Pm(f)%field_multiply_real_scalar_fast(lhs=Pm(f), rhs=0.5_R8P)
         enddo
         do f=1, 2
            if (i==0  .and.f==1) cycle
            if (i==N+1.and.f==2) cycle
            LPm(:, :, f) = Pm(f)%left_eigenvectors(eos=eos)
            RPm(:, :, f) = Pm(f)%right_eigenvectors(eos=eos)
         enddo
         do j=i+1-gc, i-1+gc
            do f=1, 2
               if (i==0  .and.f==1) cycle
               if (i==N+1.and.f==2) cycle
               do v=1, 3
                  C(f, j-i, v) = dot_product(LPm(v, :, f), [P(j)%density,                      &
                                                            P(j)%velocity .dot. normal(i+f-1), &
                                                            P(j)%pressure])
               enddo
            enddo
         enddo

         ! compute WENO reconstruction of pseudo charteristic variables
         do v=1, 3
            call solver%interpolator%interpolate(stencil=C(:, :, v), interpolation=CR(:, v))
         enddo

         ! trasform back reconstructed pseudo charteristic variables to primitive ones
         do f=1, 2
            if (i==0  .and.f==1) cycle
            if (i==N+1.and.f==2) cycle
            do v=1, 3
               buffer(v) = dot_product(RPm(v, :, f), CR(f, :))
            enddo
            PR(f, i)%density  = buffer(1)
            PR(f, i)%velocity = buffer(2) * normal(i+f-1) + tangential(f,i)
            PR(f, i)%pressure = buffer(3)
         enddo
      enddo

      ! compute reconstructed conservative variables
      do i=0, N+1
         UR(1, i) = primitive_to_conservative_compressible(primitive=PR(1, i), eos=eos)
         UR(2, i) = primitive_to_conservative_compressible(primitive=PR(2, i), eos=eos)
      enddo
   endselect
   endsubroutine reconstruct_interfaces_characteristic
endmodule off_block_object
