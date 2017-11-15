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

use, intrinsic :: iso_c_binding, only : c_ptr
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_block_signature_object, only : block_signature_object
use off_cell_object, only : cell_object
use off_error_object, only : error_object
use off_face_object, only : face_object
use off_node_object, only : node_object
use cgal_polyhedra, only : cgal_polyhedron_closest, cgal_polyhedron_finalize, cgal_polyhedron_inside, cgal_polyhedron_read
use penf, only : FR8P, FI4P, I1P, I4P, I8P, R8P
use vecfor, only : vector, ex, ey, ez
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
   contains
      ! public methods
      procedure, pass(self) :: cells_number              !< Return the number of cells.
      procedure, pass(self) :: compute_space_operator    !< Compute space operator.
      procedure, pass(self) :: create_linspace           !< Create a Cartesian block with linearly spaced nodes.
      procedure, pass(self) :: destroy                   !< Destroy block.
      procedure, pass(self) :: immerge_off_geometry      !< "Immerge" geometry (described into a OFF file) into the block grid.
      procedure, pass(self) :: interpolate_at_nodes      !< Interpolate cell-centered variable at nodes.
      procedure, pass(self) :: initialize                !< Initialize block.
      procedure, pass(self) :: load_nodes_from_file      !< Load nodes from file.
      procedure, pass(self) :: nodes_number              !< Return the number of nodes.
      procedure, pass(self) :: save_file_grid            !< Save gird file.
      procedure, pass(self) :: save_nodes_into_file      !< Save nodes into file.
      procedure, pass(self) :: update_level_set_distance !< Update level set distance.
      ! operators
      generic :: assignment(=) => block_assign_block !< Overload `=`.
      ! private methods
      procedure, pass(lhs),  private :: block_assign_block    !< Operator `=`.
      procedure, pass(self), private :: compute_cells_center  !< Compute cells center.
      procedure, pass(self), private :: compute_extents       !< Compute block extents.
      procedure, pass(self), private :: compute_faces_metrics !< Compute block faces metrics.
      procedure, pass(self), private :: compute_metrics       !< Compute block metrics.
      procedure, pass(self), private :: compute_volumes       !< Compute block volumes.
      procedure, pass(self), private :: correct_metrics       !< Correct block metrics.
      procedure, pass(self), private :: node_to_center        !< Compute cell centers coordinates from cell nodes.
      procedure, pass(self), private :: nullify_normals       !< Nullify normals for 2D or 1D domains.
      procedure, pass(self), private :: save_file_grid_tec    !< Save grid file in Tecplot format.
      procedure, pass(self), private :: save_file_grid_vtk    !< Save grid file in VTK format.
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

   elemental subroutine destroy(self)
   !< Destroy block.
   class(block_object), intent(inout) :: self  !< Block.
   type(block_object)                 :: fresh !< Fresh instance of block object.

   self = fresh
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
   endsubroutine destroy

   subroutine immerge_off_geometry(self, file_name, n, outside_reference)
   !< "Immerge" geometry (described into a OFF file) into the block grid.
   !<
   !< @note OFF file format (Object File Format) is a geometry definition file format containing the description
   !< of the composing polygons of the 3D object used by CGAL.
   class(block_object), intent(inout) :: self              !< Block.
   character(*),        intent(in)    :: file_name         !< Name of OFF file.
   integer(I4P),        intent(in)    :: n                 !< Number of geometry in the global numbering.
   type(vector),        intent(in)    :: outside_reference !< A reference point outside the body.
   type(c_ptr)                        :: geometry_ptr      !< Geometry tree pointer.
   type(vector)                       :: closest           !< Closest point coordinates.
   logical                            :: is_inside         !< Logical dummy.
   integer(I4P)                       :: i, j, k           !< Counter.

   if (allocated(self%cell)) then
      call cgal_polyhedron_read(ptree=geometry_ptr, fname=trim(adjustl(file_name)))
      do k=1, self%signature%nk
         do j=1, self%signature%nj
            do i=1, self%signature%ni
               call cgal_polyhedron_closest(ptree=geometry_ptr,           &
                                            xq=self%cell(i,j,k)%center%x, &
                                            yq=self%cell(i,j,k)%center%y, &
                                            zq=self%cell(i,j,k)%center%z, &
                                            xn=closest%x, yn=closest%y, zn=closest%z)
               self%cell(i,j,k)%level_set%distances(n) = sqrt((self%cell(i,j,k)%center%x - closest%x)**2 + &
                                                              (self%cell(i,j,k)%center%y - closest%y)**2 + &
                                                              (self%cell(i,j,k)%center%z - closest%z)**2)
               is_inside = cgal_polyhedron_inside(ptree=geometry_ptr,           &
                                                  xq=self%cell(i,j,k)%center%x, &
                                                  yq=self%cell(i,j,k)%center%y, &
                                                  zq=self%cell(i,j,k)%center%z, &
                                                  xr=outside_reference%x,       &
                                                  yr=outside_reference%y,       &
                                                  zr=outside_reference%z)
               if (is_inside) self%cell(i,j,k)%level_set%distances(n) = -self%cell(i,j,k)%level_set%distances(n)
            enddo
         enddo
      enddo
   endif
   endsubroutine immerge_off_geometry

   subroutine initialize(self, signature,                                           &
                         id, level, gc, ni, nj, nk,                                 &
                         emin, emax, is_cartesian, is_null_x, is_null_y, is_null_z, &
                         interfaces_number, distances)
   !< Initialize block.
   !<
   !< Assign block signature, allocate dynamic memory and set block features.
   class(block_object),          intent(inout)        :: self              !< Block.
   type(block_signature_object), intent(in), optional :: signature         !< Signature, namely id, level, dimensions, etc...
   integer(I8P),                 intent(in), optional :: id                !< Unique (Morton) identification code.
   integer(I4P),                 intent(in), optional :: level             !< Grid refinement level.
   integer(I4P),                 intent(in), optional :: gc(1:)            !< Number of ghost cells along each frame.
   integer(I4P),                 intent(in), optional :: ni                !< Number of cells in I direction.
   integer(I4P),                 intent(in), optional :: nj                !< Number of cells in J direction.
   integer(I4P),                 intent(in), optional :: nk                !< Number of cells in K direction.
   type(vector),                 intent(in), optional :: emin              !< Coordinates of minimum abscissa of the block.
   type(vector),                 intent(in), optional :: emax              !< Coordinates of maximum abscissa of the block.
   logical,                      intent(in), optional :: is_cartesian      !< Flag for checking if the block is Cartesian.
   logical,                      intent(in), optional :: is_null_x         !< Nullify X direction (2D yz, 1D y or z domain).
   logical,                      intent(in), optional :: is_null_y         !< Nullify Y direction (2D xy, 1D x or y domain).
   logical,                      intent(in), optional :: is_null_z         !< Nullify Z direction (2D xy, 1D x or y domain).
   integer(I4P),                 intent(in), optional :: interfaces_number !< Number of different interfaces.
   real(R8P),                    intent(in), optional :: distances(:)      !< Distance from all interfaces.
   integer(I4P)                                       :: i                 !< Counter.
   integer(I4P)                                       :: j                 !< Counter.
   integer(I4P)                                       :: k                 !< Counter.

   self%error%status = ERROR_BLOCK_CREATE_FAILED

   if ((.not.(present(signature))).and.&
       (.not.(present(id).and.present(level).and.present(gc).and.present(ni).and.present(nj).and.present(nk)))) then
      error stop 'error: either signature or (id, level, gc, ni, nj, nk) tuple must be passed to a block initialized'
   endif

   call self%destroy

   call self%signature%initialize(signature=signature,                             &
                                  id=id, level=level, gc=gc, ni=ni, nj=nj, nk=nk,  &
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
               call self%cell(i,j,k)%initialize(interfaces_number=interfaces_number, distances=distances)
            enddo
         enddo
      enddo
   endassociate

   call self%face_i%initialize
   call self%face_j%initialize
   call self%face_k%initialize

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

   subroutine save_file_grid(self, file_name, ascii, metrics, tecplot, vtk)
   !< Save grid file file.
   class(block_object), intent(inout)        :: self         !< Block.
   character(*),        intent(in)           :: file_name    !< File name.
   logical,             intent(in), optional :: ascii        !< Ascii/binary output.
   logical,             intent(in), optional :: metrics      !< Save also metrics data.
   logical,             intent(in), optional :: tecplot      !< Tecplot output format sentinel.
   logical,             intent(in), optional :: vtk          !< VTK output format sentinel.
   logical                                   :: tecplot_     !< Tecplot format sentinel, local variable.
   logical                                   :: vtk_         !< VTK format sentinel, local variable.

   tecplot_ = .false. ; if (present(tecplot)) tecplot_ = tecplot
   vtk_     = .false. ; if (present(vtk    )) vtk_     = vtk

   if (vtk_) call self%save_file_grid_vtk(file_name=file_name, ascii=ascii, metrics=metrics)
   endsubroutine save_file_grid

   subroutine save_nodes_into_file(self, file_unit, pos)
   !< Save nodes into file.
   class(block_object), intent(inout) :: self      !< Block.
   integer(I4P),        intent(in)    :: file_unit !< File unit.
   integer(I4P),        intent(in)    :: pos       !< Position to start the loading.

   write(file_unit, pos=pos, iostat=self%error%status) self%node%vertex%x, self%node%vertex%y, self%node%vertex%z
   endsubroutine save_nodes_into_file

   elemental subroutine update_level_set_distance(self)
   !< Update level set distance.
   class(block_object), intent(inout) :: self    !< Block.
   integer(I4P)                       :: i, j, k !< Counter.

   if (allocated(self%cell)) then
      do k=1, self%signature%nk
         do j=1, self%signature%nj
            do i=1, self%signature%ni
               call self%cell(i,j,k)%level_set%update_distance
            enddo
         enddo
      enddo
   endif
   endsubroutine update_level_set_distance

   ! private methods
   pure subroutine block_assign_block(lhs, rhs)
   !< Operator `=`.
   class(block_object), intent(inout) :: lhs !< Left hand side.
   type(block_object),  intent(in)    :: rhs !< Right hand side.

   lhs%error        = rhs%error
   lhs%signature    = rhs%signature
   if (allocated(rhs%cell)) then
      if (allocated(lhs%cell)) then
         call lhs%cell%destroy
         deallocate(lhs%cell)
      endif
      lhs%cell = rhs%cell
   endif
   if (allocated(rhs%face_i)) then
      if (allocated(lhs%face_i)) then
         call lhs%face_i%destroy
         deallocate(lhs%face_i)
      endif
      lhs%face_i = rhs%face_i
   endif
   if (allocated(rhs%face_j)) then
      if (allocated(lhs%face_j)) then
         call lhs%face_j%destroy
         deallocate(lhs%face_j)
      endif
      lhs%face_j = rhs%face_j
   endif
   if (allocated(rhs%face_k)) then
      if (allocated(lhs%face_k)) then
         call lhs%face_k%destroy
         deallocate(lhs%face_k)
      endif
      lhs%face_k = rhs%face_k
   endif
   if (allocated(rhs%node)) then
      if (allocated(lhs%node)) then
         call lhs%node%destroy
         deallocate(lhs%node)
      endif
      lhs%node = rhs%node
   endif
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

   associate(node=>self%node, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk)
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

      do k=1, nk
         do j=1, nj
            do i=0, ni
               call self%face_i(i, j, k)%compute_metrics(pt1 = node(i, j-1, k-1)%vertex, &
                                                         pt2 = node(i, j  , k-1)%vertex, &
                                                         pt3 = node(i, j  , k  )%vertex, &
                                                         pt4 = node(i, j-1, k  )%vertex, signd=signi)
            enddo
         enddo
      enddo

      do k=1, nk
         do j=0, nj
            do i=1, ni
               call self%face_j(i, j, k)%compute_metrics(pt1 = node(i-1, j, k-1)%vertex, &
                                                         pt2 = node(i-1, j, k  )%vertex, &
                                                         pt3 = node(i  , j, k  )%vertex, &
                                                         pt4 = node(i  , j, k-1)%vertex, signd=signj)
            enddo
         enddo
      enddo

      do k=0, nk
         do j=1, nj
            do i=1, ni
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

   associate(node=>self%node, ni=>self%signature%ni, nj=>self%signature%nj, nk=>self%signature%nk)
      do k=1, nk
         do j=1, nj
            do i=1, ni
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
   !< Save grid file in Tecplot format.
   !<
   !< @TODO implement Tecplot output.
   class(block_object), intent(inout)        :: self      !< Block.
   character(*),        intent(in)           :: file_name !< Output file name.
   logical,             intent(in), optional :: ascii     !< Ascii/binary output.
   logical,             intent(in), optional :: metrics   !< Save also metrics data.

   error stop 'error: block mesh Tecplot output to be implemented'
   endsubroutine save_file_grid_tec

   subroutine save_file_grid_vtk(self, file_name, ascii, metrics)
   !< Save mesh data into VTK file.
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
         self%error%status = vtk%xml_writer%write_dataarray(location='cell', action='close')
      endif
      self%error%status = vtk%xml_writer%write_piece()
      self%error%status = vtk%finalize()
   endassociate
   endsubroutine save_file_grid_vtk
endmodule off_block_object
