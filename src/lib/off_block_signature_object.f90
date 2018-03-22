!< OFF block signature object definition and implementation.

module off_block_signature_object
!< OFF block signature object definition and implementation.

use foreseer, only : conservative_compressible, primitive_compressible
use penf, only : I4P, I8P, str
use vecfor, only : vector

implicit none
private
public :: block_signature_object

type :: block_signature_object
   !< Block signature object class.
   !<
   !< Define the block dimensions, id and level.
   integer(I8P)   :: id=0                  !< Unique (Morton) identification code.
   integer(I4P)   :: level=0               !< block refinement level.
   integer(I4P)   :: gc(1:6)=[0,0,0,0,0,0] !< Number of ghost cells along each frame.
   integer(I4P)   :: ni=0                  !< Number of cells in I direction.
   integer(I4P)   :: nj=0                  !< Number of cells in J direction.
   integer(I4P)   :: nk=0                  !< Number of cells in K direction.
   integer(I4P)   :: nc=0                  !< Number of conservative variables.
   integer(I4P)   :: np=0                  !< Number of primitive variables.
   integer(I4P)   :: interfaces_number=0   !< Number of different interfaces (level set).
   type(vector)   :: emin                  !< Coordinates of minimum abscissa (extent) of the block.
   type(vector)   :: emax                  !< Coordinates of maximum abscissa (extent) of the block.
   character(999) :: faces_bc(6)           !< Faces boundary conditions.
   logical        :: is_cartesian=.false.  !< Flag for checking if the block is Cartesian.
   logical        :: is_null_x=.false.     !< Nullify X direction (2D yz, 1D y or z domain).
   logical        :: is_null_y=.false.     !< Nullify Y direction (2D xy, 1D x or y domain).
   logical        :: is_null_z=.false.     !< Nullify Z direction (2D xy, 1D x or y domain).
   contains
      ! public methods
      procedure, pass(self) :: cells_number   !< Return the number of cells.
      procedure, pass(self) :: description    !< Return a pretty-formatted description of block signature.
      procedure, pass(self) :: destroy        !< Destroy block signature.
      procedure, pass(self) :: initialize     !< Initialize block signature.
      procedure, pass(self) :: iolength       !< Return the IO length storage.
      procedure, pass(self) :: load_from_file !< Load block signature from file.
      procedure, pass(self) :: nodes_number   !< Return the number of nodes.
      procedure, pass(self) :: save_into_file !< Save block signature into file.
      ! operators
      generic :: assignment(=) => block_d_assign_block_d !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: block_d_assign_block_d !< Operator `=`.
endtype block_signature_object

contains
   ! public methods
   elemental function cells_number(self, with_ghosts) result(cells_number_)
   !< Return the number of cells.
   class(block_signature_object), intent(in)           :: self          !< Block.
   logical,                       intent(in), optional :: with_ghosts   !< Take into account ghost cells.
   integer(I4P)                                        :: cells_number_ !< Number of cells.
   logical                                             :: with_ghosts_  !< Take into account ghost cells, local variable.

   with_ghosts_ = .true. ; if (present(with_ghosts)) with_ghosts_ = with_ghosts
   if (with_ghosts_) then
      cells_number_ = (self%ni + self%gc(1) + self%gc(2)) &
                    * (self%nj + self%gc(3) + self%gc(4)) &
                    * (self%nk + self%gc(5) + self%gc(6))
   else
      cells_number_ = (self%ni) &
                    * (self%nj) &
                    * (self%nk)
   endif
   endfunction cells_number

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the block signature.
   class(block_signature_object), intent(in)           :: self             !< Block signature object.
   character(*),                  intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                       :: desc             !< Description.
   character(len=:), allocatable                       :: prefix_          !< Prefixing string, local variable.
   integer                                             :: i                !< Counter.
   character(len=1), parameter                         :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//prefix_//'id                 : '//trim(str(self%id,                                 no_sign=.true.))//NL
   desc = desc//prefix_//'level              : '//trim(str(self%level,                              no_sign=.true.))//NL
   desc = desc//prefix_//'gc                 : '//trim(str(self%gc,                                 no_sign=.true.))//NL
   desc = desc//prefix_//'ni                 : '//trim(str(self%ni,                                 no_sign=.true.))//NL
   desc = desc//prefix_//'nj                 : '//trim(str(self%nj,                                 no_sign=.true.))//NL
   desc = desc//prefix_//'nk                 : '//trim(str(self%nk,                                 no_sign=.true.))//NL
   desc = desc//prefix_//'nc                 : '//trim(str(self%nc,                                 no_sign=.true.))//NL
   desc = desc//prefix_//'interfaces number  : '//trim(str(self%interfaces_number,                  no_sign=.true.))//NL
   desc = desc//prefix_//'emin               : '//trim(str([self%emin%x, self%emin%y, self%emin%z]                ))//NL
   desc = desc//prefix_//'emax               : '//trim(str([self%emax%x, self%emax%y, self%emax%z]                ))//NL
   desc = desc//prefix_//'boundary conditions: '//trim(self%faces_bc(1))//' '//trim(self%faces_bc(2))//' '//&
                                                  trim(self%faces_bc(3))//' '//trim(self%faces_bc(4))//' '//&
                                                  trim(self%faces_bc(5))//' '//trim(self%faces_bc(6))//NL
   desc = desc//prefix_//'is cartesian       : '//trim(str(self%is_cartesian                                      ))//NL
   desc = desc//prefix_//'is null X          : '//trim(str(self%is_null_x                                         ))//NL
   desc = desc//prefix_//'is null Y          : '//trim(str(self%is_null_y                                         ))//NL
   desc = desc//prefix_//'is null Z          : '//trim(str(self%is_null_z                                         ))
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy block signature.
   class(block_signature_object), intent(inout) :: self  !< Block signature object.
   type(block_signature_object)                 :: fresh !< Fresh instance of block signature object.

   self = fresh
   endsubroutine destroy

   pure subroutine initialize(self, signature,           &
                              id, level, gc, ni, nj, nk, &
                              interfaces_number,         &
                              emin, emax, is_cartesian, is_null_x, is_null_y, is_null_z, U0, P0)
   !< Initialize block signature.
   !<
   !< @note If both whole `signature` and single components like `id, level, gc...` are passed, the values of
   !< `signature%id, signature%level, ...` are overridden.
   class(block_signature_object),   intent(inout)        :: self              !< Block signature object.
   type(block_signature_object),    intent(in), optional :: signature         !< Block signature input.
   integer(I8P),                    intent(in), optional :: id                !< Unique (Morton) identification code.
   integer(I4P),                    intent(in), optional :: level             !< Grid refinement level.
   integer(I4P),                    intent(in), optional :: gc(1:)            !< Number of ghost cells along each frame.
   integer(I4P),                    intent(in), optional :: ni                !< Number of cells in I direction.
   integer(I4P),                    intent(in), optional :: nj                !< Number of cells in J direction.
   integer(I4P),                    intent(in), optional :: nk                !< Number of cells in K direction.
   integer(I4P),                    intent(in), optional :: interfaces_number !< Number of different interfaces.
   type(vector),                    intent(in), optional :: emin              !< Coordinates of minimum abscissa of the block.
   type(vector),                    intent(in), optional :: emax              !< Coordinates of maximum abscissa of the block.
   logical,                         intent(in), optional :: is_cartesian      !< Flag for checking if the block is Cartesian.
   logical,                         intent(in), optional :: is_null_x         !< Nullify X direction (2D yz, 1D y or z domain).
   logical,                         intent(in), optional :: is_null_y         !< Nullify Y direction (2D xy, 1D x or y domain).
   logical,                         intent(in), optional :: is_null_z         !< Nullify Z direction (2D xy, 1D x or y domain).
   type(conservative_compressible), intent(in), optional :: U0                !< Initial state of conservative variables.
   type(primitive_compressible),    intent(in), optional :: P0                !< Initial state of primitive variables.

   call self%destroy
   if (present(signature        )) self                   = signature
   if (present(id               )) self%id                = id
   if (present(level            )) self%level             = level
   if (present(gc               )) self%gc                = gc
   if (present(ni               )) self%ni                = ni
   if (present(nj               )) self%nj                = nj
   if (present(nk               )) self%nk                = nk
   if (present(interfaces_number)) self%interfaces_number = interfaces_number
   if (present(emin             )) self%emin              = emin
   if (present(emax             )) self%emax              = emax
   if (present(is_cartesian     )) self%is_cartesian      = is_cartesian
   if (present(is_null_x        )) self%is_null_x         = is_null_x
   if (present(is_null_y        )) self%is_null_y         = is_null_y
   if (present(is_null_z        )) self%is_null_z         = is_null_z
   if (present(U0               )) self%nc                = size(U0%array(), dim=1)
   if (present(P0               )) self%np                = size(P0%array(), dim=1)
   endsubroutine initialize

   function iolength(self)
   !< Return the IO length storage.
   class(block_signature_object), intent(in) :: self     !< Block signature object.
   integer(I4P)                              :: iolength !< IO length storage.

   inquire(iolength=iolength) self%id,                               &
                              self%level,                            &
                              self%gc,                               &
                              self%ni,                               &
                              self%nj,                               &
                              self%nk,                               &
                              self%nc,                               &
                              self%np,                               &
                              self%interfaces_number,                &
                              self%emin%x, self%emin%y, self%emin%z, &
                              self%emax%x, self%emax%y, self%emax%z, &
                              self%is_cartesian,                     &
                              self%is_null_x , self%is_null_y , self%is_null_z
   endfunction iolength

   subroutine load_from_file(self, file_unit)
   !< Load block signature from file.
   class(block_signature_object), intent(inout) :: self      !< Block signature object.
   integer(I4P),                  intent(in)    :: file_unit !< File unit.

   read(unit=file_unit) self%id,                               &
                        self%level,                            &
                        self%gc,                               &
                        self%ni,                               &
                        self%nj,                               &
                        self%nk,                               &
                        self%nc,                               &
                        self%np,                               &
                        self%interfaces_number,                &
                        self%emin%x, self%emin%y, self%emin%z, &
                        self%emax%x, self%emax%y, self%emax%z, &
                        self%is_cartesian,                     &
                        self%is_null_x , self%is_null_y , self%is_null_z
   endsubroutine load_from_file

   elemental function nodes_number(self, with_ghosts) result(nodes_number_)
   !< Return the number of nodes.
   class(block_signature_object), intent(in)           :: self          !< Block.
   logical,                       intent(in), optional :: with_ghosts   !< Take into account ghost cells.
   integer(I4P)                                        :: nodes_number_ !< Number of nodes.
   logical                                             :: with_ghosts_  !< Take into account ghost cells, local variable.

   with_ghosts_ = .true. ; if (present(with_ghosts)) with_ghosts_ = with_ghosts
   if (with_ghosts_) then
      nodes_number_ = (self%ni + self%gc(1) + self%gc(2) + 1) &
                    * (self%nj + self%gc(3) + self%gc(4) + 1) &
                    * (self%nk + self%gc(5) + self%gc(6) + 1)
   else
      nodes_number_ = (self%ni + 1) &
                    * (self%nj + 1) &
                    * (self%nk + 1)
   endif
   endfunction nodes_number

   subroutine save_into_file(self, file_unit)
   !< Load the block signature of all blocks from file.
   class(block_signature_object), intent(in) :: self      !< Block signature object.
   integer(I4P),                  intent(in) :: file_unit !< File unit.

   write(unit=file_unit) self%id,                               &
                         self%level,                            &
                         self%gc,                               &
                         self%ni,                               &
                         self%nj,                               &
                         self%nk,                               &
                         self%nc,                               &
                         self%np,                               &
                         self%interfaces_number,                &
                         self%emin%x, self%emin%y, self%emin%z, &
                         self%emax%x, self%emax%y, self%emax%z, &
                         self%is_cartesian,                     &
                         self%is_null_x , self%is_null_y , self%is_null_z
   endsubroutine save_into_file

   ! private methods
   pure subroutine block_d_assign_block_d(lhs, rhs)
   !< Operator `=`.
   class(block_signature_object), intent(inout) :: lhs !< Left hand side.
   type(block_signature_object),  intent(in)    :: rhs !< Right hand side.

   lhs%id                = rhs%id
   lhs%level             = rhs%level
   lhs%gc                = rhs%gc
   lhs%ni                = rhs%ni
   lhs%nj                = rhs%nj
   lhs%nk                = rhs%nk
   lhs%nc                = rhs%nc
   lhs%np                = rhs%np
   lhs%interfaces_number = rhs%interfaces_number
   lhs%emin              = rhs%emin
   lhs%emax              = rhs%emax
   lhs%faces_bc          = rhs%faces_bc
   lhs%is_cartesian      = rhs%is_cartesian
   lhs%is_null_x         = rhs%is_null_x
   lhs%is_null_y         = rhs%is_null_y
   lhs%is_null_z         = rhs%is_null_z
   endsubroutine block_d_assign_block_d
endmodule off_block_signature_object
