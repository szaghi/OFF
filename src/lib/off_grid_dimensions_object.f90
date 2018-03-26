!< OFF grid dimensions object definition and implementation.

module off_grid_dimensions_object
!< OFF grid dimensions object definition and implementation.

use off_block_signature_object, only : block_signature_object
use flow, only : conservative_compressible
use penf, only : I4P, I8P
use vecfor, only : vector

implicit none
private
public :: grid_dimensions_object

type :: grid_dimensions_object
   !< Grid dimensions object class.
   integer(I4P)                              :: blocks_number=0    !< Number of blocks, blobal (whole) number on all process/image.
   type(block_signature_object), allocatable :: block_signature(:) !< Signature of each block.
   contains
      ! public methods
      procedure, pass(self) :: alloc                     !< Allocate blocks.
      procedure, pass(self) :: description               !< Return a pretty-formatted description of grid dimensions.
      procedure, pass(self) :: destroy                   !< Destroy grid dimensions.
      procedure, pass(self) :: initialize                !< Initialize grid dimensions.
      procedure, pass(self) :: iolength                  !< Return the IO length storage.
      procedure, pass(self) :: iopos_block_nodes         !< Return the IO position where nodes of block b-th are stored.
      procedure, pass(self) :: iopos_block_conservatives !< Return the IO position where conservatives of block b-th are stored.
      procedure, pass(self) :: load_from_file            !< Load grid dimensions from file.
      procedure, pass(self) :: save_into_file            !< Save grid dimensions into file.
      ! operators
      generic :: assignment(=) => grid_d_assign_grid_d !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: grid_d_assign_grid_d !< Operator `=`.
endtype grid_dimensions_object

contains
   ! public methods
   elemental subroutine alloc(self, blocks_number)
   !< Allocate blocks.
   class(grid_dimensions_object), intent(inout)        :: self          !< Grid dimensions object.
   integer(I4P),                  intent(in), optional :: blocks_number !< Number of blocks.

   if (.not.allocated(self%block_signature)) then
      if (present(blocks_number)) self%blocks_number = blocks_number
      if (self%blocks_number > 0) allocate(self%block_signature(1:self%blocks_number))
   endif
   endsubroutine alloc

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the grid dimensions.
   class(grid_dimensions_object), intent(in)           :: self             !< Grid dimensions object.
   character(*),                  intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                       :: desc             !< Description.
   character(len=1), parameter                         :: NL=new_line('a') !< New line character.
   integer(I4P)                                        :: b                !< Counter.

   desc = ''
   if (self%blocks_number > 0) then
      do b=1, self%blocks_number - 1
         desc = desc//self%block_signature(b)%description(prefix=prefix)//NL//NL
      enddo
      desc = desc//self%block_signature(self%blocks_number)%description(prefix=prefix)
   endif
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy grid dimensions.
   class(grid_dimensions_object), intent(inout) :: self  !< Grid dimensions object.
   type(grid_dimensions_object)                 :: fresh !< Fresh instance of grid dimensions object.

   self = fresh
   if (allocated(self%block_signature)) then
      call self%block_signature%destroy
      deallocate(self%block_signature)
   endif
   endsubroutine destroy

   pure subroutine initialize(self, block_signature)
   !< Initialize grid dimensions.
   class(grid_dimensions_object), intent(inout)        :: self                !< Grid dimensions object.
   type(block_signature_object),  intent(in), optional :: block_signature(1:) !< Dimensions of each block.

   call self%destroy
   if (present(block_signature)) then
      self%blocks_number = size(block_signature, dim=1)
      allocate(self%block_signature(1:self%blocks_number), source=block_signature)
   endif
   endsubroutine initialize

   function iolength(self)
   !< Return the IO length storage.
   class(grid_dimensions_object), intent(in) :: self     !< Grid dimensions object.
   integer(I4P)                              :: iolength !< IO length storage.

   inquire(iolength=iolength) self%blocks_number
   if (self%blocks_number > 0) iolength = iolength + self%blocks_number * self%block_signature(1)%iolength()
   endfunction iolength

   function iopos_block_nodes(self, b)
   !< Return the IO position where nodes of block b-th are stored.
   class(grid_dimensions_object), intent(in) :: self              !< Grid dimensions object.
   integer(I4P),                  intent(in) :: b                 !< Block index.
   integer(I4P)                              :: iopos_block_nodes !< IO position where nodes of block b-th are stored.
   type(vector)                              :: node              !< A node coordinate prototype.
   integer(I4P)                              :: node_iolength     !< Node IO length storage.
   integer(I4P)                              :: bb                !< Counter.

   inquire(iolength=node_iolength) node%x, node%y, node%z
   iopos_block_nodes = self%iolength() ! file header length
   if (self%blocks_number > 0) then
      if (b > 1.and. b <= self%blocks_number) then
         ! b-th block, there are the file header and [1:b-1] blocks before its nodes
         do bb=1, b-1
            iopos_block_nodes = iopos_block_nodes + self%block_signature(bb)%nodes_number() * node_iolength
         enddo
         iopos_block_nodes = iopos_block_nodes + 1
      else
         ! first block, there is only the file header before its nodes
         iopos_block_nodes = iopos_block_nodes + 1
      endif
   endif
   endfunction iopos_block_nodes

   function iopos_block_conservatives(self, b) result(iopos)
   !< Return the IO position where conservatives of block b-th are stored.
   class(grid_dimensions_object), intent(in) :: self          !< Grid dimensions object.
   integer(I4P),                  intent(in) :: b             !< Block index.
   integer(I4P)                              :: iopos         !< IO position where conservatives of block b-th are stored.
   type(conservative_compressible)           :: cons          !< A conservatives prototype.
   integer(I4P)                              :: cons_iolength !< Node IO length storage.
   integer(I4P)                              :: bb            !< Counter.

   inquire(iolength=cons_iolength) cons%array()
   iopos = self%iolength() ! file header length
   if (self%blocks_number > 0) then
      if (b > 1.and. b <= self%blocks_number) then
         ! b-th block, there are the file header and [1:b-1] blocks before its nodes
         do bb=1, b-1
            iopos = iopos + self%block_signature(bb)%cells_number() * cons_iolength
         enddo
         iopos = iopos + 1
      else
         ! first block, there is only the file header before its nodes
         iopos = iopos + 1
      endif
   endif
   endfunction iopos_block_conservatives

   subroutine load_from_file(self, file_unit)
   !< Load grid dimensions from file.
   class(grid_dimensions_object), intent(inout) :: self      !< Grid dimensions object.
   integer(I4P),                  intent(in)    :: file_unit !< File unit.
   integer(I4P)                                 :: b         !< Counter.

   call self%destroy
   read(unit=file_unit) self%blocks_number
   if (self%blocks_number > 0) then
      call self%alloc
      do b=1, self%blocks_number
         call self%block_signature(b)%load_from_file(file_unit=file_unit)
      enddo
   endif
   endsubroutine load_from_file

   subroutine save_into_file(self, file_unit)
   !< Load the grid dimensions of all blocks from file.
   class(grid_dimensions_object), intent(in) :: self      !< Grid dimensions object.
   integer(I4P),                  intent(in) :: file_unit !< File unit.
   integer(I4P)                              :: b         !< Counter.

   if (self%blocks_number > 0) then
      write(unit=file_unit) self%blocks_number
      do b=1, self%blocks_number
         call self%block_signature(b)%save_into_file(file_unit=file_unit)
      enddo
   endif
   endsubroutine save_into_file

   ! private methods
   pure subroutine grid_d_assign_grid_d(lhs, rhs)
   !< Operator `=`.
   class(grid_dimensions_object), intent(inout) :: lhs !< Left hand side.
   type(grid_dimensions_object),  intent(in)    :: rhs !< Right hand side.

   if (allocated(rhs%block_signature)) then
      call lhs%destroy
      lhs%blocks_number = rhs%blocks_number
      allocate(lhs%block_signature, source=rhs%block_signature)
   endif
   endsubroutine grid_d_assign_grid_d
endmodule off_grid_dimensions_object
