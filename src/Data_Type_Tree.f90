!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_TreeDerivedType Data_Type_Tree
!> Module definition of Type_Tree
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_TreePublicProcedure Data_Type_Tree
!> Module definition of Type_Tree
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_TreePrivateProcedure Data_Type_Tree
!> Module definition of Type_Tree
!> @}

!> @brief Module Data_Type_Tree contains the definition of Type_Tree type and useful procedures for its
!> handling.
!> Type_Tree contains generic data stored as a dynamic Hierarchical Data Structure, namely Octree or Quadtree. The data are stored by
!> means of a generic Hash Table structure. To retrieve a specific data (identified by a unique key, ID) a hash function is used.
!> The unique ID is computed by means of Morton filling curve, namely the Z-ordering.
!> In order to resolve the "keys collisions" the "chaining" (based on single linked list) technique is used.
!> @note The hierarchical tree data structure is a dynamic n-tree (typically an octree in 3D or quadtree in 2D). In the following
!> a quadtree is used as an example. Let us assume to have initially 3 nodes (number of ancestor nodes), \f$Na=3\f$. The tree
!> refinement ratio is 4, \f$ref_ratio=4\f$, a quadtree being considered. Let us use the Morton order (or Z-order) to apply
!> a unique ID-key to each node of each level, ID=1,2,3 being the ancestor nodes. Consequently the ID of each child of any nodes can
!> be computed by the equation \f$ ID_{child} = ref_ratio*ID_{node}+i+Na-ref_{ratio}\f$ where \f$i\f$ is the index of child node
!> considered in the children numbering \f$i \in [1,ref_ratio]\f$. Conversely, the parent of a node is computed by the equation
!> \f$ID_{parent} = int(\frac{ID_{node}-1-(Na-ref_{ratio})}{ref_{ratio}})\f$.
!> In the following sketch the tree represented has 3 ancestor nodes and
!> 3 levels of refinement. Node 1 and 3 at level are refined at level 2 where nodes 4 and 13 are refined again at level 3.
!> @code
!>
!>                     +---+          +---+             +---+
!>  Ancestor nodes =>  | 1 |          | 2 |             | 3 |             LEVEL 1: first-last IDs=1-3
!>                     +---+          +---+             +---+
!>                   _/ | | \_                        _/ | | \_
!>                 _/   | |   \_                    _/   | |   \_
!>                /    /   \    \                  /    /   \    \
!>            +---+ +---+ +---+ +---+          +---+ +---+ +---+ +---+
!>            | 4 | | 5 | | 6 | | 7 |          | 12| | 13| | 14| | 15|    LEVEL 2: first-last IDs=4-15
!>            +---+ +---+ +---+ +---+          +---+ +---+ +---+ +---+
!>          _/ | | \_                              _/ | | \_
!>        _/   | |   \_                          _/   | |   \_
!>       /    /   \    \                        /    /   \    \
!>   +---+ +---+ +---+ +---+                +---+ +---+ +---+ +---+
!>   | 16| | 17| | 18| | 19|                | 52| | 53| | 54| | 55|       LEVEL 3: first-last IDs=16-63
!>   +---+ +---+ +---+ +---+                +---+ +---+ +---+ +---+
!> @endcode
!> The previous representation is equivalent to the following quadtree sketch where only leaf nodes are represented.
!> @code
!>   +-------------------+
!>   |         |         |
!>   |         |         |
!>   |   14    |   15    |
!>   |         |         |
!>   |         |         |
!>   |---------3---------|
!>   |         | 54 | 55 |
!>   |         |    |    |
!>   |   12    |---13----|
!>   |         | 52 | 53 |
!>   |         |    |    |
!>   +-------------------+-------------------+
!>   |         |         |                   |
!>   |         |         |                   |
!>   |   6     |   7     |                   |
!>   |         |         |                   |
!>   |         |         |                   |
!>   |---------1---------|        2          |
!>   | 18 | 19 |         |                   |
!>   |    |    |         |                   |
!>   |----4----|   5     |                   |
!>   | 16 | 17 |         |                   |
!>   |    |    |         |                   |
!>   +-------------------+-------------------+
!> @endcode
module Data_Type_Tree
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
integer(I4P), parameter:: ht_leng_def   = 9973_I4P !< Default length of hash table.
integer(I4P), parameter:: ref_ratio_def = 8_I4P    !< Default refinement ratio.
integer(I4P), parameter:: Na_def        = 8_I4P    !< Default ancestor nodes.
!> @brief Derived type containing the node of the Tree. It implements a generic single linked list.
!> @ingroup Data_Type_TreeDerivedType
type, public:: Type_Tree_Node
  private
  integer(I8P),         allocatable:: ID             !< ID (unique) of the current node.
  class(*),             pointer::     d    => null() !< Node data.
  type(Type_Tree_Node), pointer::     next => null() !< Pointer to the next node of the list.
  contains
    procedure:: free         => free_node         ! Procedure for freeing (destroying) the list.
    procedure:: node         => node_node         ! Procedure for returning the ID-th node pointer of the list.
    procedure:: dat          => dat_node          ! Procedure for returning the data pointer of the ID-th node of the list.
    procedure:: put          => put_node          ! Procedure for inserting data into the ID-th node of the list.
    procedure:: get          => get_node          ! Procedure for getting data from the ID-th node of the list.
    procedure:: del          => del_node          ! Procedure for deleting the ID-th node of the list.
    procedure:: length       => length_node       ! Procedure for computing the length of the list.
    procedure:: getIDs       => getIDs_node       ! Procedure for getting the list of actually stored IDs.
    procedure:: print        => print_node        ! Procedure for printing IDs list with a pretty format.
    procedure:: parent_ID    => parent_ID_node    ! Procedure for computing the parent ID.
    procedure:: child_ID     => child_ID_node     ! Procedure for computing the children IDs.
    procedure:: siblings_IDs => siblings_IDs_node ! Procedure for computing the IDs of siblings nodes.
    procedure:: child_number => child_number_node ! Procedure for computing the local number in children numbering [1,ref_ratio].
    procedure:: ref_level    => ref_level_node    ! Procedure for computing the refinement level.
    procedure:: path_IDs     => path_IDs_node     ! Procedure for computing path-IDs from node to its ancestor node (root level).
    procedure:: max_level    => max_level_node    ! Procedure for computing the maximum refinement level of the list.
    final::     finalize_node                     ! Procedure for freeing dynamic memory when finalizing.
    ! private procedures
    procedure:: assign_tree_node                  ! Procedure for assignment.
endtype Type_Tree_Node
!> @brief Derived type defining a generic dynamic Hierarchical Data Structure, namely Octree or Quadtree. The data are stored by
!> means of a generic Hash Table structure. To retrieve a specific data (identified by a unique key, ID) a hash function is used.
!> In order to resolve the "keys collisions" the "chaining" (based on single linked list) technique is used.
!> @ingroup Data_Type_TreeDerivedType
type, public:: Type_Tree
  private
  type(Type_Tree_Node), allocatable:: ht(:)             !< Hash table.
  integer(I4P)::                      leng      = 0_I4P !< Length of the hash table.
  integer(I4P)::                      Na        = 0_I4P !< Number of ancestor nodes (refinement level 1).
  integer(I4P)::                      ref_ratio = 0_I4P !< Refinement ratio.
  integer(I4P)::                      ref_max   = 0_I4P !< Maximum refinement level.
  integer(I8P), allocatable, public:: IDs(:)            !< List of actually stored IDs.
  contains
    procedure:: hash                                ! Procedure for performing the hashing of the ID (unique key).
    procedure:: init          => init_tree          ! Procedure for initializing the tree.
    procedure:: free          => free_tree          ! Procedure for freeing (destroying) the tree.
    procedure:: node          => node_tree          ! Procedure for returning the ID-th node pointer of the tree.
    procedure:: put           => put_tree           ! Procedure for inserting data into node ID-th of the tree.
    procedure:: get           => get_tree           ! Procedure for getting data from node ID-th of the tree.
    procedure:: dat           => dat_tree           ! Procedure for getting data pointer from node ID-th of the tree.
    procedure:: del           => del_tree           ! Procedure for deleting node ID-th of the tree.
    procedure:: length        => length_tree        ! Procedure for computing the length of the tree.
    procedure:: getIDs        => getIDs_tree        ! Procedure for getting the list of actually stored IDs.
    procedure:: update        => update_tree        ! Procedure for updating tree data.
    procedure:: loopID        => loopID_tree        ! Procedure for performing a while loop returning the ID.
    procedure:: print         => print_tree         ! Procedure for printing IDs list with a pretty format.
    procedure:: parent_ID     => parent_ID_tree     ! Procedure for computing the parent ID.
    procedure:: child_ID      => child_ID_tree      ! Procedure for computing the children IDs.
    procedure:: siblings_IDs  => siblings_IDs_tree  ! Procedure for computing the IDs of siblings nodes.
    procedure:: child_number  => child_number_tree  ! Procedure for computing the local number in children numbering [1,ref_ratio].
    procedure:: ref_level     => ref_level_tree     ! Procedure for computing the refinement level.
    procedure:: path_IDs      => path_IDs_tree      ! Procedure for computing path-IDs from node to its ancestor node (root level).
    procedure:: first_ID      => first_ID_tree      ! Procedure for computing first ID of a given level.
    procedure:: last_ID       => last_ID_tree       ! Procedure for computing last  ID of a given level.
    procedure:: get_max_level => get_max_level_tree ! Procedure for computing the maximum refinement level of the list.
    procedure:: str_Na_ref_ratio                    ! Procedure for casting Na and ref_ratio to string.
    final::     finalize_tree                       ! Procedure for freeing dynamic memory when finalizing.
    ! private procedures
    procedure:: assign_tree                         ! Procedure for assignment.
endtype Type_Tree
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_TreePrivateProcedure
  !> @{
  !> @brief Procedure for freeing (destroying) the list.
  recursive subroutine free_node_rec(list)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(INOUT):: list !< List.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(list%next)) then
    call free_node(list=list%next)
    deallocate(list%next)
  endif
  if (allocated( list%ID)) deallocate(list%ID)
  if (associated(list%d )) deallocate(list%d )
  list%d    => null()
  list%next => null()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_node_rec

  !> @brief Procedure for freeing (destroying) the list.
  elemental subroutine free_node(list)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(INOUT):: list !< List.
  class(Type_Tree_Node), pointer::       n,c  !< Pointers for scanning the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do while(associated(list%next))
    c => list%next
    if (associated(c%next)) then
      ! list has at least 3 nodes
      n => c%next
      do while(associated(n))
        c => n
        n => n%next
        do while(associated(n%next))
          c => n
          n => n%next
        enddo
        if (allocated( n%ID)) deallocate(n%ID)
        if (associated(n%d )) deallocate(n%d )
        deallocate(n)
        c%next => null()
      enddo
    else
      ! list has 2 nodes
      if (allocated( c%ID)) deallocate(c%ID)
      if (associated(c%d )) deallocate(c%d )
      deallocate(c)
      list%next => null()
    endif
  enddo
  if (allocated( list%ID)) deallocate(list%ID)
  if (associated(list%d )) deallocate(list%d )
  list%d    => null()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_node

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize_node(list)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tree_Node), intent(INOUT):: list !< List.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call list%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize_node

  !> @brief Procedure for returning the ID-th node pointer of the list.
  !> @note If ID key is not present a null pointer is returned.
  function node_node(list,ID) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), target, intent(IN):: list !< List.
  integer(I8P),                  intent(IN):: ID   !< Unique key of the node of the list to be found.
  type(Type_Tree_Node), pointer::             n    !< Pointer to "ID-th" node of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  n => list
  scan_list: do
    if (allocated(n%ID)) then
      if (n%ID==ID) exit scan_list
    elseif (associated(n%next)) then
      n => n%next
    else
      n => null()
      exit scan_list
    endif
  enddo scan_list
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction node_node

  !> @brief Procedure for returning the data pointer of the ID-th node of the list.
  !> @note If ID key is not present a null pointer is returned.
  function dat_node(list,ID) result(d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), target, intent(IN):: list  !< List.
  integer(I8P),                  intent(IN):: ID    !< Unique key of the node of the list to be found.
  class(*), pointer::                         d     !< Pointer to the data of the "ID-th" node of the list.
  type(Type_Tree_Node), pointer::             n     !< Pointer to "ID-th" node of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  n=>list%node(ID=ID)
  if (associated(n)) then
    d=>n%d
  else
    d=>null()
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dat_node

  !> @brief Procedure for inserting data into the ID-th node of the list. If a node with the provided ID is not present into the
  !> list a new one is created.
  recursive subroutine put_node(list,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(INOUT):: list !< List.
  integer(I8P),          intent(IN)::    ID   !< ID (unique) of the current node.
  class(*),              intent(IN)::    d    !< Data of the ID-th node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(list%ID)) then
    allocate(list%ID,source=ID)
    if (associated(list%d)) deallocate(list%d) ; allocate(list%d,source=d)
    return
  endif
  if (list%ID==ID) then
    if (associated(list%d)) deallocate(list%d) ; allocate(list%d,source=d)
    return
  endif
  ! tail recursion
  if (.not.associated(list%next)) allocate(list%next)
  call put_node(list=list%next,ID=ID,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_node

  !> @brief Procedure for getting data from the ID-th node of the list. If the requested node ID is not present, the output data
  !> is returned in deallocated form.
  recursive subroutine get_node(list,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(IN)::  list !< Llist.
  integer(I8P),          intent(IN)::  ID   !< ID (unique) of the current node.
  class(*), allocatable, intent(OUT):: d    !< Data of the ID-th node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(d)) deallocate(d)
  if (allocated(list%ID)) then
    if (list%ID==ID) then
      allocate(d,source=list%d)
      return
    endif
  endif
  if (associated(list%next)) then
    call get_node(list=list%next,ID=ID,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_node

  !> @brief Procedure for deleting the ID-th node of the list.
  subroutine del_node(list,ID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), target, intent(INOUT):: list !< List.
  integer(I8P),                  intent(IN)::    ID   !< ID (unique) of the node.
  type(Type_Tree_Node), pointer::                c,n  !< Pointers for scanning the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  c => list
  scan_list: do while(associated(c))
    n => c%next
    if (allocated(c%ID)) then
      if (c%ID==ID) then
        if (associated(c%d)) deallocate(c%d)
        if (associated(n)) then
          if (allocated(n%ID)) then
            c%ID = n%ID
          else
            deallocate(c%ID)
          endif
          if (associated(n%d)) then
            allocate(c%d,source=n%d)
          endif
          c%next => n%next
          deallocate(n)
        else
          deallocate(c%ID)
          c%next => null()
        endif
        exit scan_list
      endif
    endif
    c => n
  enddo scan_list
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine del_node

  !> @brief Procedure for computing the length of the list.
  function length_node(list) result(Nn)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(IN):: list !< List.
  integer(I4P)::                      Nn   !< Nodes number.
  type(Type_Tree_Node), pointer::     n    !< Pointer for scanning the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nn = 0 ; if (associated(list%d)) Nn = 1
  n => list%next
  do while (associated(n))
    if (associated(n%d)) then
      Nn = Nn + 1
    endif
    n => n%next
  enddo
  n => null()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction length_node

  !> @brief Procedure for getting the list of actually stored IDs.
  pure subroutine getIDs_node(list,IDs)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(INOUT):: list    !< List.
  integer(I8P),          intent(INOUT):: IDs(1:) !< List of actually stored IDs.
  type(Type_Tree_Node), pointer::        n       !< Pointer for list scanning.
  integer(I4P)::                         i       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(list%ID)) IDs(1)=list%ID
  i=1
  n=>list%next
  do while(associated(n))
    if (allocated(n%ID)) then
      i=i+1
      IDs(i)=n%ID
      n=>n%next
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine getIDs_node

  !> @brief Procedure for printing IDs list with a pretty format.
  subroutine print_node(list,Na,ref_ratio,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), target, intent(IN)::  list               !< List.
  integer(I4P),                  intent(IN)::  Na                 !< Number of ancestors.
  integer(I4P),                  intent(IN)::  ref_ratio          !< Refinement ratio.
  character(*), optional,        intent(IN)::  pref               !< Prefixing string.
  integer(I4P), optional,        intent(OUT):: iostat             !< IO error.
  character(*), optional,        intent(OUT):: iomsg              !< IO error message.
  integer(I4P),                  intent(IN)::  unit               !< Logic unit.
  character(len=:), allocatable::              prefd              !< Prefixing string.
  integer(I4P)::                               iostatd            !< IO error.
  character(500)::                             iomsgd             !< Temporary variable for IO error message.
  integer(I8P), allocatable::                  path(:)            !< Path-IDs from node to root.
  character(len=:), allocatable::              string             !< String containing path-IDs.
  integer(I8P)::                               sib(1:ref_ratio-1) !< Nodes siblings.
  integer(I4P)::                               i                  !< Counter.
  type(Type_Tree_Node), pointer::              n                  !< Pointer for scanning the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  n => list
  do while(allocated(n%ID))
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' ID: '//trim(str(.true.,n%ID))
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Level:            '//&
      trim(str(.true.,n%ref_level(Na=Na,ref_ratio=ref_ratio)))
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Local index:      '//&
      trim(str(.true.,n%child_number(Na=Na,ref_ratio=ref_ratio)))
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Parent-ID:        '//&
      trim(str(.true.,n%parent_ID(Na=Na,ref_ratio=ref_ratio)))
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Children-IDs:     '//&
      trim(str(.true.,n%child_ID(Na=Na,ref_ratio=ref_ratio,i=1)))//'-'//&
      trim(str(.true.,n%child_ID(Na=Na,ref_ratio=ref_ratio,i=ref_ratio)))
    path = n%path_IDs(Na=Na,ref_ratio=ref_ratio)
    string = trim(str(.true.,path(1)))
    if (size(path)>1) then
      do i=2,size(path)
        string = trim(string)//'-'//trim(str(.true.,path(i)))
      enddo
    endif
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Path-IDs-to-root: '//trim(string)
    sib = n%siblings_IDs(Na=Na,ref_ratio=ref_ratio)
    string = trim(str(.true.,sib(1)))
    if (size(sib)>1) then
      do i=2,size(sib)
        string = trim(string)//'-'//trim(str(.true.,sib(i)))
      enddo
    endif
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Siblings-IDs:     '//trim(string)
    if (associated(n%next)) then
      n => n%next
    else
      exit
    endif
  enddo
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_node

  !> @brief Procedure for computing the parent ID.
  elemental function parent_ID_node(node,Na,ref_ratio) result(p)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(IN):: node      !< Node.
  integer(I4P),          intent(IN):: Na        !< Number of ancestors.
  integer(I4P),          intent(IN):: ref_ratio !< Refinement ratio.
  integer(I8P)::                      p         !< Block-ID of block's parent.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  p = int(real(node%ID-1_I8P-(Na-ref_ratio))/ref_ratio,kind=I8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction parent_ID_node

  !> @brief Procedure for computing the children IDs.
  elemental function child_ID_node(node,Na,ref_ratio,i) result(c)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(IN):: node      !< Node.
  integer(I4P),          intent(IN):: Na        !< Number of ancestors.
  integer(I4P),          intent(IN):: ref_ratio !< Refinement ratio.
  integer(I4P),          intent(IN):: i         !< Index of i-th child [1,ref_ratio].
  integer(I8P)::                      c         !< ID of node's i-th child.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  c = ref_ratio*node%ID+i+Na-ref_ratio
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction child_ID_node

  !> @brief Procedure for computing the IDs of siblings nodes.
  pure function siblings_IDs_node(node,Na,ref_ratio) result(sib)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(IN):: node               !< Node.
  integer(I4P),          intent(IN):: Na                 !< Number of ancestors.
  integer(I4P),          intent(IN):: ref_ratio          !< Refinement ratio.
  integer(I8P)::                      sib(1:ref_ratio-1) !< Siblings IDs.
  integer(I4P)::                      local,start        !< Local number of node [1,ref_ratio].
  integer(I4P)::                      i,s                !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  local = node%child_number(Na=Na,ref_ratio=ref_ratio)
  start = node%ID-local+1
  s = 0
  do i = 1,ref_ratio
    if (i/=local) then
      s = s + 1
      sib(s) = start + i - 1
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction siblings_IDs_node

  !> @brief Procedure for computing the local number in the children numbering [1,ref_ratio].
  elemental function child_number_node(node,Na,ref_ratio) result(c)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(IN):: node      !< Node.
  integer(I4P),          intent(IN):: Na        !< Number of ancestors.
  integer(I4P),          intent(IN):: ref_ratio !< Refinement ratio.
  integer(I4P)::                      c         !< Index of child [1,ref_ratio].
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (node%ID<=Na) then
    c = int(node%ID,I4P)
  else
    c = int(node%ID-(Na-ref_ratio) - ((node%ID-1_I8P-(Na-ref_ratio))/ref_ratio)*ref_ratio,kind=I4P)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction child_number_node

  !> @brief Procedure for computing the refinement level.
  elemental function ref_level_node(node,Na,ref_ratio) result(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(IN):: node      !< Node.
  integer(I4P),          intent(IN):: Na        !< Number of ancestors.
  integer(I4P),          intent(IN):: ref_ratio !< Refinement ratio.
  integer(I4P)::                      l         !< Refinement level.
  integer(I8P)::                      i         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  l = 0
  i = node%ID
  do while (i>0)
    l = l + 1
    i = (i-1-(Na-ref_ratio))/ref_ratio
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ref_level_node

  !> @brief Procedure for computing path-IDs through the list from node to its ancestor node (root level).
  pure function path_IDs_node(node,Na,ref_ratio) result(path)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(IN):: node      !< Node.
  integer(I4P),          intent(IN):: Na        !< Number of ancestors.
  integer(I4P),          intent(IN):: ref_ratio !< Refinement ratio.
  integer(I8P), allocatable::         path(:)   !< Path-IDs from node to root.
  integer(I8P), allocatable::         temp(:)   !< Temporary Path-IDs.
  type(Type_Tree_Node)::              n         !< Node counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(node%ID)) then
    n%ID = node%ID
    path = [node%ID]
    do while(n%ref_level(Na=Na,ref_ratio=ref_ratio)>1)
      allocate(temp(1:size(path)+1))
      temp(1:size(path)) = path
      temp(size(path)+1) = n%parent_ID(Na=Na,ref_ratio=ref_ratio)
      call move_alloc(from=temp,to=path)
      n%ID = n%parent_ID(Na=Na,ref_ratio=ref_ratio)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction path_IDs_node

  !> @brief Procedure for computing the maximum refinement level of the list.
  elemental subroutine max_level_node(list,Na,ref_ratio,ref_max)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(INOUT):: list      !< List.
  integer(I4P),          intent(IN)::    Na        !< Number of ancestors.
  integer(I4P),          intent(IN)::    ref_ratio !< Refinement ratio.
  integer(I4P),          intent(OUT)::   ref_max   !< Maximum refinement level.
  type(Type_Tree_Node), pointer::        n         !< Pointer for scanning the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ref_max = 0 ; if (allocated(list%ID)) ref_max = list%ref_level(Na=Na,ref_ratio=ref_ratio)
  n => list%next
  do while (associated(n))
    if (allocated(n%ID)) then
      ref_max = max(ref_max,n%ref_level(Na=Na,ref_ratio=ref_ratio))
    endif
    n => n%next
  enddo
  n => null()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine max_level_node

  ! Assignment (=)
  !> @brief Procedure for assignment between two tree nodes variables.
  elemental subroutine assign_tree_node(node1,node2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree_Node), intent(INOUT):: node1 !< LHS.
  type(Type_Tree_Node),  intent(INOUT):: node2 !< RHS.
  type(Type_Tree_Node), pointer::        n1,n2 !< Pointers for scanning the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call node1%free
  if (allocated( node2%ID)) allocate(node1%ID,source=node2%ID)
  if (associated(node2%d )) allocate(node1%d ,source=node2%d )
  if (associated(node2%next)) then
    allocate(node1%next)
    n1 => node1%next
    n2 => node2%next
    do while (associated(n2))
      if (allocated( n2%ID)) allocate(n1%ID,source=n2%ID)
      if (associated(n2%d )) allocate(n1%d ,source=n2%d )
      n2 => n2%next
      if (associated(n2)) then
        allocate(n1%next)
        n1 => n1%next
      endif
    enddo
  endif
  n1 => null()
  n2 => null()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_tree_node

  !> @brief Procedure for performing the hashing of the ID (unique key).
  elemental function hash(tree,ID) result(bucket)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN):: tree   !< Tree.
  integer(I8P),     intent(IN):: ID     !< ID (unique) of the node.
  integer(I4P)::                 bucket !< Bucket index of the node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bucket = mod(int(ID,I4P),tree%leng)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction hash

  !> @brief Procedure for initializing the Tree.
  elemental subroutine init_tree(tree,source,leng,Na,ref_ratio)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree),          intent(INOUT):: tree      !< Tree.
  type(Type_Tree), optional, intent(IN)::    source    !< Source prototype for initialization.
  integer(I4P),    optional, intent(IN)::    leng      !< Length of the Tree.
  integer(I4P),    optional, intent(IN)::    Na        !< Number of ancestor nodes.
  integer(I4P),    optional, intent(IN)::    ref_ratio !< Refinement ratio.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tree%free
  if (present(source)) then
    tree%leng      = source%leng
    tree%Na        = source%Na
    tree%ref_ratio = source%ref_ratio
  else
    tree%leng      = ht_leng_def   ; if (present(leng     )) tree%leng      = leng
    tree%Na        = Na_def        ; if (present(Na       )) tree%Na        = Na
    tree%ref_ratio = ref_ratio_def ; if (present(ref_ratio)) tree%ref_ratio = ref_ratio
  endif
  allocate(tree%ht(0:tree%leng-1))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_tree

  !> @brief Procedure for freeing (destroying) the Tree.
  elemental subroutine free_tree(tree)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(INOUT):: tree !< Tree.
  integer(I4P)::                    b    !< Bucket counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(tree%ht)) THEN
    do b=lbound(tree%ht,dim=1),ubound(tree%ht,dim=1)
      call tree%ht(b)%free
    enddo
    deallocate(tree%ht)
  endif
  tree%leng = 0_I4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_tree

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize_tree(tree)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tree), intent(INOUT):: tree !< Tree.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tree%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize_tree

  !> @brief Procedure for returning the ID-th node pointer of the tree.
  !> @note If ID key is not present a null pointer is returned.
  function node_tree(tree,ID) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN)::  tree !< Tree.
  integer(I8P),     intent(IN)::  ID   !< Unique key of the node of the tree to be found.
  type(Type_Tree_Node), pointer:: n    !< Pointer to "ID-th" node of the tree.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  n => tree%ht(tree%hash(ID=ID))%node(ID=ID)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction node_tree

  !> @brief Procedure for inserting data into node ID-th of the tree.
  subroutine put_tree(tree,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(INOUT):: tree !< Tree.
  integer(I8P),     intent(IN)::    ID   !< ID (unique) of the node.
  class(*),         intent(IN)::    d    !< Data of the node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tree%ht(tree%hash(ID=ID))%put(ID=ID,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_tree

  !> @brief Procedure for getting data from node ID-th of the tree.
  subroutine get_tree(tree,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree),      intent(IN)::  tree !< Tree.
  integer(I8P),          intent(IN)::  ID   !< ID (unique) of the node.
  class(*), allocatable, intent(OUT):: d    !< Data of the node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tree%ht(tree%hash(ID=ID))%get(ID=ID,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_tree

  !> @brief Procedure for getting data pointer from node ID-th of the tree.
  function dat_tree(tree,ID) result(d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN):: tree !< Tree.
  integer(I8P),     intent(IN):: ID   !< ID (unique) of the node.
  class(*), pointer::            d    !< Data of the node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d => tree%ht(tree%hash(ID=ID))%dat(ID=ID)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dat_tree

  !> @brief Procedure for deleting node ID-th of the tree.
  subroutine del_tree(tree,ID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(INOUT):: tree !< Tree.
  integer(I8P),     intent(IN)::    ID   !< ID (unique) of the node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tree%ht(tree%hash(ID=ID))%del(ID=ID)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine del_tree

  !> @brief Procedure for computing the length of the tree.
  function length_tree(tree) result(Nn)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN):: tree !< Tree.
  integer(I4P)::                 Nn   !< Nodes number.
  integer(I4P)::                 b    !< Bucket counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nn = 0
  if (allocated(tree%ht)) THEN
    do b=lbound(tree%ht,dim=1),ubound(tree%ht,dim=1)
       Nn = Nn + tree%ht(b)%length()
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction length_tree

  !> @brief Procedure for getting the list of actually stored IDs.
  pure subroutine getIDs_tree(tree)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(INOUT):: tree   !< Tree.
  integer(I4P)::                    length !< Actual tree length.
  integer(I4P)::                    b,i    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(tree%ht)) THEN
    length=tree%length()
    if (length>0) then
      if (allocated(tree%IDs)) deallocate(tree%IDs) ; allocate(tree%IDs(1:length))
      i=1
      do b=lbound(tree%ht,dim=1),ubound(tree%ht,dim=1)
        length=tree%ht(b)%length()
        call tree%ht(b)%getIDs(IDs=tree%IDs(i:i+length-1))
        i=i+length
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine getIDs_tree

  !> @brief Procedure for updating tree data.
  pure subroutine update_tree(tree)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(INOUT):: tree !< Tree.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tree%getIDs
  call tree%get_max_level
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update_tree

  !> @brief Procedure for performing a while loop returning the ID value defined into IDs (for ID looping).
  function loopID_tree(tree,ID) result(again)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN)::  tree  !< Tree.
  integer(I8P),     intent(OUT):: ID    !< ID value.
  logical::                       again !< Flag continuing the loop.
  integer(I4P), save::            i = 0 !< ID counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  again = .false.
  if (allocated(tree%IDs)) then
    if (i==0) then
      i = lbound(tree%IDs,dim=1)
      ID = tree%IDs(i)
      again = .true.
    elseif (i<ubound(tree%IDs,dim=1)) then
      i = i + 1
      ID = tree%IDs(i)
      again = .true.
    else
      i = 0
      ID = 0
      again = .false.
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction loopID_tree

  !> @brief Procedure for printing IDs list with a pretty format.
  subroutine print_tree(tree,verbose,objname,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree),       intent(IN)::  tree     !< Tree.
  logical,      optional, intent(IN)::  verbose  !< Flag for activating verbose printing.
  character(*), optional, intent(IN)::  objname  !< Object name (default set to 'generic' nodes).
  character(*), optional, intent(IN)::  pref     !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat   !< IO error.
  character(*), optional, intent(OUT):: iomsg    !< IO error message.
  integer(I4P),           intent(IN)::  unit     !< Logic unit.
  logical::                             verb     !< Flag for activating verbose printing.
  character(len=:), allocatable::       objnamed !< Object name (default set to 'generic' nodes).
  character(len=:), allocatable::       prefd    !< Prefixing string.
  integer(I4P)::                        iostatd  !< IO error.
  character(500)::                      iomsgd   !< Temporary variable for IO error message.
  character(len=:), allocatable::       string   !< Temporary string.
  integer(I4P)::                        b,l      !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  verb = .false. ; if (present(verbose)) verb = verbose
  objnamed = 'nodes' ; if (present(objname)) objnamed = objname
  prefd = '' ; if (present(pref)) prefd = pref
  if (allocated(tree%ht)) then
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Number of ancestor '//trim(objnamed)//' '//&
      trim(str(.true.,tree%Na))
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Refinement ratio '//trim(str(.true.,tree%ref_ratio))
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' MAX refinement level '//trim(str(.true.,tree%ref_max))
    if (tree%ref_max>1) then
      write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' First-Last IDs of each refinement level:'
      do l=1,tree%ref_max
        write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Level('//trim(str(.true.,l))//')='//&
          trim(str(.true.,tree%first_ID(level=l)))//'-'//trim(str(.true.,tree%last_ID(level=l)))
      enddo
    endif
    if (allocated(tree%IDs)) then
      string = ''
      do b=1,size(tree%IDs)-1
        string = string//trim(str(.true.,tree%IDs(b)))//'-'
      enddo
      string = string//trim(str(.true.,tree%IDs(size(tree%IDs))))
      write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Stored '//trim(objnamed)//' IDs '//trim(string)
    endif
    if (verb) then
      write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Stored '//trim(objnamed)//' topology'
      do b=lbound(tree%ht,dim=1),ubound(tree%ht,dim=1)
        call tree%ht(b)%print(Na=tree%Na,ref_ratio=tree%ref_ratio,pref=prefd//"    ",iostat=iostatd,iomsg=iomsgd,unit=unit)
      enddo
    endif
  else
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' The tree of '//trim(objnamed)//' is empty'
  endif
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_tree

  !> @brief Procedure for computing the parent ID.
  elemental function parent_ID_tree(tree,ID) result(p)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN)::  tree      !< Tree.
  integer(I8P),     intent(IN)::  ID        !< ID (unique) of the node.
  integer(I8P)::                  p         !< Block-ID of block's parent.
  type(Type_Tree_Node), pointer:: n         !< Pointer to "ID-th" node of the tree.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  n => tree%node(ID=ID)
  p = -1 ; if (associated(n)) p = n%parent_ID(Na=tree%Na,ref_ratio=tree%ref_ratio)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction parent_ID_tree

  !> @brief Procedure for computing the children IDs.
  elemental function child_ID_tree(tree,ID,i) result(c)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN)::  tree      !< Tree.
  integer(I8P),     intent(IN)::  ID        !< ID (unique) of the node.
  integer(I4P),     intent(IN)::  i         !< Index of i-th child [1,ref_ratio].
  integer(I8P)::                  c         !< ID of node's i-th child.
  type(Type_Tree_Node), pointer:: n         !< Pointer to "ID-th" node of the tree.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  n => tree%node(ID=ID)
  c = -1 ; if (associated(n)) c = n%child_ID(Na=tree%Na,ref_ratio=tree%ref_ratio,i=i)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction child_ID_tree

  !> @brief Procedure for computing the IDs of siblings nodes.
  pure function siblings_IDs_tree(tree,ID) result(sib)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN)::  tree                    !< Tree.
  integer(I8P),     intent(IN)::  ID                      !< ID (unique) of the node.
  integer(I8P)::                  sib(1:tree%ref_ratio-1) !< Siblings IDs.
  type(Type_Tree_Node), pointer:: n                       !< Pointer to "ID-th" node of the tree.
  !----------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  n => tree%node(ID=ID)
  sib = -1 ; if (associated(n)) sib = n%siblings_IDs(Na=tree%Na,ref_ratio=tree%ref_ratio)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction siblings_IDs_tree

  !> @brief Procedure for computing the local number in the children numbering [1,ref_ratio].
  elemental function child_number_tree(tree,ID) result(c)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN)::  tree      !< Tree.
  integer(I8P),     intent(IN)::  ID        !< ID (unique) of the node.
  integer(I4P)::                  c         !< Index of child [1,ref_ratio].
  type(Type_Tree_Node), pointer:: n         !< Pointer to "ID-th" node of the tree.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  n => tree%node(ID=ID)
  c = -1 ; if (associated(n)) c = n%child_number(Na=tree%Na,ref_ratio=tree%ref_ratio)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction child_number_tree

  !> @brief Procedure for computing the refinement level.
  elemental function ref_level_tree(tree,ID) result(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN)::  tree      !< Tree.
  integer(I8P),     intent(IN)::  ID        !< ID (unique) of the node.
  integer(I4P)::                  l         !< Refinement level.
  type(Type_Tree_Node), pointer:: n         !< Pointer to "ID-th" node of the tree.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  n => tree%node(ID=ID)
  l = -1 ; if (associated(n)) l = n%ref_level(Na=tree%Na,ref_ratio=tree%ref_ratio)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ref_level_tree

  !> @brief Procedure for computing path-IDs through the tree from node to its ancestor node (root level).
  pure function path_IDs_tree(tree,ID) result(path)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN)::  tree      !< Tree.
  integer(I8P),     intent(IN)::  ID        !< ID (unique) of the node.
  integer(I8P), allocatable::     path(:)   !< Path-IDs from node to root.
  type(Type_Tree_Node), pointer:: n         !< Pointer to "ID-th" node of the tree.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  n => tree%node(ID=ID)
  if (associated(n)) path = n%path_IDs(Na=tree%Na,ref_ratio=tree%ref_ratio)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction path_IDs_tree

  !> @brief Procedure for computing first ID of a given level.
  elemental function first_ID_tree(tree,level) result(firstID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN):: tree    !< Tree.
  integer(I4P),     intent(IN):: level   !< Refinement level considered.
  integer(I8P)::                 firstID !< First ID of level.
  type(Type_Tree_Node)::         node    !< Temporary node.
  integer(I4P)::                 l       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  firstID = 1_I8P
  if (level>1) then
    do l=2,level
      node%ID = firstID
      firstID = node%child_ID(Na=tree%Na,ref_ratio=tree%ref_ratio,i=1)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction first_ID_tree

  !> @brief Procedure for computing last ID of a given level.
  elemental function last_ID_tree(tree,level) result(lastID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN):: tree   !< Tree.
  integer(I4P),     intent(IN):: level  !< Refinement level considered.
  integer(I8P)::                 lastID !< Last ID of level.
  type(Type_Tree_Node)::         node   !< Temporary node.
  integer(I4P)::                 l      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  lastID = tree%Na
  if (level>1) then
    do l=2,level
      node%ID = lastID
      lastID = node%child_ID(Na=tree%Na,ref_ratio=tree%ref_ratio,i=tree%ref_ratio)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction last_ID_tree

  !> @brief Procedure for computing the maximum refinement level of the tree.
  elemental subroutine get_max_level_tree(tree)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(INOUT):: tree !< Tree.
  integer(I4P)::                    rm   !< Maximum refinement.
  integer(I4P)::                    b    !< Bucket counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  tree%ref_max = 0
  if (allocated(tree%ht)) THEN
    do b=lbound(tree%ht,dim=1),ubound(tree%ht,dim=1)
       call tree%ht(b)%max_level(Na=tree%Na,ref_ratio=tree%ref_ratio,ref_max=rm)
       tree%ref_max = max(tree%ref_max,rm)
       !tree%ref_max = max(tree%ref_max,tree%ht(b)%max_level(Na=tree%Na,ref_ratio=tree%ref_ratio))
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_max_level_tree

  !> @brief Procedure for casting Na and ref_ratio to string.
  elemental function str_Na_ref_ratio(tree) result(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(IN)::     tree   !< Tree.
  character(len=18+DI4P+20+DI4P+1):: string !< Output string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'ancestors_number="'//trim(str(n=tree%Na))//'" refinement_ratio="'//trim(str(n=tree%ref_ratio))//'"'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_Na_ref_ratio

  ! Assignment (=)
  !> @brief Procedure for assignment between two trees variables.
  elemental subroutine assign_tree(tree1,tree2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tree), intent(INOUT):: tree1 !< LHS.
  type(Type_Tree),  intent(INOUT)::    tree2 !< RHS.
  integer(I4P)::                    b     !< Bucket counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tree1%free
  if (allocated(tree2%ht)) then
    allocate(tree1%ht(lbound(tree2%ht,dim=1):ubound(tree2%ht,dim=1)))
    do b=lbound(tree2%ht,dim=1),ubound(tree2%ht,dim=1)
      tree1%ht(b) = tree2%ht(b)
    enddo
  endif
  tree1%leng      = tree2%leng
  tree1%Na        = tree2%Na
  tree1%ref_ratio = tree2%ref_ratio
  tree1%ref_max   = tree2%ref_max
  if (allocated(tree2%IDs)) allocate(tree1%IDs,source=tree2%IDs)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_tree
  !> @}
endmodule Data_Type_Tree
