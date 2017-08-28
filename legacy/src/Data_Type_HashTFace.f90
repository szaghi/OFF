!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_HashTFaceDerivedType Data_Type_HashTFace
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_HashTFacePublicProcedure Data_Type_HashTFace
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_HashTFacePrivateProcedure Data_Type_HashTFace
!> @}

!> @brief Module Data_Type_HashTFace contains the definition of Type_HashTFace type and useful procedures for its handling.
!> Type_HashTFace contains the (inter)faces data stored as a dynamic Hierachical Structure. The links data are stored by
!> means of Hash Table. To retrive a specific link data (identified by a unique key, ID) a provided hash function must be used.
!> In order to resolve the "keys collisions" the "chaining" (based on single linked list) technique is used.
module Data_Type_HashTFace
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision     !< Integers and reals precision definition.
USE Data_Type_Face   !< Definition of Type_Face.
USE Data_Type_HashID !< Definition of Type_HashID.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing single linked list data of cell one. The Single Linked List data structure is used to
!> resolve collisions into the hash table of the cells data.
!> @ingroup Data_Type_HashTFaceDerivedType
type, public:: Type_SLL
  type(Type_SLL), pointer::        next => NULL() !< Pointer to the next link of the list.
  type(Type_HashID), allocatable:: ID             !< ID (unique) of the current link node.
  type(Type_Face),   allocatable:: d              !< Link data.
  contains
    procedure:: put =>  put_sll  ! Procedure for inserting a link into the list.
    procedure:: get =>  get_sll  ! Procedure for getting a link from the list.
    procedure:: free => free_sll ! Procedure for freeing (destroyng) the list.
endtype Type_SLL

!> @brief Derived type containing the hash table of cells data. The the collisions of different (unique) IDs into  the same
!> hash buckets are resolved by the "chaining" technique by means of a single linked list data structure.
!> @ingroup Data_Type_HashTFaceDerivedType
type, public:: Type_HashTFace
  type(Type_SLL), allocatable:: ht(:)        !< Hash table.
  integer(I4P)::                leng = 0_I4P !< Lenght of the hash table.
  contains
    procedure:: init => init_hasht ! Procedure for initializing the table.
    procedure:: put =>  put_hasht  ! Procedure for inserting a node into the table.
    procedure:: get =>  get_hasht  ! Procedure for getting a node from the table.
    procedure:: free => free_hasht ! Procedure for freeing (destroyng) the table.
endtype Type_HashTFace
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_HashTFacePrivateProcedure
  !> @{
  !> @brief Recursive subroutine for inserting a link data into the list.
  recursive subroutine put_sll(list,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SLL),   intent(INOUT):: list !< Single linked list.
  type(Type_HashID), intent(IN)::    ID   !< ID (unique) of the current link.
  type(Type_Face),   intent(IN)::    d    !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(list%ID)) then
    if (list%ID/=ID) then
      if (.not.associated(list%next) ) allocate(list%next)
      call put_sll(list=list%next,ID=ID,d=d)
    endif
  else
    if (.not.allocated(list%ID)) allocate(list%ID) ; list%ID = ID
    if (.not.allocated(list%d )) allocate(list%d ) ; list%d  = d
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_sll

  !> @brief Recursive subroutine for getting a link data from the list.
  recursive subroutine get_sll(list,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SLL),              intent(IN)::  list !< Single linked list.
  type(Type_HashID),            intent(IN)::  ID   !< ID (unique) of the current link.
  type(Type_Face), allocatable, intent(OUT):: d    !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(list%ID).and.(list%ID==ID)) then
    if (.not.allocated(d)) allocate(d) ; d = list%d
  elseif(associated(list%next)) then
    call get_sll(list=list%next,ID=ID,d=d)
  else
    if (allocated(d)) deallocate(d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_sll

  !> @brief Recursive subroutine for freeing (destroyng) the list.
  recursive subroutine free_sll(list)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SLL), intent(INOUT) :: list !< Single linked list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(list%next)) then
    call free_sll(list%next)
    deallocate(list%next)
  endif
  list%next => null()
  if (allocated(list%ID)) deallocate(list%ID)
  if (allocated(list%d )) deallocate(list%d )
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_sll

  !> @brief Subroutine for initializing the hash table.
  elemental subroutine init_hasht(htable,leng)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_HashTFace),  intent(INOUT):: htable !< Hash table.
  integer(I4P), optional, intent(IN)::    leng   !< Lenght of the hash table.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(htable%ht)) deallocate(htable%ht)
  if (present(leng)) then
    htable%leng = leng
  else
    htable%leng = ht_def_leng
  endif
  allocate(htable%ht(0:htable%leng-1))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_hasht

  !> @brief Subroutine for inserting a link data into the hash table.
  subroutine put_hasht(htable,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_HashTFace), intent(INOUT):: htable !< Hash table.
  type(Type_HashID),     intent(IN)::    ID     !< ID (unique) of the current link.
  type(Type_Face),       intent(IN)::    d      !< Data of the current link.
  integer(I4P)::                         bucket !< Bucket index of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bucket = ID%hash(tleng=htable%leng)
  call htable%ht(bucket)%put(ID=ID,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_hasht

  !> @brief Subroutine for getting a link data from hash table.
  subroutine get_hasht(htable,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_HashTFace),        intent(IN)::  htable !< Hash table.
  type(Type_HashID),            intent(IN)::  ID     !< ID (unique) of the current link.
  type(Type_Face), allocatable, intent(OUT):: d      !< Data of the current link.
  integer(I4P)::                              bucket !< Bucket index of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bucket = ID%hash(tleng=htable%leng)
  call htable%ht(bucket)%get(ID=ID,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_hasht

  !> @brief Subroutine for freeing (destroyng) the hash table.
  subroutine free_hasht(htable)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_HashTFace), intent(INOUT):: htable  !< Hash table.
  integer(I4P)::                         b,b1,b2 !< Buckets counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(htable%ht)) THEN
    b1 = lbound(htable%ht,dim=1)
    b2 = ubound(htable%ht,dim=1)
    do b=b1,b2
       call htable%ht(b)%free
    enddo
    deallocate(htable%ht)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_hasht
  !> @}
endmodule Data_Type_HashTFace
