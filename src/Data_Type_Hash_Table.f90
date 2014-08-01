!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Hash_TableDerivedType Data_Type_Hash_Table
!> Module definition of Type_Hash_Table
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_Hash_TablePublicProcedure Data_Type_Hash_Table
!> Module definition of Type_Hash_Table
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Hash_TablePrivateProcedure Data_Type_Hash_Table
!> Module definition of Type_Hash_Table
!> @}

!> @brief Module Data_Type_Hash_Table contains the definition of Type_Hash_Table type and useful procedures for its
!> handling.
!> Type_Hash_Table contains generic data stored as a dynamic Hierachical Structure. The data are stored by
!> means of a generic Hash Table structure. To retrive a specific data (identified by a unique key, ID) a hash function is used.
!> In order to resolve the "keys collisions" the "chaining" (based on single linked list) technique is used.
module Data_Type_Hash_Table
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                    ! Integers and reals precision definition.
USE Data_Type_Generic_List, only: Type_Generic_List ! Definition of Type_Generic_List.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
integer(I8P), parameter:: ht_def_leng = 9973_I4P !< Default lenght of hash table.
!> @brief Derived type defining a generic hash table data structure.vThe the collisions of different (unique) IDs into the same
!> hash buckets are resolved by the "chaining" technique by means of a single linked list data structure.
!> @ingroup Data_Type_Hash_TableDerivedType
type, public:: Type_Hash_Table
  private
  type(Type_Generic_List), allocatable:: ht(:)        !< Hash table.
  integer(I4P)::                         leng = 0_I4P !< Lenght of the hash table.
  integer(I8P), allocatable::            IDs(:)       !< List of actually stored IDs.
  contains
    procedure:: hash     ! Procedure defining the hash function.
    procedure:: init     ! Procedure for initializing the table.
    procedure:: free     ! Procedure for freeing (destroying) the table.
    procedure:: put      ! Procedure for inserting data into node ID-th of the table.
    procedure:: get      ! Procedure for getting data from node ID-th of the table.
    procedure:: dat      ! Procedure for getting data pointer from node ID-th of the table.
    procedure:: del      ! Procedure for deleting node ID-th of the table.
    procedure:: length   ! Procedure for computing the length of the table.
    procedure:: getIDs   ! Procedure for getting the list of actually stored IDs.
    procedure:: print    ! Procedure for printing IDs list with a pretty format.
    final::     finalize ! Procedure for freeing dynamic memory when finalizing.
endtype Type_Hash_Table
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Hash_TablePrivateProcedure
  !> @{
  !> @brief
  elemental function hash(table,ID) result(bucket)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Hash_Table), intent(IN):: table  !< Hash table.
  integer(I8P),           intent(IN):: ID     !< ID (unique) of the node.
  integer(I4P)::                       bucket !< Bucket index of the node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bucket = mod(int(ID,I4P),table%leng)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction hash

  !> @brief Procedure for initializing the hash table.
  subroutine init(table,leng)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Hash_Table), intent(INOUT):: table !< Hash table.
  integer(I4P), optional, intent(IN)::    leng  !< Lenght of the hash table.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call table%free
  if (present(leng)) then
    table%leng = leng
  else
    table%leng = ht_def_leng
  endif
  allocate(table%ht(0:table%leng-1))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  !> @brief Procedure for freeing (destroying) the hash table.
  subroutine free(table)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Hash_Table), intent(INOUT):: table !< Hash table.
  integer(I4P)::                          b     !< Bucket counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(table%ht)) THEN
    do b=lbound(table%ht,dim=1),ubound(table%ht,dim=1)
      call table%ht(b)%free
    enddo
    deallocate(table%ht)
  endif
  table%leng = 0_I4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Procedure for freeing dynamic memory when finalizing.
  subroutine finalize(table)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Hash_Table), intent(INOUT):: table !< Hash table.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call table%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for inserting data into node ID-th of the table.
  subroutine put(table,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Hash_Table), intent(INOUT):: table !< Hash table.
  integer(I8P),           intent(IN)::    ID    !< ID (unique) of the node.
  class(*),               intent(IN)::    d     !< Data of the node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call table%ht(table%hash(ID=ID))%put(ID=ID,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put

  !> @brief Procedure for getting data from node ID-th of the table.
  subroutine get(table,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Hash_Table), intent(IN)::  table !< Hash table.
  integer(I8P),           intent(IN)::  ID    !< ID (unique) of the node.
  class(*), allocatable,  intent(OUT):: d     !< Data of the node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call table%ht(table%hash(ID=ID))%get(ID=ID,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get

  !> @brief Procedure for getting data pointer from node ID-th of the table.
  function dat(table,ID) result(d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Hash_Table), intent(IN)::  table !< Hash table.
  integer(I8P),           intent(IN)::  ID    !< ID (unique) of the node.
  class(*), pointer::                   d     !< Data of the node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d => table%ht(table%hash(ID=ID))%dat(ID=ID)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dat

  !> @brief Procedure for deleting node ID-th of the table.
  subroutine del(table,ID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Hash_Table), intent(INOUT):: table !< Hash table.
  integer(I8P),           intent(IN)::    ID    !< ID (unique) of the node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call table%ht(table%hash(ID=ID))%del(ID=ID)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine del

  !> @brief Procedure for computing the length of the table.
  function length(table) result(Nn)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Hash_Table), intent(IN):: table !< Hash table.
  integer(I4P)::                       Nn    !< Nodes number.
  integer(I4P)::                       b     !< Bucket counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nn = 0
  if (allocated(table%ht)) THEN
    do b=lbound(table%ht,dim=1),ubound(table%ht,dim=1)
       Nn = Nn + table%ht(b)%length()
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction length

  !> @brief Procedure for getting the list of actually stored IDs.
  subroutine getIDs(table)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Hash_Table), intent(INOUT):: table  !< Hash table.
  integer(I4P)::                          length !< Actual table length.
  integer(I4P)::                          b,i    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(table%ht)) THEN
    length=table%length()
    if (length>0) then
      if (allocated(table%IDs)) deallocate(table%IDs) ; allocate(table%IDs(1:length))
      i=1
      do b=lbound(table%ht,dim=1),ubound(table%ht,dim=1)
        length=table%ht(b)%length()
        call table%ht(b)%getIDs(IDs=table%IDs(i:i+length-1))
        i=i+length
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine getIDs

  !> @brief Procedure for printing IDs list with a pretty format.
  subroutine print(table,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Hash_Table), intent(IN)::  table   !< Hash table.
  character(*), optional, intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::       prefd   !< Prefixing string.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  integer(I4P)::                        b       !< Buckets counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  if (allocated(table%ht)) THEN
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' The hash table contains the following IDs'
    do b=lbound(table%ht,dim=1),ubound(table%ht,dim=1)
      call table%ht(b)%print(pref=prefd//"  ",iostat=iostatd,iomsg=iomsgd,unit=unit)
    enddo
  endif
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print
  !> @}
endmodule Data_Type_Hash_Table
