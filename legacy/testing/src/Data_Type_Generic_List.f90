!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Generic_ListDerivedType Data_Type_Generic_List
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_Generic_ListPublicProcedure Data_Type_Generic_List
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Generic_ListPrivateProcedure Data_Type_Generic_List
!> @}

!> @brief Module Data_Type_Generic_List contains the definition of Type_Generic_List type and useful procedures for its handling.
module Data_Type_Generic_List
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision !< Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the node of generic list thus defining the list itself.
!> @ingroup Data_Type_Generic_ListDerivedType
type, public:: Type_Generic_List
  private
  integer(I8P),            allocatable:: ID             !< ID (unique) of the current node.
  class(*),                pointer::     d    => null() !< Node data.
  type(Type_Generic_List), pointer::     next => null() !< Pointer to the next node of the list.
  contains
    procedure:: free     ! Procedure for freeing (destroying) the list.
    procedure:: node     ! Procedure for returning the ID-th node pointer of the list.
    procedure:: dat      ! Procedure for returning the data pointer of the ID-th node of the list.
    procedure:: put      ! Procedure for inserting data into the ID-th node of the list.
    procedure:: get      ! Procedure for getting data from the ID-th node of the list.
    procedure:: del      ! Procedure for deleting the ID-th node of the list.
    procedure:: length   ! Procedure for computing the length of the list.
    procedure:: getIDs   ! Procedure for getting the list of actually stored IDs.
    procedure:: print    ! Procedure for printing IDs list with a pretty format.
    final::     finalize ! Procedure for freeing dynamic memory when finalizing.
endtype Type_Generic_List
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Generic_ListPrivateProcedure
  !> @{
  !> @brief Procedure for freeing (destroying) the list.
  recursive subroutine free(list)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Generic_List), intent(INOUT):: list !< List.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(list%next)) then
    call free(list%next)
    deallocate(list%next)
  endif
  if (allocated( list%ID)) deallocate(list%ID)
  if (associated(list%d )) deallocate(list%d )
  list%d    => null()
  list%next => null()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Procedure for freeing dynamic memory when finalizing.
  subroutine finalize(list)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Generic_List), intent(INOUT):: list !< List.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call list%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for returning the ID-th node pointer of the list.
  !> @note If ID key is not present a null pointer is returned.
  function node(list,ID) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Generic_List), target, intent(IN):: list !< List.
  integer(I8P),                     intent(IN):: ID   !< Unique key of the node of the list to be found.
  type(Type_Generic_List), pointer::             n    !< Pointer to "ID-th" node of the list.
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
  endfunction node

  !> @brief Procedure for returning the data pointer of the ID-th node of the list.
  !> @note If ID key is not present a null pointer is returned.
  function dat(list,ID) result(d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Generic_List), target, intent(IN):: list  !< List.
  integer(I8P),                     intent(IN):: ID    !< Unique key of the node of the list to be found.
  class(*), pointer::                            d     !< Pointer to the data of the "ID-th" node of the list.
  type(Type_Generic_List), pointer::             n     !< Pointer to "ID-th" node of the list.
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
  endfunction dat

  !> @brief Procedure for inserting data into the ID-th node of the list. If a node with the provided ID is not present into the
  !> list a new one is created.
  recursive subroutine put(list,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Generic_List), intent(INOUT):: list !< List.
  integer(I8P),             intent(IN)::    ID   !< ID (unique) of the current node.
  class(*),                 intent(IN)::    d    !< Data of the ID-th node.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(list%ID)) then
    if (list%ID/=ID) then
      if (.not.associated(list%next)) allocate(list%next)
      call put(list=list%next,ID=ID,d=d)
    else
      if (associated(list%d)) deallocate(list%d) ; allocate(list%d,source=d)
    endif
  else
    allocate(list%ID,source=ID)
    if (associated(list%d)) deallocate(list%d) ; allocate(list%d,source=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put

  !> @brief Procedure for getting data from the ID-th node of the list. If the requested node ID is not present, the output data
  !> is returned in deallocated form.
  recursive subroutine get(list,ID,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Generic_List), intent(IN)::  list !< Llist.
  integer(I8P),             intent(IN)::  ID   !< ID (unique) of the current node.
  class(*), allocatable,    intent(OUT):: d    !< Data of the ID-th node.
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
    call get(list=list%next,ID=ID,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get

  !> @brief Procedure for deleting the ID-th node of the list.
  subroutine del(list,ID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Generic_List), target, intent(INOUT):: list !< List.
  integer(I8P),                     intent(IN)::    ID   !< ID (unique) of the node.
  type(Type_Generic_List), pointer::                c,n  !< Pointers for scanning the list.
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
  endsubroutine del

  !> @brief Procedure for computing the length of the list.
  function length(list) result(Nl)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Generic_List), intent(IN):: list !< List.
  integer(I4P)::                         Nl   !< Elements number.
  type(Type_Generic_List), pointer::     n    !< Pointer for scanning the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = 0 ; if (associated(list%d)) Nl = 1
  n => list%next
  do while (associated(n))
    if (associated(n%d)) then
      Nl = Nl + 1
    endif
    n => n%next
  enddo
  n => null()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction length

  !> @brief Procedure for getting the list of actually stored IDs.
  subroutine getIDs(list,IDs)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Generic_List), intent(IN)::    list    !< List.
  integer(I8P),             intent(INOUT):: IDs(1:) !< List of actually stored IDs.
  type(Type_Generic_List), pointer::        n       !< Pointer for list scanning.
  integer(I4P)::                            i       !< Counter.
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
  endsubroutine getIDs

  !> @brief Procedure for printing IDs list with a pretty format.
  subroutine print(list,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Generic_List), target, intent(IN)::  list    !< List.
  character(*), optional,           intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional,           intent(OUT):: iostat  !< IO error.
  character(*), optional,           intent(OUT):: iomsg   !< IO error message.
  integer(I4P),                     intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::                 prefd   !< Prefixing string.
  integer(I4P)::                                  iostatd !< IO error.
  character(500)::                                iomsgd  !< Temporary variable for IO error message.
  type(Type_Generic_List), pointer::              n       !< Pointer for scanning the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  n => list
  do while(allocated(n%ID))
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' ID='//str(n=n%ID)
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
  endsubroutine print
  !> @}
endmodule Data_Type_Generic_List
