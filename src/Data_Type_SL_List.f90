!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_SL_ListDerivedType Data_Type_SL_List
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_SL_ListPublicProcedure Data_Type_SL_List
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_SL_ListPrivateProcedure Data_Type_SL_List
!> @}

!> @brief Module Data_Type_SL_List contains the definition of Type_SL_List type and useful procedures for its handling.
module Data_Type_SL_List
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision !< Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: SL_List_Mold
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Parametrized types ID for identifying list data.
integer(I1P), parameter:: tp_NULL = 0  !< Undefined type.
integer(I1P), parameter:: tp_SI1  = 1  !< Scalar integer(I1P) type.
integer(I1P), parameter:: tp_SI2  = 2  !< Scalar integer(I2P) type.
integer(I1P), parameter:: tp_SI4  = 3  !< Scalar integer(I4P) type.
integer(I1P), parameter:: tp_SI8  = 4  !< Scalar integer(I8P) type.
integer(I1P), parameter:: tp_SR4  = 5  !< Scalar real(R4P) type.
integer(I1P), parameter:: tp_SR8  = 6  !< Scalar real(R8P) type.
integer(I1P), parameter:: tp_SR16 = 7  !< Scalar real(R16P) type.
integer(I1P), parameter:: tp_SCH  = 8  !< Scalar character type.
integer(I1P), parameter:: tp_SLOG = 9  !< Scalar character type.
integer(I1P), parameter:: tp_AI1  = 10 !< Array of integer(I1P) type.
integer(I1P), parameter:: tp_AI2  = 11 !< Array of integer(I2P) type.
integer(I1P), parameter:: tp_AI4  = 12 !< Array of integer(I4P) type.
integer(I1P), parameter:: tp_AI8  = 13 !< Array of integer(I8P) type.
integer(I1P), parameter:: tp_AR4  = 14 !< Array of real(R4P) type.
integer(I1P), parameter:: tp_AR8  = 15 !< Array of real(R8P) type.
integer(I1P), parameter:: tp_AR16 = 16 !< Array of real(R16P) type.
integer(I1P), parameter:: tp_ACH  = 17 !< Array of character type.
integer(I1P), parameter:: tp_ALOG = 18 !< Array of character type.
!> @brief Derived type containing link data of a single linked list.
!> @ingroup Data_Type_SL_ListDerivedType
type:: Type_SL_Link
  private
  integer(I1p)::                tp=tp_NULL   !< Type of data.
  integer(I1P),       pointer:: d(:)         !< Link data.
  type(Type_SL_Link), pointer:: n => null()  !< Pointer to the next link of the list.
  contains
    procedure, non_overridable:: free => free_link        ! Procedure for freeing (destroying) the list.
    procedure, non_overridable:: leng => leng_link        ! Procedure for computing the length of the list.
    procedure, non_overridable:: homogeneous => homo_link ! Procedure for checking if the list is homogeneous.
    procedure, non_overridable:: link                     ! Procedure for returning the n-th link pointer of the list.
    procedure, non_overridable:: del => del_link          ! Procedure for deleting the n-th link of the list.
    procedure, non_overridable:: delh => delh_link        ! Procedure for deleting head link of the list.
    procedure, non_overridable:: delt => delt_link        ! Procedure for deleting tail link of the list.
    generic::                    get =>  &         ! Procedure for getting data from the n-th link of the list (generic interface).
                                 get_link_R8,get_link_R4, &
                                 get_link_I8,get_link_I4,get_link_I2,get_link_I1,get_link_ch,get_link_log
    generic::                    geth =>  &        ! Procedure for getting data from the list head (generic interface).
#ifdef r16p
                                 geth_link_R16,&
#endif
                                 geth_link_R8,geth_link_R4, &
                                 geth_link_I8,geth_link_I4,geth_link_I2,geth_link_I1,geth_link_ch,geth_link_log
    generic::                    gett =>  &        ! Procedure for getting a link data from the list tail (generic interface).
#ifdef r16p
                                 gett_link_R16,&
#endif
                                 gett_link_R8,gett_link_R4, &
                                 gett_link_I8,gett_link_I4,gett_link_I2,gett_link_I1,gett_link_ch,gett_link_log
    generic::                    put =>  &         ! Procedure for inserting a link into the n-th link of list (generic interface).
#ifdef r16p
                                 put_link_R16,&
#endif
                                 put_link_R8,put_link_R4, &
                                 put_link_I8,put_link_I4,put_link_I2,put_link_I1,put_link_ch,put_link_log
    generic::                    puth =>  &        ! Procedure for inserting a link into the list head (generic interface).
#ifdef r16p
                                 puth_link_R16,&
#endif
                                 puth_link_R8,puth_link_R4, &
                                 puth_link_I8,puth_link_I4,puth_link_I2,puth_link_I1,puth_link_ch,puth_link_log
    generic::                    putt =>  &        ! Procedure for inserting a link into the list tail (generic interface).
#ifdef r16p
                                 putt_link_R16,&
#endif
                                 putt_link_R8,putt_link_R4, &
                                 putt_link_I8,putt_link_I4,putt_link_I2,putt_link_I1,putt_link_ch,putt_link_log
    generic::                    array =>  &       ! Procedure for converting list to array (generic interface).
#ifdef r16p
                                 array_link_R16,&
#endif
                                 array_link_R8,array_link_R4, &
                                 array_link_I8,array_link_I4,array_link_I2,array_link_I1,array_link_ch,array_link_log
#ifdef r16p
    procedure, non_overridable:: get_link_R16      ! Procedure for getting data from the n-th link of the list (R16P).
#endif
    procedure, non_overridable:: get_link_R8       ! Procedure for getting data from the n-th link of the list (R8P).
    procedure, non_overridable:: get_link_R4       ! Procedure for getting data from the n-th link of the list (R4P).
    procedure, non_overridable:: get_link_I8       ! Procedure for getting data from the n-th link of the list (I8P).
    procedure, non_overridable:: get_link_I4       ! Procedure for getting data from the n-th link of the list (I4P).
    procedure, non_overridable:: get_link_I2       ! Procedure for getting data from the n-th link of the list (I2P).
    procedure, non_overridable:: get_link_I1       ! Procedure for getting data from the n-th link of the list (I1P).
    procedure, non_overridable:: get_link_ch       ! Procedure for getting data from the n-th link of the list (character).
    procedure, non_overridable:: get_link_log      ! Procedure for getting data from the n-th link of the list (logical).
#ifdef r16p
    procedure, non_overridable:: geth_link_R16     ! Procedure for getting a link data from the list head (R16P).
#endif
    procedure, non_overridable:: geth_link_R8      ! Procedure for getting a link data from the list head (R8P).
    procedure, non_overridable:: geth_link_R4      ! Procedure for getting a link data from the list head (R4P).
    procedure, non_overridable:: geth_link_I8      ! Procedure for getting a link data from the list head (I4P).
    procedure, non_overridable:: geth_link_I4      ! Procedure for getting a link data from the list head (I4P).
    procedure, non_overridable:: geth_link_I2      ! Procedure for getting a link data from the list head (I4P).
    procedure, non_overridable:: geth_link_I1      ! Procedure for getting a link data from the list head (I4P).
    procedure, non_overridable:: geth_link_ch      ! Procedure for getting a link data from the list head (character).
    procedure, non_overridable:: geth_link_log     ! Procedure for getting a link data from the list head (logical).
#ifdef r16p
    procedure, non_overridable:: gett_link_R16     ! Procedure for getting a link data from the list tail (R16P).
#endif
    procedure, non_overridable:: gett_link_R8      ! Procedure for getting a link data from the list tail (R8P).
    procedure, non_overridable:: gett_link_R4      ! Procedure for getting a link data from the list tail (R4P).
    procedure, non_overridable:: gett_link_I8      ! Procedure for getting a link data from the list tail (I4P).
    procedure, non_overridable:: gett_link_I4      ! Procedure for getting a link data from the list tail (I4P).
    procedure, non_overridable:: gett_link_I2      ! Procedure for getting a link data from the list tail (I4P).
    procedure, non_overridable:: gett_link_I1      ! Procedure for getting a link data from the list tail (I4P).
    procedure, non_overridable:: gett_link_ch      ! Procedure for getting a link data from the list tail (character).
    procedure, non_overridable:: gett_link_log     ! Procedure for getting a link data from the list tail (logical).
#ifdef r16p
    procedure, non_overridable:: put_link_R16      ! Procedure for inserting a link into the n-th link of list (R16P).
#endif
    procedure, non_overridable:: put_link_R8       ! Procedure for inserting a link into the n-th link of list (R8P).
    procedure, non_overridable:: put_link_R4       ! Procedure for inserting a link into the n-th link of list (R4P).
    procedure, non_overridable:: put_link_I8       ! Procedure for inserting a link into the n-th link of list (I8P).
    procedure, non_overridable:: put_link_I4       ! Procedure for inserting a link into the n-th link of list (I4P).
    procedure, non_overridable:: put_link_I2       ! Procedure for inserting a link into the n-th link of list (I2P).
    procedure, non_overridable:: put_link_I1       ! Procedure for inserting a link into the n-th link of list (I1P).
    procedure, non_overridable:: put_link_ch       ! Procedure for inserting a link into the n-th link of list (character).
    procedure, non_overridable:: put_link_log      ! Procedure for inserting a link into the n-th link of list (logical).
#ifdef r16p
    procedure, non_overridable:: puth_link_R16     ! Procedure for inserting a link into the list head (R16P).
#endif
    procedure, non_overridable:: puth_link_R8      ! Procedure for inserting a link into the list head (R8P).
    procedure, non_overridable:: puth_link_R4      ! Procedure for inserting a link into the list head (R4P).
    procedure, non_overridable:: puth_link_I8      ! Procedure for inserting a link into the list head (I8P).
    procedure, non_overridable:: puth_link_I4      ! Procedure for inserting a link into the list head (I4P).
    procedure, non_overridable:: puth_link_I2      ! Procedure for inserting a link into the list head (I2P).
    procedure, non_overridable:: puth_link_I1      ! Procedure for inserting a link into the list head (I1P).
    procedure, non_overridable:: puth_link_ch      ! Procedure for inserting a link into the list head (character).
    procedure, non_overridable:: puth_link_log     ! Procedure for inserting a link into the list head (logical).
#ifdef r16p
    procedure, non_overridable:: putt_link_R16     ! Procedure for inserting a link into the list tail (R16P).
#endif
    procedure, non_overridable:: putt_link_R8      ! Procedure for inserting a link into the list tail (R8P).
    procedure, non_overridable:: putt_link_R4      ! Procedure for inserting a link into the list tail (R4P).
    procedure, non_overridable:: putt_link_I8      ! Procedure for inserting a link into the list tail (I8P).
    procedure, non_overridable:: putt_link_I4      ! Procedure for inserting a link into the list tail (I4P).
    procedure, non_overridable:: putt_link_I2      ! Procedure for inserting a link into the list tail (I2P).
    procedure, non_overridable:: putt_link_I1      ! Procedure for inserting a link into the list tail (I1P).
    procedure, non_overridable:: putt_link_ch      ! Procedure for inserting a link into the list tail (character).
    procedure, non_overridable:: putt_link_log     ! Procedure for inserting a link into the list tail (logical).
#ifdef r16p
    procedure, non_overridable:: array_link_R16    ! Procedure for converting list to array (R16P).
#endif
    procedure, non_overridable:: array_link_R8     ! Procedure for converting list to array (R8P).
    procedure, non_overridable:: array_link_R4     ! Procedure for converting list to array (R4P).
    procedure, non_overridable:: array_link_I8     ! Procedure for converting list to array (I8P).
    procedure, non_overridable:: array_link_I4     ! Procedure for converting list to array (I4P).
    procedure, non_overridable:: array_link_I2     ! Procedure for converting list to array (I2P).
    procedure, non_overridable:: array_link_I1     ! Procedure for converting list to array (I1P).
    procedure, non_overridable:: array_link_ch     ! Procedure for converting list to array (character).
    procedure, non_overridable:: array_link_log    ! Procedure for converting list to array (logical).
endtype Type_SL_Link
!> @brief Derived type containing single linked list.
!> Type_SL_List is a \b generic Single Linked List. The term \b generic means that the data stored in each node (link) of the list
!> is \b generic: at present, the stored data can be integer (of any kinds defined in IR_Precision module), real (of any kinds
!> defined in IR_Precision module), characters (of any length) and logical. In order to insert and retrieve data from the list nodes
!> (and in general for list manipulation) the type bound methods must be used, because the data is not directly accessible. This is
!> due to the internal representation of the data. In order to allow generic data, the intrinsic \b transfer function is used to
!> encode all user data into the internal integer(I1P) array representation. As a consequence, if the node is directed accessed the
!> user can obtain unpredictable results due to the encoded internal representation.
!> The provide methods for list handling are the following:
!> @code
!> ...
!> type(Type_SL_List):: list
!> ...
!> @endcode
!> - insert data into list head:
!>   @code
!>   call list%puth(d=#user_data)
!>   @endcode
!> - insert data into list tail:
!>   @code
!>   call list%putt(d=#user_data)
!>   @endcode
!> - retrieve data from list at some node "n":
!>   @code
!>   call list%get(n=n,d=#user_data)
!>   @endcode
!> - retrieve data from list head:
!>   @code
!>   call list%geth(d=#user_data)
!>   @endcode
!> - retrieve data from list tail:
!>   @code
!>   call list%gett(d=#user_data)
!>   @endcode
!> - deleting list head link:
!>   @code
!>   call list%delh
!>   @endcode
!> - deleting list tail link:
!>   @code
!>   call list%delt
!>   @endcode
!> - freeing (destroying) list:
!>   @code
!>   call list%free
!>   @endcode
!> - converting list to array:
!>   @code
!>   call list%array(a=#user_array)
!>   @endcode
!> @note For the character data, the methods get and array need a supplementary arguments, namely \b Nc being the length of the
!> character variable (array element for array method).
!> @note Array method converts the list into an array thus the list must have \b homogeneous data. Even if the list can have nodes
!> with different data type the conversion to a standard array can be completed successfully only if all data nodes are of the same
!> kind (e.g. all real or all integer or all character etc...).
!> @ingroup Data_Type_SL_ListDerivedType
type, public:: Type_SL_List
  private
  integer(I4P), public::        l = 0         !< List length.
  logical, public::             homo = .true. !< Flag for checking if the list is homogeneous.
  type(Type_SL_Link), pointer:: f => null()   !< Pointer to the first link of the list.
  contains
    procedure, non_overridable:: free => free_list        ! Procedure for freeing (destroying) the list.
    procedure, non_overridable:: leng => leng_list        ! Procedure for computing the length of the list.
    procedure, non_overridable:: homogeneous => homo_list ! Procedure for checking if the list is homogeneous.
    procedure, non_overridable:: del => del_list          ! Procedure for deleting the n-th link of the list.
    procedure, non_overridable:: delh => delh_list        ! Procedure for deleting head link of the list.
    procedure, non_overridable:: delt => delt_list        ! Procedure for deleting tail link of the list.
    generic::                    get =>  &         ! Procedure for getting data from the n-th link of the list (generic interface).
#ifdef r16p
                                 get_list_R16,&
#endif
                                 get_list_R8,get_list_R4, &
                                 get_list_I8,get_list_I4,get_list_I2,get_list_I1,get_list_ch,get_list_log
    generic::                    geth =>  &        ! Procedure for getting a link data from the list head (generic interface).
#ifdef r16p
                                 geth_list_R16,&
#endif
                                 geth_list_R8,geth_list_R4, &
                                 geth_list_I8,geth_list_I4,geth_list_I2,geth_list_I1,geth_list_ch,geth_list_log
    generic::                    gett =>  &        ! Procedure for getting a link data from the list tail (generic interface).
#ifdef r16p
                                 gett_list_R16,&
#endif
                                 gett_list_R8,gett_list_R4, &
                                 gett_list_I8,gett_list_I4,gett_list_I2,gett_list_I1,gett_list_ch,gett_list_log
    generic::                    put =>  &         ! Procedure for inserting a link into the n-th link of list (generic interface).
#ifdef r16p
                                 put_list_R16,&
#endif
                                 put_list_R8,put_list_R4, &
                                 put_list_I8,put_list_I4,put_list_I2,put_list_I1,put_list_ch,put_list_log
    generic::                    puth =>  &        ! Procedure for inserting a link into the list head (generic interface).
#ifdef r16p
                                 puth_list_R16,&
#endif
                                 puth_list_R8,puth_list_R4, &
                                 puth_list_I8,puth_list_I4,puth_list_I2,puth_list_I1,puth_list_ch,puth_list_log
    generic::                    putt =>  &        ! Procedure for inserting a link into the list tail (generic interface).
#ifdef r16p
                                 putt_list_R16,&
#endif
                                 putt_list_R8,putt_list_R4, &
                                 putt_list_I8,putt_list_I4,putt_list_I2,putt_list_I1,putt_list_ch,putt_list_log
    generic::                    array =>  &       ! Procedure for converting list to array (generic interface).
#ifdef r16p
                                 array_list_R16,&
#endif
                                 array_list_R8,array_list_R4, &
                                 array_list_I8,array_list_I4,array_list_I2,array_list_I1,array_list_ch,array_list_log
#ifdef r16p
    procedure, non_overridable:: get_list_R16      ! Procedure for getting data from the n-th link of the list (R16P).
#endif
    procedure, non_overridable:: get_list_R8       ! Procedure for getting data from the n-th link of the list (R8P).
    procedure, non_overridable:: get_list_R4       ! Procedure for getting data from the n-th link of the list (R4P).
    procedure, non_overridable:: get_list_I8       ! Procedure for getting data from the n-th link of the list (I8P).
    procedure, non_overridable:: get_list_I4       ! Procedure for getting data from the n-th link of the list (I4P).
    procedure, non_overridable:: get_list_I2       ! Procedure for getting data from the n-th link of the list (I2P).
    procedure, non_overridable:: get_list_I1       ! Procedure for getting data from the n-th link of the list (I1P).
    procedure, non_overridable:: get_list_ch       ! Procedure for getting data from the n-th link of the list (character).
    procedure, non_overridable:: get_list_log      ! Procedure for getting data from the n-th link of the list (logical).
#ifdef r16p
    procedure, non_overridable:: geth_list_R16     ! Procedure for getting a link data from the list head (R16P).
#endif
    procedure, non_overridable:: geth_list_R8      ! Procedure for getting a link data from the list head (R8P).
    procedure, non_overridable:: geth_list_R4      ! Procedure for getting a link data from the list head (R4P).
    procedure, non_overridable:: geth_list_I8      ! Procedure for getting a link data from the list head (I4P).
    procedure, non_overridable:: geth_list_I4      ! Procedure for getting a link data from the list head (I4P).
    procedure, non_overridable:: geth_list_I2      ! Procedure for getting a link data from the list head (I4P).
    procedure, non_overridable:: geth_list_I1      ! Procedure for getting a link data from the list head (I4P).
    procedure, non_overridable:: geth_list_ch      ! Procedure for getting a link data from the list head (character).
    procedure, non_overridable:: geth_list_log     ! Procedure for getting a link data from the list head (logical).
#ifdef r16p
    procedure, non_overridable:: gett_list_R16     ! Procedure for getting a link data from the list tail (R16P).
#endif
    procedure, non_overridable:: gett_list_R8      ! Procedure for getting a link data from the list tail (R8P).
    procedure, non_overridable:: gett_list_R4      ! Procedure for getting a link data from the list tail (R4P).
    procedure, non_overridable:: gett_list_I8      ! Procedure for getting a link data from the list tail (I4P).
    procedure, non_overridable:: gett_list_I4      ! Procedure for getting a link data from the list tail (I4P).
    procedure, non_overridable:: gett_list_I2      ! Procedure for getting a link data from the list tail (I4P).
    procedure, non_overridable:: gett_list_I1      ! Procedure for getting a link data from the list tail (I4P).
    procedure, non_overridable:: gett_list_ch      ! Procedure for getting a link data from the list tail (character).
    procedure, non_overridable:: gett_list_log     ! Procedure for getting a link data from the list tail (logical).
#ifdef r16p
    procedure, non_overridable:: put_list_R16      ! Procedure for inserting a link into the n-th link of list (R16P).
#endif
    procedure, non_overridable:: put_list_R8       ! Procedure for inserting a link into the n-th link of list (R8P).
    procedure, non_overridable:: put_list_R4       ! Procedure for inserting a link into the n-th link of list (R4P).
    procedure, non_overridable:: put_list_I8       ! Procedure for inserting a link into the n-th link of list (I8P).
    procedure, non_overridable:: put_list_I4       ! Procedure for inserting a link into the n-th link of list (I4P).
    procedure, non_overridable:: put_list_I2       ! Procedure for inserting a link into the n-th link of list (I2P).
    procedure, non_overridable:: put_list_I1       ! Procedure for inserting a link into the n-th link of list (I1P).
    procedure, non_overridable:: put_list_ch       ! Procedure for inserting a link into the n-th link of list (character).
    procedure, non_overridable:: put_list_log      ! Procedure for inserting a link into the n-th link of list (logical).
#ifdef r16p
    procedure, non_overridable:: puth_list_R16     ! Procedure for inserting a link into the list head (R16P).
#endif
    procedure, non_overridable:: puth_list_R8      ! Procedure for inserting a link into the list head (R8P).
    procedure, non_overridable:: puth_list_R4      ! Procedure for inserting a link into the list head (R4P).
    procedure, non_overridable:: puth_list_I8      ! Procedure for inserting a link into the list head (I8P).
    procedure, non_overridable:: puth_list_I4      ! Procedure for inserting a link into the list head (I4P).
    procedure, non_overridable:: puth_list_I2      ! Procedure for inserting a link into the list head (I2P).
    procedure, non_overridable:: puth_list_I1      ! Procedure for inserting a link into the list head (I1P).
    procedure, non_overridable:: puth_list_ch      ! Procedure for inserting a link into the list head (character).
    procedure, non_overridable:: puth_list_log     ! Procedure for inserting a link into the list head (logical).
#ifdef r16p
    procedure, non_overridable:: putt_list_R16     ! Procedure for inserting a link into the list tail (R16P).
#endif
    procedure, non_overridable:: putt_list_R8      ! Procedure for inserting a link into the list tail (R8P).
    procedure, non_overridable:: putt_list_R4      ! Procedure for inserting a link into the list tail (R4P).
    procedure, non_overridable:: putt_list_I8      ! Procedure for inserting a link into the list tail (I8P).
    procedure, non_overridable:: putt_list_I4      ! Procedure for inserting a link into the list tail (I4P).
    procedure, non_overridable:: putt_list_I2      ! Procedure for inserting a link into the list tail (I2P).
    procedure, non_overridable:: putt_list_I1      ! Procedure for inserting a link into the list tail (I1P).
    procedure, non_overridable:: putt_list_ch      ! Procedure for inserting a link into the list tail (character).
    procedure, non_overridable:: putt_list_log     ! Procedure for inserting a link into the list tail (logical).
#ifdef r16p
    procedure, non_overridable:: array_list_R16    ! Procedure for converting list to array (R16P).
#endif
    procedure, non_overridable:: array_list_R8     ! Procedure for converting list to array (R8P).
    procedure, non_overridable:: array_list_R4     ! Procedure for converting list to array (R4P).
    procedure, non_overridable:: array_list_I8     ! Procedure for converting list to array (I8P).
    procedure, non_overridable:: array_list_I4     ! Procedure for converting list to array (I4P).
    procedure, non_overridable:: array_list_I2     ! Procedure for converting list to array (I2P).
    procedure, non_overridable:: array_list_I1     ! Procedure for converting list to array (I1P).
    procedure, non_overridable:: array_list_ch     ! Procedure for converting list to array (character).
    procedure, non_overridable:: array_list_log    ! Procedure for converting list to array (logical).
endtype Type_SL_List
!> @ingroup Data_Type_SL_ListGlobalVarPar
!> @{
integer(I1P), pointer:: SL_List_Mold(:)=>null() !< Mold prototype for list data.
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_SL_ListPrivateProcedure
  !> @{
  !> @brief Recursive subroutine for freeing (destroying) the list links.
  recursive subroutine free_link(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%n)) then
    call free_link(l%n)
    deallocate(l%n)
  endif
  l%n => null()
  if (associated(l%d)) deallocate(l%d)
  l%d => null()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_link

  !> @brief Function for computing the length of the list links.
  function leng_link(l) result(Nl)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN):: l    !< List.
  integer(I4P)::                    Nl   !< Elements number.
  type(Type_SL_Link), pointer::     n    !< Pointer for scanning the list.
  integer(I1P),       pointer::     d(:) !< Pointer for scanning the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = 0
  d => l%d
  n => l%n
  do while (associated(d))
    Nl = Nl + 1
    if (associated(n)) then
      d => n%d
      n => n%n
    else
      d => null()
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction leng_link

  !> @brief Function for checking if the list is homogeneous.
  function homo_link(l) result(homo)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN):: l    !< List.
  type(Type_SL_Link), pointer::     n    !< Pointer for scanning the list.
  integer(I1P)::                    tp   !< Type counter.
  logical::                         homo !< Flag for checking if the list is homogeneous.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  homo = .true.
  tp = l%tp
  n => l%n
  do while (associated(n))
    homo = (homo.and.tp==n%tp)
    n => n%n
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction homo_link

  !> @brief Function for returning the n-th link pointer of the list.
  !> @note If n<1 then the first link of the list is returned, whereas if n> list's length the last link is returned.
  function link(l,n) result(ln)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(IN):: l  !< List.
  integer(I4P),                intent(IN):: n  !< Element of the list to be found.
  type(Type_SL_Link), pointer::             ln !< Pointer to "n-th" link of the list.
  integer(I4P)::                            e  !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ln => l
  if (n>1) then
    do e=2,n
      if (associated(ln%n)) then
        ln => ln%n
      else
        exit
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction link

  !> @brief Subroutine for deleting n-th link of the list.
  subroutine del_link(l,n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l  !< List.
  integer(I4P),                intent(IN)::    n  !< Element of the list to be deleted.
  type(Type_SL_Link), pointer::                ln !< Pointer to "n-th" link of the list.
  type(Type_SL_Link), pointer::                tmp !< Temporary link of the list.
  integer(I4P)::                               e  !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (n<=1) then
    call l%delh
  else
    ln => l
    do e=2,n
      if (associated(ln%n)) then
        ln => ln%n
      else
        exit
      endif
    enddo
    if (associated(ln%d)) then
      if (associated(ln%n)) then
        tmp => ln%n
        ln%tp = tmp%tp
        deallocate(ln%d) ; allocate(ln%d(1:size(tmp%d))) ; ln%d = tmp%d
        ln%n => tmp%n
        deallocate(tmp)
      else
        ln%tp = tp_NULL
        deallocate(ln%d)
        ln%d => null()
      endif
    endif
  endif

  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine del_link

  !> @brief Subroutine for deleting head link of the list.
  subroutine delh_link(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  type(Type_SL_Link), pointer::        n !< Pointer to next link of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%d)) then
    if (associated(l%n)) then
      n => l%n
      l%tp = n%tp
      deallocate(l%d) ; allocate(l%d(1:size(n%d))) ; l%d = n%d
      l%n => n%n
      deallocate(n)
    else
      l%tp = tp_NULL
      deallocate(l%d)
      l%d => null()
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine delh_link

  !> @brief Subroutine for deleting tail link of the list.
  subroutine delt_link(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l     !< List.
  type(Type_SL_Link), pointer::                p,c,n !< Pointers for scanning the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  p => l
  c => l%n
  if (associated(c%n)) then
    n => c%n
  else
    n => c
  endif
  do while (associated(n))
    p => c
    c => n
    n => n%n
  enddo
  if (associated(p%n)) then
    deallocate(p%n)
    p%n => null()
  else
    deallocate(l%n)
    l%n => null()
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine delt_link

#ifdef r16p
  !> @brief Subroutine for getting data from the n-th link of list (R16P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_link_R16(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(IN)::  l            !< List.
  integer(I4P),                intent(IN)::  n            !< Element of the list to be found.
  real(R16P),                  intent(OUT):: d            !< Link data.
  type(Type_SL_Link), pointer::              ln => null() !< Pointer to "n-th" link of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR16P
  ln => link(l=l,n=n)
  if (associated(ln%d)) d = transfer(ln%d,1._R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_link_R16
#endif

  !> @brief Subroutine for getting data from the n-th link of list (R8P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_link_R8(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(IN)::  l            !< List.
  integer(I4P),                intent(IN)::  n            !< Element of the list to be found.
  real(R8P),                   intent(OUT):: d            !< Link data.
  type(Type_SL_Link), pointer::              ln => null() !< Pointer to "n-th" link of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR8P
  ln => link(l=l,n=n)
  if (associated(ln%d)) d = transfer(ln%d,1._R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_link_R8

  !> @brief Subroutine for getting data from the n-th link of list (R4P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_link_R4(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(IN)::  l            !< List.
  integer(I4P),                intent(IN)::  n            !< Element of the list to be found.
  real(R4P),                   intent(OUT):: d            !< Link data.
  type(Type_SL_Link), pointer::              ln => null() !< Pointer to "n-th" link of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR4P
  ln => link(l=l,n=n)
  if (associated(ln%d)) d = transfer(ln%d,1._R4P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_link_R4

  !> @brief Subroutine for getting data from the n-th link of list (I8P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_link_I8(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(IN)::  l            !< List.
  integer(I4P),                intent(IN)::  n            !< Element of the list to be found.
  integer(I8P),                intent(OUT):: d            !< Link data.
  type(Type_SL_Link), pointer::              ln => null() !< Pointer to "n-th" link of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI8P
  ln => link(l=l,n=n)
  if (associated(ln%d)) d = transfer(ln%d,1_I8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_link_I8

  !> @brief Subroutine for getting data from the n-th link of list (I4P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_link_I4(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(IN)::  l            !< List.
  integer(I4P),                intent(IN)::  n            !< Element of the list to be found.
  integer(I4P),                intent(OUT):: d            !< Link data.
  type(Type_SL_Link), pointer::              ln => null() !< Pointer to "n-th" link of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI4P
  ln => link(l=l,n=n)
  if (associated(ln%d)) d = transfer(ln%d,1_I4P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_link_I4

  !> @brief Subroutine for getting data from the n-th link of list (I2P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_link_I2(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(IN)::  l            !< List.
  integer(I4P),                intent(IN)::  n            !< Element of the list to be found.
  integer(I2P),                intent(OUT):: d            !< Link data.
  type(Type_SL_Link), pointer::              ln => null() !< Pointer to "n-th" link of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI2P
  ln => link(l=l,n=n)
  if (associated(ln%d)) d = transfer(ln%d,1_I2P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_link_I2

  !> @brief Subroutine for getting data from the n-th link of list (I1P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_link_I1(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(IN)::  l            !< List.
  integer(I4P),                intent(IN)::  n            !< Element of the list to be found.
  integer(I1P),                intent(OUT):: d            !< Link data.
  type(Type_SL_Link), pointer::              ln => null() !< Pointer to "n-th" link of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI1P
  ln => link(l=l,n=n)
  if (associated(ln%d)) d = transfer(ln%d,1_I1P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_link_I1

  !> @brief Subroutine for getting data from the n-th link of list (character).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_link_ch(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(IN)::  l            !< List.
  integer(I4P),                intent(IN)::  n            !< Element of the list to be found.
  character(*),                intent(OUT):: d            !< Link data.
  character(len(d))::                        mold         !< Mold prototype.
  type(Type_SL_Link), pointer::              ln => null() !< Pointer to "n-th" link of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = ''
  ln => link(l=l,n=n)
  if (associated(ln%d)) d = transfer(ln%d,mold(1:min(len(d),size(ln%d,dim=1))))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_link_ch

  !> @brief Subroutine for getting data from the n-th link of list (logical).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_link_log(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(IN)::  l            !< List.
  integer(I4P),                intent(IN)::  n            !< Element of the list to be found.
  logical,                     intent(OUT):: d            !< Link data.
  type(Type_SL_Link), pointer::              ln => null() !< Pointer to "n-th" link of the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = .false.
  ln => link(l=l,n=n)
  if (associated(ln%d)) d = transfer(ln%d,.false.)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_link_log

#ifdef r16p
  !> @brief Subroutine for getting a link data from the list head (R16P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_link_R16(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  real(R16P),          intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR16P
  if (associated(l%d)) d = transfer(l%d,1._R16P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_link_R16
#endif

  !> @brief Subroutine for getting a link data from the list head (R8P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_link_R8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  real(R8P),           intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR8P
  if (associated(l%d)) d = transfer(l%d,1._R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_link_R8

  !> @brief Subroutine for getting a link data from the list head (R4P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_link_R4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  real(R4P),           intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR4P
  if (associated(l%d)) d = transfer(l%d,1._R4P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_link_R4

  !> @brief Subroutine for getting a link data from the list head (I8P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_link_I8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  integer(I8P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI8P
  if (associated(l%d)) d = transfer(l%d,1_I8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_link_I8

  !> @brief Subroutine for getting a link data from the list head (I4P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_link_I4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  integer(I4P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI4P
  if (associated(l%d)) d = transfer(l%d,1_I4P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_link_I4

  !> @brief Subroutine for getting a link data from the list head (I2P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_link_I2(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  integer(I2P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI2P
  if (associated(l%d)) d = transfer(l%d,1_I2P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_link_I2

  !> @brief Subroutine for getting a link data from the list head (I1P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_link_I1(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  integer(I1P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI1P
  if (associated(l%d)) d = transfer(l%d,1_I1P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_link_I1

  !> @brief Subroutine for getting a link data from the list head (character).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_link_ch(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l    !< List.
  character(*),        intent(OUT):: d    !< Link data.
  character(len(d))::                mold !< Mold prototype.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = ''
  if (associated(l%d)) d = transfer(l%d,mold(1:min(len(d),size(l%d,dim=1))))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_link_ch

  !> @brief Subroutine for getting a link data from the list head (logical).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_link_log(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  logical,             intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = .false.
  if (associated(l%d)) d = transfer(l%d,.false.)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_link_log

#ifdef r16p
  !> @brief Subroutine for getting a link data from the list tail (R16P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  recursive subroutine gett_link_R16(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  real(R16P),          intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%n)) then
    d = MinR16P
    if (associated(l%d)) d = transfer(l%d,1._R16P)
  else
    call gett_link_R16(l=l%n,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_link_R16
#endif

  !> @brief Subroutine for getting a link data from the list tail (R8P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  recursive subroutine gett_link_R8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  real(R8P),           intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%n)) then
    d = MinR8P
    if (associated(l%d)) d = transfer(l%d,1._R8P)
  else
    call gett_link_R8(l=l%n,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_link_R8

  !> @brief Subroutine for getting a link data from the list tail (R4P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  recursive subroutine gett_link_R4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  real(R4P),           intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%n)) then
    d = MinR4P
    if (associated(l%d)) d = transfer(l%d,1._R4P)
  else
    call gett_link_R4(l=l%n,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_link_R4

  !> @brief Subroutine for getting a link data from the list tail (I8P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  recursive subroutine gett_link_I8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  integer(I8P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%n)) then
    d = MinI8P
    if (associated(l%d)) d = transfer(l%d,1_I8P)
  else
    call gett_link_I8(l=l%n,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_link_I8

  !> @brief Subroutine for getting a link data from the list tail (I4P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  recursive subroutine gett_link_I4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  integer(I4P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%n)) then
    d = MinI4P
    if (associated(l%d)) d = transfer(l%d,1_I4P)
  else
    call gett_link_I4(l=l%n,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_link_I4

  !> @brief Subroutine for getting a link data from the list tail (I2P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  recursive subroutine gett_link_I2(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  integer(I2P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%n)) then
    d = MinI2P
    if (associated(l%d)) d = transfer(l%d,1_I2P)
  else
    call gett_link_I2(l=l%n,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_link_I2

  !> @brief Subroutine for getting a link data from the list tail (I1P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  recursive subroutine gett_link_I1(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  integer(I1P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%n)) then
    d = MinI1P
    if (associated(l%d)) d = transfer(l%d,1_I1P)
  else
    call gett_link_I1(l=l%n,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_link_I1

  !> @brief Subroutine for getting a link data from the list tail (character).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  recursive subroutine gett_link_ch(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l    !< List.
  character(*),        intent(OUT):: d    !< Link data.
  character(len(d))::                mold !< Mold prototype.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%n)) then
    d = ''
    if (associated(l%d)) d = transfer(l%d,mold(1:min(len(d),size(l%d,dim=1))))
  else
    call gett_link_ch(l=l%n,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_link_ch

  !> @brief Subroutine for getting a link data from the list tail (logical).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  recursive subroutine gett_link_log(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(IN)::  l !< List.
  logical,             intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%n)) then
    d = .false.
    if (associated(l%d)) d = transfer(l%d,.false.)
  else
    call gett_link_log(l=l%n,d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_link_log

#ifdef r16p
  !> @brief Subroutine for inserting a link into the n-th link of list (R16P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_link_R16(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l                                !< List.
  integer(I4P),                intent(IN)::    n                                !< Element of the list to be inserted.
  real(R16P),                  intent(IN)::    d                                !< Data of the current link.
  type(Type_SL_Link), pointer::                lp=>null(),lc=>null(),ln=>null() !< Previous, current and new links.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(ln) ; allocate(ln%d(1:bit_size(d)/8)) ; ln%d = transfer(d,SL_List_Mold) ; ln%tp = tp_SR16
  lp => l%link(n=n-1)
  if (associated(lp)) then
    lc => lp%n
    lp%n => ln
    ln%n => lc
  else
    call l%putt(d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_link_R16
#endif

  !> @brief Subroutine for inserting a link into the n-th link of list (R8P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_link_R8(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l                                !< List.
  integer(I4P),                intent(IN)::    n                                !< Element of the list to be inserted.
  real(R8P),                   intent(IN)::    d                                !< Data of the current link.
  type(Type_SL_Link), pointer::                lp=>null(),lc=>null(),ln=>null() !< Previous, current and new links.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(ln) ; allocate(ln%d(1:bit_size(d)/8)) ; ln%d = transfer(d,SL_List_Mold) ; ln%tp = tp_SR8
  lp => l%link(n=n-1)
  if (associated(lp)) then
    lc => lp%n
    lp%n => ln
    ln%n => lc
  else
    call l%putt(d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_link_R8

  !> @brief Subroutine for inserting a link into the n-th link of list (R4P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_link_R4(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l                                !< List.
  integer(I4P),                intent(IN)::    n                                !< Element of the list to be inserted.
  real(R4P),                   intent(IN)::    d                                !< Data of the current link.
  type(Type_SL_Link), pointer::                lp=>null(),lc=>null(),ln=>null() !< Previous, current and new links.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(ln) ; allocate(ln%d(1:bit_size(d)/8)) ; ln%d = transfer(d,SL_List_Mold) ; ln%tp = tp_SR4
  lp => l%link(n=n-1)
  if (associated(lp)) then
    lc => lp%n
    lp%n => ln
    ln%n => lc
  else
    call l%putt(d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_link_R4

  !> @brief Subroutine for inserting a link into the n-th link of list (I8P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_link_I8(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l                                !< List.
  integer(I4P),                intent(IN)::    n                                !< Element of the list to be inserted.
  integer(I8P),                intent(IN)::    d                                !< Data of the current link.
  type(Type_SL_Link), pointer::                lp=>null(),lc=>null(),ln=>null() !< Previous, current and new links.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(ln) ; allocate(ln%d(1:bit_size(d)/8)) ; ln%d = transfer(d,SL_List_Mold) ; ln%tp = tp_SI8
  lp => l%link(n=n-1)
  if (associated(lp)) then
    lc => lp%n
    lp%n => ln
    ln%n => lc
  else
    call l%putt(d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_link_I8

  !> @brief Subroutine for inserting a link into the n-th link of list (I4P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_link_I4(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l                                !< List.
  integer(I4P),                intent(IN)::    n                                !< Element of the list to be inserted.
  integer(I4P),                intent(IN)::    d                                !< Data of the current link.
  type(Type_SL_Link), pointer::                lp=>null(),lc=>null(),ln=>null() !< Previous, current and new links.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(ln) ; allocate(ln%d(1:bit_size(d)/8)) ; ln%d = transfer(d,SL_List_Mold) ; ln%tp = tp_SI4
  lp => l%link(n=n-1)
  if (associated(lp)) then
    lc => lp%n
    lp%n => ln
    ln%n => lc
  else
    call l%putt(d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_link_I4

  !> @brief Subroutine for inserting a link into the n-th link of list (I2P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_link_I2(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l                                !< List.
  integer(I4P),                intent(IN)::    n                                !< Element of the list to be inserted.
  integer(I2P),                intent(IN)::    d                                !< Data of the current link.
  type(Type_SL_Link), pointer::                lp=>null(),lc=>null(),ln=>null() !< Previous, current and new links.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(ln) ; allocate(ln%d(1:bit_size(d)/8)) ; ln%d = transfer(d,SL_List_Mold) ; ln%tp = tp_SI2
  lp => l%link(n=n-1)
  if (associated(lp)) then
    lc => lp%n
    lp%n => ln
    ln%n => lc
  else
    call l%putt(d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_link_I2

  !> @brief Subroutine for inserting a link into the n-th link of list (I1P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_link_I1(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l                                !< List.
  integer(I4P),                intent(IN)::    n                                !< Element of the list to be inserted.
  integer(I1P),                intent(IN)::    d                                !< Data of the current link.
  type(Type_SL_Link), pointer::                lp=>null(),lc=>null(),ln=>null() !< Previous, current and new links.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(ln) ; allocate(ln%d(1:bit_size(d)/8)) ; ln%d = transfer(d,SL_List_Mold) ; ln%tp = tp_SI1
  lp => l%link(n=n-1)
  if (associated(lp)) then
    lc => lp%n
    lp%n => ln
    ln%n => lc
  else
    call l%putt(d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_link_I1

  !> @brief Subroutine for inserting a link into the n-th link of list (character).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_link_ch(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l                                !< List.
  integer(I4P),                intent(IN)::    n                                !< Element of the list to be inserted.
  character(*),                intent(IN)::    d                                !< Data of the current link.
  type(Type_SL_Link), pointer::                lp=>null(),lc=>null(),ln=>null() !< Previous, current and new links.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(ln) ; allocate(ln%d(1:len(d))) ; ln%d = transfer(d,SL_List_Mold) ; ln%tp = tp_SCH
  lp => l%link(n=n-1)
  if (associated(lp)) then
    lc => lp%n
    lp%n => ln
    ln%n => lc
  else
    call l%putt(d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_link_ch

  !> @brief Subroutine for inserting a link into the n-th link of list (logical).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_link_log(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), target, intent(INOUT):: l                                !< List.
  integer(I4P),                intent(IN)::    n                                !< Element of the list to be inserted.
  logical,                     intent(IN)::    d                                !< Data of the current link.
  type(Type_SL_Link), pointer::                lp=>null(),lc=>null(),ln=>null() !< Previous, current and new links.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(ln) ; allocate(ln%d(1:1)) ; ln%d = transfer(d,SL_List_Mold) ; ln%tp = tp_SLOG
  lp => l%link(n=n-1)
  if (associated(lp)) then
    lc => lp%n
    lp%n => ln
    ln%n => lc
  else
    call l%putt(d=d)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_link_log

#ifdef r16p
  !> @brief Subroutine for inserting a link data into the list head (R16P).
  subroutine puth_link_R16(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l               !< List.
  real(R16P),          intent(IN)::    d               !< Data of the current link.
  integer(I1P),       pointer::        dtmp(:)=>null() !< Temporary data.
  type(Type_SL_Link), pointer::        n=>null()       !< New link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%d)) then
    allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold)
  else
    allocate(dtmp(1:bit_size(l%d)/8)) ; dtmp = l%d ; l%d = transfer(d,SL_List_Mold)
    allocate(n) ; allocate(n%d(1:bit_size(dtmp)/8)) ; n%d = dtmp ; n%tp = l%tp ; n%n => l%n ; l%n => n
    deallocate(dtmp)
  endif
  l%tp = tp_SR16
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_link_R16
#endif

  !> @brief Subroutine for inserting a link data into the list head (R8P).
  subroutine puth_link_R8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l               !< List.
  real(R8P),           intent(IN)::    d               !< Data of the current link.
  integer(I1P),       pointer::        dtmp(:)=>null() !< Temporary data.
  type(Type_SL_Link), pointer::        n=>null()       !< New link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%d)) then
    allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold)
  else
    allocate(dtmp(1:bit_size(l%d)/8)) ; dtmp = l%d ; l%d = transfer(d,SL_List_Mold)
    allocate(n) ; allocate(n%d(1:bit_size(dtmp)/8)) ; n%d = dtmp ; n%tp = l%tp ; n%n => l%n ; l%n => n
    deallocate(dtmp)
  endif
  l%tp = tp_SR8
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_link_R8

  !> @brief Subroutine for inserting a link data into the list head (R4P).
  subroutine puth_link_R4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l               !< List.
  real(R4P),           intent(IN)::    d               !< Data of the current link.
  integer(I1P),       pointer::        dtmp(:)=>null() !< Temporary data.
  type(Type_SL_Link), pointer::        n=>null()       !< New link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%d)) then
    allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold)
  else
    allocate(dtmp(1:bit_size(l%d)/8)) ; dtmp = l%d ; l%d = transfer(d,SL_List_Mold)
    allocate(n) ; allocate(n%d(1:bit_size(dtmp)/8)) ; n%d = dtmp ; n%tp = l%tp ; n%n => l%n ; l%n => n
    deallocate(dtmp)
  endif
  l%tp = tp_SR4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_link_R4

  !> @brief Subroutine for inserting a link data into the list head (I8P).
  subroutine puth_link_I8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l               !< List.
  integer(I8P),        intent(IN)::    d               !< Data of the current link.
  integer(I1P),       pointer::        dtmp(:)=>null() !< Temporary data.
  type(Type_SL_Link), pointer::        n=>null()       !< New link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%d)) then
    allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold)
  else
    allocate(dtmp(1:bit_size(l%d)/8)) ; dtmp = l%d ; l%d = transfer(d,SL_List_Mold)
    allocate(n) ; allocate(n%d(1:bit_size(dtmp)/8)) ; n%d = dtmp ; n%tp = l%tp ; n%n => l%n ; l%n => n
    deallocate(dtmp)
  endif
  l%tp = tp_SI8
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_link_I8

  !> @brief Subroutine for inserting a link data into the list head (I4P).
  subroutine puth_link_I4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l               !< List.
  integer(I4P),        intent(IN)::    d               !< Data of the current link.
  integer(I1P),       pointer::        dtmp(:)=>null() !< Temporary data.
  type(Type_SL_Link), pointer::        n=>null()       !< New link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%d)) then
    allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold)
  else
    allocate(dtmp(1:bit_size(l%d)/8)) ; dtmp = l%d ; l%d = transfer(d,SL_List_Mold)
    allocate(n) ; allocate(n%d(1:bit_size(dtmp)/8)) ; n%d = dtmp ; n%tp = l%tp ; n%n => l%n ; l%n => n
    deallocate(dtmp)
  endif
  l%tp = tp_SI4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_link_I4

  !> @brief Subroutine for inserting a link data into the list head (I2P).
  subroutine puth_link_I2(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l               !< List.
  integer(I2P),        intent(IN)::    d               !< Data of the current link.
  integer(I1P),       pointer::        dtmp(:)=>null() !< Temporary data.
  type(Type_SL_Link), pointer::        n=>null()       !< New link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%d)) then
    allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold)
  else
    allocate(dtmp(1:bit_size(l%d)/8)) ; dtmp = l%d ; l%d = transfer(d,SL_List_Mold)
    allocate(n) ; allocate(n%d(1:bit_size(dtmp)/8)) ; n%d = dtmp ; n%tp = l%tp ; n%n => l%n ; l%n => n
    deallocate(dtmp)
  endif
  l%tp = tp_SI2
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_link_I2

  !> @brief Subroutine for inserting a link data into the list head (I1P).
  subroutine puth_link_I1(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l               !< List.
  integer(I1P),        intent(IN)::    d               !< Data of the current link.
  integer(I1P),       pointer::        dtmp(:)=>null() !< Temporary data.
  type(Type_SL_Link), pointer::        n=>null()       !< New link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%d)) then
    allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold)
  else
    allocate(dtmp(1:bit_size(l%d)/8)) ; dtmp = l%d ; l%d = transfer(d,SL_List_Mold)
    allocate(n) ; allocate(n%d(1:bit_size(dtmp)/8)) ; n%d = dtmp ; n%tp = l%tp ; n%n => l%n ; l%n => n
    deallocate(dtmp)
  endif
  l%tp = tp_SI1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_link_I1

  !> @brief Subroutine for inserting a link data into the list head (character).
  subroutine puth_link_ch(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l               !< List.
  character(*),        intent(IN)::    d               !< Data of the current link.
  integer(I1P),       pointer::        dtmp(:)=>null() !< Temporary data.
  type(Type_SL_Link), pointer::        n=>null()       !< New link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%d)) then
    allocate(l%d(1:len(d))) ; l%d = transfer(d,SL_List_Mold)
  else
    allocate(dtmp(1:bit_size(l%d)/8)) ; dtmp = l%d ; l%d = transfer(d,SL_List_Mold)
    deallocate(l%d) ; allocate(l%d(1:len(d))) ; l%d = transfer(d,SL_List_Mold)
    allocate(n) ; allocate(n%d(1:bit_size(dtmp)/8)) ; n%d = dtmp ; n%tp = l%tp ; n%n => l%n ; l%n => n
    deallocate(dtmp)
  endif
  l%tp = tp_SCH
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_link_ch

  !> @brief Subroutine for inserting a link data into the list head (logical).
  subroutine puth_link_log(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l               !< List.
  logical,             intent(IN)::    d               !< Data of the current link.
  integer(I1P),       pointer::        dtmp(:)=>null() !< Temporary data.
  type(Type_SL_Link), pointer::        n=>null()       !< New link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%d)) then
    allocate(l%d(1:1)) ; l%d = transfer(d,SL_List_Mold)
  else
    allocate(dtmp(1:bit_size(l%d)/8)) ; dtmp = l%d ; l%d = transfer(d,SL_List_Mold)
    allocate(n) ; allocate(n%d(1:bit_size(dtmp)/8)) ; n%d = dtmp ; n%tp = l%tp ; n%n => l%n ; l%n => n
    deallocate(dtmp)
  endif
  l%tp = tp_SLOG
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_link_log

#ifdef r16p
  !> @brief Recursive subroutine for inserting a link data into the list tail (R16P).
  recursive subroutine putt_link_R16(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  real(R16P),          intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%d)) then
    if (.not.associated(l%n)) allocate(l%n) ; call putt_link_R16(l=l%n,d=d)
  else
    if (.not.associated(l%d)) allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold) ; l%tp = tp_SR16
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_link_R16
#endif

  !> @brief Recursive subroutine for inserting a link data into the list tail (R8P).
  recursive subroutine putt_link_R8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  real(R8P),           intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%d)) then
    if (.not.associated(l%n)) allocate(l%n) ; call putt_link_R8(l=l%n,d=d)
  else
    if (.not.associated(l%d)) allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold) ; l%tp = tp_SR8
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_link_R8

  !> @brief Recursive subroutine for inserting a link data into the list tail (R4P).
  recursive subroutine putt_link_R4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  real(R4P),           intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%d)) then
    if (.not.associated(l%n)) allocate(l%n) ; call putt_link_R4(l=l%n,d=d)
  else
    if (.not.associated(l%d)) allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold) ; l%tp = tp_SR4
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_link_R4

  !> @brief Recursive subroutine for inserting a link data into the list tail (I18).
  recursive subroutine putt_link_I8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  integer(I8P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%d)) then
    if (.not.associated(l%n)) allocate(l%n) ; call putt_link_I8(l=l%n,d=d)
  else
    if (.not.associated(l%d)) allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold) ; l%tp = tp_SI8
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_link_I8

  !> @brief Recursive subroutine for inserting a link data into the list tail (I4P).
  recursive subroutine putt_link_I4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  integer(I4P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%d)) then
    if (.not.associated(l%n)) allocate(l%n) ; call putt_link_I4(l=l%n,d=d)
  else
    if (.not.associated(l%d)) allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold) ; l%tp = tp_SI4
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_link_I4

  !> @brief Recursive subroutine for inserting a link data into the list tail (I2P).
  recursive subroutine putt_link_I2(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  integer(I2P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%d)) then
    if (.not.associated(l%n)) allocate(l%n) ; call putt_link_I2(l=l%n,d=d)
  else
    if (.not.associated(l%d)) allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold) ; l%tp = tp_SI2
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_link_I2

  !> @brief Recursive subroutine for inserting a link data into the list tail (I1P).
  recursive subroutine putt_link_I1(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  integer(I1P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%d)) then
    if (.not.associated(l%n)) allocate(l%n) ; call putt_link_I1(l=l%n,d=d)
  else
    if (.not.associated(l%d)) allocate(l%d(1:bit_size(d)/8)) ; l%d = transfer(d,SL_List_Mold) ; l%tp = tp_SI1
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_link_I1

  !> @brief Recursive subroutine for inserting a link data into the list tail (character).
  recursive subroutine putt_link_ch(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  character(*),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%d)) then
    if (.not.associated(l%n)) allocate(l%n) ; call putt_link_ch(l=l%n,d=d)
  else
    if (.not.associated(l%d)) allocate(l%d(1:len(d))) ; l%d = transfer(d,SL_List_Mold) ; l%tp = tp_SCH
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_link_ch

  !> @brief Recursive subroutine for inserting a link data into the list tail (logical).
  recursive subroutine putt_link_log(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link), intent(INOUT):: l !< List.
  logical,             intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%d)) then
    if (.not.associated(l%n)) allocate(l%n) ; call putt_link_log(l=l%n,d=d)
  else
    if (.not.associated(l%d)) allocate(l%d(1:1)) ; l%d = transfer(d,SL_List_Mold) ; l%tp = tp_SLOG
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_link_log

#ifdef r16p
  !> @brief Subroutine for converting list to array (R16P).
  subroutine array_link_R16(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link),     intent(IN)::  l    ! List.
  real(R16P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  type(Type_SL_Link), pointer::          n    ! Pointer for scanning the list.
  integer(I1P),       pointer::          d(:) ! Pointer for scanning the list.
  integer(I4P)::                         Nl   ! Elements number.
  integer(I4P)::                         e    ! Elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = l%leng()
  if (allocated(a)) deallocate(a)
  if (Nl>0) then
    allocate(a(1:Nl))
    e = 0
    d => l%d
    n => l%n
    do while (associated(d))
      e = e + 1
      a(e) = transfer(d,1._R16P)
      if (associated(n)) then
        d => n%d
        n => n%n
      else
        d => null()
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_link_R16
#endif

  !> @brief Subroutine for converting list to array (R8P).
  subroutine array_link_R8(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link),    intent(IN)::  l    ! List.
  real(R8P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  type(Type_SL_Link), pointer::         n    ! Pointer for scanning the list.
  integer(I1P),       pointer::         d(:) ! Pointer for scanning the list.
  integer(I4P)::                        Nl   ! Elements number.
  integer(I4P)::                        e    ! Elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = l%leng()
  if (allocated(a)) deallocate(a)
  if (Nl>0) then
    allocate(a(1:Nl))
    e = 0
    d => l%d
    n => l%n
    do while (associated(d))
      e = e + 1
      a(e) = transfer(d,1._R8P)
      if (associated(n)) then
        d => n%d
        n => n%n
      else
        d => null()
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_link_R8

  !> @brief Subroutine for converting list to array (R4P).
  subroutine array_link_R4(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link),    intent(IN)::  l    ! List.
  real(R4P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  type(Type_SL_Link), pointer::         n    ! Pointer for scanning the list.
  integer(I1P),       pointer::         d(:) ! Pointer for scanning the list.
  integer(I4P)::                        Nl   ! Elements number.
  integer(I4P)::                        e    ! Elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = l%leng()
  if (allocated(a)) deallocate(a)
  if (Nl>0) then
    allocate(a(1:Nl))
    e = 0
    d => l%d
    n => l%n
    do while (associated(d))
      e = e + 1
      a(e) = transfer(d,1._R4P)
      if (associated(n)) then
        d => n%d
        n => n%n
      else
        d => null()
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_link_R4

  !> @brief Subroutine for converting list to array (I8P).
  subroutine array_link_I8(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link),       intent(IN)::  l    ! List.
  integer(I8P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  type(Type_SL_Link), pointer::            n    ! Pointer for scanning the list.
  integer(I1P),       pointer::            d(:) ! Pointer for scanning the list.
  integer(I4P)::                           Nl   ! Elements number.
  integer(I4P)::                           e    ! Elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = l%leng()
  if (allocated(a)) deallocate(a)
  if (Nl>0) then
    allocate(a(1:Nl))
    e = 0
    d => l%d
    n => l%n
    do while (associated(d))
      e = e + 1
      a(e) = transfer(d,1_I8P)
      if (associated(n)) then
        d => n%d
        n => n%n
      else
        d => null()
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_link_I8

  !> @brief Subroutine for converting list to array (I4P).
  subroutine array_link_I4(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link),       intent(IN)::  l    ! List.
  integer(I4P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  type(Type_SL_Link), pointer::            n    ! Pointer for scanning the list.
  integer(I1P),       pointer::            d(:) ! Pointer for scanning the list.
  integer(I4P)::                           Nl   ! Elements number.
  integer(I4P)::                           e    ! Elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = l%leng()
  if (allocated(a)) deallocate(a)
  if (Nl>0) then
    allocate(a(1:Nl))
    e = 0
    d => l%d
    n => l%n
    do while (associated(d))
      e = e + 1
      a(e) = transfer(d,1_I4P)
      if (associated(n)) then
        d => n%d
        n => n%n
      else
        d => null()
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_link_I4

  !> @brief Subroutine for converting list to array (I2P).
  subroutine array_link_I2(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link),       intent(IN)::  l    ! List.
  integer(I2P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  type(Type_SL_Link), pointer::            n    ! Pointer for scanning the list.
  integer(I1P),       pointer::            d(:) ! Pointer for scanning the list.
  integer(I4P)::                           Nl   ! Elements number.
  integer(I4P)::                           e    ! Elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = l%leng()
  if (allocated(a)) deallocate(a)
  if (Nl>0) then
    allocate(a(1:Nl))
    e = 0
    d => l%d
    n => l%n
    do while (associated(d))
      e = e + 1
      a(e) = transfer(d,1_I2P)
      if (associated(n)) then
        d => n%d
        n => n%n
      else
        d => null()
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_link_I2

  !> @brief Subroutine for converting list to array (I1P).
  subroutine array_link_I1(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link),       intent(IN)::  l    ! List.
  integer(I1P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  type(Type_SL_Link), pointer::            n    ! Pointer for scanning the list.
  integer(I1P),       pointer::            d(:) ! Pointer for scanning the list.
  integer(I4P)::                           Nl   ! Elements number.
  integer(I4P)::                           e    ! Elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = l%leng()
  if (allocated(a)) deallocate(a)
  if (Nl>0) then
    allocate(a(1:Nl))
    e = 0
    d => l%d
    n => l%n
    do while (associated(d))
      e = e + 1
      a(e) = transfer(d,1_I1P)
      if (associated(n)) then
        d => n%d
        n => n%n
      else
        d => null()
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_link_I1

  !> @brief Subroutine for converting list to array (character).
  subroutine array_link_ch(l,Nc,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link),        intent(IN)::  l    ! List.
  integer(I4P),               intent(IN)::  Nc   ! Number of character of each list element.
  character(Nc), allocatable, intent(OUT):: a(:) ! Array containing the list.
  character(Nc)::                           mold ! Mold prototype.
  type(Type_SL_Link), pointer::             n    ! Pointer for scanning the list.
  integer(I1P),       pointer::             d(:) ! Pointer for scanning the list.
  integer(I4P)::                            Nl   ! Elements number.
  integer(I4P)::                            e    ! Elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = l%leng()
  if (allocated(a)) deallocate(a)
  if (Nl>0) then
    allocate(a(1:Nl))
    e = 0
    d => l%d
    n => l%n
    do while (associated(d))
      e = e + 1
      a(e) = transfer(d,mold)
      if (associated(n)) then
        d => n%d
        n => n%n
      else
        d => null()
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_link_ch

  !> @brief Subroutine for converting list to array (logical).
  subroutine array_link_log(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Link),  intent(IN)::  l    ! List.
  logical, allocatable, intent(OUT):: a(:) ! Array containing the list.
  type(Type_SL_Link), pointer::       n    ! Pointer for scanning the list.
  integer(I1P),       pointer::       d(:) ! Pointer for scanning the list.
  integer(I4P)::                      Nl   ! Elements number.
  integer(I4P)::                      e    ! Elements counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = l%leng()
  if (allocated(a)) deallocate(a)
  if (Nl>0) then
    allocate(a(1:Nl))
    e = 0
    d => l%d
    n => l%n
    do while (associated(d))
      e = e + 1
      a(e) = transfer(d,.false.)
      if (associated(n)) then
        d => n%d
        n => n%n
      else
        d => null()
      endif
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_link_log

  !> @brief Subroutine for freeing (destroying) the list.
  subroutine free_list(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  l%l = 0
  l%homo = .true.
  if (associated(l%f)) then
    call l%f%free
    deallocate(l%f)
    l%f => null()
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_list

  !> @brief Function for computing the length of the list.
  function leng_list(l) result(Nl)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN):: l    ! List.
  integer(I4P)::                    Nl   ! Elements number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nl = 0
  if (associated(l%f)) Nl = l%f%leng()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction leng_list

  !> @brief Recursive subroutine for checking if the list is homogeneous.
  function homo_list(l) result(homo)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN):: l    !< List.
  logical::                         homo !< Flag for checking if the list is homogeneous.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  homo = .true.
  if (associated(l%f)) homo = l%f%homogeneous()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction homo_list

  !> @brief Subroutine for deleting n-th link of the list.
  subroutine del_list(l,n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  integer(I4P),        intent(IN)::    n !< Element of the list to be deleted.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    call l%f%del(n=n)
    l%l = l%l - 1
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine del_list

  !> @brief Subroutine for deleting head link of the list.
  subroutine delh_list(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    call l%f%delh
    l%l = l%l - 1
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine delh_list

  !> @brief Subroutine for deleting tail link of the list.
  subroutine delt_list(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    call l%f%delt
    l%l = l%l - 1
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine delt_list

#ifdef r16p
  !> @brief Subroutine for getting data from the n-th link of list (R16P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_list_R16(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(IN)::  n !< Element of the list to be found.
  real(R16P),          intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR16P
  if (associated(l%f)) call l%f%get(n=n,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_list_R16
#endif

  !> @brief Subroutine for getting data from the n-th link of list (R8P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_list_R8(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(IN)::  n !< Element of the list to be found.
  real(R8P),           intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR8P
  if (associated(l%f)) call l%f%get(n=n,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_list_R8

  !> @brief Subroutine for getting data from the n-th link of list (R4P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_list_R4(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(IN)::  n !< Element of the list to be found.
  real(R4P),           intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR4P
  if (associated(l%f)) call l%f%get(n=n,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_list_R4

  !> @brief Subroutine for getting data from the n-th link of list (I8P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_list_I8(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(IN)::  n !< Element of the list to be found.
  integer(I8P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI8P
  if (associated(l%f)) call l%f%get(n=n,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_list_I8

  !> @brief Subroutine for getting data from the n-th link of list (I4P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_list_I4(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(IN)::  n !< Element of the list to be found.
  integer(I4P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI4P
  if (associated(l%f)) call l%f%get(n=n,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_list_I4

  !> @brief Subroutine for getting data from the n-th link of list (I2P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_list_I2(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(IN)::  n !< Element of the list to be found.
  integer(I2P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI2P
  if (associated(l%f)) call l%f%get(n=n,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_list_I2

  !> @brief Subroutine for getting data from the n-th link of list (I1P).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_list_I1(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(IN)::  n !< Element of the list to be found.
  integer(I1P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI1P
  if (associated(l%f)) call l%f%get(n=n,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_list_I1

  !> @brief Subroutine for getting data from the n-th link of list (character).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_list_ch(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(IN)::  n !< Element of the list to be found.
  character(*),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = ''
  if (associated(l%f)) call l%f%get(n=n,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_list_ch

  !> @brief Subroutine for getting data from the n-th link of list (logical).
  !> @note If n<1 then the first link data is returned, whereas if n> list's length the last data is returned.
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine get_list_log(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(IN)::  n !< Element of the list to be found.
  logical,             intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = .false.
  if (associated(l%f)) call l%f%get(n=n,d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_list_log

#ifdef r16p
  !> @brief Subroutine for getting a link data from the list head (R16P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_list_R16(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Liist, intent(IN)::  l !< List.
  real(R16P),          intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR16P
  if (associated(l%f)) call l%f%geth(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_list_R16
#endif

  !> @brief Subroutine for getting a link data from the list head (R8P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_list_R8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  real(R8P),           intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR8P
  if (associated(l%f)) call l%f%geth(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_list_R8

  !> @brief Subroutine for getting a link data from the list head (R4P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_list_R4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  real(R4P),           intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR4P
  if (associated(l%f)) call l%f%geth(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_list_R4

  !> @brief Subroutine for getting a link data from the list head (I8P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_list_I8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I8P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI8P
  if (associated(l%f)) call l%f%geth(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_list_I8

  !> @brief Subroutine for getting a link data from the list head (I4P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_list_I4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI4P
  if (associated(l%f)) call l%f%geth(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_list_I4

  !> @brief Subroutine for getting a link data from the list head (I2P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_list_I2(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I2P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI2P
  if (associated(l%f)) call l%f%geth(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_list_I2

  !> @brief Subroutine for getting a link data from the list head (I1P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_list_I1(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I1P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI1P
  if (associated(l%f)) call l%f%geth(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_list_I1

  !> @brief Subroutine for getting a link data from the list head (character).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_list_ch(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  character(*),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = ''
  if (associated(l%f)) call l%f%geth(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_list_ch

  !> @brief Subroutine for getting a link data from the list head (logical).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine geth_list_log(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  logical,             intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = .false.
  if (associated(l%f)) call l%f%geth(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine geth_list_log

#ifdef r16p
  !> @brief Subroutine for getting a link data from the list tail (R16P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine gett_list_R16(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_Liist, intent(IN)::  l !< List.
  real(R16P),          intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR16P
  if (associated(l%f)) call l%f%gett(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_list_R16
#endif

  !> @brief Subroutine for getting a link data from the list tail (R8P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine gett_list_R8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  real(R8P),           intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR8P
  if (associated(l%f)) call l%f%gett(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_list_R8

  !> @brief Subroutine for getting a link data from the list tail (R4P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine gett_list_R4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  real(R4P),           intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinR4P
  if (associated(l%f)) call l%f%gett(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_list_R4

  !> @brief Subroutine for getting a link data from the list tail (I8P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine gett_list_I8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I8P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI8P
  if (associated(l%f)) call l%f%gett(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_list_I8

  !> @brief Subroutine for getting a link data from the list tail (I4P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine gett_list_I4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I4P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI4P
  if (associated(l%f)) call l%f%gett(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_list_I4

  !> @brief Subroutine for getting a link data from the list tail (I2P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine gett_list_I2(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I2P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI2P
  if (associated(l%f)) call l%f%gett(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_list_I2

  !> @brief Subroutine for getting a link data from the list tail (I1P).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine gett_list_I1(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  integer(I1P),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = MinI1P
  if (associated(l%f)) call l%f%gett(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_list_I1

  !> @brief Subroutine for getting a link data from the list tail (character).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine gett_list_ch(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  character(*),        intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = ''
  if (associated(l%f)) call l%f%gett(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_list_ch

  !> @brief Subroutine for getting a link data from the list tail (logical).
  !> @note If data is not present the minimum representable value is returned (for error checking).
  subroutine gett_list_log(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(IN)::  l !< List.
  logical,             intent(OUT):: d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = .false.
  if (associated(l%f)) call l%f%gett(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine gett_list_log

#ifdef r16p
  !> @brief Subroutine for inserting a link into the n-th link of list (R16P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_list_R16(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT)::  l !< List.
  integer(I4P),        intent(IN)::     n !< Element of the list to be inserted.
  real(R16P),          intent(IN)::     d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%put(n=n,d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_list_R16
#endif

  !> @brief Subroutine for inserting a link into the n-th link of list (R8P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_list_R8(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT)::  l !< List.
  integer(I4P),        intent(IN)::     n !< Element of the list to be inserted.
  real(R8P),           intent(IN)::     d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%put(n=n,d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_list_R8

  !> @brief Subroutine for inserting a link into the n-th link of list (R4P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_list_R4(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT)::  l !< List.
  integer(I4P),        intent(IN)::     n !< Element of the list to be inserted.
  real(R4P),           intent(IN)::     d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%put(n=n,d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_list_R4

  !> @brief Subroutine for inserting a link into the n-th link of list (I8P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_list_I8(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT)::  l !< List.
  integer(I4P),        intent(IN)::     n !< Element of the list to be inserted.
  integer(I8P),        intent(IN)::     d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%put(n=n,d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_list_I8

  !> @brief Subroutine for inserting a link into the n-th link of list (I4P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_list_I4(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT)::  l !< List.
  integer(I4P),        intent(IN)::     n !< Element of the list to be inserted.
  integer(I4P),        intent(IN)::     d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%put(n=n,d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_list_I4

  !> @brief Subroutine for inserting a link into the n-th link of list (I2P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_list_I2(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT)::  l !< List.
  integer(I4P),        intent(IN)::     n !< Element of the list to be inserted.
  integer(I2P),        intent(IN)::     d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%put(n=n,d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_list_I2

  !> @brief Subroutine for inserting a link into the n-th link of list (I1P).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_list_I1(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT)::  l !< List.
  integer(I4P),        intent(IN)::     n !< Element of the list to be inserted.
  integer(I1P),        intent(IN)::     d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%put(n=n,d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_list_I1

  !> @brief Subroutine for inserting a link into the n-th link of list (character).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_list_ch(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT)::  l !< List.
  integer(I4P),        intent(IN)::     n !< Element of the list to be inserted.
  character(*),        intent(IN)::     d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%put(n=n,d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_list_ch

  !> @brief Subroutine for inserting a link into the n-th link of list (logical).
  !> @note If n<1 then the data is inserted at first link , whereas if n> list's length the data is inserted at last.
  subroutine put_list_log(l,n,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT)::  l !< List.
  integer(I4P),        intent(IN)::     n !< Element of the list to be inserted.
  logical,             intent(IN)::     d !< Link data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%put(n=n,d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine put_list_log

#ifdef r16p
  !> @brief Subroutine for inserting a link data into the list head (R16P).
  subroutine puth_list_R16(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  real(R16P),          intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%puth(d=d)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_list_R16
#endif

  !> @brief Subroutine for inserting a link data into the list head (R8P).
  subroutine puth_list_R8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  real(R8P),           intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%puth(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_list_R8

  !> @brief Subroutine for inserting a link data into the list head (R4P).
  subroutine puth_list_R4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  real(R4P),           intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%puth(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_list_R4

  !> @brief Subroutine for inserting a link data into the list head (I8P).
  subroutine puth_list_I8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  integer(I8P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%puth(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_list_I8

  !> @brief Subroutine for inserting a link data into the list head (I4P).
  subroutine puth_list_I4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  integer(I4P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%puth(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_list_I4

  !> @brief Subroutine for inserting a link data into the list head (I2P).
  subroutine puth_list_I2(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  integer(I2P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%puth(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_list_I2

  !> @brief Subroutine for inserting a link data into the list head (I1P).
  subroutine puth_list_I1(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  integer(I1P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%puth(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_list_I1

  !> @brief Subroutine for inserting a link data into the list head (character).
  subroutine puth_list_ch(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  character(*),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%puth(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_list_ch

  !> @brief Subroutine for inserting a link data into the list head (logical).
  subroutine puth_list_log(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  logical,             intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%puth(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine puth_list_log

#ifdef r16p
  !> @brief Subroutine for inserting a link data into the list tail (R16P).
  subroutine putt_list_R16(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  real(R16P),          intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%putt(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_list_R16
#endif

  !> @brief Subroutine for inserting a link data into the list tail (R8P).
  subroutine putt_list_R8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  real(R8P),           intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%putt(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_list_R8

  !> @brief Subroutine for inserting a link data into the list tail (R4P).
  subroutine putt_list_R4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  real(R4P),           intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%putt(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_list_R4

  !> @brief Subroutine for inserting a link data into the list tail (I8P).
  subroutine putt_list_I8(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  integer(I8P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%putt(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_list_I8

  !> @brief Subroutine for inserting a link data into the list tail (I4P).
  subroutine putt_list_I4(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  integer(I4P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%putt(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_list_I4

  !> @brief Subroutine for inserting a link data into the list tail (I2P).
  subroutine putt_list_I2(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  integer(I2P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%putt(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_list_I2

  !> @brief Subroutine for inserting a link data into the list tail (I1P).
  subroutine putt_list_I1(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  integer(I1P),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%putt(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_list_I1

  !> @brief Subroutine for inserting a link data into the list tail (character).
  subroutine putt_list_ch(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  character(*),        intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%putt(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_list_ch

  !> @brief Subroutine for inserting a link data into the list tail (logical).
  subroutine putt_list_log(l,d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List), intent(INOUT):: l !< List.
  logical,             intent(IN)::    d !< Data of the current link.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(l%f)) allocate(l%f)
  call l%f%putt(d=d)
  l%l = l%l + 1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine putt_list_log

#ifdef r16p
  !> @brief Subroutine for converting list to array (R16P).
  subroutine array_list_R16(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List),     intent(IN)::  l    ! List.
  real(R16P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    if (l%f%homogeneous()) call l%f%array(a=a)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_list_R16
#endif

  !> @brief Subroutine for converting list to array (R8P).
  subroutine array_list_R8(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List),    intent(IN)::  l    ! List.
  real(R8P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    if (l%f%homogeneous()) call l%f%array(a=a)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_list_R8

  !> @brief Subroutine for converting list to array (R4P).
  subroutine array_list_R4(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List),    intent(IN)::  l    ! List.
  real(R4P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    if (l%f%homogeneous()) call l%f%array(a=a)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_list_R4

  !> @brief Subroutine for converting list to array (I8P).
  subroutine array_list_I8(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List),       intent(IN)::  l    ! List.
  integer(I8P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    if (l%f%homogeneous()) call l%f%array(a=a)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_list_I8

  !> @brief Subroutine for converting list to array (I4P).
  subroutine array_list_I4(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List),       intent(IN)::  l    ! List.
  integer(I4P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    if (l%f%homogeneous()) call l%f%array(a=a)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_list_I4

  !> @brief Subroutine for converting list to array (I2P).
  subroutine array_list_I2(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List),       intent(IN)::  l    ! List.
  integer(I2P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    if (l%f%homogeneous()) call l%f%array(a=a)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_list_I2

  !> @brief Subroutine for converting list to array (I1P).
  subroutine array_list_I1(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List),       intent(IN)::  l    ! List.
  integer(I1P), allocatable, intent(OUT):: a(:) ! Array containing the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    if (l%f%homogeneous()) call l%f%array(a=a)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_list_I1

  !> @brief Subroutine for converting list to array (character).
  subroutine array_list_ch(l,Nc,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List),        intent(IN)::  l    ! List.
  integer(I4P),               intent(IN)::  Nc   ! Number of character of each list element.
  character(Nc), allocatable, intent(OUT):: a(:) ! Array containing the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    if (l%f%homogeneous()) call l%f%array(Nc=Nc,a=a)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_list_ch

  !> @brief Subroutine for converting list to array (logical).
  subroutine array_list_log(l,a)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SL_List),  intent(IN)::  l    ! List.
  logical, allocatable, intent(OUT):: a(:) ! Array containing the list.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(l%f)) then
    if (l%f%homogeneous()) call l%f%array(a=a)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine array_list_log
  !> @}
endmodule Data_Type_SL_List
