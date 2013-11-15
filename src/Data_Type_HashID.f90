!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_HashIDDerivedType Data_Type_HashID
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_HashIDInterface Data_Type_HashID
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_HashIDGlobalVarPar Data_Type_HashID
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup Data_Type_HashIDPrivateVarPar Data_Type_HashID
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_HashIDPublicProcedure Data_Type_HashID
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_HashIDPrivateProcedure Data_Type_HashID
!> @}

!> @brief Module Data_Type_HashID contains the definition of Type_HashID type and useful procedures for its handling.
!> Type_HashID contains the definition of the unique key, ID, used within  the hash tables.
module Data_Type_HashID
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
USE Lib_Morton   ! Procedure for Morton's encoding.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: ht_def_leng
public:: lvl_max
public:: treedim
public:: operator (/=)
public:: operator (==)
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Data_Type_HashIDGlobalVarPar
!> @{
integer(I4P), parameter:: ht_def_leng = 9973_I4P !< Default lenght of hash table.
integer(I1P)::            lvl_max_u = 0_I1P      !< Maximum refinement level used (maximum possible is lvl_max).
!> @}
!> @ingroup Data_Type_HashIDPrivateVarPar
integer(I1P), parameter:: lvl_max = 16_I1P       !< Maximum refinement level, arbitrarily chosen.
!> @brief Derived type containing the definition of the unique key, ID, used within the hash tables.
!> @note The identifier is a function of the form ID=(blk,bcl,lvl,pth) where:
!> @code
!> 'blk' are the bits identifying block;
!> 'bcl' are the bits identifying the (parent) base cell;
!> 'lvl' are the bits identifying the level;
!> 'pth' are the bits describing the path (of children).
!> @endcode
!> The range of ID keys can be computed as following. Each base cell can
!> contain a maximum of \f$ Nchild=\sum_{l=1}^{lvl_{max_u}}{2^{3 \cdot l}}\f$ children cells. As a consequence the maximum number of
!> cells are \f$ N_b \cdot N_{bc} \cdot Nchild \f$ where Nb is the (maximum) number of blocks, Nbc is the (maximum) number of the
!> base cells. Assuming (arbitrarily) \f$ N_b = 100 \f$, \f$ N_{bc} = 1 \cdot 10^6 \f$ and \f$ lvl_{max_u} = 5 \f$ the maximum
!> number of keys (cells) is of the order of \f$ 3 \cdot 10^{12}\f$, a very large value.
!> @code
!>                  +----+----+
!>                 /|   /|   /|
!>                / |  6 |  7 |
!>               /  +----*----+
!>              /  /|/  /|/  /|
!>             /  / |  4 |  5 |
!>            /  / /+----+----+
!>           /  / // / // /  /
!>          +----+----+/ /  /
!>         /| / /| / /| /  /
!>        / |/ //|/ //|/  /
!>       /  +----*----+  /
!>      /  /|// /|// /| /
!>     /  / |/ / |/ / |/
!>    /  / /+----+----+
!>   /  / // / // /  /
!>  +----+----+/ /  /       y(j)  z(k)
!>  | /  | /  | /  /        ^    ^
!>  |/ 2/|/ 3/|/  /         |   /
!>  +----*----+  /          |  /
!>  | /  | /  | /           | /
!>  |/ 0 |/ 1 |/            |/
!>  +----+----+             +------->x(i)
!> @endcode
!> @ingroup Data_Type_HashIDDerivedType
type, public:: Type_HashID
  integer(I2P):: blk = 0_I2P !< Block 16 bits: 65536 maximum blocks.
  integer(I8P):: bcl = 0_I4P !< Base cells 48(3x16) bits: 281474976710656 maximum cells (Morton encoded).
  integer(I1P):: lvl = 0_I1P !< Level 8 bits: 256 maximum levels, but arbitrarily fixed to 16 (due to path memory limit).
  integer(I8P):: pth = 0_I8P !< Path (of children) 64 bits: (4 bits, values in [0,1,2...,7])x(16 max levels), Morton encoded.
  contains
    procedure:: build ! Procedure for building ID from indexes information.
    procedure:: hash  ! Procedure for "hashing" ID.
endtype Type_HashID
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Not-equal-to boolean operator (/=) overloading.
!> @ingroup Data_Type_HashIDInterface
interface operator (/=)
  module procedure id_not_eq_id
endinterface
!> @brief Equal-to boolean operator (==) overloading.
!> @ingroup Data_Type_HashIDInterface
interface operator (==)
  module procedure id_eq_id
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_HashIDPublicProcedure
  !> @{
  !> @brief Function for computing the dimension of a tree.
  !> @note A tree can contain a maximum of \f$ Nchild=\sum_{l=1}^{lvl_{max_u}}{2^{s \cdot l}}\f$ children cells, where \f$s\f$ is
  !> the space dimension (1 in 1D, 2 in 2D, 3 in 3D ecc...).
  elemental function treedim(s,lm) result(d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN):: s  !< Space dimension (1=>1D, 2=>2D, 3=>3D...).
  integer(I1P), intent(IN):: lm !< Maximum level of the tree.
  integer(I8P)::             d  !< Dimension of the (complete) tree of dimension 's' and level 'lm'.
  integer(I1P)::             l  !< Levels counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = 0_I8P
  do l=1,lm
    d = d + (2**(s*l))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction treedim
  !> @}

  !> @ingroup Data_Type_HashIDPrivateProcedure
  !> @{
  !> @brief Subroutine for building ID from block,bcell-ijk,level informations.
  elemental subroutine build(ID,b,i,j,k,l,p)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_HashID), intent(INOUT):: ID !< Resulting ID.
  integer(I2P),       intent(IN)::    b  !< Block number (in the global ordering).
  integer(I2P),       intent(IN)::    i  !< Base cell i index (in the local-block ordering).
  integer(I2P),       intent(IN)::    j  !< Base cell j index (in the local-block ordering).
  integer(I2P),       intent(IN)::    k  !< Base cell k index (in the local-block ordering).
  integer(I1P),       intent(IN)::    l  !< Level (>0 if the key is of a child).
  integer(I8P),       intent(IN)::    p  !< Path (>0 if the key is of a child).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ID%blk = b
  ID%bcl = morton3(i=i,j=j,k=k)
  ID%lvl = l
  ID%pth = p
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine build

  !> @brief Function for computing the hash function of a given ID and hash table length.
  elemental function hash(id,tleng) result(bucket)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_HashID), intent(IN):: id     !< ID.
  integer(I4P),       intent(IN):: tleng  !< Hash table length.
  integer(I4P)::                   bucket !< Bucket corresponding to the ID.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  bucket = mod(int(id%blk,I4P),tleng)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction hash

  !> @brief This function returns .true. if id1 is /= id2, .false. otherwise.
  elemental function id_not_eq_id(id1,id2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_HashID), intent(IN):: id1     ! First ID.
  type(Type_HashID), intent(IN):: id2     ! Second ID.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = ((id1%blk/=id2%blk).OR.(id1%blk/=id2%blk))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction id_not_eq_id

  !> @brief This function returns .true. if id1 is 0= id2, .false. otherwise.
  elemental function id_eq_id(id1,id2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_HashID), intent(IN):: id1     ! First ID.
  type(Type_HashID), intent(IN):: id2     ! Second ID.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = ((id1%blk==id2%blk).AND.(id1%blk==id2%blk))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction id_eq_id
  !> @}
endmodule Data_Type_HashID
