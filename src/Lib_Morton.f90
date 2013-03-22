!> @ingroup Library
!> @{
!> @defgroup Lib_MortonLibrary Lib_Morton
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Lib_MortonInterface Lib_Morton
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup Lib_MortonPrivateVarPar Lib_Morton
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_MortonPrivateProcedure Lib_Morton
!> @}

!> This module contains procedures and data for encoding (and deconding) integers into Morton's (or z-order) ordering. Morton
!> (Z-ordering) is a scheme to map multi-dimensional arrays onto to a linear with a great deal of spatial locality (see A Computer
!> Oriented Geodetic Data Base and a New Technique in File Sequencing, Morton G.M., technical report, IBM, 1966). Two main
!> procedures are available: one for encoding (decoding) 2 integer(I1P/I2P/I4P) into 1 integer(I8P) and the second for encoding
!> (decoding) 3 integer(I1P/I2P) into 1 integer(I8P). The econding is performed by means of the Morton's algorithm: e.g. if two
!> 8 bits integers, 'x' and 'y', must be encoded into 1 integer(I8P) the Morton's key is 'y7x7y6x6y5x5y4x4y3x3y2x2y1x2y0x0' where
!> the numbers indicate the bits. This order, called Z-order or Morton order, preserves important proprieties like the spatial
!> locality of multi-refined hierarchical mesh.
!> This is a library module.
!> @ingroup Lib_MortonLibrary
module Lib_Morton
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition. module Lib_Morton
!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: dil,ctc
public:: morton2,demorton2,morton3,demorton3
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Lib_MortonPrivateVarPar
!> @{
! Binary masks: used into the bits dilating and contracting algorithms.
integer(I8P), parameter:: mask32_32 = Z'FFFFFFFF'         !< 0000000000000000000000000000000011111111111111111111111111111111.
integer(I8P), parameter:: mask16_48 = Z'FFFF'             !< 0000000000000000000000000000000000000000000000001111111111111111.
integer(I8P), parameter:: mask16_32 = Z'FFFF00000000FFFF' !< 1111111111111111000000000000000000000000000000001111111111111111.
integer(I8P), parameter:: mask16_16 = Z'FFFF0000FFFF'     !< 0000000000000000111111111111111100000000000000001111111111111111.
integer(I8P), parameter:: mask8_56  = Z'FF'               !< 0000000000000000000000000000000000000000000000000000000011111111.
integer(I8P), parameter:: mask8_16  = Z'FF0000FF0000FF'   !< 0000000011111111000000000000000011111111000000000000000011111111.
integer(I8P), parameter:: mask8_8   = Z'FF00FF00FF00FF'   !< 0000000011111111000000001111111100000000111111110000000011111111.
integer(I8P), parameter:: mask4_60  = Z'F'                !< 0000000000000000000000000000000000000000000000000000000000001111.
integer(I8P), parameter:: mask4_8   = Z'F00F00F00F00F00F' !< 1111000000001111000000001111000000001111000000001111000000001111.
integer(I8P), parameter:: mask4_4   = Z'F0F0F0F0F0F0F0F'  !< 0000111100001111000011110000111100001111000011110000111100001111.
integer(I8P), parameter:: mask2_62  = Z'3'                !< 0000000000000000000000000000000000000000000000000000000000000011.
integer(I8P), parameter:: mask2_4   = z'30C30C30C30C30C3' !< 0011000011000011000011000011000011000011000011000011000011000011.
integer(I8P), parameter:: mask2_2   = Z'3333333333333333' !< 0011001100110011001100110011001100110011001100110011001100110011.
integer(I8P), parameter:: mask1_2   = Z'9249249249249249' !< 1001001001001001001001001001001001001001001001001001001001001001.
integer(I8P), parameter:: mask1_1   = Z'5555555555555555' !< 0101010101010101010101010101010101010101010101010101010101010101.
!> Binary mask for selecting significant bits.
integer(I8P), parameter:: signif(1:5) = (/mask2_62,mask4_60,mask8_56,mask16_48,mask32_32/)
!> Binary mask for perfoming significant bits shifting.
integer(I8P), parameter:: mask(1:6,1:2) = reshape((/mask1_1,mask2_2,mask4_4,mask8_8, mask16_16,mask32_32,   & ! 1 bit interleaving.
                                                    mask1_2,mask2_4,mask4_8,mask8_16,mask16_32,mask32_32/), & ! 2 bits interleaving.
                                                   (/6,2/))
!> Shift number array.
integer(I1P), parameter:: shft(1:6,1:2) = reshape((/1,2,4,8, 16,32,&          ! 1 bit interleaving.
                                                    2,4,8,16,32,64/),(/6,2/)) ! 2 bits interleaving.
!> Real parameter for computing the number of shifts (Ns). The number of shifts is computed by \f$2^{Ns}=b\f$ where 'b' is the
!> significant bits. As a consequence \f$Ns=\frac{log(b)}{log2}\f$, therefore it is convenient to store the value of
!> \f$\frac{1}{log2}\f$.
real(R4P), parameter:: log10_2_inv = 1._R4P/log10(2._R4P)
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Function for dilating 1 integer(I1P/I2P/I4P) (ref On Spatial Orders and Location Codes, Stocco, LJ and Schrack, G,
!> IEEE Transaction on Computers, vol 58, n 3, March 2009). The resulting integer has 64 bits; it has only 'b' significant
!> bits interleaved by 'z' zeros: bb/zx0/bb-1/zx0../b1/zx0/b0; e.g. for (b=4,z=1): b3/b2/b1/b0 => b3/0/b2/0/b1/0/b0.
!> @ingroup Lib_MortonInterface
interface dil
  module procedure dil_I1P,dil_I2P,dil_I4P
endinterface
!> @brief Subroutine for contracting 1 integer of 64 bits (ref On Spatial Orders and Location Codes, Stocco, LJ and Schrack,
!> G, IEEE Transaction on Computers, vol 58, n 3, March 2009). The resulting integer(I1P/I2P/I4P) has only 'b' significant bits
!> obtained by the following contraction: if bb/zx0/bb-1/zx0../b1/zx0/b0 => bb/bb-1/.../b1/b0; e.g. for (b=4,z=1):
!> b3/0/b2/0/b1/0/b0 => b3/b2/b1/b0.
!> @ingroup Lib_MortonInterface
interface ctc
  module procedure ctc_I1P,ctc_I2P,ctc_I4P
endinterface
!> @brief Function for encoding 2 integer(I1P/I2P/I4P) indexes into 1 integer(I8P) Morton's key.
!> @ingroup Lib_MortonInterface
interface morton2
  module procedure morton2_I1P,morton2_I2P,morton2_I4P
endinterface
!> @brief Subroutine for decoding 1 integer(I8P) Morton's key into 2 integer(I1P/I2P/I4P) indexes.
!> @ingroup Lib_MortonInterface
interface demorton2
  module procedure demorton2_I1P,demorton2_I2P,demorton2_I4P
endinterface
!> @brief Function for encoding 3 integer(I1P/I2P) indexes into 1 integer(I8P) Morton's key.
!> @ingroup Lib_MortonInterface
interface morton3
  module procedure morton3_I1P,morton3_I2P
endinterface
!> @brief Subroutine for decoding 1 integer(I8P) Morton's key into 3 integer(I1P/I2P) indexes.
!> @ingroup Lib_MortonInterface
interface demorton3
  module procedure demorton3_I1P,demorton3_I2P
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_MortonPrivateProcedure
  !> @{
  !> @brief Function for dilating integers (ref On Spatial Orders and Location Codes, Stocco, LJ and Schrack, G,
  !> IEEE Transaction on Computers, vol 58, n 3, March 2009). The resulting integer has 64 bits; it has only 'b' significant
  !> bits interleaved by 'z' zeros: bb/zx0/bb-1/zx0../b1/zx0/b0; e.g. for (b=4,z=1): b3/b2/b1/b0 => b3/0/b2/0/b1/0/b0.
  elemental function dil_I1P(i,b,z) result(d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN):: i   !< Input integer.
  integer(I2P), intent(IN):: b   !< Number of significant bits of 'i' (2/4/8).
  integer(I1P), intent(IN):: z   !< Number of zero 'i' (1/2).
  integer(I8P)::             d   !< Dilated integer.
  integer(I1P)::             l,m !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  l = int(log10(b*1._R4P)*log10_2_inv,I1P)
  d = int(i,I8P)
  d = iand(d,signif(l))
  do m=l,1,-1 !(3/2/1,1,-1)
    d = iand(ior(d,ishft(d,shft(m,z))),mask(m,z))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dil_I1P

  !> @brief Function for dilating integers (ref On Spatial Orders and Location Codes, Stocco, LJ and Schrack, G,
  !> IEEE Transaction on Computers, vol 58, n 3, March 2009). The resulting integer has 64 bits; it has only 'b' significant
  !> bits interleaved by 'z' zeros: bb/zx0/bb-1/zx0../b1/zx0/b0; e.g. for (b=4,z=1): b3/b2/b1/b0 => b3/0/b2/0/b1/0/b0.
  elemental function dil_I2P(i,b,z) result(d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P), intent(IN):: i   !< Input integer.
  integer(I2P), intent(IN):: b   !< Number of significant bits of 'i' (2/4/8/16).
  integer(I1P), intent(IN):: z   !< Number of zero 'i' (1/2).
  integer(I8P)::             d   !< Dilated integer.
  integer(I1P)::             l,m !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  l = int(log10(b*1._R4P)*log10_2_inv,I1P)
  d = int(i,I8P)
  d = iand(d,signif(l))
  do m=l,1,-1 !(4/3/2/1,1,-1)
    d = iand(ior(d,ishft(d,shft(m,z))),mask(m,z))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dil_I2P

  !> @brief Function for dilating integers (ref On Spatial Orders and Location Codes, Stocco, LJ and Schrack, G,
  !> IEEE Transaction on Computers, vol 58, n 3, March 2009). The resulting integer has 64 bits; it has only 'b' significant
  !> bits interleaved by 'z' zeros: bb/zx0/bb-1/zx0../b1/zx0/b0; e.g. for (b=4,z=1): b3/b2/b1/b0 => b3/0/b2/0/b1/0/b0.
  elemental function dil_I4P(i,b,z) result(d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: i   !< Input integer.
  integer(I2P), intent(IN):: b   !< Number of significant bits of 'i' (2/4/8/16/32).
  integer(I1P), intent(IN):: z   !< Number of zero 'i' (1/2).
  integer(I8P)::             d   !< Dilated integer.
  integer(I1P)::             l,m !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  l = int(log10(b*1._R4P)*log10_2_inv,I1P)
  d = int(i,I8P)
  d = iand(d,signif(l))
  do m=l,1,-1 !(5/4/3/2/1,1,-1)
    d = iand(ior(d,ishft(d,shft(m,z))),mask(m,z))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dil_I4P

  !> @brief Subroutine for contracting integer of 64 bits (ref On Spatial Orders and Location Codes, Stocco, LJ and Schrack,
  !> G, IEEE Transaction on Computers, vol 58, n 3, March 2009). The resulting integer has only 'b' significant bits obtained by
  !> the following contraction: if bb/zx0/bb-1/zx0../b1/zx0/b0 => bb/bb-1/.../b1/b0; e.g. for (b=4,z=1):
  !> b3/0/b2/0/b1/0/b0 => b3/b2/b1/b0.
  elemental subroutine ctc_I1P(i,b,z,c)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN)::  i !< Input integer.
  integer(I2P), intent(IN)::  b !< Number of significant bits of 'i' (2/4/8).
  integer(I1P), intent(IN)::  z !< Number of zero 'i' (1/2).
  integer(I1P), intent(OUT):: c !< Contracted integer.
  integer(I8P)::              d !< Temporary dilated integer.
  integer(I1P)::              m !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = iand(i,mask(1,z))
  do m=1,int(log10(b*1._R4P)*log10_2_inv,I1P),1 !(1,3/2/1)
    d = iand(ieor(d,ishft(d,-shft(m,z))),mask(m+1,z))
  enddo
  c = int(d,I1P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine ctc_I1P

  !> @brief Subroutine for contracting integer of 64 bits (ref On Spatial Orders and Location Codes, Stocco, LJ and Schrack,
  !> G, IEEE Transaction on Computers, vol 58, n 3, March 2009). The resulting integer has only 'b' significant bits obtained by
  !> the following contraction: if bb/zx0/bb-1/zx0../b1/zx0/b0 => bb/bb-1/.../b1/b0; e.g. for (b=4,z=1):
  !> b3/0/b2/0/b1/0/b0 => b3/b2/b1/b0.
  elemental subroutine ctc_I2P(i,b,z,c)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN)::  i !< Input integer.
  integer(I2P), intent(IN)::  b !< Number of significant bits of 'i' (2/4/8/16).
  integer(I1P), intent(IN)::  z !< Number of zero 'i' (1/2).
  integer(I2P), intent(OUT):: c !< Contracted integer.
  integer(I8P)::              d !< Temporary dilated integer.
  integer(I1P)::              m !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = iand(i,mask(1,z))
  do m=1,int(log10(b*1._R4P)*log10_2_inv,I1P),1 !(1,4/3/2/1)
    d = iand(ieor(d,ishft(d,-shft(m,z))),mask(m+1,z))
  enddo
  c = int(d,I2P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine ctc_I2P

  !> @brief Subroutine for contracting integer of 64 bits (ref On Spatial Orders and Location Codes, Stocco, LJ and Schrack,
  !> G, IEEE Transaction on Computers, vol 58, n 3, March 2009). The resulting integer has only 'b' significant bits obtained by
  !> the following contraction: if bb/zx0/bb-1/zx0../b1/zx0/b0 => bb/bb-1/.../b1/b0; e.g. for (b=4,z=1):
  !> b3/0/b2/0/b1/0/b0 => b3/b2/b1/b0.
  elemental subroutine ctc_I4P(i,b,z,c)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN)::  i !< Input integer.
  integer(I2P), intent(IN)::  b !< Number of significant bits of 'i' (2/4/8/16/32).
  integer(I1P), intent(IN)::  z !< Number of zero 'i' (1/2).
  integer(I4P), intent(OUT):: c !< Contracted integer.
  integer(I8P)::              d !< Temporary dilated integer.
  integer(I1P)::              m !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d = iand(i,mask(1,z))
  do m=1,int(log10(b*1._R4P)*log10_2_inv,I1P),1 !(1,5/4/3/2/1)
    d = iand(ieor(d,ishft(d,-shft(m,z))),mask(m+1,z))
  enddo
  c = int(d,I4P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine ctc_I4P

  !> @brief Function for encoding 2 integer(I1P) indexes into 1 integer(I8P) Morton's key.
  elemental function morton2_I1P(i,j,b) result(key)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN)::           i     !< I index.
  integer(I1P), intent(IN)::           j     !< J index.
  integer(I2P), intent(IN), optional:: b     !< Number of significant bits of 'i' (2/4/8).
  integer(I8P)::                       key   !  Morton's key.
  integer(I8P)::                       di,dj !< Dilated indexes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    di = dil_I1P(i=i,b=b,z=1_I1P)
    dj = dil_I1P(i=j,b=b,z=1_I1P)
  else
    di = dil_I1P(i=i,b=8_I2P,z=1_I1P)
    dj = dil_I1P(i=j,b=8_I2P,z=1_I1P)
  endif
  key = ishft(dj,1) + di
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction morton2_I1P

  !> @brief Function for encoding 2 integer(I2P) indexes into 1 integer(I8P) Morton's key.
  elemental function morton2_I2P(i,j,b) result(key)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P), intent(IN)::           i     !< I index.
  integer(I2P), intent(IN)::           j     !< J index.
  integer(I2P), intent(IN), optional:: b     !< Number of significant bits of 'i' (2/4/8/16).
  integer(I8P)::                       key   !  Morton's key.
  integer(I8P)::                       di,dj !< Dilated indexes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    di = dil_I2P(i=i,b=b,z=1_I1P)
    dj = dil_I2P(i=j,b=b,z=1_I1P)
  else
    di = dil_I2P(i=i,b=16_I2P,z=1_I1P)
    dj = dil_I2P(i=j,b=16_I2P,z=1_I1P)
  endif
  key = ishft(dj,1) + di
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction morton2_I2P

  !> @brief Function for encoding 2 integer(I4P) indexes into 1 integer(I8P) Morton's key.
  elemental function morton2_I4P(i,j,b) result(key)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN)::           i     !< I index.
  integer(I4P), intent(IN)::           j     !< J index.
  integer(I2P), intent(IN), optional:: b     !< Number of significant bits of 'i' (2/4/8/16/32).
  integer(I8P)::                       key   !  Morton's key.
  integer(I8P)::                       di,dj !< Dilated indexes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    di = dil_I4P(i=i,b=b,z=1_I1P)
    dj = dil_I4P(i=j,b=b,z=1_I1P)
  else
    di = dil_I4P(i=i,b=32_I2P,z=1_I1P)
    dj = dil_I4P(i=j,b=32_I2P,z=1_I1P)
  endif
  key = ishft(dj,1) + di
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction morton2_I4P

  !> @brief Subroutine for decoding 1 integer(I8P) Morton's key into 2 integer(I1P) indexes.
  elemental subroutine demorton2_I1P(key,i,j,b)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN)::           key !< Morton's key.
  integer(I1P), intent(INOUT)::        i   !< I index.
  integer(I1P), intent(INOUT)::        j   !< J index.
  integer(I2P), intent(IN), optional:: b   !< Number of significant bits of 'i' (2/4/8).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    call ctc(i=key,          b=b,z=1_I1P,c=i)
    call ctc(i=ishft(key,-1),b=b,z=1_I1P,c=j)
  else
    call ctc(i=key,          b=8_I2P,z=1_I1P,c=i)
    call ctc(i=ishft(key,-1),b=8_I2P,z=1_I1P,c=j)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine demorton2_I1P

  !> @brief Subroutine for decoding 1 integer(I8P) Morton's key into 2 integer(I2P) indexes.
  elemental subroutine demorton2_I2P(key,i,j,b)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN)::           key !< Morton's key.
  integer(I2P), intent(INOUT)::        i   !< I index.
  integer(I2P), intent(INOUT)::        j   !< J index.
  integer(I2P), intent(IN), optional:: b   !< Number of significant bits of 'i' (2/4/8/16).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    call ctc(i=key,          b=b,z=1_I1P,c=i)
    call ctc(i=ishft(key,-1),b=b,z=1_I1P,c=j)
  else
    call ctc(i=key,          b=16_I2P,z=1_I1P,c=i)
    call ctc(i=ishft(key,-1),b=16_I2P,z=1_I1P,c=j)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine demorton2_I2P

  !> @brief Subroutine for decoding 1 integer(I8P) Morton's key into 2 integer(I4P) indexes.
  elemental subroutine demorton2_I4P(key,i,j,b)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN)::           key !< Morton's key.
  integer(I4P), intent(INOUT)::        i   !< I index.
  integer(I4P), intent(INOUT)::        j   !< J index.
  integer(I2P), intent(IN), optional:: b   !< Number of significant bits of 'i' (2/4/8/16/32).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    call ctc(i=key,          b=b,z=1_I1P,c=i)
    call ctc(i=ishft(key,-1),b=b,z=1_I1P,c=j)
  else
    call ctc(i=key,          b=32_I2P,z=1_I1P,c=i)
    call ctc(i=ishft(key,-1),b=32_I2P,z=1_I1P,c=j)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine demorton2_I4P

  !> @brief Function for encoding 3 integer(I1P) indexes into 1 integer(I8P) Morton's key.
  elemental function morton3_I1P(i,j,k,b) result(key)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN)::           i        !< I index.
  integer(I1P), intent(IN)::           j        !< J index.
  integer(I1P), intent(IN)::           k        !< K index.
  integer(I2P), intent(IN), optional:: b        !< Number of significant bits of 'i' (2/4/8).
  integer(I8P)::                       key      !  Morton's key.
  integer(I8P)::                       di,dj,dk !< Dilated indexes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    di = dil_I1P(i=i,b=b,z=2_I1P)
    dj = dil_I1P(i=j,b=b,z=2_I1P)
    dk = dil_I1P(i=k,b=b,z=2_I1P)
  else
    di = dil_I1P(i=i,b=8_I2P,z=2_I1P)
    dj = dil_I1P(i=j,b=8_I2P,z=2_I1P)
    dk = dil_I1P(i=k,b=8_I2P,z=2_I1P)
  endif
  !key = ishft(dk,2) + ishft(dj,1) + di
  key = ishft(dj,2) + ishft(dk,1) + di
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction morton3_I1P

  !> @brief Function for encoding 3 integer(I2P) indexes into 1 integer(I8P) Morton's key.
  elemental function morton3_I2P(i,j,k,b) result(key)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P), intent(IN)::           i        !< I index.
  integer(I2P), intent(IN)::           j        !< J index.
  integer(I2P), intent(IN)::           k        !< K index.
  integer(I2P), intent(IN), optional:: b        !< Number of significant bits of 'i' (2/4/8/16).
  integer(I8P)::                       key      !  Morton's key.
  integer(I8P)::                       di,dj,dk !< Dilated indexes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    di = dil_I2P(i=i,b=b,z=2_I1P)
    dj = dil_I2P(i=j,b=b,z=2_I1P)
    dk = dil_I2P(i=k,b=b,z=2_I1P)
  else
    di = dil_I2P(i=i,b=16_I2P,z=2_I1P)
    dj = dil_I2P(i=j,b=16_I2P,z=2_I1P)
    dk = dil_I2P(i=k,b=16_I2P,z=2_I1P)
  endif
  !key = ishft(dk,2) + ishft(dj,1) + di
  key = ishft(dj,2) + ishft(dk,1) + di
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction morton3_I2P

  !> @brief Subroutine for decoding 1 integer(I8P) Morton's key into 3 integer(I1P) indexes.
  elemental subroutine demorton3_I1P(key,i,j,k,b)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN)::           key !< Morton's key.
  integer(I1P), intent(INOUT)::        i   !< I index.
  integer(I1P), intent(INOUT)::        j   !< J index.
  integer(I1P), intent(INOUT)::        k   !< K index.
  integer(I2P), intent(IN), optional:: b   !< Number of significant bits of 'i' (2/4/8).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    !call ctc(i=key,          b=b,z=2_I1P,c=i)
    !call ctc(i=ishft(key,-1),b=b,z=2_I1P,c=j)
    !call ctc(i=ishft(key,-2),b=b,z=2_I1P,c=k)
    call ctc(i=key,          b=b,z=2_I1P,c=i)
    call ctc(i=ishft(key,-1),b=b,z=2_I1P,c=k)
    call ctc(i=ishft(key,-2),b=b,z=2_I1P,c=j)
  else
    !call ctc(i=key,          b=8_I2P,z=2_I1P,c=i)
    !call ctc(i=ishft(key,-1),b=8_I2P,z=2_I1P,c=j)
    !call ctc(i=ishft(key,-2),b=8_I2P,z=2_I1P,c=k)
    call ctc(i=key,          b=8_I2P,z=2_I1P,c=i)
    call ctc(i=ishft(key,-1),b=8_I2P,z=2_I1P,c=k)
    call ctc(i=ishft(key,-2),b=8_I2P,z=2_I1P,c=j)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine demorton3_I1P

  !> @brief Subroutine for decoding 1 integer(I8P) Morton's key into 3 integer(I2P) indexes.
  elemental subroutine demorton3_I2P(key,i,j,k,b)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN)::           key !< Morton's key.
  integer(I2P), intent(INOUT)::        i   !< I index.
  integer(I2P), intent(INOUT)::        j   !< J index.
  integer(I2P), intent(INOUT)::        k   !< K index.
  integer(I2P), intent(IN), optional:: b   !< Number of significant bits of 'i' (2/4/8/16).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) then
    !call ctc(i=key,          b=b,z=2_I1P,c=i)
    !call ctc(i=ishft(key,-1),b=b,z=2_I1P,c=j)
    !call ctc(i=ishft(key,-2),b=b,z=2_I1P,c=k)
    call ctc(i=key,          b=b,z=2_I1P,c=i)
    call ctc(i=ishft(key,-1),b=b,z=2_I1P,c=k)
    call ctc(i=ishft(key,-2),b=b,z=2_I1P,c=j)
  else
    !call ctc(i=key,          b=16_I2P,z=2_I1P,c=i)
    !call ctc(i=ishft(key,-1),b=16_I2P,z=2_I1P,c=j)
    !call ctc(i=ishft(key,-2),b=16_I2P,z=2_I1P,c=k)
    call ctc(i=key,          b=16_I2P,z=2_I1P,c=i)
    call ctc(i=ishft(key,-1),b=16_I2P,z=2_I1P,c=k)
    call ctc(i=ishft(key,-2),b=16_I2P,z=2_I1P,c=j)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine demorton3_I2P
  !> @}
endmodule Lib_Morton
