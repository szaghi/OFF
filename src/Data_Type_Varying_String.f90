!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Varying_StringDerivedType Data_Type_Varying_String
!> Module definition of Type_Varying_String
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Varying_StringPrivateProcedure Data_Type_Varying_String
!> Module definition of Type_Varying_String
!> @}

!> @brief Module Data_Type_Varying_String contains the definition of Type_Varying_String and useful procedures for its handling.
module Data_Type_Varying_String
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: len,len_trim,adjustl,adjustr,trim
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type defining Type_Varying_String, a useful type for implementing ISO_VARYING_STRING by means of fortran 2003
!> features.
!> @ingroup Data_Type_Varying_StringDerivedType
type, public:: Type_Varying_String
  character(len=:), allocatable:: vs !< Deferred (vaying) length string.
  contains
    procedure:: free     ! Procedure for freeing dynamic memory.
    final::     finalize ! Procedure for freeing dynamic memory when finalizing.
    ! operators overloading
    generic:: read(formatted)    => read_fmt_self
    generic:: read(unformatted)  => read_ufmt_self
    generic:: write(formatted)   => write_fmt_self
    generic:: write(unformatted) => write_ufmt_self
    generic:: len                => len_self
    generic:: len_trim           => len_trim_self
    generic:: adjustl            => adjustl_self
    generic:: adjustr            => adjustr_self
    generic:: trim               => trim_self
    generic:: assignment(=)      => assign_self,assign_character,character_assign
    generic:: operator(//)       => concat_self,self_concat_character,character_concat_self
    ! private procedures
    procedure, pass(self ), private:: read_fmt_self
    procedure, pass(self ), private:: read_ufmt_self
    procedure, pass(self ), private:: write_fmt_self
    procedure, pass(self ), private:: write_ufmt_self
    procedure, pass(self ), private:: len_self
    procedure, pass(self ), private:: len_trim_self
    procedure, pass(self ), private:: adjustl_self
    procedure, pass(self ), private:: adjustr_self
    procedure, pass(self ), private:: trim_self
    procedure, pass(self1), private:: assign_self
    procedure, pass(self ), private:: assign_character
    procedure, pass(self ), private:: character_assign
    procedure, pass(self1), private:: concat_self
    procedure, pass(self ), private:: self_concat_character
    procedure, pass(self ), private:: character_concat_self
endtype Type_Varying_String
interface len
    module procedure len_self
endinterface len
interface len_trim
    module procedure len_trim_self
endinterface len_trim
interface adjustl
    module procedure adjustl_self
endinterface adjustl
interface adjustr
    module procedure adjustr_self
endinterface adjustr
interface trim
    module procedure trim_self
endinterface trim
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_XML_TagPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free(vstring)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(INOUT):: vstring !< Varying (lenght) string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(vstring%vs)) deallocate(vstring%vs)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(vstring)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Varying_String), intent(INOUT):: vstring !< XML tag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call vstring%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  ! IO Procedures
  !> @brief Procedure for formatted reading self.
  subroutine read_fmt_self(self,unit,iotype,vlist,iostat,iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(INOUT):: self
  integer,                    intent(IN)::    unit
  character(*),               intent(IN)::    iotype
  integer,                    intent(IN)::    vlist(:)
  integer,                    intent(OUT)::   iostat
  character(*),               intent(INOUT):: iomsg
  character(len=:), allocatable::             fm
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  iomsg=iotype
  if (size(vlist)>0) then
    self%vs = repeat(' ',vlist(1))
    fm = '(A'//trim(str(.true.,vlist(1)))//')'
  else
    if (.not.allocated(self%vs)) then
      iostat=1
      iomsg=' Error: self%vs in not allocated'
      return
    endif
    fm = '(A)'
  endif
  read(unit=unit,fmt=fm,iostat=iostat,iomsg=iomsg)self%vs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine read_fmt_self

  !> @brief Procedure for unformatted reading self.
  subroutine read_ufmt_self(self,unit,iostat,iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(INOUT):: self
  integer,                    intent(IN)::    unit
  integer,                    intent(OUT)::   iostat
  character(*),               intent(INOUT):: iomsg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit=unit,iostat=iostat,iomsg=iomsg)self%vs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine read_ufmt_self

  !> @brief Procedure for formatted writing self.
  subroutine write_fmt_self(self,unit,iotype,vlist,iostat,iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(IN)::    self
  integer,                    intent(IN)::    unit
  character(*),               intent(IN)::    iotype
  integer,                    intent(IN)::    vlist(:)
  integer,                    intent(OUT)::   iostat
  character(*),               intent(INOUT):: iomsg
  character(len=:), allocatable::             fm
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  iomsg=iotype
  fm ='(A)' ; if (size(vlist)>0) fm = '(A'//trim(str(.true.,vlist(1)))//')'
  write(unit=unit,fmt=fm,iostat=iostat,iomsg=iomsg)self%vs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine write_fmt_self

  !> @brief Procedure for unformatted writing self.
  subroutine write_ufmt_self(self,unit,iostat,iomsg)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(IN)::    self
  integer,                    intent(IN)::    unit
  integer,                    intent(OUT)::   iostat
  character(*),               intent(INOUT):: iomsg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit=unit,iostat=iostat,iomsg=iomsg)self%vs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine write_ufmt_self

  !> @brief Procedure for computing the length of self.
  elemental function len_self(self) result(leng)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(IN):: self
  integer::                                leng
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  leng = 0 ; if (allocated(self%vs)) leng = len(self%vs)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction len_self

  !> @brief Procedure for computing the length of trimmed self.
  elemental function len_trim_self(self) result(leng)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(IN):: self
  integer::                                leng
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  leng = 0 ; if (allocated(self%vs)) leng = len_trim(self%vs)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction len_trim_self

  !> @brief Procedure for adjusting self to the left, removing leading blanks and inserting trailing blanks.
  elemental function adjustl_self(self) result(adj)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(IN):: self
  type(Type_Varying_String)::              adj
  integer::                                leng
  integer::                                fws
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  leng = len(self)
  if (leng>0) then
    adj%vs = self%vs
    fws = index(self%vs,' ')
    if (fws>0) adj%vs = self%vs(fws+1:)//repeat(' ',fws)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction adjustl_self

  !> @brief Procedure for adjusting self to the right, removing trailing blanks and inserting leading blanks.
  elemental function adjustr_self(self) result(adj)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(IN):: self
  type(Type_Varying_String)::              adj
  integer::                                leng
  integer::                                lengt
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  leng = len(self)
  if (leng>0) then
    adj%vs = self%vs
    lengt = len_trim(self)
    if (leng>lengt) adj%vs = repeat(' ',leng-lengt)//self%vs(1:lengt)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction adjustr_self

  !> @brief Procedure for removing trailing blanks.
  elemental function trim_self(self) result(selft)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(IN):: self
  type(Type_Varying_String)::              selft
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%vs)) selft%vs = self%vs(1:len_trim(self%vs))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction trim_self

  ! Assignment (=)
  !> @brief Procedure for assignment between two selfs.
  elemental subroutine assign_self(self1,self2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(INOUT):: self1
  type(Type_Varying_String),  intent(IN)::    self2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self2%vs)) self1%vs = self2%vs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_self

  !> @brief Procedure for assignment between self and character string.
  elemental subroutine assign_character(self,car)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(INOUT):: self
  character(*),               intent(IN)::    car
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%vs = car
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_character

  !> @brief Procedure for assignment between character string and self.
  elemental subroutine character_assign(car,self)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),               intent(INOUT):: car
  class(Type_Varying_String), intent(IN)::    self
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  car = self%vs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine character_assign

  ! Concatenation (//)
  !> @brief Procedure for concatenating two selfs.
  elemental function concat_self(self1,self2) result(concat)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(IN):: self1
  type(Type_Varying_String),  intent(IN):: self2
  type(Type_Varying_String)::              concat
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  concat%vs = self1%vs//self2%vs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction concat_self

  !> @brief Procedure for concatenating self and character string.
  elemental function self_concat_character(self,car) result(concat)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(IN):: self
  character(*),               intent(IN):: car
  type(Type_Varying_String)::              concat
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  concat%vs = self%vs//car
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_concat_character

  !> @brief Procedure for concatenating character string and self.
  elemental function character_concat_self(car,self) result(concat)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),               intent(IN):: car
  class(Type_Varying_String), intent(IN):: self
  type(Type_Varying_String)::              concat
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  concat%vs = car//self%vs
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction character_concat_self
  !> @}
endmodule Data_Type_Varying_String
