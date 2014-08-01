!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_ProcmapDerivedType Data_Type_File_Procmap
!> Module definition of Type_File_Procmap
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_File_ProcmapInterface Data_Type_File_Procmap
!> Module definition of Type_File_Procmap
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_File_ProcmapPublicProcedure Data_Type_File_Procmap
!> Module definition of Type_File_Procmap
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_ProcmapPrivateProcedure Data_Type_File_Procmap
!> Module definition of Type_File_Procmap
!> @}

!> @brief Module Data_Type_File_Procmap contains the definition of Type_File_Procmap, that is the processes/blocks map.
!> @ingroup Data_Type_File_ProcmapDerivedType
module Data_Type_File_Procmap
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                      ! Integers and reals precision definition.
USE Data_Type_File_Base,       only: Type_File_Base,err_Nproc_unmatch ! Definition of Type_File_Base.
USE Data_Type_Mesh_Dimensions, only: Type_Mesh_Dimensions             ! Definition of Type_Mesh_Dimensions.
USE Data_Type_Parallel,        only: Type_Parallel                    ! Definition of Type_Parallel.
USE Data_Type_XML_Tag,         only: Type_XML_Tag                     ! Definition of Type_XML_Tag.
USE Lib_IO_Misc,               only: iostat_eor                       ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Procmap.
!> @ingroup Data_Type_File_ProcmapDerivedType
type, public, extends(Type_File_Base):: Type_File_Procmap
  contains
    procedure:: load  => load_procmap ! Procedure for loading procmap file.
    procedure:: save  => save_procmap ! Procedure for saving procmap file.
endtype Type_File_Procmap
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_ProcmapPrivateProcedure
  !> @{
  !> @brief Procedure for loading the processes/blocks map and local/global blocks map.
  !> @note This options file is a very simple XML file. This file contains the following data:
  !>   - mesh_dim%Nb;
  !>   - parallel%Nproc (that must be equals to the one defined by the MPI environment if used);
  !>   - parallel%procmap;
  !>   - parallel%blockmap.
  !> @note This options file is formatted by means of XML syntax. A validated file should be formatted as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <block_process b="1"      p="..."/>
  !>   <block_process b="2"      p="..."/>
  !>   ...
  !>   <block_process b="Nb_tot" p="..."/>
  !> @endcode
  !> The main tag is 'block_process' that contains the blocks map, namely the association of each block (in global numeration) to
  !> the processes handling it. The processes numeration must start from "0".
  !> The above 'block_process' tags can appear in any order as well as their attributes. The total number of blocks, "Nb_tot" is
  !> computed parsing the whole file as well as the totale number of processes "Nproc" that must be equal to the one defined by the
  !> MPI envinroment. It is worth nothing that the syntax of tag names and attributes is case sensitive, whereas the syntax of
  !> their values is case insensitive.
  subroutine load_procmap(file_d,mesh_dims,parallel)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Procmap),   intent(INOUT):: file_d    !< File data.
  type(Type_Mesh_Dimensions), intent(INOUT):: mesh_dims !< File data.
  type(Type_Parallel),        intent(INOUT):: parallel  !< Parallel data.
  type(Type_XML_Tag)::                        tag       !< Dummy XML tag for parsing file.
  integer(I4P)::                              Nproc     !< Number of MPI processes.
  integer(I4P)::                              Nb_tot    !< Global number of blocks (all processes).
  integer(I4P)::                              b,bb      !< Blocks counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%open(ascii=.true.,action='READ') ; if (file_d%iostat/=0) return
  allocate(character(len=1000):: tag%string%vs)
  call tag%set(tag_name='block_process',att_name=['b','p'])
  associate(Nb=>mesh_dims%Nb,unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg)
    Nb_tot = 0
    Blocks_Count: do
      read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=10)tag%string%vs ; if (iostat/=0.and.iostat/=iostat_eor) return
      if (index(string=tag%string%vs,substring='<block_process')>0) Nb_tot = Nb_tot + 1
    enddo Blocks_Count
    10 rewind(unit=unit)
    call parallel%alloc(Nb_tot = Nb_tot) ! allocating procmap
    Read_Loop: do
      read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=20)tag%string%vs ; if (iostat/=0.and.iostat/=iostat_eor) return
      if (index(string=tag%string%vs,substring='<block_process')>0) then
        call tag%get_attributes
        b                   = cton(str=tag%att_val(1)%vs,knd=1_I4P)
        parallel%procmap(b) = cton(str=tag%att_val(2)%vs,knd=1_I4P)
      endif
    enddo Read_Loop
    20 call file_d%close
    Nproc = maxval(parallel%procmap,dim=1)+1
    if (Nproc/=parallel%Nproc) then
      ! the number of processes defined into the file is different from the one defined by the MPI envinroment
      call file_d%raise_error(errtype=err_Nproc_unmatch)
    endif
    ! computing the local/global blocks map
    if (parallel%Nproc==1_I4P) then
      ! there is no MPI environment thus all blocks are loaded by process 0
      Nb = Nb_tot
      call parallel%alloc(Nb = Nb) ! allocating blockmap
      do b=1,Nb
        parallel%blockmap(b) = b ! the blocks map is identity
      enddo
    else
      ! computing the local (of myrank) number of blocks
      Nb = 0
      do b=1,Nb_tot
        if (parallel%procmap(b)==parallel%myrank) Nb = Nb + 1
      enddo
      call parallel%alloc(Nb = Nb) ! allocating blockmap
      bb = 0
      do b=1,Nb_tot
        if (parallel%procmap(b)==parallel%myrank) then
          bb = bb + 1
          parallel%blockmap(bb) = b
        endif
      enddo
    endif
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_procmap

  !> @brief Procedure for saving the processes/blocks map and local/global blocks map.
  subroutine save_procmap(file_d,parallel)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Procmap), intent(INOUT):: file_d   !< File data.
  type(Type_Parallel),      intent(IN)::    parallel !< Parallel data.
  integer(I4P)::                            b        !< Block counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(Nb_tot=>size(parallel%procmap),procmap=>parallel%procmap,unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg)
    call file_d%open(ascii=.true.,replace=.true.,action='WRITE') ; if (file_d%iostat/=0) return
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<?xml version="1.0"?>'
    do b=1,Nb_tot
      write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<block_process b="'//trim(str(n=b))//'" p="'//trim(str(n=procmap(b)))//&
        '"/>'
    enddo
    call file_d%close()
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_procmap
  !> @}
endmodule Data_Type_File_Procmap
