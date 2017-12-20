!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_GNUDerivedType Data_Type_File_GNU
!> Module definition of Gnuplot post-processed file, Type_File_GNU
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_GNUPrivateProcedure Data_Type_File_GNU
!> Module definition of Gnuplot post-processed file, Type_File_GNU
!> @}

!> @brief Module Data_Type_File_GNU contains the definition of Type_File_GNU, that is the post-processed output file in Gnuplot
!> format.
module Data_Type_File_GNU
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                               ! Integers and reals precision definition.
USE Data_Type_File_Base,   only: Type_File_Base,err_gnu_binary ! Definition of Type_File_Base.
USE Data_Type_Global,      only: Type_Global                   ! Definition of Type_Global.
USE Data_Type_PostProcess, only: Type_PostProcess              ! Definition of Type_PostProcess.
USE Data_Type_Primitive,   only: Type_Primitive                ! Definition of Type_Primitive.
USE Data_Type_SBlock,      only: Type_SBlock                   ! Definition of Type_SBlock.
USE Lib_IO_Misc,           only: Get_Unit,set_extension        ! Procedures for IO and strings operations.
USE Lib_Math,              only: digit                         ! Procedure for computing the significant digits of a number.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_GNU.
!> @ingroup Data_Type_File_GNUDerivedType
type, public, extends(Type_File_Base):: Type_File_GNU
  type(Type_PostProcess):: pp !< Post-processing options.
  contains
    procedure:: save_block => save_block_gnu ! Procedure for saving one block.
    procedure:: save       => save_file_gnu  ! Procedure for saving Gnuplot file.
endtype Type_File_GNU
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_GNUPrivateProcedure
  !> @{
  !> @brief Procedure for saving one block.
  subroutine save_block_gnu(file_d,filename,global,ID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_GNU), intent(INOUT):: file_d                  !< File data.
  character(*),         intent(IN)::    filename                !< File name.
  type(Type_Global),    intent(IN)::    global                  !< Global-level data.
  integer(I8P),         intent(IN)::    ID                      !< ID-key of block.
  type(Type_SBlock), pointer::          block                   !< Pointer for scanning global%block tree.
  type(Type_Primitive), allocatable::   P(:,:,:)                !< Prim variables intepolated at nodes.
  real(R8P),            allocatable::   r(:,:,:,:)              !< Array containing species densities.
  integer(I4P)::                        ni1,ni2,nj1,nj2,nk1,nk2 !< Bounds of of node-centered data.
  integer(I4P)::                        ci1,ci2,cj1,cj2,ck1,ck2 !< Bounds of of cell-centered data.
  integer(I4P)::                        nnode,ncell             !< Number of nodes and cells.
  integer(I4P)::                        i,j,k,s                 !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg,pp=>file_d%pp)
    associate(n=>global%time_step%n,t=>global%time_step%t)
      block => global%block%dat(ID=ID)
      associate(Ni=>block%dims%Ni,Nj=>block%dims%Nj,Nk=>block%dims%Nk,Ns=>block%dims%Ns)
        allocate(r(1:Ns,0:Ni,0:Nj,0:Nk))
        if (pp%node) call block%interpolate_primitive(primN=P) ! tri-linear interpolation of cell-centered values at nodes
        ! initialize the zone dimensions
        call pp%compute_dimensions(block=block,                                    &
                                   ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,&
                                   ci1=ci1,ci2=ci2,cj1=cj1,cj2=cj2,ck1=ck1,ck2=ck2)
        nnode = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1)
        ncell = (ci2-ci1+1)*(cj2-cj1+1)*(ck2-ck1+1)
        ! storing species densities into an array for avoiding problems with nonzero rank pointer
        if (pp%node) then
          do k=0,Nk
            do j=0,Nj
              do i=0,Ni
                do s=1,Ns
                  r(s,i,j,k) = P(i,j,k)%r(s)
                enddo
              enddo
            enddo
          enddo
        else
          do k=1,Nk
            do j=1,Nj
              do i=1,Ni
                do s=1,Ns
                  r(s,i,j,k) = block%C(i,j,k)%P%r(s)
                enddo
              enddo
            enddo
          enddo
        endif
        ! initializing file
        if (pp%binary) then
          call file_d%raise_error(errtype=err_gnu_binary)
        else
          open(unit=Get_Unit(unit),file=trim(filename))
          if (pp%node) then
            if (pp%meshonly) then
              do k=nk1,nk2
                do j=nj1,nj2
                  do i=ni1,ni2
                    write(unit,'(3('//FR8P//',1X))',iostat=iostat) block%node(i,j,k)%x,block%node(i,j,k)%y,block%node(i,j,k)%z
                  enddo
                  write(unit,*,iostat=iostat)
                enddo
                write(unit,*,iostat=iostat)
              enddo
            else
              do k=nk1,nk2
                do j=nj1,nj2
                  do i=ni1,ni2
                    write(unit,'(9('//FR8P//',1X))',iostat=iostat) block%node(i,j,k)%x, &
                                                                block%node(i,j,k)%y, &
                                                                block%node(i,j,k)%z, &
                                                                P(i,j,k)%d,          &
                                                                P(i,j,k)%v%x,        &
                                                                P(i,j,k)%v%y,        &
                                                                P(i,j,k)%v%z,        &
                                                                P(i,j,k)%p,          &
                                                                P(i,j,k)%g
                  enddo
                  write(unit,*,iostat=iostat)
                enddo
                write(unit,*,iostat=iostat)
              enddo
            endif
          else
            if (pp%meshonly) then
              do k=ck1,ck2
                do j=cj1,cj2
                  do i=ci1,ci2
                    write(unit,'(3('//FR8P//',1X))',iostat=iostat) block%C(i,j,k)%cent%x,block%C(i,j,k)%cent%y,block%C(i,j,k)%cent%z
                  enddo
                  write(unit,*,iostat=iostat)
                enddo
                write(unit,*,iostat=iostat)
              enddo
            else
              do k=ck1,ck2
                do j=cj1,cj2
                  do i=ci1,ci2
                    write(unit,'(9('//FR8P//',1X))',iostat=iostat) block%C(i,j,k)%cent%x, &
                                                                block%C(i,j,k)%cent%y, &
                                                                block%C(i,j,k)%cent%z, &
                                                                block%C(i,j,k)%P%d,    &
                                                                block%C(i,j,k)%P%v%x,  &
                                                                block%C(i,j,k)%P%v%y,  &
                                                                block%C(i,j,k)%P%v%z,  &
                                                                block%C(i,j,k)%P%p,    &
                                                                block%C(i,j,k)%P%g
                  enddo
                  write(unit,*,iostat=iostat)
                enddo
                write(unit,*,iostat=iostat)
              enddo
            endif
          endif
          close(unit)
        endif
      endassociate
    endassociate
  endassociate
  call file_d%fallback
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_block_gnu

  !> @brief Procedure for saving Gnuplot file.
  subroutine save_file_gnu(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_GNU), intent(INOUT):: file_d !< File data.
  type(Type_Global),    intent(IN)::    global !< Global-level data.
  integer(I8P)::                        ID     !< Counter.
  integer(I8P)::                        Nb_tot !< Total number of blocks.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nb_tot = global%block%length()
  associate(iostat=>file_d%iostat,pp=>file_d%pp)
    ! writing data blocks
    do while(global%block%loopID(ID=ID))
      call file_d%save_block(filename=set_extension(filename=file_d%name//'-ID'//trim(strz(digit(Nb_tot),ID)),extension='gnu'),&
                             global=global,ID=ID)
    enddo
  endassociate
  call file_d%fallback
  return
  !---------------------------------------------------------------------------------------------------------------------------------

  endsubroutine save_file_gnu
  !> @}
endmodule Data_Type_File_GNU
