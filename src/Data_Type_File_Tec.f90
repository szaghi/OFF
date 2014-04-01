!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_TecDerivedType Data_Type_File_Tec
!> Module definition of Tecplot post-processed file, Type_File_Tec
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_TecPrivateProcedure Data_Type_File_Tec
!> Module definition of Tecplot post-processed file, Type_File_Tec
!> @}

!> @brief Module Data_Type_File_Tec contains the definition of Type_File_Tec, that is the post-processed output file in Tecplot
!> format.
module Data_Type_File_Tec
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                              ! Integers and reals precision definition.
USE Data_Type_File_Base,   only: Type_File_Base,err_not_tecio ! Definition of Type_File_Base.
USE Data_Type_Global,      only: Type_Global                  ! Definition of Type_Global.
USE Data_Type_PostProcess, only: Type_PostProcess             ! Definition of Type_PostProcess.
USE Data_Type_Primitive,   only: Type_Primitive               ! Definition of Type_Primitive.
USE Data_Type_SBlock,      only: Type_SBlock                  ! Definition of Type_SBlock.
USE Lib_IO_Misc,           only: set_extension                ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Tec.
!> @ingroup Data_Type_File_TecDerivedType
type, public, extends(Type_File_Base):: Type_File_Tec
  type(Type_PostProcess)::        pp                  !< Post-processing options.
  character(1)::                  tecendrec = char(0) !< End-character for binary-record end.
  character(len=:), allocatable:: tecvarname          !< Variables name for tecplot header file.
  integer(I4P),     allocatable:: tecvarloc(:)        !< Tecplot array of variables location.
  character(len=:), allocatable:: tecvarlocstr        !< Tecplot string of variables location.
  integer(I4P),     allocatable:: tecnull(:)          !< Tecplot null array.
  integer(I4P)::                  nvar                !< Number of variables saved.
  contains
    procedure:: free       => free_file_tec  ! Procedure for freeing dynamic memory.
    procedure:: init       => init_file_tec  ! Procedure for initializing file data.
    procedure:: save_block => save_block_tec ! Procedure for saving one block.
    procedure:: save       => save_file_tec  ! Procedure for saving Tecplot file.
    final::     finalize                     ! Procedure for freeing dynamic memory when finalizing.
endtype Type_File_Tec
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_TecPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_file_tec(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Tec), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free_base
  if (allocated(file_d%tecvarname   )) deallocate(file_d%tecvarname   )
  if (allocated(file_d%tecvarloc    )) deallocate(file_d%tecvarloc    )
  if (allocated(file_d%tecvarlocstr )) deallocate(file_d%tecvarlocstr )
  if (allocated(file_d%tecnull      )) deallocate(file_d%tecnull      )
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_file_tec

  !> @brief Procedure for freeing dynamic memory when finalizng.
  elemental subroutine finalize(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_File_Tec), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for initializing file data.
  subroutine init_file_tec(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Tec), intent(INOUT):: file_d !< File data.
  type(Type_Global),    intent(IN)::    global !< Global-level data.
  integer(I4P)::                        s      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(file_d%tecvarloc(1:3+global%species0%Np))
  allocate(file_d%tecnull(  1:3+global%species0%Np))
  associate(pp=>file_d%pp,nvar=>file_d%nvar,Ns=>global%species0%Ns,Np=>global%species0%Np)
    if (pp%binary) then
      file_d%tecvarname = 'x y z'
      if (.not.pp%meshonly) then
        do s=1,Ns
          file_d%tecvarname = trim(file_d%tecvarname)//' r_'//trim(str(.true.,s))
        enddo
        file_d%tecvarname = trim(file_d%tecvarname)//' r u v w p g'
      endif
      if (pp%node) then
        file_d%tecvarloc = 1
      else
        file_d%tecvarloc(1:3) = 1 ; file_d%tecvarloc(4:3+Np)= 0
      endif
      file_d%tecnull = 0
      if (.not.pp%meshonly) then
        nvar = 3 + Np
      else
        nvar = 3
      endif
    else
      file_d%tecvarname = ' VARIABLES ="x" "y" "z"'
      if (.not.pp%meshonly) then
        do s=1,Ns
          file_d%tecvarname = trim(file_d%tecvarname)//' "r_'//trim(str(.true.,s))//'"'
        enddo
        file_d%tecvarname = trim(file_d%tecvarname)//' "r" "u" "v" "w" "p" "g"'
      endif
      if (.not.pp%meshonly) then
        if (pp%node) then
          file_d%tecvarlocstr = ', VARLOCATION=([1-'//trim(str(.true.,Np))//']=NODAL)'
        else
          file_d%tecvarlocstr = ', VARLOCATION=([1-3]=NODAL,[4-'//trim(str(.true.,Np))//']=CELLCENTERED)'
        endif
      else
        file_d%tecvarlocstr = ', VARLOCATION=([1-3]=NODAL)'
      endif
    endif
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_file_tec

  !> @brief Procedure for saving one block.
  subroutine save_block_tec(file_d,global,b,l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Tec), intent(INOUT):: file_d                  !< File data.
  type(Type_Global),    intent(IN)::    global                  !< Global-level data.
  integer(I4P),         intent(IN)::    b                       !< Block number.
  integer(I4P),         intent(IN)::    l                       !< Grid level.
  type(Type_SBlock)::                   blockmir                !< Mirrored block.
  character(len=:), allocatable::       zname                   !< Zone name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call get_zname(prefix='')
  call save_block(zname=zname,blk=global%block(b,l),t=global%time_step%t)
  if (file_d%pp%mirrorX) then
    call get_zname(prefix='mirX-')
    blockmir = global%block(b,l)%mirror(mirrorX=.true.)
    call save_block(zname=zname,blk=blockmir,t=global%time_step%t)
    call blockmir%free
  endif
  if (file_d%pp%mirrorY) then
    call get_zname(prefix='mirY-')
    blockmir = global%block(b,l)%mirror(mirrorY=.true.)
    call save_block(zname=zname,blk=blockmir,t=global%time_step%t)
    call blockmir%free
  endif
  if (file_d%pp%mirrorZ) then
    call get_zname(prefix='mirZ-')
    blockmir = global%block(b,l)%mirror(mirrorZ=.true.)
    call save_block(zname=zname,blk=blockmir,t=global%time_step%t)
    call blockmir%free
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    subroutine get_zname(prefix)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    character(*), intent(IN):: prefix !< Zone name prefix.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (file_d%pp%binary) then
      zname = trim(adjustl(prefix))//'n_'//trim(strz(nz_pad=10,n=global%time_step%n))//'-b_'//&
              trim(strz(nz_pad=global%mesh_dims%Nb,n=b))//file_d%tecendrec
    else
      zname = ' ZONE  T = "'//trim(adjustl(prefix))//'n_'//trim(strz(nz_pad=10,n=global%time_step%n))//'-b_'//&
              trim(strz(nz_pad=global%mesh_dims%Nb,n=b))//'"'
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine get_zname

    subroutine save_block(zname,blk,t)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    character(*),      intent(IN)::     zname                   !< Zone name.
    type(Type_SBlock), intent(IN)::     blk                     !< Block-level data.
    real(R8P),         intent(IN)::     t                       !< Time.
    type(Type_Primitive), allocatable:: P(:,:,:)                !< Prim variables intepolated at nodes.
    real(R8P),            allocatable:: r(:,:,:,:)              !< Array containing species densities.
    integer(I4P)::                      ni1,ni2,nj1,nj2,nk1,nk2 !< Bounds of of node-centered data.
    integer(I4P)::                      ci1,ci2,cj1,cj2,ck1,ck2 !< Bounds of of cell-centered data.
    integer(I4P)::                      nnode,ncell             !< Number of nodes and cells.
    integer(I4P)::                      i,j,k,s                 !< Counters.
    character(len=:), allocatable::     teczoneheader           !< Tecplot string of zone header.
#ifdef TECIO
    integer(I4P), external::            teczne112,tecdat112     !< Tecplot 'tecio' external functions.
#endif
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg,&
              pp=>file_d%pp,                                              &
              tecvarloc=>file_d%tecvarloc,                                &
              tecnull=>file_d%tecnull,                                    &
              nvar=>file_d%nvar)
      associate(Ni=>blk%dims%Ni,Nj=>blk%dims%Nj,Nk=>blk%dims%Nk,Ns=>blk%dims%Ns)
        allocate(r(1:Ns,0:Ni,0:Nj,0:Nk))
        if (pp%node) call blk%interpolate_primitive(primN=P) ! tri-linear interpolation of cell-centered values at nodes
        ! initialize the zone dimensions
        call pp%compute_dimensions(block=blk,                                      &
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
                  r(s,i,j,k) = blk%C(i,j,k)%P%r(s)
                enddo
              enddo
            enddo
          enddo
        endif
        ! writing the block data
        if (pp%binary) then
#ifdef TECIO
          iostat=teczne112(trim(adjustl(zname)),                                                         &
                           0,                                                                            &
                           ni2-ni1+1,                                                                    &
                           nj2-nj1+1,                                                                    &
                           nk2-nk1+1,                                                                    &
                           0,                                                                            &
                           0,                                                                            &
                           0,                                                                            &
                           t,                                                                            &
                           0,                                                                            &
                           0,                                                                            &
                           1,                                                                            &
                           0,                                                                            &
                           0,                                                                            &
                           0,                                                                            &
                           0,                                                                            &
                           0,                                                                            &
                           tecnull(1:nvar),                                                              &
                           tecvarloc(1:nvar),                                                            &
                           tecnull(1:nvar),                                                              &
                           0)
          iostat=tecdat112(nnode,blk%node(ni1:ni2,nj1:nj2,nk1:nk2)%x,1)
          iostat=tecdat112(nnode,blk%node(ni1:ni2,nj1:nj2,nk1:nk2)%y,1)
          iostat=tecdat112(nnode,blk%node(ni1:ni2,nj1:nj2,nk1:nk2)%z,1)
          if (.not.pp%meshonly) then
            if (pp%node) then
              do s=1,Ns
                iostat=tecdat112(ncell,r(s,ni1:ni2,nj1:nj2,nk1:nk2),1)
              enddo
              iostat=tecdat112(nnode,P(ni1:ni2,nj1:nj2,nk1:nk2)%d  ,1)
              iostat=tecdat112(nnode,P(ni1:ni2,nj1:nj2,nk1:nk2)%v%x,1)
              iostat=tecdat112(nnode,P(ni1:ni2,nj1:nj2,nk1:nk2)%v%y,1)
              iostat=tecdat112(nnode,P(ni1:ni2,nj1:nj2,nk1:nk2)%v%z,1)
              iostat=tecdat112(nnode,P(ni1:ni2,nj1:nj2,nk1:nk2)%p  ,1)
              iostat=tecdat112(nnode,P(ni1:ni2,nj1:nj2,nk1:nk2)%g  ,1)
            else
              do s=1,Ns
                file_d%iostat=tecdat112(ncell,r(s,ci1:ci2,cj1:cj2,ck1:ck2),1)
              enddo
              iostat=tecdat112(ncell,blk%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%d  ,1)
              iostat=tecdat112(ncell,blk%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%x,1)
              iostat=tecdat112(ncell,blk%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%y,1)
              iostat=tecdat112(ncell,blk%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%z,1)
              iostat=tecdat112(ncell,blk%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%p  ,1)
              iostat=tecdat112(ncell,blk%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%g  ,1)
            endif
          endif
#else
          call file_d%raise_error(errtype=err_not_tecio) ; return
#endif
        else
          teczoneheader = ' ZONE  T = "'//trim(adjustl(zname))//'"'//    &
                          ', I='//trim(str(no_sign=.true.,n=ni2-ni1+1))//&
                          ', J='//trim(str(no_sign=.true.,n=nj2-nj1+1))//&
                          ', K='//trim(str(no_sign=.true.,n=nk2-nk1+1))//&
                          ', DATAPACKING=BLOCK'//file_d%tecvarlocstr
          write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)teczoneheader
          write(unit=unit,fmt=FR8P, iostat=iostat,iomsg=iomsg)(((blk%node(i,j,k)%x,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(unit=unit,fmt=FR8P, iostat=iostat,iomsg=iomsg)(((blk%node(i,j,k)%y,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(unit=unit,fmt=FR8P, iostat=iostat,iomsg=iomsg)(((blk%node(i,j,k)%z,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          if (.not.pp%meshonly) then
            if (pp%node) then
              do s=1,Ns
                write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((P(i,j,k)%r(s),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
              enddo
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((P(i,j,k)%d  ,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((P(i,j,k)%v%x,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((P(i,j,k)%v%y,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((P(i,j,k)%v%z,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((P(i,j,k)%p  ,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((P(i,j,k)%g  ,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
            else
              do s=1,Ns
                write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((blk%C(i,j,k)%P%r(s),i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
              enddo
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((blk%C(i,j,k)%P%d  ,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((blk%C(i,j,k)%P%v%x,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((blk%C(i,j,k)%P%v%y,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((blk%C(i,j,k)%P%v%z,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((blk%C(i,j,k)%P%p  ,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
              write(unit=unit,fmt=FR8P,iostat=iostat,iomsg=iomsg)(((blk%C(i,j,k)%P%g  ,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
            endif
          endif
        endif
        deallocate(r)
      endassociate
    endassociate
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine save_block
  endsubroutine save_block_tec

  !> @brief Procedure for saving Tecplot file.
  subroutine save_file_tec(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Tec), intent(INOUT):: file_d                           !< File data.
  type(Type_Global),    intent(IN)::    global                           !< Global-level data.
  integer(I4P)::                        b,l                              !< Counters.
#ifdef TECIO
  integer(I4P), external::              tecini112,tecauxstr112,tecend112 !< Tecplot 'tecio' external functions.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%init(global=global)
  associate(pp=>file_d%pp,tecendrec=>file_d%tecendrec,tecvarname=>file_d%tecvarname)
    file_d%name = set_extension(filename=file_d%name,extension='plt')
    if (pp%binary) then
#ifdef TECIO
      file_d%iostat = tecini112(tecendrec,trim(tecvarname)//tecendrec,trim(file_d%name)//tecendrec,'.'//tecendrec,0,0,1)
      file_d%iostat = tecauxstr112("Time"//tecendrec,trim(str(n=global%time_step%t))//tecendrec)
#else
      call file_d%raise_error(errtype=err_not_tecio) ; return
#endif
    else
      call file_d%open(ascii=.true.,replace=.true.,action='WRITE') ; if (file_d%iostat/=0) return
      write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat)tecvarname
    endif
    ! writing data blocks
    do l=1,global%mesh_dims%Nl ; do b=1,global%mesh_dims%Nb
      call file_d%save_block(global=global,b=b,l=l)
    enddo ; enddo
    ! finalizing tecplot file
    if (pp%binary) then
#ifdef TECIO
      file_d%iostat = tecend112()
#else
      call file_d%raise_error(errtype=err_not_tecio) ; return
#endif
    else
      call file_d%close()
    endif
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------

  endsubroutine save_file_tec
  !> @}
endmodule Data_Type_File_Tec
