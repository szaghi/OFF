!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_VTKDerivedType Data_Type_File_VTK
!> Module definition of VTK post-processed file, Type_File_VTK
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_VTKPrivateProcedure Data_Type_File_VTK
!> Module definition of VTK post-processed file, Type_File_VTK
!> @}

!> @brief Module Data_Type_File_VTK contains the definition of Type_File_VTK, that is the post-processed output file in VTK format.
module Data_Type_File_VTK
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                  ! Integers and reals precision definition.
USE Data_Type_File_Base,   only: Type_File_Base   ! Definition of Type_File_Base.
USE Data_Type_Global,      only: Type_Global      ! Definition of Type_Global.
USE Data_Type_PostProcess, only: Type_PostProcess ! Definition of Type_PostProcess.
USE Data_Type_Primitive,   only: Type_Primitive   ! Definition of Type_Primitive.
USE Lib_IO_Misc,           only: set_extension    ! Procedures for imposing extension to a file name.
USE Lib_Math,              only: digit            ! Procedure for computing the significant digits of a number.
USE Lib_VTK_IO                                    ! Library for IO VTK files.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_VTK.
!> @ingroup Data_Type_File_VTKDerivedType
type, public, extends(Type_File_Base):: Type_File_VTK
  type(Type_PostProcess):: pp !< Post-processing options.
  contains
    procedure:: save_block => save_block_vtk ! Procedure for saving one block.
    procedure:: save       => save_file_vtk  ! Procedure for saving VTK file.
endtype Type_File_VTK
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_VTKPrivateProcedure
  !> @{
  !> @brief Procedure for saving one block.
  subroutine save_block_vtk(file_d,filename,global,b,l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_VTK), intent(INOUT):: file_d                  !< File data.
  character(*),         intent(IN)::    filename                !< File name.
  type(Type_Global),    intent(IN)::    global                  !< Global-level data.
  integer(I4P),         intent(IN)::    b                       !< Block number.
  integer(I4P),         intent(IN)::    l                       !< Grid level.
  type(Type_Primitive), allocatable::   P(:,:,:)                !< Prim variables intepolated at nodes.
  real(R8P),            allocatable::   r(:,:,:,:)              !< Array containing species densities.
  integer(I4P)::                        ni1,ni2,nj1,nj2,nk1,nk2 !< Bounds of of node-centered data.
  integer(I4P)::                        ci1,ci2,cj1,cj2,ck1,ck2 !< Bounds of of cell-centered data.
  integer(I4P)::                        nnode,ncell             !< Number of nodes and cells.
  integer(I4P)::                        i,j,k,s                 !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg,pp=>file_d%pp)
    associate(Nb=>global%mesh_dims%Nb,n=>global%time_step%n,t=>global%time_step%t,block=>global%block(b,l))
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
          iostat = VTK_INI_XML(output_format = 'binary',         &
                               filename      = trim(filename),   &
                               mesh_topology = 'StructuredGrid', &
                               nx1 = ni1, nx2 = ni2,             &
                               ny1 = nj1, ny2 = nj2,             &
                               nz1 = nk1, nz2 = nk2)
        else
          iostat = VTK_INI_XML(output_format = 'ascii',          &
                               filename      = trim(filename),   &
                               mesh_topology = 'StructuredGrid', &
                               nx1 = ni1, nx2 = ni2,             &
                               ny1 = nj1, ny2 = nj2,             &
                               nz1 = nk1, nz2 = nk2)
        endif
        ! saving auxiliary data (time and time step)
        iostat = VTK_FLD_XML(fld_action='open')
        iostat = VTK_FLD_XML(fld=global%time_step%t,fname='TIME')
        iostat = VTK_FLD_XML(fld=global%time_step%n,fname='CYCLE')
        iostat = VTK_FLD_XML(fld_action='close')
        ! saving the geometry
        iostat = VTK_GEO_XML(nx1 = ni1, nx2 = ni2,                    &
                             ny1 = nj1, ny2 = nj2,                    &
                             nz1 = nk1, nz2 = nk2,                    &
                             NN = nnode,                              &
                             X=block%node(ni1:ni2,nj1:nj2,nk1:nk2)%x, &
                             Y=block%node(ni1:ni2,nj1:nj2,nk1:nk2)%y, &
                             Z=block%node(ni1:ni2,nj1:nj2,nk1:nk2)%z)
        if (.not.pp%meshonly) then
          ! saving dependent variables
          if (pp%node) then
            iostat = VTK_DAT_XML(var_location = 'node', var_block_action = 'open')
            do s=1,block%dims%Ns
              iostat=VTK_VAR_XML(NC_NN=nnode,varname='r('//trim(str(.true.,s))//')',var=r(s,ni1:ni2,nj1:nj2,nk1:nk2))
            enddo
            iostat=VTK_VAR_XML(NC_NN=nnode,varname='r',var=P(ni1:ni2,nj1:nj2,nk1:nk2)%d  )
            iostat=VTK_VAR_XML(NC_NN=nnode,varname='u',var=P(ni1:ni2,nj1:nj2,nk1:nk2)%v%x)
            iostat=VTK_VAR_XML(NC_NN=nnode,varname='v',var=P(ni1:ni2,nj1:nj2,nk1:nk2)%v%y)
            iostat=VTK_VAR_XML(NC_NN=nnode,varname='w',var=P(ni1:ni2,nj1:nj2,nk1:nk2)%v%z)
            iostat=VTK_VAR_XML(NC_NN=nnode,varname='p',var=P(ni1:ni2,nj1:nj2,nk1:nk2)%p  )
            iostat=VTK_VAR_XML(NC_NN=nnode,varname='g',var=P(ni1:ni2,nj1:nj2,nk1:nk2)%g  )
            iostat=VTK_DAT_XML(var_location ='node',var_block_action = 'close')
          else
            iostat = VTK_DAT_XML(var_location = 'cell', var_block_action = 'open')
            do s=1,block%dims%Ns
              iostat=VTK_VAR_XML(NC_NN=ncell,varname='r('//trim(str(.true.,s))//')',var=r(s,ci1:ci2,cj1:cj2,ck1:ck2))
            enddo
            iostat=VTK_VAR_XML(NC_NN=ncell,varname='r',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%d  )
            iostat=VTK_VAR_XML(NC_NN=ncell,varname='u',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%x)
            iostat=VTK_VAR_XML(NC_NN=ncell,varname='v',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%y)
            iostat=VTK_VAR_XML(NC_NN=ncell,varname='w',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%z)
            iostat=VTK_VAR_XML(NC_NN=ncell,varname='p',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%p  )
            iostat=VTK_VAR_XML(NC_NN=ncell,varname='g',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%g  )
            iostat=VTK_DAT_XML(var_location ='cell',var_block_action = 'close')
          endif
        endif
        ! closing VTK file
        iostat = VTK_GEO_XML()
        iostat = VTK_END_XML()
      endassociate
    endassociate
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_block_vtk

  !> @brief Procedure for saving VTK file.
  subroutine save_file_vtk(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_VTK), intent(INOUT):: file_d !< File data.
  type(Type_Global),    intent(IN)::    global !< Global-level data.
  integer(I4P)::                        b,l    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(iostat=>file_d%iostat,pp=>file_d%pp,Nb_tot=>global%mesh_dims%Nb_tot,Nb=>global%mesh_dims%Nb,Nl=>global%mesh_dims%Nl)
    ! writing data blocks
    do l=1,Nl ; do b=1,Nb
    call file_d%save_block(filename=set_extension(filename=file_d%name//                                                    &
                                                           '-b'//trim(strz(digit(Nb_tot),b))//'-l'//trim(strz(digit(Nl),l)),&
                                                  extension='vts'),global=global,b=b,l=l)
    enddo ; enddo
    ! writing the VTM wrapper for multiblock dataset
    iostat = VTM_INI_XML(filename=set_extension(filename=file_d%name,extension='vtm'))
    iostat = VTM_BLK_XML(block_action='open')
    iostat = VTM_WRF_XML(flist=[((global%OS%basename(&
                                  set_extension(filename=file_d%name//                                                    &
                                                         '-b'//trim(strz(digit(Nb_tot),b))//'-l'//trim(strz(digit(Nl),l)),&
                                                 extension='vts')),b=1,Nb),l=1,Nl)])
    iostat = VTM_BLK_XML(block_action='close')
    iostat = VTM_END_XML()
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------

  endsubroutine save_file_vtk
  !> @}
endmodule Data_Type_File_VTK
