!> @ingroup Library
!> @{
!> @defgroup Lib_PostProcessingLibrary Lib_PostProcessing
!> @}

!> @ingroup DerivedType
!> @{
!> @defgroup Lib_PostProcessingDerivedType Lib_PostProcessing
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Lib_PostProcessingGlobalVarPar Lib_PostProcessing
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_PostProcessingPublicProcedure Lib_PostProcessing
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_PostProcessingPrivateProcedure Lib_PostProcessing
!> @}

!> This module contains the definition of procedures and variables useful for post-process the OFF data.
!> @ingroup Lib_PostProcessingLibrary
module Lib_PostProcessing
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision        ! Integers and reals precision definition.
USE Data_Type_Global    ! Definition of Type_Global.
USE Data_Type_Primitive ! Definition of Type_Primitive.
USE Data_Type_SBlock    ! Definition of Type_SBlock.
USE Lib_IO_Misc         ! Procedures for IO and strings operations.
USE Lib_VTK_IO          ! Library of functions and subroutines for I/O VTK files.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: pp_format
public:: tec_output
public:: vtk_output
public:: gnu_output
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing the post-processing options.
!> @ingroup Lib_PostProcessingDerivedType
type, public:: Type_PP_Format
  logical:: binary = .true.  !< Binary or ascii post-process file.
  logical:: node   = .true.  !< Node or cell data location.
  logical:: bc     = .false. !< Saving or not boundary conditions cells.
  logical:: tec    = .true.  !< Tecplot file.
  logical:: vtk    = .false. !< VTK file.
  logical:: gnu    = .false. !< Gnuplot file.
endtype Type_PP_Format
!> @ingroup Lib_PostProcessingGlobalVarPar
!> @{
type(Type_PP_Format):: pp_format !< Post-processing format options (see \ref Lib_PostProcessing::Type_PP_Format "definition").
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_PostProcessingPrivateProcedure
  !> @{
  subroutine interpolate_primitive(block,P)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for interpolating primitive variables at nodes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_SBlock),    intent(INOUT):: block                               ! Block-level data.
  type(Type_Primitive), intent(INOUT):: P(0:block%Ni,0:block%Nj,0:block%Nk) ! Primitive variables intepolated at nodes.
  real(R_P)::                           mf                                  ! Mean factor.
  integer(I_P)::                        i,j,k                               ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call P%init(Ns=block%global%Ns)
#if !defined NULi && !defined NULj && !defined NULk
  ! 3D data
  mf = 0.125_R_P
#elif defined NULi
  block%C(0         ,:,:)%P = 0._R_P
  block%C(block%Ni+1,:,:)%P = 0._R_P
#if !defined NULj && !defined NULk
  ! 2D data
  mf = 0.25_R_P
#elif defined NULj
  ! 1D data
  mf = 0.5_R_P
  block%C(:,0         ,:)%P = 0._R_P
  block%C(:,block%Nj+1,:)%P = 0._R_P
#elif defined NULk
  ! 1D data
  mf = 0.5_R_P
  block%C(:,:,0         )%P = 0._R_P
  block%C(:,:,block%Nk+1)%P = 0._R_P
#endif
#elif defined NULj
  block%C(:,0         ,:)%P = 0._R_P
  block%C(:,block%Nj+1,:)%P = 0._R_P
#if !defined NULi && !defined NULk
  ! 2D data
  mf = 0.25_R_P
#elif defined NULi
  ! 1D data
  mf = 0.5_R_P
  block%C(0         ,:,:)%P = 0._R_P
  block%C(block%Ni+1,:,:)%P = 0._R_P
#elif defined NULk
  ! 1D data
  mf = 0.5_R_P
  block%C(:,:,0         )%P = 0._R_P
  block%C(:,:,block%Nk+1)%P = 0._R_P
#endif
#elif defined NULk
  block%C(:,:,0         )%P = 0._R_P
  block%C(:,:,block%Nk+1)%P = 0._R_P
#if !defined NULi && !defined NULj
  ! 2D data
  mf = 0.25_R_P
#elif defined NULi
  ! 1D data
  mf = 0.5_R_P
  block%C(0         ,:,:)%P = 0._R_P
  block%C(block%Ni+1,:,:)%P = 0._R_P
#elif defined NULj
  ! 1D data
  mf = 0.5_R_P
  block%C(:,0         ,:)%P = 0._R_P
  block%C(:,block%Nj+1,:)%P = 0._R_P
#endif
#endif
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(block,P,mf)
  !$OMP DO
  do k=0,block%Nk
    do j=0,block%Nj
      do i=0,block%Ni
          P(i,j,k) = block%C(i+1,j+1,k+1)%P + block%C(i,j+1,k+1)%P &
                   + block%C(i+1,j  ,k+1)%P + block%C(i,j,  k+1)%P &
                   + block%C(i+1,j+1,k  )%P + block%C(i,j+1,k  )%P &
                   + block%C(i+1,j  ,k  )%P + block%C(i,j  ,k  )%P
          P(i,j,k) = mf*P(i,j,k)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine interpolate_primitive

  subroutine compute_dimensions(block,ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing the dimensions of the domain to be post-processed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_SBlock), intent(IN)::  block                   ! Block-level data.
  integer(I_P),      intent(OUT):: ni1,ni2,nj1,nj2,nk1,nk2 ! Bounds of dimensions of node-centered data.
  integer(I_P),      intent(OUT):: ci1,ci2,cj1,cj2,ck1,ck2 ! Bounds of dimensions of cell-centered data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (pp_format%node) then
    ni1 = 0 ; ni2 = block%Ni
    nj1 = 0 ; nj2 = block%Nj
    nk1 = 0 ; nk2 = block%Nk

    ci1 = 0 ; ci2 = block%Ni
    cj1 = 0 ; cj2 = block%Nj
    ck1 = 0 ; ck2 = block%Nk
  else
    if (pp_format%bc) then
      ni1 = 0 - block%gc(1) ; ni2 = block%Ni + block%gc(2)
      nj1 = 0 - block%gc(3) ; nj2 = block%Nj + block%gc(4)
      nk1 = 0 - block%gc(5) ; nk2 = block%Nk + block%gc(6)

      ci1 = 1 - block%gc(1) ; ci2 = block%Ni + block%gc(2)
      cj1 = 1 - block%gc(3) ; cj2 = block%Nj + block%gc(4)
      ck1 = 1 - block%gc(5) ; ck2 = block%Nk + block%gc(6)
    else
      ni1 = 0               ; ni2 = block%Ni
      nj1 = 0               ; nj2 = block%Nj
      nk1 = 0               ; nk2 = block%Nk

      ci1 = 1               ; ci2 = block%Ni
      cj1 = 1               ; cj2 = block%Nj
      ck1 = 1               ; ck2 = block%Nk
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_dimensions
  !> @}

  !> @ingroup Lib_PostProcessingPublicProcedure
  !> @{
  !> Function for writing OFF block data to Tecplot file.
  function tec_output(meshonly,block,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,           intent(IN)::    meshonly            !< Flag for post-process only mesh.
  type(Type_SBlock), intent(INOUT):: block(1:)           !< Block-level data.
  character(*),      intent(IN)::    filename            !< File name of the output file.
  integer(I_P)::                     err                 !< Error trapping flag: 0 no errors, >0 error occurs.
  type(Type_Global), pointer::       global              !< Global-level data.
#ifdef TECIO
  integer(I_P), external::           tecini112,    &     ! |
                                     tecauxstr112, &     ! |
                                     teczne112,    &     ! | Tecplot external functions.
                                     tecdat112,    &     ! |
                                     tecend112           ! |
#endif
  character(1), parameter::          tecendrec = char(0) !< End-character for binary-record end.
  character(500)::                   tecvarname          !< Variables name for tecplot header file.
  character(500)::                   teczoneheader       !< Tecplot string of zone header.
  character(500)::                   tecvarform          !< Format for variables for tecplot file.
  integer(I_P), allocatable::        tecvarloc(:)        !< Tecplot array of variables location.
  character(500)::                   tecvarlocstr        !< Tecplot string of variables location.
  integer(I_P), allocatable::        tecnull(:)          !< Tecplot null array.
  integer(I_P)::                     tecunit             !< Free logic unit of tecplot file.
  integer(I_P)::                     nvar                !< Number of variables saved.
  integer(I_P)::                     b                   !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  global => block(1)%global
  ! allocating dynamic arrays
  allocate(tecvarloc(1:3+global%Np))
  allocate(tecnull(1:3+global%Np))
  ! initializing tecplot variables
  call tec_init()
  ! initializing tecplot file
  if (pp_format%binary) then
#ifdef TECIO
    if (filename(len_trim(filename)-4:len_trim(filename))/=".plt") then
      err = tecini112(tecendrec,trim(tecvarname)//tecendrec,trim(filename)//".plt"//tecendrec,'.'//tecendrec,0,0,1)
    else
      err = tecini112(tecendrec,trim(tecvarname)//tecendrec,trim(filename)//tecendrec,'.'//tecendrec,0,0,1)
    endif
    err = tecauxstr112("Time"//tecendrec,trim(str(n=global%t))//tecendrec)
#else
    write(stderr,'(A)') 'Error: your are trying to save binary tecplot file without compiling against the Tecplot library.'
    stop
#endif
  else
    if (filename(len_trim(filename)-4:len_trim(filename))/=".dat") then
      open(unit=Get_Unit(tecunit),file=trim(filename)//".dat")
    else
      open(unit=Get_Unit(tecunit),file=trim(filename))
    endif
    write(tecunit,'(A)',iostat=err)trim(tecvarname)
  endif
  ! writing data blocks
  do b=1,global%Nb
    err = tec_blk_data(b = b, global = global, block = block(b))
  enddo
  ! finalizing tecplot file
  if (pp_format%binary) then
#ifdef TECIO
    err = tecend112()
#else
    write(stderr,'(A)') 'Error: your are trying to save binary tecplot file without compiling against the Tecplot library.'
    stop
#endif
  else
    close(tecunit)
  endif
  ! deallocating dynamic arrays
  deallocate(tecvarloc)
  deallocate(tecnull)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    subroutine tec_init()
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Function for initializing Tecplot specific variables.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P):: s ! Counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (pp_format%binary) then
      ! header variables names
      tecvarname = 'x y z'
      if (.not.meshonly) then
        do s=1,global%Ns
          tecvarname = trim(tecvarname)//' r_'//trim(str(.true.,s))
        enddo
        tecvarname = trim(tecvarname)//' r u v w p g'
      endif
      ! variables location
      if (pp_format%node) then
        tecvarloc = 1
      else
        tecvarloc(1:3) = 1 ; tecvarloc(4:3+global%Np)= 0
      endif
      ! null array
      tecnull = 0
      if (.not.meshonly) then
        nvar = 3 + global%Np
      else
        nvar = 3
      endif
    else
      ! header variables names
      tecvarname = ' VARIABLES ="x" "y" "z"'
      if (.not.meshonly) then
        do s=1,global%Ns
          tecvarname = trim(tecvarname)//' "r_'//trim(str(.true.,s))//'"'
        enddo
        tecvarname = trim(tecvarname)//' "r" "u" "v" "w" "p" "g"'
      endif
      ! variables output format
      if (.not.meshonly) then
        write(tecvarform,'(A)')'('//trim(str(no_sign=.true.,n=global%Np))//'('//FR_P//',1X))'
      else
        write(tecvarform,'(A)')'('//trim(str(no_sign=.true.,n=3))//'('//FR_P//',1X))'
      endif
      ! variables location
      if (.not.meshonly) then
        if (pp_format%node) then
          tecvarlocstr = ', VARLOCATION=([1-'//trim(str(.true.,global%Np))//']=NODAL)'
        else
          tecvarlocstr = ', VARLOCATION=([1-3]=NODAL,[4-'//trim(str(.true.,global%Np))//']=CELLCENTERED)'
        endif
      else
        tecvarlocstr = ', VARLOCATION=([1-3]=NODAL)'
      endif
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine tec_init

    function tec_blk_data(b,global,block) result(err)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Function for writing block data.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),      intent(IN)::    b      ! Block number.
    type(Type_Global), intent(IN)::    global ! Global-level data.
    type(Type_SBlock), intent(INOUT):: block  ! Block-level data.
    integer(I_P)::                     err    ! Error trapping flag: 0 no errors, >0 error occurs.
    type(Type_Primitive):: P(            0:block%Ni,0:block%Nj,0:block%Nk) ! Primitive variables intepolated at nodes.
    real(R_P)::            r(1:global%Ns,0:block%Ni,0:block%Nj,0:block%Nk) ! Array containing species densities.
    integer(I_P)::         ni1,ni2,nj1,nj2,nk1,nk2                         ! Bounds of dimensions of node-centered data.
    integer(I_P)::         ci1,ci2,cj1,cj2,ck1,ck2                         ! Bounds of dimensions of cell-centered data.
    integer(I_P)::         nnode,ncell                                     ! Number of nodes and cells.
    integer(I_P)::         i,j,k,s                                         ! Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (pp_format%node) then
      ! tri-linear interpolation of cell-centered values at nodes
      call interpolate_primitive(block=block,P=P)
    endif
    ! initialize the zone dimensions
    call compute_dimensions(block=block,                                     &
                            ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2, &
                            ci1=ci1,ci2=ci2,cj1=cj1,cj2=cj2,ck1=ck1,ck2=ck2)
    nnode = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1)
    ncell = (ci2-ci1+1)*(cj2-cj1+1)*(ck2-ck1+1)
    ! storing species densities into an array for avoiding problems with nonzero rank pointer
    if (pp_format%node) then
      do k=0,block%Nk
        do j=0,block%Nj
          do i=0,block%Ni
            do s=1,global%Ns
              r(s,i,j,k) = P(i,j,k)%r(s)
            enddo
          enddo
        enddo
      enddo
    else
      do k=1,block%Nk
        do j=1,block%Nj
          do i=1,block%Ni
            do s=1,global%Ns
              r(s,i,j,k) = block%C(i,j,k)%P%r(s)
            enddo
          enddo
        enddo
      enddo
    endif
    ! writing the block data
    if (pp_format%binary) then
#ifdef TECIO
      err = teczne112('n_'//trim(strz(nz_pad=10,n=global%n))//'-b_'//trim(strz(nz_pad=4,n=b))//tecendrec, &
                      0,                                                                                  &
                      ni2-ni1+1,                                                                          &
                      nj2-nj1+1,                                                                          &
                      nk2-nk1+1,                                                                          &
                      0,                                                                                  &
                      0,                                                                                  &
                      0,                                                                                  &
                      global%t,                                                                           &
                      0,                                                                                  &
                      0,                                                                                  &
                      1,                                                                                  & !1=>block,0=>point
                      0,                                                                                  &
                      0,                                                                                  &
                      0,                                                                                  &
                      0,                                                                                  &
                      0,                                                                                  &
                      tecnull(1:nvar),                                                                    &
                      tecvarloc(1:nvar),                                                                  &
                      tecnull(1:nvar),                                                                    &
                      0)
#else
      write(stderr,'(A)') 'Error: your are trying to save binary tecplot file without compiling against the Tecplot library.'
      stop
#endif
      err=tec_dat(N=nnode,dat=block%node(ni1:ni2,nj1:nj2,nk1:nk2)%x)
      err=tec_dat(N=nnode,dat=block%node(ni1:ni2,nj1:nj2,nk1:nk2)%y)
      err=tec_dat(N=nnode,dat=block%node(ni1:ni2,nj1:nj2,nk1:nk2)%z)
      if (.not.meshonly) then
        if (pp_format%node) then
          do s=1,global%Ns
            err=tec_dat(N=ncell,dat=r(s,ni1:ni2,nj1:nj2,nk1:nk2))
          enddo
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%d  )
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%v%x)
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%v%y)
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%v%z)
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%p  )
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%g  )
        else
          do s=1,global%Ns
            err=tec_dat(N=ncell,dat=r(s,ci1:ci2,cj1:cj2,ck1:ck2))
          enddo
          err=tec_dat(N=ncell,dat=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%d  )
          err=tec_dat(N=ncell,dat=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%x)
          err=tec_dat(N=ncell,dat=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%y)
          err=tec_dat(N=ncell,dat=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%z)
          err=tec_dat(N=ncell,dat=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%p  )
          err=tec_dat(N=ncell,dat=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%g  )
        endif
      endif
    else
      ! tecplot zone header
      teczoneheader = ' ZONE  T = "n_'//trim(strz(nz_pad=10,n=global%n))// &
                                 '-b_'//trim(strz(nz_pad=4,n=b))//'"'//    &
                      ', I='//trim(str(no_sign=.true.,n=ni2-ni1+1))//      &
                      ', J='//trim(str(no_sign=.true.,n=nj2-nj1+1))//      &
                      ', K='//trim(str(no_sign=.true.,n=nk2-nk1+1))//      &
                      ', DATAPACKING=BLOCK'//adjustl(trim(tecvarlocstr))
      write(tecunit,'(A)',iostat=err)trim(teczoneheader)
      write(tecunit,FR_P,iostat=err)(((block%node(i,j,k)%x,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
      write(tecunit,FR_P,iostat=err)(((block%node(i,j,k)%y,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
      write(tecunit,FR_P,iostat=err)(((block%node(i,j,k)%z,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
      if (.not.meshonly) then
        if (pp_format%node) then
          do s=1,global%Ns
            write(tecunit,FR_P,iostat=err)(((P(i,j,k)%r(s),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          enddo
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%d  ,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%v%x,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%v%y,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%v%z,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%p  ,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%g  ,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
        else
          do s=1,global%Ns
            write(tecunit,FR_P,iostat=err)(((block%C(i,j,k)%P%r(s),i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          enddo
          write(tecunit,FR_P,iostat=err)(((block%C(i,j,k)%P%d  ,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          write(tecunit,FR_P,iostat=err)(((block%C(i,j,k)%P%v%x,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          write(tecunit,FR_P,iostat=err)(((block%C(i,j,k)%P%v%y,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          write(tecunit,FR_P,iostat=err)(((block%C(i,j,k)%P%v%z,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          write(tecunit,FR_P,iostat=err)(((block%C(i,j,k)%P%p  ,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          write(tecunit,FR_P,iostat=err)(((block%C(i,j,k)%P%g  ,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
        endif
      endif
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction tec_blk_data

    function tec_dat(N,dat) result(err)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Function interface for using "tecdat" function.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P), intent(IN):: N        ! Number of data to save.
    real(R_P),    intent(IN):: dat(1:N) ! Data to save.
    integer(I_P)::             err      ! Error trapping flag: 0 no errors, >0 error occurs.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
#ifdef TECIO
    err = tecdat112(N,dat,1)
#else
    write(stderr,'(A)') 'Error: your are trying to save binary tecplot file without compiling against the Tecplot library.'
    stop
#endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction tec_dat
  endfunction tec_output

  !> Function for writing OFF block data to VTK file.
  function vtk_output(meshonly,block,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,           intent(IN)::    meshonly  !< Flag for post-process only mesh.
  type(Type_SBlock), intent(INOUT):: block(1:) !< Block-level data.
  character(*),      intent(IN)::    filename  !< File name of the output file.
  integer(I_P)::                     err       !< Error trapping flag: 0 no errors, >0 error occurs.
  type(Type_Global), pointer::       global    !< Global-level data.
  integer(I_P)::                     b         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  global => block(1)%global
  ! writing data blocks
  do b=1,global%Nb
    err = vtk_blk_data( global = global, block = block(b), filename = trim(filename)//'_b'//trim(strz(4,b))//'.vts')
  enddo
  !! writing the VTM wrapper for multiblock dataset
  !err = VTM_INI_XML(filename=trim(filename)//'.vtm')
  !err = VTM_BLK_XML(block_action='open')
  !err = VTM_WRF_XML(vtk_xml_file_list=(/(trim(filename)//'_b'//trim(strz(4,b))//'.vts',b=1,global%Nb)/))
  !err = VTM_BLK_XML(block_action='close')
  !err = VTM_END_XML()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    function vtk_blk_data(global,block,filename) result(err)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Function for writing block data.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    type(Type_Global), intent(IN)::    global   ! Global-level data.
    type(Type_SBlock), intent(INOUT):: block    ! Block-level data.
    character(*),      intent(IN)::    filename ! File name of the output block file.
    integer(I_P)::                     err      ! Error trapping flag: 0 no errors, >0 error occurs.
    type(Type_Primitive):: P(            0:block%Ni,0:block%Nj,0:block%Nk) ! Primitive variables intepolated at nodes.
    real(R_P)::            r(1:global%Ns,0:block%Ni,0:block%Nj,0:block%Nk) ! Array containing species densities.
    integer(I_P)::         ni1,ni2,nj1,nj2,nk1,nk2                         ! Bounds of dimensions of node-centered data.
    integer(I_P)::         ci1,ci2,cj1,cj2,ck1,ck2                         ! Bounds of dimensions of cell-centered data.
    integer(I_P)::         nnode,ncell                                     ! Number of nodes and cells.
    integer(I_P)::         i,j,k,s                                         ! Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (pp_format%node) then
      ! tri-linear interpolation of cell-centered values at nodes
      call interpolate_primitive(block=block,P=P)
    endif
    ! initialize the block dimensions
    call compute_dimensions(block=block,                                     &
                            ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2, &
                            ci1=ci1,ci2=ci2,cj1=cj1,cj2=cj2,ck1=ck1,ck2=ck2)
    nnode = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1)
    ncell = (ci2-ci1+1)*(cj2-cj1+1)*(ck2-ck1+1)
    ! storing species densities into an array for avoiding problems with nonzero rank pointer
    if (pp_format%node) then
      do k=0,block%Nk
        do j=0,block%Nj
          do i=0,block%Ni
            do s=1,global%Ns
              r(s,i,j,k) = P(i,j,k)%r(s)
            enddo
          enddo
        enddo
      enddo
    else
      do k=1,block%Nk
        do j=1,block%Nj
          do i=1,block%Ni
            do s=1,global%Ns
              r(s,i,j,k) = block%C(i,j,k)%P%r(s)
            enddo
          enddo
        enddo
      enddo
    endif
    ! initializing VTK file
    if (pp_format%binary) then
      err = VTK_INI_XML(output_format = 'binary',         &
                        filename      = trim(filename),   &
                        mesh_topology = 'StructuredGrid', &
                        nx1 = ni1, nx2 = ni2,             &
                        ny1 = nj1, ny2 = nj2,             &
                        nz1 = nk1, nz2 = nk2)
    else
      err = VTK_INI_XML(output_format = 'ascii',          &
                        filename      = trim(filename),   &
                        mesh_topology = 'StructuredGrid', &
                        nx1 = ni1, nx2 = ni2,             &
                        ny1 = nj1, ny2 = nj2,             &
                        nz1 = nk1, nz2 = nk2)
    endif
    ! saving auxiliary data (time and time step)
    err = VTK_FLD_XML(fld_action='open')
    err = VTK_FLD_XML(fld=block%global%t,fname='TIME')
    err = VTK_FLD_XML(fld=block%global%n,fname='CYCLE')
    err = VTK_FLD_XML(fld_action='close')
    ! saving the geometry
    err = VTK_GEO_XML(nx1 = ni1, nx2 = ni2,                                            &
                      ny1 = nj1, ny2 = nj2,                                            &
                      nz1 = nk1, nz2 = nk2,                                            &
                      NN = nnode,                                                      &
                      X=reshape(block%node(ni1:ni2,nj1:nj2,nk1:nk2)%x,(/nnode/)), &
                      Y=reshape(block%node(ni1:ni2,nj1:nj2,nk1:nk2)%y,(/nnode/)), &
                      Z=reshape(block%node(ni1:ni2,nj1:nj2,nk1:nk2)%z,(/nnode/)))
    if (.not.meshonly) then
      ! saving dependent variables
      if (pp_format%node) then
        err = VTK_DAT_XML(var_location = 'node', var_block_action = 'open')
        do s=1,global%Ns
          err=VTK_VAR_XML(NC_NN=nnode,varname='r('//trim(str(.true.,s))//')',var=reshape(r(s,ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        enddo
        err=VTK_VAR_XML(NC_NN=nnode,varname='r',var=reshape(P(ni1:ni2,nj1:nj2,nk1:nk2)%d,  (/nnode/)))
        err=VTK_VAR_XML(NC_NN=nnode,varname='u',var=reshape(P(ni1:ni2,nj1:nj2,nk1:nk2)%v%x,(/nnode/)))
        err=VTK_VAR_XML(NC_NN=nnode,varname='v',var=reshape(P(ni1:ni2,nj1:nj2,nk1:nk2)%v%y,(/nnode/)))
        err=VTK_VAR_XML(NC_NN=nnode,varname='w',var=reshape(P(ni1:ni2,nj1:nj2,nk1:nk2)%v%z,(/nnode/)))
        err=VTK_VAR_XML(NC_NN=nnode,varname='p',var=reshape(P(ni1:ni2,nj1:nj2,nk1:nk2)%p,  (/nnode/)))
        err=VTK_VAR_XML(NC_NN=nnode,varname='g',var=reshape(P(ni1:ni2,nj1:nj2,nk1:nk2)%g,  (/nnode/)))
        err=VTK_DAT_XML(var_location ='node',var_block_action = 'close')
      else
        err = VTK_DAT_XML(var_location = 'cell', var_block_action = 'open')
        do s=1,global%Ns
          err=VTK_VAR_XML(NC_NN=ncell,varname='r('//trim(str(.true.,s))//')',var=reshape(r(s,ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
        enddo
        err=VTK_VAR_XML(NC_NN=ncell,varname='r',var=reshape(block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%d,  (/ncell/)))
        err=VTK_VAR_XML(NC_NN=ncell,varname='u',var=reshape(block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%x,(/ncell/)))
        err=VTK_VAR_XML(NC_NN=ncell,varname='v',var=reshape(block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%y,(/ncell/)))
        err=VTK_VAR_XML(NC_NN=ncell,varname='w',var=reshape(block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%z,(/ncell/)))
        err=VTK_VAR_XML(NC_NN=ncell,varname='p',var=reshape(block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%p,  (/ncell/)))
        err=VTK_VAR_XML(NC_NN=ncell,varname='g',var=reshape(block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%g,  (/ncell/)))
        err=VTK_DAT_XML(var_location ='cell',var_block_action = 'close')
      endif
    endif
    ! closing VTK file
    err = VTK_GEO_XML()
    err = VTK_END_XML()
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction vtk_blk_data
  endfunction vtk_output

  !> Function for writing OFF block data to Gnuplot file.
  function gnu_output(meshonly,block,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,           intent(IN):: meshonly  !< Flag for post-process only mesh.
  type(Type_SBlock), intent(IN):: block(1:) !< Block-level data.
  character(*),      intent(IN):: filename  !< File name of the output file.
  integer(I_P)::                  err       !< Error trapping flag: 0 no errors, >0 error occurs.
  type(Type_Global), pointer::    global    !< Global-level data.
  integer(I_P)::                  b         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  global => block(1)%global
  ! writing data blocks
  do b=1,global%Nb
    err = gnu_blk_data(block = block(b), filename = trim(filename)//'_b'//trim(strz(4,b))//'.gnu.dat')
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    function gnu_blk_data(block,filename) result(err)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Function for writing block data.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    type(Type_SBlock), intent(IN):: block    ! Block-level data.
    character(*),      intent(IN):: filename ! File name of the output block file.
    integer(I_P)::                  err      ! Error trapping flag: 0 no errors, >0 error occurs.
    integer(I_P)::                  i,j,k    ! Counters.
    integer(I_P)::                  ugnu     ! Logic unit.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    open(unit=Get_Unit(ugnu),file=trim(filename))
    if (meshonly) then
      do k=1,block%Nk
        do j=1,block%Nj
          do i=1,block%Ni
            write(ugnu,'(3('//FR8P//',1X))',iostat=err) block%C(i,j,k)%cent%x,block%C(i,j,k)%cent%y,block%C(i,j,k)%cent%z
          enddo
        enddo
      enddo
    else
      do k=1,block%Nk
        do j=1,block%Nj
          do i=1,block%Ni
            write(ugnu,'(9('//FR8P//',1X))',iostat=err) block%C(i,j,k)%cent%x, &
                                                        block%C(i,j,k)%cent%y, &
                                                        block%C(i,j,k)%cent%z, &
                                                        block%C(i,j,k)%P%d,    &
                                                        block%C(i,j,k)%P%v%x,  &
                                                        block%C(i,j,k)%P%v%y,  &
                                                        block%C(i,j,k)%P%v%z,  &
                                                        block%C(i,j,k)%P%p,    &
                                                        block%C(i,j,k)%P%g
          enddo
        enddo
      enddo
    endif
    close(ugnu)
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction gnu_blk_data
  endfunction gnu_output
  !> @}
endmodule Lib_PostProcessing
