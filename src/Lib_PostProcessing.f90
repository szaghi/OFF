!> @ingroup GlobalVarPar
!> @{
!> @defgroup Lib_PostProcessing Lib_PostProcessing
!> @}

!> This module contains the definition of procedures and variables useful for post-process the OFF data.
!> @ingroup Library
module Lib_PostProcessing
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision        ! Integers and reals precision definition.
USE Data_Type_Globals   ! Definition of Type_Global and Type_Block.
USE Data_Type_Primitive ! Definition of Type_Primitive.
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
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing the post-processing options.
!> @ingroup DerivedType
type, public:: Type_PP_Format
  logical:: binary = .true.  !< Binary or ascii post-process file.
  logical:: node   = .true.  !< Node or cell data location.
  logical:: bc     = .false. !< Saving or not boundary conditions cells.
  logical:: tec    = .true.  !< Tecplot file.
  logical:: vtk    = .false. !< VTK file.
endtype Type_PP_Format
!> @ingroup Lib_PostProcessing
!> @{
type(Type_PP_Format):: pp_format !< Post-processing format options (see \ref Lib_PostProcessing::Type_PP_Format "definition").
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine interpolate_primitive(global,block,P)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for interpolating primitive variables at nodes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global),    intent(IN)::    global               ! Global-level data.
  type(Type_Block),     intent(INOUT):: block                ! Block-level data.
  type(Type_Primitive), intent(INOUT):: P(0:block%mesh%Ni, &
                                          0:block%mesh%Nj, &
                                          0:block%mesh%Nk)   ! Primitive variables intepolated at nodes.
  real(R_P)::                           mf                   ! Mean factor.
  integer(I_P)::                        i,j,k,s              ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call init(Ns=global%fluid%Ns,prim=P)
#if !defined NULi && !defined NULj && !defined NULk
  ! 3D data
  mf = 0.125_R_P
#elif defined NULi
  block%fluid%P(0              ,:,:) = 0._R_P
  block%fluid%P(block%mesh%Ni+1,:,:) = 0._R_P
#if !defined NULj && !defined NULk
  ! 2D data
  mf = 0.25_R_P
#elif defined NULj
  ! 1D data
  mf = 0.5_R_P
  block%fluid%P(:,0              ,:) = 0._R_P
  block%fluid%P(:,block%mesh%Nj+1,:) = 0._R_P
#elif defined NULk
  ! 1D data
  mf = 0.5_R_P
  block%fluid%P(:,:,0              ) = 0._R_P
  block%fluid%P(:,:,block%mesh%Nk+1) = 0._R_P
#endif
#elif defined NULj
  block%fluid%P(:,0              ,:) = 0._R_P
  block%fluid%P(:,block%mesh%Nj+1,:) = 0._R_P
#if !defined NULi && !defined NULk
  ! 2D data
  mf = 0.25_R_P
#elif defined NULi
  ! 1D data
  mf = 0.5_R_P
  block%fluid%P(0              ,:,:) = 0._R_P
  block%fluid%P(block%mesh%Ni+1,:,:) = 0._R_P
#elif defined NULk
  ! 1D data
  mf = 0.5_R_P
  block%fluid%P(:,:,0              ) = 0._R_P
  block%fluid%P(:,:,block%mesh%Nk+1) = 0._R_P
#endif
#elif defined NULk
  block%fluid%P(:,:,0              ) = 0._R_P
  block%fluid%P(:,:,block%mesh%Nk+1) = 0._R_P
#if !defined NULi && !defined NULj
  ! 2D data
  mf = 0.25_R_P
#elif defined NULi
  ! 1D data
  mf = 0.5_R_P
  block%fluid%P(0              ,:,:) = 0._R_P
  block%fluid%P(block%mesh%Ni+1,:,:) = 0._R_P
#elif defined NULj
  ! 1D data
  mf = 0.5_R_P
  block%fluid%P(:,0              ,:) = 0._R_P
  block%fluid%P(:,block%mesh%Nj+1,:) = 0._R_P
#endif
#endif
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k,s)       &
  !$OMP SHARED(global,block,P,mf)
  !$OMP DO
  do k=0,block%mesh%Nk
    do j=0,block%mesh%Nj
      do i=0,block%mesh%Ni
          P(i,j,k) = block%fluid%P(i+1,j+1,k+1) + block%fluid%P(i,j+1,k+1) &
                   + block%fluid%P(i+1,j  ,k+1) + block%fluid%P(i,j,  k+1) &
                   + block%fluid%P(i+1,j+1,k  ) + block%fluid%P(i,j+1,k  ) &
                   + block%fluid%P(i+1,j  ,k  ) + block%fluid%P(i,j  ,k  )
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
  type(Type_Block), intent(IN)::  block                   ! Block-level data.
  integer(I_P),     intent(OUT):: ni1,ni2,nj1,nj2,nk1,nk2 ! Bounds of dimensions of node-centered data.
  integer(I_P),     intent(OUT):: ci1,ci2,cj1,cj2,ck1,ck2 ! Bounds of dimensions of cell-centered data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (pp_format%node) then
    ni1 = 0 ; ni2 = block%mesh%Ni
    nj1 = 0 ; nj2 = block%mesh%Nj
    nk1 = 0 ; nk2 = block%mesh%Nk

    ci1 = 0 ; ci2 = block%mesh%Ni
    cj1 = 0 ; cj2 = block%mesh%Nj
    ck1 = 0 ; ck2 = block%mesh%Nk
  else
    if (pp_format%bc) then
      ni1 = 0 - block%mesh%gc(1) ; ni2 = block%mesh%Ni + block%mesh%gc(2)
      nj1 = 0 - block%mesh%gc(3) ; nj2 = block%mesh%Nj + block%mesh%gc(4)
      nk1 = 0 - block%mesh%gc(5) ; nk2 = block%mesh%Nk + block%mesh%gc(6)

      ci1 = 1 - block%mesh%gc(1) ; ci2 = block%mesh%Ni + block%mesh%gc(2)
      cj1 = 1 - block%mesh%gc(3) ; cj2 = block%mesh%Nj + block%mesh%gc(4)
      ck1 = 1 - block%mesh%gc(5) ; ck2 = block%mesh%Nk + block%mesh%gc(6)
    else
      ni1 = 0                    ; ni2 = block%mesh%Ni
      nj1 = 0                    ; nj2 = block%mesh%Nj
      nk1 = 0                    ; nk2 = block%mesh%Nk

      ci1 = 1                    ; ci2 = block%mesh%Ni
      cj1 = 1                    ; cj2 = block%mesh%Nj
      ck1 = 1                    ; ck2 = block%mesh%Nk
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_dimensions

  !> Function for writing OFF block data to Tecplot file.
  function tec_output(meshonly,global,block,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,           intent(IN)::    meshonly                           !< Flag for post-process only mesh.
  type(Type_Global), intent(IN)::    global                             !< Global-level data.
  type(Type_Block),  intent(INOUT):: block(1:global%mesh%Nb)            !< Block-level data.
  character(*),      intent(IN)::    filename                           !< File name of the output file.
  integer(I_P)::                     err                                !< Error trapping flag: 0 no errors, >0 error occurs.
#ifdef TECIO
  integer(I_P), external::           tecini112,    &                    ! |
                                     tecauxstr112, &                    ! |
                                     teczne112,    &                    ! | Tecplot external functions.
                                     tecdat112,    &                    ! |
                                     tecend112                          ! |
#endif
  character(1), parameter::          tecendrec = char(0)                !< End-character for binary-record end.
  character(500)::                   tecvarname                         !< Variables name for tecplot header file.
  character(500)::                   teczoneheader                      !< Tecplot string of zone header.
  character(500)::                   tecvarform                         !< Format for variables for tecplot file.
  integer(I_P)::                     tecvarloc(1:3+global%fluid%Np)     !< Tecplot array of variables location.
  character(500)::                   tecvarlocstr                       !< Tecplot string of variables location.
  integer(I_P)::                     tecnull(1:3+global%fluid%Np)       !< Tecplot null array.
  integer(I_P)::                     tecunit                            !< Free logic unit of tecplot file.
  integer(I_P)::                     nvar                               !< Number of variables saved.
  integer(I_P)::                     b                                  !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
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
    err = tecauxstr112("Time"//tecendrec,trim(str(n=global%fluid%t))//tecendrec)
#else
    write(stderr,'(A)') 'Error: your are trying to save binary tecplot file without compiling against the Tecplot library.'
    stop
#endif
  else
    tecunit = Get_Unit()
    if (filename(len_trim(filename)-4:len_trim(filename))/=".dat") then
      open(unit=tecunit,file=trim(filename)//".dat")
    else
      open(unit=tecunit,file=trim(filename))
    endif
    write(tecunit,'(A)',iostat=err)trim(tecvarname)
  endif
  ! writing data blocks
  do b=1,global%mesh%Nb
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
        do s=1,global%fluid%Ns
          tecvarname = trim(tecvarname)//' r_'//trim(str(.true.,s))
        enddo
        tecvarname = trim(tecvarname)//' r u v w p g'
      endif
      ! variables location
      if (pp_format%node) then
        tecvarloc = 1
      else
        tecvarloc(1:3) = 1 ; tecvarloc(4:3+global%fluid%Np)= 0
      endif
      ! null array
      tecnull = 0
      if (.not.meshonly) then
        nvar = 3 + global%fluid%Np
      else
        nvar = 3
      endif
    else
      ! header variables names
      tecvarname = ' VARIABLES ="x" "y" "z"'
      if (.not.meshonly) then
        do s=1,global%fluid%Ns
          tecvarname = trim(tecvarname)//' "r_'//trim(str(.true.,s))//'"'
        enddo
        tecvarname = trim(tecvarname)//' "r" "u" "v" "w" "p" "g"'
      endif
      ! variables output format
      if (.not.meshonly) then
        write(tecvarform,'(A)')'('//trim(str(no_sign=.true.,n=global%fluid%Np))//'('//FR_P//',1X))'
      else
        write(tecvarform,'(A)')'('//trim(str(no_sign=.true.,n=3))//'('//FR_P//',1X))'
      endif
      ! variables location
      if (.not.meshonly) then
        if (pp_format%node) then
          tecvarlocstr = ', VARLOCATION=([1-'//trim(str(.true.,global%fluid%Np))//']=NODAL)'
        else
          tecvarlocstr = ', VARLOCATION=([1-3]=NODAL,[4-'//trim(str(.true.,global%fluid%Np))//']=CELLCENTERED)'
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
    integer(I_P),      intent(IN)::    b                       ! Block number.
    type(Type_Global), intent(IN)::    global                  ! Global-level data.
    type(Type_Block),  intent(INOUT):: block                   ! Block-level data.
    integer(I_P)::                     err                     ! Error trapping flag: 0 no errors, >0 error occurs.
    type(Type_Primitive)::             P(0:block%mesh%Ni, &    ! |
                                         0:block%mesh%Nj, &    ! | Primitive variables intepolated at nodes.
                                         0:block%mesh%Nk)      ! |
    real(R_P)::                        r(1:global%fluid%Ns, &  ! |
                                         0:block%mesh%Ni,   &  ! | Array containing species densities.
                                         0:block%mesh%Nj,   &  ! |
                                         0:block%mesh%Nk)      ! |
    integer(I_P)::                     ni1,ni2,nj1,nj2,nk1,nk2 ! Bounds of dimensions of node-centered data.
    integer(I_P)::                     ci1,ci2,cj1,cj2,ck1,ck2 ! Bounds of dimensions of cell-centered data.
    integer(I_P)::                     nnode,ncell             ! Number of nodes and cells.
    integer(I_P)::                     i,j,k,s                 ! Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (pp_format%node) then
      ! tri-linear interpolation of cell-centered values at nodes
      call interpolate_primitive(global=global,block=block,P=P)
    endif
    ! initialize the zone dimensions
    call compute_dimensions(block=block,                                     &
                            ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2, &
                            ci1=ci1,ci2=ci2,cj1=cj1,cj2=cj2,ck1=ck1,ck2=ck2)
    nnode = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1)
    ncell = (ci2-ci1+1)*(cj2-cj1+1)*(ck2-ck1+1)
    ! storing species densities into an array for avoiding problems with nonzero rank pointer
    if (pp_format%node) then
      do k=0,block%mesh%Nk
        do j=0,block%mesh%Nj
          do i=0,block%mesh%Ni
            do s=1,global%fluid%Ns
              r(s,i,j,k) = P(i,j,k)%r(s)
            enddo
          enddo
        enddo
      enddo
    else
      do k=1,block%mesh%Nk
        do j=1,block%mesh%Nj
          do i=1,block%mesh%Ni
            do s=1,global%fluid%Ns
              r(s,i,j,k) = block%fluid%P(i,j,k)%r(s)
            enddo
          enddo
        enddo
      enddo
    endif
    ! writing the block data
    if (pp_format%binary) then
#ifdef TECIO
      err = teczne112('n_'//trim(strz(nz_pad=10,n=global%fluid%n))//'-b_'//trim(strz(nz_pad=4,n=b))//tecendrec, &
                      0,                                                                                        &
                      ni2-ni1+1,                                                                                &
                      nj2-nj1+1,                                                                                &
                      nk2-nk1+1,                                                                                &
                      0,                                                                                        &
                      0,                                                                                        &
                      0,                                                                                        &
                      global%fluid%t,                                                                           &
                      0,                                                                                        &
                      0,                                                                                        &
                      1,                                                                                        & !1=>block,0=>point
                      0,                                                                                        &
                      0,                                                                                        &
                      0,                                                                                        &
                      0,                                                                                        &
                      0,                                                                                        &
                      tecnull(1:nvar),                                                                          &
                      tecvarloc(1:nvar),                                                                        &
                      tecnull(1:nvar),                                                                          &
                      0)
#else
      write(stderr,'(A)') 'Error: your are trying to save binary tecplot file without compiling against the Tecplot library.'
      stop
#endif
      err=tec_dat(N=nnode,dat=block%mesh%node(ni1:ni2,nj1:nj2,nk1:nk2)%x)
      err=tec_dat(N=nnode,dat=block%mesh%node(ni1:ni2,nj1:nj2,nk1:nk2)%y)
      err=tec_dat(N=nnode,dat=block%mesh%node(ni1:ni2,nj1:nj2,nk1:nk2)%z)
      if (.not.meshonly) then
        if (pp_format%node) then
          do s=1,global%fluid%Ns
            err=tec_dat(N=ncell,dat=r(s,ni1:ni2,nj1:nj2,nk1:nk2))
          enddo
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%d  )
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%v%x)
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%v%y)
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%v%z)
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%p  )
          err=tec_dat(N=nnode,dat=P(ni1:ni2,nj1:nj2,nk1:nk2)%g  )
        else
          do s=1,global%fluid%Ns
            err=tec_dat(N=ncell,dat=r(s,ci1:ci2,cj1:cj2,ck1:ck2))
          enddo
          err=tec_dat(N=ncell,dat=block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%d  )
          err=tec_dat(N=ncell,dat=block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%v%x)
          err=tec_dat(N=ncell,dat=block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%v%y)
          err=tec_dat(N=ncell,dat=block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%v%z)
          err=tec_dat(N=ncell,dat=block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%p  )
          err=tec_dat(N=ncell,dat=block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%g  )
        endif
      endif
    else
      ! tecplot zone header
      teczoneheader = ' ZONE  T = "n_'//trim(strz(nz_pad=10,n=global%fluid%n))// &
                                 '-b_'//trim(strz(nz_pad=4,n=b))//'"'//          &
                      ', I='//trim(str(no_sign=.true.,n=ni2-ni1+1))//            &
                      ', J='//trim(str(no_sign=.true.,n=nj2-nj1+1))//            &
                      ', K='//trim(str(no_sign=.true.,n=nk2-nk1+1))//            &
                      ', DATAPACKING=BLOCK'//adjustl(trim(tecvarlocstr))
      write(tecunit,'(A)',iostat=err)trim(teczoneheader)
      write(tecunit,FR_P,iostat=err)(((block%mesh%node(i,j,k)%x,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
      write(tecunit,FR_P,iostat=err)(((block%mesh%node(i,j,k)%y,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
      write(tecunit,FR_P,iostat=err)(((block%mesh%node(i,j,k)%z,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
      if (.not.meshonly) then
        if (pp_format%node) then
          do s=1,global%fluid%Ns
            write(tecunit,FR_P,iostat=err)(((P(i,j,k)%r(s),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          enddo
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%d  ,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%v%x,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%v%y,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%v%z,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%p  ,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
          write(tecunit,FR_P,iostat=err)(((P(i,j,k)%g  ,i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
        else
          do s=1,global%fluid%Ns
            write(tecunit,FR_P,iostat=err)(((block%fluid%P(i,j,k)%r(s),i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          enddo
          write(tecunit,FR_P,iostat=err)(((block%fluid%P(i,j,k)%d  ,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          write(tecunit,FR_P,iostat=err)(((block%fluid%P(i,j,k)%v%x,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          write(tecunit,FR_P,iostat=err)(((block%fluid%P(i,j,k)%v%y,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          write(tecunit,FR_P,iostat=err)(((block%fluid%P(i,j,k)%v%z,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          write(tecunit,FR_P,iostat=err)(((block%fluid%P(i,j,k)%p  ,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
          write(tecunit,FR_P,iostat=err)(((block%fluid%P(i,j,k)%g  ,i=ci1,ci2),j=cj1,cj2),k=ck1,ck2)
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
  function vtk_output(meshonly,global,block,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,           intent(IN)::    meshonly                !< Flag for post-process only mesh.
  type(Type_Global), intent(IN)::    global                  !< Global-level data.
  type(Type_Block),  intent(INOUT):: block(1:global%mesh%Nb) !< Block-level data.
  character(*),      intent(IN)::    filename                !< File name of the output file.
  integer(I_P)::                     err                     !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     b                       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! writing data blocks
  do b=1,global%mesh%Nb
    err = vtk_blk_data( global = global, block = block(b), filename = trim(filename)//'_b'//trim(strz(4,b))//'.vts')
  enddo
  !! writing the VTM wrapper for multiblock dataset
  !err = VTM_INI_XML(filename=trim(filename)//'.vtm')
  !err = VTM_BLK_XML(block_action='open')
  !err = VTM_WRF_XML(vtk_xml_file_list=(/(trim(filename)//'_b'//trim(strz(4,b))//'.vts',b=1,global%mesh%Nb)/))
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
    type(Type_Global), intent(IN)::    global                  ! Global-level data.
    type(Type_Block),  intent(INOUT):: block                   ! Block-level data.
    character(*),      intent(IN)::    filename                ! File name of the output block file.
    integer(I_P)::                     err                     ! Error trapping flag: 0 no errors, >0 error occurs.
    type(Type_Primitive)::             P(0:block%mesh%Ni,   &  ! |
                                         0:block%mesh%Nj,   &  ! | Primitive variables intepolated at nodes.
                                         0:block%mesh%Nk)      ! |
    real(R_P)::                        r(1:global%fluid%Ns, &  ! |
                                         0:block%mesh%Ni,   &  ! | Array containing species densities.
                                         0:block%mesh%Nj,   &  ! |
                                         0:block%mesh%Nk)      ! |
    integer(I_P)::                     ni1,ni2,nj1,nj2,nk1,nk2 ! Bounds of dimensions of node-centered data.
    integer(I_P)::                     ci1,ci2,cj1,cj2,ck1,ck2 ! Bounds of dimensions of cell-centered data.
    integer(I_P)::                     nnode,ncell             ! Number of nodes and cells.
    integer(I_P)::                     i,j,k,s                 ! Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (pp_format%node) then
      ! tri-linear interpolation of cell-centered values at nodes
      call interpolate_primitive(global=global,block=block,P=P)
    endif
    ! initialize the block dimensions
    call compute_dimensions(block=block,                                     &
                            ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2, &
                            ci1=ci1,ci2=ci2,cj1=cj1,cj2=cj2,ck1=ck1,ck2=ck2)
    nnode = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1)
    ncell = (ci2-ci1+1)*(cj2-cj1+1)*(ck2-ck1+1)
    ! storing species densities into an array for avoiding problems with nonzero rank pointer
    if (pp_format%node) then
      do k=0,block%mesh%Nk
        do j=0,block%mesh%Nj
          do i=0,block%mesh%Ni
            do s=1,global%fluid%Ns
              r(s,i,j,k) = P(i,j,k)%r(s)
            enddo
          enddo
        enddo
      enddo
    else
      do k=1,block%mesh%Nk
        do j=1,block%mesh%Nj
          do i=1,block%mesh%Ni
            do s=1,global%fluid%Ns
              r(s,i,j,k) = block%fluid%P(i,j,k)%r(s)
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
    ! saving the geometry
    err = VTK_GEO_XML(nx1 = ni1, nx2 = ni2,                                            &
                      ny1 = nj1, ny2 = nj2,                                            &
                      nz1 = nk1, nz2 = nk2,                                            &
                      NN = nnode,                                                      &
                      X=reshape(block%mesh%node(ni1:ni2,nj1:nj2,nk1:nk2)%x,(/nnode/)), &
                      Y=reshape(block%mesh%node(ni1:ni2,nj1:nj2,nk1:nk2)%y,(/nnode/)), &
                      Z=reshape(block%mesh%node(ni1:ni2,nj1:nj2,nk1:nk2)%z,(/nnode/)))
    if (.not.meshonly) then
      ! saving dependent variables
      if (pp_format%node) then
        err = VTK_DAT_XML(var_location = 'node', var_block_action = 'open')
        do s=1,global%fluid%Ns
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
        do s=1,global%fluid%Ns
          err=VTK_VAR_XML(NC_NN=ncell,varname='r('//trim(str(.true.,s))//')',var=reshape(r(s,ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
        enddo
        err=VTK_VAR_XML(NC_NN=ncell,varname='r',var=reshape(block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%d,  (/ncell/)))
        err=VTK_VAR_XML(NC_NN=ncell,varname='u',var=reshape(block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%v%x,(/ncell/)))
        err=VTK_VAR_XML(NC_NN=ncell,varname='v',var=reshape(block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%v%y,(/ncell/)))
        err=VTK_VAR_XML(NC_NN=ncell,varname='w',var=reshape(block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%v%z,(/ncell/)))
        err=VTK_VAR_XML(NC_NN=ncell,varname='p',var=reshape(block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%p,  (/ncell/)))
        err=VTK_VAR_XML(NC_NN=ncell,varname='g',var=reshape(block%fluid%P(ci1:ci2,cj1:cj2,ck1:ck2)%g,  (/ncell/)))
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
endmodule Lib_PostProcessing
