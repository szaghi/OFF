program test
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision
USE Data_Type_BC
!USE Data_Type_BC_Adjacent
USE Data_Type_Files
USE Data_Type_Global
USE Lib_IO_Misc
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!call test_file_ibm_opts
!call test_file_mesh
call test_ibm
!call test_file_fluid
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine test_ibm()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Files)::         IOFile
  type(Type_Global),target:: global
  integer(I4P)::             b,l
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------

  associate(OS=>global%OS,mesh_dims=>global%mesh_dims,ibm_opts=>IOFile%ibm_opts)
    call ibm_opts%set(name='ibm_options.opt',path_in='./input/')
    call ibm_opts%load(OS=OS,Nl=mesh_dims%Nl) ; if (ibm_opts%iostat/=0) write(stderr,'(A)')trim(ibm_opts%iomsg)

    call IOFile%spec%set(         name=ibm_opts%fn_spec,          path_in =ibm_opts%path_in )
    call IOFile%ibm_blks_desc%set(name=ibm_opts%fn_blk_desc(1)%vs,path_in =ibm_opts%path_in )
    call IOFile%mesh%set(         name=ibm_opts%fn_mesh,          path_out=ibm_opts%path_out)
    call IOFile%bc%set(           name=ibm_opts%fn_bc  ,          path_out=ibm_opts%path_out)
    call IOFile%init%set(         name=ibm_opts%fn_init,          path_out=ibm_opts%path_out)
  endassociate

  call IOFile%spec%load(species=global%species0) ; if (IOFile%spec%iostat/=0) write(stderr,'(A)')trim(IOFile%spec%iomsg)

  associate(Nl=>global%mesh_dims%Nl,ibm_blks_desc=>IOFile%ibm_blks_desc)
    call ibm_blks_desc%load ; if (ibm_blks_desc%iostat/=0) write(stderr,'(A)')trim(ibm_blks_desc%iomsg)
  endassociate

  global%mesh_dims%Nb_tot = IOFile%ibm_blks_desc%mesh_dims%Nb

  associate(mesh_dims=>global%mesh_dims,parallel=>global%parallel,proc=>IOFile%proc)
    call proc%set(name='procmap-no_mpi.dat',path_in='./input/',path_out=IOFile%ibm_opts%path_out)
    call proc%load(mesh_dims=mesh_dims,parallel=parallel) ; if (proc%iostat/=0) write(stderr,'(A)')trim(proc%iomsg)
    call global%parallel%print(unit=stdout)
    call proc%save(parallel=parallel) ; if (proc%iostat/=0) write(stderr,'(A)')trim(proc%iomsg)
  endassociate

  call global%parallel%set_serial(Nb_tot=global%mesh_dims%Nb_tot)

  ! computing the mesh dimensions
  call global%mesh_dims%alloc
  do l=1,global%mesh_dims%Nl ; do b=1,global%mesh_dims%Nb_tot
    if     (mod(IOFile%ibm_blks_desc%block(b)%dims%Ni,(2**(l-1)))/=0) then
      write(stderr,'(A)')' Attention the number of grid levels used is not consistent with the number of cells'
      write(stderr,'(A)')' Impossible to compute grid level '//trim(str(.true.,l))//' of block '//trim(str(.true.,b))
      write(stderr,'(A)')' Inconsistent direction i, Ni '//trim(str(.true.,IOFile%ibm_blks_desc%block(b)%dims%Ni))
      stop
    elseif (mod(IOFile%ibm_blks_desc%block(b)%dims%Nj,(2**(l-1)))/=0) then
      write(stderr,'(A)')' Attention the number of grid levels used is not consistent with the number of cells'
      write(stderr,'(A)')' Impossible to compute grid level '//trim(str(.true.,l))//' of block '//trim(str(.true.,b))
      write(stderr,'(A)')' Inconsistent direction j, Nj '//trim(str(.true.,IOFile%ibm_blks_desc%block(b)%dims%Nj))
      stop
    elseif (mod(IOFile%ibm_blks_desc%block(b)%dims%Nk,(2**(l-1)))/=0) then
      write(stderr,'(A)')' Attention the number of grid levels used is not consistent with the number of cells'
      write(stderr,'(A)')' Impossible to compute grid level '//trim(str(.true.,l))//' of block '//trim(str(.true.,b))
      write(stderr,'(A)')' Inconsistent direction k, Nk '//trim(str(.true.,IOFile%ibm_blks_desc%block(b)%dims%Nk))
      stop
    endif
    global%mesh_dims%block_dims(b,l)%gc = IOFile%ibm_blks_desc%block(b)%dims%gc
    global%mesh_dims%block_dims(b,l)%Ni = IOFile%ibm_blks_desc%block(b)%dims%Ni/(2**(l-1))
    global%mesh_dims%block_dims(b,l)%Nj = IOFile%ibm_blks_desc%block(b)%dims%Nj/(2**(l-1))
    global%mesh_dims%block_dims(b,l)%Nk = IOFile%ibm_blks_desc%block(b)%dims%Nk/(2**(l-1))
  enddo ; enddo
  call global%set_Ns_from_species0
  global%mesh_dims%set = .true.
  call global%alloc
  call global%block%alloc
  call global%block%alloc(members=.true.)
  ! computing node coordinates of finest level
  do b=1,global%mesh_dims%Nb_tot
    associate(block => global%block(b,1), extents => IOFile%ibm_blks_desc%block(b)%exts)
      call block%create_uniform_grid(extents=extents)
    endassociate
  enddo
  ! computing the node of other grid levels if used
  if (global%mesh_dims%Nl>1) then
    write(stderr,'(A)')' Attention: not yet tested'
    do l=2,global%mesh_dims%Nl ; do b=1,global%mesh_dims%Nb_tot
      associate(block => global%block(b,l), block_f => global%block(b,l-1))
        call block%create_grid_from_finer(block_f=block_f)
      endassociate
    enddo ; enddo
  endif
  ! boundary conditions setting
  associate(mesh_dims => global%mesh_dims)
    do l=1,mesh_dims%Nl ; do b=1,mesh_dims%Nb_tot
      associate(block => global%block(b,l), block_bc => IOFile%ibm_blks_desc%block(b)%bc%bc)
        block%F%bc=block_bc ; call block%set_cells_bc(mesh_dims=mesh_dims,l=l)
      endassociate
    enddo ; enddo
  endassociate
  ! initial conditions setting
  do l=1,global%mesh_dims%Nl ; do b=1,global%mesh_dims%Nb_tot
    global%block(b,l)%C%P = IOFile%ibm_blks_desc%block(b)%prim
  enddo ; enddo
  call IOFile%mesh%save(global=global) ; if (IOFile%mesh%iostat/=0) write(stderr,'(A)')trim(IOFile%mesh%iomsg)
  call IOFile%bc%save(  global=global) ; if (IOFile%bc%iostat  /=0) write(stderr,'(A)')trim(IOFile%bc%iomsg  )
  call IOFile%init%save(global=global) ; if (IOFile%init%iostat/=0) write(stderr,'(A)')trim(IOFile%init%iomsg)
  ! printing block infos
  write(stdout,'(A)')' Blocks infos'
  do l=1,global%mesh_dims%Nl ; do b=1,global%mesh_dims%Nb_tot
    write(stdout,'(A)')
    write(stdout,'(A)')'+-> Block b='//trim(str(n=b))//' level l='//trim(str(n=l))
    call global%block(b,l)%print(unit=stdout,pref='|->')
  enddo ; enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine test_ibm

 !subroutine test_file_ibm_opts()
 !!---------------------------------------------------------------------------------------------------------------------------------
 !implicit none
 !type(Type_Files)::         IOFile
 !type(Type_Global),target:: global
 !!---------------------------------------------------------------------------------------------------------------------------------
 !
 !!---------------------------------------------------------------------------------------------------------------------------------
 !write(stdout,'(A)')repeat('-',132)
 !write(stdout,'(A)')' Testing File IBM Options'
 !call IOFile%ibm_opts%set(name='ibm_options.opt',path_in='./input/',global=global)
 !call IOFile%ibm_opts%load
 !call IOFile%ibm_opts%save
 !write(stdout,'(A)')' path_in  '//IOFile%ibm_opts%path_in
 !write(stdout,'(A)')' path_out '//IOFile%ibm_opts%path_out
 !write(stdout,'(A)')' fn_mesh  '//IOFile%ibm_opts%fn_mesh
 !write(stdout,'(A)')' Nl       '//trim(str(n=IOFile%ibm_opts%global%mesh_dims%Nl))
 !write(stdout,'(A)')' fn_block '//IOFile%ibm_opts%fn_blk_desc(1)%vs
 !write(stdout,'(A)')repeat('-',132)
 !return
 !!---------------------------------------------------------------------------------------------------------------------------------
 !endsubroutine test_file_ibm_opts
 !
 !subroutine test_file_mesh()
 !!---------------------------------------------------------------------------------------------------------------------------------
 !implicit none
 !type(Type_Files)::         IOFile
 !type(Type_Global),target:: global
 !!---------------------------------------------------------------------------------------------------------------------------------
 !
 !!---------------------------------------------------------------------------------------------------------------------------------
 !write(stdout,'(A)')repeat('-',132)
 !write(stdout,'(A)')' Testing File Mesh'
 !call IOFile%mesh%set(name='sod.geo',path_in='./input/',path_out='./output/')
 !call IOFile%mesh%load(global=global)
 !call global%mesh_dims%print(unit=stdout)
 !call IOFile%mesh%save(global=global)
 !write(stdout,'(A)')repeat('-',132)
 !return
 !!---------------------------------------------------------------------------------------------------------------------------------
 !endsubroutine test_file_mesh
 !
 !subroutine test_file_fluid()
 !!---------------------------------------------------------------------------------------------------------------------------------
 !implicit none
 !type(Type_Files)::         IOFile
 !type(Type_Global),target:: global
 !!---------------------------------------------------------------------------------------------------------------------------------
 !
 !!---------------------------------------------------------------------------------------------------------------------------------
 !write(stdout,'(A)')repeat('-',132)
 !write(stdout,'(A)')' Testing File Fluid'
 !call IOFile%mesh%set(name='sod.geo',path_in='./input/',path_out='./output/',global=global)
 !call IOFile%mesh%load(global=global)
 !call global%mesh_dims%print(unit=stdout)
 !call IOFile%init%set(name='sod.itc',path_in='./input/',path_out='./output/',global=global)
 !call IOFile%init%load(global=global)
 !call global%species0%print(unit=stdout)
 !call IOFile%init%save(global=global)
 !call IOFile%bc%set(name='sod.bco',path_in='./input/',path_out='./output/',global=global)
 !!call IOFile%bc%load
 !call global%block(1,1)%C(10,1,1)%print(unit=stdout)
 !call global%block(1,1)%node(10,1,1)%print(unit=stdout)
 !call global%block(2,1)%C(10,1,1)%print(unit=stdout)
 !call global%block(2,1)%node(10,1,1)%print(unit=stdout)
 !
 !call global%block(1,1)%Fi(48,1,1)%BC%print(unit=stdout)
 !call global%block(1,1)%Fi(50,1,1)%BC%print(unit=stdout)
 !call global%block(2,1)%Fi(0,1,1)%BC%print(unit=stdout)
 !write(stdout,'(A)')repeat('-',132)
 !return
 !!---------------------------------------------------------------------------------------------------------------------------------
 !endsubroutine test_file_fluid

 !subroutine test_file_species()
 !!---------------------------------------------------------------------------------------------------------------------------------
 !implicit none
 !type(Type_File_Species):: file_spec
 !
 !!---------------------------------------------------------------------------------------------------------------------------------
 !call file_spec%set(name='initial_species.opt',path_in='./input/',path_out='./output/')
 !call file_spec%load
 !call file_spec%save
 !write(stdout,'(A)')' specie 1 '//trim(str(n=file_spec%species%heats(1)%cp))//' '//trim(str(n=file_spec%species%heats(1)%cv))
 !write(stdout,'(A)')' specie 2 '//trim(str(n=file_spec%species%heats(2)%cp))//' '//trim(str(n=file_spec%species%heats(2)%cv))
 !write(stdout,'(A)')' specie 3 '//trim(str(n=file_spec%species%heats(3)%cp))//' '//trim(str(n=file_spec%species%heats(3)%cv))
 !write(stdout,'(A)')' specie 4 '//trim(str(n=file_spec%species%heats(4)%cp))//' '//trim(str(n=file_spec%species%heats(4)%cv))
 !return
 !!---------------------------------------------------------------------------------------------------------------------------------
 !endsubroutine test_file_species

 !subroutine test_file_ibm_blk_desc()
 !!---------------------------------------------------------------------------------------------------------------------------------
 !implicit none
 !type(Type_File_Blocks_Cartesian):: file_ibm_blk_desc
 !!---------------------------------------------------------------------------------------------------------------------------------
 !
 !call file_ibm_blk_desc%set(name='blocks_description.opt',path_in='./input/')
 !call file_ibm_blk_desc%load
 !return
 !!---------------------------------------------------------------------------------------------------------------------------------
 !endsubroutine test_file_ibm_blk_desc

endprogram test
