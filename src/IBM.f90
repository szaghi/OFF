!> @brief IBM, Initial and Boundary conditions, Mesh generator for @off (Open Finite volume Fluid dynamics code).
!> This is an auxiliary tool useful for building proper inputs for @off code. ICG can build Initial Conditions files, Mesh files and
!> Boundary Conditions files. It accepts two kinds of inputs: \n
!> - Direct Blocks Description: this is the simplest available input. The initial and boundary descriptions as well as the geometry
!>   are directly described by means of simple ascii files. This kind of inputs can describe only simple Cartesian grids.
!> - Ansys (http://www.ansys.com) IcemCFD Multiblock INFO importer: this is a more complex (but more flexible) input. The initial
!>   conditions are described by means of simple ascii files similar to the Direct Block Description input, but the boundary
!>   conditions and the geometry are loaded by Ansys IcemCFD Multiblock INFO files. These files can describe more complex scenario
!>   with general curvilinear grids.
!> @todo \b DocImprove: Improve the documentation
!> @ingroup Program
program IBM
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                           ! Integers and reals precision definition.
USE Data_Type_BC, init_bc=>init, set_bc=>set               ! Definition of Type_BC.
USE Data_Type_Cell, init_cell=>init, set_cell=>set         ! Definition of Type_Cell.
USE Data_Type_Conservative, init_cons=>init, set_cons=>set ! Definition of Type_Conservative.
USE Data_Type_Globals                                      ! Definition of Type_Global and Type_Block.
USE Data_Type_OS, init_OS=>init, set_OS=>set               ! Definition of Type_OS.
USE Data_Type_Primitive, init_prim=>init, set_prim=>set    ! Definition of Type_Primitive.
USE Data_Type_Time, only: Get_Date_String                  ! Function for getting actual date.
USE Data_Type_Vector, set_Vector=>set                      ! Definition of Type_Vector.
USE Lib_IO_Misc                                            ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Global)::             global         ! Global-level data.
type(Type_Block), allocatable:: block(:,:)     ! Block-level data [1:Nb,1:Nl].
!> Derived type containing blocks informations.
type:: Type_Blocks
  integer(I1P)::         gc(1:6) !< Ghost cells.
  integer(I_P)::         Ni      !< Number of cells in i direction.
  integer(I_P)::         Nj      !< Number of cells in j direction.
  integer(I_P)::         Nk      !< Number of cells in k direction.
  real(R_P)::            xmin    !< Minimum value of x abscissa.
  real(R_P)::            ymin    !< Minimum value of y abscissa.
  real(R_P)::            zmin    !< Minimum value of z abscissa.
  real(R_P)::            xmax    !< Maximum value of x abscissa.
  real(R_P)::            ymax    !< Maximum value of y abscissa.
  real(R_P)::            zmax    !< Maximum value of z abscissa.
  type(Type_BC)::        bc(1:6) !< Block faces boundary conditions.
  type(Type_Primitive):: P       !< Blocks primitive variables.
endtype Type_Blocks
type(Type_Blocks), allocatable:: blocks(:)           ! Blocks data [1:Nb].
character(7)::                   In_type = 'ICEMCFD' ! Input type: 'BLOCKS' use direct blocks description, 'ICEMCFD' geo/topo files.
integer(I_P)::                   err                 ! Error trapping flag: 0 no errors, >0 error occurs.
character(20)::                  date                ! Actual date.
integer(I_P),      allocatable:: UnitScratch(:,:,:)  ! Free logic units for scratch files [1:3,1:Nb,1:Nl].
character(60)::                  File_Blocks         ! Blocks option file name.
integer(I_P)::                   Nf_Icem             ! Number of icemcfd files.
character(60),     allocatable:: File_Icem(:)        ! Icemcfd files name [1:Nf_Icem].
integer(I_P)::                   myrank = 0_I_P      ! Actual rank process.
integer(I_P)::                   b,l                 ! Counters for blocks and grid levels.
integer(I_P)::                   i,j,k               ! Counters for i, j, k directions.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
write(stdout,'(A)')' ICG started on'
date = Get_Date_String()
write(stdout,'(A)')' '//date
write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
write(stdout,*)

! initializing the grid generation
call icg_init(global)

! loading initial species
write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
write(stdout,'(A)',iostat=err)' Loading '//adjustl(trim(global%file%Path_InPut))//trim(global%file%File_Spec)
err = load_gfluid_0species(myrank=myrank,filename=adjustl(trim(global%file%Path_InPut))//trim(global%file%File_Spec),global=global)

! computing the mesh, boundary and initial conditions and storing in scratch files
select case(trim(In_type))
case('BLOCKS')
  ! loading blocks data
  write(stdout,'(A)',iostat=err)' Loading '//adjustl(trim(global%file%Path_InPut))//trim(File_Blocks)
  err = load_blocks(filename = adjustl(trim(global%file%Path_InPut))//trim(File_Blocks))
case('ICEMCFD')
  ! loading icemcfd data
  write(stdout,'(A)',iostat=err)' Loading Icemcfd data'
  err = load_icemcfd(Nf = Nf_Icem, filenames = adjustl(trim(global%file%Path_InPut))//File_Icem(1:Nf_Icem))
endselect
write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
! the number of global blocks is artificially set to 1
global%mesh%Nb = 1
! allocating block-level data
if (allocated(block)) then
  do l=lbound(block,dim=2),ubound(block,dim=2)
    do b=lbound(block,dim=1),ubound(block,dim=1)
      call free_block(block=block(b,l))
    enddo
  enddo
  deallocate(block)
endif
allocate(block(1:global%mesh%Nb,1:global%mesh%Nl))
do b=1,global%mesh%Nb_tot
  write(stdout,'(A)',iostat=err)' Saving output files of block '//str(.true.,b)
  ! setting mesh dimensions for grid level 1
  block(1,1)%mesh%Ni = blocks(b)%Ni
  block(1,1)%mesh%Nj = blocks(b)%Nj
  block(1,1)%mesh%Nk = blocks(b)%Nk
  block(1,1)%mesh%gc = blocks(b)%gc
  ! computing mesh dimensions for other grid levels
  if (global%mesh%Nl>1) then
    do l=2,global%mesh%Nl
      ! computing the number of cells of coarser grid levels
      if     (mod(block(1,l-1)%mesh%Ni,2)/=0) then
        write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
        write(stderr,'(A,I4)')' Impossible to compute grid level ',l
        write(stderr,'(A,I4)')' Inconsistent direction i, Ni ',block(1,l-1)%mesh%Ni
        write(stderr,'(A,I4)')' level ',l-1
        write(stderr,'(A,I4)')' block ',b
        stop
      elseif (mod(block(1,l-1)%mesh%Nj,2)/=0) then
        write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
        write(stderr,'(A,I4)')' Impossible to compute grid level ',l
        write(stderr,'(A,I4)')' Inconsistent direction j, Nj ',block(1,l-1)%mesh%Nj
        write(stderr,'(A,I4)')' level ',l-1
        write(stderr,'(A,I4)')' block ',b
        stop
      elseif (mod(block(1,l-1)%mesh%Nk,2)/=0) then
        write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
        write(stderr,'(A,I4)')' Impossible to compute grid level ',l
        write(stderr,'(A,I4)')' Inconsistent direction k, Nk ',block(1,l-1)%mesh%Nk
        write(stderr,'(A,I4)')' level ',l-1
        write(stderr,'(A,I4)')' block ',b
        stop
      endif
      block(1,l)%mesh%Ni = block(1,l-1)%mesh%Ni/2
      block(1,l)%mesh%Nj = block(1,l-1)%mesh%Nj/2
      block(1,l)%mesh%Nk = block(1,l-1)%mesh%Nk/2
      block(1,l)%mesh%gc = block(1,l-1)%mesh%gc
    enddo
  endif
    ! allocating block
  do l=1,global%mesh%Nl
    call alloc_block(global=global,block=block(1,l))
  enddo
  do l=1,global%mesh%Nl
    write(stdout,'(A,I3)',iostat=err)'   Grid level ',l
    ! mesh files
    err = read(UnitScratch(1,b,l),block(1,l)%mesh%node)
    err = read(UnitScratch(1,b,l),block(1,l)%mesh%cell)
    close(UnitScratch(1,b,l))
    err = save_bmesh(filename=file_name(basename = trim(global%file%Path_OutPut)//trim(global%file%File_Mesh), &
                                        suffix   = '.geo',                                                     &
                                        blk      = b,                                                          &
                                        grl      = l),                                                         &
                     block=block(1,l))

    ! boundary condition files
    err = read(UnitScratch(2,b,l),block(1,l)%bc%BCi)
    err = read(UnitScratch(2,b,l),block(1,l)%bc%BCj)
    err = read(UnitScratch(2,b,l),block(1,l)%bc%BCk)
    close(UnitScratch(2,b,l))
    err = save_bbc(filename=file_name(basename = trim(global%file%Path_OutPut)//trim(global%file%File_BC), &
                                      suffix   = '.bco',                                                   &
                                      blk      = b,                                                        &
                                      grl      = l),                                                       &
                   block=block(1,l))

    ! initial condition files
    read(UnitScratch(3,b,l),iostat=err)block(1,l)%fluid%Dt
    err = read(UnitScratch(3,b,l),block(1,l)%fluid%P)
    close(UnitScratch(3,b,l))
    err = save_bfluid(filename=file_name(basename = trim(global%file%Path_OutPut)//trim(global%file%File_Mesh), &
                                         suffix   = '.itc',                                                     &
                                         blk      = b,                                                          &
                                         grl      = l),                                                         &
                      global=global,                                                                            &
                      block=block(1,l))
  enddo
  write(stdout,*)
enddo
write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
write(stdout,*)
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine print_usage()
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for printing the correct use of the program.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(stderr,'(A)')' A valid file name of the options file must be provided as command line argument'
  write(stderr,'(A)')' No argument has been passed to command line'
  write(stderr,'(A)')' Correct use is:'
  write(stderr,*)
  write(stderr,'(A)')' ICG "valid_option_file_name"'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_usage

  subroutine icg_init(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for initializing the initial conditions according to the input options.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(OUT):: global      ! Global-level data.
  integer(I_P)::                   b           ! Counter.
  integer(I_P)::                   Nca = 0     ! Number of command line arguments.
  character(60)::                  File_Option ! Options file name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! parsing command line for getting global option file name
  Nca = command_argument_count ()
  if (Nca==0) then
    call print_usage
    stop
  else
    call get_command_argument (1, File_Option)
    File_Option = string_OS_sep(File_Option) ; File_Option = adjustl(trim(File_Option))
  endif
  write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
  write(stdout,'(A)') ' Loading '//trim(File_Option)
  err = load_option_file(File_Option)
  write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
  write(stdout,*)
  write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
  write(stdout,'(A)')     ' Input Files'
  write(stdout,'(A)')     '  Species file:             '//adjustl(trim(global%file%Path_InPut))//trim(global%file%File_Spec)
  select case(trim(In_type))
  case('BLOCKS')
    write(stdout,'(A)')   '  Blocks file:              '//adjustl(trim(global%file%Path_InPut))//trim(File_Blocks)
  case('ICEMCFD')
    write(stdout,'(A)')   '  ICEMCFD files:'
    do b=1,Nf_Icem
      write(stdout,'(A)') '                            '//adjustl(trim(global%file%Path_InPut))//trim(File_Icem(b))//'.geo/.topo'
    enddo
  endselect
  write(stdout,'(A)')     ' Output Files'
  write(stdout,'(A)')     '  Mesh file:                '//adjustl(trim(global%file%Path_OutPut))//trim(global%file%File_Mesh)
  write(stdout,'(A)')     '  Boundary conditions file: '//adjustl(trim(global%file%Path_OutPut))//trim(global%file%File_BC)
  write(stdout,'(A)')     '  Initial conditions file:  '//adjustl(trim(global%file%Path_OutPut))//trim(global%file%File_Init)
  write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
  write(stdout,*)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine icg_init

  function load_option_file(filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for loading option file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: filename ! Name of file where initial species are saved.
  integer(I_P)::             err      ! Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::             UnitFree ! Free logic unit.
  logical::                  is_file  ! Flag for inquiring the presence of option file.
  character(3)::             os_type  ! Type operating system.
  integer(I_P)::             f        ! Files counter.
  !--------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_option')
  UnitFree = Get_Unit()
  open(unit = UnitFree, file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*)
  read(UnitFree,*,iostat=err)!               OS
  read(UnitFree,*,iostat=err)os_type ; os_type = Upper_Case(os_type) ; OS = init_OS(c_id=os_type)
  read(UnitFree,*)!       INPUT FILES
  read(UnitFree,*,iostat=err)global%file%Path_InPut
  read(UnitFree,*,iostat=err)global%file%File_Spec
  global%file%Path_InPut = string_OS_sep(global%file%Path_InPut) ; global%file%Path_InPut = adjustl(trim(global%file%Path_InPut))
  global%file%File_Spec  = string_OS_sep(global%file%File_Spec ) ; global%file%File_Spec  = adjustl(trim(global%file%File_Spec ))
  read(UnitFree,*,iostat=err)In_type ; In_type = Upper_Case(In_type)
  select case(trim(In_type))
  case('BLOCKS')
    read(UnitFree,*,iostat=err)File_Blocks ; File_Blocks = string_OS_sep(File_Blocks) ; File_Blocks = adjustl(trim(File_Blocks))
  case('ICEMCFD')
    read(UnitFree,*,iostat=err)Nf_Icem
    if (allocated(File_Icem)) deallocate(File_Icem) ; allocate(File_Icem(1:Nf_Icem))
    do f=1,Nf_Icem
      read(UnitFree,*,iostat=err)File_Icem(f) ; File_Icem(f)=string_OS_sep(File_Icem(f)) ; File_Icem(f)=adjustl(trim(File_Icem(f)))
    enddo
  endselect
  read(UnitFree,*)!       GRID OPTIONS
  read(UnitFree,*,iostat=err)global%mesh%Nl
  read(UnitFree,*)!       OUTPUT FILES
  read(UnitFree,*,iostat=err)global%file%Path_OutPut
  read(UnitFree,*,iostat=err)global%file%File_Mesh
  read(UnitFree,*,iostat=err)global%file%File_BC
  read(UnitFree,*,iostat=err)global%file%File_Init
  global%file%Path_OutPut=string_OS_sep(global%file%Path_OutPut) ; global%file%Path_OutPut=adjustl(trim(global%file%Path_OutPut))
  global%file%File_Mesh  =string_OS_sep(global%file%File_Mesh  ) ; global%file%File_Mesh  =adjustl(trim(global%file%File_Mesh  ))
  global%file%File_BC    =string_OS_sep(global%file%File_BC    ) ; global%file%File_BC    =adjustl(trim(global%file%File_BC    ))
  global%file%File_Init  =string_OS_sep(global%file%File_Init  ) ; global%file%File_Init  =adjustl(trim(global%file%File_Init  ))
  close(UnitFree)
  ! creating the output directory
  err = make_dir(directory=global%file%Path_OutPut)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_option_file

  function load_blocks(filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for loading blocks file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: filename ! Name of file where blocks informations are saved.
  integer(I_P)::             err      ! Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::             UnitFree ! Free logic unit.
  logical::                  is_file  ! Flag for inquiring the presence of blocks file.
  character(3)::             bc_str   ! Type of boundary condition.
  real(R_P)::                Di,Dj,Dk ! Space steps in i, j, k directions.
  integer(I_P)::             b,f,l    ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_blocks')
  UnitFree = Get_Unit()
  open(unit = UnitFree, file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*)
  read(UnitFree,*)
  read(UnitFree,*) f
  if (f/=global%fluid%Ns) then
    write(stderr,'(A)')   ' The number of initial species (Ns) of the blocks file'
    write(stderr,'(A)')   ' '//adjustl(trim(filename))
    write(stderr,'(A)')   ' is different from one of the initial species file'
    write(stderr,'(A,I4)')' Blocks file          = ',f
    write(stderr,'(A,I4)')' Initial species file = ',global%fluid%Ns
    stop
  endif
  read(UnitFree,*) global%mesh%Nb_tot
  ! allocating array for blocks informations
  if (allocated(blocks))  then
    do b=lbound(blocks,dim=1),ubound(blocks,dim=1)
      err = free(blocks(b)%P)
      err = free(blocks(b)%bc)
    enddo
    deallocate(blocks)
  endif
  allocate(blocks(1:global%mesh%Nb_tot))
  call init_prim(Ns=global%fluid%Ns,prim=blocks%P)
  if (allocated(UnitScratch)) deallocate(UnitScratch) ; allocate(UnitScratch(1:3,1:global%mesh%Nb_tot,1:global%mesh%Nl))
  ! reading informations of each block
  do b=1,global%mesh%Nb_tot
    read(UnitFree,*)
    read(UnitFree,*) (blocks(b)%gc(f),f=1,6)
    read(UnitFree,*) blocks(b)%Ni,blocks(b)%Nj,blocks(b)%Nk
    read(UnitFree,*) blocks(b)%xmin,blocks(b)%ymin,blocks(b)%zmin
    read(UnitFree,*) blocks(b)%xmax,blocks(b)%ymax,blocks(b)%zmax
    do f=1,6
      read(UnitFree,*) bc_str,l
      blocks(b)%bc(f)%tp = get_bc_id(myrank,bc_str)
      if (blocks(b)%bc(f)%tp==bc_adj) then
        !blocks(b)%bc(f) = set_bc(adj=Type_Adj(b=l,i=0,j=0,k=0))
        call set_bc(adj=Type_Adj(b=l,i=0,j=0,k=0),bc=blocks(b)%bc(f))
      elseif (blocks(b)%bc(f)%tp==bc_in1.or.blocks(b)%bc(f)%tp==bc_in2) then
        !blocks(b)%bc(f) = set_bc(inf=l)
        call set_bc(inf=l,bc=blocks(b)%bc(f))
      endif
    enddo
    do f=1,global%fluid%Ns
      read(UnitFree,*) blocks(b)%P%r(f)
    enddo
    read(UnitFree,*) blocks(b)%P%v%x
    read(UnitFree,*) blocks(b)%P%v%y
    read(UnitFree,*) blocks(b)%P%v%z
    read(UnitFree,*) blocks(b)%P%p
    read(UnitFree,*) blocks(b)%P%d
    read(UnitFree,*) blocks(b)%P%g
  enddo
  close(UnitFree)
  ! for memory efficiency each block is generated alone and then stored in scratch file
  ! opening the scratch file where the current block is temporarily stored
  do l=1,global%mesh%Nl
    do b=1,global%mesh%Nb_tot
      UnitScratch(1,b,l)=Get_Unit() ; open(unit=UnitScratch(1,b,l),form='UNFORMATTED',status='SCRATCH',iostat=err) ! geo file
      UnitScratch(2,b,l)=Get_Unit() ; open(unit=UnitScratch(2,b,l),form='UNFORMATTED',status='SCRATCH',iostat=err) ! bco file
      UnitScratch(3,b,l)=Get_Unit() ; open(unit=UnitScratch(3,b,l),form='UNFORMATTED',status='SCRATCH',iostat=err) ! itc file
    enddo
  enddo
  ! the number of global blocks is artificially set to 1
  global%mesh%Nb = 1
  ! allocating block-level data
  if (allocated(block)) then
    do l=lbound(block,dim=2),ubound(block,dim=2)
      do b=lbound(block,dim=1),ubound(block,dim=1)
        call free_block(block=block(b,l))
      enddo
    enddo
    deallocate(block)
  endif
  allocate(block(1:global%mesh%Nb,1:global%mesh%Nl))
  ! generating mesh, boundary and initial conditions of each block and storing in scratch files
  do b=1,global%mesh%Nb_tot
    ! setting mesh dimensions for grid level 1
    block(1,1)%mesh%Ni = blocks(b)%Ni
    block(1,1)%mesh%Nj = blocks(b)%Nj
    block(1,1)%mesh%Nk = blocks(b)%Nk
    block(1,1)%mesh%gc = blocks(b)%gc
    ! computing mesh dimensions for other grid levels
    if (global%mesh%Nl>1) then
      do l=2,global%mesh%Nl
        ! computing the number of cells of coarser grid levels
        if     (mod(block(1,l-1)%mesh%Ni,2)/=0) then
          write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
          write(stderr,'(A,I4)')' Impossible to compute grid level ',l
          write(stderr,'(A,I4)')' Inconsistent direction i, Ni ',block(1,l-1)%mesh%Ni
          write(stderr,'(A,I4)')' level ',l-1
          write(stderr,'(A,I4)')' block ',b
          stop
        elseif (mod(block(1,l-1)%mesh%Nj,2)/=0) then
          write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
          write(stderr,'(A,I4)')' Impossible to compute grid level ',l
          write(stderr,'(A,I4)')' Inconsistent direction j, Nj ',block(1,l-1)%mesh%Nj
          write(stderr,'(A,I4)')' level ',l-1
          write(stderr,'(A,I4)')' block ',b
          stop
        elseif (mod(block(1,l-1)%mesh%Nk,2)/=0) then
          write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
          write(stderr,'(A,I4)')' Impossible to compute grid level ',l
          write(stderr,'(A,I4)')' Inconsistent direction k, Nk ',block(1,l-1)%mesh%Nk
          write(stderr,'(A,I4)')' level ',l-1
          write(stderr,'(A,I4)')' block ',b
          stop
        endif
        block(1,l)%mesh%Ni = block(1,l-1)%mesh%Ni/2
        block(1,l)%mesh%Nj = block(1,l-1)%mesh%Nj/2
        block(1,l)%mesh%Nk = block(1,l-1)%mesh%Nk/2
        block(1,l)%mesh%gc = block(1,l-1)%mesh%gc
      enddo
    endif
      ! allocating block
    do l=1,global%mesh%Nl
      call alloc_block(global=global,block=block(1,l))
    enddo
    ! node coordinates computing
    write(stdout,'(A)')'  Computing the nodes coordinates'
    Di = (blocks(b)%xmax-blocks(b)%xmin)/real(block(1,1)%mesh%Ni,R_P)
    Dj = (blocks(b)%ymax-blocks(b)%ymin)/real(block(1,1)%mesh%Nj,R_P)
    Dk = (blocks(b)%zmax-blocks(b)%zmin)/real(block(1,1)%mesh%Nk,R_P)
    do k=0-block(1,1)%mesh%gc(5),block(1,1)%mesh%Nk+block(1,1)%mesh%gc(6)
      do j=0-block(1,1)%mesh%gc(3),block(1,1)%mesh%Nj+block(1,1)%mesh%gc(4)
        do i=0-block(1,1)%mesh%gc(1),block(1,1)%mesh%Ni+block(1,1)%mesh%gc(2)
          !block(1,1)%mesh%node(i,j,k) = set_Vector(x = blocks(b)%xmin + real(i,R_P)*Di, &
          !                                         y = blocks(b)%ymin + real(j,R_P)*Dj, &
          !                                         z = blocks(b)%zmin + real(k,R_P)*Dk)
          call set_Vector(x = blocks(b)%xmin + real(i,R_P)*Di, &
                          y = blocks(b)%ymin + real(j,R_P)*Dj, &
                          z = blocks(b)%zmin + real(k,R_P)*Dk, &
                          vec=block(1,1)%mesh%node(i,j,k))
        enddo
      enddo
    enddo
    ! computing the node of other grid levels if used
    if (global%mesh%Nl>1) then
      ! inner nodes
      do l=2,global%mesh%Nl
        do k=0,block(1,l)%mesh%Nk
          do j=0,block(1,l)%mesh%Nj
            do i=0,block(1,l)%mesh%Ni
              block(1,l)%mesh%node(i,j,k) = block(1,l-1)%mesh%node(i*2,j*2,k*2)
            enddo
          enddo
        enddo
        ! ghost cells nodes
        ! left i
        do k=0-block(1,l)%mesh%gc(5),block(1,l)%mesh%Nk+block(1,l)%mesh%gc(6)
          do j=0-block(1,l)%mesh%gc(3),block(1,l)%mesh%Nj+block(1,l)%mesh%gc(4)
            do i=0-block(1,l)%mesh%gc(1),-1
              !block(1,l)%mesh%node(i,j,k) = set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
              !                                         y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
              !                                         z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk)
              call set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
                              y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
                              z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk, &
                              vec=block(1,l)%mesh%node(i,j,k))
            enddo
          enddo
        enddo
        ! right i
        do k=0-block(1,l)%mesh%gc(5),block(1,l)%mesh%Nk+block(1,l)%mesh%gc(6)
          do j=0-block(1,l)%mesh%gc(3),block(1,l)%mesh%Nj+block(1,l)%mesh%gc(4)
            do i=block(1,l)%mesh%Ni+1,block(1,l)%mesh%Ni+block(1,l)%mesh%gc(2)
              !block(1,l)%mesh%node(i,j,k) = set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
              !                                         y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
              !                                         z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk)
              call set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
                              y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
                              z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk, &
                              vec=block(1,l)%mesh%node(i,j,k))
            enddo
          enddo
        enddo
        ! left j
        do k=0-block(1,l)%mesh%gc(5),block(1,l)%mesh%Nk+block(1,l)%mesh%gc(6)
          do j=0-block(1,l)%mesh%gc(3),-1
            do i=0-block(1,l)%mesh%gc(1),block(1,l)%mesh%Ni+block(1,l)%mesh%gc(2)
              !block(1,l)%mesh%node(i,j,k) = set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
              !                                         y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
              !                                         z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk)
              call set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
                              y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
                              z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk, &
                              vec=block(1,l)%mesh%node(i,j,k))
            enddo
          enddo
        enddo
        ! right j
        do k=0-block(1,l)%mesh%gc(5),block(1,l)%mesh%Nk+block(1,l)%mesh%gc(6)
          do j=block(1,l)%mesh%Nj+1,block(1,l)%mesh%Nj+block(1,l)%mesh%gc(4)
            do i=0-block(1,l)%mesh%gc(1),block(1,l)%mesh%Ni+block(1,l)%mesh%gc(2)
              !block(1,l)%mesh%node(i,j,k) = set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
              !                                         y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
              !                                         z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk)
              call set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
                              y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
                              z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk, &
                              vec=block(1,l)%mesh%node(i,j,k))
            enddo
          enddo
        enddo
        ! left k
        do k=0-block(1,l)%mesh%gc(5),-1
          do j=0-block(1,l)%mesh%gc(3),block(1,l)%mesh%Nj+block(1,l)%mesh%gc(4)
            do i=0-block(1,l)%mesh%gc(1),block(1,l)%mesh%Ni+block(1,l)%mesh%gc(2)
              !block(1,l)%mesh%node(i,j,k) = set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
              !                                         y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
              !                                         z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk)
              call set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
                              y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
                              z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk, &
                              vec=block(1,l)%mesh%node(i,j,k))
            enddo
          enddo
        enddo
        ! right k
        do k=block(1,l)%mesh%Nk+1,block(1,l)%mesh%Nk+block(1,l)%mesh%gc(6)
          do j=0-block(1,l)%mesh%gc(3),block(1,l)%mesh%Nj+block(1,l)%mesh%gc(4)
            do i=0-block(1,l)%mesh%gc(1),block(1,l)%mesh%Ni+block(1,l)%mesh%gc(2)
              !block(1,l)%mesh%node(i,j,k) = set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
              !                                         y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
              !                                         z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk)
              call set_Vector(x = blocks(b)%xmin + real(i*2*(l-1),R_P)*Di, &
                              y = blocks(b)%ymin + real(j*2*(l-1),R_P)*Dj, &
                              z = blocks(b)%zmin + real(k*2*(l-1),R_P)*Dk, &
                              vec=block(1,l)%mesh%node(i,j,k))
            enddo
          enddo
        enddo
      enddo
    endif
    ! boundary conditions setting
    write(stdout,'(A)')'  Setting boundary conditions'
    do l=1,global%mesh%Nl
      ! left i
      do k=1,block(1,l)%mesh%Nk
        do j=1,block(1,l)%mesh%Nj
          do i=1-block(1,l)%mesh%gc(1),0
            block(1,l)%bc%BCi(i,j,k)%tp = blocks(b)%bc(1)%tp
            select case(blocks(b)%bc(1)%tp)
            case(bc_adj)
              !block(1,l)%bc%BCi(i,j,k) = set_bc(adj=Type_Adj(b = blocks(b)%bc(1)%adj%b,                           &
              !                                               i = blocks(blocks(b)%bc(1)%adj%b)%Ni/(2**(l-1)) + i, &
              !                                               j = j,                                               &
              !                                               k = k))
              call set_bc(adj=Type_Adj(b = blocks(b)%bc(1)%adj%b,                           &
                                       i = blocks(blocks(b)%bc(1)%adj%b)%Ni/(2**(l-1)) + i, &
                                       j = j,                                               &
                                       k = k),&
                          bc=block(1,l)%bc%BCi(i,j,k))
            case(bc_in1)
              !block(1,l)%bc%BCi(i,j,k) = set_bc(inf=blocks(b)%bc(1)%inf)
              call set_bc(inf=blocks(b)%bc(1)%inf,bc=block(1,l)%bc%BCi(i,j,k))
            endselect
          enddo
        enddo
      enddo
      ! right i
      do k=1,block(1,l)%mesh%Nk
        do j=1,block(1,l)%mesh%Nj
          do i=block(1,l)%mesh%Ni+1,block(1,l)%mesh%Ni+block(1,l)%mesh%gc(2)
            block(1,l)%bc%BCi(i-1,j,k)%tp = blocks(b)%bc(2)%tp
            select case(blocks(b)%bc(2)%tp)
            case(bc_adj)
              !block(1,l)%bc%BCi(i-1,j,k) = set_bc(adj=Type_Adj(b = blocks(b)%bc(2)%adj%b,             &
              !                                                 i = i - block(1,l)%mesh%Ni/(2**(l-1)), &
              !                                                 j = j,                                 &
              !                                                 k = k))
              call set_bc(adj=Type_Adj(b = blocks(b)%bc(2)%adj%b,             &
                                       i = i - block(1,l)%mesh%Ni/(2**(l-1)), &
                                       j = j,                                 &
                                       k = k), &
                          bc=block(1,l)%bc%BCi(i-1,j,k))
            case(bc_in1)
              !block(1,l)%bc%BCi(i-1,j,k) = set_bc(inf=blocks(b)%bc(2)%inf)
              call set_bc(inf=blocks(b)%bc(2)%inf,bc=block(1,l)%bc%BCi(i-1,j,k))
            endselect
          enddo
        enddo
      enddo
      ! left j
      do k=1,block(1,l)%mesh%Nk
        do i=1,block(1,l)%mesh%Ni
          do j=1-block(1,l)%mesh%gc(3),0
            block(1,l)%bc%BCj(i,j,k)%tp = blocks(b)%bc(3)%tp
            select case(blocks(b)%bc(3)%tp)
            case(bc_adj)
              !block(1,l)%bc%BCj(i,j,k) = set_bc(adj=Type_Adj(b = blocks(b)%bc(3)%adj%b, &
              !                                               i = i, &
              !                                               j = blocks(blocks(b)%bc(3)%adj%b)%Nj/(2**(l-1)) + j, &
              !                                               k = k))
              call set_bc(adj=Type_Adj(b = blocks(b)%bc(3)%adj%b, &
                                                             i = i, &
                                                             j = blocks(blocks(b)%bc(3)%adj%b)%Nj/(2**(l-1)) + j, &
                                                             k = k), &
                          bc=block(1,l)%bc%BCj(i,j,k))
            case(bc_in1)
              !block(1,l)%bc%BCj(i,j,k) = set_bc(inf=blocks(b)%bc(3)%inf)
              call set_bc(inf=blocks(b)%bc(3)%inf,bc=block(1,l)%bc%BCj(i,j,k))
            endselect
          enddo
        enddo
      enddo
      ! right j
      do k=1,block(1,l)%mesh%Nk
        do i=1,block(1,l)%mesh%Ni
          do j=block(1,l)%mesh%Nj+1,block(1,l)%mesh%Nj+block(1,l)%mesh%gc(4)
            block(1,l)%bc%BCj(i,j-1,k)%tp = blocks(b)%bc(4)%tp
            select case(blocks(b)%bc(4)%tp)
            case(bc_adj)
              !block(1,l)%bc%BCj(i,j-1,k) = set_bc(adj=Type_Adj(b = blocks(b)%bc(4)%adj%b,             &
              !                                                 i = i,                                 &
              !                                                 j = j - block(1,l)%mesh%Nj/(2**(l-1)), &
              !                                                 k = k))
              call set_bc(adj=Type_Adj(b = blocks(b)%bc(4)%adj%b,             &
                                                               i = i,                                 &
                                                               j = j - block(1,l)%mesh%Nj/(2**(l-1)), &
                                                               k = k), &
                          bc=block(1,l)%bc%BCj(i,j-1,k))
            case(bc_in1)
              !block(1,l)%bc%BCj(i,j-1,k) = set_bc(inf=blocks(b)%bc(4)%inf)
              call set_bc(inf=blocks(b)%bc(4)%inf,bc=block(1,l)%bc%BCj(i,j-1,k))
            endselect
          enddo
        enddo
      enddo
      ! left k
      do j=1,block(1,l)%mesh%Nj
        do i=1,block(1,l)%mesh%Ni
          do k=1-block(1,l)%mesh%gc(5),0
            block(1,l)%bc%BCk(i,j,k)%tp = blocks(b)%bc(5)%tp
            select case(blocks(b)%bc(5)%tp)
            case(bc_adj)
              !block(1,l)%bc%BCk(i,j,k) = set_bc(adj=Type_Adj(b = blocks(b)%bc(5)%adj%b, &
              !                                               i = i,                     &
              !                                               j = j,                     &
              !                                               k = blocks(blocks(b)%bc(5)%adj%b)%Nk/(2**(l-1)) + k))
              call set_bc(adj=Type_Adj(b = blocks(b)%bc(5)%adj%b, &
                                                             i = i,                     &
                                                             j = j,                     &
                                                             k = blocks(blocks(b)%bc(5)%adj%b)%Nk/(2**(l-1)) + k), &
                          bc=block(1,l)%bc%BCk(i,j,k))
            case(bc_in1)
              !block(1,l)%bc%BCk(i,j,k) = set_bc(inf=blocks(b)%bc(5)%inf)
              call set_bc(inf=blocks(b)%bc(5)%inf,bc=block(1,l)%bc%BCk(i,j,k))
            endselect
          enddo
        enddo
      enddo
      ! right k
      do j=1,block(1,l)%mesh%Nj
        do i=1,block(1,l)%mesh%Ni
          do k=block(1,l)%mesh%Nk+1,block(1,l)%mesh%Nk+block(1,l)%mesh%gc(6)
            block(1,l)%bc%BCk(i,j,k-1)%tp = blocks(b)%bc(6)%tp
            select case(blocks(b)%bc(6)%tp)
            case(bc_adj)
              !block(1,l)%bc%BCk(i,j,k-1) = set_bc(adj=Type_Adj(b = blocks(b)%bc(6)%adj%b            , &
              !                                                 i = i                                , &
              !                                                 j = j                                , &
              !                                                 k = k - block(1,l)%mesh%Nk/(2**(l-1))))
              call set_bc(adj=Type_Adj(b = blocks(b)%bc(6)%adj%b            , &
                                                               i = i                                , &
                                                               j = j                                , &
                                                               k = k - block(1,l)%mesh%Nk/(2**(l-1))), &
                          bc=block(1,l)%bc%BCk(i,j,k-1))
            case(bc_in1)
              !block(1,l)%bc%BCk(i,j,k-1) = set_bc(inf=blocks(b)%bc(6)%inf)
              call set_bc(inf=blocks(b)%bc(6)%inf,bc=block(1,l)%bc%BCk(i,j,k-1))
            endselect
          enddo
        enddo
      enddo
    enddo
    ! initial conditions setting
    write(stdout,'(A)')'  Setting initial conditions'//str(.true.,b)
    do l=1,global%mesh%Nl
      !call fluid_set(l = l, b = 1, P = blocks(b)%P)
      block(1,l)%fluid%P = blocks(b)%P
    enddo

    ! storing the mesh, boundary and initial conditions in the scratch files
    do l=1,global%mesh%Nl
      ! mesh data
      err = write(UnitScratch(1,b,l),block(1,l)%mesh%node)
      err = write(UnitScratch(1,b,l),block(1,l)%mesh%cell)
      ! boundary conditions data
      err = write(UnitScratch(2,b,l),block(1,l)%bc%BCi)
      err = write(UnitScratch(2,b,l),block(1,l)%bc%BCj)
      err = write(UnitScratch(2,b,l),block(1,l)%bc%BCk)
      ! initial conditions data
      write(UnitScratch(3,b,l),iostat=err)block(1,l)%fluid%Dt
      err = write(UnitScratch(3,b,l),block(1,l)%fluid%P)
      ! rewinding scratch files
      rewind(UnitScratch(1,b,l))
      rewind(UnitScratch(2,b,l))
      rewind(UnitScratch(3,b,l))
    enddo
  enddo
  global%mesh%Nb = global%mesh%Nb_tot
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_blocks

  function load_icemcfd(Nf,filenames) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for loading icemcfd fils.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: Nf               ! Number of input files.
  character(*), intent(IN):: filenames(1:Nf)  ! Base name of icemcfd files.
  integer(I_P)::             err              ! Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::             UnitFree         ! Free logic unit.
  integer(I_P)::             Unit_itc         ! Free logic unit for initial conditions files.
  integer(I_P)::             Unit_gc          ! Free logic unit for ghost cells definition.
  logical::                  is_file          ! Flag for inquiring the presence of file.
  integer(I_P)::             b,f,v,b1,b2,l    ! Counters.
  integer(I_P)::             bc               ! Boundary conditions counters.
  integer(I_P)::             Nb_l             ! Number of local blocks.
  character(500)::           line,line1,line2 ! Dummy strings for parsing stuff.
  integer(I_P), parameter::  No = 48          ! Number of possible orientations.
  integer(I_P)::             o,o1,o2          ! Orientation counters.
  integer(I_P)::             d1_min(1:3),d1_max(1:3)
  integer(I_P)::             d2_min(1:3),d2_max(1:3)
  integer(I_P)::             i1_min,i1_max,j1_min,j1_max,k1_min,k1_max
  integer(I_P)::             i2_min,i2_max,j2_min,j2_max,k2_min,k2_max
  integer(I_P)::             i1_sgn,j1_sgn,k1_sgn
  integer(I_P)::             i2_sgn,j2_sgn,k2_sgn
  character(6), parameter::  orientation(1:No) = (/'-i j k', '-i j-k', '-i k j', '-i k-j', '-i-j k', '-i-j-k', & ! Orientation list.
                                                   '-i-k j', '-i-k-j', '-j i k', '-j i-k', '-j k i', '-j k-i', &
                                                   '-j-i k', '-j-i-k', '-j-k i', '-j-k-i', '-k i j', '-k i-j', &
                                                   '-k j i', '-k j-i', '-k-i j', '-k-i-j', '-k-j i', '-k-j-i', &
                                                   ' i j k', ' i j-k', ' i k j', ' i k-j', ' i-j k', ' i-j-k', &
                                                   ' i-k j', ' i-k-j', ' j i k', ' j i-k', ' j k i', ' j k-i', &
                                                   ' j-i k', ' j-i-k', ' j-k i', ' j-k-i', ' k i j', ' k i-j', &
                                                   ' k j i', ' k j-i', ' k-i j', ' k-i-j', ' k-j i', ' k-j-i'/)
  character(1), parameter::  tab = char(9)
  integer(I_P)::             blk_map(1:Nf)    ! Map of blocks over icem files.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! verifying the presence of input files
  do f=1,Nf
    inquire(file=adjustl(trim(filenames(f)))//'.geo',exist=is_file,iostat=err)
    if (.NOT.is_file) call File_Not_Found(myrank,adjustl(trim(filenames(f)))//'.geo','load_icemcfd')
    inquire(file=adjustl(trim(filenames(f)))//'.topo',exist=is_file,iostat=err)
    if (.NOT.is_file) call File_Not_Found(myrank,adjustl(trim(filenames(f)))//'.topo','load_blocks')
  enddo
  ! computing the number of blocks
  global%mesh%Nb_tot = 0
  do f=1,Nf
    ! opening geometry file
    UnitFree = Get_Unit()
    open(unit = UnitFree, file = adjustl(trim(filenames(f)))//'.geo', status = 'OLD', action = 'READ', form = 'FORMATTED')
    blk_map(f)=0
    do
      read(UnitFree,'(A)',iostat=err) line
      if (err /= 0) exit
      if (index(line,'domain')>0) blk_map(f) = blk_map(f) + 1
    enddo
    close(UnitFree)                                      ! close file
    global%mesh%Nb_tot = global%mesh%Nb_tot + blk_map(f)
  enddo
  ! allocating array for blocks informations
  if (allocated(blocks))  then
    do b=lbound(blocks,dim=1),ubound(blocks,dim=1)
      err = free(blocks(b)%P)
      err = free(blocks(b)%bc)
    enddo
    deallocate(blocks)
  endif
  allocate(blocks(1:global%mesh%Nb_tot))
  call init_prim(Ns=global%fluid%Ns,prim=blocks%P)
  if (allocated(UnitScratch)) deallocate(UnitScratch) ; allocate(UnitScratch(0:3,1:global%mesh%Nb_tot,1:global%mesh%Nl))
  ! reading the number of cells for each blocks contained into the file
  b = 0
  do f=1,Nf
    ! opening geometry file
    UnitFree = Get_Unit()
    open(unit = UnitFree, file = adjustl(trim(filenames(f)))//'.geo', status = 'OLD', action = 'READ', form = 'FORMATTED')
    do
      read(UnitFree,'(A)',iostat=err) line
      if (err /= 0) exit
      if (index(line,'domain')>0)  then
        b = b + 1 ! updating the block counter
        read(line(index(line,' ')+1:),*)blocks(b)%Ni,blocks(b)%Nj,blocks(b)%Nk
        blocks(b)%Ni = blocks(b)%Ni - 1_I_P
        blocks(b)%Nj = blocks(b)%Nj - 1_I_P
        blocks(b)%Nk = blocks(b)%Nk - 1_I_P
        Unit_gc = Get_Unit()
        open(unit = Unit_gc, file = 'BLK'//trim(str(.true.,b))//'.gc', status = 'OLD', action = 'READ', form = 'FORMATTED')
        read(Unit_gc,*)(blocks(b)%gc(i),i=1,6)
        close(Unit_gc)
        write(stdout,*)
        write(stdout,'(A,'//FI_P//')') ' Ghost cells of block: ',b
        write(stdout,'(A,'//FI_P//')') ' Ghost cells:          ',blocks(b)%gc
        write(stdout,*)
      endif
    enddo
    close(UnitFree)                                      ! close file
  enddo
  ! for efficiency each block is generated alone and then stored in scratch file
  ! opening the scratch file where the current block is temporarily stored
  do l=1,global%mesh%Nl
    do b=1,global%mesh%Nb_tot
      if (l==1) then
        UnitScratch(0,b,l)=Get_Unit() ; open(unit=UnitScratch(0,b,l),form=  'FORMATTED',status='SCRATCH',iostat=err) ! topo file
      endif
        UnitScratch(1,b,l)=Get_Unit() ; open(unit=UnitScratch(1,b,l),form='UNFORMATTED',status='SCRATCH',iostat=err) ! geo file
        UnitScratch(2,b,l)=Get_Unit() ; open(unit=UnitScratch(2,b,l),form='UNFORMATTED',status='SCRATCH',iostat=err) ! bco file
        UnitScratch(3,b,l)=Get_Unit() ; open(unit=UnitScratch(3,b,l),form='UNFORMATTED',status='SCRATCH',iostat=err) ! itc file
    enddo
  enddo
  ! the number of global blocks is artificially set to 1
  global%mesh%Nb = 1
  ! allocating block-level data
  if (allocated(block)) then
    do l=lbound(block,dim=2),ubound(block,dim=2)
      do b=lbound(block,dim=1),ubound(block,dim=1)
        call free_block(block=block(b,l))
      enddo
    enddo
    deallocate(block)
  endif
  allocate(block(1:global%mesh%Nb,1:global%mesh%Nl))
  ! reading mesh, boundary and initial conditions of each block and storing in scratch files
  write(stdout,'(A)')'  Reading nodes coordinates from icemcfd files'
  b = 0
  do f=1,Nf
    ! opening geometry file
    UnitFree = Get_Unit()
    open(unit = UnitFree, file = adjustl(trim(filenames(f)))//'.geo', status = 'OLD', action = 'READ', form = 'FORMATTED')
    do
      read(UnitFree,'(A)',iostat=err) line
      if (err /= 0) exit
      if (index(line,'domain')>0)  then
        b = b + 1 ! updating the block counter
        ! setting mesh dimensions for grid level 1
        block(1,1)%mesh%Ni = blocks(b)%Ni
        block(1,1)%mesh%Nj = blocks(b)%Nj
        block(1,1)%mesh%Nk = blocks(b)%Nk
        block(1,1)%mesh%gc = blocks(b)%gc
        ! computing mesh dimensions for other grid levels
        if (global%mesh%Nl>1) then
          do l=2,global%mesh%Nl
            ! computing the number of cells of coarser grid levels
            if     (mod(block(1,l-1)%mesh%Ni,2)/=0) then
              write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
              write(stderr,'(A,I4)')' Impossible to compute grid level ',l
              write(stderr,'(A,I4)')' Inconsistent direction i, Ni ',block(1,l-1)%mesh%Ni
              write(stderr,'(A,I4)')' level ',l-1
              write(stderr,'(A,I4)')' block ',b
              stop
            elseif (mod(block(1,l-1)%mesh%Nj,2)/=0) then
              write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
              write(stderr,'(A,I4)')' Impossible to compute grid level ',l
              write(stderr,'(A,I4)')' Inconsistent direction j, Nj ',block(1,l-1)%mesh%Nj
              write(stderr,'(A,I4)')' level ',l-1
              write(stderr,'(A,I4)')' block ',b
              stop
            elseif (mod(block(1,l-1)%mesh%Nk,2)/=0) then
              write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
              write(stderr,'(A,I4)')' Impossible to compute grid level ',l
              write(stderr,'(A,I4)')' Inconsistent direction k, Nk ',block(1,l-1)%mesh%Nk
              write(stderr,'(A,I4)')' level ',l-1
              write(stderr,'(A,I4)')' block ',b
              stop
            endif
            block(1,l)%mesh%Ni = block(1,l-1)%mesh%Ni/2
            block(1,l)%mesh%Nj = block(1,l-1)%mesh%Nj/2
            block(1,l)%mesh%Nk = block(1,l-1)%mesh%Nk/2
            block(1,l)%mesh%gc = block(1,l-1)%mesh%gc
          enddo
        endif
          ! allocating block
        do l=1,global%mesh%Nl
          call alloc_block(global=global,block=block(1,l))
        enddo
        ! reading node coordinates
        do k=0,block(1,1)%mesh%Nk
          do j=0,block(1,1)%mesh%Nj
            do i=0,block(1,1)%mesh%Ni
              read(UnitFree,*)block(1,1)%mesh%node(i,j,k)%x,block(1,1)%mesh%node(i,j,k)%y,block(1,1)%mesh%node(i,j,k)%z
            enddo
          enddo
        enddo
        ! computing the node of other grid levels if used
        if (global%mesh%Nl>1) then
          do l=2,global%mesh%Nl
            do k=0,block(1,l)%mesh%Nk
              do j=0,block(1,l)%mesh%Nj
                do i=0,block(1,l)%mesh%Ni
                  block(1,l)%mesh%node(i,j,k) = block(1,l-1)%mesh%node(i*2,j*2,k*2)
                enddo
              enddo
            enddo
          enddo
        endif
        ! storing the mesh in the scratch files
        do l=1,global%mesh%Nl
          err = write(UnitScratch(1,b,l),block(1,l)%mesh%node)
          err = write(UnitScratch(1,b,l),block(1,l)%mesh%cell)
          ! rewinding scratch file
          rewind(UnitScratch(1,b,l))
        enddo
      endif
    enddo
    close(UnitFree)                                      ! close file
  enddo
  write(stdout,'(A)')'  Reading boundary and initial conditions from icemcfd files'
  ! the topological files are pre-processed and arranged one for block in the scratch files
  b = 0
  do f=1,Nf
    Nb_l = 0 ; if (f>1) Nb_l = sum(blk_map(1:f-1))
    ! opening geometry file
    UnitFree = Get_Unit()
    open(unit = UnitFree, file = adjustl(trim(filenames(f)))//'.topo', status = 'OLD', action = 'READ', form = 'FORMATTED')
    do
      read(UnitFree,'(A)',iostat=err) line
      if (err /= 0) exit
      if (index(line,'# Connectivity for domain.')>0)  then
        ! found a block connectivity information
        read(line(index(line,'# Connectivity for domain.')+26:),*)b ; b = b + Nb_l ! current block
        write(UnitScratch(0,b,1),'(A)')'# Connectivity for domain.'//trim(str(.true.,b))
        do
          read(UnitFree,'(A)',iostat=err) line1
          if (err /= 0) exit
          if (index(line1,tab//'f')>0) then
            read(UnitFree,'(A)',iostat=err) line2
            if (err /= 0) exit
            read(line1(index(line1,'domain.')+7:),*)b1 ; b1 = b1 + Nb_l ! reading the number of current block
            read(line2(index(line2,'domain.')+7:),*)b2 ; b2 = b2 + Nb_l ! reading the number of connected block
            write(UnitScratch(0,b,1),'(A)')line1(:index(line1,' domain.')+7)//trim(str(.true.,b1))//trim(line1(index(line1,tab)-1:))
            write(UnitScratch(0,b,1),'(A)')line2(:index(line2,' domain.')+7)//trim(str(.true.,b2))//trim(line2(index(line2,tab)-1:))
          elseif (index(line1,tab//'b')>0) then
            cycle
          elseif (index(line1,tab//'e')>0) then
            cycle
          elseif (index(line1,tab//'v')>0) then
            cycle
          else
            write(UnitScratch(0,b,1),*)
            exit
          endif
        enddo
      elseif (index(line,'# Boundary conditions and/or properties for domain.')>0)  then
        ! found a block boundary or initial conditions information
        read(line(index(line,'# Boundary conditions and/or properties for domain.')+51:),*)b ; b = b + Nb_l ! current block
        write(UnitScratch(0,b,1),'(A)')'# Boundary conditions and/or properties for domain.'//trim(str(.true.,b))
        do
          read(UnitFree,'(A)',iostat=err) line1
          if (err /= 0) exit
          if (index(line1,tab//'f')>0) then
            write(UnitScratch(0,b,1),'(A)')trim(line1)
          elseif (index(line1,tab//'b')>0) then
            write(UnitScratch(0,b,1),'(A)')trim(line1)
          elseif (index(line1,tab//'e')>0) then
            cycle
          elseif (index(line1,tab//'v')>0) then
            cycle
          else
            write(UnitScratch(0,b,1),*)
            exit
          endif
        enddo
      endif
    enddo
    close(UnitFree)
  enddo
  do b=1,global%mesh%Nb_tot
    rewind(UnitScratch(0,b,1))
  enddo
  ! reading boundary and initial conditions from scratch topological files
  do b=1,global%mesh%Nb_tot
    ! setting mesh dimensions for grid level 1
    block(1,1)%mesh%Ni = blocks(b)%Ni
    block(1,1)%mesh%Nj = blocks(b)%Nj
    block(1,1)%mesh%Nk = blocks(b)%Nk
    block(1,1)%mesh%gc = blocks(b)%gc
    ! computing mesh dimensions for other grid levels
    if (global%mesh%Nl>1) then
      do l=2,global%mesh%Nl
        ! computing the number of cells of coarser grid levels
        if     (mod(block(1,l-1)%mesh%Ni,2)/=0) then
          write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
          write(stderr,'(A,I4)')' Impossible to compute grid level ',l
          write(stderr,'(A,I4)')' Inconsistent direction i, Ni ',block(1,l-1)%mesh%Ni
          write(stderr,'(A,I4)')' level ',l-1
          write(stderr,'(A,I4)')' block ',b
          stop
        elseif (mod(block(1,l-1)%mesh%Nj,2)/=0) then
          write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
          write(stderr,'(A,I4)')' Impossible to compute grid level ',l
          write(stderr,'(A,I4)')' Inconsistent direction j, Nj ',block(1,l-1)%mesh%Nj
          write(stderr,'(A,I4)')' level ',l-1
          write(stderr,'(A,I4)')' block ',b
          stop
        elseif (mod(block(1,l-1)%mesh%Nk,2)/=0) then
          write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
          write(stderr,'(A,I4)')' Impossible to compute grid level ',l
          write(stderr,'(A,I4)')' Inconsistent direction k, Nk ',block(1,l-1)%mesh%Nk
          write(stderr,'(A,I4)')' level ',l-1
          write(stderr,'(A,I4)')' block ',b
          stop
        endif
        block(1,l)%mesh%Ni = block(1,l-1)%mesh%Ni/2
        block(1,l)%mesh%Nj = block(1,l-1)%mesh%Nj/2
        block(1,l)%mesh%Nk = block(1,l-1)%mesh%Nk/2
        block(1,l)%mesh%gc = block(1,l-1)%mesh%gc
      enddo
    endif
      ! allocating block
    do l=1,global%mesh%Nl
      call alloc_block(global=global,block=block(1,l))
    enddo
    do
      read(UnitScratch(0,b,1),'(A)',iostat=err) line
      if (err /= 0) exit
      if (index(line,'# Connectivity for domain.')>0)  then
        ! found a block connectivity information
        do
          read(UnitScratch(0,b,1),'(A)',iostat=err) line1
          if (err /= 0) exit
          if (index(line1,tab//'f')>0) then
            read(UnitScratch(0,b,1),'(A)',iostat=err) line2
            if (err /= 0) exit
            read(line1(index(line1,'domain.')+7:),*)b1 ! reading the number of current block
            read(line2(index(line2,'domain.')+7:),*)b2 ! reading the number of connected block
            read(line1(index(line1,tab//'f',BACK=.true.)+2:),*)d1_min(1),d1_min(2),d1_min(3),d1_max(1),d1_max(2),d1_max(3)
            read(line2(index(line2,tab//'f',BACK=.true.)+2:),*)d2_min(1),d2_min(2),d2_min(3),d2_max(1),d2_max(2),d2_max(3)
            ! reading the orientation of the 2 blocks
            do o=1,No
              if (index(line1,orientation(o))>0) o1 = o
              if (index(line2,orientation(o))>0) o2 = o
            enddo
            ! mapping the block orientation in i,j,k direction
            ! block 1
            select case(orientation(o1)(2:2))
            case('i')
              if (orientation(o1)(1:1)==' ') then
                i1_min = d1_min(1) - 1 ; i1_max = d1_max(1) - 1
                i1_sgn =  1
              else
                i1_min = d1_max(1) - 1 ; i1_max = d1_min(1) - 1
                i1_sgn = -1
              endif
            case('j')
              if (orientation(o1)(1:1)==' ') then
                j1_min = d1_min(1) - 1 ; j1_max = d1_max(1) - 1
                j1_sgn =  1
              else
                j1_min = d1_max(1) - 1 ; j1_max = d1_min(1) - 1
                j1_sgn = -1
              endif
            case('k')
              if (orientation(o1)(1:1)==' ') then
                k1_min = d1_min(1) - 1 ; k1_max = d1_max(1) - 1
                k1_sgn =  1
              else
                k1_min = d1_max(1) - 1 ; k1_max = d1_max(1) - 1
                k1_sgn = -1
              endif
            endselect
            select case(orientation(o1)(4:4))
            case('i')
              if (orientation(o1)(3:3)==' ') then
                i1_min = d1_min(2) - 1 ; i1_max = d1_max(2) - 1
                i1_sgn =  1
              else
                i1_min = d1_max(2) - 1 ; i1_max = d1_min(2) - 1
                i1_sgn = -1
              endif
            case('j')
              if (orientation(o1)(3:3)==' ') then
                j1_min = d1_min(2) - 1 ; j1_max = d1_max(2) - 1
                j1_sgn =  1
              else
                j1_min = d1_max(2) - 1 ; j1_max = d1_min(2) - 1
                j1_sgn = -1
              endif
            case('k')
              if (orientation(o1)(3:3)==' ') then
                k1_min = d1_min(2) - 1 ; k1_max = d1_max(2) - 1
                k1_sgn =  1
              else
                k1_min = d1_max(2) - 1 ; k1_max = d1_min(2) - 1
                k1_sgn = -1
              endif
            endselect
            select case(orientation(o1)(6:6))
            case('i')
              if (orientation(o1)(5:5)==' ') then
                i1_min = d1_min(3) - 1 ; i1_max = d1_max(3) - 1
                i1_sgn =  1
              else
                i1_min = d1_max(3) - 1 ; i1_max = d1_min(3) - 1
                i1_sgn = -1
              endif
            case('j')
              if (orientation(o1)(5:5)==' ') then
                j1_min = d1_min(3) - 1 ; j1_max = d1_max(3) - 1
                j1_sgn =  1
              else
                j1_min = d1_max(3) - 1 ; j1_max = d1_min(3) - 1
                j1_sgn = -1
              endif
            case('k')
              if (orientation(o1)(5:5)==' ') then
                k1_min = d1_min(3) - 1 ; k1_max = d1_max(3) - 1
                k1_sgn =  1
              else
                k1_min = d1_max(3) - 1 ; k1_max = d1_min(3) - 1
                k1_sgn = -1
              endif
            endselect
            ! block 2
            select case(orientation(o2)(2:2))
            case('i')
              if (orientation(o2)(1:1)==' ') then
                i2_min = d2_min(1) - 1 ; i2_max = d2_max(1) - 1
                i2_sgn =  1
              else
                i2_min = d2_max(1) - 1 ; i2_max = d2_min(1) - 1
                i2_sgn = -1
              endif
            case('j')
              if (orientation(o2)(1:1)==' ') then
                j2_min = d2_min(1) - 1 ; j2_max = d2_max(1) - 1
                j2_sgn =  1
              else
                j2_min = d2_max(1) - 1 ; j2_max = d2_min(1) - 1
                j2_sgn = -1
              endif
            case('k')
              if (orientation(o2)(1:1)==' ') then
                k2_min = d2_min(1) - 1 ; k2_max = d2_max(1) - 1
                k2_sgn =  1
              else
                k2_min = d2_max(1) - 1 ; k2_max = d2_min(1) - 1
                k2_sgn = -1
              endif
            endselect
            select case(orientation(o2)(4:4))
            case('i')
              if (orientation(o2)(3:3)==' ') then
                i2_min = d2_min(2) - 1 ; i2_max = d2_max(2) - 1
                i2_sgn =  1
              else
                i2_min = d2_max(2) - 1 ; i2_max = d2_min(2) - 1
                i2_sgn = -1
              endif
            case('j')
              if (orientation(o2)(3:3)==' ') then
                j2_min = d2_min(2) - 1 ; j2_max = d2_max(2) - 1
                j2_sgn =  1
              else
                j2_min = d2_max(2) - 1 ; j2_max = d2_min(2) - 1
                j2_sgn = -1
              endif
            case('k')
              if (orientation(o2)(3:3)==' ') then
                k2_min = d2_min(2) - 1 ; k2_max = d2_max(2) - 1
                k2_sgn =  1
              else
                k2_min = d2_max(2) - 1 ; k2_max = d2_min(2) - 1
                k2_sgn = -1
              endif
            endselect
            select case(orientation(o2)(6:6))
            case('i')
              if (orientation(o2)(5:5)==' ') then
                i2_min = d2_min(3) - 1 ; i2_max = d2_max(3) - 1
                i2_sgn =  1
              else
                i2_min = d2_max(3) - 1 ; i2_max = d2_min(3) - 1
                i2_sgn = -1
              endif
            case('j')
              if (orientation(o2)(5:5)==' ') then
                j2_min = d2_min(3) - 1 ; j2_max = d2_max(3) - 1
                j2_sgn =  1
              else
                j2_min = d2_max(3) - 1 ; j2_max = d2_min(3) - 1
                j2_sgn = -1
              endif
            case('k')
              if (orientation(o2)(5:5)==' ') then
                k2_min = d2_min(3) - 1 ; k2_max = d2_max(3) - 1
                k2_sgn =  1
              else
                k2_min = d2_max(3) - 1 ; k2_max = d2_min(3) - 1
                k2_sgn = -1
              endif
            endselect
            ! imposing the adjacent boundary conditions for all grid levels
            do l=1,global%mesh%Nl
              select case(orientation(o1)(6:6))
              case('i')
                ! the adjacent face is i-face for block 1
                if (i1_min==0) then
                  i1_min = 1-block(1,l)%mesh%gc(1)
                  i1_max = 0
                else
                  i1_min = blocks(b1)%Ni/(2**(l-1))
                  i1_max = blocks(b1)%Ni/(2**(l-1))+block(1,l)%mesh%gc(2)-1
                endif
                do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
                  do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
                    do i=i1_min,i1_max
                      block(1,l)%bc%BCi(i,j,k)%tp    = bc_adj
                      if (allocated(block(1,l)%bc%BCi(i,j,k)%adj)) deallocate(block(1,l)%bc%BCi(i,j,k)%adj)
                      allocate(block(1,l)%bc%BCi(i,j,k)%adj)
                      block(1,l)%bc%BCi(i,j,k)%adj%b = b2
                      do o=2,4,2
                        select case(orientation(o1)(o:o))
                        case('j')
                          select case(orientation(o2)(o:o))
                          case('i')
                            if (j1_sgn*i2_sgn<0) then
                              block(1,l)%bc%BCi(i,j,k)%adj%i = j1_max/2**(l-1) - j
                            else
                              block(1,l)%bc%BCi(i,j,k)%adj%i = j
                            endif
                          case('j')
                            if (j1_sgn*j2_sgn<0) then
                              block(1,l)%bc%BCi(i,j,k)%adj%j = j1_max/2**(l-1) - j
                            else
                              block(1,l)%bc%BCi(i,j,k)%adj%j = j
                            endif
                          case('k')
                            if (j1_sgn*k2_sgn<0) then
                              block(1,l)%bc%BCi(i,j,k)%adj%k = j1_max/2**(l-1) - j
                            else
                              block(1,l)%bc%BCi(i,j,k)%adj%k = j
                            endif
                          endselect
                        case('k')
                          select case(orientation(o2)(o:o))
                          case('i')
                            if (k1_sgn*i2_sgn<0) then
                              block(1,l)%bc%BCi(i,j,k)%adj%i = k1_max/2**(l-1) - k
                            else
                              block(1,l)%bc%BCi(i,j,k)%adj%i = k
                            endif
                          case('j')
                            if (k1_sgn*j2_sgn<0) then
                              block(1,l)%bc%BCi(i,j,k)%adj%j = k1_max/2**(l-1) - k
                            else
                              block(1,l)%bc%BCi(i,j,k)%adj%j = k
                            endif
                          case('k')
                            if (k1_sgn*k2_sgn<0) then
                              block(1,l)%bc%BCi(i,j,k)%adj%k = k1_max/2**(l-1) - k
                            else
                              block(1,l)%bc%BCi(i,j,k)%adj%k = k
                            endif
                          endselect
                        endselect
                      enddo
                      select case(orientation(o2)(6:6))
                      case('i')
                        if (i1_max==0) then
                          block(1,l)%bc%BCi(i,j,k)%adj%i = (blocks(b2)%Ni/(2**(l-1))+i)*(1-i1_sgn*i2_sgn)/2 + &
                                                           (-i+1                      )*(1+i1_sgn*i2_sgn)/2
                        else
                          block(1,l)%bc%BCi(i,j,k)%adj%i =                                                  &
                            (i-blocks(b1)%Ni/(2**(l-1))+1                           )*(1-i1_sgn*i2_sgn)/2 + &
                            (  blocks(b2)%Ni/(2**(l-1))-(i-blocks(b1)%Ni/(2**(l-1))))*(1+i1_sgn*i2_sgn)/2
                        endif
                      case('j')
                        if (i1_max==0) then
                          block(1,l)%bc%BCi(i,j,k)%adj%j = (blocks(b2)%Nj/(2**(l-1))+i)*(1-i1_sgn*j2_sgn)/2 + &
                                                           (-i+1                      )*(1+i1_sgn*j2_sgn)/2
                        else
                          block(1,l)%bc%BCi(i,j,k)%adj%j =                                                  &
                            (i-blocks(b1)%Ni/(2**(l-1))+1                           )*(1-i1_sgn*j2_sgn)/2 + &
                            (  blocks(b2)%Nj/(2**(l-1))-(i-blocks(b1)%Ni/(2**(l-1))))*(1+i1_sgn*j2_sgn)/2
                        endif
                      case('k')
                        if (i1_max==0) then
                          block(1,l)%bc%BCi(i,j,k)%adj%k = (blocks(b2)%Nk/(2**(l-1))+i)*(1-i1_sgn*k2_sgn)/2 + &
                                                           (-i+1                      )*(1+i1_sgn*k2_sgn)/2
                        else
                          block(1,l)%bc%BCi(i,j,k)%adj%k =                                                  &
                            (i-blocks(b1)%Ni/(2**(l-1))+1                           )*(1-i1_sgn*k2_sgn)/2 + &
                            (  blocks(b2)%Nk/(2**(l-1))-(i-blocks(b1)%Ni/(2**(l-1))))*(1+i1_sgn*k2_sgn)/2
                        endif
                      endselect
                    enddo
                  enddo
                enddo
              case('j')
                ! the adjacent face is j-face for block 1
                if (j1_min==0) then
                  j1_min = 1-block(1,l)%mesh%gc(3)
                  j1_max = 0
                else
                  j1_min = blocks(b1)%Nj/(2**(l-1))
                  j1_max = blocks(b1)%Nj/(2**(l-1))+block(1,l)%mesh%gc(4)-1
                endif
                do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
                  do j=j1_min,j1_max
                    do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
                      block(1,l)%bc%BCj(i,j,k)%tp    = bc_adj
                      if (allocated(block(1,l)%bc%BCj(i,j,k)%adj)) deallocate(block(1,l)%bc%BCj(i,j,k)%adj)
                      allocate(block(1,l)%bc%BCj(i,j,k)%adj)
                      block(1,l)%bc%BCj(i,j,k)%adj%b = b2
                      do o=2,4,2
                        select case(orientation(o1)(o:o))
                        case('i')
                          select case(orientation(o2)(o:o))
                          case('i')
                            if (i1_sgn*i2_sgn<0) then
                              block(1,l)%bc%BCj(i,j,k)%adj%i = i1_max/2**(l-1) - i
                            else
                              block(1,l)%bc%BCj(i,j,k)%adj%i = i
                            endif
                          case('j')
                            if (i1_sgn*j2_sgn<0) then
                              block(1,l)%bc%BCj(i,j,k)%adj%j = i1_max/2**(l-1) - i
                            else
                              block(1,l)%bc%BCj(i,j,k)%adj%j = i
                            endif
                          case('k')
                            if (i1_sgn*k2_sgn<0) then
                              block(1,l)%bc%BCj(i,j,k)%adj%k = i1_max/2**(l-1) - i
                            else
                              block(1,l)%bc%BCj(i,j,k)%adj%k = i
                            endif
                          endselect
                        case('k')
                          select case(orientation(o2)(o:o))
                          case('i')
                            if (k1_sgn*i2_sgn<0) then
                              block(1,l)%bc%BCj(i,j,k)%adj%i = k1_max/2**(l-1) - k
                            else
                              block(1,l)%bc%BCj(i,j,k)%adj%i = k
                            endif
                          case('j')
                            if (k1_sgn*j2_sgn<0) then
                              block(1,l)%bc%BCj(i,j,k)%adj%j = k1_max/2**(l-1) - k
                            else
                              block(1,l)%bc%BCj(i,j,k)%adj%j = k
                            endif
                          case('k')
                            if (k1_sgn*k2_sgn<0) then
                              block(1,l)%bc%BCj(i,j,k)%adj%k = k1_max/2**(l-1) - k
                            else
                              block(1,l)%bc%BCj(i,j,k)%adj%k = k
                            endif
                          endselect
                        endselect
                      enddo
                      select case(orientation(o2)(6:6))
                      case('i')
                        if (j1_max==0) then
                          block(1,l)%bc%BCj(i,j,k)%adj%i = (blocks(b2)%Ni/(2**(l-1))+j)*(1-j1_sgn*i2_sgn)/2 + &
                                                           (-j+1                      )*(1+j1_sgn*i2_sgn)/2
                        else
                          block(1,l)%bc%BCj(i,j,k)%adj%i =                                                  &
                            (j-blocks(b1)%Nj/(2**(l-1))+1                           )*(1-j1_sgn*i2_sgn)/2 + &
                            (  blocks(b2)%Ni/(2**(l-1))-(j-blocks(b1)%Nj/(2**(l-1))))*(1+j1_sgn*i2_sgn)/2
                        endif
                      case('j')
                        if (j1_max==0) then
                          block(1,l)%bc%BCj(i,j,k)%adj%j = (blocks(b2)%Nj/(2**(l-1))+j)*(1-j1_sgn*j2_sgn)/2 + &
                                                           (-j+1                      )*(1+j1_sgn*j2_sgn)/2
                        else
                          block(1,l)%bc%BCj(i,j,k)%adj%j =                                                  &
                            (j-blocks(b1)%Nj/(2**(l-1))+1                           )*(1-j1_sgn*j2_sgn)/2 + &
                            (  blocks(b2)%Nj/(2**(l-1))-(j-blocks(b1)%Nj/(2**(l-1))))*(1+j1_sgn*j2_sgn)/2
                        endif
                      case('k')
                        if (j1_max==0) then
                          block(1,l)%bc%BCj(i,j,k)%adj%k = (blocks(b2)%Nk/(2**(l-1))+j)*(1-j1_sgn*k2_sgn)/2 + &
                                                           (-j+1                      )*(1+j1_sgn*k2_sgn)/2
                        else
                          block(1,l)%bc%BCj(i,j,k)%adj%k =                                                  &
                            (j-blocks(b1)%Nj/(2**(l-1))+1                           )*(1-j1_sgn*k2_sgn)/2 + &
                            (  blocks(b2)%Nk/(2**(l-1))-(j-blocks(b1)%Nj/(2**(l-1))))*(1+j1_sgn*k2_sgn)/2
                        endif
                      endselect
                    enddo
                  enddo
                enddo
              case('k')
                ! the adjacent face is k-face for block 1
                if (k1_min==0) then
                  k1_min = 1-block(1,l)%mesh%gc(5)
                  k1_max = 0
                else
                  k1_min = blocks(b1)%Nk/(2**(l-1))
                  k1_max = blocks(b1)%Nk/(2**(l-1))+block(1,l)%mesh%gc(6)-1
                endif
                do k=k1_min,k1_max
                  do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
                    do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
                      block(1,l)%bc%BCk(i,j,k)%tp    = bc_adj
                      if (allocated(block(1,l)%bc%BCk(i,j,k)%adj)) deallocate(block(1,l)%bc%BCk(i,j,k)%adj)
                      allocate(block(1,l)%bc%BCk(i,j,k)%adj)
                      block(1,l)%bc%BCk(i,j,k)%adj%b = b2
                      do o=2,4,2
                        select case(orientation(o1)(o:o))
                        case('i')
                          select case(orientation(o2)(o:o))
                          case('i')
                            if (i1_sgn*i2_sgn<0) then
                              block(1,l)%bc%BCk(i,j,k)%adj%i = i1_max/2**(l-1) - i
                            else
                              block(1,l)%bc%BCk(i,j,k)%adj%i = i
                            endif
                          case('j')
                            if (i1_sgn*j2_sgn<0) then
                              block(1,l)%bc%BCk(i,j,k)%adj%j = i1_max/2**(l-1) - i
                            else
                              block(1,l)%bc%BCk(i,j,k)%adj%j = i
                            endif
                          case('k')
                            if (i1_sgn*k2_sgn<0) then
                              block(1,l)%bc%BCk(i,j,k)%adj%k = i1_max/2**(l-1) - i
                            else
                              block(1,l)%bc%BCk(i,j,k)%adj%k = i
                            endif
                          endselect
                        case('j')
                          select case(orientation(o2)(o:o))
                          case('i')
                            if (j1_sgn*i2_sgn<0) then
                              block(1,l)%bc%BCk(i,j,k)%adj%i = j1_max/2**(l-1) - j
                            else
                              block(1,l)%bc%BCk(i,j,k)%adj%i = j
                            endif
                          case('j')
                            if (j1_sgn*j2_sgn<0) then
                              block(1,l)%bc%BCk(i,j,k)%adj%j = j1_max/2**(l-1) - j
                            else
                              block(1,l)%bc%BCk(i,j,k)%adj%j = j
                            endif
                          case('k')
                            if (j1_sgn*k2_sgn<0) then
                              block(1,l)%bc%BCk(i,j,k)%adj%k = j1_max/2**(l-1) - j
                            else
                              block(1,l)%bc%BCk(i,j,k)%adj%k = j
                            endif
                          endselect
                        endselect
                      enddo
                      select case(orientation(o2)(6:6))
                      case('i')
                        if (k1_max==0) then
                          block(1,l)%bc%BCk(i,j,k)%adj%i = (blocks(b2)%Ni/(2**(l-1))+k)*(1-k1_sgn*i2_sgn)/2 + &
                                                           (-k+1                      )*(1+k1_sgn*i2_sgn)/2
                        else
                          block(1,l)%bc%BCk(i,j,k)%adj%i =                                                  &
                            (k-blocks(b1)%Nk/(2**(l-1))+1                           )*(1-k1_sgn*i2_sgn)/2 + &
                            (  blocks(b2)%Nk/(2**(l-1))-(k-blocks(b1)%Ni/(2**(l-1))))*(1+k1_sgn*i2_sgn)/2
                        endif
                      case('j')
                        if (k1_max==0) then
                          block(1,l)%bc%BCk(i,j,k)%adj%j = (blocks(b2)%Nj/(2**(l-1))+k)*(1-k1_sgn*j2_sgn)/2 + &
                                                           (-k+1                      )*(1+k1_sgn*j2_sgn)/2
                        else
                          block(1,l)%bc%BCk(i,j,k)%adj%j =                                                  &
                            (k-blocks(b1)%Nk/(2**(l-1))+1                           )*(1-k1_sgn*j2_sgn)/2 + &
                            (  blocks(b2)%Nj/(2**(l-1))-(k-blocks(b1)%Nk/(2**(l-1))))*(1+k1_sgn*j2_sgn)/2
                        endif
                      case('k')
                        if (k1_max==0) then
                          block(1,l)%bc%BCk(i,j,k)%adj%k = (blocks(b2)%Nk/(2**(l-1))+k)*(1-k1_sgn*k2_sgn)/2 + &
                                                           (-k+1                      )*(1+k1_sgn*k2_sgn)/2
                        else
                          block(1,l)%bc%BCk(i,j,k)%adj%k =                                                  &
                            (k-blocks(b1)%Nk/(2**(l-1))+1                           )*(1-k1_sgn*k2_sgn)/2 + &
                            (  blocks(b2)%Nk/(2**(l-1))-(k-blocks(b1)%Nk/(2**(l-1))))*(1+k1_sgn*k2_sgn)/2
                        endif
                      endselect
                    enddo
                  enddo
                enddo
              endselect
            enddo
          elseif (index(line1,tab//'b')>0) then
            cycle
          elseif (index(line1,tab//'e')>0) then
            cycle
          elseif (index(line1,tab//'v')>0) then
            cycle
          else
            exit
          endif
        enddo
      elseif (index(line,'# Boundary conditions and/or properties for domain.')>0)  then
        ! found a block boundary or initial conditions information
        do
          read(UnitScratch(0,b,1),'(A)',iostat=err) line1
          if (err /= 0) exit
          if (index(line1,tab//'f')>0) then
            do bc=1,Nbc
              if (index(line1,bc_list_str(bc))>0) then
                read(line1(index(line1,tab//'f',BACK=.true.)+2:),*)i1_min,j1_min,k1_min,i1_max,j1_max,k1_max
                i1_min = i1_min - 1 ; i1_max = i1_max - 1
                j1_min = j1_min - 1 ; j1_max = j1_max - 1
                k1_min = k1_min - 1 ; k1_max = k1_max - 1
                do l=1,global%mesh%Nl
                  if     (i1_min==i1_max) then
                    ! i face
                    if (i1_min==0) then
                      do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
                        do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
                          do i=1-block(1,l)%mesh%gc(1),0
                            block(1,l)%bc%BCi(i,j,k)%tp = bc_list(bc)
                            if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
                              if (allocated(block(1,l)%bc%BCi(i,j,k)%inf)) deallocate(block(1,l)%bc%BCi(i,j,k)%inf)
                              allocate(block(1,l)%bc%BCi(i,j,k)%inf)
                              block(1,l)%bc%BCi(i,j,k)%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
                                                                        index(line1,tab//'f')-1),1_I_P)
                            endif
                          enddo
                        enddo
                      enddo
                    else
                      do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
                        do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
                          do i=i1_max/2**(l-1)+1,i1_max/2**(l-1)+block(1,l)%mesh%gc(2)
                            block(1,l)%bc%BCi(i-1,j,k)%tp = bc_list(bc)
                            if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
                              if (allocated(block(1,l)%bc%BCi(i-1,j,k)%inf)) deallocate(block(1,l)%bc%BCi(i-1,j,k)%inf)
                              allocate(block(1,l)%bc%BCi(i-1,j,k)%inf)
                              block(1,l)%bc%BCi(i-1,j,k)%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
                                                                          index(line1,tab//'f')-1),1_I_P)
                            endif
                          enddo
                        enddo
                      enddo
                    endif
                  elseif (j1_min==j1_max) then
                    ! j face
                    if (j1_min==0) then
                      do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
                        do j=1-block(1,l)%mesh%gc(3),0
                          do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
                            block(1,l)%bc%BCj(i,j,k)%tp = bc_list(bc)
                            if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
                              if (allocated(block(1,l)%bc%BCj(i,j,k)%inf)) deallocate(block(1,l)%bc%BCj(i,j,k)%inf)
                              allocate(block(1,l)%bc%BCj(i,j,k)%inf)
                              block(1,l)%bc%BCj(i,j,k)%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
                                                                        index(line1,tab//'f')-1),1_I_P)
                            endif
                          enddo
                        enddo
                      enddo
                    else
                      do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
                        do j=j1_max/2**(l-1)+1,j1_max/2**(l-1)+block(1,l)%mesh%gc(4)
                          do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
                            block(1,l)%bc%BCj(i,j-1,k)%tp = bc_list(bc)
                            if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
                              if (allocated(block(1,l)%bc%BCj(i,j-1,k)%inf)) deallocate(block(1,l)%bc%BCj(i,j-1,k)%inf)
                              allocate(block(1,l)%bc%BCj(i,j-1,k)%inf)
                              block(1,l)%bc%BCj(i,j-1,k)%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
                                                                          index(line1,tab//'f')-1),1_I_P)
                            endif
                          enddo
                        enddo
                      enddo
                    endif
                  elseif (k1_min==k1_max) then
                    ! k face
                    if (k1_min==0) then
                      do k=1-block(1,l)%mesh%gc(5),0
                        do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
                          do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
                            block(1,l)%bc%BCk(i,j,k)%tp = bc_list(bc)
                            if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
                              if (allocated(block(1,l)%bc%BCk(i,j,k)%inf)) deallocate(block(1,l)%bc%BCk(i,j,k)%inf)
                              allocate(block(1,l)%bc%BCk(i,j,k)%inf)
                              block(1,l)%bc%BCk(i,j,k)%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
                                                                        index(line1,tab//'f')-1),1_I_P)
                            endif
                          enddo
                        enddo
                      enddo
                    else
                      do k=k1_max/2**(l-1)+1,k1_max/2**(l-1)+block(1,l)%mesh%gc(6)
                        do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
                          do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
                            block(1,l)%bc%BCk(i,j,k-1)%tp = bc_list(bc)
                            if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
                              if (allocated(block(1,l)%bc%BCk(i,j,k-1)%inf)) deallocate(block(1,l)%bc%BCk(i,j,k-1)%inf)
                              allocate(block(1,l)%bc%BCk(i,j,k-1)%inf)
                              block(1,l)%bc%BCk(i,j,k-1)%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
                                                                          index(line1,tab//'f')-1),1_I_P)
                            endif
                          enddo
                        enddo
                      enddo
                    endif
                  endif
                enddo
              endif
            enddo
          elseif (index(line1,tab//'b')>0) then
            if (index(line1,'BLK')>0) then
              inquire(file=adjustl(trim(line1(1:index(line1,tab)-1)))//'.itc',exist=is_file)
              if (.NOT.is_file) then
                write(stderr,'(A)')' File'
                write(stderr,'(A)')' '//adjustl(trim(line1(1:index(line1,tab)-1)))//'.itc'
                write(stderr,'(A)')' File Not Found'
                stop
              endif
              Unit_itc = Get_Unit()
              open(unit=Unit_itc,file=adjustl(trim(line1(1:index(line1,tab)-1)))//'.itc')
              do v=1,global%fluid%Ns
                read(Unit_itc,*) blocks(b)%P%r(v)
              enddo
              read(Unit_itc,*) blocks(b)%P%v%x
              read(Unit_itc,*) blocks(b)%P%v%y
              read(Unit_itc,*) blocks(b)%P%v%z
              read(Unit_itc,*) blocks(b)%P%p
              read(Unit_itc,*) blocks(b)%P%d
              read(Unit_itc,*) blocks(b)%P%g
              close(Unit_itc)
              do l=1,global%mesh%Nl
                block(1,l)%fluid%P = blocks(b)%P
              enddo
            endif
          elseif (index(line1,tab//'e')>0) then
            cycle
          elseif (index(line1,tab//'v')>0) then
            cycle
          else
            exit
          endif
        enddo
      endif
    enddo
    close(UnitScratch(0,b,1))
    ! storing the boundary and initial conditions in the scratch files
    do l=1,global%mesh%Nl
      ! boundary conditions data
      err = write(UnitScratch(2,b,l),block(1,l)%bc%BCi)
      err = write(UnitScratch(2,b,l),block(1,l)%bc%BCj)
      err = write(UnitScratch(2,b,l),block(1,l)%bc%BCk)
      ! initial conditions data
      write(UnitScratch(3,b,l),iostat=err)block(1,l)%fluid%Dt
      err = write(UnitScratch(3,b,l),block(1,l)%fluid%P)
      ! rewinding scratch files
      rewind(UnitScratch(2,b,l))
      rewind(UnitScratch(3,b,l))
    enddo
  enddo
  global%mesh%Nb = global%mesh%Nb_tot
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_icemcfd
endprogram IBM
