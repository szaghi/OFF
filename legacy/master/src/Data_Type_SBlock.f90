!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_SBlockDerivedType Data_Type_SBlock
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_SBlockPrivateProcedure Data_Type_SBlock
!> @}

!> @brief Module Data_Type_SBlock contains the definition of Type_SBlock (structured block) type and useful procedures for its
!> handling.
module Data_Type_SBlock
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision           !< Integers and reals precision definition.
USE Data_Type_BC           !< Definition of Type_BC.
USE Data_Type_Cell         !< Definition of Type_Cell.
!USE Data_Type_Conservative !< Definition of Type_Conservative.
USE Data_Type_Face         !< Definition of Type_Face.
USE Data_Type_Global       !< Definition of Type_Global.
USE Data_Type_Primitive    !< Definition of Type_Primitive.
USE Data_Type_Vector       !< Definition of Type_Vector.
USE Lib_IO_Misc            !< Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------

!> @brief Derived type containing structured block-level data.
!> Structured block-level type contains data (mesh, boundary conditions and fluid dynamic data) of structured (implicit
!> connectivity) numerical grid.
!> @ingroup Data_Type_SBlockDerivedType
type, public:: Type_SBlock
  type(Type_Global), pointer::     global                 !< Pointer to global-level data.
  integer(I1P)::                   gc(1:6)=&
                                          [1_I1P,1_I1P, & ! gc(1) => left i, gc(2) => right i.
                                           1_I1P,1_I1P, & ! gc(3) => left j, gc(4) => right j.
                                           1_I1P,1_I1P  & ! gc(5) => left k, gc(6) => right k.
                                           ]              !< Number of ghost cells for the 6 faces of the block.
  integer(I4P)::                   Ni     = 0_I4P         !< Number of cells in i direction.
  integer(I4P)::                   Nj     = 0_I4P         !< Number of cells in j direction.
  integer(I4P)::                   Nk     = 0_I4P         !< Number of cells in k direction.
  type(Type_Vector), allocatable:: node(:,:,:)            !< Nodes coord.  [0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)].
  type(Type_Face),   allocatable:: Fi(:,:,:)              !< Faces i data  [0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Face),   allocatable:: Fj(:,:,:)              !< Faces j data  [1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Face),   allocatable:: Fk(:,:,:)              !< Faces k data  [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)].
  type(Type_Cell),   allocatable:: C(:,:,:)               !< Cells data    [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  contains
    procedure:: alloc => alloc_block                  ! Procedure for allocating memory.
    procedure:: free => free_block                    ! Procedure for freeing memory.
    procedure:: save_mesh => save_bmesh               ! Procedure for saving the mesh data.
    procedure:: load_mesh_dims => load_bmesh_dims     ! Procedure for loading the mesh data dimensions.
    procedure:: load_mesh => load_bmesh               ! Procedure for loading the mesh data.
    procedure:: print_info_mesh => print_info_bmesh   ! Procedure for printing mesh data.
    procedure:: save_bc => save_bbc                   ! Procedure for saving the bc data.
    procedure:: load_bc => load_bbc                   ! Procedure for loading the bc data.
    procedure:: save_fluid => save_bfluid             ! Procedure for saving the fluid dynamic data.
    procedure:: load_fluid => load_bfluid             ! Procedure for loading the fluid dynamic data.
    procedure:: print_info_fluid => print_info_bfluid ! Procedure for printing fluid dynamic data.
    procedure:: metrics                               ! Procedure for computing block metrics.
    procedure:: metrics_correction                    ! Procedure for correcting block metrics of bc cells.
    procedure:: node2center                           ! Procedure for computing cell center coo from cell nodes.
endtype Type_SBlock
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_SBlockPrivateProcedure
  !> @{
  !> Subroutine for allocating dynamic data of Type_SBlock variables.
  subroutine alloc_block(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block data.
  integer(I1P)::                      gc(1:6)  !< Temporary variable  for storing block ghost cells number.
  integer(I_P)::                      Ni,Nj,Nk !< Temporary variables for storing block dimensions.
  integer(I_P)::                      Ns       !< Temporary variable  for storing number of species.
  integer(I1P)::                      rk_ord   !< Temporary variable  for storing rk_ord.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! free dynamic data of block if previously allocated
  call block%free
  ! storing block dimensions into temporary variables to simplify the code
  gc(1:6) = block%gc(1:6)
  Ni      = block%Ni
  Nj      = block%Nj
  Nk      = block%Nk
  Ns      = block%global%Ns
  rk_ord  = block%global%rk_ord
  allocate(block%node(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  allocate(block%Fi  (0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; call block%Fi%init
  allocate(block%Fj  (1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; call block%Fj%init
  allocate(block%Fk  (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; call block%Fk%init
  allocate(block%C   (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; call block%C%init(Ns=Ns,rk_ord=rk_ord)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_block

  !> Subroutine for freeing dynamic data of Type_SBlock variables.
  subroutine free_block(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block !< Block data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(block%node)) deallocate(block%node)
  if (allocated(block%Fi)) deallocate(block%Fi)
  if (allocated(block%Fj)) deallocate(block%Fj)
  if (allocated(block%Fk)) deallocate(block%Fk)
  if (allocated(block%C)) deallocate(block%C)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_block

  !> Function for saving block mesh file.
  !> @return \b err integer(I_P) variable.
  function save_bmesh(block,ascii,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(IN):: block    !< Block level data.
  logical, optional,  intent(IN):: ascii    !< Flag for ascii file.
  character(*),       intent(IN):: filename !< Name of file where mesh variables are saved.
  integer(I_P)::                   err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                   UnitFree !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(ascii)) then
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), action = 'WRITE', form = 'FORMATTED')
    write(UnitFree,'(6('//FI1P//',1X))',iostat=err)block%gc
    write(UnitFree,'(3('//FI_P//',1X))',iostat=err)block%Ni,block%Nj,block%Nk
    err = write_vector(array3D=block%node,format='*',unit=UnitFree)
    !err = write_cell(array3D=block%C%cell,format='*',unit=UnitFree)
    close(UnitFree)
  else
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), action = 'WRITE', form = 'UNFORMATTED')
    write(UnitFree,iostat=err)block%gc
    write(UnitFree,iostat=err)block%Ni,block%Nj,block%Nk
    err = write_vector(array3D=block%node,unit=UnitFree)
    !err = write_cell(array3D=block%C%cell,unit=UnitFree)
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bmesh

  !> Function for loading the mesh data dimensions of block from the mesh file "filename".
  !> @return \b err integer(I_P) variable.
  function load_bmesh_dims(block,ascii,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block level data.
  logical, optional,  intent(IN)::    ascii    !< Flag for ascii file.
  character(*),       intent(IN)::    filename !< Name of file where mesh variables are saved.
  integer(I_P)::                      UnitFree !< Free logic unit.
  logical::                           is_file  !< Flag for inquiring the presence of mesh file.
  integer(I_P)::                      err      !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=trim(filename),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(block%global%myrank,filename,'load_bmesh_dims')
  if (present(ascii)) then
    open(unit = Get_Unit(UnitFree), file = trim(filename), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err)block%gc
    read(UnitFree,*,iostat=err)block%Ni,block%Nj,block%Nk
    close(UnitFree)
  else
    open(unit = Get_Unit(UnitFree ), file = trim(filename), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)block%gc
    read(UnitFree,iostat=err)block%Ni,block%Nj,block%Nk
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bmesh_dims

  !> Function for loading block mesh file.
  !> @return \b err integer(I_P) variable.
  function load_bmesh(block,ascii,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block level data.
  logical, optional,  intent(IN)::    ascii    !< Flag for ascii file.
  character(*),       intent(IN)::    filename !< Name of file where mesh variables are saved.
  integer(I_P)::                      err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                      UnitFree !< Free logic unit.
  logical::                           is_file  !< Flag for inquiring the presence of mesh file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(block%global%myrank,filename,'load_bmesh')
  if (present(ascii)) then
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err)block%gc
    read(UnitFree,*,iostat=err)block%Ni,block%Nj,block%Nk
    err = read_vector(array3D=block%node,format='*',unit=UnitFree)
    !err = read_cell(array3D=block%C%cell,format='*',unit=UnitFree)
    close(UnitFree)
  else
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)block%gc
    read(UnitFree,iostat=err)block%Ni,block%Nj,block%Nk
    err = read_vector(array3D=block%node,unit=UnitFree)
    !err = read_cell(array3D=block%C%cell,unit=UnitFree)
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bmesh

  !> Function for printing to standard output info of block mesh data.
  !> @return \b err integer(I_P) variable.
  function print_info_bmesh(block,blk,grl) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(IN):: block       !< Block level data.
  integer(I_P),       intent(IN):: blk         !< Actual block number.
  integer(I_P),       intent(IN):: grl         !< Actual grid level number.
  integer(I_P)::                   err         !< Error trapping flag: 0 no errors, >0 error occurs.
  character(DR_P)::                vmax,vmin   !< String for printing max and min of variables.
  character(DR_P)::                xmm,ymm,zmm !< String for printing max and min of variables.
  type(Type_Vector), allocatable:: NFi(:,:,:)  !< |
  type(Type_Vector), allocatable:: NFj(:,:,:)  !< |
  type(Type_Vector), allocatable:: NFk(:,:,:)  !< |
  real(R_P),         allocatable:: Si (:,:,:)  !< | Dummy variables for printing only internal cells info.
  real(R_P),         allocatable:: Sj (:,:,:)  !< |
  real(R_P),         allocatable:: Sk (:,:,:)  !< |
  real(R_P),         allocatable:: V  (:,:,:)  !< |
  integer(I_P)::                   Ni,Nj,Nk    !< Temporary variables for storing block dimensions.
  character(DI_P)::                rks         !< String containing myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,block%global%myrank))
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)trim(rks)//' Mesh infos of myrank '//trim(rks)
  write(stdout,'(A)',iostat=err)trim(rks)//'  Global number of blocks '//trim(str(.true.,block%global%Nb_tot))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Block '//trim(str(.true.,blk))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Grid level '//trim(str(.true.,grl))
  Ni = block%Ni
  Nj = block%Nj
  Nk = block%Nk
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: left  i '//trim(str(.true.,block%gc(1)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: right i '//trim(str(.true.,block%gc(2)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: left  j '//trim(str(.true.,block%gc(3)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: right j '//trim(str(.true.,block%gc(4)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: left  k '//trim(str(.true.,block%gc(5)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: right k '//trim(str(.true.,block%gc(6)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Number of cells Ni '//trim(str(.true.,Ni))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Number of cells Nj '//trim(str(.true.,Nj))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Number of cells Nk '//trim(str(.true.,Nk))
  if (allocated(NFi)) deallocate(NFi) ; allocate(NFi(0:Ni,1:Nj,1:Nk)) ; NFi(0:Ni,1:Nj,1:Nk) = block%Fi(0:Ni,1:Nj,1:Nk)%N
  if (allocated(NFj)) deallocate(NFj) ; allocate(NFj(1:Ni,0:Nj,1:Nk)) ; NFj(1:Ni,0:Nj,1:Nk) = block%Fj(1:Ni,0:Nj,1:Nk)%N
  if (allocated(NFk)) deallocate(NFk) ; allocate(NFk(1:Ni,1:Nj,0:Nk)) ; NFk(1:Ni,1:Nj,0:Nk) = block%Fk(1:Ni,1:Nj,0:Nk)%N
  if (allocated(Si )) deallocate(Si ) ; allocate(Si (0:Ni,1:Nj,1:Nk)) ; Si (0:Ni,1:Nj,1:Nk) = block%Fi(0:Ni,1:Nj,1:Nk)%S
  if (allocated(Sj )) deallocate(Sj ) ; allocate(Sj (1:Ni,0:Nj,1:Nk)) ; Sj (1:Ni,0:Nj,1:Nk) = block%Fj(1:Ni,0:Nj,1:Nk)%S
  if (allocated(Sk )) deallocate(Sk ) ; allocate(Sk (1:Ni,1:Nj,0:Nk)) ; Sk (1:Ni,1:Nj,0:Nk) = block%Fk(1:Ni,1:Nj,0:Nk)%S
  if (allocated(V  )) deallocate(V  ) ; allocate(V  (1:Ni,1:Nj,1:Nk)) ; V  (1:Ni,1:Nj,1:Nk) = block%C (1:Ni,1:Nj,1:Nk)%V
  xmm=trim(str(n=maxval(NFi(:,:,:)%x))) ; ymm=trim(str(n=maxval(NFi(:,:,:)%y))) ; zmm=trim(str(n=maxval(NFi(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFi-x-max '//trim(xmm)//' NFi-y-max '//trim(ymm)//' NFi-z-max '//trim(zmm)
  xmm=trim(str(n=minval(NFi(:,:,:)%x))) ; ymm=trim(str(n=minval(NFi(:,:,:)%y))) ; zmm=trim(str(n=minval(NFi(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFi-x-min '//trim(xmm)//' NFi-y-min '//trim(ymm)//' NFi-z-min '//trim(zmm)
  xmm=trim(str(n=maxval(NFj(:,:,:)%x))) ; ymm=trim(str(n=maxval(NFj(:,:,:)%y))) ; zmm=trim(str(n=maxval(NFj(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFj-x-max '//trim(xmm)//' NFj-y-max '//trim(ymm)//' NFj-z-max '//trim(zmm)
  xmm=trim(str(n=minval(NFj(:,:,:)%x))) ; ymm=trim(str(n=minval(NFj(:,:,:)%y))) ; zmm=trim(str(n=minval(NFj(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFj-x-min '//trim(xmm)//' NFj-y-min '//trim(ymm)//' NFj-z-min '//trim(zmm)
  xmm=trim(str(n=maxval(NFk(:,:,:)%x))) ; ymm=trim(str(n=maxval(NFk(:,:,:)%y))) ; zmm=trim(str(n=maxval(NFk(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFk-x-max '//trim(xmm)//' NFk-y-max '//trim(ymm)//' NFk-z-max '//trim(zmm)
  xmm=trim(str(n=minval(NFk(:,:,:)%x))) ; ymm=trim(str(n=minval(NFk(:,:,:)%y))) ; zmm=trim(str(n=minval(NFk(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFk-x-min '//trim(xmm)//' NFk-y-min '//trim(ymm)//' NFk-z-min '//trim(zmm)
  vmax=trim(str(n=maxval(Si(:,:,:)))) ; vmin=trim(str(n=minval(Si(:,:,:))))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Si-max '//trim(vmax)//' Si-min '//trim(vmin)
  vmax=trim(str(n=maxval(Sj(:,:,:)))) ; vmin=trim(str(n=minval(Sj(:,:,:))))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Sj-max '//trim(vmax)//' Sj-min '//trim(vmin)
  vmax=trim(str(n=maxval(Sk(:,:,:)))) ; vmin=trim(str(n=minval(Sk(:,:,:))))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Sk-max '//trim(vmax)//' Sk-min '//trim(vmin)
  vmax=trim(str(n=maxval(V(:,:,:)))) ; vmin=trim(str(n=minval(V(:,:,:))))
  write(stdout,'(A)',iostat=err)trim(rks)//'  V-max '//trim(vmax)//' V-min '//trim(vmin)
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,*)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction print_info_bmesh

  !> Function for saving block boundary conditions file.
  !> @return \b err integer(I_P) variable.
  function save_bbc(block,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(IN):: block    !< Block level data.
  character(*),       intent(IN):: filename !< Name of file where mesh variables are saved.
  integer(I_P)::                   err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                   UnitFree !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  open(unit=Get_Unit(UnitFree), file=adjustl(trim(filename)), action='WRITE', form='UNFORMATTED')
  err = write_bc(array3D=block%Fi%BC,unit=UnitFree)
  err = write_bc(array3D=block%Fj%BC,unit=UnitFree)
  err = write_bc(array3D=block%Fk%BC,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bbc

  !> Function for loading block boundary conditions file.
  !> @return \b err integer(I_P) variable.
  function load_bbc(block,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block level data.
  character(*),       intent(IN)::    filename !< Name of file where mesh variables are saved.
  integer(I_P)::                      err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                      UnitFree !< Free logic unit.
  logical::                           is_file  !< Flag for inquiring the presence of boundary conditions file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(block%global%myrank,filename,'load_bbc')
  open(unit=Get_Unit(UnitFree), file=adjustl(trim(filename)), status='OLD', action='READ', form='UNFORMATTED')
  err = read_bc(array3D=block%Fi%BC,unit=UnitFree)
  err = read_bc(array3D=block%Fj%BC,unit=UnitFree)
  err = read_bc(array3D=block%Fk%BC,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bbc

  !> Function for saving fluid dynamic file.
  !> @return \b err integer(I_P) variable.
  function save_bfluid(block,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(IN):: block    !< Block level data.
  character(*),       intent(IN):: filename !< Name of file where mesh variables are saved.
  integer(I_P)::                   err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                   UnitFree !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), action = 'WRITE', form = 'UNFORMATTED')
  write(UnitFree,iostat=err)block%global%Ns
  write(UnitFree,iostat=err)block%global%cp0,block%global%cv0
  write(UnitFree,iostat=err)block%global%n
  write(UnitFree,iostat=err)block%global%t
  write(UnitFree,iostat=err)block%C%Dt
  err = write_primitive(array3D=block%C%P,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bfluid

  !> Function for loading fluid dynamic file.
  !> @return \b err integer(I_P) variable.
  function load_bfluid(block,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block level data.
  character(*),       intent(IN)::    filename !< Name of file where mesh variables are saved.
  integer(I_P)::                      err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                      UnitFree !< Free logic unit.
  logical::                           is_file  !< Flag for inquiring the presence of fluid dynamic file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(block%global%myrank,filename,'load_bfluid')
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
  read(UnitFree,iostat=err)block%global%Ns
  read(UnitFree,iostat=err)block%global%cp0,block%global%cv0
  read(UnitFree,iostat=err)block%global%n
  read(UnitFree,iostat=err)block%global%t
  read(UnitFree,iostat=err)block%C%Dt
  err = read_primitive(array3D=block%C%P,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bfluid

  !> Function for printing to standard output info of block fluid data.
  !> @return \b err integer(I_P) variable.
  function print_info_bfluid(block,blk,grl) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(IN):: block          !< Block level data.
  integer(I_P),       intent(IN):: blk            !< Actual block number.
  integer(I_P),       intent(IN):: grl            !< Actual grid level number.
  integer(I_P)::                   err            !< Error trapping flag: 0 no errors, >0 error occurs.
  character(DR_P)::                vmax,vmin      !< String for printing max and min of variables.
  real(R_P), allocatable::         dummy(:,:,:,:) !<  Dummy variables for printing only internal cells info.
  integer(I_P)::                   Ni,Nj,Nk,Np,Nc !< Temporary variables for storing block dimensions.
  character(DI_P)::                rks            !< String containing myrank.
  integer(I_P)::                   i,j,k,s        !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,block%global%myrank))
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)trim(rks)//' Fluid infos of myrank '//trim(rks)
  write(stdout,'(A)',iostat=err)trim(rks)//'  Global number of blocks '//trim(str(.true.,block%global%Nb_tot))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Number of initial species '//trim(str(.true.,block%global%Ns))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Block '//trim(str(.true.,blk))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Grid level '//trim(str(.true.,grl))
  Ni = block%Ni
  Nj = block%Nj
  Nk = block%Nk
  Np = block%global%Np
  Nc = block%global%Nc
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Nc,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(:,i,j,k) = block%C(i,j,k)%U%cons2array()
      enddo
    enddo
  enddo
  write(stdout,'(A)',iostat=err)trim(rks)//'  Conservative variables'
  do s=1,block%global%Nc
    vmax=trim(str(n=maxval(dummy(s,:,:,:)))) ; vmin=trim(str(n=minval(dummy(s,:,:,:))))
    write(stdout,'(A)',iostat=err)trim(rks)//'  U('//trim(str(.true.,s))//'): max '//trim(vmax)//' min '//trim(vmin)
  enddo
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Np,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(:,i,j,k) = block%C(i,j,k)%P%prim2array()
      enddo
    enddo
  enddo
  write(stdout,'(A)',iostat=err)trim(rks)//'  Primitive variables'
  do s=1,block%global%Np
    vmax=trim(str(n=maxval(dummy(s,:,:,:)))) ; vmin=trim(str(n=minval(dummy(s,:,:,:))))
    write(stdout,'(A)',iostat=err)trim(rks)//'  P('//trim(str(.true.,s))//'): max '//trim(vmax)//' min '//trim(vmin)
  enddo
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(1,i,j,k) = block%C(i,j,k)%Dt
      enddo
    enddo
  enddo
  vmax=trim(str(n=maxval(dummy(1,:,:,:)))) ; vmin=trim(str(n=minval(dummy(1,:,:,:))))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Time steps Dt: max '//trim(vmax)//' min '//trim(vmin)
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,*)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction print_info_bfluid

  !> Subroutine for computing the metrics of structured block.
  subroutine metrics(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block             !< Block-level data.
  type(Type_Vector)::                 NFS,s1,s2,nd,db   !< Dummy vector variables.
  real(R_P)::                         signi,signj,signk !< Dummy variables for checking the directions of normals.
  real(R_P)::                         Vx,Vy,Vz          !< Dummy variables for computing volume.
  real(R_P)::                         xp,yp,zp          !< Dummy variables for computing face coordinates.
  real(R_P)::                         xm,ym,zm          !< Dummy variables for computing face coordinates.
  integer(I_P)::                      i,j,k             !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing faces normals
  ! positioning at the middle of the block
  i = max(1,block%Ni)
  j = max(1,block%Nj)
  k = max(1,block%Nk)
  ! checking the direction of i normals
  s1 = block%node(i,j  ,k) - block%node(i,j-1,k-1)
  s2 = block%node(i,j-1,k) - block%node(i,j,  k-1)
  nd = s1.cross.s2
  s1 = 0.25_R_P*(block%node(i,  j,k)+block%node(i,  j-1,k)+block%node(i,  j,k-1)+block%node(i,  j-1,k-1))
  s2 = 0.25_R_P*(block%node(i-1,j,k)+block%node(i-1,j-1,k)+block%node(i-1,j,k-1)+block%node(i-1,j-1,k-1))
  db = s1 - s2
  signi = sign(1._R_P,(nd.dot.db))
  ! checking the direction of j normals
  s1 = block%node(i,j,k  ) - block%node(i-1,j,k-1)
  s2 = block%node(i,j,k-1) - block%node(i-1,j,k  )
  nd = s1.cross.s2
  s1 = 0.25_R_P*(block%node(i,j,  k)+block%node(i-1,j,  k)+block%node(i,j,  k-1)+block%node(i-1,j,  k-1))
  s2 = 0.25_R_P*(block%node(i,j-1,k)+block%node(i-1,j-1,k)+block%node(i,j-1,k-1)+block%node(i-1,j-1,k-1))
  db = s1 - s2
  signj = sign(1._R_P,(nd.dot.db))
  ! checking the direction of k normals
  s1 = block%node(i,  j,k) - block%node(i-1,j-1,k)
  s2 = block%node(i-1,j,k) - block%node(i,  j-1,k)
  nd = s1.cross.s2
  s1 = 0.25_R_P*(block%node(i,j,k  )+block%node(i-1,j,k  )+block%node(i,j-1,k  )+block%node(i-1,j-1,k  ))
  s2 = 0.25_R_P*(block%node(i,j,k-1)+block%node(i-1,j,k-1)+block%node(i,j-1,k-1)+block%node(i-1,j-1,k-1))
  db = s1 - s2
  signk = sign(1._R_P,(nd.dot.db))
  !!$OMP PARALLEL DEFAULT(NONE)                        &
  !!$OMP PRIVATE(i,j,k,NFS,Vx,Vy,Vz,xp,yp,zp,xm,ym,zm) &
  !!$OMP SHARED(block,signi,signj,signk)
  !!$OMP DO
  do k=1,block%Nk
    do j=1,block%Nj
      do i=0,block%Ni
        NFS = face_normal4(pt1 = block%node(i,j-1,k-1), &
                           pt2 = block%node(i,j  ,k-1), &
                           pt3 = block%node(i,j  ,k  ), &
                           pt4 = block%node(i,j-1,k  ))
        NFS = NFS*signi
        block%Fi(i,j,k)%N = normalize(NFS)
        block%Fi(i,j,k)%S =    normL2(NFS)
      enddo
    enddo
  enddo
  !!$OMP DO
  do k=1,block%Nk
    do j=0,block%Nj
      do i=1,block%Ni
        NFS = face_normal4(pt1 = block%node(i-1,j,k-1), &
                           pt2 = block%node(i-1,j,k  ), &
                           pt3 = block%node(i  ,j,k  ), &
                           pt4 = block%node(i  ,j,k-1))
        NFS = NFS*signj
        block%Fj(i,j,k)%N = normalize(NFS)
        block%Fj(i,j,k)%S =    normL2(NFS)
      enddo
    enddo
  enddo
  !!$OMP DO
  do k=0,block%Nk
    do j=1,block%Nj
      do i=1,block%Ni
        NFS = face_normal4(pt1 = block%node(i-1,j-1,k), &
                           pt2 = block%node(i  ,j-1,k), &
                           pt3 = block%node(i  ,j  ,k), &
                           pt4 = block%node(i-1,j  ,k))
        NFS = NFS*signk
        block%Fk(i,j,k)%N = normalize(NFS)
        block%Fk(i,j,k)%S =    normL2(NFS)
      enddo
    enddo
  enddo
  ! computing finte volumes
  !!$OMP DO
  do k=1,block%Nk
    do j=1,block%Nj
      do i=1,block%Ni
        Vx = 0._R_P
        Vy = 0._R_P
        Vz = 0._R_P

        xp = 0.25_R_P*(block%node(i  ,j  ,k  )%x + block%node(i  ,j  ,k-1)%x + &
                       block%node(i  ,j-1,k  )%x + block%node(i  ,j-1,k-1)%x)
        yp = 0.25_R_P*(block%node(i  ,j  ,k  )%y + block%node(i  ,j  ,k-1)%y + &
                       block%node(i  ,j-1,k  )%y + block%node(i  ,j-1,k-1)%y)
        zp = 0.25_R_P*(block%node(i  ,j  ,k  )%z + block%node(i  ,j  ,k-1)%z + &
                       block%node(i  ,j-1,k  )%z + block%node(i  ,j-1,k-1)%z)
        xm = 0.25_R_P*(block%node(i-1,j  ,k  )%x + block%node(i-1,j  ,k-1)%x + &
                       block%node(i-1,j-1,k  )%x + block%node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(block%node(i-1,j  ,k  )%y + block%node(i-1,j  ,k-1)%y + &
                       block%node(i-1,j-1,k  )%y + block%node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(block%node(i-1,j  ,k  )%z + block%node(i-1,j  ,k-1)%z + &
                       block%node(i-1,j-1,k  )%z + block%node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*block%Fi(i,j,k)%N%x*block%Fi(i,j,k)%S - xm*block%Fi(i-1,j,k)%N%x*block%Fi(i-1,j,k)%S
        Vy = Vy + yp*block%Fi(i,j,k)%N%y*block%Fi(i,j,k)%S - ym*block%Fi(i-1,j,k)%N%y*block%Fi(i-1,j,k)%S
        Vz = Vz + zp*block%Fi(i,j,k)%N%z*block%Fi(i,j,k)%S - zm*block%Fi(i-1,j,k)%N%z*block%Fi(i-1,j,k)%S

        xp = 0.25_R_P*(block%node(i  ,j  ,k  )%x + block%node(i  ,j  ,k-1)%x + &
                       block%node(i-1,j  ,k  )%x + block%node(i-1,j  ,k-1)%x)
        yp = 0.25_R_P*(block%node(i  ,j  ,k  )%y + block%node(i  ,j  ,k-1)%y + &
                       block%node(i-1,j  ,k  )%y + block%node(i-1,j  ,k-1)%y)
        zp = 0.25_R_P*(block%node(i  ,j  ,k  )%z + block%node(i  ,j  ,k-1)%z + &
                       block%node(i-1,j  ,k  )%z + block%node(i-1,j  ,k-1)%z)
        xm = 0.25_R_P*(block%node(i  ,j-1,k  )%x + block%node(i  ,j-1,k-1)%x + &
                       block%node(i-1,j-1,k  )%x + block%node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(block%node(i  ,j-1,k  )%y + block%node(i  ,j-1,k-1)%y + &
                       block%node(i-1,j-1,k  )%y + block%node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(block%node(i  ,j-1,k  )%z + block%node(i  ,j-1,k-1)%z + &
                       block%node(i-1,j-1,k  )%z + block%node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*block%Fj(i,j,k)%N%x*block%Fj(i,j,k)%S - xm*block%Fj(i,j-1,k)%N%x*block%Fj(i,j-1,k)%S
        Vy = Vy + yp*block%Fj(i,j,k)%N%y*block%Fj(i,j,k)%S - ym*block%Fj(i,j-1,k)%N%y*block%Fj(i,j-1,k)%S
        Vz = Vz + zp*block%Fj(i,j,k)%N%z*block%Fj(i,j,k)%S - zm*block%Fj(i,j-1,k)%N%z*block%Fj(i,j-1,k)%S

        xp = 0.25_R_P*(block%node(i  ,j  ,k  )%x + block%node(i  ,j-1,k  )%x + &
                       block%node(i-1,j  ,k  )%x + block%node(i-1,j-1,k  )%x)
        yp = 0.25_R_P*(block%node(i  ,j  ,k  )%y + block%node(i  ,j-1,k  )%y + &
                       block%node(i-1,j  ,k  )%y + block%node(i-1,j-1,k  )%y)
        zp = 0.25_R_P*(block%node(i  ,j  ,k  )%z + block%node(i  ,j-1,k  )%z + &
                       block%node(i-1,j  ,k  )%z + block%node(i-1,j-1,k  )%z)
        xm = 0.25_R_P*(block%node(i  ,j  ,k-1)%x + block%node(i  ,j-1,k-1)%x + &
                       block%node(i-1,j  ,k-1)%x + block%node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(block%node(i  ,j  ,k-1)%y + block%node(i  ,j-1,k-1)%y + &
                       block%node(i-1,j  ,k-1)%y + block%node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(block%node(i  ,j  ,k-1)%z + block%node(i  ,j-1,k-1)%z + &
                       block%node(i-1,j  ,k-1)%z + block%node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*block%Fk(i,j,k)%N%x*block%Fk(i,j,k)%S - xm*block%Fk(i,j,k-1)%N%x*block%Fk(i,j,k-1)%S
        Vy = Vy + yp*block%Fk(i,j,k)%N%y*block%Fk(i,j,k)%S - ym*block%Fk(i,j,k-1)%N%y*block%Fk(i,j,k-1)%S
        Vz = Vz + zp*block%Fk(i,j,k)%N%z*block%Fk(i,j,k)%S - zm*block%Fk(i,j,k-1)%N%z*block%Fk(i,j,k-1)%S

        block%C(i,j,k)%V = max(Vx,Vy,Vz)
      enddo
    enddo
  enddo
  !!$OMP END PARALLEL
#ifdef NULi
  block%Fj(:,:,:)%N%x=0._R_P
  block%Fk(:,:,:)%N%x=0._R_P
  block%Fi(:,:,:)%N%y=0._R_P
  block%Fi(:,:,:)%N%z=0._R_P
#endif
#ifdef NULj
  block%Fi(:,:,:)%N%y=0._R_P
  block%Fk(:,:,:)%N%y=0._R_P
  block%Fj(:,:,:)%N%x=0._R_P
  block%Fj(:,:,:)%N%z=0._R_P
#endif
#ifdef NULk
  block%Fi(:,:,:)%N%z=0._R_P
  block%Fj(:,:,:)%N%z=0._R_P
  block%Fk(:,:,:)%N%x=0._R_P
  block%Fk(:,:,:)%N%y=0._R_P
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine metrics

  !> Subroutine for correcting the metrics of natural (and negative volume) boundary conditions cells of structured block.
  subroutine metrics_correction(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block-level data.
  logical::                           correct  !< Flag for inquiring if metrics must be corrected.
  logical::                           wall     !< Flag for inquiring if bc is "wall-type": different corrections must be used.
  real(R_P)::                         tm       !< Tangential metrics parameter (-1 for wall-type bc).
  real(R_P)::                         sn       !< Normal metrics coefficient correction.
  integer(I_P)::                      Ni,Nj,Nk !< Dimensions of the block.
  integer(I_P)::                      i,j,k    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = block%Ni
  Nj = block%Nj
  Nk = block%Nk
  !!$OMP PARALLEL DEFAULT(NONE)            &
  !!$OMP PRIVATE(i,j,k,correct,wall,tm,sn) &
  !!$OMP SHARED(Ni,Nj,Nk,block)
  ! left i
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      correct = ((block%Fi(0,j,k)%BC%tp==bc_ext).OR. &
                 (block%Fi(0,j,k)%BC%tp==bc_ref).OR. &
                 (block%Fi(0,j,k)%BC%tp==bc_in1).OR. &
                 (block%Fi(0,j,k)%BC%tp==bc_in2).OR. &
                 (block%C(  0,j,k)%V<(0.2_R_P*block%C(1,j,k)%V)))
      wall    = ((block%Fi(0,j,k)%BC%tp==bc_ext).OR.(block%Fi(0,j,k)%BC%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%Fi(1,j,k)%N*block%Fi(1,j,k)%S).dot.block%Fi(0,j,k)%N)
         block%Fi(-1,j,  k  )%N = -(block%Fi(1,j,k)%N*block%Fi(1,j,k)%S) + sn*block%Fi(0,j,k)%N
         block%Fi(-1,j,  k  )%S = normL2(   block%Fi(-1,j,k)%N)
         block%Fi(-1,j,  k  )%N = normalize(block%Fi(-1,j,k)%N)
         ! tangential metrics
         block%Fj( 0,j  ,k  )%N = tm*block%Fj(1,j  ,k  )%N
         block%Fj( 0,j-1,k  )%N = tm*block%Fj(1,j-1,k  )%N
         block%Fj( 0,j  ,k-1)%N = tm*block%Fj(1,j  ,k-1)%N
         block%Fj( 0,j-1,k-1)%N = tm*block%Fj(1,j-1,k-1)%N
         block%Fj( 0,j  ,k  )%S = tm*block%Fj(1,j  ,k  )%S
         block%Fj( 0,j-1,k  )%S = tm*block%Fj(1,j-1,k  )%S
         block%Fj( 0,j  ,k-1)%S = tm*block%Fj(1,j  ,k-1)%S
         block%Fj( 0,j-1,k-1)%S = tm*block%Fj(1,j-1,k-1)%S

         block%Fk( 0,j  ,k  )%N = tm*block%Fk(1,j  ,k  )%N
         block%Fk( 0,j-1,k  )%N = tm*block%Fk(1,j-1,k  )%N
         block%Fk( 0,j  ,k-1)%N = tm*block%Fk(1,j  ,k-1)%N
         block%Fk( 0,j-1,k-1)%N = tm*block%Fk(1,j-1,k-1)%N
         block%Fk( 0,j  ,k  )%S = tm*block%Fk(1,j  ,k  )%S
         block%Fk( 0,j-1,k  )%S = tm*block%Fk(1,j-1,k  )%S
         block%Fk( 0,j  ,k-1)%S = tm*block%Fk(1,j  ,k-1)%S
         block%Fk( 0,j-1,k-1)%S = tm*block%Fk(1,j-1,k-1)%S
         ! volume
         block%C(  0,j,  k  )%V = block%C(    1,j,  k  )%V
      end if
    enddo
  enddo
  ! right i
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      correct = ((block%Fi(Ni+1,j,k)%BC%tp==bc_ext).OR. &
                 (block%Fi(Ni+1,j,k)%BC%tp==bc_ref).OR. &
                 (block%Fi(Ni+1,j,k)%BC%tp==bc_in1).OR. &
                 (block%Fi(Ni+1,j,k)%BC%tp==bc_in2).OR. &
                 (block%C(Ni+1,j,k)%V<(0.2_R_P*block%C(Ni,j,k)%V)))
      wall    = ((block%Fi(Ni+1,j,k)%BC%tp==bc_ext).OR.(block%Fi(Ni+1,j,k)%BC%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%Fi(Ni-1,j,k)%N*block%Fi(Ni-1,j,k)%S).dot.block%Fi(Ni,j,k)%N)
         block%Fi(Ni+1,j,  k  )%N = -(block%Fi(Ni-1,j,k)%N*block%Fi(Ni-1,j,k)%S) + sn*block%Fi(Ni,j,k)%N
         block%Fi(Ni+1,j,  k  )%S = normL2(   block%Fi(Ni+1,j,k)%N)
         block%Fi(Ni+1,j,  k  )%N = normalize(block%Fi(Ni+1,j,k)%N)
         ! tangential metrics
         block%Fj(Ni+1,j  ,k  )%N = tm*block%Fj(Ni,j  ,k  )%N
         block%Fj(Ni+1,j-1,k  )%N = tm*block%Fj(Ni,j-1,k  )%N
         block%Fj(Ni+1,j  ,k-1)%N = tm*block%Fj(Ni,j  ,k-1)%N
         block%Fj(Ni+1,j-1,k-1)%N = tm*block%Fj(Ni,j-1,k-1)%N
         block%Fj(Ni+1,j  ,k  )%S = tm*block%Fj(Ni,j  ,k  )%S
         block%Fj(Ni+1,j-1,k  )%S = tm*block%Fj(Ni,j-1,k  )%S
         block%Fj(Ni+1,j  ,k-1)%S = tm*block%Fj(Ni,j  ,k-1)%S
         block%Fj(Ni+1,j-1,k-1)%S = tm*block%Fj(Ni,j-1,k-1)%S

         block%Fk(Ni+1,j  ,k  )%N = tm*block%Fk(Ni,j  ,k  )%N
         block%Fk(Ni+1,j-1,k  )%N = tm*block%Fk(Ni,j-1,k  )%N
         block%Fk(Ni+1,j  ,k-1)%N = tm*block%Fk(Ni,j  ,k-1)%N
         block%Fk(Ni+1,j-1,k-1)%N = tm*block%Fk(Ni,j-1,k-1)%N
         block%Fk(Ni+1,j  ,k  )%S = tm*block%Fk(Ni,j  ,k  )%S
         block%Fk(Ni+1,j-1,k  )%S = tm*block%Fk(Ni,j-1,k  )%S
         block%Fk(Ni+1,j  ,k-1)%S = tm*block%Fk(Ni,j  ,k-1)%S
         block%Fk(Ni+1,j-1,k-1)%S = tm*block%Fk(Ni,j-1,k-1)%S
         ! volume
         block%C( Ni+1,j,  k  )%V = block%C(    Ni,j,  k  )%V
      end if
    enddo
  enddo
  ! left j
  !!$OMP DO
  do k=1,Nk
    do i=1,Ni
      correct = ((block%Fj(i,0,k)%BC%tp==bc_ext).OR. &
                 (block%Fj(i,0,k)%BC%tp==bc_ref).OR. &
                 (block%Fj(i,0,k)%BC%tp==bc_in1).OR. &
                 (block%Fj(i,0,k)%BC%tp==bc_in2).OR. &
                 (block%C(i,0,k)%V<(0.2_R_P*block%C(i,1,k)%V)))
      wall    = ((block%Fj(i,0,k)%BC%tp==bc_ext).OR.(block%Fj(i,0,k)%BC%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%Fj(i,1,k)%N*block%Fj(i,1,k)%S).dot.block%Fj(i,0,k)%N)
         block%Fj(i, -1,k  )%N = -(block%Fj(i,1,k)%N*block%Fj(i,1,k)%S) + sn*block%Fj(i,0,k)%N
         block%Fj(i, -1,k  )%S = normL2(   block%Fj(i,-1,k)%N)
         block%Fj(i, -1,k  )%N = normalize(block%Fj(i,-1,k)%N)
         ! tangential metrics
         block%Fi(i  ,0,k  )%N = tm*block%Fi(i  ,1,k  )%N
         block%Fi(i-1,0,k  )%N = tm*block%Fi(i-1,1,k  )%N
         block%Fi(i  ,0,k-1)%N = tm*block%Fi(i  ,1,k-1)%N
         block%Fi(i-1,0,k-1)%N = tm*block%Fi(i-1,1,k-1)%N
         block%Fi(i  ,0,k  )%S = tm*block%Fi(i  ,1,k  )%S
         block%Fi(i-1,0,k  )%S = tm*block%Fi(i-1,1,k  )%S
         block%Fi(i  ,0,k-1)%S = tm*block%Fi(i  ,1,k-1)%S
         block%Fi(i-1,0,k-1)%S = tm*block%Fi(i-1,1,k-1)%S

         block%Fk(i  ,0,k  )%N = tm*block%Fk(i  ,1,k  )%N
         block%Fk(i-1,0,k  )%N = tm*block%Fk(i-1,1,k  )%N
         block%Fk(i  ,0,k-1)%N = tm*block%Fk(i  ,1,k-1)%N
         block%Fk(i-1,0,k-1)%N = tm*block%Fk(i-1,1,k-1)%N
         block%Fk(i  ,0,k  )%S = tm*block%Fk(i  ,1,k  )%S
         block%Fk(i-1,0,k  )%S = tm*block%Fk(i-1,1,k  )%S
         block%Fk(i  ,0,k-1)%S = tm*block%Fk(i  ,1,k-1)%S
         block%Fk(i-1,0,k-1)%S = tm*block%Fk(i-1,1,k-1)%S
         ! volume
         block%C( i,  0,k  )%V = block%C(    i,  1,k  )%V
      end if
    enddo
  enddo
  ! right j
  !!$OMP DO
  do k=1,Nk
    do i=1,Ni
      correct = ((block%Fj(i,Nj+1,k)%BC%tp==bc_ext).OR. &
                 (block%Fj(i,Nj+1,k)%BC%tp==bc_ref).OR. &
                 (block%Fj(i,Nj+1,k)%BC%tp==bc_in1).OR. &
                 (block%Fj(i,Nj+1,k)%BC%tp==bc_in2).OR. &
                 (block%C(i,Nj+1,k)%V<(0.2_R_P*block%C(i,Nj,k)%V)))
      wall    = ((block%Fj(i,Nj+1,k)%BC%tp==bc_ext).OR.(block%Fj(i,Nj+1,k)%BC%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%Fj(i,Nj-1,k)%N*block%Fj(i,Nj-1,k)%S).dot.block%Fj(i,Nj,k)%N)
         block%Fj(i,Nj+1,  k  )%N = -(block%Fj(i,Nj-1,k)%N*block%Fj(i,Nj-1,k)%S) + sn*block%Fj(i,Nj,k)%N
         block%Fj(i,Nj+1,  k  )%S = normL2(   block%Fj(i,Nj+1,k)%N)
         block%Fj(i,Nj+1,  k  )%N = normalize(block%Fj(i,Nj+1,k)%N)
         ! tangential metrics
         block%Fi(i  ,Nj+1,k  )%N = tm*block%Fi(i  ,Nj,k  )%N
         block%Fi(i-1,Nj+1,k  )%N = tm*block%Fi(i-1,Nj,k  )%N
         block%Fi(i  ,Nj+1,k-1)%N = tm*block%Fi(i  ,Nj,k-1)%N
         block%Fi(i-1,Nj+1,k-1)%N = tm*block%Fi(i-1,Nj,k-1)%N
         block%Fi(i  ,Nj+1,k  )%S = tm*block%Fi(i  ,Nj,k  )%S
         block%Fi(i-1,Nj+1,k  )%S = tm*block%Fi(i-1,Nj,k  )%S
         block%Fi(i  ,Nj+1,k-1)%S = tm*block%Fi(i  ,Nj,k-1)%S
         block%Fi(i-1,Nj+1,k-1)%S = tm*block%Fi(i-1,Nj,k-1)%S

         block%Fk(i  ,Nj+1,k  )%N = tm*block%Fk(i  ,Nj,k  )%N
         block%Fk(i-1,Nj+1,k  )%N = tm*block%Fk(i-1,Nj,k  )%N
         block%Fk(i  ,Nj+1,k-1)%N = tm*block%Fk(i  ,Nj,k-1)%N
         block%Fk(i-1,Nj+1,k-1)%N = tm*block%Fk(i-1,Nj,k-1)%N
         block%Fk(i  ,Nj+1,k  )%S = tm*block%Fk(i  ,Nj,k  )%S
         block%Fk(i-1,Nj+1,k  )%S = tm*block%Fk(i-1,Nj,k  )%S
         block%Fk(i  ,Nj+1,k-1)%S = tm*block%Fk(i  ,Nj,k-1)%S
         block%Fk(i-1,Nj+1,k-1)%S = tm*block%Fk(i-1,Nj,k-1)%S
         ! volume
         block%C( i,  Nj+1,k  )%V = block%C(    i,  Nj,k  )%V
      end if
    enddo
  enddo
  ! left k
  !!$OMP DO
  do j=1,Nj
    do i=1,Ni
      correct = ((block%Fk(i,j,0)%BC%tp==bc_ext).OR. &
                 (block%Fk(i,j,0)%BC%tp==bc_ref).OR. &
                 (block%Fk(i,j,0)%BC%tp==bc_in1).OR. &
                 (block%Fk(i,j,0)%BC%tp==bc_in2).OR. &
                 (block%C(i,j,0)%V<(0.2_R_P*block%C(i,j,1)%V)))
      wall    = ((block%Fk(i,j,0)%BC%tp==bc_ext).OR.(block%Fk(i,j,0)%BC%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%Fk(i,j,1)%N*block%Fk(i,j,1)%S).dot.block%Fk(i,j,0)%N)
         block%Fk(i,  j, -1)%N = -(block%Fk(i,j,1)%N*block%Fk(i,j,1)%S) + sn*block%Fk(i,j,0)%N
         block%Fk(i,  j, -1)%S = normL2(   block%Fk(i,j,-1)%N)
         block%Fk(i,  j, -1)%N = normalize(block%Fk(i,j,-1)%N)
         ! tangential metrics
         block%Fi(i  ,j  ,0)%N = tm*block%Fi(i  ,j  ,1)%N
         block%Fi(i-1,j  ,0)%N = tm*block%Fi(i-1,j  ,1)%N
         block%Fi(i  ,j-1,0)%N = tm*block%Fi(i  ,j-1,1)%N
         block%Fi(i-1,j-1,0)%N = tm*block%Fi(i-1,j-1,1)%N
         block%Fi(i  ,j  ,0)%S = tm*block%Fi(i  ,j  ,1)%S
         block%Fi(i-1,j  ,0)%S = tm*block%Fi(i-1,j  ,1)%S
         block%Fi(i  ,j-1,0)%S = tm*block%Fi(i  ,j-1,1)%S
         block%Fi(i-1,j-1,0)%S = tm*block%Fi(i-1,j-1,1)%S

         block%Fj(i  ,j  ,0)%N = tm*block%Fj(i  ,j  ,1)%N
         block%Fj(i-1,j  ,0)%N = tm*block%Fj(i-1,j  ,1)%N
         block%Fj(i  ,j-1,0)%N = tm*block%Fj(i  ,j-1,1)%N
         block%Fj(i-1,j-1,0)%N = tm*block%Fj(i-1,j-1,1)%N
         block%Fj(i  ,j  ,0)%S = tm*block%Fj(i  ,j  ,1)%S
         block%Fj(i-1,j  ,0)%S = tm*block%Fj(i-1,j  ,1)%S
         block%Fj(i  ,j-1,0)%S = tm*block%Fj(i  ,j-1,1)%S
         block%Fj(i-1,j-1,0)%S = tm*block%Fj(i-1,j-1,1)%S
         ! volume
         block%C( i,  j,  0)%V = block%C(    i,  j,  1)%V
      end if
    enddo
  enddo
  ! right k
  !!$OMP DO
  do j=1,Nj
    do i=1,Ni
      correct = ((block%Fk(i,j,Nk+1)%BC%tp==bc_ext).OR. &
                 (block%Fk(i,j,Nk+1)%BC%tp==bc_ref).OR. &
                 (block%Fk(i,j,Nk+1)%BC%tp==bc_in1).OR. &
                 (block%Fk(i,j,Nk+1)%BC%tp==bc_in2).OR. &
                 (block%C(i,j,Nk+1)%V<(0.2_R_P*block%C(i,j,Nk)%V)))
      wall    = ((block%Fk(i,j,Nk+1)%BC%tp==bc_ext).OR.(block%Fk(i,j,Nk+1)%BC%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%Fk(i,j,Nk-1)%N*block%Fk(i,j,Nk-1)%S).dot.block%Fk(i,j,Nk)%N)
         block%Fk(i,  j,  Nk+1)%N = -(block%Fk(i,j,Nk-1)%N*block%Fk(i,j,Nk-1)%S) + sn*block%Fk(i,j,Nk)%N
         block%Fk(i,  j,  Nk+1)%S = normL2(   block%Fk(i,j,Nk+1)%N)
         block%Fk(i,  j,  Nk+1)%N = normalize(block%Fk(i,j,Nk+1)%N)
         ! tangential metrics
         block%Fi(i  ,j  ,Nk+1)%N = tm*block%Fi(i  ,j  ,Nk)%N
         block%Fi(i-1,j  ,Nk+1)%N = tm*block%Fi(i-1,j  ,Nk)%N
         block%Fi(i  ,j-1,Nk+1)%N = tm*block%Fi(i  ,j-1,Nk)%N
         block%Fi(i-1,j-1,Nk+1)%N = tm*block%Fi(i-1,j-1,Nk)%N
         block%Fi(i  ,j  ,Nk+1)%S = tm*block%Fi(i  ,j  ,Nk)%S
         block%Fi(i-1,j  ,Nk+1)%S = tm*block%Fi(i-1,j  ,Nk)%S
         block%Fi(i  ,j-1,Nk+1)%S = tm*block%Fi(i  ,j-1,Nk)%S
         block%Fi(i-1,j-1,Nk+1)%S = tm*block%Fi(i-1,j-1,Nk)%S

         block%Fj(i  ,j  ,Nk+1)%N = tm*block%Fj(i  ,j  ,Nk)%N
         block%Fj(i-1,j  ,Nk+1)%N = tm*block%Fj(i-1,j  ,Nk)%N
         block%Fj(i  ,j-1,Nk+1)%N = tm*block%Fj(i  ,j-1,Nk)%N
         block%Fj(i-1,j-1,Nk+1)%N = tm*block%Fj(i-1,j-1,Nk)%N
         block%Fj(i  ,j  ,Nk+1)%S = tm*block%Fj(i  ,j  ,Nk)%S
         block%Fj(i-1,j  ,Nk+1)%S = tm*block%Fj(i-1,j  ,Nk)%S
         block%Fj(i  ,j-1,Nk+1)%S = tm*block%Fj(i  ,j-1,Nk)%S
         block%Fj(i-1,j-1,Nk+1)%S = tm*block%Fj(i-1,j-1,Nk)%S
         ! volume
         block%C( i,  j,  Nk+1)%V = block%C(    i,  j,  Nk)%V
      end if
    enddo
  enddo
  !!$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine metrics_correction

  !> Subroutine for computing cell center coordinates from cell nodes ones.
  subroutine node2center(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block  !< Block-level data.
  integer(I_P)::                      i,j,k !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(block)
  !$OMP DO
  do k=1-block%gc(5),block%Nk+block%gc(6)
    do j=1-block%gc(3),block%Nj+block%gc(4)
      do i=1-block%gc(1),block%Ni+block%gc(2)
        block%C(i,j,k)%cent%x = (block%node(i,  j,  k  )%x + &
                                 block%node(i-1,j,  k  )%x + &
                                 block%node(i  ,j-1,k  )%x + &
                                 block%node(i  ,j  ,k-1)%x + &
                                 block%node(i-1,j-1,k-1)%x + &
                                 block%node(i  ,j-1,k-1)%x + &
                                 block%node(i-1,j  ,k-1)%x + &
                                 block%node(i-1,j-1,k  )%x)*0.125_R_P
        block%C(i,j,k)%cent%y = (block%node(i,  j,  k  )%y + &
                                 block%node(i-1,j,  k  )%y + &
                                 block%node(i  ,j-1,k  )%y + &
                                 block%node(i  ,j  ,k-1)%y + &
                                 block%node(i-1,j-1,k-1)%y + &
                                 block%node(i  ,j-1,k-1)%y + &
                                 block%node(i-1,j  ,k-1)%y + &
                                 block%node(i-1,j-1,k  )%y)*0.125_R_P
        block%C(i,j,k)%cent%z = (block%node(i,  j,  k  )%z + &
                                 block%node(i-1,j,  k  )%z + &
                                 block%node(i  ,j-1,k  )%z + &
                                 block%node(i  ,j  ,k-1)%z + &
                                 block%node(i-1,j-1,k-1)%z + &
                                 block%node(i  ,j-1,k-1)%z + &
                                 block%node(i-1,j  ,k-1)%z + &
                                 block%node(i-1,j-1,k  )%z)*0.125_R_P
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine node2center
  !> @}
endmodule Data_Type_SBlock
