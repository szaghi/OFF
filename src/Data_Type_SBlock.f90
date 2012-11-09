!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_SBlockPublicProcedure Data_Type_SBlock
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
USE Data_Type_Conservative !< Definition of Type_Conservative.
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
!> @brief Derived type containing the structured block-level data.
!> Structured block-level type contains data (mesh, boundary conditions and fluid dynamic data) of structured (implicit
!> connectivity) numerical grid.
!> @ingroup DerivedType
type, public:: Type_SBlock
  integer(I_P):: myrank = 0_I_P !< Rank of the process which block belongs to.
  ! mesh data
  integer(I1P)::                   gc(1:6)=&
                                          (/1_I1P, & ! gc(1) => left  i.
                                            1_I1P, & ! gc(2) => right i.
                                            1_I1P, & ! gc(3) => left  j.
                                            1_I1P, & ! gc(4) => right j.
                                            1_I1P, & ! gc(5) => left  k.
                                            1_I1P  & ! gc(6) => right k.
                                            /)       !< Number of ghost cells for the 6 faces of the block.
  integer(I_P)::                   Ni     = 0_I_P    !< Number of cells in i direction.
  integer(I_P)::                   Nj     = 0_I_P    !< Number of cells in j direction.
  integer(I_P)::                   Nk     = 0_I_P    !< Number of cells in k direction.
  type(Type_Vector), allocatable:: node(:,:,:)       !< Cell nodes coordinates [0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)].
  type(Type_Vector), allocatable:: NFi(:,:,:)        !< Face i normals, versor [0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Vector), allocatable:: NFj(:,:,:)        !< Face j normals, versor [1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Vector), allocatable:: NFk(:,:,:)        !< Face k normals, versor [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)].
  real(R_P),         allocatable:: Si(:,:,:)         !< Face i area            [0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  real(R_P),         allocatable:: Sj(:,:,:)         !< Face j area            [1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  real(R_P),         allocatable:: Sk(:,:,:)         !< Face k area            [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)].
  real(R_P),         allocatable:: V(:,:,:)          !< Cell volume            [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Cell),   allocatable:: cell(:,:,:)       !< Cell data informations [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Vector), allocatable:: cent(:,:,:)       !< Cell center coordinates[1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  ! boundary conditions data
  type(Type_BC), allocatable:: BCi(:,:,:) !< Boundary conditions of i faces [0-gc(1):1+gc(2),1-gc(3):1+gc(4),1-gc(5):1+gc(6)].
  type(Type_BC), allocatable:: BCj(:,:,:) !< Boundary conditions of j faces [1-gc(1):1+gc(2),0-gc(3):1+gc(4),1-gc(5):1+gc(6)].
  type(Type_BC), allocatable:: BCk(:,:,:) !< Boundary conditions of k faces [1-gc(1):1+gc(2),1-gc(3):1+gc(4),0-gc(5):1+gc(6)].
  ! fluid dynamic data
  real(R_P),               allocatable:: Dt(:,:,:)   !< Local time step.
  type(Type_Primitive),    allocatable:: P (:,:,:)   !< Primitive variables.
  type(Type_Conservative), allocatable:: U (:,:,:)   !< Conservative variables.
  type(Type_Conservative), allocatable:: KS(:,:,:,:) !< Runge-Kutta stages of conservative variables [1:rk_ord].
  contains
    procedure                 :: free => free_block                    ! Procedure for freeing memory.
    procedure                 :: alloc => alloc_block                  ! Procedure for allocating memory.
    procedure, non_overridable:: save_mesh => save_bmesh               ! Procedure for saving the mesh data.
    procedure, non_overridable:: load_mesh_dims => load_bmesh_dims     ! Procedure for loading the mesh data dimensions.
    procedure, non_overridable:: load_mesh => load_bmesh               ! Procedure for loading the mesh data.
    procedure, non_overridable:: print_info_mesh => print_info_bmesh   ! Procedure for printing of the mesh data.
    procedure, non_overridable:: save_bc => save_bbc                   ! Procedure for saving the bc data.
    procedure, non_overridable:: load_bc => load_bbc                   ! Procedure for loading the bc data.
    procedure, non_overridable:: save_fluid => save_bfluid             ! Procedure for saving the fluidynamic data.
    procedure, non_overridable:: load_fluid => load_bfluid             ! Procedure for loading the fluidynamic data.
    procedure, non_overridable:: print_info_fluid => print_info_bfluid ! Procedure for printing of the mesh data.
endtype Type_SBlock
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_SBlockPrivateProcedure
  !> @{
  !> Subroutine for freeing dynamic data of Type_SBlock variables.
  subroutine free_block(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block !< Block data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! mesh data
  if (allocated(block%node)) deallocate(block%node)
  if (allocated(block%NFi )) deallocate(block%NFi )
  if (allocated(block%NFj )) deallocate(block%NFj )
  if (allocated(block%NFk )) deallocate(block%NFk )
  if (allocated(block%Si  )) deallocate(block%Si  )
  if (allocated(block%Sj  )) deallocate(block%Sj  )
  if (allocated(block%Sk  )) deallocate(block%Sk  )
  if (allocated(block%V   )) deallocate(block%V   )
  if (allocated(block%cent)) deallocate(block%cent)
  if (allocated(block%cell)) deallocate(block%cell)
  ! boundary conditions data
  if (allocated(block%BCi)) then
    call block%BCi%free ; deallocate(block%BCi)
  endif
  if (allocated(block%BCj)) then
    call block%BCj%free ; deallocate(block%BCj)
  endif
  if (allocated(block%BCk)) then
    call block%BCk%free ; deallocate(block%BCk)
  endif
  ! fluid dynamic data
  if (allocated(block%Dt)) deallocate(block%Dt)
  if (allocated(block%P)) then
    call block%P%free ; deallocate(block%P)
  endif
  if (allocated(block%U)) then
    call block%U%free ; deallocate(block%U)
  endif
  if (allocated(block%KS)) then
    call block%KS%free ; deallocate(block%KS)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_block

  !> Subroutine for allocating dynamic data of Type_SBlock variables.
  subroutine alloc_block(block,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block data.
  type(Type_Global),  intent(IN)::    global   !< Global data.
  integer(I_P)::                      gc(1:6)  !< Temporary variable  for storing block ghost cells number.
  integer(I_P)::                      Ni,Nj,Nk !< Temporary variables for storing block dimensions.
  integer(I_P)::                      Ns       !< Temporary variable  for storing number of species.
  integer(I_P)::                      rk_ord   !< Temporary variable  for storing rk_ord.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! free dynamic data of block if previously allocated
  call block%free
  ! storing block dimensions into temporary variables to simplify the code
  gc(1:6) = block%gc(1:6)
  Ni      = block%Ni
  Nj      = block%Nj
  Nk      = block%Nk
  Ns      = global%Ns
  rk_ord  = global%rk_ord
  ! mesh data
  allocate(block%node(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  allocate(block%NFi (0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%NFj (1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%NFk (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  allocate(block%Si  (0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%Si = 0._R_P
  allocate(block%Sj  (1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%Sj = 0._R_P
  allocate(block%Sk  (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; block%Sk = 0._R_P
  allocate(block%V   (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%V  = 0._R_P
  allocate(block%cell(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%cent(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  ! boundary conditions data
  allocate(block%BCi(0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%BCj(1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%BCk(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  ! fluid dynamic data
  allocate(block%Dt(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; block%Dt = 0._R_P
  allocate(block%P (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; call block%P%init(Ns=Ns)
  allocate(block%U (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; call block%U%init(Ns=Ns)
  allocate(block%KS(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6),1:rk_ord)) ; call block%KS%init(Ns=Ns)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_block

  !> Function for saving block mesh file.
  !> @return \b err integer(I4P) variable.
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
    err = write_cell(array3D=block%cell,format='*',unit=UnitFree)
    close(UnitFree)
  else
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), action = 'WRITE', form = 'UNFORMATTED')
    write(UnitFree,iostat=err)block%gc
    write(UnitFree,iostat=err)block%Ni,block%Nj,block%Nk
    err = write_vector(array3D=block%node,unit=UnitFree)
    err = write_cell(array3D=block%cell,unit=UnitFree)
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bmesh

  !> Function for loading the mesh data dimensions of block from the mesh file "filename".
  !> @return \b err integer(I4P) variable.
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
  if (.NOT.is_file) call File_Not_Found(block%myrank,filename,'load_bmesh_dims')
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
  !> @return \b err integer(I4P) variable.
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
  if (.NOT.is_file) call File_Not_Found(block%myrank,filename,'load_bmesh')
  if (present(ascii)) then
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err)block%gc
    read(UnitFree,*,iostat=err)block%Ni,block%Nj,block%Nk
    err = read_vector(array3D=block%node,format='*',unit=UnitFree)
    err = read_cell(array3D=block%cell,format='*',unit=UnitFree)
    close(UnitFree)
  else
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)block%gc
    read(UnitFree,iostat=err)block%Ni,block%Nj,block%Nk
    err = read_vector(array3D=block%node,unit=UnitFree)
    err = read_cell(array3D=block%cell,unit=UnitFree)
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bmesh

  !> Function for printing to standard output info of block mesh data.
  !> @return \b err integer(I4P) variable.
  function print_info_bmesh(block,blk,grl,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(IN):: block       !< Block level data.
  integer(I_P),       intent(IN):: blk         !< Actual block number.
  integer(I_P),       intent(IN):: grl         !< Actual grid level number.
  type(Type_Global),  intent(IN):: global      !< Global level data.
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
  rks = 'rank'//trim(str(.true.,block%myrank))
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)trim(rks)//' Mesh infos of myrank '//trim(rks)
  write(stdout,'(A)',iostat=err)trim(rks)//'  Global number of blocks '//trim(str(.true.,global%Nb_tot))
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
  if (allocated(NFi)) deallocate(NFi) ; allocate(NFi(0:Ni,1:Nj,1:Nk)) ; NFi(0:Ni,1:Nj,1:Nk) = block%NFi(0:Ni,1:Nj,1:Nk)
  if (allocated(NFj)) deallocate(NFj) ; allocate(NFj(1:Ni,0:Nj,1:Nk)) ; NFj(1:Ni,0:Nj,1:Nk) = block%NFj(1:Ni,0:Nj,1:Nk)
  if (allocated(NFk)) deallocate(NFk) ; allocate(NFk(1:Ni,1:Nj,0:Nk)) ; NFk(1:Ni,1:Nj,0:Nk) = block%NFk(1:Ni,1:Nj,0:Nk)
  if (allocated(Si )) deallocate(Si ) ; allocate(Si (0:Ni,1:Nj,1:Nk)) ; Si (0:Ni,1:Nj,1:Nk) = block%Si (0:Ni,1:Nj,1:Nk)
  if (allocated(Sj )) deallocate(Sj ) ; allocate(Sj (1:Ni,0:Nj,1:Nk)) ; Sj (1:Ni,0:Nj,1:Nk) = block%Sj (1:Ni,0:Nj,1:Nk)
  if (allocated(Sk )) deallocate(Sk ) ; allocate(Sk (1:Ni,1:Nj,0:Nk)) ; Sk (1:Ni,1:Nj,0:Nk) = block%Sk (1:Ni,1:Nj,0:Nk)
  if (allocated(V  )) deallocate(V  ) ; allocate(V  (1:Ni,1:Nj,1:Nk)) ; V  (1:Ni,1:Nj,1:Nk) = block%V  (1:Ni,1:Nj,1:Nk)
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
  !> @return \b err integer(I4P) variable.
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
  err = write_bc(array3D=block%BCi,unit=UnitFree)
  err = write_bc(array3D=block%BCj,unit=UnitFree)
  err = write_bc(array3D=block%BCk,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bbc

  !> Function for loading block boundary conditions file.
  !> @return \b err integer(I4P) variable.
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
  if (.NOT.is_file) call File_Not_Found(block%myrank,filename,'load_bbc')
  open(unit=Get_Unit(UnitFree), file=adjustl(trim(filename)), status='OLD', action='READ', form='UNFORMATTED')
  err = read_bc(array3D=block%BCi,unit=UnitFree)
  err = read_bc(array3D=block%BCj,unit=UnitFree)
  err = read_bc(array3D=block%BCk,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bbc

  !> Function for saving fluid dynamic file.
  !> @return \b err integer(I4P) variable.
  function save_bfluid(block,filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(IN):: block    !< Block level data.
  character(*),       intent(IN):: filename !< Name of file where mesh variables are saved.
  type(Type_Global),  intent(IN):: global   !< Global level data.
  integer(I_P)::                   err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                   UnitFree !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), action = 'WRITE', form = 'UNFORMATTED')
  write(UnitFree,iostat=err)global%Ns
  write(UnitFree,iostat=err)global%cp0,global%cv0
  write(UnitFree,iostat=err)global%n
  write(UnitFree,iostat=err)global%t
  write(UnitFree,iostat=err)block%Dt
  err = write_primitive(array3D=block%P,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bfluid

  !> Function for loading fluid dynamic file.
  !> @return \b err integer(I4P) variable.
  function load_bfluid(block,filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block level data.
  character(*),       intent(IN)::    filename !< Name of file where mesh variables are saved.
  type(Type_Global),  intent(INOUT):: global   !< Global level data.
  integer(I_P)::                      err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                      UnitFree !< Free logic unit.
  logical::                           is_file  !< Flag for inquiring the presence of fluid dynamic file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(block%myrank,filename,'load_bfluid')
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
  read(UnitFree,iostat=err)global%Ns
  read(UnitFree,iostat=err)global%cp0,global%cv0
  read(UnitFree,iostat=err)global%n
  read(UnitFree,iostat=err)global%t
  read(UnitFree,iostat=err)block%Dt
  err = read_primitive(array3D=block%P,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bfluid

  !> Function for printing to standard output info of block fluid data.
  !> @return \b err integer(I4P) variable.
  function print_info_bfluid(block,blk,grl,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(IN):: block          !< Block level data.
  integer(I_P),       intent(IN):: blk            !< Actual block number.
  integer(I_P),       intent(IN):: grl            !< Actual grid level number.
  type(Type_Global),  intent(IN):: global         !< Global level data.
  integer(I_P)::                   err            !< Error trapping flag: 0 no errors, >0 error occurs.
  character(DR_P)::                vmax,vmin      !< String for printing max and min of variables.
  real(R_P), allocatable::         dummy(:,:,:,:) !<  Dummy variables for printing only internal cells info.
  integer(I_P)::                   Ni,Nj,Nk,Np,Nc !< Temporary variables for storing block dimensions.
  character(DI_P)::                rks            !< String containing myrank.
  integer(I_P)::                   i,j,k,s        !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,block%myrank))
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)trim(rks)//' Fluid infos of myrank '//trim(rks)
  write(stdout,'(A)',iostat=err)trim(rks)//'  Global number of blocks '//trim(str(.true.,global%Nb_tot))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Number of initial species '//trim(str(.true.,global%Ns))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Block '//trim(str(.true.,blk))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Grid level '//trim(str(.true.,grl))
  Ni = block%Ni
  Nj = block%Nj
  Nk = block%Nk
  Np = global%Np
  Nc = global%Nc
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Nc,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(:,i,j,k) = block%U(i,j,k)%cons2array()
      enddo
    enddo
  enddo
  write(stdout,'(A)',iostat=err)trim(rks)//'  Conservative variables'
  do s=1,global%Nc
    vmax=trim(str(n=maxval(dummy(s,:,:,:)))) ; vmin=trim(str(n=minval(dummy(s,:,:,:))))
    write(stdout,'(A)',iostat=err)trim(rks)//'  U('//trim(str(.true.,s))//'): max '//trim(vmax)//' min '//trim(vmin)
  enddo
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Np,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(:,i,j,k) = block%P(i,j,k)%prim2array()
      enddo
    enddo
  enddo
  write(stdout,'(A)',iostat=err)trim(rks)//'  Primitive variables'
  do s=1,global%Np
    vmax=trim(str(n=maxval(dummy(s,:,:,:)))) ; vmin=trim(str(n=minval(dummy(s,:,:,:))))
    write(stdout,'(A)',iostat=err)trim(rks)//'  P('//trim(str(.true.,s))//'): max '//trim(vmax)//' min '//trim(vmin)
  enddo
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(1,i,j,k) = block%Dt(i,j,k)
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
  !> @}
endmodule Data_Type_SBlock
