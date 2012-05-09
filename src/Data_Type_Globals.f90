!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_GlobalsPrivateProcedure Data_Type_Globals
!> @}

!> @brief Module Data_Type_Globals contains the definition of global types variables and useful procedures for their handling.
!> Two main derived type are defined: \n
!> - Type_Global: derived type containing global-level data.
!> - Type_Block: derived type containing block-level data.
module Data_Type_Globals
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                            !< Integers and reals precision definition.
USE Data_Type_BC, init_bc=>init, set_bc=>set                !< Definition of Type_BC.
USE Data_Type_Cell, init_cell=>init, set_cell=>set          !< Definition of Type_Cell.
USE Data_Type_Conservative, init_cons=>init, set_cons=>set  !< Definition of Type_Conservative.
USE Data_Type_Primitive, init_prim=>init, set_prim=>set     !< Definition of Type_Primitive.
USE Data_Type_Vector, set_vec=>set                          !< Definition of Type_Vector.
USE Lib_IO_Misc                                             !< Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: file_name
public:: alloc_global_bc,alloc_global_fluid
public:: free_block,alloc_block
public:: load_bmesh_dims,load_bmesh,save_bmesh,print_info_bmesh
public:: load_gbc_in1,load_bbc,save_bbc
public:: load_gfluid_soption,load_gfluid_Ns,load_gfluid_0species,load_bfluid,save_bfluid,print_info_bfluid
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the global-level file data.
!> Global-level file data are referred to those informations concerning with files of global interest.
!> @ingroup DerivedType
type, public:: Type_File_Global
  character(60)::  Path_InPut           !< Path of input files.
  character(60)::  Path_OutPut          !< Path of output files.
  character(60)::  File_BC     ='unset' !< Base name of boundary conditions file.
  character(60)::  File_Init   ='unset' !< Base name of initial conditions file.
  character(60)::  File_Mesh   ='unset' !< Base name of mesh file.
  character(60)::  File_Spec   ='unset' !< Name of initial species file.
  character(60)::  File_Sol    ='unset' !< Base name of solution file.
  character(60)::  File_Solver ='unset' !< Name of solver options file.
  character(60)::  File_Pout   ='unset' !< Post-processed output file name.
  character(500):: varform_res          !< Gnuplot residuals writing format.
  integer(I8P)::   screen_out  = 1      !< Console refresh frequency.
  integer(I8P)::   sol_out     = 0      !< Actual solution writing frequency (if 0 only restart solution is saved).
  integer(I8P)::   restart_out = 1      !< Restart Solution writing frequency.
  integer(I8P)::   probe_out   = 1      !< Probes writing frequency.
  integer(I_P)::   unit_res             !< Logical unit of gnuplot log file of residuals.
endtype Type_File_Global

!> @brief Derived type containing the global-level mesh data.
!> Global-level mesh data are referred to those informations concerning with mesh's details of global interest.
!> @ingroup DerivedType
type, public:: Type_Mesh_Global
  integer(I_P):: Nl     = 1_I_P !< Number of grid levels.
  integer(I_P):: Nb     = 0_I_P !< Number of blocks.
  integer(I_P):: Nb_tot = 0_I_P !< Number of total blocks (sum over each process).
  integer(I1P):: gco    = 1_I_P !< Number of ghost cells necessary to achieve the space reconstruction order.
endtype Type_Mesh_Global
!> @brief Derived type containing the block-level mesh data.
!> Block-level mesh data are referred to those informations concerning with mesh's details of blocks interest.
type, public:: Type_Mesh_Block
  integer(I_P)::                   Ni     = 0_I_P    !< Number of cells in i direction.
  integer(I_P)::                   Nj     = 0_I_P    !< Number of cells in j direction.
  integer(I_P)::                   Nk     = 0_I_P    !< Number of cells in k direction.
  integer(I1P)::                   gc(1:6)=&
                                          (/1_I1P, & ! gc(1) => left  i.
                                            1_I1P, & ! gc(2) => right i.
                                            1_I1P, & ! gc(3) => left  j.
                                            1_I1P, & ! gc(4) => right j.
                                            1_I1P, & ! gc(5) => left  k.
                                            1_I1P  & ! gc(6) => right k.
                                            /)       !< Number of ghost cells for the 6 faces of the block.
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
  ! Adaptive Mesh Refinement
  integer(I1P)::                   octree = 0_I_P    !< Octree refinement level.
  type(Type_Mesh_Block), pointer:: amr(:,:,:)=>null()!< Mesh data of refined cells.
endtype Type_Mesh_Block

!> @brief Derived type containing the global-level boundary conditions data.
!> Global-level boundary conditions data are referred to those informations concerning with bc's details of global interest.
!> @ingroup DerivedType
type, public:: Type_BC_Global
  integer(I_P)::                      Nin1 = 0_I_P !< Number of inflow 1 boundary conditions.
  type(Type_Primitive), allocatable:: in1(:)       !< Inflow 1 boundary conditions primitive variables [1:Nin1].
endtype Type_BC_Global
!> @brief Derived type containing the block-level boundary conditions data.
!> Block-level boundary conditions data are referred to those informations concerning with bc's details of blocks interest.
!> @ingroup DerivedType
type, public:: Type_BC_Block
  type(Type_BC), allocatable:: BCi(:,:,:) !< Boundary conditions of i faces [0-gc(1):1+gc(2),1-gc(3):1+gc(4),1-gc(5):1+gc(6)].
  type(Type_BC), allocatable:: BCj(:,:,:) !< Boundary conditions of j faces [1-gc(1):1+gc(2),0-gc(3):1+gc(4),1-gc(5):1+gc(6)].
  type(Type_BC), allocatable:: BCk(:,:,:) !< Boundary conditions of k faces [1-gc(1):1+gc(2),1-gc(3):1+gc(4),0-gc(5):1+gc(6)].
endtype Type_BC_Block

!> Derived type containing the global-level fluid dynamic data.
!> Global-level fluid dynamic data are referred to those informations concerning with fluids details of global interest.
!> @ingroup DerivedType
type, public:: Type_Fluid_Global
  integer(I8P)::           n             = 0_I8P    !< Time steps counter.
  real(R_P)::              t             = 0._R_P   !< Time.
  integer(I_P)::           Ns            = 1_I_P    !< Number of species.
  integer(I_P)::           Np            = 7_I_P    !< Number of primitive variables    (Np = Ns + 6).
  integer(I_P)::           Nc            = 5_I_P    !< Number of conservative variables (Nc = Ns + 4).
  logical::                unsteady      = .true.   !< Type of simulation: unsteady or not.
  integer(I8P)::           Nmax          = 0_I8P    !< Max number of iterations.
  real(R_P)::              Tmax          = 0._R_P   !< Max time, ignored if Nmax>0.
  integer(I1P)::           sp_ord        = 1_I_P    !< Order of space convergence (number of ghost cells).
  integer(I1P)::           rk_ord        = 1_I_P    !< Order of time convergence (number of Runge-Kutta stages).
  real(R_P)::              CFL           = 0.3_R_P  !< Value of stability coefficient.
  real(R_P)::              residual_toll = 0.01_R_P !< Tolerance for residuals vanishing evaluation.
  logical::                residual_stop = .false.  !< Sentinel for stopping steady simulation when residuals vanish.
  real(R_P), allocatable:: cp0(:)                   !< Initial specific heat cp for each specie [1:Ns].
  real(R_P), allocatable:: cv0(:)                   !< Initial specific heat cv for each specie [1:Ns].
endtype Type_Fluid_Global
!> Derived type containing the block-level fluid dynamic data.
!> Block-level fluid dynamic data are referred to those informations concerning with fluids details of blocks interest.
!> @note Dimensions of 3D array are [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
!> @ingroup DerivedType
type, public:: Type_Fluid_Block
  real(R_P),               allocatable:: Dt(:,:,:)   !< Local time step.
  type(Type_Primitive),    allocatable:: P (:,:,:)   !< Primitive variables.
  type(Type_Conservative), allocatable:: U (:,:,:)   !< Conservative variables.
  type(Type_Conservative), allocatable:: KS(:,:,:,:) !< Runge-Kutta stages of conservative variables [1:rk_ord].
endtype Type_Fluid_Block

!> @brief Derived type containing the global-level data.
!> @ingroup DerivedType
type, public:: Type_Global
  type(Type_File_Global)::  file  !< File data.
  type(Type_Mesh_Global)::  mesh  !< Mesh data.
  type(Type_BC_Global)::    bc    !< Boundary conditions data.
  type(Type_Fluid_Global):: fluid !< Fluid dynamic data.
endtype Type_Global
!> Derived type containing the block-level data.
!> @ingroup DerivedType
type, public:: Type_Block
  type(Type_Mesh_Block)::  mesh  !< Mesh data.
  type(Type_BC_Block)::    bc    !< Boundary conditions data.
  type(Type_Fluid_Block):: fluid !< Fluid dynamic data.
endtype Type_Block
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Generic procedure for building file names.
!> Three different types of output file name can be build by means of this generic interface: \n
!> - File name of \b block-clean data; this type of file name is referred to those files containing block-level data; calling
!>   signature is of the type:
!>   @code fname = file_name(basename,suffix,blk,grl) @endcode
!> - File name of \b block-time_step data; this type of file name is referred to those files containing block-level data that
!>   varying with time; calling signature is of the type:
!>   @code fname = file_name(basename,suffix,blk,grl,n) @endcode
!> - File name of \b block-flip_flop data; this type of file name is referred to those files containing block-level data stored as
!>   backup flip/flop file; calling signature is of the type:
!>   @code fname = file_name(basename,suffix,blk,grl,flip) @endcode
!> @ingroup Interface
interface file_name
  module procedure Block_File_Name,Block_Step_File_Name,Block_Flip_File_Name
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_GlobalsPrivateProcedure
  !> @{
  !> Function for building block-clean file name.
  !> @return \b filename character(len_trim(basename)+4+5+len_trim(suffix)) variable.
  function Block_File_Name(basename,suffix,blk,grl) result(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN)::                           basename !< Base name of file name.
  character(*), intent(IN)::                           suffix   !< Suffix   of file name.
  integer(I_P), intent(IN)::                           blk      !< Block number.
  integer(I_P), intent(IN)::                           grl      !< Grid refinement level.
  character(len_trim(basename)+4+5+len_trim(suffix)):: filename !< Output file name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  filename = trim(basename)//'.g'//trim(strz(2,grl))//'.b'//trim(strz(3,blk))//trim(suffix)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Block_File_Name

  !> Function for building block-time_step file name.
  !> @return \b filename character(len_trim(basename)+4+5+13+len_trim(suffix)) variable.
  function Block_Step_File_Name(basename,suffix,blk,grl,n) result(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN)::                              basename !< Base name of file name.
  character(*), intent(IN)::                              suffix   !< Suffix   of file name.
  integer(I_P), intent(IN)::                              blk      !< Block number.
  integer(I_P), intent(IN)::                              grl      !< Grid refinement level.
  integer(I8P), intent(IN)::                              n        !< Time step number.
  character(len_trim(basename)+4+5+13+len_trim(suffix)):: filename !< Output file name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  filename = trim(basename)//'.g'//trim(strz(2,grl))//'.b'//trim(strz(3,blk))//'-N_'//trim(strz(10,n))//trim(suffix)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Block_Step_File_Name

  !> Function for building block-flip_flop file name.
  !> @return \b filename character(len_trim(basename)+3+4+5+len_trim(suffix)) variable.
  function Block_Flip_File_Name(basename,suffix,blk,grl,flip) result(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN)::                             basename !< Base name of file name.
  character(*), intent(IN)::                             suffix   !< Suffix   of file name.
  integer(I_P), intent(IN)::                             blk      !< Block number.
  integer(I_P), intent(IN)::                             grl      !< Grid refinement level.
  integer(I1P), intent(IN)::                             flip     !< Flip-flop number.
  character(len_trim(basename)+3+4+5+len_trim(suffix)):: filename !< Output file name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  filename = trim(basename)//'-f'//trim(strz(1,flip))//'.g'//trim(strz(2,grl))//'.b'//trim(strz(3,blk))//trim(suffix)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Block_Flip_File_Name
  !> @}

  !> Subroutine for allocating dynamic data of Type_Global boundary conditions variables.
  subroutine alloc_global_bc(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(INOUT):: global !< Global data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(global%bc%in1)) deallocate(global%bc%in1) ; if (global%bc%Nin1>0) allocate(global%bc%in1(1:global%bc%Nin1))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_global_bc

  !> Subroutine for allocating dynamic data of Type_Global fluid dynamic variables.
  subroutine alloc_global_fluid(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(INOUT):: global !< Global data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(global%fluid%cp0)) deallocate(global%fluid%cp0) ; allocate(global%fluid%cp0(1:global%fluid%Ns))
  if (allocated(global%fluid%cv0)) deallocate(global%fluid%cv0) ; allocate(global%fluid%cv0(1:global%fluid%Ns))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_global_fluid

  !> Subroutine for free dynamic data of Type_Block variables.
  subroutine free_block(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_block), intent(INOUT):: block !< Block data.
  integer(I_P)::                    err   !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Mesh data
  if (allocated(block%mesh%node)) deallocate(block%mesh%node)
  if (allocated(block%mesh%NFi )) deallocate(block%mesh%NFi )
  if (allocated(block%mesh%NFj )) deallocate(block%mesh%NFj )
  if (allocated(block%mesh%NFk )) deallocate(block%mesh%NFk )
  if (allocated(block%mesh%Si  )) deallocate(block%mesh%Si  )
  if (allocated(block%mesh%Sj  )) deallocate(block%mesh%Sj  )
  if (allocated(block%mesh%Sk  )) deallocate(block%mesh%Sk  )
  if (allocated(block%mesh%V   )) deallocate(block%mesh%V   )
  if (allocated(block%mesh%cent)) deallocate(block%mesh%cent)
  if (allocated(block%mesh%cell)) deallocate(block%mesh%cell)
  ! Boundary conditions data
  if (allocated(block%bc%BCi)) then
    err = free(block%bc%BCi) ! free the dynamic memory of Type_BC data
    deallocate(block%bc%BCi) ! free the block dynamic data
  endif
  if (allocated(block%bc%BCj)) then
    err = free(block%bc%BCj) ! free the dynamic memory of Type_BC data
    deallocate(block%bc%BCj) ! free the block dynamic data
  endif
  if (allocated(block%bc%BCk)) then
    err = free(block%bc%BCk) ! free the dynamic memory of Type_BC data
    deallocate(block%bc%BCk) ! free the block dynamic data
  endif
  ! Fluid dynamic data
  if (allocated(block%fluid%Dt)) deallocate(block%fluid%Dt)
  if (allocated(block%fluid%P)) then
    err = free(block%fluid%P) ! free the dynamic memory of Type_Primitive data
    deallocate(block%fluid%P) ! free the block dynamic data
  endif
  if (allocated(block%fluid%U)) then
    err = free(block%fluid%U) ! free the dynamic memory of Type_Conservative data
    deallocate(block%fluid%U) ! free the block dynamic data
  endif
  if (allocated(block%fluid%KS)) then
    err = free(block%fluid%KS) ! free the dynamic memory of Type_Conservative data
    deallocate(block%fluid%KS) ! free the block dynamic data
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_block

  !> Subroutine for allocating dynamic data of Type_Block variables.
  subroutine alloc_block(global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(IN)::    global   !< Global data.
  type(Type_block),  intent(INOUT):: block    !< Block data.
  integer(I_P)::                     Ni,Nj,Nk !< Temporary variables for storing blocks dimensions.
  integer(I_P)::                     gc(1:6)  !< Temporary variable  for storing blocks ghost cells number.
  integer(I_P)::                     Ns       !< Temporary variable  for storing number of species.
  integer(I_P)::                     rk_ord   !< Temporary variable  for storing rk_ord.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Free dynamic data block if previously allocated
  call free_block(block)
  ! Storing block dimensions into temporary variables to simplify the code
  Ni      = block%mesh%Ni
  Nj      = block%mesh%Nj
  Nk      = block%mesh%Nk
  gc(1:6) = block%mesh%gc(1:6)
  Ns      = global%fluid%Ns
  rk_ord  = global%fluid%rk_ord
  ! Mesh data
  allocate(block%mesh%node(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; block%mesh%node=0._R_P
  allocate(block%mesh%NFi (0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%mesh%NFi =0._R_P
  allocate(block%mesh%NFj (1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%mesh%NFj =0._R_P
  allocate(block%mesh%NFk (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; block%mesh%NFk =0._R_P
  allocate(block%mesh%Si  (0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%mesh%Si  =0._R_P
  allocate(block%mesh%Sj  (1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%mesh%Sj  =0._R_P
  allocate(block%mesh%Sk  (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; block%mesh%Sk  =0._R_P
  allocate(block%mesh%V   (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%mesh%V   =0._R_P
  allocate(block%mesh%cell(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%mesh%cell=init_cell()
  allocate(block%mesh%cent(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%mesh%cent=0._R_P
  ! Boundary conditions data
  allocate(block%bc%BCi(0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%bc%BCi=init_bc()
  allocate(block%bc%BCj(1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%bc%BCj=init_bc()
  allocate(block%bc%BCk(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; block%bc%BCk=init_bc()
  ! Fluid dynamic data
  allocate(block%fluid%Dt(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; block%fluid%Dt = 0._R_P
  allocate(block%fluid%P (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; call init_prim(Ns=Ns,prim=block%fluid%P)
  allocate(block%fluid%U (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; call init_cons(Ns=Ns,cons=block%fluid%U)
  allocate(block%fluid%KS(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6),1:rk_ord)) ; call init_cons(Ns=Ns,cons=block%fluid%KS)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_block

  !> Function for loading the mesh data dimensions of block from the mesh file "filename".
  !> @return \b err integer(I4P) variable.
  function load_bmesh_dims(ascii,myrank,filename,block) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical, optional, intent(IN)::    ascii    !< Flag for ascii file.
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< Name of file where mesh variables are saved.
  type(Type_Block),  intent(INOUT):: block    !< Block level data.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of mesh file.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=trim(filename),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_bmesh_dims')
  if (present(ascii)) then
    UnitFree = Get_Unit()
    open(unit = UnitFree, file = trim(filename), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err)block%mesh%gc
    read(UnitFree,*,iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    close(UnitFree)
  else
    UnitFree = Get_Unit()
    open(unit = UnitFree, file = trim(filename), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)block%mesh%gc
    read(UnitFree,iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bmesh_dims

  !> Function for loading block mesh file.
  !> @return \b err integer(I4P) variable.
  function load_bmesh(ascii,myrank,filename,block) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical, optional, intent(IN)::    ascii    !< Flag for ascii file.
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< Name of file where mesh variables are saved.
  type(Type_Block),  intent(INOUT):: block    !< Block level data.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of mesh file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_bmesh')
  if (present(ascii)) then
    UnitFree = Get_Unit()
    open(unit = UnitFree, file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err)block%mesh%gc
    read(UnitFree,*,iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    err = read(UnitFree,'*',block%mesh%node)
    err = read(UnitFree,'*',block%mesh%cell)
    close(UnitFree)
  else
    UnitFree = Get_Unit()
    open(unit = UnitFree, file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)block%mesh%gc
    read(UnitFree,iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    err = read(UnitFree,block%mesh%node)
    err = read(UnitFree,block%mesh%cell)
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bmesh

  !> Function for saving block mesh file.
  !> @return \b err integer(I4P) variable.
  function save_bmesh(ascii,filename,block) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical, optional, intent(IN):: ascii    !< Flag for ascii file.
  character(*),      intent(IN):: filename !< Name of file where mesh variables are saved.
  type(Type_Block),  intent(IN):: block    !< Block level data.
  integer(I_P)::                  err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                  UnitFree !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(ascii)) then
    UnitFree = Get_Unit()
    open(unit = UnitFree, file = adjustl(trim(filename)), action = 'WRITE', form = 'FORMATTED')
    write(UnitFree,'(6('//FI1P//',1X))',iostat=err)block%mesh%gc
    write(UnitFree,'(3('//FI_P//',1X))',iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    err = write(UnitFree,'*',block%mesh%node)
    err = write(UnitFree,'*',block%mesh%cell)
    close(UnitFree)
  else
    UnitFree = Get_Unit()
    open(unit = UnitFree, file = adjustl(trim(filename)), action = 'WRITE', form = 'UNFORMATTED')
    write(UnitFree,iostat=err)block%mesh%gc
    write(UnitFree,iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    err = write(UnitFree,block%mesh%node)
    err = write(UnitFree,block%mesh%cell)
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bmesh

  !> Function for printing to standard output info of block mesh data.
  !> @return \b err integer(I4P) variable.
  function print_info_bmesh(myrank,blk,grl,global,block) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::  myrank       !< Rank process identification for MPI communications.
  integer(I_P),      intent(IN)::  blk          !< Actual block number.
  integer(I_P),      intent(IN)::  grl          !< Actual grid level number.
  type(Type_Global), intent(IN)::  global       !< Global level data.
  type(Type_Block),  intent(IN)::  block        !< Block level data.
  integer(I_P)::                   err          !< Error trapping flag: 0 no errors, >0 error occurs.
  character(DR_P)::                vmax,vmin    !< String for printing max and min of variables.
  character(DR_P)::                xmm,ymm,zmm  !< String for printing max and min of variables.
  type(Type_Vector), allocatable:: NFi(:,:,:)   !< |
  type(Type_Vector), allocatable:: NFj(:,:,:)   !< |
  type(Type_Vector), allocatable:: NFk(:,:,:)   !< |
  real(R_P),         allocatable:: Si (:,:,:)   !< | Dummy variables for printing only internal cells info.
  real(R_P),         allocatable:: Sj (:,:,:)   !< |
  real(R_P),         allocatable:: Sk (:,:,:)   !< |
  real(R_P),         allocatable:: V  (:,:,:)   !< |
  integer(I_P)::                   Ni,Nj,Nk     !< Temporary variables for storing block dimensions.
  character(DI_P)::                rks          !< String containing myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,myrank))
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)trim(rks)//' Mesh infos of myrank '//trim(rks)
  write(stdout,'(A)',iostat=err)trim(rks)//'  Global number of blocks '//trim(str(.true.,global%mesh%Nb_tot))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Block '//trim(str(.true.,blk))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Grid level '//trim(str(.true.,grl))
  Ni = block%mesh%Ni
  Nj = block%mesh%Nj
  Nk = block%mesh%Nk
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: left  i '//trim(str(.true.,block%mesh%gc(1)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: right i '//trim(str(.true.,block%mesh%gc(2)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: left  j '//trim(str(.true.,block%mesh%gc(3)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: right j '//trim(str(.true.,block%mesh%gc(4)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: left  k '//trim(str(.true.,block%mesh%gc(5)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: right k '//trim(str(.true.,block%mesh%gc(6)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Number of cells Ni '//trim(str(.true.,Ni))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Number of cells Nj '//trim(str(.true.,Nj))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Number of cells Nk '//trim(str(.true.,Nk))
  if (allocated(NFi)) deallocate(NFi) ; allocate(NFi(0:Ni,1:Nj,1:Nk)) ; NFi(0:Ni,1:Nj,1:Nk) = block%mesh%NFi(0:Ni,1:Nj,1:Nk)
  if (allocated(NFj)) deallocate(NFj) ; allocate(NFj(1:Ni,0:Nj,1:Nk)) ; NFj(1:Ni,0:Nj,1:Nk) = block%mesh%NFj(1:Ni,0:Nj,1:Nk)
  if (allocated(NFk)) deallocate(NFk) ; allocate(NFk(1:Ni,1:Nj,0:Nk)) ; NFk(1:Ni,1:Nj,0:Nk) = block%mesh%NFk(1:Ni,1:Nj,0:Nk)
  if (allocated(Si )) deallocate(Si ) ; allocate(Si (0:Ni,1:Nj,1:Nk)) ; Si (0:Ni,1:Nj,1:Nk) = block%mesh%Si (0:Ni,1:Nj,1:Nk)
  if (allocated(Sj )) deallocate(Sj ) ; allocate(Sj (1:Ni,0:Nj,1:Nk)) ; Sj (1:Ni,0:Nj,1:Nk) = block%mesh%Sj (1:Ni,0:Nj,1:Nk)
  if (allocated(Sk )) deallocate(Sk ) ; allocate(Sk (1:Ni,1:Nj,0:Nk)) ; Sk (1:Ni,1:Nj,0:Nk) = block%mesh%Sk (1:Ni,1:Nj,0:Nk)
  if (allocated(V  )) deallocate(V  ) ; allocate(V  (1:Ni,1:Nj,1:Nk)) ; V  (1:Ni,1:Nj,1:Nk) = block%mesh%V  (1:Ni,1:Nj,1:Nk)
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

  !> Function for loading inflow 1 boundary conditions.
  !> @return \b err integer(I4P) variable.
  function load_gbc_in1(myrank,filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< Name of file where boundary conditions variables are saved.
  type(Type_Global), intent(INOUT):: global   !< Global level data.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of boundary conditions file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_gbc_in1')
  UnitFree = Get_Unit()
  open(unit = UnitFree, file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  err = read(UnitFree,'*',global%bc%in1)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gbc_in1

  !> Function for loading block boundary conditions file.
  !> @return \b err integer(I4P) variable.
  function load_bbc(myrank,filename,block) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),     intent(IN)::    myrank   !< Actual rank process.
  character(*),     intent(IN)::    filename !< Name of file where mesh variables are saved.
  type(Type_Block), intent(INOUT):: block    !< Block level data.
  integer(I_P)::                    err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                    UnitFree !< Free logic unit.
  logical::                         is_file  !< Flag for inquiring the presence of boundary conditions file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_bbc')
  UnitFree = Get_Unit()
  open(unit=UnitFree, file=adjustl(trim(filename)), status='OLD', action='READ', form='UNFORMATTED')
  err = read(UnitFree,block%bc%BCi)
  err = read(UnitFree,block%bc%BCj)
  err = read(UnitFree,block%bc%BCk)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bbc

  !> Function for saving block boundary conditions file.
  !> @return \b err integer(I4P) variable.
  function save_bbc(filename,block) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),     intent(IN):: filename !< Name of file where mesh variables are saved.
  type(Type_Block), intent(IN):: block    !< Block level data.
  integer(I_P)::                 err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                 UnitFree !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  UnitFree = Get_Unit()
  open(unit=UnitFree, file=adjustl(trim(filename)), action='WRITE', form='UNFORMATTED')
  err = write(UnitFree,block%bc%BCi)
  err = write(UnitFree,block%bc%BCj)
  err = write(UnitFree,block%bc%BCk)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bbc

  !> Function for loading solver option variables.
  !> @return \b err integer(I4P) variable.
  function load_gfluid_soption(myrank,filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    myrank   !< Rank process identification for MPI communications.
  character(*),      intent(IN)::    filename !< Name of file where option variables are saved.
  type(Type_Global), intent(INOUT):: global   !< Global level data.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of option file.
  character(8)::                     timing   !< timing = 'unsteady' simulation otherwise steady one.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 1
  inquire(file=adjustl(trim(filename)),exist=is_file)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_gfluid_soption')
  UnitFree = Get_Unit()
  open(unit = UnitFree, file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*,iostat=err) ! record skipped because unnecessary
  read(UnitFree,*,iostat=err) ! record skipped because unnecessary
  read(UnitFree,*,iostat=err)timing ; timing = Upper_Case(timing) ; global%fluid%unsteady = (adjustl(trim(timing))=='UNSTEADY')
  read(UnitFree,*,iostat=err)global%fluid%Nmax
  read(UnitFree,*,iostat=err)global%fluid%Tmax
  read(UnitFree,*,iostat=err)global%fluid%sp_ord
  read(UnitFree,*,iostat=err)global%fluid%rk_ord
  read(UnitFree,*,iostat=err)global%fluid%CFL
  read(UnitFree,*,iostat=err)global%fluid%residual_toll
  close(UnitFree)
  if (global%fluid%unsteady) then
    if (global%fluid%Nmax>0_I8P) then
      global%fluid%Tmax=-1._R_P ! the value of Tmax is ignored because Nmax>0
    endif
  else
    global%fluid%Tmax=-1._R_P ! the value of Tmax is ignored because steady simulation
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gfluid_soption

  !> Function for loading the number of species from fluid dynamic file "filename".
  !> @return \b err integer(I4P) variable.
  function load_gfluid_Ns(binary,myrank,filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical,           intent(IN), optional:: binary   !< Flag for binary of ascii input file.
  integer(I_P),      intent(IN)::           myrank   !< Actual rank process.
  character(*),      intent(IN)::           filename !< Name of file where mesh variables are saved.
  type(Type_Global), intent(INOUT)::        global   !< Global level data.
  integer(I_P)::                            err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                            UnitFree !< Free logic unit.
  logical::                                 is_file  !< Flag for inquiring the presence of fluid dynamic file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=trim(filename),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_gfluid_Ns')
  UnitFree = Get_Unit()
  if (present(binary)) then
    open(unit = UnitFree, file = trim(filename), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)global%fluid%Ns
  else
    open(unit = UnitFree, file = trim(filename), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err) ! record skipped
    read(UnitFree,*,iostat=err) ! record skipped
    read(UnitFree,*,iostat=err)global%fluid%Ns
  endif
  close(UnitFree)
  global%fluid%Np = global%fluid%Ns + 6
  global%fluid%Nc = global%fluid%Ns + 4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gfluid_Ns

  !> Function for loading initial species from ascii file.
  !> @return \b err integer(I4P) variable.
  function load_gfluid_0species(myrank,filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< Name of file where initial species are saved.
  type(Type_Global), intent(INOUT):: global   !< Global level data.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of initial species file.
  integer(I_P)::                     s        !< Species counter.
  real(R_P)::                        g        !< Specific heats ratio (cp/cv).
  real(R_P)::                        R        !< Specific heats difference (cp-cv).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_gfluid_0species')
  UnitFree = Get_Unit()
  open(unit = UnitFree, file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*,iostat=err) ! record skipped
  read(UnitFree,*,iostat=err) ! record skipped
  read(UnitFree,*,iostat=err)global%fluid%Ns ; call alloc_global_fluid(global=global)
  do s=1,global%fluid%Ns
    read(UnitFree,*,iostat=err) ! record skipped
    read(UnitFree,*,iostat=err)g
    read(UnitFree,*,iostat=err)R
    read(UnitFree,*,iostat=err)global%fluid%cp0(s)
    read(UnitFree,*,iostat=err)global%fluid%cv0(s)
    if ((g>0.0_R_P).AND.(R>0.0_R_P)) then
      global%fluid%cv0(s)=R/(g-1._R_P)
      global%fluid%cp0(s)=g*global%fluid%cv0(s)
    elseif ((g>0.0_R_P).AND.(global%fluid%cp0(s)>0.0_R_P)) then
      global%fluid%cv0(s)=global%fluid%cp0(s)/g
    elseif ((g>0.0_R_P).AND.(global%fluid%cv0(s)>0.0_R_P)) then
      global%fluid%cp0(s)=global%fluid%cv0(s)*g
    elseif ((R>0.0_R_P).AND.(global%fluid%cp0(s)>0.0_R_P)) then
      global%fluid%cv0(s)=global%fluid%cp0(s)-R
    elseif ((R>0.0_R_P).AND.(global%fluid%cv0(s)>0.0_R_P)) then
      global%fluid%cp0(s)=global%fluid%cv0(s)+R
    endif
  enddo
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gfluid_0species

  !> Function for loading fluid dynamic file.
  !> @return \b err integer(I4P) variable.
  function load_bfluid(myrank,filename,global,block) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< Name of file where mesh variables are saved.
  type(Type_Global), intent(INOUT):: global   !< Global level data.
  type(Type_Block),  intent(INOUT):: block    !< Block level data.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of fluid dynamic file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_bfluid')
  UnitFree = Get_Unit()
  open(unit = UnitFree, file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
  read(UnitFree,iostat=err)global%fluid%Ns
  read(UnitFree,iostat=err)global%fluid%cp0,global%fluid%cv0
  read(UnitFree,iostat=err)global%fluid%n
  read(UnitFree,iostat=err)global%fluid%t
  read(UnitFree,iostat=err)block%fluid%Dt
  err = read(UnitFree,block%fluid%P)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bfluid

  !> Function for saving fluid dynamic file.
  !> @return \b err integer(I4P) variable.
  function save_bfluid(filename,global,block) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),      intent(IN):: filename !< Name of file where mesh variables are saved.
  type(Type_Global), intent(IN):: global   !< Global level data.
  type(Type_Block),  intent(IN):: block    !< Block level data.
  integer(I_P)::                  err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                  UnitFree !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  UnitFree = Get_Unit()
  open(unit = UnitFree, file = adjustl(trim(filename)), action = 'WRITE', form = 'UNFORMATTED')
  write(UnitFree,iostat=err)global%fluid%Ns
  write(UnitFree,iostat=err)global%fluid%cp0,global%fluid%cv0
  write(UnitFree,iostat=err)global%fluid%n
  write(UnitFree,iostat=err)global%fluid%t
  write(UnitFree,iostat=err)block%fluid%Dt
  err = write(UnitFree,block%fluid%P)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bfluid

  !> Function for printing to standard output info of block fluid data.
  !> @return \b err integer(I4P) variable.
  function print_info_bfluid(myrank,blk,grl,global,block) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN):: myrank         !< Rank process identification for MPI communications.
  integer(I_P),      intent(IN):: blk            !< Actual block number.
  integer(I_P),      intent(IN):: grl            !< Actual grid level number.
  type(Type_Global), intent(IN):: global         !< Global level data.
  type(Type_Block),  intent(IN):: block          !< Block level data.
  integer(I_P)::                  err            !< Error trapping flag: 0 no errors, >0 error occurs.
  character(DR_P)::               vmax,vmin      !< String for printing max and min of variables.
  real(R_P), allocatable::        dummy(:,:,:,:) !<  Dummy variables for printing only internal cells info.
  integer(I_P)::                  Ni,Nj,Nk,Np,Nc !< Temporary variables for storing block dimensions.
  character(DI_P)::               rks            !< String containing myrank.
  integer(I_P)::                  i,j,k,s        !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,myrank))
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)trim(rks)//' Fluid infos of myrank '//trim(rks)
  write(stdout,'(A)',iostat=err)trim(rks)//'  Global number of blocks '//trim(str(.true.,global%mesh%Nb_tot))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Number of initial species '//trim(str(.true.,global%fluid%Ns))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Block '//trim(str(.true.,blk))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Grid level '//trim(str(.true.,grl))
  Ni = block%mesh%Ni
  Nj = block%mesh%Nj
  Nk = block%mesh%Nk
  Np = global%fluid%Np
  Nc = global%fluid%Nc
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Nc,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(:,i,j,k) = cons2array(block%fluid%U(i,j,k))
      enddo
    enddo
  enddo
  write(stdout,'(A)',iostat=err)trim(rks)//'  Conservative variables'
  do s=1,global%fluid%Nc
    vmax=trim(str(n=maxval(dummy(s,:,:,:)))) ; vmin=trim(str(n=minval(dummy(s,:,:,:))))
    write(stdout,'(A)',iostat=err)trim(rks)//'  U('//trim(str(.true.,s))//'): max '//trim(vmax)//' min '//trim(vmin)
  enddo
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Np,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(:,i,j,k) = prim2array(block%fluid%P(i,j,k))
      enddo
    enddo
  enddo
  write(stdout,'(A)',iostat=err)trim(rks)//'  Primitive variables'
  do s=1,global%fluid%Np
    vmax=trim(str(n=maxval(dummy(s,:,:,:)))) ; vmin=trim(str(n=minval(dummy(s,:,:,:))))
    write(stdout,'(A)',iostat=err)trim(rks)//'  P('//trim(str(.true.,s))//'): max '//trim(vmax)//' min '//trim(vmin)
  enddo
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(1,i,j,k) = block%fluid%Dt(i,j,k)
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
endmodule Data_Type_Globals
