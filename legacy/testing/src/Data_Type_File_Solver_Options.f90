!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_Solver_OptionsDerivedType Data_Type_File_Solver_Options
!> Module definition of Type_File_Solver_Options
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_File_Solver_OptionsInterface Data_Type_File_Solver_Options
!> Module definition of Type_File_Solver_Options
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_File_Solver_OptionsPublicProcedure Data_Type_File_Solver_Options
!> Module definition of Type_File_Solver_Options
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_Solver_OptionsPrivateProcedure Data_Type_File_Solver_Options
!> Module definition of Type_File_Solver_Options
!> @}

!> @brief Module Data_Type_File_Solver_Options contains the definition of Type_File_Solver_Options, that defines the main fluid
!> dynamic solver options.
module Data_Type_File_Solver_Options
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                    ! Integers and reals precision definition.
USE Data_Type_Adimensional, only: Type_Adimensional ! Definition of Type_Adimensional.
USE Data_Type_File_Base,    only: Type_File_Base    ! Definition of Type_File_Base.
USE Data_Type_Space_Step,   only: Type_Space_Step   ! Definition of Type_Space_Step.
USE Data_Type_Time_Step,    only: Type_Time_Step    ! Definition of Type_Time_Step.
USE Data_Type_XML_Tag,      only: Type_XML_Tag      ! Definition of Type_XML_Tag.
USE Lib_IO_Misc,            only: iostat_eor        ! Procedures for IO and strings operations.
USE Lib_Strings,            only: Upper_Case        ! Library for strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Solver_Options.
!> @note This options file is a very simple XML file that drives the running of OFF program.
!> @note This options file is formatted by means of XML syntax. A validated file should be formatted as following:
!> @code
!>   <?xml version="1.0"?>
!>   <Solver_Options>
!>     <timing> unsteady(or steady) </timing>
!>     <max_time_iterations> ... </max_time_iterations>
!>     <final_time> ... </final_time>
!>     <time_integration_accuracy> 1(or 2,3,4,5) </time_integration_accuracy>
!>     <space_integration_accuracy> 1(or 3,5,7) </space_integration_accuracy>
!>     <CFL> ... </CFL>
!>     <residuals_tolerance> ... </residuals_tolerance>
!>     <riemann_solver> ... </riemann_solver>
!>     <reconstruction_type> ... </reconstruction_type>
!>   </Solver_Options>
!>   <Non_Dimensional_Numbers>
!>     <Re> ... </Re>
!>     <Fr> ... </Fr>
!>     <Pr> ... </Pr>
!>   </Non_Dimensional_Numbers>
!>   <Reference_Values>
!>     <L0> ... </L0>
!>     <r0> ... </r0>
!>     <v0> ... </v0>
!>     <c0> ... </c0>
!>   </Reference_Values>
!> @endcode
!> There are 3 main tags, namely 'Solver_Options', that contains the main fluid dynamic options, 'Non_Dimensional_Numbers',
!> that contains the value of Reynolds, Froude and Prandtl Numbers, and 'Reference_Values', that contains the reference length
!> (L0), density (r0), velocity (v0) and specific heat (c0). Note that the final time value is ignored if the max number of
!> iterations provided is greater than zero.
!> The solver options file contains the following data:
!>   - file_d%time_step%unsteady;
!>   - file_d%time_step%Nmax;
!>   - file_d%time_step%Tmax;
!>   - file_d%time_step%rk_ord;
!>   - file_d%space_step%sp_ord;
!>   - file_d%space_step%CFL;
!>   - file_d%adim%Re;
!>   - file_d%adim%Fr;
!>   - file_d%adim%Pr;
!>   - file_d%adim%L0;
!>   - file_d%adim%r0;
!>   - file_d%adim%v0;
!>   - file_d%adim%c0.
!> It is worth nothing that the syntax of tag names and attributes is case sensitive, whereas the syntax of their values is case
!> insensitive.
!> @ingroup Data_Type_File_Solver_OptionsDerivedType
type, public, extends(Type_File_Base):: Type_File_Solver_Options
  type(Type_XML_Tag)::            Solver_Options          !< Solver options XML tag.
  type(Type_XML_Tag)::            Non_Dimensional_Numbers !< Non dimensional numbers Main XML tag.
  type(Type_XML_Tag)::            Reference_Values        !< Reference values XML tag.
  type(Type_Time_Step)::          time_step               !< Time step data.
  type(Type_Space_Step)::         space_step              !< Space step data.
  type(Type_Adimensional)::       adim                    !< Non-dimensionalization data.
  character(len=:), allocatable:: RSU                     !< Riemann solver used.
  character(len=:), allocatable:: recon_tp                !< Reconstruction type.
  contains
    procedure:: free  => free_fs_options  ! Procedure for freeing dynamic memory.
    procedure:: load  => load_fs_options  ! Procedure for loading fluid dynamic solver options.
    procedure:: save  => save_fs_options  ! Procedure for saving  fluid dynamic solver options.
    procedure:: print => print_fs_options ! Procedure for printing solver options with a pretty format.
    final::     finalize_fs_options       ! Procedure for freeing dynamic memory when finalizing.
endtype Type_File_Solver_Options
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_Solver_OptionsPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_fs_options(file_d,also_base)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Solver_Options), intent(INOUT):: file_d    !< File data.
  logical, optional,               intent(IN)::    also_base !< Flag for freeing also file_base dynamic memory.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(also_base)) then
    if (also_base) call file_d%free_base
  endif
  if (allocated(file_d%RSU     )) deallocate(file_d%RSU     )
  if (allocated(file_d%recon_tp)) deallocate(file_d%recon_tp)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_fs_options

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize_fs_options(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_File_Solver_Options), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free(also_base=.true.)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize_fs_options

  !> @brief Procedure for loading solver option variables.
  subroutine load_fs_options(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Solver_Options), intent(INOUT):: file_d  !< File data.
  type(Type_XML_Tag)::                             tag     !< Generic XML tag.
  character(len=:), allocatable::                  stream  !< String containing the file data as a single stream.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%get_stream(stream=stream)

  call file_d%Solver_Options%set(          tag_name = 'Solver_Options'          , search_into_string = stream)
  call file_d%Non_Dimensional_Numbers%set( tag_name = 'Non_Dimensional_Numbers' , search_into_string = stream)
  call file_d%Reference_Values%set(        tag_name = 'Reference_Values'        , search_into_string = stream)

  call tag%set(tag_name='timing',search_into_other_tag=file_d%Solver_Options)
  file_d%time_step%unsteady = (Upper_Case(tag%tag_val%vs)=='UNSTEADY')
  call tag%set(tag_name='max_time_iterations',search_into_other_tag=file_d%Solver_Options)
  file_d%time_step%Nmax = cton(str=tag%tag_val%vs,knd=1_I8P)
  call tag%set(tag_name='final_time',search_into_other_tag=file_d%Solver_Options)
  file_d%time_step%Tmax = cton(str=tag%tag_val%vs,knd=1._R8P)
  call tag%set(tag_name='time_integration_accuracy',search_into_other_tag=file_d%Solver_Options)
  file_d%time_step%rk_ord = cton(str=tag%tag_val%vs,knd=1_I1P)
  call tag%set(tag_name='space_integration_accuracy',search_into_other_tag=file_d%Solver_Options)
  file_d%space_step%sp_ord = cton(str=tag%tag_val%vs,knd=1_I1P)
  call tag%set(tag_name='CFL',search_into_other_tag=file_d%Solver_Options)
  file_d%time_step%CFL = cton(str=tag%tag_val%vs,knd=1._R8P)
  call tag%set(tag_name='residuals_tolerance',search_into_other_tag=file_d%Solver_Options)
  file_d%space_step%residual_toll = cton(str=tag%tag_val%vs,knd=1._R8P)
  call tag%set(tag_name='riemann_solver',search_into_other_tag=file_d%Solver_Options)
  file_d%RSU = tag%tag_val%vs
  call tag%set(tag_name='reconstruction_type',search_into_other_tag=file_d%Solver_Options)
  file_d%recon_tp = tag%tag_val%vs

  call tag%set(tag_name='Re',search_into_other_tag=file_d%Non_Dimensional_Numbers)
  file_d%adim%Re = cton(str=tag%tag_val%vs,knd=1._R8P)
  call tag%set(tag_name='Fr',search_into_other_tag=file_d%Non_Dimensional_Numbers)
  file_d%adim%Fr = cton(str=tag%tag_val%vs,knd=1._R8P)
  call tag%set(tag_name='Pr',search_into_other_tag=file_d%Non_Dimensional_Numbers)
  file_d%adim%Pr = cton(str=tag%tag_val%vs,knd=1._R8P)

  call tag%set(tag_name='L0',search_into_other_tag=file_d%Reference_Values)
  file_d%adim%L0 = cton(str=tag%tag_val%vs,knd=1._R8P)
  call tag%set(tag_name='r0',search_into_other_tag=file_d%Reference_Values)
  file_d%adim%r0 = cton(str=tag%tag_val%vs,knd=1._R8P)
  call tag%set(tag_name='v0',search_into_other_tag=file_d%Reference_Values)
  file_d%adim%v0 = cton(str=tag%tag_val%vs,knd=1._R8P)
  call tag%set(tag_name='c0',search_into_other_tag=file_d%Reference_Values)
  file_d%adim%c0 = cton(str=tag%tag_val%vs,knd=1._R8P)

  call file_d%fallback

  call file_d%space_step%set_gc    ! setting the number of ghost cells used
  call file_d%time_step%set_limits ! setting limits for simulation stop condition
  call file_d%adim%compute_values0 ! computing the reference values for non-dimensional quantities
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_fs_options

  !> @brief Procedure for saving solver options file.
  subroutine save_fs_options(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Solver_Options), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg,&
            time_step=>file_d%time_step,space_step=>file_d%space_step,adim=>file_d%adim,RSU=>file_d%RSU,recon_tp=>file_d%recon_tp)
    call file_d%open(ascii=.true.,action='WRITE') ; if (iostat/=0) return
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<?xml version="1.0"?>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<Solver_Options>'
    if (time_step%unsteady) then
      write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <timing>unsteady</timing>'
    else
      write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <timing>steady</timing>'
    endif
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <max_time_iterations>'//trim(str(n=time_step%Nmax))//&
                                                         '</max_time_iterations>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <final_time>'//trim(str(n=time_step%Tmax))//'</final_time>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <time_integration_accuracy>'//trim(str(n=time_step%rk_ord)) //&
                                                         '</time_integration_accuracy>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <CFL>'//trim(str(n=time_step%CFL))//'</CFL>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <space_integration_accuracy>'//trim(str(n=space_step%sp_ord))//&
                                                         '</space_integration_accuracy>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <residuals_tolerance>'//trim(str(n=space_step%residual_toll))//&
                                                         '</residuals_tolerance>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <riemann_solver>'//RSU//'</riemann_solver>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <reconstruction_type>'//recon_tp//'</reconstruction_type>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'</Solver_Options>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<Non_Dimensional_Numbers>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <Re>'//trim(str(n=adim%Re))//'</Re>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <Fr>'//trim(str(n=adim%Fr))//'</Fr>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <Pr>'//trim(str(n=adim%Pr))//'</Pr>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'</Non_Dimensional_Numbers>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<Reference_Values>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <L0>'//trim(str(n=adim%L0))//'</L0>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <r0>'//trim(str(n=adim%r0))//'</r0>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <v0>'//trim(str(n=adim%v0))//'</v0>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <c0>'//trim(str(n=adim%c0))//'</c0>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'</Reference_Values>'
    call file_d%close
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_fs_options

  !> @brief Procedure for printing solver options with a pretty format.
  subroutine print_fs_options(file_d,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Solver_Options), intent(INOUT):: file_d  !< File data.
  character(*), optional,          intent(IN)::    pref    !< Prefixing string.
  integer(I4P), optional,          intent(OUT)::   iostat  !< IO error.
  character(*), optional,          intent(OUT)::   iomsg   !< IO error message.
  integer(I4P),                    intent(IN)::    unit    !< Logic unit.
  character(len=:), allocatable::                  prefd   !< Prefixing string.
  integer(I4P)::                                   iostatd !< IO error.
  character(500)::                                 iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Solver options'
  call file_d%time_step%print( unit=unit,iostat=iostatd,iomsg=iomsgd,pref=prefd//'  ')
  call file_d%space_step%print(unit=unit,iostat=iostatd,iomsg=iomsgd,pref=prefd//'  ')
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Riemann solver used: '//file_d%RSU
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Variable reconstruction type: '//file_d%recon_tp
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Non dimensional numbers and reference values'
  call file_d%adim%print(unit=unit,iostat=iostatd,iomsg=iomsgd,pref=prefd//'  ')
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_fs_options
  !> @}
endmodule Data_Type_File_Solver_Options
