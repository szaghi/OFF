!> @ingroup Library
!> @{
!> @defgroup Lib_FluidynamicLibrary Lib_Fluidynamic
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_FluidynamicPublicProcedure Lib_Fluidynamic
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_FluidynamicPrivateProcedure Lib_Fluidynamic
!> @}

!> This module contains the definition of fluid dynamic procedures.
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_FluidynamicLibrary
module Lib_Fluidynamic
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                   ! Integers and reals precision definition.
USE Data_Type_BC                                   ! Definition of Type_BC.
USE Data_Type_Cell                                 ! Definition of Type_Cell.
USE Data_Type_Conservative                         ! Definition of Type_Conservative.
USE Data_Type_Face                                 ! Definition of Type_Face.
USE Data_Type_Files,  only: Type_Files                  ! Definition of Type_Files.
USE Data_Type_Global                               ! Definition of Type_Global.
USE Data_Type_Primitive                            ! Definition of Type_Primitive.
USE Data_Type_SBlock                               ! Definition of Type_SBlock.
USE Data_Type_Time                                 ! Definition of Type_Time.
USE Data_Type_Vector                               ! Definition of Type_Vector.
USE Lib_Fluxes_Convective, only: fluxes_convective ! Subroutine for computing convective fluxes.
!USE Lib_Fluxes_Diffusive, only: fluxes_diffusive   ! Subroutine for computing diffusive fluxes.
USE Lib_IO_Misc                                    ! Procedures for IO and strings operations.
USE Lib_Math, only: interpolate1                   ! Function for computing linear interpolation.
#ifdef PROFILING
USE Lib_Profiling                                  ! Procedures for profiling the code.
#endif
USE Lib_Runge_Kutta                                ! Runge-Kutta time integration library.
USE Lib_Thermodynamic_Laws_Ideal, only: a          ! Function for computing speed of sound.
#ifdef _MPI
USE Lib_Parallel, only: procmap, &                 ! Proc/blocks map.
                        prim_sendrecv              ! Library for send/receive data for parallel (MPI) operations.
USE MPI                                            ! MPI runtime library.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
integer(I1P):: flip = 0_I1P !< Flip-Flop flag for restart solution file.
!-----------------------------------------------------------------------------------------------------------------------------------
contains

  !> Subroutine for computing Runge-Kutta one time step integration.
  !> @ingroup Lib_FluidynamicPrivateProcedure
  subroutine rk_time_integration(block,RU)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_SBlock), intent(INOUT):: block               !< Block-level data.
  real(R_P),         intent(OUT)::   RU(1:)              !< NormL2 of residuals of conservative variables.
  real(R_P)                          R(1:size(RU,dim=1)) !< Residuals of conservative variables.
  integer(I4P)::                     i,j,k               !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  RU = 0._R_P
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k,R)       &
  !$OMP SHARED(block)          &
  !$OMP REDUCTION(+: RU)
  !$OMP DO
  do k=1,block%dims%Nk
    do j=1,block%dims%Nj
      do i=1,block%dims%Ni
        call rk_time_integration_backend(Ns=block%dims%Ns,rk_ord=block%dims%Nrk,&
                                         Dt=block%C(i,j,k)%Dt,KS=block%C(i,j,k)%KS,U=block%C(i,j,k)%U,R=R)
        RU = RU + (R*R)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  RU = sqrt(RU)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    pure subroutine rk_time_integration_backend(Ns,rk_ord,Dt,KS,U,R)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),            intent(IN)::    Ns                            !< Number of species.
    integer(I1P),            intent(IN)::    rk_ord                        !< Number of Runge-Kutta stages.
    real(R_P),               intent(IN)::    Dt                            !< Current time step.
    type(Type_Conservative), intent(IN)::    KS(1:)                        !< Runge-Kutta stages.
    type(Type_Conservative), intent(INOUT):: U                             !< Conservative variables.
    real(R_P),               intent(OUT)::   R(1:)                         !< Residuals of conservative variables.
    real(R_P)::                              KSa(1:size(R,dim=1),1:rk_ord) !< Dummy Runge-Kutta stages.
    integer(I1P)::                           s1                            !< Counter.
    integer(I4P)::                           s                             !< Counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    do s1=1_I1P,rk_ord
      KSa(:,s1) = KS(s1)%cons2array()
    enddo
    do s=1_I4P,Ns
      call rk_time_integ(Dt=Dt,KS=KSa(s,:),R=R(s),Unp1=U%rs(s))
    enddo
    call rk_time_integ(Dt=Dt,KS=KSa(Ns+1_I4P,:),R=R(Ns+1_I4P),Unp1=U%rv%x)
    call rk_time_integ(Dt=Dt,KS=KSa(Ns+2_I4P,:),R=R(Ns+2_I4P),Unp1=U%rv%y)
    call rk_time_integ(Dt=Dt,KS=KSa(Ns+3_I4P,:),R=R(Ns+3_I4P),Unp1=U%rv%z)
    call rk_time_integ(Dt=Dt,KS=KSa(Ns+4_I4P,:),R=R(Ns+4_I4P),Unp1=U%re  )
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine rk_time_integration_backend
  endsubroutine rk_time_integration
endmodule Lib_Fluidynamic
