!> @ingroup Library
!> @{
!> @defgroup Lib_Fluxes_DiffusiveLibrary Lib_Fluxes_Diffusive
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Fluxes_DiffusivePublicProcedure Lib_Fluxes_Diffusive
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_Fluxes_DiffusivePrivateProcedure Lib_Fluxes_Diffusive
!> @}

!> This module contains the definition of procedures for computing diffusive fluxes (parabolic operator).
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_Fluxes_DiffusiveLibrary
module Lib_Fluxes_Diffusive
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                             ! Integers and reals precision definition.
USE Data_Type_Conservative                   ! Definition of Type_Conservative.
USE Data_Type_Primitive                      ! Definition of Type_Primitive.
USE Data_Type_SBlock                         ! Definition of Type_SBlock.
USE Data_Type_Tensor, sq_norm_ten => sq_norm ! Definition of Type_Tensor.
USE Data_Type_Vector, sq_norm_vec => sq_norm ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: fluxes_diffusive
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> Subroutine for computing interfaces diffusive fluxes (parabolic operator).
  subroutine fluxes_diffusive(block,i,j,k,dir,F)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_SBlock),       intent(IN)::    block  !< Block-level data.
  integer(I_P),            intent(IN)::    i      !< Interface i indexe.
  integer(I_P),            intent(IN)::    j      !< Interface j indexe.
  integer(I_P),            intent(IN)::    k      !< Interface k indexe.
  character(1),            intent(IN)::    dir    !< Direction of fluxes ('i','j','k').
  type(Type_Conservative), intent(INOUT):: F      !< Diffusive fluxes at the interface.
  type(Type_Vector):: NFS(1:6) !< Normals, with interfaces area as module, of the 6 hexahedron's faces of finite
                               !< volume centered at interface i,j,k.
  type(Type_Vector):: vel(1:6) !< Velocity vector at interfaces of finite volume centered at interface i,j,k.
  type(Type_Tensor):: tau      !< Shear stress tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing metrics of finite volume centered at interface i,j,k
  call fv_interface_metrics(block=block)
  ! computing velocity vector at the interfaces of the finite volume centered at interface i,j,k
  call fv_interface_velocity(block=block)
  ! computing the shear stress tensor at interface i,j,k
  call shear_stress()
  ! computing the diffusive fluxes
  select case(dir)
  case('i')
    F%rs = 0._R_P
    F%rv = -(tau.dot.block%Fi(i,j,k)%N)
    F%re = -((tau.dot.(0.5_R_P*(block%C(i,j,k)%P%v+block%C(i+1,j,k)%P%v))).dot.block%Fi(i,j,k)%N)
  case('j')
    F%rs = 0._R_P
    F%rv = -(tau.dot.block%Fj(i,j,k)%N)
    F%re = -((tau.dot.(0.5_R_P*(block%C(i,j,k)%P%v+block%C(i,j+1,k)%P%v))).dot.block%Fj(i,j,k)%N)
  case('k')
    F%rs = 0._R_P
    F%rv = -(tau.dot.block%Fk(i,j,k)%N)
    F%re = -((tau.dot.(0.5_R_P*(block%C(i,j,k)%P%v+block%C(i,j,k+1)%P%v))).dot.block%Fk(i,j,k)%N)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> Subroutine for computing the metrics (normals and interfaces area) of finite volume centered at interface i,j,k.
    !> Such a 'shifted' finite volume is used in computing partial derivative by means of finite volume approach (using Green's
    !> theorem).
    !> @note The computed metrics is stored as normals to the 6 hexahedron faces of the shifted finite volume. Each normal has
    !> module equal to its face area.
    !> @return \b NFS(1:6) Type_Vector array.
    subroutine fv_interface_metrics(block)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    type(Type_SBlock), intent(IN):: block !< Block-level data.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    select case(dir)
    case('i')
      NFS(1) = 0.5_R_P*(block%Fi(i,j  ,k  )%S*block%Fi(i,j  ,k  )%N + block%Fi(i-1,j  ,k  )%S*block%Fi(i-1,j  ,k  )%N)
      NFS(2) = 0.5_R_P*(block%Fi(i,j  ,k  )%S*block%Fi(i,j  ,k  )%N + block%Fi(i+1,j  ,k  )%S*block%Fi(i+1,j  ,k  )%N)
      NFS(3) = 0.5_R_P*(block%Fj(i,j-1,k  )%S*block%Fj(i,j-1,k  )%N + block%Fj(i+1,j-1,k  )%S*block%Fj(i+1,j-1,k  )%N)
      NFS(4) = 0.5_R_P*(block%Fj(i,j  ,k  )%S*block%Fj(i,j  ,k  )%N + block%Fj(i+1,j  ,k  )%S*block%Fj(i+1,j  ,k  )%N)
      NFS(5) = 0.5_R_P*(block%Fk(i,j  ,k-1)%S*block%Fk(i,j  ,k-1)%N + block%Fk(i+1,j  ,k-1)%S*block%Fk(i+1,j  ,k-1)%N)
      NFS(6) = 0.5_R_P*(block%Fk(i,j  ,k  )%S*block%Fk(i,j  ,k  )%N + block%Fk(i+1,j  ,k  )%S*block%Fk(i+1,j  ,k  )%N)
    case('j')
      NFS(1) = 0.5_R_P*(block%Fj(i  ,j,k  )%S*block%Fj(i  ,j,k  )%N + block%Fj(i  ,j-1,k  )%S*block%Fj(i  ,j-1,k  )%N)
      NFS(2) = 0.5_R_P*(block%Fj(i  ,j,k  )%S*block%Fj(i  ,j,k  )%N + block%Fj(i  ,j+1,k  )%S*block%Fj(i  ,j+1,k  )%N)
      NFS(3) = 0.5_R_P*(block%Fk(i  ,j,k-1)%S*block%Fk(i  ,j,k-1)%N + block%Fk(i  ,j+1,k-1)%S*block%Fk(i  ,j+1,k-1)%N)
      NFS(4) = 0.5_R_P*(block%Fk(i  ,j,k  )%S*block%Fk(i  ,j,k  )%N + block%Fk(i  ,j+1,k  )%S*block%Fk(i  ,j+1,k  )%N)
      NFS(5) = 0.5_R_P*(block%Fi(i-1,j,k  )%S*block%Fi(i-1,j,k  )%N + block%Fi(i-1,j+1,k  )%S*block%Fi(i-1,j+1,k  )%N)
      NFS(6) = 0.5_R_P*(block%Fi(i  ,j,k  )%S*block%Fi(i  ,j,k  )%N + block%Fi(i  ,j+1,k  )%S*block%Fi(i  ,j+1,k  )%N)
    case('k')
      NFS(1) = 0.5_R_P*(block%Fk(i  ,j  ,k)%S*block%Fk(i  ,j  ,k)%N + block%Fk(i  ,j  ,k-1)%S*block%Fk(i  ,j  ,k-1)%N)
      NFS(2) = 0.5_R_P*(block%Fk(i  ,j  ,k)%S*block%Fk(i  ,j  ,k)%N + block%Fk(i  ,j  ,k+1)%S*block%Fk(i  ,j  ,k+1)%N)
      NFS(3) = 0.5_R_P*(block%Fi(i-1,j  ,k)%S*block%Fi(i-1,j  ,k)%N + block%Fi(i-1,j  ,k+1)%S*block%Fi(i-1,j  ,k+1)%N)
      NFS(4) = 0.5_R_P*(block%Fi(i  ,j  ,k)%S*block%Fi(i  ,j  ,k)%N + block%Fi(i  ,j  ,k+1)%S*block%Fi(i  ,j  ,k+1)%N)
      NFS(5) = 0.5_R_P*(block%Fj(i  ,j-1,k)%S*block%Fj(i  ,j-1,k)%N + block%Fj(i  ,j-1,k+1)%S*block%Fj(i  ,j-1,k+1)%N)
      NFS(6) = 0.5_R_P*(block%Fj(i  ,j  ,k)%S*block%Fj(i  ,j  ,k)%N + block%Fj(i  ,j  ,k+1)%S*block%Fj(i  ,j  ,k+1)%N)
    endselect
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine fv_interface_metrics

    !> Subroutine for computing the velocity vector at interfaces of finite volume centered at interface i,j,k.
    !> @return \b vel(1:6) Type_Vector array.
    subroutine fv_interface_velocity(block)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    type(Type_SBlock), intent(IN):: block !< Block-level data.
    type(Type_Vector)::             mean  !< Mean velocity vector across the interface.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    select case(dir)
    case('i')
      mean = block%C(i,j,k)%P%v + block%C(i+1,j,k)%P%v

      vel(1) = block%C(i  ,j,k)%P%v
      vel(2) = block%C(i+1,j,k)%P%v
      vel(3) = 0.25_R_P*(mean + block%C(i,j-1,k  )%P%v + block%C(i+1,j-1,k  )%P%v)
      vel(4) = 0.25_R_P*(mean + block%C(i,j+1,k  )%P%v + block%C(i+1,j+1,k  )%P%v)
      vel(5) = 0.25_R_P*(mean + block%C(i,j  ,k-1)%P%v + block%C(i+1,j  ,k-1)%P%v)
      vel(6) = 0.25_R_P*(mean + block%C(i,j  ,k+1)%P%v + block%C(i+1,j  ,k+1)%P%v)
    case('j')
      mean = block%C(i,j,k)%P%v + block%C(i,j+1,k)%P%v

      vel(1) = block%C(i,j  ,k)%P%v
      vel(2) = block%C(i,j+1,k)%P%v
      vel(3) = 0.25_R_P*(mean + block%C(i  ,j,k-1)%P%v + block%C(i  ,j+1,k-1)%P%v)
      vel(4) = 0.25_R_P*(mean + block%C(i  ,j,k+1)%P%v + block%C(i  ,j+1,k+1)%P%v)
      vel(5) = 0.25_R_P*(mean + block%C(i-1,j,k  )%P%v + block%C(i-1,j+1,k  )%P%v)
      vel(6) = 0.25_R_P*(mean + block%C(i+1,j,k  )%P%v + block%C(i+1,j+1,k  )%P%v)
    case('k')
      mean = block%C(i,j,k)%P%v + block%C(i,j,k+1)%P%v

      vel(1) = block%C(i,j,k  )%P%v
      vel(2) = block%C(i,j,k+1)%P%v
      vel(3) = 0.25_R_P*(mean + block%C(i-1,j  ,k)%P%v + block%C(i-1,j  ,k+1)%P%v)
      vel(4) = 0.25_R_P*(mean + block%C(i+1,j  ,k)%P%v + block%C(i+1,j  ,k+1)%P%v)
      vel(5) = 0.25_R_P*(mean + block%C(i  ,j-1,k)%P%v + block%C(i  ,j-1,k+1)%P%v)
      vel(6) = 0.25_R_P*(mean + block%C(i  ,j+1,k)%P%v + block%C(i  ,j+1,k+1)%P%v)
    endselect
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine fv_interface_velocity

    !> Function for computing first partial derivative by means of Finite Volumes approach (Green's theorem).
    !> Using the Green's theorem the first partial derivative of \f$u\f$ in \f$\vec i\f$ direction could be written as
    !> \f$\frac{{\partial u}}{{\partial i}} =\frac{1}{V}\sum\limits_{f = 1}^6 {{u_f}\overrightarrow {{n_f}} \cdot \vec i\,{S_f}}\f$
    !> being \f$V\f$ the value of finite volume, \f$\overrightarrow {{n_f}}\f$ the outward unit normal of \f$f^{th}\f$ face which
    !> area is \f$S_f\f$.
    !> @note It is assumed that the finite volume is discretized by means of a hexahedron.
    !> @return \b fpd real(R_P) variable.
    function dudi_FV(u,nsi,v) result(fpd)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(R_P), intent(IN):: u(1:6)   !< Values of variable to be differentiated at each of 6 interfaces surrounding finite volume.
    real(R_P), intent(IN):: nsi(1:6) !< Area of 6 interfaces surrounding finite volume multiplied by normals projected along 'i'.
    real(R_P), intent(IN):: v        !< Value of finite volume.
    real(R_P)::             fpd      !< First partial derivative of 'u' in 'i' direction.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    fpd = (u(2)*nsi(2) - u(1)*nsi(1) + u(4)*nsi(4) - u(3)*nsi(3) + u(6)*nsi(6) - u(5)*nsi(5))/v
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction dudi_FV

    !> Subroutine for computing the shear stress tensor at center of a finite volume.
    subroutine shear_stress()
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(R_P)::         vol !< Value of finite volume.
    type(Type_Tensor):: rst !< Rate of strain tensor.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    ! computing the value of the finite volume
    select case(dir)
    case('i')
      vol=0.5_R_P*(block%C(i,j,k)%V + block%C(i+1,j,k)%V)
    case('j')
      vol=0.5_R_P*(block%C(i,j,k)%V + block%C(i,j+1,k)%V)
    case('k')
      vol=0.5_R_P*(block%C(i,j,k)%V + block%C(i,j,k+1)%V)
    endselect
    ! computing the gradient of velocity vector
    rst%x%x = dudi_FV(u=vel%x,nsi=NFS%x,v=vol)
    rst%x%y = dudi_FV(u=vel%x,nsi=NFS%y,v=vol)
    rst%x%z = dudi_FV(u=vel%x,nsi=NFS%z,v=vol)
    rst%y%x = dudi_FV(u=vel%y,nsi=NFS%x,v=vol)
    rst%y%y = dudi_FV(u=vel%y,nsi=NFS%y,v=vol)
    rst%y%z = dudi_FV(u=vel%y,nsi=NFS%z,v=vol)
    rst%z%x = dudi_FV(u=vel%z,nsi=NFS%x,v=vol)
    rst%z%y = dudi_FV(u=vel%z,nsi=NFS%y,v=vol)
    rst%z%z = dudi_FV(u=vel%z,nsi=NFS%z,v=vol)
    ! computing the shear stress tensor
    tau%x%x = block%global%adim%Re_inv*rst%x%x*4._R_P/3._R_P
    tau%x%y = block%global%adim%Re_inv*(rst%x%y+rst%y%x)
    tau%x%z = block%global%adim%Re_inv*(rst%x%z+rst%z%x)
    tau%y%x = tau%x%y
    tau%y%y = block%global%adim%Re_inv*rst%y%y*4._R_P/3._R_P
    tau%y%z = block%global%adim%Re_inv*(rst%y%z+rst%z%y)
    tau%z%x = tau%x%z
    tau%z%y = tau%y%z
    tau%z%z = block%global%adim%Re_inv*rst%z%z*4._R_P/3._R_P
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine shear_stress
endsubroutine fluxes_diffusive
endmodule Lib_Fluxes_Diffusive
