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
!> @ingroup Library
module Lib_Fluidynamic
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                   ! Integers and reals precision definition.
USE Data_Type_AMRBlock                             ! Definition of Type_AMRBlock.
USE Data_Type_BC                                   ! Definition of Type_BC.
USE Data_Type_Conservative                         ! Definition of Type_Conservative.
USE Data_Type_Global                               ! Definition of Type_Global.
USE Data_Type_Primitive                            ! Definition of Type_Primitive.
USE Data_Type_SBlock                               ! Definition of Type_SBlock.
USE Data_Type_Time                                 ! Definition of Type_Time.
USE Data_Type_Vector                               ! Definition of Type_Vector.
USE Lib_Fluxes_Convective, only: fluxes_convective ! Subroutine for computing convective fluxes.
USE Lib_Fluxes_Diffusive, only: fluxes_diffusive   ! Subroutine for computing diffusive fluxes.
USE Lib_IO_Misc                                    ! Procedures for IO and strings operations.
USE Lib_Math, only: interpolate1                   ! Function for computing linear interpolation.
USE Lib_Parallel, only: blockmap                   ! Local/global blocks map.
#ifdef PROFILING
USE Lib_Profiling                                  ! Procedures for profiling the code.
#endif
USE Lib_Runge_Kutta                                ! Runge-Kutta time integration library.
USE Lib_Thermodynamic_Laws_Ideal, only: a          ! Function for computing speed of sound.
#ifdef MPI2
USE Lib_Parallel, only: procmap, &                 ! Proc/blocks map.
                        Psendrecv                  ! Subroutine for send/receive Primitive variables.
USE MPI                                            ! MPI runtime library.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: prim2cons,primitive2conservative,cons2prim,conservative2primitive
public:: boundary_conditions
public:: solve_grl
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
integer(I1P):: flip = 0_I_P !< Flip-Flop flag for restart solution file.
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> Subroutine for converting primitive variables to conservative variables.
  !> The conversion laws are: \n
  !> - Partial densities\f$|_{conservative} =\left| \rho_s\right|_{primitive}\f$
  !> - Momentum\f$|_{conservative} = \left| \rho \vec V \right|_{primitive}\f$
  !> - Energy\f$|_{conservative}=\left|\frac{p}{\gamma-1}+\frac{1}{2}\rho
  !>   \left(V_x^2+V_y^2+V_z^2\right)\right|_{primitive}\f$
  !> @ingroup Lib_FluidynamicPublicProcedure
  pure subroutine prim2cons(prim,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive),    intent(IN)::    prim !< Primitive variables (see \ref Data_Type_Primitive::Type_Primitive
                                                !< "Type_Primitive" definition).
  type(Type_Conservative), intent(INOUT):: cons !< Conservative variables (see \ref Data_Type_Conservative::Type_Conservative
                                                !< "Type_Conservative" definition).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons%rs = prim%r
  cons%rv = prim%d*prim%v
  cons%re = prim%p/(prim%g-1._R_P) + 0.5_R_P*prim%d*sq_norm(prim%v)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prim2cons

  !> Subroutine for converting conservative variables to primitive variables.
  !> The conversion laws are: \n
  !> - Partial densities\f$|_{primitive} =\left| \rho_s\right|_{conservative}\f$
  !> - Density\f$|_{primitive} =\left|\sum_{s=1}^{Ns}\rho_s\right|_{conservative}\f$
  !> - Specific heats ratio\f$|_{primitive}=\left|\frac{\sum_{s=1}^{Ns}{\frac{\rho_s}{\rho}cp_s^0}}
  !>   {\sum_{s=1}^{Ns}{\frac{\rho_s}{\rho}cv_s^0}}\right|_{conservative}\f$
  !> - Velocity\f$|_{primitive} =\left|\frac{\vec momentum}{\rho}\right|_{conservative}\f$
  !> - Pressure\f$|_{primitive}=\left|(\gamma-1)\left[energy-\frac{1}{2}\rho
  !>   \left(V_x^2+V_y^2+V_z^2\right)\right]\right|_{conservative}\f$ \n
  !> where Ns is the number of initial species and \f$cp_s^0,cv_s^0\f$ are the initial specific heats.
  !> @ingroup Lib_FluidynamicPublicProcedure
  pure subroutine cons2prim(cp0,cv0,cons,prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN)::    cp0(:)                  !< Specific heat at constant p of initial species.
  real(R_P),               intent(IN)::    cv0(:)                  !< Specific heat at constant v of initial species.
  type(Type_Conservative), intent(IN)::    cons                    !< Conservative variables (see
                                                                   !< \ref Data_Type_Conservative::Type_Conservative
                                                                   !< "Type_Conservative" definition).
  type(Type_Primitive),    intent(INOUT):: prim                    !< Primitive variables (see
                                                                   !< \ref Data_Type_Primitive::Type_Primitive "Type_Primitive"
                                                                   !< definition).
  real(R_P)::                              c(1:size(cp0(:),dim=1)) !< Species concentration.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%r     = cons%rs
  prim%d     = sum(cons%rs(:))
  c          = prim%r/prim%d
  prim%g     = dot_product(c,cp0)/dot_product(c,cv0)
  prim%v     = cons%rv/prim%d
  prim%p     = (prim%g - 1._R_P)*(cons%re - 0.5_R_P*prim%d*(sq_norm(prim%v)))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine cons2prim

  !> Subroutine for converting primitive variables to conservative variables of all cells of a block.
  !> @note Only the inner cells of the block are converted.
  !> @ingroup Lib_FluidynamicPublicProcedure
  subroutine primitive2conservative(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block !< Block-level data (see \ref Data_Type_SBlock::Type_SBlock "Type_SBlock" definition).
  integer(I_P)::                      i,j,k !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !!$OMP PARALLEL DEFAULT(NONE) &
  !!$OMP PRIVATE(i,j,k)         &
  !!$OMP SHARED(block)
  !!$OMP DO
  do k=1,block%Nk
    do j=1,block%Nj
      do i=1,block%Ni
        call prim2cons(prim = block%P(i,j,k), cons = block%U(i,j,k))
      enddo
    enddo
  enddo
  !!$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine primitive2conservative

  !> Subroutine for converting conservative variables to primitive variables of all cells of a block.
  !> @note Only the inner cells of the block are converted.
  !> @ingroup Lib_FluidynamicPublicProcedure
  subroutine conservative2primitive(global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global),  intent(IN)::    global !< Global-level data (see \ref Data_Type_Global::Type_Global "Type_Global" definition).
  class(Type_SBlock), intent(INOUT):: block  !< Block-level data (see \ref Data_Type_SBlock::Type_SBlock "Type_SBlock" definition).
  integer(I_P)::                      i,j,k  !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !!$OMP PARALLEL DEFAULT(NONE) &
  !!$OMP PRIVATE(i,j,k)         &
  !!$OMP SHARED(global,block)
  !!$OMP DO
  do k=1,block%Nk
    do j=1,block%Nj
      do i=1,block%Ni
        call cons2prim(cp0 = global%cp0, cv0 = global%cv0, cons = block%U(i,j,k), prim = block%P(i,j,k))
      enddo
    enddo
  enddo
  !!$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine conservative2primitive

  !> Function for evaluating the local and global time step value by CFL condition.
  !> @ingroup Lib_FluidynamicPrivateProcedure
  subroutine compute_time(global,block,Dtmin)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global),  intent(IN)::    global                        !< Global-level data.
  class(Type_SBlock), intent(INOUT):: block                         !< Block-level data.
  real(R_P),          intent(OUT)::   Dtmin                         !< Minimum Dt.
  real(R_P)::                         vmax                          !< Maximum speed of waves.
  real(R_P)::                         ss                            !< Speed of sound.
  real(R_P)::                         vmiL,vmiR,vmjL,vmjR,vmkL,vmkR !< Dummy velocities.
  type(Type_Vector)::                 vm                            !< Dummy vectorial velocities.
  integer(I_P)::                      Ni,Nj,Nk,gc(1:6)              !< Temp var for storing block dimensions.
  integer(I_P)::                      i,j,k                         !< Space counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  gc = block%gc
  Ni = block%Ni
  Nj = block%Nj
  Nk = block%Nk
  ! computing the minimum Dt into the inner cells
  !!$OMP PARALLEL DEFAULT(NONE)                                   &
  !!$OMP PRIVATE(i,j,k,vmax,ss,vmiL,vmiR,vmjL,vmjR,vmkL,vmkR,vm)  &
  !!$OMP SHARED(gc,Ni,Nj,Nk,global,block,Dtmin)
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        ! computing the local speed of sound
        ss = a(p=block%P(i,j,k)%p,r=block%P(i,j,k)%d,g=block%P(i,j,k)%g)
        ! evaluating the maximum propagation speed of acoustic segnals multiplied for face area
        ! left i
        vm   = 0.5_R_P*(block%P(i-1,j,k)%v+block%P(i,j,k)%v)
        vmiL = (vm.dot.block%NFi(i-1,j,k))*block%Si(i-1,j,k)
        vmiL = abs(vmiL) + ss
        ! right i
        vm   = 0.5_R_P*(block%P(i,j,k)%v+block%P(i+1,j,k)%v)
        vmiR = (vm.dot.block%NFi(i,j,k))*block%Si(i,j,k)
        vmiR = abs(vmiR) + ss
        ! left j
        vm   = 0.5_R_P*(block%P(i,j-1,k)%v+block%P(i,j,k)%v)
        vmjL = (vm.dot.block%NFj(i,j-1,k))*block%Sj(i,j-1,k)
        vmjL = abs(vmjL) + ss
        ! right j
        vm   = 0.5_R_P*(block%P(i,j,k)%v+block%P(i,j+1,k)%v)
        vmjR = (vm.dot.block%NFj(i,j,k))*block%Sj(i,j,k)
        vmjR = abs(vmjR) + ss
        ! left k
        vm   = 0.5_R_P*(block%P(i,j,k-1)%v+block%P(i,j,k)%v)
        vmkL = (vm.dot.block%NFk(i,j,k-1))*block%Sk(i,j,k-1)
        vmkL = abs(vmkL) + ss
        ! right k
        vm   = 0.5_R_P*(block%P(i,j,k)%v+block%P(i,j,k+1)%v)
        vmkR = (vm.dot.block%NFk(i,j,k))*block%Sk(i,j,k)
        vmkR = abs(vmkR) + ss
        ! vmax
        vmax = max(vmiL,vmiR,vmjL,vmjR,vmkL,vmkR)
        block%Dt(i,j,k) = block%V(i,j,k)/vmax*global%CFL
      enddo
    enddo
  enddo
  ! computing minimum Dt
  !!$OMP SINGLE
  Dtmin = minval(block%Dt(1:Ni,1:Nj,1:Nk))
  !!$OMP END SINGLE
  ! ghost cells estrapolation: imposing the minum value of Dt
  ! left i frame
  !!$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=1-gc(1),0
        block%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  ! right i frame
  !!$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=Ni+1,Ni+gc(2)
        block%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  ! left j frame
  !!$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=1-gc(3),0
      do i=1-gc(1),Ni+gc(2)
        block%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  ! right j frame
  !!$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=Nj+1,Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
        block%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  ! left k frame
  !!$OMP DO
  do k=1-gc(5),0
    do j=1-gc(3),Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
        block%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  ! right k frame
  !!$OMP DO
  do k=Nk+1,Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
        block%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  !!$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_time

  !> Subroutine for computing the residuals. This the space operator. The residuals are stored in block%KS(s1) conservative
  !> variables.
  !> @ingroup Lib_FluidynamicPrivateProcedure
  subroutine residuals(s1,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),       intent(IN)::    s1                                    !< Current Runge-kutta stage.
  type(Type_Global),  intent(IN)::    global                                !< Global-level data.
  class(Type_SBlock), intent(INOUT):: block                                 !< Block-level data.
  type(Type_Conservative)::           Fic(0:block%Ni,1:block%Nj,1:block%Nk) !< I convective fluxes.
  type(Type_Conservative)::           Fjc(1:block%Ni,0:block%Nj,1:block%Nk) !< J convective fluxes.
  type(Type_Conservative)::           Fkc(1:block%Ni,1:block%Nj,0:block%Nk) !< K convective fluxes.
  type(Type_Conservative)::           Fid(0:block%Ni,1:block%Nj,1:block%Nk) !< I diffusive fluxes.
  type(Type_Conservative)::           Fjd(1:block%Ni,0:block%Nj,1:block%Nk) !< J diffusive fluxes.
  type(Type_Conservative)::           Fkd(1:block%Ni,1:block%Nj,0:block%Nk) !< K diffusive fluxes.
  integer(I1P)::                      gcu                                   !< Number of ghost cells used.
  integer(I_P)::                      Ni,Nj,Nk,Ns                           !< Temp var for storing block dims.
  integer(I1P)::                      gc(1:6)                               !< Temp var for storing ghost cells number.
  integer(I_P)::                      i,j,k,s                               !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  gc = block%gc
  Ni = block%Ni
  Nj = block%Nj
  Nk = block%Nk
  Ns = global%Ns
  call Fic%init(Ns=Ns)
  call Fid%init(Ns=Ns)
  call Fjc%init(Ns=Ns)
  call Fjd%init(Ns=Ns)
  call Fkc%init(Ns=Ns)
  call Fkd%init(Ns=Ns)
  ! computing convective fluxes
  !!$OMP PARALLEL DEFAULT(NONE) &
  !!$OMP PRIVATE(i,j,k,s)       &
  !!$OMP SHARED(s1,gc,gcu,Ni,Nj,Nk,Ns,global,block,Fic,Fjc,Fkc,Fid,Fjd,Fkd)
#ifndef NULi
  ! i direction
  !!$OMP SINGLE
  gcu = min(global%gco,gc(1),gc(2))
  !!$OMP END SINGLE
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      call fluxes_convective(gc  = gcu,                         &
                             N   = Ni,                          &
                             Ns  = Ns,                          &
                             cp0 = global%cp0,                  &
                             cv0 = global%cv0,                  &
                             NF  = block%NFi(0-gcu:Ni+gcu,j,k), &
                             P   = block%P  (1-gcu:Ni+gcu,j,k), &
                             F   = Fic      (    0:Ni    ,j,k))
    enddo
  enddo
#endif
#ifndef NULj
  ! j direction
  !!$OMP SINGLE
  gcu = min(global%gco,gc(3),gc(4))
  !!$OMP END SINGLE
  !!$OMP DO
  do k=1,Nk
    do i=1,Ni
      call fluxes_convective(gc  = gcu,                         &
                             N   = Nj,                          &
                             Ns  = Ns,                          &
                             cp0 = global%cp0,                  &
                             cv0 = global%cv0,                  &
                             NF  = block%NFj(i,0-gcu:Nj+gcu,k), &
                             P   = block%P  (i,1-gcu:Nj+gcu,k), &
                             F   = Fjc      (i,    0:Nj    ,k))
    enddo
  enddo
#endif
#ifndef NULk
  ! k direction
  !!$OMP SINGLE
  gcu = min(global%gco,gc(5),gc(6))
  !!$OMP END SINGLE
  !!$OMP DO
  do j=1,Nj
    do i=1,Ni
      call fluxes_convective(gc  = gcu,                         &
                             N   = Nk,                          &
                             Ns  = Ns,                          &
                             cp0 = global%cp0,                  &
                             cv0 = global%cv0,                  &
                             NF  = block%NFk(i,j,0-gcu:Nk+gcu), &
                             P   = block%P  (i,j,1-gcu:Nk+gcu), &
                             F   = Fkc      (i,j,    0:Nk    ))
    enddo
  enddo
#endif
  ! computing diffusive fluxes
  if (.not.global%inviscid) then
#ifndef NULi
    ! i direction
    !!$OMP DO
    do k=1,Nk
      do j=1,Nj
        do i=0,Ni
          call fluxes_diffusive(global=global,block=block,i=i,j=j,k=k,dir='i',F=Fid(i,j,k))
        enddo
      enddo
    enddo
#endif
#ifndef NULj
    ! j direction
    !!$OMP DO
    do k=1,Nk
      do j=0,Nj
        do i=1,Ni
          call fluxes_diffusive(global=global,block=block,i=i,j=j,k=k,dir='j',F=Fjd(i,j,k))
        enddo
      enddo
    enddo
#endif
#ifndef NULk
    ! k direction
    !!$OMP DO
    do k=0,Nk
      do j=1,Nj
        do i=1,Ni
          call fluxes_diffusive(global=global,block=block,i=i,j=j,k=k,dir='k',F=Fkd(i,j,k))
        enddo
      enddo
    enddo
#endif
  endif

  ! computing the residuals
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        ! overloaded operators form: not efficient!
        !block%KS(i,j,k,s1) =                                                                                     &
        !  (block%Si(i-1,j,  k  )*(Fic(i-1,j,  k  )+Fid(i-1,j,  k  )) - block%Si(i,j,k)*(Fic(i,j,k)+Fid(i,j,k)) + &
        !   block%Sj(i,  j-1,k  )*(Fjc(i,  j-1,k  )+Fjd(i,  j-1,k  )) - block%Sj(i,j,k)*(Fjc(i,j,k)+Fjd(i,j,k)) + &
        !   block%Sk(i,  j,  k-1)*(Fkc(i,  j,  k-1)+Fkd(i,  j,  k-1)) - block%Sk(i,j,k)*(Fkc(i,j,k)+Fkd(i,j,k))   &
        !  )/block%V(i,j,k)
        do s=1,Ns
          block%KS(i,j,k,s1)%rs(s) =                                             &
          (block%Si(i-1,j,  k  )*(Fic(i-1,j,  k  )%rs(s)+Fid(i-1,j,  k  )%rs(s))-&
           block%Si(i,  j,  k  )*(Fic(i,  j,  k  )%rs(s)+Fid(i,  j,  k  )%rs(s))+&
           block%Sj(i,  j-1,k  )*(Fjc(i,  j-1,k  )%rs(s)+Fjd(i,  j-1,k  )%rs(s))-&
           block%Sj(i,  j,  k  )*(Fjc(i,  j,  k  )%rs(s)+Fjd(i,  j,  k  )%rs(s))+&
           block%Sk(i,  j,  k-1)*(Fkc(i,  j,  k-1)%rs(s)+Fkd(i,  j,  k-1)%rs(s))-&
           block%Sk(i,  j,  k  )*(Fkc(i,  j,  k  )%rs(s)+Fkd(i,  j,  k  )%rs(s)) &
          )/block%V(i,j,k)
        enddo
        block%KS(i,j,k,s1)%rv =                                                                                        &
        (block%Si(i-1,j,  k  )*(Fic(i-1,j,  k  )%rv+Fid(i-1,j,  k  )%rv)-block%Si(i,j,k)*(Fic(i,j,k)%rv+Fid(i,j,k)%rv)+&
         block%Sj(i,  j-1,k  )*(Fjc(i,  j-1,k  )%rv+Fjd(i,  j-1,k  )%rv)-block%Sj(i,j,k)*(Fjc(i,j,k)%rv+Fjd(i,j,k)%rv)+&
         block%Sk(i,  j,  k-1)*(Fkc(i,  j,  k-1)%rv+Fkd(i,  j,  k-1)%rv)-block%Sk(i,j,k)*(Fkc(i,j,k)%rv+Fkd(i,j,k)%rv) &
        )/block%V(i,j,k)
        block%KS(i,j,k,s1)%re =                                                                                        &
        (block%Si(i-1,j,  k  )*(Fic(i-1,j,  k  )%re+Fid(i-1,j,  k  )%re)-block%Si(i,j,k)*(Fic(i,j,k)%re+Fid(i,j,k)%re)+&
         block%Sj(i,  j-1,k  )*(Fjc(i,  j-1,k  )%re+Fjd(i,  j-1,k  )%re)-block%Sj(i,j,k)*(Fjc(i,j,k)%re+Fjd(i,j,k)%re)+&
         block%Sk(i,  j,  k-1)*(Fkc(i,  j,  k-1)%re+Fkd(i,  j,  k-1)%re)-block%Sk(i,j,k)*(Fkc(i,j,k)%re+Fkd(i,j,k)%re) &
        )/block%V(i,j,k)
      enddo
    enddo
  enddo
#ifdef NULi
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        block%KS(i,j,k,s1)%rv%x = 0._R_P
      enddo
    enddo
  enddo
#endif
#ifdef NULj
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        block%KS(i,j,k,s1)%rv%y = 0._R_P
      enddo
    enddo
  enddo
#endif
#ifdef NULk
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        block%KS(i,j,k,s1)%rv%z = 0._R_P
      enddo
    enddo
  enddo
#endif
  !!$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine residuals

  !> Subroutine for imposing the boundary conditions of blocks of grid level "l" updating the ghost cells.
  !> @note Considering the ghost cell \f$ c^0 \f$ and the inner one \f$ c^1 \f$ (being \f$ c^N \f$ the other bound)
  !> the available boundary conditions are:
  !> - \b REF: reflective boundary condition \f$ \begin{array}{*{20}{c}} P^0 = P^1 \\ P_{\vec v\cdot \vec n}^0 =
  !>           - P_{\vec v\cdot \vec n}^1\end{array} \f$;
  !> - \b EXT: extrapolation boundary condition \f$ P^0 = P^1 \f$;
  !> - \b PER: periodic boundary condition \f$ P^0 = P^N \f$;
  !> - \b ADJ: adjacent (cell) boundary condition \f$ P^0 = P^a \f$ where "a" is the adjacent cell (b,i,j,k indexes must be
  !>           specified);
  !> - \b IN1: supersonic inflow steady boundary condition \f$ P^0 = P^{in1} \f$. \n
  !> where \f$P\f$ are the primitive variables.
  !> @ingroup Lib_FluidynamicPublicProcedure
  subroutine boundary_conditions(l,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),       intent(IN)::    l                  !< Current grid level.
  type(Type_Global),  intent(IN)::    global             !< Global-level data.
  class(Type_SBlock), intent(INOUT):: block(1:global%Nb) !< Block-level data.
  integer(I_P)::                      Ni,Nj,Nk,gc(1:6)   !< Temporary var for storing block dimensions.
  integer(I_P)::                      b,i,j,k            !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef MPI2
  ! doing the multi-processes communications if necessary
  call Psendrecv(l=l,global=global,block=block)
#endif
  do b=1,global%Nb
    gc = block(b)%gc
    Ni = block(b)%Ni
    Nj = block(b)%Nj
    Nk = block(b)%Nk
    !!$OMP PARALLEL DEFAULT(NONE) &
    !!$OMP PRIVATE(i,j,k)         &
    !!$OMP SHARED(b,l,Ni,Nj,Nk,gc,global,block)
    !!$OMP DO
    do k=1,Nk
      do j=1,Nj
        ! left i
        select case(block(b)%BCi(0,j,k)%tp)
        case(bc_ext)
          call set_ext_l(gc=gc(1),ic=0,N=Ni,P=block(b)%P(1-gc(1):0+gc(1),j,k))
        case(bc_ref)
          call set_ref_l(gc=gc(1),ic=0,N=Ni,NF=block(b)%NFi(0,j,k),P=block(b)%P(1-gc(1):0+gc(1),j,k))
        case(bc_per)
          call set_per(gc=gc(1),ic=0,N=Ni,boundary='l',P=block(b)%P(1-gc(1):Ni+gc(1),j,k))
        case(bc_adj)
          call set_adj(gc=gc(1),bc=block(b)%BCi(1-gc(1):0,j,k),P=block(b)%P(1-gc(1):0,j,k))
        case(bc_in1)
          call set_in1(gc=gc(1),bc=block(b)%BCi(1-gc(1):0,j,k),P=block(b)%P(1-gc(1):0,j,k))
        endselect
        ! right i
        select case(block(b)%BCi(Ni,j,k)%tp)
        case(bc_ext)
          call set_ext_r(gc=gc(2),ic=0,N=Ni,P=block(b)%P(Ni-gc(2):Ni+gc(2),j,k))
        case(bc_ref)
          call set_ref_r(gc=gc(2),ic=0,N=Ni,NF=block(b)%NFi(Ni,j,k),P=block(b)%P(Ni-gc(2):Ni+gc(2),j,k))
        case(bc_per)
          call set_per(gc=gc(2),ic=0,N=Ni,boundary='r',P=block(b)%P(1-gc(2):Ni+gc(2),j,k))
        case(bc_adj)
          call set_adj(gc=gc(2),bc=block(b)%BCi(Ni:Ni+gc(2)-1,j,k),P=block(b)%P(Ni+1:Ni+gc(2),j,k))
        case(bc_in1)
          call set_in1(gc=gc(2),bc=block(b)%BCi(Ni:Ni+gc(2)-1,j,k),P=block(b)%P(Ni+1:Ni+gc(2),j,k))
        endselect
      enddo
    enddo
    !!$OMP DO
    do k=1,Nk
      do i=1,Ni
        ! left j
        select case(block(b)%BCj(i,0,k)%tp)
        case(bc_ext)
          call set_ext_l(gc=gc(3),ic=0,N=Nj,P=block(b)%P(i,1-gc(3):0+gc(3),k))
        case(bc_ref)
          call set_ref_l(gc=gc(3),ic=0,N=Nj,NF=block(b)%NFj(i,0,k),P=block(b)%P(i,1-gc(3):0+gc(3),k))
        case(bc_per)
          call set_per(gc=gc(3),ic=0,N=Nj,boundary='l',P=block(b)%P(i,1-gc(3):Nj+gc(3),k))
        case(bc_adj)
          call set_adj(gc=gc(3),bc=block(b)%BCj(i,1-gc(3):0,k),P=block(b)%P(i,1-gc(3):0,k))
        case(bc_in1)
          call set_in1(gc=gc(3),bc=block(b)%BCj(i,1-gc(3):0,k),P=block(b)%P(i,1-gc(3):0,k))
        endselect
        ! right j
        select case(block(b)%BCj(i,Nj,k)%tp)
        case(bc_ext)
          call set_ext_r(gc=gc(4),ic=0,N=Nj,P=block(b)%P(i,Nj-gc(4):Nj+gc(4),k))
        case(bc_ref)
          call set_ref_r(gc=gc(4),ic=0,N=Nj,NF=block(b)%NFj(i,Nj,k),P=block(b)%P(i,Nj-gc(4):Nj+gc(4),k))
        case(bc_per)
          call set_per(gc=gc(4),ic=0,N=Nj,boundary='r',P=block(b)%P(i,1-gc(4):Nj+gc(4),k))
        case(bc_adj)
          call set_adj(gc=gc(4),bc=block(b)%BCj(i,Nj:Nj+gc(4)-1,k),P=block(b)%P(i,Nj+1:Nj+gc(4),k))
        case(bc_in1)
          call set_in1(gc=gc(4),bc=block(b)%BCj(i,Nj:Nj+gc(4)-1,k),P=block(b)%P(i,Nj+1:Nj+gc(4),k))
        endselect
      enddo
    enddo
    !!$OMP DO
    do j=1,Nj
      do i=1,Ni
        ! left k
        select case(block(b)%BCk(i,j,0)%tp)
        case(bc_ext)
          call set_ext_l(gc=gc(5),ic=0,N=Nk,P=block(b)%P(i,j,1-gc(5):0+gc(5)))
        case(bc_ref)
          call set_ref_l(gc=gc(5),ic=0,N=Nk,NF=block(b)%NFk(i,j,0),P=block(b)%P(i,j,1-gc(5):0+gc(5)))
        case(bc_per)
          call set_per(gc=gc(5),ic=0,N=Nk,boundary='l',P=block(b)%P(i,j,1-gc(5):Nk+gc(5)))
        case(bc_adj)
          call set_adj(gc=gc(5),bc=block(b)%BCk(i,j,1-gc(5):0),P=block(b)%P(i,j,1-gc(5):0))
        case(bc_in1)
          call set_in1(gc=gc(5),bc=block(b)%BCk(i,j,1-gc(5):0),P=block(b)%P(i,j,1-gc(5):0))
        endselect
        ! right k
        select case(block(b)%BCk(i,j,Nk)%tp)
        case(bc_ext)
          call set_ext_r(gc=gc(6),ic=0,N=Nk,P=block(b)%P(i,j,Nk-gc(6):Nk+gc(6)))
        case(bc_ref)
          call set_ref_r(gc=gc(6),ic=0,N=Nk,NF=block(b)%NFk(i,j,Nk),P=block(b)%P(i,j,Nk-gc(6):Nk+gc(6)))
        case(bc_per)
          call set_per(gc=gc(6),ic=0,N=Nk,boundary='r',P=block(b)%P(i,j,1-gc(6):Nk+gc(6)))
        case(bc_adj)
          call set_adj(gc=gc(6),bc=block(b)%BCk(i,j,Nk:Nk+gc(6)-1),P=block(b)%P(i,j,Nk+1:Nk+gc(6)))
        case(bc_in1)
          call set_in1(gc=gc(6),bc=block(b)%BCk(i,j,Nk:Nk+gc(6)-1),P=block(b)%P(i,j,Nk+1:Nk+gc(6)))
        endselect
      enddo
    enddo
    !!$OMP END PARALLEL
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Subroutine for imposing extrapolation of ghost cells from internal ones (left boundary).
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the array \b P is declared as assumed-shape
    !> with only the lower bound defined. Its extentions is: P [1-gc:0+gc].
    pure subroutine set_ext_l(gc,ic,N,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc       !< Number of ghost cells.
    integer(I_P),         intent(IN)::    ic       !< Number of internal cells used for extrapolation (1 or gc).
    integer(I_P),         intent(IN)::    N        !< Number of internal cells.
    type(Type_Primitive), intent(INOUT):: P(1-gc:) !< Primitive variables [1-gc:0+gc].
    integer(I_P)::                        i        !< Cell counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (ic==1.OR.N<gc) then
      ! extrapolation using only the cell 1
      do i=1-gc,0
        P(i) = P(1)
      enddo
    else
      ! extrapolation using the cells 1,2,...,gc
      do i=1-gc,0
        P(i) = P(-i+1)
      enddo
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_ext_l

    !> @brief Subroutine for imposing extrapolation of ghost cells from internal ones (right boundary).
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the array \b P is declared as assumed-shape
    !> with only the lower bound defined. Its extentions is: P [N-gc:N+gc].
    pure subroutine set_ext_r(gc,ic,N,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc       !< Number of ghost cells.
    integer(I_P),         intent(IN)::    ic       !< Number of internal cells used for extrapolation (1 or gc).
    integer(I_P),         intent(IN)::    N        !< Number of internal cells.
    type(Type_Primitive), intent(INOUT):: P(N-gc:) !< Primitive variables [N-gc:N+gc].
    integer(I_P)::                        i        !< Cell counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (ic==1.OR.N<gc) then
      ! extrapolation using only the cell N
      do i=N+1,N+gc
        P(i) = P(N)
      enddo
    else
      ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
      do i=N+1,N+gc
        P(i) = P(N+1-(i-N))
      enddo
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_ext_r

    !> @brief Subroutine for imposing reflective boundary conditions (left boundary).
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the array \b P is declared as assumed-shape
    !> with only the lower bound defined. Its extentions is: P [1-gc:0+gc].
    pure subroutine set_ref_l(gc,ic,N,NF,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc       !< Number of ghost cells.
    integer(I_P),         intent(IN)::    ic       !< Number of internal cells used for extrapolation (1 or gc).
    integer(I_P),         intent(IN)::    N        !< Number of internal cells.
    type(Type_Vector),    intent(IN)::    NF       !< Left face normal.
    type(Type_Primitive), intent(INOUT):: P(1-gc:) !< Left section of primitive variables [1-gc:0+gc].
    integer(I_P)::                        i        !< Cell counter.
    type(Type_Vector)::                   vr       !< Reflected velocity vector.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (ic==1.OR.N<gc) then
      ! reflection using only the cell 1
      vr = P(1)%v - (2._R_P*(P(1)%v.paral.NF)) ! reflected velocity
      do i=1-gc,0
        P(i)%r = P(1)%r
        P(i)%v = vr
        P(i)%p = P(1)%p
        P(i)%d = P(1)%d
        P(i)%g = P(1)%g
      enddo
    else
      ! reflection using the cells 1,2,...,gc
      do i=1-gc,0
        vr = P(-i+1)%v - (2._R_P*(P(-i+1)%v.paral.NF)) ! reflected velocity
        P(i)%r = P(-i+1)%r
        P(i)%v = vr
        P(i)%p = P(-i+1)%p
        P(i)%d = P(-i+1)%d
        P(i)%g = P(-i+1)%g
      enddo
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_ref_l

    !> @brief Subroutine for imposing reflective boundary conditions (right boundary).
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the array \b P is declared as assumed-shape
    !> with only the lower bound defined. Its extentions is: P [N-gc:N+gc].
    pure subroutine set_ref_r(gc,ic,N,NF,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc       !< Number of ghost cells.
    integer(I_P),         intent(IN)::    ic       !< Number of internal cells used for extrapolation (1 or gc).
    integer(I_P),         intent(IN)::    N        !< Number of internal cells.
    type(Type_Vector),    intent(IN)::    NF       !< Right face normal.
    type(Type_Primitive), intent(INOUT):: P(N-gc:) !< Right section of primitive variables [N-gc:N+gc].
    integer(I_P)::                        i        !< Cell counter.
    type(Type_Vector)::                   vr       !< Reflected velocity vector.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (ic==1.OR.N<gc) then
      ! reflection using only the cell N
      vr = P(N)%v - (2._R_P*(P(N)%v.paral.NF)) ! reflected velocity
      do i=N+1,N+gc
        P(i)%r = P(N)%r
        P(i)%v = vr
        P(i)%p = P(N)%p
        P(i)%d = P(N)%d
        P(i)%g = P(N)%g
      enddo
    else
      ! reflection using the cells N-gc,N-gc+1,N-gc+2,...,N
      do i=N+1,N+gc
        vr = P(N+1-(i-N))%v - (2._R_P*(P(N+1-(i-N))%v.paral.NF)) ! reflected velocity
        P(i)%r = P(N+1-(i-N))%r
        P(i)%v = vr
        P(i)%p = P(N+1-(i-N))%p
        P(i)%d = P(N+1-(i-N))%d
        P(i)%g = P(N+1-(i-N))%g
      enddo
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_ref_r

    !> @brief Subroutine for imposing periodic boundary conditions.
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the array \b P is declared as assumed-shape
    !> with only the lower bound defined. Its extentions is: P [1-gc:N+gc].
    pure subroutine set_per(gc,ic,N,boundary,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc       !< Number of ghost cells.
    integer(I_P),         intent(IN)::    ic       !< Number of internal cells used for extrapolation (1 or gc).
    integer(I_P),         intent(IN)::    N        !< Number of internal cells.
    character(1),         intent(IN)::    boundary !< Boundary left ('l') or right ('r').
    type(Type_Primitive), intent(INOUT):: P(1-gc:) !< Left section of primitive variables [1-gc:N+gc].
    integer(I_P)::                        i        !< Cell counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (boundary=='l') then
      if (ic==1.OR.N<gc) then
        ! extrapolation using only the cell N
        do i=1-gc,0
          P(i) = P(N)
        enddo
      else
        ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
        do i=1-gc,0
          P(i) = P(i+N)
        enddo
      endif
    endif
    if (boundary=='r') then
      if (ic==1.OR.N<gc) then
        ! extrapolation using only the cell 1
        do i=N+1,N+gc
          P(i) = P(1)
        enddo
      else
        ! extrapolation using the cells 1,2,...,gc
        do i=N+1,N+gc
          P(i) = P(i-N)
        enddo
      endif
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_per

    !> @brief Subroutine for imposing adjacent boundary conditions.
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b bc and \b P are declared as
    !> assumed-shape with only the lower bound defined. Their extentions are: bc [1-gc:0], P [1-gc:0].
    !> @note When this subroutine is called for a 'right' (Ni,Nj,Nk) boundary the section of arrays BC and P must be
    !> properly remapped: in the section [1-gc:0] must be passed the actual section [N:N+gc-1] for BC and [N+1:N+gc] for P.
    pure subroutine set_adj(gc,bc,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc        !< Number of ghost cells.
    type(Type_BC),        intent(IN)::    bc(1-gc:) !< Boundary conditions infos           [1-gc:0].
    type(Type_Primitive), intent(INOUT):: P (1-gc:) !< Left section of primitive variables [1-gc:0].
    integer(I_P)::                        i,b       !< Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
#ifdef MPI2
    ! for multi-processes simulation it is possible that the adjacent cell is in other processes than the actual and in case the
    ! data have been already exchanged by the MPI subroutine Psendrecv
    if (procmap(bc(0)%adj%b)/=global%myrank) return
#endif
    do i=1-gc,0
      b = minloc(array=blockmap,dim=1,mask=blockmap==bc(i)%adj%b)
      P(i) = block(b)%P(bc(i)%adj%i,bc(i)%adj%j,bc(i)%adj%k)
      !call get_cell_lbijk(P = P(i),        &
      !                    l = l,           &
      !                    b = b,           &
      !                    i = BC(i)%adj%i, &
      !                    j = BC(i)%adj%j, &
      !                    k = BC(i)%adj%k)
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_adj

    !> @brief Subroutine for imposing inflow 1 boundary conditions.
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b bc and \b P are declared as
    !> assumed-shape with only the lower bound defined. Their extentions are: bc [1-gc:0], P [1-gc:0].
    !> @note When this subroutine is called for a 'right' (Ni,Nj,Nk) boundary the section of arrays BC and P must be
    !> properly remapped: in the section [1-gc:0] must be passed the actual section [N:N+gc-1] for BC and [N+1:N+gc] for P.
    pure subroutine set_in1(gc,bc,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc        !< Number of ghost cells.
    type(Type_BC),        intent(IN)::    bc(1-gc:) !< Boundary conditions                 [1-gc:0].
    type(Type_Primitive), intent(INOUT):: P (1-gc:) !< Left section of primitive variables [1-gc:0].
    integer(I_P)::                        i         !< Cell counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    do i=1-gc,0
      P(i) = global%in1(bc(i)%inf)
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_in1
  endsubroutine boundary_conditions

  !> Subroutine for summing Runge-Kutta stages for updating primitive variables (block%P).
  !> @ingroup Lib_FluidynamicPrivateProcedure
  subroutine rk_stages_sum(s1,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),       intent(IN)::    s1     !< Current Runge-Kutta stage.
  type(Type_Global),  intent(IN)::    global !< Global-level data.
  class(Type_SBlock), intent(INOUT):: block  !< Block-level data.
  type(Type_Conservative)::           Ud     !< Dummy conservative variables.
  integer(I8P)::                      i,j,k  !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !!$OMP PARALLEL DEFAULT(NONE) &
  !!$OMP PRIVATE(i,j,k,Ud)      &
  !!$OMP SHARED(global,block,s1)
  call Ud%init(Ns = global%Ns)
  !!$OMP DO
  do k=1,block%Nk
    do j=1,block%Nj
      do i=1,block%Ni
        call rk_stage(s1=s1,Dt=block%Dt(i,j,k),Un=block%U(i,j,k),KS=block%KS(i,j,k,1:s1-1),KS1=Ud)
        call cons2prim(cp0 = global%cp0, cv0 = global%cv0, cons = Ud, prim = block%P(i,j,k))
      enddo
    enddo
  enddo
  !!$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rk_stages_sum

  !> Subroutine for computing Runge-Kutta one time step integration.
  !> @ingroup Lib_FluidynamicPrivateProcedure
  subroutine rk_time_integration(global,block,RU)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global),  intent(IN)::    global          !< Global-level data.
  class(Type_SBlock), intent(INOUT):: block           !< Block-level data.
  real(R_P),          intent(OUT)::   RU(1:global%Nc) !< NormL2 of residuals of conservative variables.
  type(Type_Conservative)::           Ud,R            !< Dummy conservative variables.
  integer(I8P)::                      i,j,k           !< counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  RU = 0._R_P
  !!$OMP PARALLEL DEFAULT(NONE) &
  !!$OMP PRIVATE(i,j,k,Ud,R)    &
  !!$OMP SHARED(global,block)   &
  !!$OMP REDUCTION(+: RU)
  call Ud%init(Ns = global%Ns)
  call R%init( Ns = global%Ns)
  !!$OMP DO
  do k=1,block%Nk
    do j=1,block%Nj
      do i=1,block%Ni
        call rk_time_integ(Dt=block%Dt(i,j,k),Un=block%U(i,j,k),KS=block%KS(i,j,k,1:global%rk_ord),Unp1=Ud)
        R%rs = (Ud%rs - block%U(i,j,k)%rs)/block%Dt(i,j,k) ; R%rs = R%rs*R%rs
        R%rv = (Ud%rv - block%U(i,j,k)%rv)/block%Dt(i,j,k) ; R%rv = R%rv*R%rv
        R%re = (Ud%re - block%U(i,j,k)%re)/block%Dt(i,j,k) ; R%re = R%re*R%re
        RU = RU + R%cons2array()
        block%U(i,j,k) = Ud
      enddo
    enddo
  enddo
  !!$OMP DO
  do i=1,global%Nc
   RU(i) = sqrt(RU(i))
  enddo
  !!$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rk_time_integration

  !> @ingroup Lib_FluidynamicPublicProcedure
  !> Subroutine for solving (performing one time step integration) the conservation equations for grid level "l".
  subroutine solve_grl(l,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),       intent(IN)::    l                  !< Current grid level.
  type(Type_Global),  intent(INOUT):: global             !< Global-level data (see \ref Data_Type_Global::Type_Global
                                                         !< "Type_Global" definition).
  class(Type_SBlock), intent(INOUT):: block(1:global%Nb) !< Block-level data (see \ref Data_Type_SBlock::Type_SBlock
                                                         !< "Type_SBlock" definition).
  real(R_P)::       Dtmin(1:global%Nb)                   !< Min t step of actual process for each blk.
  real(R_P)::       DtminL                               !< Min t step of actual process over all blks.
  real(R_P)::       gDtmin                               !< Global (all processes/all blks) min t step.
  real(R_P)::       RU  (1:global%Nc,1:global%Nb)        !< NormL2 of conservartive residuals.
  real(R_P)::       mRU (1:global%Nc)                    !< Maximum of RU of actual process.
  real(R_P)::       gmRU(1:global%Nc)                    !< Global (all processes) maximum of RU.
  integer(I_P)::    err                                  !< Error trapping flag: 0 no errors, >0 errors.
  integer(I_P)::    b                                    !< Blocks counter.
  integer(I_P)::    s1                                   !< Runge-Kutta stages counters.
  real(R_P)::       sec_elp                              !< Seconds elapsed from the simulation start.
  real(R_P)::       sec_res                              !< Seconds residual from the simulation end.
  type(Type_Time):: time_elp                             !< Time elapsed (in days,hours,min... format).
  type(Type_Time):: time_res                             !< Time residual (in days,hours,min... format).
  real(R_P)::       progress                             !< Status (%) of simulation progress.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! converting conservative variables to primitive ones
#ifdef PROFILING
  call profile(p=2,pstart=.true.,myrank=global%myrank)
#endif
  do b=1,global%Nb
    call conservative2primitive(global = global, block = block(b))
  enddo
#ifdef PROFILING
  call profile(p=2,pstop=.true.,myrank=global%myrank)
#endif

  ! imposing the boundary conditions
#ifdef PROFILING
  call profile(p=3,pstart=.true.,myrank=global%myrank)
#endif
  call boundary_conditions(l = l, global = global, block = block)
#ifdef PROFILING
  call profile(p=3,pstop=.true.,myrank=global%myrank)
#endif

  ! saving the actual solution
  if (global%file%sol_out>0) then
    if ((mod(global%n,global%file%sol_out)==0).OR.(global%t==global%Tmax).OR.(global%n==global%Nmax)) then
      do b=1,global%Nb
        err= block(b)%save_fluid(filename=file_name(basename=trim(global%file%Path_OutPut)//global%file%File_Sol,&
                                                    suffix='.sol',blk=blockmap(b),grl=l,n=global%n),             &
                                 global=global)
      enddo
    endif
  endif

  ! saving the restart solution file
  if ((mod(global%n,global%file%restart_out)==0).OR.(global%t==global%Tmax).OR.(global%n==global%Nmax)) then
    flip = 1_I1P - flip
    do b=1,global%Nb
      err= block(b)%save_fluid(filename=file_name(basename=trim(global%file%Path_OutPut)//global%file%File_Sol,&
                                                  suffix='.sol',blk=blockmap(b),grl=l,flip=flip),              &
                                                  global=global)
    enddo
  endif

  ! updating time varying variables: Dt,Dtmin
  global%n = global%n + 1_I8P
#ifdef PROFILING
  call profile(p=4,pstart=.true.,myrank=global%myrank)
#endif
  do b=1,global%Nb
    call compute_time(global=global,block=block(b),Dtmin=Dtmin(b))
  enddo
#ifdef PROFILING
  call profile(p=4,pstop=.true.,myrank=global%myrank)
#endif
  DtminL = minval(Dtmin)
#ifdef MPI2
  ! for multi-processes simulation all processes must exchange their DtminL for computing the global variables
  call MPI_ALLREDUCE(DtminL,gDtmin,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,err)
#else
  ! for single processes DtminL are already global variables
  gDtmin = DtminL
#endif
  if (global%unsteady) then ! for an unsteady accurate simulation each cell is updated by means of global minimum time step
    ! control for the last iterate
    if (global%Nmax<=0) then
      if ((global%t+gDtmin)>global%Tmax) then
        ! the global minimum time step is so high that the last iteration will go over Tmax
        ! it is decreased both for in order to achieve exactly Tmax
        gDtmin=abs(global%Tmax-global%t)
      endif
    endif
    global%t = global%t + gDtmin
    do b=1,global%Nb
      block(b)%Dt = gDtmin
    enddo
  endif

  ! updating console
  if (global%myrank==0) then
    if ((mod(global%n,global%file%screen_out)==0).OR.(global%t==global%Tmax).OR.(global%n==global%Nmax).OR.(global%n==1)) then
      sec_elp=Crono(instant0=.true.)
      if (global%Nmax>0) then
        progress = global%n*100/(global%Nmax*1._R_P)
        sec_res  = global%Nmax*sec_elp/global%n - sec_elp
      elseif (global%Tmax>0._R_P) then
        progress = 100*global%t/global%Tmax
        sec_res  = global%Tmax*sec_elp/global%t - sec_elp
      else
        progress = 0._R_P
        sec_res  = 0._R_P
      endif
      time_elp = Seconds_To_Time(sec_elp)
      time_res = Seconds_To_Time(sec_res)
      write(stdout,'(A)',                      iostat=err)'----------------------------------------------------------------------'
      write(stdout,'(A39,I30)',                iostat=err)' Current grid level                 l: ',l
      write(stdout,'(A39,23X,F6.2,A)',         iostat=err)' Simulation progress                p: ',progress,'%'
      write(stdout,'(A39,I30)',                iostat=err)' Current step number                n: ',global%n
      write(stdout,'(A39,ES30.12)',            iostat=err)' Current simulation time            t: ',global%t
      write(stdout,'(A39,ES30.12)',            iostat=err)' Current time step value          gDt: ',gDtmin
      write(stdout,*)
      write(stdout,'(10X,A15,20X,A15)',        iostat=err)' Elapsed Time','Residual Time'
      write(stdout,'(10X,A9,I6,20X,A9,I6)',    iostat=err)' Days    =',time_elp%Days,   'Days    =',time_res%Days
      write(stdout,'(10X,A9,I6,20X,A9,I6)',    iostat=err)' Hours   =',time_elp%Hours,  'Hours   =',time_res%Hours
      write(stdout,'(10X,A9,I6,20X,A9,I6)',    iostat=err)' Minutes =',time_elp%Minutes,'Minutes =',time_res%Minutes
      write(stdout,'(10X,A9,F6.1,20X,A9,F6.1)',iostat=err)' Seconds =',time_elp%Seconds,'Seconds =',time_res%Seconds
      write(stdout,'(A)',                      iostat=err)'----------------------------------------------------------------------'
      write(stdout,*)
    endif
  endif

  ! evaluating the Runge-Kutta stages
  ! Runge-Kutta stages initialization
  do b=1,global%Nb
    block(b)%KS = 0._R_P
  enddo
  do s1=1,global%rk_ord
    if (s1>1) then
      ! summing the stages up to s1-1 for update P
      do b=1,global%Nb
          call rk_stages_sum(s1=s1,global=global,block=block(b))
      enddo
      ! imposing the boundary conditions
      call boundary_conditions(l = l, global = global, block = block)
    endif
    ! computing the s1-th Runge-Kutta stage: K_s1=R(Un+sum_s2=1^s1-1(Dt*rk_c2(s1,s2)*K_s2))
#ifdef PROFILING
    call profile(p=5,pstart=.true.,myrank=global%myrank)
#endif
    do b=1,global%Nb
      call residuals(s1=s1,global=global,block=block(b))
    enddo
#ifdef PROFILING
    call profile(p=5,pstop=.true.,myrank=global%myrank)
#endif
  enddo
  ! Runge-Kutta time integration
#ifdef PROFILING
    call profile(p=6,pstart=.true.,myrank=global%myrank)
#endif
  do b=1,global%Nb
    call rk_time_integration(global=global,block=block(b),RU=RU(:,b))
  enddo
#ifdef PROFILING
    call profile(p=6,pstop=.true.,myrank=global%myrank)
#endif
  ! finding the maximum value of residuals of actual process
  mRU = maxval(RU,dim=2)
#ifdef MPI2
  ! for multi-processes simulation all processes must exchange their mRU for computing the global gmRU
  call MPI_ALLREDUCE(mRU,gmRU,global%Nc,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,err)
#else
  ! for single processes mRU is already global gmRU
  gmRU = mRU
#endif

  ! updating the log file of residuals
  if (global%myrank==0) then
    if ((mod(global%n,global%file%probe_out)==0).OR.(global%t==global%Tmax).OR.(global%n==global%Nmax).OR.(global%n==1)) then
      write(global%file%unit_res,trim(global%file%varform_res))global%n,global%t,(gmRU(b),b=1,global%Nc)
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine solve_grl
endmodule Lib_Fluidynamic
