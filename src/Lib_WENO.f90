!> @ingroup PrivateVarPar
!> @{
!> @defgroup Lib_WENOPrivateVarPar Lib_WENO
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_WENOPublicProcedure Lib_WENO
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_WENOPrivateProcedure Lib_WENO
!> @}

!> This module contains the definition of procedures for computing WENO reconstruction with a user-defined \f$P^{th}\f$ order.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Library
module Lib_WENO
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision       ! Integers and reals precision definition.
USE Data_Type_AMRBlock ! Definition of Type_AMRBlock.
USE Data_Type_Global   ! Definition of Type_Global.
USE Data_Type_SBlock   ! Definition of Type_SBlock.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: weno_init
public:: noweno_central
public:: weno_optimal
public:: weno
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! WENO coefficients (S = number of stencils used)
!> @ingroup Lib_WENOPrivateVarPar
!> @{
real(R_P), allocatable:: weno_c(:,:)   !< Central difference coefficients    [1:2,1:2*S].
real(R_P), allocatable:: weno_a(:,:)   !< Optimal weights                    [1:2,0:S-1].
real(R_P), allocatable:: weno_p(:,:,:) !< Polynomials coefficients           [1:2,0:S-1,0:S-1].
real(R_P), allocatable:: weno_d(:,:,:) !< Smoothness indicators coefficients [0:S-1,0:S-1,0:S-1].
real(R_P)::              weno_eps      !< Parameter for avoiding divided by zero when computing smoothness indicators.
integer(I_P)::           weno_odd      !< Constant for distinguishing between odd and even number of stencils (mod(S,2)).
integer(I_P)::           weno_exp      !< Exponent for growing the diffusive part of weights.
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_WENOPublicProcedure
  !> @{
  !> @brief Subroutine for initialization of WENO coefficients.
  subroutine weno_init(global,block,S)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global),  intent(IN):: global             !< Global-level data.
  class(Type_SBlock), intent(IN):: block(1:global%Nb) !< Block-level data.
  integer(I1P),       intent(IN):: S                  !< Number of stencils used.
  integer(I_P)::                   b                  !< Blocks counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initialize weno_exp
  weno_exp = S
  if (S>4) weno_exp = S - 1
  ! computing weno_odd
  weno_odd = mod(S,2_I1P)
  ! computing a reasonable value of weno_eps according to the (finest) grid spacing
  weno_eps = MaxR_P
  do b=1,global%Nb
    call get_min_sstep(block = block(b), ss = weno_eps)
  enddo
  weno_eps = weno_eps**(3_I1P*S-4_I1P)
  if (global%myrank==0) then
    write(*,'(A)')           '----------------------------------------------------------------------'
    write(*,'(A,'//FR_P//')')' The "eps" WENO parameter used is: ',weno_eps
    write(*,'(A,'//FI_P//')')' The "exp" WENO parameter used is: ',weno_exp
    write(*,'(A,'//FI1P//')')' The "S"   WENO parameter used is: ',S
    write(*,'(A)')           '----------------------------------------------------------------------'
    write(*,*)
  endif
  ! allocating variables
  if (allocated(weno_c)) deallocate(weno_c) ; allocate(weno_c(1:2,1:2*S))
  if (allocated(weno_a)) deallocate(weno_a) ; allocate(weno_a(1:2,0:S-1))
  if (allocated(weno_p)) deallocate(weno_p) ; allocate(weno_p(1:2,0:S-1,0:S-1))
  if (allocated(weno_d)) deallocate(weno_d) ; allocate(weno_d(0:S-1,0:S-1,0:S-1))
  ! inizializing the coefficients
  select case(S)
  case(2_I1P) ! 3rd order WENO reconstruction
    ! central difference coefficients
    ! 1 => left interface (i-1/2)
    weno_c(1,1) = -1._R_P/12._R_P ! cell -2
    weno_c(1,2) =  7._R_P/12._R_P ! cell -1
    weno_c(1,3) =  7._R_P/12._R_P ! cell  0
    weno_c(1,4) = -1._R_P/12._R_P ! cell  1
    ! 2 => right interface (i+1/2)
    weno_c(2,1) = -1._R_P/12._R_P ! cell -1
    weno_c(2,2) =  7._R_P/12._R_P ! cell  0
    weno_c(2,3) =  7._R_P/12._R_P ! cell  1
    weno_c(2,4) = -1._R_P/12._R_P ! cell  2

    ! optimal weights
    ! 1 => left interface (i-1/2)
    weno_a(1,0) = 2._R_P/3._R_P ! stencil 0
    weno_a(1,1) = 1._R_P/3._R_P ! stencil 1
    ! 2 => right interface (i+1/2)
    weno_a(2,0) = 1._R_P/3._R_P ! stencil 0
    weno_a(2,1) = 2._R_P/3._R_P ! stencil 1

    ! polinomials coefficients
    ! 1 => left interface (i-1/2)
    !  cell  0               ;    cell  1
    weno_p(1,0,0) =  0.5_R_P ; weno_p(1,1,0) =  0.5_R_P ! stencil 0
    weno_p(1,0,1) = -0.5_R_P ; weno_p(1,1,1) =  1.5_R_P ! stencil 1
    ! 2 => right interface (i+1/2)
    !  cell  0               ;    cell  1
    weno_p(2,0,0) =  1.5_R_P ; weno_p(2,1,0) = -0.5_R_P ! stencil 0
    weno_p(2,0,1) =  0.5_R_P ; weno_p(2,1,1) =  0.5_R_P ! stencil 1

    ! smoothness indicators coefficients
    ! stencil 0
    !      i*i             ;       (i-1)*i
    weno_d(0,0,0) = 1._R_P ; weno_d(1,0,0) =-2._R_P
    !      /               ;       (i-1)*(i-1)
    weno_d(0,1,0) = 0._R_P ; weno_d(1,1,0) = 1._R_P
    ! stencil 1
    !     (i+1)*(i+1)      ;       (i+1)*i
    weno_d(0,0,1) = 1._R_P ; weno_d(1,0,1) =-2._R_P
    !      /               ;        i*i
    weno_d(0,1,1) = 0._R_P ; weno_d(1,1,1) = 1._R_P
  case(3_I1P) ! 5th order WENO reconstruction
    ! central difference coefficients
    ! 1 => left interface (i-1/2)
    weno_c(1,1) =  1._R_P/60._R_P ! cell -3
    weno_c(1,2) = -7.5_R_P        ! cell -2
    weno_c(1,3) = 37._R_P/60._R_P ! cell -1
    weno_c(1,4) = 37._R_P/60._R_P ! cell  0
    weno_c(1,5) = -7.5_R_P        ! cell  1
    weno_c(1,6) =  1._R_P/60._R_P ! cell  2
    ! 2 => right interface (i+1/2)
    weno_c(1,1) =  1._R_P/60._R_P ! cell -2
    weno_c(1,2) = -7.5_R_P        ! cell -1
    weno_c(1,3) = 37._R_P/60._R_P ! cell  0
    weno_c(1,4) = 37._R_P/60._R_P ! cell  1
    weno_c(1,5) = -7.5_R_P        ! cell  2
    weno_c(1,6) =  1._R_P/60._R_P ! cell  3

    ! optimal weights
    ! 1 => left interface (i-1/2)
    weno_a(1,0) = 0.3_R_P ! stencil 0
    weno_a(1,1) = 0.6_R_P ! stencil 1
    weno_a(1,2) = 0.1_R_P ! stencil 2
    ! 2 => right interface (i+1/2)
    weno_a(2,0) = 0.1_R_P ! stencil 0
    weno_a(2,1) = 0.6_R_P ! stencil 1
    weno_a(2,2) = 0.3_R_P ! stencil 2

    ! polinomials coefficients
    ! 1 => left interface (i-1/2)
    !  cell  0                     ;    cell  1                     ;    cell  2
    weno_p(1,0,0) =  1._R_P/3._R_P ; weno_p(1,1,0) =  5._R_P/6._R_P ; weno_p(1,2,0) = -1._R_P/6._R_P ! stencil 0
    weno_p(1,0,1) = -1._R_P/6._R_P ; weno_p(1,1,1) =  5._R_P/6._R_P ; weno_p(1,2,1) =  1._R_P/3._R_P ! stencil 1
    weno_p(1,0,2) =  1._R_P/3._R_P ; weno_p(1,1,2) = -7._R_P/6._R_P ; weno_p(1,2,2) = 11._R_P/6._R_P ! stencil 2
    ! 2 => right interface (i+1/2)
    !  cell  0                     ;    cell  1                     ;    cell  2
    weno_p(2,0,0) = 11._R_P/6._R_P ; weno_p(2,1,0) = -7._R_P/6._R_P ; weno_p(2,2,0) =  1._R_P/3._R_P ! stencil 0
    weno_p(2,0,1) =  1._R_P/3._R_P ; weno_p(2,1,1) =  5._R_P/6._R_P ; weno_p(2,2,1) = -1._R_P/6._R_P ! stencil 1
    weno_p(2,0,2) = -1._R_P/6._R_P ; weno_p(2,1,2) =  5._R_P/6._R_P ; weno_p(2,2,2) =  1._R_P/3._R_P ! stencil 2

    ! smoothness indicators coefficients
    ! stencil 0
    !      i*i                      ;       (i-1)*i                   ;       (i-2)*i
    weno_d(0,0,0) =  10._R_P/3._R_P ; weno_d(1,0,0) = -31._R_P/3._R_P ; weno_d(2,0,0) =  11._R_P/3._R_P
    !      /                        ;       (i-1)*(i-1)               ;       (i-2)*(i-1)
    weno_d(0,1,0) =   0._R_P        ; weno_d(1,1,0) =  25._R_P/3._R_P ; weno_d(2,1,0) = -19._R_P/3._R_P
    !      /                        ;        /                        ;       (i-2)*(i-2)
    weno_d(0,2,0) =   0._R_P        ; weno_d(1,2,0) =   0._R_P        ; weno_d(2,2,0) =   4._R_P/3._R_P
    ! stencil 1
    !     (i+1)*(i+1)               ;        i*(i+1)                  ;       (i-1)*(i+1)
    weno_d(0,0,1) =   4._R_P/3._R_P ; weno_d(1,0,1) = -13._R_P/3._R_P ; weno_d(2,0,1) =   5._R_P/3._R_P
    !      /                        ;        i*i                      ;       (i-1)*i
    weno_d(0,1,1) =   0._R_P        ; weno_d(1,1,1) =  13._R_P/3._R_P ; weno_d(2,1,1) = -13._R_P/3._R_P
    !      /                        ;        /                        ;       (i-1)*(i-1)
    weno_d(0,2,1) =   0._R_P        ; weno_d(1,2,1) =   0._R_P        ; weno_d(2,2,1) =   4._R_P/3._R_P
    ! stencil 2
    !     (i+2)*(i+2)               ;       (i+1)*(i+2)               ;        i*(i+2)
    weno_d(0,0,2) =   4._R_P/3._R_P ; weno_d(1,0,2) = -19._R_P/3._R_P ; weno_d(2,0,2) =  11._R_P/3._R_P
    !      /                        ;       (i+1)*(i+1)               ;        i*(i+1)
    weno_d(0,1,2) =   0._R_P        ; weno_d(1,1,2) =  25._R_P/3._R_P ; weno_d(2,1,2) = -31._R_P/3._R_P
    !      /                        ;        /                        ;        i*i
    weno_d(0,2,2) =   0._R_P        ; weno_d(1,2,2) =   0._R_P        ; weno_d(2,2,2) =  10._R_P/3._R_P
  case(4_I1P) ! 7th order WENO reconstruction
    ! optimal weights
    ! 1 => left interface (i-1/2)
    weno_a(1,0) =  4._R_P/35._R_P ! stencil 0
    weno_a(1,1) = 18._R_P/35._R_P ! stencil 1
    weno_a(1,2) = 12._R_P/35._R_P ! stencil 2
    weno_a(1,3) =  1._R_P/35._R_P ! stencil 3
    ! 2 => right interface (i+1/2)
    weno_a(2,0) =  1._R_P/35._R_P ! stencil 0
    weno_a(2,1) = 12._R_P/35._R_P ! stencil 1
    weno_a(2,2) = 18._R_P/35._R_P ! stencil 2
    weno_a(2,3) =  4._R_P/35._R_P ! stencil 3

    ! polinomials coefficients
    ! 1 => left interface (i-1/2)
    !  cell  0                   ;   cell  1                   ;   cell  2                    ;   cell  3
    weno_p(1,0,0)= 1._R_P/4._R_P ;weno_p(1,1,0)=13._R_P/12._R_P;weno_p(1,2,0)= -5._R_P/12._R_P;weno_p(1,3,0)= 1._R_P/12._R_P! sten 0
    weno_p(1,0,1)=-1._R_P/12._R_P;weno_p(1,1,1)= 7._R_P/12._R_P;weno_p(1,2,1)=  7._R_P/12._R_P;weno_p(1,3,1)=-1._R_P/12._R_P! sten 1
    weno_p(1,0,2)= 1._R_P/12._R_P;weno_p(1,1,2)=-5._R_P/12._R_P;weno_p(1,2,2)= 13._R_P/12._R_P;weno_p(1,3,2)= 1._R_P/4._R_P ! sten 2
    weno_p(1,0,3)=-1._R_P/4._R_P ;weno_p(1,1,3)=13._R_P/12._R_P;weno_p(1,2,3)=-23._R_P/12._R_P;weno_p(1,3,3)=25._R_P/12._R_P! sten 3
    ! 2 => right interface (i+1/2)
    !  cell  0                   ;   cell  1                    ;   cell  2                   ;   cell  3
    weno_p(2,0,0)=25._R_P/12._R_P;weno_p(2,1,0)=-23._R_P/12._R_P;weno_p(2,2,0)=13._R_P/12._R_P;weno_p(2,3,0)=-1._R_P/4._R_P ! sten 0
    weno_p(2,0,1)= 1._R_P/4._R_P ;weno_p(2,1,1)= 13._R_P/12._R_P;weno_p(2,2,1)=-5._R_P/12._R_P;weno_p(2,3,1)= 1._R_P/12._R_P! sten 1
    weno_p(2,0,2)=-1._R_P/12._R_P;weno_p(2,1,2)=  7._R_P/12._R_P;weno_p(2,2,2)= 7._R_P/12._R_P;weno_p(2,3,2)=-1._R_P/12._R_P! sten 2
    weno_p(2,0,3)= 1._R_P/12._R_P;weno_p(2,1,3)= -5._R_P/12._R_P;weno_p(2,2,3)=13._R_P/12._R_P;weno_p(2,3,3)= 1._R_P/4._R_P ! sten 3

    ! smoothness indicators coefficients
    ! stencil 0
    !      i*i                ;       (i-1)*i             ;       (i-2)*i              ;       (i-3)*i
    weno_d(0,0,0) = 2107._R_P ; weno_d(1,0,0) =-9402._R_P ; weno_d(2,0,0) = 7042._R_P  ; weno_d(3,0,0) = -1854._R_P
    !      /                  ;       (i-1)*(i-1)         ;       (i-2)*(i-1)          ;       (i-3)*(i-1)
    weno_d(0,1,0) =   0._R_P  ; weno_d(1,1,0) =11003._R_P ; weno_d(2,1,0) =-17246._R_P ; weno_d(3,1,0) =  4642._R_P
    !      /                  ;        /                  ;       (i-2)*(i-2)          ;       (i-3)*(i-2)
    weno_d(0,2,0) =   0._R_P  ; weno_d(1,2,0) =   0._R_P  ; weno_d(2,2,0) = 7043._R_P  ; weno_d(3,2,0) = -3882._R_P
    !      /                  ;        /                  ;        /                   ;       (i-3)*(i-3)
    weno_d(0,3,0) =   0._R_P  ; weno_d(1,3,0) =   0._R_P  ; weno_d(2,3,0) =   0._R_P   ; weno_d(3,3,0) = 547._R_P
    ! stencil 1
    !     (i+1)*(i+1)         ;        i*(i+1)            ;       (i-1)*(i+1)          ;       (i-2)*(i+1)
    weno_d(0,0,1) =  547._R_P ; weno_d(1,0,1) =-2522._R_P ; weno_d(2,0,1) = 1922._R_P  ; weno_d(3,0,1) = -494._R_P
    !      /                  ;        i*i                ;       (i-1)*i              ;       (i-2)*i
    weno_d(0,1,1) =   0._R_P  ; weno_d(1,1,1) = 3443._R_P ; weno_d(2,1,1) = -5966._R_P ; weno_d(3,1,1) =  1602._R_P
    !      /                  ;        /                  ;       (i-1)*(i-1)          ;       (i-2)*(i-1)
    weno_d(0,2,1) =   0._R_P  ; weno_d(1,2,1) =   0._R_P  ; weno_d(2,2,1) = 2843._R_P  ; weno_d(3,2,1) = -1642._R_P
    !      /                  ;        /                  ;        /                   ;       (i-2)*(i-2)
    weno_d(0,3,1) =   0._R_P  ; weno_d(1,3,1) =   0._R_P  ; weno_d(2,3,1) =   0._R_P   ; weno_d(3,3,1) = 267._R_P
    ! stencil 2
    !     (i+2)*(i+2)         ;       (i+1)*(i+2)         ;        i*(i+2)             ;       (i-1)*(i+2)
    weno_d(0,0,2) =  267._R_P ; weno_d(1,0,2) =-1642._R_P ; weno_d(2,0,2) = 1602._R_P  ; weno_d(3,0,2) = -494._R_P
    !      /                  ;       (i+1)*(i+1)         ;        i*(i+1)             ;       (i-1)*(i+1)
    weno_d(0,1,2) =   0._R_P  ; weno_d(1,1,2) = 2843._R_P ; weno_d(2,1,2) = -5966._R_P ; weno_d(3,1,2) =  1922._R_P
    !      /                  ;        /                  ;        i*i                 ;       (i-1)*i
    weno_d(0,2,2) =   0._R_P  ; weno_d(1,2,2) =   0._R_P  ; weno_d(2,2,2) = 3443._R_P  ; weno_d(3,2,2) = -2522._R_P
    !      /                  ;        /                  ;        /                   ;       (i-1)*(i-1)
    weno_d(0,3,2) =   0._R_P  ; weno_d(1,3,2) =   0._R_P  ; weno_d(2,3,2) =   0._R_P   ; weno_d(3,3,2) = 547._R_P
    ! stencil 3
    !     (i+3)*(i+3)         ;       (i+2)*(i+3)         ;       (i+1)*(i+3)          ;        i*(i+3)
    weno_d(0,0,3) =  547._R_P ; weno_d(1,0,3) =-3882._R_P ; weno_d(2,0,3) = 4642._R_P  ; weno_d(3,0,3) = -1854._R_P
    !      /                  ;       (i+2)*(i+2)         ;       (i+1)*(i+2)          ;        i*(i+2)
    weno_d(0,1,3) =   0._R_P  ; weno_d(1,1,3) = 7043._R_P ; weno_d(2,1,3) =-17246._R_P ; weno_d(3,1,3) =  7042._R_P
    !      /                  ;        /                  ;       (i+1)*(i+1)          ;        i*(i+1)
    weno_d(0,2,3) =   0._R_P  ; weno_d(1,2,3) =   0._R_P  ; weno_d(2,2,3) =11003._R_P  ; weno_d(3,2,3) = -9402._R_P
    !      /                  ;        /                  ;        /                   ;        i*i
    weno_d(0,3,3) =   0._R_P  ; weno_d(1,3,3) =   0._R_P  ; weno_d(2,3,3) =   0._R_P   ; weno_d(3,3,3) = 2107._R_P
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    pure subroutine get_min_sstep(block,ss)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Function for evaluating minimum space step.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    class(Type_SBlock), intent(IN)::    block ! Block-level data.
    real(R_P),          intent(INOUT):: ss    ! Minimum space step.
    integer(I_P)::                      i,j,k ! Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    do k=1,block%Nk
      do j=1,block%Nj
        do i=1,block%Ni
          ss = min(ss,                                                                 &
                   (0.25_R_P*(block%node(i  ,j  ,k  )%x + block%node(i  ,j-1,k  )%x +  &
                              block%node(i  ,j-1,k-1)%x + block%node(i  ,j  ,k-1)%x)-  &
                    0.25_R_P*(block%node(i-1,j  ,k  )%x + block%node(i-1,j-1,k  )%x +  &
                              block%node(i-1,j-1,k-1)%x + block%node(i-1,j  ,k-1)%x)), &
                   (0.25_R_P*(block%node(i  ,j  ,k  )%y + block%node(i-1,j  ,k  )%y +  &
                              block%node(i-1,j  ,k-1)%y + block%node(i  ,j  ,k-1)%y)-  &
                    0.25_R_P*(block%node(i  ,j-1,k  )%y + block%node(i-1,j-1,k  )%y +  &
                              block%node(i-1,j-1,k-1)%y + block%node(i  ,j-1,k-1)%y)), &
                   (0.25_R_P*(block%node(i  ,j  ,k  )%z + block%node(i-1,j-1,k  )%z +  &
                              block%node(i-1,j  ,k  )%z + block%node(i  ,j-1,k  )%z)-  &
                    0.25_R_P*(block%node(i  ,j  ,k-1)%z + block%node(i-1,j-1,k-1)%z +  &
                              block%node(i-1,j  ,k-1)%z + block%node(i  ,j-1,k-1)%z)))

        enddo
      enddo
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine get_min_sstep
  endsubroutine weno_init

  !> @brief Subroutine for computing central difference reconstruction of \f$2S^{th}\f$ order instead of WENO \f$2S^{th}-1\f$ one.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b V and \b VR are declared as
  !> assumed-shape with only the lower bound defined. Their extentions are: V [1:2,0:2*S], VR [1:2].
  pure subroutine noweno_central(S,V,VR)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  S         !< Number of stencils used.
  real(R_P),    intent(IN)::  V (1:,0:) !< Variable to be reconstructed.                           [1:2,0:2*S].
  real(R_P),    intent(OUT):: VR(1:)    !< Left and right (1,2) interface value of reconstructed V [1:2      ].
  integer(I_P)::              k,f       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing the central difference reconstruction
  VR = 0._R_P
  do k=1,2*S
    do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      VR(f) = VR(f) + weno_c(f,k)*V(f,k+f-2)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine noweno_central

  !> @brief Subroutine for computing WENO reconstruction with optimal weights (without smoothness indicators computations) of
  !> \f$2S^{th}-1\f$ order.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b V and \b VR are declared as
  !> assumed-shape with only the lower bound defined. Their extentions are: V [1:2,1-S:-1+S], VR [1:2].
  pure subroutine weno_optimal(S,V,VR)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  S             !< Number of stencils used.
  real(R_P),    intent(IN)::  V (1:,1-S:)   !< Variable to be reconstructed                            [1:2,1-S:-1+S].
  real(R_P),    intent(OUT):: VR(1:)        !< Left and right (1,2) interface value of reconstructed V [1:2         ].
  real(R_P)::                 VP(1:2,0:S-1) !< Polynomial reconstructions.
  real(R_P)::                 w (1:2,0:S-1) !< Weights of the stencils.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing the polynomials
  call weno_polynomials( S = S, V = V, VP = VP)
  ! using optimal weights
  w(1:2,0:S-1) = weno_a(1:2,0:S-1)
  ! computing the convultion of reconstructing plynomials
  call weno_convolution(S = S, VP = VP, w = w, VR = VR)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine weno_optimal

  !> @brief Subroutine for computing WENO reconstruction of \f$2S^{th}-1\f$ order.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b V and \b VR are declared as
  !> assumed-shape with only the lower bound defined. Their extentions are: V [1:2,1-S:-1+S], VR [1:2].
  pure subroutine weno(S,V,VR)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  S             !< Number of stencils used.
  real(R_P),    intent(IN)::  V (1:,1-S:)   !< Variable to be reconstructed                            [1:2,1-S:-1+S].
  real(R_P),    intent(OUT):: VR(1:)        !< Left and right (1,2) interface value of reconstructed V [1:2         ].
  real(R_P)::                 VP(1:2,0:S-1) !< Polynomial reconstructions.
  real(R_P)::                 w (1:2,0:S-1) !< Weights of the stencils.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing the polynomials
  call weno_polynomials( S = S, V = V, VP = VP)
  ! computing the weights associated to the polynomials
  call weno_weights(S = S, V = V, w = w)
  ! computing the convultion of reconstructing plynomials
  call weno_convolution(S = S, VP = VP, w = w, VR = VR)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine weno
  !> @}

  !> @ingroup Lib_WENOPrivateProcedure
  !> @{
  !> @brief Subroutine for computing WENO polynomials.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b V and \b VP are declared as
  !> assumed-shape with only the lower bound defined. Their extentions are: V [1:2,1-S:-1+S], VP [1:2,0:S-1].
  pure subroutine weno_polynomials(S,V,VP)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  S           !< Number of stencils used.
  real(R_P),    intent(IN)::  V (1:,1-S:) !< Variable to be reconstructed [1:2,1-S:-1+S].
  real(R_P),    intent(OUT):: VP(1:,0:)   !< Polynomial reconstructions   [1:2,  0:S-1 ].
  integer(I_P)::              s1,s2,f     !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing the polynomials
  VP = 0._R_P
  do s1=0,S-1 ! stencil counter
    do s2=0,S-1 ! cell counter counter
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        VP(f,s1) = VP(f,s1) + weno_p(f,s2,s1)*V(f,-s2+s1)
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine weno_polynomials

  !> @brief Subroutine for computing WENO weights of the polynomial reconstructions.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b V and \b w are declared as
  !> assumed-shape with only the lower bound defined. Their extentions are: V [1:2,1-S:-1+S], w [1:2,0:S-1].
  pure subroutine weno_weights(S,V,w)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  S                !< Number of stencils used.
  real(R_P),    intent(IN)::  V    (1:,1-S:)   !< Variable to be reconstructed [1:2,1-S:-1+S].
  real(R_P),    intent(OUT):: w    (1:,0:)     !< Weights of the stencils      [1:2,  0:S-1 ].
  real(R_P)::                 IS   (1:2,0:S-1) !< Smoothness indicators of the stencils.
  real(R_P)::                 a    (1:2,0:S-1) !< Alpha coifficients for the weights.
  real(R_P)::                 a_tot(1:2)       !< Summ of the alpha coefficients.
#ifdef WENOZ
  real(R_P)::                 tao  (1:2)       !< Normalization factor.
#endif
  integer(I_P)::              s1,s2,s3,f       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing smoothness indicators
  do s1=0,S-1 ! stencil counter
    do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      IS(f,s1) = 0._R_P
      do s2=0,S-1
        do s3=0,S-1
          IS(f,s1) = IS(f,s1) + weno_d(s3,s2,s1)*V(f,s1-s3)*V(f,s1-s2)
        enddo
      enddo
    enddo
  enddo
#ifdef WENOZ
  ! computing smooth-difference tao
  if (S>2) then
    do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      tao(f) = abs(IS(f,0) - (1-weno_odd)*IS(f,1) - (1-weno_odd)*IS(f,S-2) + (1-2*weno_odd)*IS(f,S-1))
    enddo
  else
    do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      tao(f) = udiff(p=2*S-2,v=V(f,1-S:S-1))
      tao(f) = tao(f)*tao(f)
    enddo
  endif
#endif
  ! computing alfa coefficients
  a_tot = 0._R_P
  do s1=0,S-1
    do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
#ifdef WENOZ
      if (tao(f)>minval(IS(f,:))) then
        a(f,s1) = weno_a(f,s1)*(1._R_P+(tao(f)/(weno_eps+IS(f,s1)))**(weno_exp)) ; a_tot(f) = a_tot(f) + a(f,s1)
      else
        a(f,s1) = weno_a(f,s1)                                                   ; a_tot(f) = a_tot(f) + a(f,s1)
      endif
#else
      a(f,s1) = weno_a(f,s1)*(1._R_P/(weno_eps+IS(f,s1))**S)                     ; a_tot(f) = a_tot(f) + a(f,s1)
#endif
    enddo
  enddo
  ! computing the weights
  do s1=0,S-1
    do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      w(f,s1) = a(f,s1)/a_tot(f)
    enddo
  enddo
#ifdef WENOM
  ! Henrick mapping
  a_tot = 0._R_P
  do s1=0,S-1
    do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      a(f,s1)  = ( w(f,s1)*( weno_a(f,s1) + weno_a(f,s1)*weno_a(f,s1) - 3._R_P*weno_a(f,s1)*w(f,s1) + w(f,s1)*w(f,s1) ) )/ &
                 ( weno_a(f,s1)*weno_a(f,s1) + w(f,s1)*( 1._R_P - 2._R_P*weno_a(f,s1) ) )
      a_tot(f) = a_tot(f) + a(f,s1)
    enddo
  enddo
  ! computing the weights with the mapped a(f,s1)
  do s1=0,S-1
    do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      w(f,s1) = a(f,s1)/a_tot(f)
    enddo
  enddo
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine weno_weights

  !> @brief Function for computing the WENO convulution of the polynomial recontructions.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b VP \b w and \b VR are declared as
  !> assumed-shape with only the lower bound defined. Their extentions are: VP [1:2,0:S-1], w [1:2,0:S-1], V [1:2].
  pure subroutine weno_convolution(S,VP,w,VR)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  S         !< Number of stencils used.
  real(R_P),    intent(IN)::  VP(1:,0:) !< Polynomial reconstructions                              [1:2,0:S-1].
  real(R_P),    intent(IN)::  w (1:,0:) !< Weights of the stencils                                 [1:2,0:S-1].
  real(R_P),    intent(OUT):: VR(1:   ) !< Left and right (1,2) interface value of reconstructed V [1:2      ].
  integer(I_P)::              k,f       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing the convultion
  VR = 0._R_P
  do k=0,S-1
    do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      VR(f) = VR(f) + w(f,k)*VP(f,k)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine weno_convolution

  recursive pure function udiff(p,v) result(diff)
  !-------------------------------------------------------------------------------------------------------------------------------
  ! Function for computing finite undivided difference of order p of nodal function v.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: p        ! Order of finite undivided difference.
  real(R_P),    intent(IN):: v(1:p+1) ! Nodal values of function.
  real(R_P)::                diff     ! Finite undivided difference.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  if (p==1_I_P) then
    diff = v(2) - v(1)
  else
    diff = udiff(p=p-1,v=v(2:p)) - udiff(p=p-1,v=v(1:p-1))
  endif
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endfunction udiff
  !> @}
endmodule Lib_WENO
