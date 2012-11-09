!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_MeshPublicProcedure Lib_Mesh
!> @}

!> This module contains mesh procedures.
!> @ingroup Library
module Lib_Mesh
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision       ! Integers and reals precision definition.
USE Data_Type_AMRBlock ! Definition of Type_AMRBlock.
USE Data_Type_BC       ! Definition of Type_BC.
USE Data_Type_SBlock   ! Definition of Type_SBlock.
USE Data_Type_Vector   ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: node2center
public:: mesh_metrics
public:: mesh_metrics_correction
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> Subroutine for computing cell center coordinates from cell nodes ones.
  !> @ingroup Lib_MeshPublicProcedure
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
        block%cent(i,j,k)%x = (block%node(i,  j,  k  )%x + &
                               block%node(i-1,j,  k  )%x + &
                               block%node(i  ,j-1,k  )%x + &
                               block%node(i  ,j  ,k-1)%x + &
                               block%node(i-1,j-1,k-1)%x + &
                               block%node(i  ,j-1,k-1)%x + &
                               block%node(i-1,j  ,k-1)%x + &
                               block%node(i-1,j-1,k  )%x)*0.125_R_P
        block%cent(i,j,k)%y = (block%node(i,  j,  k  )%y + &
                               block%node(i-1,j,  k  )%y + &
                               block%node(i  ,j-1,k  )%y + &
                               block%node(i  ,j  ,k-1)%y + &
                               block%node(i-1,j-1,k-1)%y + &
                               block%node(i  ,j-1,k-1)%y + &
                               block%node(i-1,j  ,k-1)%y + &
                               block%node(i-1,j-1,k  )%y)*0.125_R_P
        block%cent(i,j,k)%z = (block%node(i,  j,  k  )%z + &
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

  !> Subroutine for computing the metrics of blocks.
  !> @ingroup Lib_MeshPublicProcedure
  subroutine mesh_metrics(block)
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
        block%NFi(i,j,k) = normalize(NFS)
        block%Si (i,j,k) =    normL2(NFS)
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
        block%NFj(i,j,k) = normalize(NFS)
        block%Sj (i,j,k) =    normL2(NFS)
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
        block%NFk(i,j,k) = normalize(NFS)
        block%Sk (i,j,k) =    normL2(NFS)
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

        Vx = Vx + xp*block%NFi(i,j,k)%x*block%Si(i,j,k) - xm*block%NFi(i-1,j,k)%x*block%Si(i-1,j,k)
        Vy = Vy + yp*block%NFi(i,j,k)%y*block%Si(i,j,k) - ym*block%NFi(i-1,j,k)%y*block%Si(i-1,j,k)
        Vz = Vz + zp*block%NFi(i,j,k)%z*block%Si(i,j,k) - zm*block%NFi(i-1,j,k)%z*block%Si(i-1,j,k)

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

        Vx = Vx + xp*block%NFj(i,j,k)%x*block%Sj(i,j,k) - xm*block%NFj(i,j-1,k)%x*block%Sj(i,j-1,k)
        Vy = Vy + yp*block%NFj(i,j,k)%y*block%Sj(i,j,k) - ym*block%NFj(i,j-1,k)%y*block%Sj(i,j-1,k)
        Vz = Vz + zp*block%NFj(i,j,k)%z*block%Sj(i,j,k) - zm*block%NFj(i,j-1,k)%z*block%Sj(i,j-1,k)

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

        Vx = Vx + xp*block%NFk(i,j,k)%x*block%Sk(i,j,k) - xm*block%NFk(i,j,k-1)%x*block%Sk(i,j,k-1)
        Vy = Vy + yp*block%NFk(i,j,k)%y*block%Sk(i,j,k) - ym*block%NFk(i,j,k-1)%y*block%Sk(i,j,k-1)
        Vz = Vz + zp*block%NFk(i,j,k)%z*block%Sk(i,j,k) - zm*block%NFk(i,j,k-1)%z*block%Sk(i,j,k-1)

        block%V(i,j,k) = max(Vx,Vy,Vz)
      enddo
    enddo
  enddo
  !!$OMP END PARALLEL
#ifdef NULi
  block%NFj(:,:,:)%x=0._R_P
  block%NFk(:,:,:)%x=0._R_P
  block%NFi(:,:,:)%y=0._R_P
  block%NFi(:,:,:)%z=0._R_P
#endif
#ifdef NULj
  block%NFi(:,:,:)%y=0._R_P
  block%NFk(:,:,:)%y=0._R_P
  block%NFj(:,:,:)%x=0._R_P
  block%NFj(:,:,:)%z=0._R_P
#endif
#ifdef NULk
  block%NFi(:,:,:)%z=0._R_P
  block%NFj(:,:,:)%z=0._R_P
  block%NFk(:,:,:)%x=0._R_P
  block%NFk(:,:,:)%y=0._R_P
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine mesh_metrics

  !> Subroutine for correcting the metrics of natural (and negative volume) boundary conditions cells.
  !> @ingroup Lib_MeshPublicProcedure
  subroutine mesh_metrics_correction(block)
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
      correct = ((block%BCi(0,j,k)%tp==bc_ext).OR. &
                 (block%BCi(0,j,k)%tp==bc_ref).OR. &
                 (block%BCi(0,j,k)%tp==bc_in1).OR. &
                 (block%BCi(0,j,k)%tp==bc_in2).OR. &
                 (block%V(  0,j,k)<(0.2_R_P*block%V(1,j,k))))
      wall    = ((block%BCi(0,j,k)%tp==bc_ext).OR.(block%BCi(0,j,k)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%NFi(1,j,k)*block%Si(1,j,k)).dot.block%NFi(0,j,k))
         block%NFi(-1,j,  k  ) = -(block%NFi(1,j,k)*block%Si(1,j,k)) + sn*block%NFi(0,j,k)
         block%Si (-1,j,  k  ) = normL2(   block%NFi(-1,j,k))
         block%NFi(-1,j,  k  ) = normalize(block%NFi(-1,j,k))
         ! tangential metrics
         block%NFj( 0,j  ,k  ) = tm*block%NFj(1,j  ,k  )
         block%NFj( 0,j-1,k  ) = tm*block%NFj(1,j-1,k  )
         block%NFj( 0,j  ,k-1) = tm*block%NFj(1,j  ,k-1)
         block%NFj( 0,j-1,k-1) = tm*block%NFj(1,j-1,k-1)
         block%Sj ( 0,j  ,k  ) = tm*block%Sj (1,j  ,k  )
         block%Sj ( 0,j-1,k  ) = tm*block%Sj (1,j-1,k  )
         block%Sj ( 0,j  ,k-1) = tm*block%Sj (1,j  ,k-1)
         block%Sj ( 0,j-1,k-1) = tm*block%Sj (1,j-1,k-1)

         block%NFk( 0,j  ,k  ) = tm*block%NFk(1,j  ,k  )
         block%NFk( 0,j-1,k  ) = tm*block%NFk(1,j-1,k  )
         block%NFk( 0,j  ,k-1) = tm*block%NFk(1,j  ,k-1)
         block%NFk( 0,j-1,k-1) = tm*block%NFk(1,j-1,k-1)
         block%Sk ( 0,j  ,k  ) = tm*block%Sk (1,j  ,k  )
         block%Sk ( 0,j-1,k  ) = tm*block%Sk (1,j-1,k  )
         block%Sk ( 0,j  ,k-1) = tm*block%Sk (1,j  ,k-1)
         block%Sk ( 0,j-1,k-1) = tm*block%Sk (1,j-1,k-1)
         ! volume
         block%V(   0,j,  k  ) = block%V(     1,j,  k  )
      end if
    enddo
  enddo
  ! right i
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      correct = ((block%BCi(Ni+1,j,k)%tp==bc_ext).OR. &
                 (block%BCi(Ni+1,j,k)%tp==bc_ref).OR. &
                 (block%BCi(Ni+1,j,k)%tp==bc_in1).OR. &
                 (block%BCi(Ni+1,j,k)%tp==bc_in2).OR. &
                 (block%V(Ni+1,j,k)<(0.2_R_P*block%V(Ni,j,k))))
      wall    = ((block%BCi(Ni+1,j,k)%tp==bc_ext).OR.(block%BCi(Ni+1,j,k)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%NFi(Ni-1,j,k)*block%Si(Ni-1,j,k)).dot.block%NFi(Ni,j,k))
         block%NFi(Ni+1,j,  k  ) = -(block%NFi(Ni-1,j,k)*block%Si(Ni-1,j,k)) + sn*block%NFi(Ni,j,k)
         block%Si (Ni+1,j,  k  ) = normL2(   block%NFi(Ni+1,j,k))
         block%NFi(Ni+1,j,  k  ) = normalize(block%NFi(Ni+1,j,k))
         ! tangential metrics
         block%NFj(Ni+1,j  ,k  ) = tm*block%NFj(Ni,j  ,k  )
         block%NFj(Ni+1,j-1,k  ) = tm*block%NFj(Ni,j-1,k  )
         block%NFj(Ni+1,j  ,k-1) = tm*block%NFj(Ni,j  ,k-1)
         block%NFj(Ni+1,j-1,k-1) = tm*block%NFj(Ni,j-1,k-1)
         block%Sj (Ni+1,j  ,k  ) = tm*block%Sj (Ni,j  ,k  )
         block%Sj (Ni+1,j-1,k  ) = tm*block%Sj (Ni,j-1,k  )
         block%Sj (Ni+1,j  ,k-1) = tm*block%Sj (Ni,j  ,k-1)
         block%Sj (Ni+1,j-1,k-1) = tm*block%Sj (Ni,j-1,k-1)

         block%NFk(Ni+1,j  ,k  ) = tm*block%NFk(Ni,j  ,k  )
         block%NFk(Ni+1,j-1,k  ) = tm*block%NFk(Ni,j-1,k  )
         block%NFk(Ni+1,j  ,k-1) = tm*block%NFk(Ni,j  ,k-1)
         block%NFk(Ni+1,j-1,k-1) = tm*block%NFk(Ni,j-1,k-1)
         block%Sk (Ni+1,j  ,k  ) = tm*block%Sk (Ni,j  ,k  )
         block%Sk (Ni+1,j-1,k  ) = tm*block%Sk (Ni,j-1,k  )
         block%Sk (Ni+1,j  ,k-1) = tm*block%Sk (Ni,j  ,k-1)
         block%Sk (Ni+1,j-1,k-1) = tm*block%Sk (Ni,j-1,k-1)
         ! volume
         block%V(  Ni+1,j,  k  ) = block%V(     Ni,j,  k  )
      end if
    enddo
  enddo
  ! left j
  !!$OMP DO
  do k=1,Nk
    do i=1,Ni
      correct = ((block%BCj(i,0,k)%tp==bc_ext).OR. &
                 (block%BCj(i,0,k)%tp==bc_ref).OR. &
                 (block%BCj(i,0,k)%tp==bc_in1).OR. &
                 (block%BCj(i,0,k)%tp==bc_in2).OR. &
                 (block%V(i,0,k)<(0.2_R_P*block%V(i,1,k))))
      wall    = ((block%BCj(i,0,k)%tp==bc_ext).OR.(block%BCj(i,0,k)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%NFj(i,1,k)*block%Sj(i,1,k)).dot.block%NFj(i,0,k))
         block%NFj(i, -1,k  ) = -(block%NFj(i,1,k)*block%Sj(i,1,k)) + sn*block%NFj(i,0,k)
         block%Sj (i, -1,k  ) = normL2(   block%NFj(i,-1,k))
         block%NFj(i, -1,k  ) = normalize(block%NFj(i,-1,k))
         ! tangential metrics
         block%NFi(i  ,0,k  ) = tm*block%NFi(i  ,1,k  )
         block%NFi(i-1,0,k  ) = tm*block%NFi(i-1,1,k  )
         block%NFi(i  ,0,k-1) = tm*block%NFi(i  ,1,k-1)
         block%NFi(i-1,0,k-1) = tm*block%NFi(i-1,1,k-1)
         block%Si (i  ,0,k  ) = tm*block%Si (i  ,1,k  )
         block%Si (i-1,0,k  ) = tm*block%Si (i-1,1,k  )
         block%Si (i  ,0,k-1) = tm*block%Si (i  ,1,k-1)
         block%Si (i-1,0,k-1) = tm*block%Si (i-1,1,k-1)

         block%NFk(i  ,0,k  ) = tm*block%NFk(i  ,1,k  )
         block%NFk(i-1,0,k  ) = tm*block%NFk(i-1,1,k  )
         block%NFk(i  ,0,k-1) = tm*block%NFk(i  ,1,k-1)
         block%NFk(i-1,0,k-1) = tm*block%NFk(i-1,1,k-1)
         block%Sk (i  ,0,k  ) = tm*block%Sk (i  ,1,k  )
         block%Sk (i-1,0,k  ) = tm*block%Sk (i-1,1,k  )
         block%Sk (i  ,0,k-1) = tm*block%Sk (i  ,1,k-1)
         block%Sk (i-1,0,k-1) = tm*block%Sk (i-1,1,k-1)
         ! volume
         block%V(  i,  0,k  ) = block%V(     i,  1,k  )
      end if
    enddo
  enddo
  ! right j
  !!$OMP DO
  do k=1,Nk
    do i=1,Ni
      correct = ((block%BCj(i,Nj+1,k)%tp==bc_ext).OR. &
                 (block%BCj(i,Nj+1,k)%tp==bc_ref).OR. &
                 (block%BCj(i,Nj+1,k)%tp==bc_in1).OR. &
                 (block%BCj(i,Nj+1,k)%tp==bc_in2).OR. &
                 (block%V(i,Nj+1,k)<(0.2_R_P*block%V(i,Nj,k))))
      wall    = ((block%BCj(i,Nj+1,k)%tp==bc_ext).OR.(block%BCj(i,Nj+1,k)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%NFj(i,Nj-1,k)*block%Sj(i,Nj-1,k)).dot.block%NFj(i,Nj,k))
         block%NFj(i,Nj+1,  k  ) = -(block%NFj(i,Nj-1,k)*block%Sj(i,Nj-1,k)) + sn*block%NFj(i,Nj,k)
         block%Sj (i,Nj+1,  k  ) = normL2(   block%NFj(i,Nj+1,k))
         block%NFj(i,Nj+1,  k  ) = normalize(block%NFj(i,Nj+1,k))
         ! tangential metrics
         block%NFi(i  ,Nj+1,k  ) = tm*block%NFi(i  ,Nj,k  )
         block%NFi(i-1,Nj+1,k  ) = tm*block%NFi(i-1,Nj,k  )
         block%NFi(i  ,Nj+1,k-1) = tm*block%NFi(i  ,Nj,k-1)
         block%NFi(i-1,Nj+1,k-1) = tm*block%NFi(i-1,Nj,k-1)
         block%Si (i  ,Nj+1,k  ) = tm*block%Si (i  ,Nj,k  )
         block%Si (i-1,Nj+1,k  ) = tm*block%Si (i-1,Nj,k  )
         block%Si (i  ,Nj+1,k-1) = tm*block%Si (i  ,Nj,k-1)
         block%Si (i-1,Nj+1,k-1) = tm*block%Si (i-1,Nj,k-1)

         block%NFk(i  ,Nj+1,k  ) = tm*block%NFk(i  ,Nj,k  )
         block%NFk(i-1,Nj+1,k  ) = tm*block%NFk(i-1,Nj,k  )
         block%NFk(i  ,Nj+1,k-1) = tm*block%NFk(i  ,Nj,k-1)
         block%NFk(i-1,Nj+1,k-1) = tm*block%NFk(i-1,Nj,k-1)
         block%Sk (i  ,Nj+1,k  ) = tm*block%Sk (i  ,Nj,k  )
         block%Sk (i-1,Nj+1,k  ) = tm*block%Sk (i-1,Nj,k  )
         block%Sk (i  ,Nj+1,k-1) = tm*block%Sk (i  ,Nj,k-1)
         block%Sk (i-1,Nj+1,k-1) = tm*block%Sk (i-1,Nj,k-1)
         ! volume
         block%V(  i,  Nj+1,k  ) = block%V(     i,  Nj,k  )
      end if
    enddo
  enddo
  ! left k
  !!$OMP DO
  do j=1,Nj
    do i=1,Ni
      correct = ((block%BCk(i,j,0)%tp==bc_ext).OR. &
                 (block%BCk(i,j,0)%tp==bc_ref).OR. &
                 (block%BCk(i,j,0)%tp==bc_in1).OR. &
                 (block%BCk(i,j,0)%tp==bc_in2).OR. &
                 (block%V(i,j,0)<(0.2_R_P*block%V(i,j,1))))
      wall    = ((block%BCk(i,j,0)%tp==bc_ext).OR.(block%BCk(i,j,0)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%NFk(i,j,1)*block%Sk(i,j,1)).dot.block%NFk(i,j,0))
         block%NFk(i,  j, -1) = -(block%NFk(i,j,1)*block%Sk(i,j,1)) + sn*block%NFk(i,j,0)
         block%Sk (i,  j, -1) = normL2(   block%NFk(i,j,-1))
         block%NFk(i,  j, -1) = normalize(block%NFk(i,j,-1))
         ! tangential metrics
         block%NFi(i  ,j  ,0) = tm*block%NFi(i  ,j  ,1)
         block%NFi(i-1,j  ,0) = tm*block%NFi(i-1,j  ,1)
         block%NFi(i  ,j-1,0) = tm*block%NFi(i  ,j-1,1)
         block%NFi(i-1,j-1,0) = tm*block%NFi(i-1,j-1,1)
         block%Si (i  ,j  ,0) = tm*block%Si (i  ,j  ,1)
         block%Si (i-1,j  ,0) = tm*block%Si (i-1,j  ,1)
         block%Si (i  ,j-1,0) = tm*block%Si (i  ,j-1,1)
         block%Si (i-1,j-1,0) = tm*block%Si (i-1,j-1,1)

         block%NFj(i  ,j  ,0) = tm*block%NFj(i  ,j  ,1)
         block%NFj(i-1,j  ,0) = tm*block%NFj(i-1,j  ,1)
         block%NFj(i  ,j-1,0) = tm*block%NFj(i  ,j-1,1)
         block%NFj(i-1,j-1,0) = tm*block%NFj(i-1,j-1,1)
         block%Sj (i  ,j  ,0) = tm*block%Sj (i  ,j  ,1)
         block%Sj (i-1,j  ,0) = tm*block%Sj (i-1,j  ,1)
         block%Sj (i  ,j-1,0) = tm*block%Sj (i  ,j-1,1)
         block%Sj (i-1,j-1,0) = tm*block%Sj (i-1,j-1,1)
         ! volume
         block%V(  i,  j,  0) = block%V(     i,  j,  1)
      end if
    enddo
  enddo
  ! right k
  !!$OMP DO
  do j=1,Nj
    do i=1,Ni
      correct = ((block%BCk(i,j,Nk+1)%tp==bc_ext).OR. &
                 (block%BCk(i,j,Nk+1)%tp==bc_ref).OR. &
                 (block%BCk(i,j,Nk+1)%tp==bc_in1).OR. &
                 (block%BCk(i,j,Nk+1)%tp==bc_in2).OR. &
                 (block%V(i,j,Nk+1)<(0.2_R_P*block%V(i,j,Nk))))
      wall    = ((block%BCk(i,j,Nk+1)%tp==bc_ext).OR.(block%BCk(i,j,Nk+1)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%NFk(i,j,Nk-1)*block%Sk(i,j,Nk-1)).dot.block%NFk(i,j,Nk))
         block%NFk(i,  j,  Nk+1) = -(block%NFk(i,j,Nk-1)*block%Sk(i,j,Nk-1)) + sn*block%NFk(i,j,Nk)
         block%Sk (i,  j,  Nk+1) = normL2(   block%NFk(i,j,Nk+1))
         block%NFk(i,  j,  Nk+1) = normalize(block%NFk(i,j,Nk+1))
         ! tangential metrics
         block%NFi(i  ,j  ,Nk+1) = tm*block%NFi(i  ,j  ,Nk)
         block%NFi(i-1,j  ,Nk+1) = tm*block%NFi(i-1,j  ,Nk)
         block%NFi(i  ,j-1,Nk+1) = tm*block%NFi(i  ,j-1,Nk)
         block%NFi(i-1,j-1,Nk+1) = tm*block%NFi(i-1,j-1,Nk)
         block%Si (i  ,j  ,Nk+1) = tm*block%Si (i  ,j  ,Nk)
         block%Si (i-1,j  ,Nk+1) = tm*block%Si (i-1,j  ,Nk)
         block%Si (i  ,j-1,Nk+1) = tm*block%Si (i  ,j-1,Nk)
         block%Si (i-1,j-1,Nk+1) = tm*block%Si (i-1,j-1,Nk)

         block%NFj(i  ,j  ,Nk+1) = tm*block%NFj(i  ,j  ,Nk)
         block%NFj(i-1,j  ,Nk+1) = tm*block%NFj(i-1,j  ,Nk)
         block%NFj(i  ,j-1,Nk+1) = tm*block%NFj(i  ,j-1,Nk)
         block%NFj(i-1,j-1,Nk+1) = tm*block%NFj(i-1,j-1,Nk)
         block%Sj (i  ,j  ,Nk+1) = tm*block%Sj (i  ,j  ,Nk)
         block%Sj (i-1,j  ,Nk+1) = tm*block%Sj (i-1,j  ,Nk)
         block%Sj (i  ,j-1,Nk+1) = tm*block%Sj (i  ,j-1,Nk)
         block%Sj (i-1,j-1,Nk+1) = tm*block%Sj (i-1,j-1,Nk)
         ! volume
         block%V(  i,  j,  Nk+1) = block%V(     i,  j,  Nk)
      end if
    enddo
  enddo
  !!$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine mesh_metrics_correction
endmodule Lib_Mesh
