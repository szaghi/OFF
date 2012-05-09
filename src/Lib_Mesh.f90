!> This module contains mesh procedures.
!> @ingroup Library
module Lib_Mesh
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision      ! Integers and reals precision definition.
USE Data_Type_BC      ! Definition of Type_BC.
USE Data_Type_Globals ! Definition of Type_Global and Type_Block.
USE Data_Type_Vector  ! Definition of Type_Vector.
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
  subroutine node2center(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Block), intent(INOUT):: block !< Block-level data.
  integer(I_P)::                    i,j,k !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(block)
  !$OMP DO
  do k=1-block%mesh%gc(5),block%mesh%Nk+block%mesh%gc(6)
    do j=1-block%mesh%gc(3),block%mesh%Nj+block%mesh%gc(4)
      do i=1-block%mesh%gc(1),block%mesh%Ni+block%mesh%gc(2)
        block%mesh%cent(i,j,k)%x = (block%mesh%node(i,  j,  k  )%x + &
                                    block%mesh%node(i-1,j,  k  )%x + &
                                    block%mesh%node(i  ,j-1,k  )%x + &
                                    block%mesh%node(i  ,j  ,k-1)%x + &
                                    block%mesh%node(i-1,j-1,k-1)%x + &
                                    block%mesh%node(i  ,j-1,k-1)%x + &
                                    block%mesh%node(i-1,j  ,k-1)%x + &
                                    block%mesh%node(i-1,j-1,k  )%x)*0.125_R_P
        block%mesh%cent(i,j,k)%y = (block%mesh%node(i,  j,  k  )%y + &
                                    block%mesh%node(i-1,j,  k  )%y + &
                                    block%mesh%node(i  ,j-1,k  )%y + &
                                    block%mesh%node(i  ,j  ,k-1)%y + &
                                    block%mesh%node(i-1,j-1,k-1)%y + &
                                    block%mesh%node(i  ,j-1,k-1)%y + &
                                    block%mesh%node(i-1,j  ,k-1)%y + &
                                    block%mesh%node(i-1,j-1,k  )%y)*0.125_R_P
        block%mesh%cent(i,j,k)%z = (block%mesh%node(i,  j,  k  )%z + &
                                    block%mesh%node(i-1,j,  k  )%z + &
                                    block%mesh%node(i  ,j-1,k  )%z + &
                                    block%mesh%node(i  ,j  ,k-1)%z + &
                                    block%mesh%node(i-1,j-1,k-1)%z + &
                                    block%mesh%node(i  ,j-1,k-1)%z + &
                                    block%mesh%node(i-1,j  ,k-1)%z + &
                                    block%mesh%node(i-1,j-1,k  )%z)*0.125_R_P
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine node2center

  !> Subroutine for computing the metrics of blocks.
  subroutine mesh_metrics(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Block), intent(INOUT):: block             !< Block-level data.
  type(Type_Vector)::               NFS,s1,s2,nd,db   !< Dummy vector variables.
  real(R_P)::                       signi,signj,signk !< Dummy variables for checking the directions of normals.
  real(R_P)::                       Vx,Vy,Vz          !< Dummy variables for computing volume.
  real(R_P)::                       xp,yp,zp          !< Dummy variables for computing face coordinates.
  real(R_P)::                       xm,ym,zm          !< Dummy variables for computing face coordinates.
  integer(I_P)::                    i,j,k             !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing faces normals
  ! positioning at the middle of the block
  i = max(1,block%mesh%Ni)
  j = max(1,block%mesh%Nj)
  k = max(1,block%mesh%Nk)
  ! checking the direction of i normals
  s1 = block%mesh%node(i,j  ,k) - block%mesh%node(i,j-1,k-1)
  s2 = block%mesh%node(i,j-1,k) - block%mesh%node(i,j,  k-1)
  nd = s1.cross.s2
  s1 = 0.25_R_P*(block%mesh%node(i,  j,k)+block%mesh%node(i,  j-1,k)+block%mesh%node(i,  j,k-1)+block%mesh%node(i,  j-1,k-1))
  s2 = 0.25_R_P*(block%mesh%node(i-1,j,k)+block%mesh%node(i-1,j-1,k)+block%mesh%node(i-1,j,k-1)+block%mesh%node(i-1,j-1,k-1))
  db = s1 - s2
  signi = sign(1._R_P,(nd.dot.db))
  ! checking the direction of j normals
  s1 = block%mesh%node(i,j,k  ) - block%mesh%node(i-1,j,k-1)
  s2 = block%mesh%node(i,j,k-1) - block%mesh%node(i-1,j,k  )
  nd = s1.cross.s2
  s1 = 0.25_R_P*(block%mesh%node(i,j,  k)+block%mesh%node(i-1,j,  k)+block%mesh%node(i,j,  k-1)+block%mesh%node(i-1,j,  k-1))
  s2 = 0.25_R_P*(block%mesh%node(i,j-1,k)+block%mesh%node(i-1,j-1,k)+block%mesh%node(i,j-1,k-1)+block%mesh%node(i-1,j-1,k-1))
  db = s1 - s2
  signj = sign(1._R_P,(nd.dot.db))
  ! checking the direction of k normals
  s1 = block%mesh%node(i,  j,k) - block%mesh%node(i-1,j-1,k)
  s2 = block%mesh%node(i-1,j,k) - block%mesh%node(i,  j-1,k)
  nd = s1.cross.s2
  s1 = 0.25_R_P*(block%mesh%node(i,j,k  )+block%mesh%node(i-1,j,k  )+block%mesh%node(i,j-1,k  )+block%mesh%node(i-1,j-1,k  ))
  s2 = 0.25_R_P*(block%mesh%node(i,j,k-1)+block%mesh%node(i-1,j,k-1)+block%mesh%node(i,j-1,k-1)+block%mesh%node(i-1,j-1,k-1))
  db = s1 - s2
  signk = sign(1._R_P,(nd.dot.db))
  !$OMP PARALLEL DEFAULT(NONE)                        &
  !$OMP PRIVATE(i,j,k,NFS,Vx,Vy,Vz,xp,yp,zp,xm,ym,zm) &
  !$OMP SHARED(block,signi,signj,signk)
  !$OMP DO
  do k=1,block%mesh%Nk
    do j=1,block%mesh%Nj
      do i=0,block%mesh%Ni
        NFS = face_normal4(pt1 = block%mesh%node(i,j-1,k-1), &
                           pt2 = block%mesh%node(i,j  ,k-1), &
                           pt3 = block%mesh%node(i,j  ,k  ), &
                           pt4 = block%mesh%node(i,j-1,k  ))
        NFS = NFS*signi
        block%mesh%NFi(i,j,k) = normalize(NFS)
        block%mesh%Si (i,j,k) =    normL2(NFS)
      enddo
    enddo
  enddo
  !$OMP DO
  do k=1,block%mesh%Nk
    do j=0,block%mesh%Nj
      do i=1,block%mesh%Ni
        NFS = face_normal4(pt1 = block%mesh%node(i-1,j,k-1), &
                           pt2 = block%mesh%node(i-1,j,k  ), &
                           pt3 = block%mesh%node(i  ,j,k  ), &
                           pt4 = block%mesh%node(i  ,j,k-1))
        NFS = NFS*signj
        block%mesh%NFj(i,j,k) = normalize(NFS)
        block%mesh%Sj (i,j,k) =    normL2(NFS)
      enddo
    enddo
  enddo
  !$OMP DO
  do k=0,block%mesh%Nk
    do j=1,block%mesh%Nj
      do i=1,block%mesh%Ni
        NFS = face_normal4(pt1 = block%mesh%node(i-1,j-1,k), &
                           pt2 = block%mesh%node(i  ,j-1,k), &
                           pt3 = block%mesh%node(i  ,j  ,k), &
                           pt4 = block%mesh%node(i-1,j  ,k))
        NFS = NFS*signk
        block%mesh%NFk(i,j,k) = normalize(NFS)
        block%mesh%Sk (i,j,k) =    normL2(NFS)
      enddo
    enddo
  enddo
  ! computing finte volumes
  !$OMP DO
  do k=1,block%mesh%Nk
    do j=1,block%mesh%Nj
      do i=1,block%mesh%Ni
        Vx = 0._R_P
        Vy = 0._R_P
        Vz = 0._R_P

        xp = 0.25_R_P*(block%mesh%node(i  ,j  ,k  )%x + block%mesh%node(i  ,j  ,k-1)%x + &
                       block%mesh%node(i  ,j-1,k  )%x + block%mesh%node(i  ,j-1,k-1)%x)
        yp = 0.25_R_P*(block%mesh%node(i  ,j  ,k  )%y + block%mesh%node(i  ,j  ,k-1)%y + &
                       block%mesh%node(i  ,j-1,k  )%y + block%mesh%node(i  ,j-1,k-1)%y)
        zp = 0.25_R_P*(block%mesh%node(i  ,j  ,k  )%z + block%mesh%node(i  ,j  ,k-1)%z + &
                       block%mesh%node(i  ,j-1,k  )%z + block%mesh%node(i  ,j-1,k-1)%z)
        xm = 0.25_R_P*(block%mesh%node(i-1,j  ,k  )%x + block%mesh%node(i-1,j  ,k-1)%x + &
                       block%mesh%node(i-1,j-1,k  )%x + block%mesh%node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(block%mesh%node(i-1,j  ,k  )%y + block%mesh%node(i-1,j  ,k-1)%y + &
                       block%mesh%node(i-1,j-1,k  )%y + block%mesh%node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(block%mesh%node(i-1,j  ,k  )%z + block%mesh%node(i-1,j  ,k-1)%z + &
                       block%mesh%node(i-1,j-1,k  )%z + block%mesh%node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*block%mesh%NFi(i,j,k)%x*block%mesh%Si(i,j,k) - xm*block%mesh%NFi(i-1,j,k)%x*block%mesh%Si(i-1,j,k)
        Vy = Vy + yp*block%mesh%NFi(i,j,k)%y*block%mesh%Si(i,j,k) - ym*block%mesh%NFi(i-1,j,k)%y*block%mesh%Si(i-1,j,k)
        Vz = Vz + zp*block%mesh%NFi(i,j,k)%z*block%mesh%Si(i,j,k) - zm*block%mesh%NFi(i-1,j,k)%z*block%mesh%Si(i-1,j,k)

        xp = 0.25_R_P*(block%mesh%node(i  ,j  ,k  )%x + block%mesh%node(i  ,j  ,k-1)%x + &
                       block%mesh%node(i-1,j  ,k  )%x + block%mesh%node(i-1,j  ,k-1)%x)
        yp = 0.25_R_P*(block%mesh%node(i  ,j  ,k  )%y + block%mesh%node(i  ,j  ,k-1)%y + &
                       block%mesh%node(i-1,j  ,k  )%y + block%mesh%node(i-1,j  ,k-1)%y)
        zp = 0.25_R_P*(block%mesh%node(i  ,j  ,k  )%z + block%mesh%node(i  ,j  ,k-1)%z + &
                       block%mesh%node(i-1,j  ,k  )%z + block%mesh%node(i-1,j  ,k-1)%z)
        xm = 0.25_R_P*(block%mesh%node(i  ,j-1,k  )%x + block%mesh%node(i  ,j-1,k-1)%x + &
                       block%mesh%node(i-1,j-1,k  )%x + block%mesh%node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(block%mesh%node(i  ,j-1,k  )%y + block%mesh%node(i  ,j-1,k-1)%y + &
                       block%mesh%node(i-1,j-1,k  )%y + block%mesh%node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(block%mesh%node(i  ,j-1,k  )%z + block%mesh%node(i  ,j-1,k-1)%z + &
                       block%mesh%node(i-1,j-1,k  )%z + block%mesh%node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*block%mesh%NFj(i,j,k)%x*block%mesh%Sj(i,j,k) - xm*block%mesh%NFj(i,j-1,k)%x*block%mesh%Sj(i,j-1,k)
        Vy = Vy + yp*block%mesh%NFj(i,j,k)%y*block%mesh%Sj(i,j,k) - ym*block%mesh%NFj(i,j-1,k)%y*block%mesh%Sj(i,j-1,k)
        Vz = Vz + zp*block%mesh%NFj(i,j,k)%z*block%mesh%Sj(i,j,k) - zm*block%mesh%NFj(i,j-1,k)%z*block%mesh%Sj(i,j-1,k)

        xp = 0.25_R_P*(block%mesh%node(i  ,j  ,k  )%x + block%mesh%node(i  ,j-1,k  )%x + &
                       block%mesh%node(i-1,j  ,k  )%x + block%mesh%node(i-1,j-1,k  )%x)
        yp = 0.25_R_P*(block%mesh%node(i  ,j  ,k  )%y + block%mesh%node(i  ,j-1,k  )%y + &
                       block%mesh%node(i-1,j  ,k  )%y + block%mesh%node(i-1,j-1,k  )%y)
        zp = 0.25_R_P*(block%mesh%node(i  ,j  ,k  )%z + block%mesh%node(i  ,j-1,k  )%z + &
                       block%mesh%node(i-1,j  ,k  )%z + block%mesh%node(i-1,j-1,k  )%z)
        xm = 0.25_R_P*(block%mesh%node(i  ,j  ,k-1)%x + block%mesh%node(i  ,j-1,k-1)%x + &
                       block%mesh%node(i-1,j  ,k-1)%x + block%mesh%node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(block%mesh%node(i  ,j  ,k-1)%y + block%mesh%node(i  ,j-1,k-1)%y + &
                       block%mesh%node(i-1,j  ,k-1)%y + block%mesh%node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(block%mesh%node(i  ,j  ,k-1)%z + block%mesh%node(i  ,j-1,k-1)%z + &
                       block%mesh%node(i-1,j  ,k-1)%z + block%mesh%node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*block%mesh%NFk(i,j,k)%x*block%mesh%Sk(i,j,k) - xm*block%mesh%NFk(i,j,k-1)%x*block%mesh%Sk(i,j,k-1)
        Vy = Vy + yp*block%mesh%NFk(i,j,k)%y*block%mesh%Sk(i,j,k) - ym*block%mesh%NFk(i,j,k-1)%y*block%mesh%Sk(i,j,k-1)
        Vz = Vz + zp*block%mesh%NFk(i,j,k)%z*block%mesh%Sk(i,j,k) - zm*block%mesh%NFk(i,j,k-1)%z*block%mesh%Sk(i,j,k-1)

        block%mesh%V(i,j,k) = max(Vx,Vy,Vz)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
#ifdef NULi
  block%mesh%NFj(:,:,:)%x=0._R_P
  block%mesh%NFk(:,:,:)%x=0._R_P
  block%mesh%NFi(:,:,:)%y=0._R_P
  block%mesh%NFi(:,:,:)%z=0._R_P
#endif
#ifdef NULj
  block%mesh%NFi(:,:,:)%y=0._R_P
  block%mesh%NFk(:,:,:)%y=0._R_P
  block%mesh%NFj(:,:,:)%x=0._R_P
  block%mesh%NFj(:,:,:)%z=0._R_P
#endif
#ifdef NULk
  block%mesh%NFi(:,:,:)%z=0._R_P
  block%mesh%NFj(:,:,:)%z=0._R_P
  block%mesh%NFk(:,:,:)%x=0._R_P
  block%mesh%NFk(:,:,:)%y=0._R_P
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine mesh_metrics

  !> Subroutine for correcting the metrics of natural (and negative volume) boundary conditions cells.
  subroutine mesh_metrics_correction(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Block), intent(INOUT):: block    !< Block-level data.
  logical::                         correct  !< Flag for inquiring if metrics must be corrected.
  logical::                         wall     !< Flag for inquiring if bc is "wall-type": different corrections must be used.
  real(R_P)::                       tm       !< Tangential metrics parameter (-1 for wall-type bc).
  real(R_P)::                       sn       !< Normal metrics coefficient correction.
  integer(I_P)::                    Ni,Nj,Nk !< Dimensions of the block.
  integer(I_P)::                    i,j,k    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = block%mesh%Ni
  Nj = block%mesh%Nj
  Nk = block%mesh%Nk
  !$OMP PARALLEL DEFAULT(NONE)            &
  !$OMP PRIVATE(i,j,k,correct,wall,tm,sn) &
  !$OMP SHARED(Ni,Nj,Nk,block)
  ! left i
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
      correct = ((block%bc%BCi(0,j,k)%tp==bc_ext).OR. &
                 (block%bc%BCi(0,j,k)%tp==bc_ref).OR. &
                 (block%bc%BCi(0,j,k)%tp==bc_in1).OR. &
                 (block%bc%BCi(0,j,k)%tp==bc_in2).OR. &
                 (block%mesh%V(0,j,k)<(0.2_R_P*block%mesh%V(1,j,k))))
      wall    = ((block%bc%BCi(0,j,k)%tp==bc_ext).OR.(block%bc%BCi(0,j,k)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%mesh%NFi(1,j,k)*block%mesh%Si(1,j,k)).dot.block%mesh%NFi(0,j,k))
         block%mesh%NFi(-1,j,  k  ) = -(block%mesh%NFi(1,j,k)*block%mesh%Si(1,j,k)) + sn*block%mesh%NFi(0,j,k)
         block%mesh%Si (-1,j,  k  ) = normL2(   block%mesh%NFi(-1,j,k))
         block%mesh%NFi(-1,j,  k  ) = normalize(block%mesh%NFi(-1,j,k))
         ! tangential metrics
         block%mesh%NFj( 0,j  ,k  ) = tm*block%mesh%NFj(1,j  ,k  )
         block%mesh%NFj( 0,j-1,k  ) = tm*block%mesh%NFj(1,j-1,k  )
         block%mesh%NFj( 0,j  ,k-1) = tm*block%mesh%NFj(1,j  ,k-1)
         block%mesh%NFj( 0,j-1,k-1) = tm*block%mesh%NFj(1,j-1,k-1)
         block%mesh%Sj ( 0,j  ,k  ) = tm*block%mesh%Sj (1,j  ,k  )
         block%mesh%Sj ( 0,j-1,k  ) = tm*block%mesh%Sj (1,j-1,k  )
         block%mesh%Sj ( 0,j  ,k-1) = tm*block%mesh%Sj (1,j  ,k-1)
         block%mesh%Sj ( 0,j-1,k-1) = tm*block%mesh%Sj (1,j-1,k-1)

         block%mesh%NFk( 0,j  ,k  ) = tm*block%mesh%NFk(1,j  ,k  )
         block%mesh%NFk( 0,j-1,k  ) = tm*block%mesh%NFk(1,j-1,k  )
         block%mesh%NFk( 0,j  ,k-1) = tm*block%mesh%NFk(1,j  ,k-1)
         block%mesh%NFk( 0,j-1,k-1) = tm*block%mesh%NFk(1,j-1,k-1)
         block%mesh%Sk ( 0,j  ,k  ) = tm*block%mesh%Sk (1,j  ,k  )
         block%mesh%Sk ( 0,j-1,k  ) = tm*block%mesh%Sk (1,j-1,k  )
         block%mesh%Sk ( 0,j  ,k-1) = tm*block%mesh%Sk (1,j  ,k-1)
         block%mesh%Sk ( 0,j-1,k-1) = tm*block%mesh%Sk (1,j-1,k-1)
         ! volume
         block%mesh%V(   0,j,  k  ) = block%mesh%V(     1,j,  k  )
      end if
    enddo
  enddo
  ! right i
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
      correct = ((block%bc%BCi(Ni+1,j,k)%tp==bc_ext).OR. &
                 (block%bc%BCi(Ni+1,j,k)%tp==bc_ref).OR. &
                 (block%bc%BCi(Ni+1,j,k)%tp==bc_in1).OR. &
                 (block%bc%BCi(Ni+1,j,k)%tp==bc_in2).OR. &
                 (block%mesh%V(Ni+1,j,k)<(0.2_R_P*block%mesh%V(Ni,j,k))))
      wall    = ((block%bc%BCi(Ni+1,j,k)%tp==bc_ext).OR.(block%bc%BCi(Ni+1,j,k)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%mesh%NFi(Ni-1,j,k)*block%mesh%Si(Ni-1,j,k)).dot.block%mesh%NFi(Ni,j,k))
         block%mesh%NFi(Ni+1,j,  k  ) = -(block%mesh%NFi(Ni-1,j,k)*block%mesh%Si(Ni-1,j,k)) + sn*block%mesh%NFi(Ni,j,k)
         block%mesh%Si (Ni+1,j,  k  ) = normL2(   block%mesh%NFi(Ni+1,j,k))
         block%mesh%NFi(Ni+1,j,  k  ) = normalize(block%mesh%NFi(Ni+1,j,k))
         ! tangential metrics
         block%mesh%NFj(Ni+1,j  ,k  ) = tm*block%mesh%NFj(Ni,j  ,k  )
         block%mesh%NFj(Ni+1,j-1,k  ) = tm*block%mesh%NFj(Ni,j-1,k  )
         block%mesh%NFj(Ni+1,j  ,k-1) = tm*block%mesh%NFj(Ni,j  ,k-1)
         block%mesh%NFj(Ni+1,j-1,k-1) = tm*block%mesh%NFj(Ni,j-1,k-1)
         block%mesh%Sj (Ni+1,j  ,k  ) = tm*block%mesh%Sj (Ni,j  ,k  )
         block%mesh%Sj (Ni+1,j-1,k  ) = tm*block%mesh%Sj (Ni,j-1,k  )
         block%mesh%Sj (Ni+1,j  ,k-1) = tm*block%mesh%Sj (Ni,j  ,k-1)
         block%mesh%Sj (Ni+1,j-1,k-1) = tm*block%mesh%Sj (Ni,j-1,k-1)

         block%mesh%NFk(Ni+1,j  ,k  ) = tm*block%mesh%NFk(Ni,j  ,k  )
         block%mesh%NFk(Ni+1,j-1,k  ) = tm*block%mesh%NFk(Ni,j-1,k  )
         block%mesh%NFk(Ni+1,j  ,k-1) = tm*block%mesh%NFk(Ni,j  ,k-1)
         block%mesh%NFk(Ni+1,j-1,k-1) = tm*block%mesh%NFk(Ni,j-1,k-1)
         block%mesh%Sk (Ni+1,j  ,k  ) = tm*block%mesh%Sk (Ni,j  ,k  )
         block%mesh%Sk (Ni+1,j-1,k  ) = tm*block%mesh%Sk (Ni,j-1,k  )
         block%mesh%Sk (Ni+1,j  ,k-1) = tm*block%mesh%Sk (Ni,j  ,k-1)
         block%mesh%Sk (Ni+1,j-1,k-1) = tm*block%mesh%Sk (Ni,j-1,k-1)
         ! volume
         block%mesh%V(  Ni+1,j,  k  ) = block%mesh%V(     Ni,j,  k  )
      end if
    enddo
  enddo
  ! left j
  !$OMP DO
  do k=1,Nk
    do i=1,Ni
      correct = ((block%bc%BCj(i,0,k)%tp==bc_ext).OR. &
                 (block%bc%BCj(i,0,k)%tp==bc_ref).OR. &
                 (block%bc%BCj(i,0,k)%tp==bc_in1).OR. &
                 (block%bc%BCj(i,0,k)%tp==bc_in2).OR. &
                 (block%mesh%V(i,0,k)<(0.2_R_P*block%mesh%V(i,1,k))))
      wall    = ((block%bc%BCj(i,0,k)%tp==bc_ext).OR.(block%bc%BCj(i,0,k)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%mesh%NFj(i,1,k)*block%mesh%Sj(i,1,k)).dot.block%mesh%NFj(i,0,k))
         block%mesh%NFj(i, -1,k  ) = -(block%mesh%NFj(i,1,k)*block%mesh%Sj(i,1,k)) + sn*block%mesh%NFj(i,0,k)
         block%mesh%Sj (i, -1,k  ) = normL2(   block%mesh%NFj(i,-1,k))
         block%mesh%NFj(i, -1,k  ) = normalize(block%mesh%NFj(i,-1,k))
         ! tangential metrics
         block%mesh%NFi(i  ,0,k  ) = tm*block%mesh%NFi(i  ,1,k  )
         block%mesh%NFi(i-1,0,k  ) = tm*block%mesh%NFi(i-1,1,k  )
         block%mesh%NFi(i  ,0,k-1) = tm*block%mesh%NFi(i  ,1,k-1)
         block%mesh%NFi(i-1,0,k-1) = tm*block%mesh%NFi(i-1,1,k-1)
         block%mesh%Si (i  ,0,k  ) = tm*block%mesh%Si (i  ,1,k  )
         block%mesh%Si (i-1,0,k  ) = tm*block%mesh%Si (i-1,1,k  )
         block%mesh%Si (i  ,0,k-1) = tm*block%mesh%Si (i  ,1,k-1)
         block%mesh%Si (i-1,0,k-1) = tm*block%mesh%Si (i-1,1,k-1)

         block%mesh%NFk(i  ,0,k  ) = tm*block%mesh%NFk(i  ,1,k  )
         block%mesh%NFk(i-1,0,k  ) = tm*block%mesh%NFk(i-1,1,k  )
         block%mesh%NFk(i  ,0,k-1) = tm*block%mesh%NFk(i  ,1,k-1)
         block%mesh%NFk(i-1,0,k-1) = tm*block%mesh%NFk(i-1,1,k-1)
         block%mesh%Sk (i  ,0,k  ) = tm*block%mesh%Sk (i  ,1,k  )
         block%mesh%Sk (i-1,0,k  ) = tm*block%mesh%Sk (i-1,1,k  )
         block%mesh%Sk (i  ,0,k-1) = tm*block%mesh%Sk (i  ,1,k-1)
         block%mesh%Sk (i-1,0,k-1) = tm*block%mesh%Sk (i-1,1,k-1)
         ! volume
         block%mesh%V(  i,  0,k  ) = block%mesh%V(     i,  1,k  )
      end if
    enddo
  enddo
  ! right j
  !$OMP DO
  do k=1,Nk
    do i=1,Ni
      correct = ((block%bc%BCj(i,Nj+1,k)%tp==bc_ext).OR. &
                 (block%bc%BCj(i,Nj+1,k)%tp==bc_ref).OR. &
                 (block%bc%BCj(i,Nj+1,k)%tp==bc_in1).OR. &
                 (block%bc%BCj(i,Nj+1,k)%tp==bc_in2).OR. &
                 (block%mesh%V(i,Nj+1,k)<(0.2_R_P*block%mesh%V(i,Nj,k))))
      wall    = ((block%bc%BCj(i,Nj+1,k)%tp==bc_ext).OR.(block%bc%BCj(i,Nj+1,k)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%mesh%NFj(i,Nj-1,k)*block%mesh%Sj(i,Nj-1,k)).dot.block%mesh%NFj(i,Nj,k))
         block%mesh%NFj(i,Nj+1,  k  ) = -(block%mesh%NFj(i,Nj-1,k)*block%mesh%Sj(i,Nj-1,k)) + sn*block%mesh%NFj(i,Nj,k)
         block%mesh%Sj (i,Nj+1,  k  ) = normL2(   block%mesh%NFj(i,Nj+1,k))
         block%mesh%NFj(i,Nj+1,  k  ) = normalize(block%mesh%NFj(i,Nj+1,k))
         ! tangential metrics
         block%mesh%NFi(i  ,Nj+1,k  ) = tm*block%mesh%NFi(i  ,Nj,k  )
         block%mesh%NFi(i-1,Nj+1,k  ) = tm*block%mesh%NFi(i-1,Nj,k  )
         block%mesh%NFi(i  ,Nj+1,k-1) = tm*block%mesh%NFi(i  ,Nj,k-1)
         block%mesh%NFi(i-1,Nj+1,k-1) = tm*block%mesh%NFi(i-1,Nj,k-1)
         block%mesh%Si (i  ,Nj+1,k  ) = tm*block%mesh%Si (i  ,Nj,k  )
         block%mesh%Si (i-1,Nj+1,k  ) = tm*block%mesh%Si (i-1,Nj,k  )
         block%mesh%Si (i  ,Nj+1,k-1) = tm*block%mesh%Si (i  ,Nj,k-1)
         block%mesh%Si (i-1,Nj+1,k-1) = tm*block%mesh%Si (i-1,Nj,k-1)

         block%mesh%NFk(i  ,Nj+1,k  ) = tm*block%mesh%NFk(i  ,Nj,k  )
         block%mesh%NFk(i-1,Nj+1,k  ) = tm*block%mesh%NFk(i-1,Nj,k  )
         block%mesh%NFk(i  ,Nj+1,k-1) = tm*block%mesh%NFk(i  ,Nj,k-1)
         block%mesh%NFk(i-1,Nj+1,k-1) = tm*block%mesh%NFk(i-1,Nj,k-1)
         block%mesh%Sk (i  ,Nj+1,k  ) = tm*block%mesh%Sk (i  ,Nj,k  )
         block%mesh%Sk (i-1,Nj+1,k  ) = tm*block%mesh%Sk (i-1,Nj,k  )
         block%mesh%Sk (i  ,Nj+1,k-1) = tm*block%mesh%Sk (i  ,Nj,k-1)
         block%mesh%Sk (i-1,Nj+1,k-1) = tm*block%mesh%Sk (i-1,Nj,k-1)
         ! volume
         block%mesh%V(  i,  Nj+1,k  ) = block%mesh%V(     i,  Nj,k  )
      end if
    enddo
  enddo
  ! left k
  !$OMP DO
  do j=1,Nj
    do i=1,Ni
      correct = ((block%bc%BCk(i,j,0)%tp==bc_ext).OR. &
                 (block%bc%BCk(i,j,0)%tp==bc_ref).OR. &
                 (block%bc%BCk(i,j,0)%tp==bc_in1).OR. &
                 (block%bc%BCk(i,j,0)%tp==bc_in2).OR. &
                 (block%mesh%V(i,j,0)<(0.2_R_P*block%mesh%V(i,j,1))))
      wall    = ((block%bc%BCk(i,j,0)%tp==bc_ext).OR.(block%bc%BCk(i,j,0)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%mesh%NFk(i,j,1)*block%mesh%Sk(i,j,1)).dot.block%mesh%NFk(i,j,0))
         block%mesh%NFk(i,  j, -1) = -(block%mesh%NFk(i,j,1)*block%mesh%Sk(i,j,1)) + sn*block%mesh%NFk(i,j,0)
         block%mesh%Sk (i,  j, -1) = normL2(   block%mesh%NFk(i,j,-1))
         block%mesh%NFk(i,  j, -1) = normalize(block%mesh%NFk(i,j,-1))
         ! tangential metrics
         block%mesh%NFi(i  ,j  ,0) = tm*block%mesh%NFi(i  ,j  ,1)
         block%mesh%NFi(i-1,j  ,0) = tm*block%mesh%NFi(i-1,j  ,1)
         block%mesh%NFi(i  ,j-1,0) = tm*block%mesh%NFi(i  ,j-1,1)
         block%mesh%NFi(i-1,j-1,0) = tm*block%mesh%NFi(i-1,j-1,1)
         block%mesh%Si (i  ,j  ,0) = tm*block%mesh%Si (i  ,j  ,1)
         block%mesh%Si (i-1,j  ,0) = tm*block%mesh%Si (i-1,j  ,1)
         block%mesh%Si (i  ,j-1,0) = tm*block%mesh%Si (i  ,j-1,1)
         block%mesh%Si (i-1,j-1,0) = tm*block%mesh%Si (i-1,j-1,1)

         block%mesh%NFj(i  ,j  ,0) = tm*block%mesh%NFj(i  ,j  ,1)
         block%mesh%NFj(i-1,j  ,0) = tm*block%mesh%NFj(i-1,j  ,1)
         block%mesh%NFj(i  ,j-1,0) = tm*block%mesh%NFj(i  ,j-1,1)
         block%mesh%NFj(i-1,j-1,0) = tm*block%mesh%NFj(i-1,j-1,1)
         block%mesh%Sj (i  ,j  ,0) = tm*block%mesh%Sj (i  ,j  ,1)
         block%mesh%Sj (i-1,j  ,0) = tm*block%mesh%Sj (i-1,j  ,1)
         block%mesh%Sj (i  ,j-1,0) = tm*block%mesh%Sj (i  ,j-1,1)
         block%mesh%Sj (i-1,j-1,0) = tm*block%mesh%Sj (i-1,j-1,1)
         ! volume
         block%mesh%V(  i,  j,  0) = block%mesh%V(     i,  j,  1)
      end if
    enddo
  enddo
  ! right k
  !$OMP DO
  do j=1,Nj
    do i=1,Ni
      correct = ((block%bc%BCk(i,j,Nk+1)%tp==bc_ext).OR. &
                 (block%bc%BCk(i,j,Nk+1)%tp==bc_ref).OR. &
                 (block%bc%BCk(i,j,Nk+1)%tp==bc_in1).OR. &
                 (block%bc%BCk(i,j,Nk+1)%tp==bc_in2).OR. &
                 (block%mesh%V(i,j,Nk+1)<(0.2_R_P*block%mesh%V(i,j,Nk))))
      wall    = ((block%bc%BCk(i,j,Nk+1)%tp==bc_ext).OR.(block%bc%BCk(i,j,Nk+1)%tp==bc_ref))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
         ! normal metrics
         sn=2._R_P*((block%mesh%NFk(i,j,Nk-1)*block%mesh%Sk(i,j,Nk-1)).dot.block%mesh%NFk(i,j,Nk))
         block%mesh%NFk(i,  j,  Nk+1) = -(block%mesh%NFk(i,j,Nk-1)*block%mesh%Sk(i,j,Nk-1)) + sn*block%mesh%NFk(i,j,Nk)
         block%mesh%Sk (i,  j,  Nk+1) = normL2(   block%mesh%NFk(i,j,Nk+1))
         block%mesh%NFk(i,  j,  Nk+1) = normalize(block%mesh%NFk(i,j,Nk+1))
         ! tangential metrics
         block%mesh%NFi(i  ,j  ,Nk+1) = tm*block%mesh%NFi(i  ,j  ,Nk)
         block%mesh%NFi(i-1,j  ,Nk+1) = tm*block%mesh%NFi(i-1,j  ,Nk)
         block%mesh%NFi(i  ,j-1,Nk+1) = tm*block%mesh%NFi(i  ,j-1,Nk)
         block%mesh%NFi(i-1,j-1,Nk+1) = tm*block%mesh%NFi(i-1,j-1,Nk)
         block%mesh%Si (i  ,j  ,Nk+1) = tm*block%mesh%Si (i  ,j  ,Nk)
         block%mesh%Si (i-1,j  ,Nk+1) = tm*block%mesh%Si (i-1,j  ,Nk)
         block%mesh%Si (i  ,j-1,Nk+1) = tm*block%mesh%Si (i  ,j-1,Nk)
         block%mesh%Si (i-1,j-1,Nk+1) = tm*block%mesh%Si (i-1,j-1,Nk)

         block%mesh%NFj(i  ,j  ,Nk+1) = tm*block%mesh%NFj(i  ,j  ,Nk)
         block%mesh%NFj(i-1,j  ,Nk+1) = tm*block%mesh%NFj(i-1,j  ,Nk)
         block%mesh%NFj(i  ,j-1,Nk+1) = tm*block%mesh%NFj(i  ,j-1,Nk)
         block%mesh%NFj(i-1,j-1,Nk+1) = tm*block%mesh%NFj(i-1,j-1,Nk)
         block%mesh%Sj (i  ,j  ,Nk+1) = tm*block%mesh%Sj (i  ,j  ,Nk)
         block%mesh%Sj (i-1,j  ,Nk+1) = tm*block%mesh%Sj (i-1,j  ,Nk)
         block%mesh%Sj (i  ,j-1,Nk+1) = tm*block%mesh%Sj (i  ,j-1,Nk)
         block%mesh%Sj (i-1,j-1,Nk+1) = tm*block%mesh%Sj (i-1,j-1,Nk)
         ! volume
         block%mesh%V(  i,  j,  Nk+1) = block%mesh%V(     i,  j,  Nk)
      end if
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine mesh_metrics_correction
endmodule Lib_Mesh
