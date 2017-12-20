!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_ICEMCFDDerivedType Data_Type_File_ICEMCFD
!> Module definition of Type_File_ICEMCFD
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_ICEMCFDPrivateProcedure Data_Type_File_ICEMCFD
!> Module definition of Type_File_ICEMCFD
!> @}

!> @brief Module Data_Type_File_ICEMCFD contains the definition of Type_File_ICEMCFD, that is
!> the file handle for Ansys ICEM CFD files.
module Data_Type_File_ICEMCFD
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_ICEMCFD.
!> @ingroup Data_Type_File_ICEMCFDDerivedType
type, public, extends(Type_File_Base):: Type_File_ICEMCFD
  !contains
endtype Type_File_ICEMCFD
!-----------------------------------------------------------------------------------------------------------------------------------
!contains
  !> @ingroup Data_Type_File_ICEMCFDPrivateProcedure
  !> @{
!function load_icemcfd(Nf,filenames) result(err)
!!---------------------------------------------------------------------------------------------------------------------------------
!! Function for loading icemcfd fils.
!!---------------------------------------------------------------------------------------------------------------------------------
!
!!---------------------------------------------------------------------------------------------------------------------------------
!implicit none
!integer(I_P), intent(IN):: Nf               ! Number of input files.
!character(*), intent(IN):: filenames(1:Nf)  ! Base name of icemcfd files.
!integer(I_P)::             err              ! Error trapping flag: 0 no errors, >0 error occurs.
!integer(I_P)::             UnitFree         ! Free logic unit.
!integer(I_P)::             Unit_itc         ! Free logic unit for initial conditions files.
!integer(I_P)::             Unit_gc          ! Free logic unit for ghost cells definition.
!logical::                  is_file          ! Flag for inquiring the presence of file.
!integer(I_P)::             b,f,v,b1,b2,l    ! Counters.
!integer(I_P)::             bc               ! Boundary conditions counters.
!integer(I_P)::             Nb_l             ! Number of local blocks.
!character(500)::           line,line1,line2 ! Dummy strings for parsing stuff.
!integer(I_P), parameter::  No = 48          ! Number of possible orientations.
!integer(I_P)::             o,o1,o2          ! Orientation counters.
!integer(I_P)::             d1_min(1:3),d1_max(1:3)
!integer(I_P)::             d2_min(1:3),d2_max(1:3)
!integer(I_P)::             i1_min,i1_max,j1_min,j1_max,k1_min,k1_max
!integer(I_P)::             i2_min,i2_max,j2_min,j2_max,k2_min,k2_max
!integer(I_P)::             i1_sgn,j1_sgn,k1_sgn
!integer(I_P)::             i2_sgn,j2_sgn,k2_sgn
!character(6), parameter::  orientation(1:No) = (/'-i j k', '-i j-k', '-i k j', '-i k-j', '-i-j k', '-i-j-k', & ! Orientation list.
!                                                 '-i-k j', '-i-k-j', '-j i k', '-j i-k', '-j k i', '-j k-i', &
!                                                 '-j-i k', '-j-i-k', '-j-k i', '-j-k-i', '-k i j', '-k i-j', &
!                                                 '-k j i', '-k j-i', '-k-i j', '-k-i-j', '-k-j i', '-k-j-i', &
!                                                 ' i j k', ' i j-k', ' i k j', ' i k-j', ' i-j k', ' i-j-k', &
!                                                 ' i-k j', ' i-k-j', ' j i k', ' j i-k', ' j k i', ' j k-i', &
!                                                 ' j-i k', ' j-i-k', ' j-k i', ' j-k-i', ' k i j', ' k i-j', &
!                                                 ' k j i', ' k j-i', ' k-i j', ' k-i-j', ' k-j i', ' k-j-i'/)
!character(1), parameter::  tab = char(9)
!integer(I_P)::             blk_map(1:Nf)    ! Map of blocks over icem files.
!!---------------------------------------------------------------------------------------------------------------------------------
!
!!---------------------------------------------------------------------------------------------------------------------------------
!! verifying the presence of input files
!do f=1,Nf
!  inquire(file=adjustl(trim(filenames(f)))//'.geo',exist=is_file,iostat=err)
!  if (.NOT.is_file) then
!    call File_Not_Found(filename=adjustl(trim(filenames(f)))//'.geo',cpn='load_icemcfd')
!    stop
!  endif
!  inquire(file=adjustl(trim(filenames(f)))//'.topo',exist=is_file,iostat=err)
!  if (.NOT.is_file) then
!    call File_Not_Found(filename=adjustl(trim(filenames(f)))//'.topo',cpn='load_blocks')
!    stop
!  endif
!enddo
!! computing the number of blocks
!global%mesh_dims%Nb_tot = 0
!do f=1,Nf
!  ! opening geometry file
!  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filenames(f)))//'.geo', status = 'OLD', action = 'READ', form = 'FORMATTED')
!  blk_map(f)=0
!  do
!    read(UnitFree,'(A)',iostat=err) line
!    if (err /= 0) exit
!    if (index(line,'domain')>0) blk_map(f) = blk_map(f) + 1
!  enddo
!  close(UnitFree)                                      ! close file
!  global%mesh_dims%Nb_tot = global%mesh_dims%Nb_tot + blk_map(f)
!enddo
!! allocating array for blocks informations
!if (allocated(blocks))  then
!  do b=lbound(blocks,dim=1),ubound(blocks,dim=1)
!    call blocks(b)%P%free
!    call blocks(b)%bc%free
!  enddo
!  deallocate(blocks)
!endif
!allocate(blocks(1:global%mesh_dims%Nb_tot))
!call blocks%P%init(Ns=global%Ns)
!if (allocated(global%dfile%IBM_Scratch)) deallocate(global%dfile%IBM_Scratch)
!allocate(global%dfile%IBM_Scratch(0:3,1:global%mesh_dims%Nb_tot,1:global%mesh_dims%Nl))
!! reading the number of cells for each blocks contained into the file
!b = 0
!do f=1,Nf
!  ! opening geometry file
!  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filenames(f)))//'.geo', status = 'OLD', action = 'READ', form = 'FORMATTED')
!  do
!    read(UnitFree,'(A)',iostat=err) line
!    if (err /= 0) exit
!    if (index(line,'domain')>0)  then
!      b = b + 1 ! updating the block counter
!      read(line(index(line,' ')+1:),*)blocks(b)%Ni,blocks(b)%Nj,blocks(b)%Nk
!      blocks(b)%Ni = blocks(b)%Ni - 1_I_P
!      blocks(b)%Nj = blocks(b)%Nj - 1_I_P
!      blocks(b)%Nk = blocks(b)%Nk - 1_I_P
!      open(unit=Get_Unit(Unit_gc), file='BLK'//trim(str(.true.,b))//'.gc', status='OLD', action='READ', form='FORMATTED')
!      read(Unit_gc,*)(blocks(b)%gc(i),i=1,6)
!      close(Unit_gc)
!      write(stdout,*)
!      write(stdout,'(A,'//FI_P//')') ' Ghost cells of block: ',b
!      write(stdout,'(A,'//FI_P//')') ' Ghost cells:          ',blocks(b)%gc
!      write(stdout,*)
!    endif
!  enddo
!  close(UnitFree)                                      ! close file
!enddo
!! for efficiency each block is generated alone and then stored in scratch file
!! opening the scratch file where the current block is temporarily stored
!do l=1,global%mesh_dims%Nl
!  do b=1,global%mesh_dims%Nb_tot
!    if (l==1) then
!      open(unit=Get_Unit(global%dfile%IBM_Scratch(0,b,l)),form=  'FORMATTED',status='SCRATCH',iostat=err) ! topo file
!    endif
!      open(unit=Get_Unit(global%dfile%IBM_Scratch(1,b,l)),form='UNFORMATTED',status='SCRATCH',iostat=err) ! geo file
!      open(unit=Get_Unit(global%dfile%IBM_Scratch(2,b,l)),form='UNFORMATTED',status='SCRATCH',iostat=err) ! bco file
!      open(unit=Get_Unit(global%dfile%IBM_Scratch(3,b,l)),form='UNFORMATTED',status='SCRATCH',iostat=err) ! itc file
!  enddo
!enddo
!! the number of global blocks is artificially set to 1
!global%Nb = 1
!! allocating block-level data
!if (allocated(block)) then
!  do l=lbound(block,dim=2),ubound(block,dim=2)
!    do b=lbound(block,dim=1),ubound(block,dim=1)
!      call block(b,l)%free
!    enddo
!  enddo
!  deallocate(block)
!endif
!allocate(block(1:global%Nb,1:global%mesh_dims%Nl))
!do l=1,global%mesh_dims%Nl
!  do b=1,global%Nb
!    block(b,l)%global => global
!  enddo
!enddo
!! reading mesh, boundary and initial conditions of each block and storing in scratch files
!write(stdout,'(A)')'  Reading nodes coordinates from icemcfd files'
!b = 0
!do f=1,Nf
!  ! opening geometry file
!  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filenames(f)))//'.geo', status = 'OLD', action = 'READ', form = 'FORMATTED')
!  do
!    read(UnitFree,'(A)',iostat=err) line
!    if (err /= 0) exit
!    if (index(line,'domain')>0)  then
!      b = b + 1 ! updating the block counter
!      ! setting mesh dimensions for grid level 1
!      block(1,1)%gc = blocks(b)%gc
!      block(1,1)%Ni = blocks(b)%Ni
!      block(1,1)%Nj = blocks(b)%Nj
!      block(1,1)%Nk = blocks(b)%Nk
!      ! computing mesh dimensions for other grid levels
!      if (global%mesh_dims%Nl>1) then
!        do l=2,global%mesh_dims%Nl
!          ! computing the number of cells of coarser grid levels
!          if     (mod(block(1,l-1)%Ni,2)/=0) then
!            write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
!            write(stderr,'(A,I4)')' Impossible to compute grid level ',l
!            write(stderr,'(A,I4)')' Inconsistent direction i, Ni ',block(1,l-1)%Ni
!            write(stderr,'(A,I4)')' level ',l-1
!            write(stderr,'(A,I4)')' block ',b
!            stop
!          elseif (mod(block(1,l-1)%Nj,2)/=0) then
!            write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
!            write(stderr,'(A,I4)')' Impossible to compute grid level ',l
!            write(stderr,'(A,I4)')' Inconsistent direction j, Nj ',block(1,l-1)%Nj
!            write(stderr,'(A,I4)')' level ',l-1
!            write(stderr,'(A,I4)')' block ',b
!            stop
!          elseif (mod(block(1,l-1)%Nk,2)/=0) then
!            write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
!            write(stderr,'(A,I4)')' Impossible to compute grid level ',l
!            write(stderr,'(A,I4)')' Inconsistent direction k, Nk ',block(1,l-1)%Nk
!            write(stderr,'(A,I4)')' level ',l-1
!            write(stderr,'(A,I4)')' block ',b
!            stop
!          endif
!          block(1,l)%gc = block(1,l-1)%gc
!          block(1,l)%Ni = block(1,l-1)%Ni/2
!          block(1,l)%Nj = block(1,l-1)%Nj/2
!          block(1,l)%Nk = block(1,l-1)%Nk/2
!        enddo
!      endif
!        ! allocating block
!      do l=1,global%mesh_dims%Nl
!        call block(1,l)%alloc
!      enddo
!      ! reading node coordinates
!      do k=0,block(1,1)%Nk
!        do j=0,block(1,1)%Nj
!          do i=0,block(1,1)%Ni
!            read(UnitFree,*)block(1,1)%node(i,j,k)%x,block(1,1)%node(i,j,k)%y,block(1,1)%node(i,j,k)%z
!          enddo
!        enddo
!      enddo
!      ! computing the node of other grid levels if used
!      if (global%mesh_dims%Nl>1) then
!        do l=2,global%mesh_dims%Nl
!          do k=0,block(1,l)%Nk
!            do j=0,block(1,l)%Nj
!              do i=0,block(1,l)%Ni
!                block(1,l)%node(i,j,k) = block(1,l-1)%node(i*2,j*2,k*2)
!              enddo
!            enddo
!          enddo
!        enddo
!      endif
!      ! storing the mesh in the scratch files
!      do l=1,global%mesh_dims%Nl
!        err = write_vector(array3D=block(1,l)%node,unit=global%dfile%IBM_Scratch(1,b,l))
!        !err = write_cell(  array3D=block(1,l)%cell,unit=global%dfile%IBM_Scratch(1,b,l))
!        ! rewinding scratch file
!        rewind(global%dfile%IBM_Scratch(1,b,l))
!      enddo
!    endif
!  enddo
!  close(UnitFree)                                      ! close file
!enddo
!write(stdout,'(A)')'  Reading boundary and initial conditions from icemcfd files'
!! the topological files are pre-processed and arranged one for block in the scratch files
!b = 0
!do f=1,Nf
!  Nb_l = 0 ; if (f>1) Nb_l = sum(blk_map(1:f-1))
!  ! opening geometry file
!  open(unit=Get_Unit(UnitFree), file=adjustl(trim(filenames(f)))//'.topo', status='OLD', action='READ', form='FORMATTED')
!  do
!    read(UnitFree,'(A)',iostat=err) line
!    if (err /= 0) exit
!    if (index(line,'# Connectivity for domain.')>0)  then
!      ! found a block connectivity information
!      read(line(index(line,'# Connectivity for domain.')+26:),*)b ; b = b + Nb_l ! current block
!      write(global%dfile%IBM_Scratch(0,b,1),'(A)')'# Connectivity for domain.'//trim(str(.true.,b))
!      do
!        read(UnitFree,'(A)',iostat=err) line1
!        if (err /= 0) exit
!        if (index(line1,tab//'f')>0) then
!          read(UnitFree,'(A)',iostat=err) line2
!          if (err /= 0) exit
!          read(line1(index(line1,'domain.')+7:),*)b1 ; b1 = b1 + Nb_l ! reading the number of current block
!          read(line2(index(line2,'domain.')+7:),*)b2 ; b2 = b2 + Nb_l ! reading the number of connected block
!          write(global%dfile%IBM_Scratch(0,b,1),'(A)')line1(:index(line1,' domain.')+7)//trim(str(.true.,b1))//&
!trim(line1(index(line1,tab)-1:))
!          write(global%dfile%IBM_Scratch(0,b,1),'(A)')line2(:index(line2,' domain.')+7)//trim(str(.true.,b2))//&
!trim(line2(index(line2,tab)-1:))
!        elseif (index(line1,tab//'b')>0) then
!          cycle
!        elseif (index(line1,tab//'e')>0) then
!          cycle
!        elseif (index(line1,tab//'v')>0) then
!          cycle
!        else
!          write(global%dfile%IBM_Scratch(0,b,1),*)
!          exit
!        endif
!      enddo
!    elseif (index(line,'# Boundary conditions and/or properties for domain.')>0)  then
!      ! found a block boundary or initial conditions information
!      read(line(index(line,'# Boundary conditions and/or properties for domain.')+51:),*)b ; b = b + Nb_l ! current block
!      write(global%dfile%IBM_Scratch(0,b,1),'(A)')'# Boundary conditions and/or properties for domain.'//trim(str(.true.,b))
!      do
!        read(UnitFree,'(A)',iostat=err) line1
!        if (err /= 0) exit
!        if (index(line1,tab//'f')>0) then
!          write(global%dfile%IBM_Scratch(0,b,1),'(A)')trim(line1)
!        elseif (index(line1,tab//'b')>0) then
!          write(global%dfile%IBM_Scratch(0,b,1),'(A)')trim(line1)
!        elseif (index(line1,tab//'e')>0) then
!          cycle
!        elseif (index(line1,tab//'v')>0) then
!          cycle
!        else
!          write(global%dfile%IBM_Scratch(0,b,1),*)
!          exit
!        endif
!      enddo
!    endif
!  enddo
!  close(UnitFree)
!enddo
!do b=1,global%mesh_dims%Nb_tot
!  rewind(global%dfile%IBM_Scratch(0,b,1))
!enddo
!! reading boundary and initial conditions from scratch topological files
!do b=1,global%mesh_dims%Nb_tot
!  ! setting mesh dimensions for grid level 1
!  block(1,1)%Ni = blocks(b)%Ni
!  block(1,1)%Nj = blocks(b)%Nj
!  block(1,1)%Nk = blocks(b)%Nk
!  block(1,1)%gc = blocks(b)%gc
!  ! computing mesh dimensions for other grid levels
!  if (global%mesh_dims%Nl>1) then
!    do l=2,global%mesh_dims%Nl
!      ! computing the number of cells of coarser grid levels
!      if     (mod(block(1,l-1)%Ni,2)/=0) then
!        write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
!        write(stderr,'(A,I4)')' Impossible to compute grid level ',l
!        write(stderr,'(A,I4)')' Inconsistent direction i, Ni ',block(1,l-1)%Ni
!        write(stderr,'(A,I4)')' level ',l-1
!        write(stderr,'(A,I4)')' block ',b
!        stop
!      elseif (mod(block(1,l-1)%Nj,2)/=0) then
!        write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
!        write(stderr,'(A,I4)')' Impossible to compute grid level ',l
!        write(stderr,'(A,I4)')' Inconsistent direction j, Nj ',block(1,l-1)%Nj
!        write(stderr,'(A,I4)')' level ',l-1
!        write(stderr,'(A,I4)')' block ',b
!        stop
!      elseif (mod(block(1,l-1)%Nk,2)/=0) then
!        write(stderr,'(A)')   ' Attention the number of grid levels used is not consistent with the number of cells'
!        write(stderr,'(A,I4)')' Impossible to compute grid level ',l
!        write(stderr,'(A,I4)')' Inconsistent direction k, Nk ',block(1,l-1)%Nk
!        write(stderr,'(A,I4)')' level ',l-1
!        write(stderr,'(A,I4)')' block ',b
!        stop
!      endif
!      block(1,l)%gc = block(1,l-1)%gc
!      block(1,l)%Ni = block(1,l-1)%Ni/2
!      block(1,l)%Nj = block(1,l-1)%Nj/2
!      block(1,l)%Nk = block(1,l-1)%Nk/2
!    enddo
!  endif
!    ! allocating block
!  do l=1,global%mesh_dims%Nl
!    call block(1,l)%alloc
!  enddo
!  do
!    read(global%dfile%IBM_Scratch(0,b,1),'(A)',iostat=err) line
!    if (err /= 0) exit
!    if (index(line,'# Connectivity for domain.')>0)  then
!      ! found a block connectivity information
!      do
!        read(global%dfile%IBM_Scratch(0,b,1),'(A)',iostat=err) line1
!        if (err /= 0) exit
!        if (index(line1,tab//'f')>0) then
!          read(global%dfile%IBM_Scratch(0,b,1),'(A)',iostat=err) line2
!          if (err /= 0) exit
!          read(line1(index(line1,'domain.')+7:),*)b1 ! reading the number of current block
!          read(line2(index(line2,'domain.')+7:),*)b2 ! reading the number of connected block
!          read(line1(index(line1,tab//'f',BACK=.true.)+2:),*)d1_min(1),d1_min(2),d1_min(3),d1_max(1),d1_max(2),d1_max(3)
!          read(line2(index(line2,tab//'f',BACK=.true.)+2:),*)d2_min(1),d2_min(2),d2_min(3),d2_max(1),d2_max(2),d2_max(3)
!          ! reading the orientation of the 2 blocks
!          do o=1,No
!            if (index(line1,orientation(o))>0) o1 = o
!            if (index(line2,orientation(o))>0) o2 = o
!          enddo
!          ! mapping the block orientation in i,j,k direction
!          ! block 1
!          select case(orientation(o1)(2:2))
!          case('i')
!            if (orientation(o1)(1:1)==' ') then
!              i1_min = d1_min(1) - 1 ; i1_max = d1_max(1) - 1
!              i1_sgn =  1
!            else
!              i1_min = d1_max(1) - 1 ; i1_max = d1_min(1) - 1
!              i1_sgn = -1
!            endif
!          case('j')
!            if (orientation(o1)(1:1)==' ') then
!              j1_min = d1_min(1) - 1 ; j1_max = d1_max(1) - 1
!              j1_sgn =  1
!            else
!              j1_min = d1_max(1) - 1 ; j1_max = d1_min(1) - 1
!              j1_sgn = -1
!            endif
!          case('k')
!            if (orientation(o1)(1:1)==' ') then
!              k1_min = d1_min(1) - 1 ; k1_max = d1_max(1) - 1
!              k1_sgn =  1
!            else
!              k1_min = d1_max(1) - 1 ; k1_max = d1_max(1) - 1
!              k1_sgn = -1
!            endif
!          endselect
!          select case(orientation(o1)(4:4))
!          case('i')
!            if (orientation(o1)(3:3)==' ') then
!              i1_min = d1_min(2) - 1 ; i1_max = d1_max(2) - 1
!              i1_sgn =  1
!            else
!              i1_min = d1_max(2) - 1 ; i1_max = d1_min(2) - 1
!              i1_sgn = -1
!            endif
!          case('j')
!            if (orientation(o1)(3:3)==' ') then
!              j1_min = d1_min(2) - 1 ; j1_max = d1_max(2) - 1
!              j1_sgn =  1
!            else
!              j1_min = d1_max(2) - 1 ; j1_max = d1_min(2) - 1
!              j1_sgn = -1
!            endif
!          case('k')
!            if (orientation(o1)(3:3)==' ') then
!              k1_min = d1_min(2) - 1 ; k1_max = d1_max(2) - 1
!              k1_sgn =  1
!            else
!              k1_min = d1_max(2) - 1 ; k1_max = d1_min(2) - 1
!              k1_sgn = -1
!            endif
!          endselect
!          select case(orientation(o1)(6:6))
!          case('i')
!            if (orientation(o1)(5:5)==' ') then
!              i1_min = d1_min(3) - 1 ; i1_max = d1_max(3) - 1
!              i1_sgn =  1
!            else
!              i1_min = d1_max(3) - 1 ; i1_max = d1_min(3) - 1
!              i1_sgn = -1
!            endif
!          case('j')
!            if (orientation(o1)(5:5)==' ') then
!              j1_min = d1_min(3) - 1 ; j1_max = d1_max(3) - 1
!              j1_sgn =  1
!            else
!              j1_min = d1_max(3) - 1 ; j1_max = d1_min(3) - 1
!              j1_sgn = -1
!            endif
!          case('k')
!            if (orientation(o1)(5:5)==' ') then
!              k1_min = d1_min(3) - 1 ; k1_max = d1_max(3) - 1
!              k1_sgn =  1
!            else
!              k1_min = d1_max(3) - 1 ; k1_max = d1_min(3) - 1
!              k1_sgn = -1
!            endif
!          endselect
!          ! block 2
!          select case(orientation(o2)(2:2))
!          case('i')
!            if (orientation(o2)(1:1)==' ') then
!              i2_min = d2_min(1) - 1 ; i2_max = d2_max(1) - 1
!              i2_sgn =  1
!            else
!              i2_min = d2_max(1) - 1 ; i2_max = d2_min(1) - 1
!              i2_sgn = -1
!            endif
!          case('j')
!            if (orientation(o2)(1:1)==' ') then
!              j2_min = d2_min(1) - 1 ; j2_max = d2_max(1) - 1
!              j2_sgn =  1
!            else
!              j2_min = d2_max(1) - 1 ; j2_max = d2_min(1) - 1
!              j2_sgn = -1
!            endif
!          case('k')
!            if (orientation(o2)(1:1)==' ') then
!              k2_min = d2_min(1) - 1 ; k2_max = d2_max(1) - 1
!              k2_sgn =  1
!            else
!              k2_min = d2_max(1) - 1 ; k2_max = d2_min(1) - 1
!              k2_sgn = -1
!            endif
!          endselect
!          select case(orientation(o2)(4:4))
!          case('i')
!            if (orientation(o2)(3:3)==' ') then
!              i2_min = d2_min(2) - 1 ; i2_max = d2_max(2) - 1
!              i2_sgn =  1
!            else
!              i2_min = d2_max(2) - 1 ; i2_max = d2_min(2) - 1
!              i2_sgn = -1
!            endif
!          case('j')
!            if (orientation(o2)(3:3)==' ') then
!              j2_min = d2_min(2) - 1 ; j2_max = d2_max(2) - 1
!              j2_sgn =  1
!            else
!              j2_min = d2_max(2) - 1 ; j2_max = d2_min(2) - 1
!              j2_sgn = -1
!            endif
!          case('k')
!            if (orientation(o2)(3:3)==' ') then
!              k2_min = d2_min(2) - 1 ; k2_max = d2_max(2) - 1
!              k2_sgn =  1
!            else
!              k2_min = d2_max(2) - 1 ; k2_max = d2_min(2) - 1
!              k2_sgn = -1
!            endif
!          endselect
!          select case(orientation(o2)(6:6))
!          case('i')
!            if (orientation(o2)(5:5)==' ') then
!              i2_min = d2_min(3) - 1 ; i2_max = d2_max(3) - 1
!              i2_sgn =  1
!            else
!              i2_min = d2_max(3) - 1 ; i2_max = d2_min(3) - 1
!              i2_sgn = -1
!            endif
!          case('j')
!            if (orientation(o2)(5:5)==' ') then
!              j2_min = d2_min(3) - 1 ; j2_max = d2_max(3) - 1
!              j2_sgn =  1
!            else
!              j2_min = d2_max(3) - 1 ; j2_max = d2_min(3) - 1
!              j2_sgn = -1
!            endif
!          case('k')
!            if (orientation(o2)(5:5)==' ') then
!              k2_min = d2_min(3) - 1 ; k2_max = d2_max(3) - 1
!              k2_sgn =  1
!            else
!              k2_min = d2_max(3) - 1 ; k2_max = d2_min(3) - 1
!              k2_sgn = -1
!            endif
!          endselect
!          ! imposing the adjacent boundary conditions for all grid levels
!          do l=1,global%mesh_dims%Nl
!            select case(orientation(o1)(6:6))
!            case('i')
!              ! the adjacent face is i-face for block 1
!              if (i1_min==0) then
!                i1_min = 1-block(1,l)%gc(1)
!                i1_max = 0
!              else
!                i1_min = blocks(b1)%Ni/(2**(l-1))
!                i1_max = blocks(b1)%Ni/(2**(l-1))+block(1,l)%gc(2)-1
!              endif
!              do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
!                do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
!                  do i=i1_min,i1_max
!                    block(1,l)%Fi(i,j,k)%BC%tp    = bc_adj ; call block(1,l)%Fi(i,j,k)%BC%init
!                    block(1,l)%Fi(i,j,k)%BC%adj%b = b2
!                    do o=2,4,2
!                      select case(orientation(o1)(o:o))
!                      case('j')
!                        select case(orientation(o2)(o:o))
!                        case('i')
!                          if (j1_sgn*i2_sgn<0) then
!                            block(1,l)%Fi(i,j,k)%BC%adj%i = j1_max/2**(l-1) - j
!                          else
!                            block(1,l)%Fi(i,j,k)%BC%adj%i = j
!                          endif
!                        case('j')
!                          if (j1_sgn*j2_sgn<0) then
!                            block(1,l)%Fi(i,j,k)%BC%adj%j = j1_max/2**(l-1) - j
!                          else
!                            block(1,l)%Fi(i,j,k)%BC%adj%j = j
!                          endif
!                        case('k')
!                          if (j1_sgn*k2_sgn<0) then
!                            block(1,l)%Fi(i,j,k)%BC%adj%k = j1_max/2**(l-1) - j
!                          else
!                            block(1,l)%Fi(i,j,k)%BC%adj%k = j
!                          endif
!                        endselect
!                      case('k')
!                        select case(orientation(o2)(o:o))
!                        case('i')
!                          if (k1_sgn*i2_sgn<0) then
!                            block(1,l)%Fi(i,j,k)%BC%adj%i = k1_max/2**(l-1) - k
!                          else
!                            block(1,l)%Fi(i,j,k)%BC%adj%i = k
!                          endif
!                        case('j')
!                          if (k1_sgn*j2_sgn<0) then
!                            block(1,l)%Fi(i,j,k)%BC%adj%j = k1_max/2**(l-1) - k
!                          else
!                            block(1,l)%Fi(i,j,k)%BC%adj%j = k
!                          endif
!                        case('k')
!                          if (k1_sgn*k2_sgn<0) then
!                            block(1,l)%Fi(i,j,k)%BC%adj%k = k1_max/2**(l-1) - k
!                          else
!                            block(1,l)%Fi(i,j,k)%BC%adj%k = k
!                          endif
!                        endselect
!                      endselect
!                    enddo
!                    select case(orientation(o2)(6:6))
!                    case('i')
!                      if (i1_max==0) then
!                        block(1,l)%Fi(i,j,k)%BC%adj%i = (blocks(b2)%Ni/(2**(l-1))+i)*(1-i1_sgn*i2_sgn)/2 + &
!                                                        (-i+1                      )*(1+i1_sgn*i2_sgn)/2
!                      else
!                        block(1,l)%Fi(i,j,k)%BC%adj%i =                                                   &
!                          (i-blocks(b1)%Ni/(2**(l-1))+1                           )*(1-i1_sgn*i2_sgn)/2 + &
!                          (  blocks(b2)%Ni/(2**(l-1))-(i-blocks(b1)%Ni/(2**(l-1))))*(1+i1_sgn*i2_sgn)/2
!                      endif
!                    case('j')
!                      if (i1_max==0) then
!                        block(1,l)%Fi(i,j,k)%BC%adj%j = (blocks(b2)%Nj/(2**(l-1))+i)*(1-i1_sgn*j2_sgn)/2 + &
!                                                        (-i+1                      )*(1+i1_sgn*j2_sgn)/2
!                      else
!                        block(1,l)%Fi(i,j,k)%BC%adj%j =                                                   &
!                          (i-blocks(b1)%Ni/(2**(l-1))+1                           )*(1-i1_sgn*j2_sgn)/2 + &
!                          (  blocks(b2)%Nj/(2**(l-1))-(i-blocks(b1)%Ni/(2**(l-1))))*(1+i1_sgn*j2_sgn)/2
!                      endif
!                    case('k')
!                      if (i1_max==0) then
!                        block(1,l)%Fi(i,j,k)%BC%adj%k = (blocks(b2)%Nk/(2**(l-1))+i)*(1-i1_sgn*k2_sgn)/2 + &
!                                                        (-i+1                      )*(1+i1_sgn*k2_sgn)/2
!                      else
!                        block(1,l)%Fi(i,j,k)%BC%adj%k =                                                   &
!                          (i-blocks(b1)%Ni/(2**(l-1))+1                           )*(1-i1_sgn*k2_sgn)/2 + &
!                          (  blocks(b2)%Nk/(2**(l-1))-(i-blocks(b1)%Ni/(2**(l-1))))*(1+i1_sgn*k2_sgn)/2
!                      endif
!                    endselect
!                  enddo
!                enddo
!              enddo
!            case('j')
!              ! the adjacent face is j-face for block 1
!              if (j1_min==0) then
!                j1_min = 1-block(1,l)%gc(3)
!                j1_max = 0
!              else
!                j1_min = blocks(b1)%Nj/(2**(l-1))
!                j1_max = blocks(b1)%Nj/(2**(l-1))+block(1,l)%gc(4)-1
!              endif
!              do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
!                do j=j1_min,j1_max
!                  do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
!                    block(1,l)%Fj(i,j,k)%BC%tp    = bc_adj ; call block(1,l)%Fj(i,j,k)%BC%init
!                    block(1,l)%Fj(i,j,k)%BC%adj%b = b2
!                    do o=2,4,2
!                      select case(orientation(o1)(o:o))
!                      case('i')
!                        select case(orientation(o2)(o:o))
!                        case('i')
!                          if (i1_sgn*i2_sgn<0) then
!                            block(1,l)%Fj(i,j,k)%BC%adj%i = i1_max/2**(l-1) - i
!                          else
!                            block(1,l)%Fj(i,j,k)%BC%adj%i = i
!                          endif
!                        case('j')
!                          if (i1_sgn*j2_sgn<0) then
!                            block(1,l)%Fj(i,j,k)%BC%adj%j = i1_max/2**(l-1) - i
!                          else
!                            block(1,l)%Fj(i,j,k)%BC%adj%j = i
!                          endif
!                        case('k')
!                          if (i1_sgn*k2_sgn<0) then
!                            block(1,l)%Fj(i,j,k)%BC%adj%k = i1_max/2**(l-1) - i
!                          else
!                            block(1,l)%Fj(i,j,k)%BC%adj%k = i
!                          endif
!                        endselect
!                      case('k')
!                        select case(orientation(o2)(o:o))
!                        case('i')
!                          if (k1_sgn*i2_sgn<0) then
!                            block(1,l)%Fj(i,j,k)%BC%adj%i = k1_max/2**(l-1) - k
!                          else
!                            block(1,l)%Fj(i,j,k)%BC%adj%i = k
!                          endif
!                        case('j')
!                          if (k1_sgn*j2_sgn<0) then
!                            block(1,l)%Fj(i,j,k)%BC%adj%j = k1_max/2**(l-1) - k
!                          else
!                            block(1,l)%Fj(i,j,k)%BC%adj%j = k
!                          endif
!                        case('k')
!                          if (k1_sgn*k2_sgn<0) then
!                            block(1,l)%Fj(i,j,k)%BC%adj%k = k1_max/2**(l-1) - k
!                          else
!                            block(1,l)%Fj(i,j,k)%BC%adj%k = k
!                          endif
!                        endselect
!                      endselect
!                    enddo
!                    select case(orientation(o2)(6:6))
!                    case('i')
!                      if (j1_max==0) then
!                        block(1,l)%Fj(i,j,k)%BC%adj%i = (blocks(b2)%Ni/(2**(l-1))+j)*(1-j1_sgn*i2_sgn)/2 + &
!                                                        (-j+1                      )*(1+j1_sgn*i2_sgn)/2
!                      else
!                        block(1,l)%Fj(i,j,k)%BC%adj%i =                                                   &
!                          (j-blocks(b1)%Nj/(2**(l-1))+1                           )*(1-j1_sgn*i2_sgn)/2 + &
!                          (  blocks(b2)%Ni/(2**(l-1))-(j-blocks(b1)%Nj/(2**(l-1))))*(1+j1_sgn*i2_sgn)/2
!                      endif
!                    case('j')
!                      if (j1_max==0) then
!                        block(1,l)%Fj(i,j,k)%BC%adj%j = (blocks(b2)%Nj/(2**(l-1))+j)*(1-j1_sgn*j2_sgn)/2 + &
!                                                        (-j+1                      )*(1+j1_sgn*j2_sgn)/2
!                      else
!                        block(1,l)%Fj(i,j,k)%BC%adj%j =                                                   &
!                          (j-blocks(b1)%Nj/(2**(l-1))+1                           )*(1-j1_sgn*j2_sgn)/2 + &
!                          (  blocks(b2)%Nj/(2**(l-1))-(j-blocks(b1)%Nj/(2**(l-1))))*(1+j1_sgn*j2_sgn)/2
!                      endif
!                    case('k')
!                      if (j1_max==0) then
!                        block(1,l)%Fj(i,j,k)%BC%adj%k = (blocks(b2)%Nk/(2**(l-1))+j)*(1-j1_sgn*k2_sgn)/2 + &
!                                                        (-j+1                      )*(1+j1_sgn*k2_sgn)/2
!                      else
!                        block(1,l)%Fj(i,j,k)%BC%adj%k =                                                   &
!                          (j-blocks(b1)%Nj/(2**(l-1))+1                           )*(1-j1_sgn*k2_sgn)/2 + &
!                          (  blocks(b2)%Nk/(2**(l-1))-(j-blocks(b1)%Nj/(2**(l-1))))*(1+j1_sgn*k2_sgn)/2
!                      endif
!                    endselect
!                  enddo
!                enddo
!              enddo
!            case('k')
!              ! the adjacent face is k-face for block 1
!              if (k1_min==0) then
!                k1_min = 1-block(1,l)%gc(5)
!                k1_max = 0
!              else
!                k1_min = blocks(b1)%Nk/(2**(l-1))
!                k1_max = blocks(b1)%Nk/(2**(l-1))+block(1,l)%gc(6)-1
!              endif
!              do k=k1_min,k1_max
!                do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
!                  do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
!                    block(1,l)%Fk(i,j,k)%BC%tp    = bc_adj ; call block(1,l)%Fk(i,j,k)%BC%init
!                    block(1,l)%Fk(i,j,k)%BC%adj%b = b2
!                    do o=2,4,2
!                      select case(orientation(o1)(o:o))
!                      case('i')
!                        select case(orientation(o2)(o:o))
!                        case('i')
!                          if (i1_sgn*i2_sgn<0) then
!                            block(1,l)%Fk(i,j,k)%BC%adj%i = i1_max/2**(l-1) - i
!                          else
!                            block(1,l)%Fk(i,j,k)%BC%adj%i = i
!                          endif
!                        case('j')
!                          if (i1_sgn*j2_sgn<0) then
!                            block(1,l)%Fk(i,j,k)%BC%adj%j = i1_max/2**(l-1) - i
!                          else
!                            block(1,l)%Fk(i,j,k)%BC%adj%j = i
!                          endif
!                        case('k')
!                          if (i1_sgn*k2_sgn<0) then
!                            block(1,l)%Fk(i,j,k)%BC%adj%k = i1_max/2**(l-1) - i
!                          else
!                            block(1,l)%Fk(i,j,k)%BC%adj%k = i
!                          endif
!                        endselect
!                      case('j')
!                        select case(orientation(o2)(o:o))
!                        case('i')
!                          if (j1_sgn*i2_sgn<0) then
!                            block(1,l)%Fk(i,j,k)%BC%adj%i = j1_max/2**(l-1) - j
!                          else
!                            block(1,l)%Fk(i,j,k)%BC%adj%i = j
!                          endif
!                        case('j')
!                          if (j1_sgn*j2_sgn<0) then
!                            block(1,l)%Fk(i,j,k)%BC%adj%j = j1_max/2**(l-1) - j
!                          else
!                            block(1,l)%Fk(i,j,k)%BC%adj%j = j
!                          endif
!                        case('k')
!                          if (j1_sgn*k2_sgn<0) then
!                            block(1,l)%Fk(i,j,k)%BC%adj%k = j1_max/2**(l-1) - j
!                          else
!                            block(1,l)%Fk(i,j,k)%BC%adj%k = j
!                          endif
!                        endselect
!                      endselect
!                    enddo
!                    select case(orientation(o2)(6:6))
!                    case('i')
!                      if (k1_max==0) then
!                        block(1,l)%Fk(i,j,k)%BC%adj%i = (blocks(b2)%Ni/(2**(l-1))+k)*(1-k1_sgn*i2_sgn)/2 + &
!                                                        (-k+1                      )*(1+k1_sgn*i2_sgn)/2
!                      else
!                        block(1,l)%Fk(i,j,k)%BC%adj%i =                                                   &
!                          (k-blocks(b1)%Nk/(2**(l-1))+1                           )*(1-k1_sgn*i2_sgn)/2 + &
!                          (  blocks(b2)%Nk/(2**(l-1))-(k-blocks(b1)%Ni/(2**(l-1))))*(1+k1_sgn*i2_sgn)/2
!                      endif
!                    case('j')
!                      if (k1_max==0) then
!                        block(1,l)%Fk(i,j,k)%BC%adj%j = (blocks(b2)%Nj/(2**(l-1))+k)*(1-k1_sgn*j2_sgn)/2 + &
!                                                        (-k+1                      )*(1+k1_sgn*j2_sgn)/2
!                      else
!                        block(1,l)%Fk(i,j,k)%BC%adj%j =                                                   &
!                          (k-blocks(b1)%Nk/(2**(l-1))+1                           )*(1-k1_sgn*j2_sgn)/2 + &
!                          (  blocks(b2)%Nj/(2**(l-1))-(k-blocks(b1)%Nk/(2**(l-1))))*(1+k1_sgn*j2_sgn)/2
!                      endif
!                    case('k')
!                      if (k1_max==0) then
!                        block(1,l)%Fk(i,j,k)%BC%adj%k = (blocks(b2)%Nk/(2**(l-1))+k)*(1-k1_sgn*k2_sgn)/2 + &
!                                                        (-k+1                      )*(1+k1_sgn*k2_sgn)/2
!                      else
!                        block(1,l)%Fk(i,j,k)%BC%adj%k =                                                   &
!                          (k-blocks(b1)%Nk/(2**(l-1))+1                           )*(1-k1_sgn*k2_sgn)/2 + &
!                          (  blocks(b2)%Nk/(2**(l-1))-(k-blocks(b1)%Nk/(2**(l-1))))*(1+k1_sgn*k2_sgn)/2
!                      endif
!                    endselect
!                  enddo
!                enddo
!              enddo
!            endselect
!          enddo
!        elseif (index(line1,tab//'b')>0) then
!          cycle
!        elseif (index(line1,tab//'e')>0) then
!          cycle
!        elseif (index(line1,tab//'v')>0) then
!          cycle
!        else
!          exit
!        endif
!      enddo
!    elseif (index(line,'# Boundary conditions and/or properties for domain.')>0)  then
!      ! found a block boundary or initial conditions information
!      do
!        read(global%dfile%IBM_Scratch(0,b,1),'(A)',iostat=err) line1
!        if (err /= 0) exit
!        if (index(line1,tab//'f')>0) then
!          do bc=1,Nbc
!            if (index(line1,bc_list_str(bc))>0) then
!              read(line1(index(line1,tab//'f',BACK=.true.)+2:),*)i1_min,j1_min,k1_min,i1_max,j1_max,k1_max
!              i1_min = i1_min - 1 ; i1_max = i1_max - 1
!              j1_min = j1_min - 1 ; j1_max = j1_max - 1
!              k1_min = k1_min - 1 ; k1_max = k1_max - 1
!              do l=1,global%mesh_dims%Nl
!                if     (i1_min==i1_max) then
!                  ! i face
!                  if (i1_min==0) then
!                    do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
!                      do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
!                        do i=1-block(1,l)%gc(1),0
!                          block(1,l)%Fi(i,j,k)%BC%tp = bc_list(bc) ; call block(1,l)%Fi(i,j,k)%BC%init
!                          if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
!                            block(1,l)%Fi(i,j,k)%BC%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
!                                                                     index(line1,tab//'f')-1),1_I_P)
!                          endif
!                        enddo
!                      enddo
!                    enddo
!                  else
!                    do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
!                      do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
!                        do i=i1_max/2**(l-1)+1,i1_max/2**(l-1)+block(1,l)%gc(2)
!                          block(1,l)%Fi(i-1,j,k)%BC%tp = bc_list(bc) ; call block(1,l)%Fi(i-1,j,k)%BC%init
!                          if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
!                            block(1,l)%Fi(i-1,j,k)%BC%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
!                                                                       index(line1,tab//'f')-1),1_I_P)
!                          endif
!                        enddo
!                      enddo
!                    enddo
!                  endif
!                elseif (j1_min==j1_max) then
!                  ! j face
!                  if (j1_min==0) then
!                    do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
!                      do j=1-block(1,l)%gc(3),0
!                        do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
!                          block(1,l)%Fj(i,j,k)%BC%tp = bc_list(bc) ; call block(1,l)%Fj(i,j,k)%BC%init
!                          if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
!                            block(1,l)%Fj(i,j,k)%BC%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
!                                                                     index(line1,tab//'f')-1),1_I_P)
!                          endif
!                        enddo
!                      enddo
!                    enddo
!                  else
!                    do k=k1_min/2**(l-1)+1,k1_max/2**(l-1)
!                      do j=j1_max/2**(l-1)+1,j1_max/2**(l-1)+block(1,l)%gc(4)
!                        do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
!                          block(1,l)%Fj(i,j-1,k)%BC%tp = bc_list(bc) ; call block(1,l)%Fj(i,j-1,k)%BC%init
!                          if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
!                            block(1,l)%Fj(i,j-1,k)%BC%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
!                                                                       index(line1,tab//'f')-1),1_I_P)
!                          endif
!                        enddo
!                      enddo
!                    enddo
!                  endif
!                elseif (k1_min==k1_max) then
!                  ! k face
!                  if (k1_min==0) then
!                    do k=1-block(1,l)%gc(5),0
!                      do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
!                        do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
!                          block(1,l)%Fk(i,j,k)%BC%tp = bc_list(bc) ; call block(1,l)%Fk(i,j,k)%BC%init
!                          if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
!                            block(1,l)%Fk(i,j,k)%BC%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
!                                                                     index(line1,tab//'f')-1),1_I_P)
!                          endif
!                        enddo
!                      enddo
!                    enddo
!                  else
!                    do k=k1_max/2**(l-1)+1,k1_max/2**(l-1)+block(1,l)%gc(6)
!                      do j=j1_min/2**(l-1)+1,j1_max/2**(l-1)
!                        do i=i1_min/2**(l-1)+1,i1_max/2**(l-1)
!                          block(1,l)%Fk(i,j,k-1)%BC%tp = bc_list(bc) ; call block(1,l)%Fk(i,j,k-1)%BC%init
!                          if (bc_list(bc)==bc_in1.OR.bc_list(bc)==bc_in2) then
!                            block(1,l)%Fk(i,j,k-1)%BC%inf = cton(line1(index(line1,bc_list_str(bc))+3: &
!                                                                       index(line1,tab//'f')-1),1_I_P)
!                          endif
!                        enddo
!                      enddo
!                    enddo
!                  endif
!                endif
!              enddo
!            endif
!          enddo
!        elseif (index(line1,tab//'b')>0) then
!          if (index(line1,'BLK')>0) then
!            inquire(file=adjustl(trim(line1(1:index(line1,tab)-1)))//'.itc',exist=is_file)
!            if (.NOT.is_file) then
!              write(stderr,'(A)')' File'
!              write(stderr,'(A)')' '//adjustl(trim(line1(1:index(line1,tab)-1)))//'.itc'
!              write(stderr,'(A)')' File Not Found'
!              stop
!            endif
!            open(unit=Get_Unit(Unit_itc),file=adjustl(trim(line1(1:index(line1,tab)-1)))//'.itc')
!            do v=1,global%Ns
!              read(Unit_itc,*) blocks(b)%P%r(v)
!            enddo
!            read(Unit_itc,*) blocks(b)%P%v%x
!            read(Unit_itc,*) blocks(b)%P%v%y
!            read(Unit_itc,*) blocks(b)%P%v%z
!            read(Unit_itc,*) blocks(b)%P%p
!            read(Unit_itc,*) blocks(b)%P%d
!            read(Unit_itc,*) blocks(b)%P%g
!            close(Unit_itc)
!            do l=1,global%mesh_dims%Nl
!              block(1,l)%C%P = blocks(b)%P
!            enddo
!          endif
!        elseif (index(line1,tab//'e')>0) then
!          cycle
!        elseif (index(line1,tab//'v')>0) then
!          cycle
!        else
!          exit
!        endif
!      enddo
!    endif
!  enddo
!  close(global%dfile%IBM_Scratch(0,b,1))
!  ! storing the boundary and initial conditions in the scratch files
!  do l=1,global%mesh_dims%Nl
!    ! boundary conditions data
!    err = write_bc(array3D=block(1,l)%Fi%BC,unit=global%dfile%IBM_Scratch(2,b,l))
!    err = write_bc(array3D=block(1,l)%Fj%BC,unit=global%dfile%IBM_Scratch(2,b,l))
!    err = write_bc(array3D=block(1,l)%Fk%BC,unit=global%dfile%IBM_Scratch(2,b,l))
!    ! initial conditions data
!    write(global%dfile%IBM_Scratch(3,b,l),iostat=err)block(1,l)%C%Dt
!    err = write_primitive(array3D=block(1,l)%C%P,unit=global%dfile%IBM_Scratch(3,b,l))
!    ! rewinding scratch files
!    rewind(global%dfile%IBM_Scratch(2,b,l))
!    rewind(global%dfile%IBM_Scratch(3,b,l))
!  enddo
!enddo
!global%Nb = global%mesh_dims%Nb_tot
!return
!!---------------------------------------------------------------------------------------------------------------------------------
!endfunction load_icemcfd
  !> @}
endmodule Data_Type_File_ICEMCFD
