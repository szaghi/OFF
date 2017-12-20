module irprecision
!-----------------------------------------------------------------------------------------------------------------------------------
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
public
interface str
  module procedure str_R8P,str_I4P
endinterface
integer, parameter:: R8P  = selected_real_kind(15,307)  !< 15  digits, range \f$[10^{-307} , 10^{+307}  - 1]\f$; 64 bits.
integer, parameter:: I4P  = selected_int_kind(9)        !< Range \f$[-2^{31},+2^{31} - 1]\f$, 10 digits plus sign; 32 bits.
integer, parameter:: DR8P = 23                          !< Number of digits of output format FR8P.
integer, parameter:: DI4P = 11                          !< Number of digits of output format I4P.
character(10), parameter:: FR8P = '(E23.15E3)'          !< Output format for kind=R8P variable.
character(5),  parameter:: FI4P = '(I11)'               !< Output format for kind=I4P variable.
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  elemental function str_I4P(n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: n   !< Integer to be converted.
  character(DI4P)::          str !< Returned string containing input number plus padding zeros.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI4P) n        ! Casting of n to string.
  str = adjustl(trim(str)) ! Removing white spaces.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_I4P

  elemental function str_R8P(n) result(str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),    intent(IN):: n   !< Real to be converted.
  character(DR8P)::          str !< Returned string containing input number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FR8P) n        ! Casting of n to string.
  str = adjustl(trim(str)) ! Removing white spaces.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction str_R8P
endmodule irprecision

program make_cart_blocks
!-----------------------------------------------------------------------------------------------------------------------------------
USE irprecision
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type:: Type_Vector
  real(R8P):: x = 0._R8P !< Cartesian component in X direction.
  real(R8P):: y = 0._R8P !< Cartesian component in Y direction.
  real(R8P):: z = 0._R8P !< Cartesian component in Z direction.
endtype Type_Vector
type:: Type_Zone
  type(Type_Vector)::      domain(1:2) !< Zone domain bounding box (whole extent).
  real(R8P), allocatable:: p(:)        !< Zone primitive variables.
endtype Type_Zone
type(Type_Vector)::             domain(1:2)     !< Domain bounding box (whole extent).
real(R8P)::                     Dx,Dy,Dz        !< Space steps.
real(R8P)::                     x,y,z           !< Current coordinates.
integer(I4P)::                  Nbi = 0_I4P     !< Number of blocks in X direction.
integer(I4P)::                  Nbj = 0_I4P     !< Number of blocks in Y direction.
integer(I4P)::                  Nbk = 0_I4P     !< Number of blocks in Z direction.
integer(I4P)::                  Nb = 0_I4P      !< Number of blocks.
integer(I4P)::                  Ni = 0_I4P      !< Number of cells in X direction.
integer(I4P)::                  Nj = 0_I4P      !< Number of cells in Y direction.
integer(I4P)::                  Nk = 0_I4P      !< Number of cells in Z direction.
integer(I4P)::                  gc(1:6)         !< Number of ghost cells.
integer(I4P)::                  Ns = 0_I4P      !< Number of species.
character(3)::                  BC(1:6)         !< Boundary conditions tags.
integer(I4P)::                  Nzone = 0_I4P   !< Number of initial zones.
type(Type_Zone), allocatable::  zone(:)         !< Initial zones [1:Nzone].
integer(I4P)::                  Nca = 0_I4P     !< Number of command line arguments.
logical::                       is_file=.false. !< Flag for inquiring the presence of input file.
character(100)::                ifile           !< Input file name.
character(100)::                ofile           !< Output file name.
integer(I4P)::                  ui = 0_I4P      !< Input file unit.
integer(I4P)::                  uo = 0_I4P      !< Output file unit.
integer(I4P)::                  i,j,k,c,s       !< Counters.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
Nca = command_argument_count()
if (Nca==0) then
  write(stderr,'(A)') ' A valid file name of the input and output files must be provided as command line argument (in this order)'
  write(stderr,'(A)') ' No argument has been passed to command line'
  write(stderr,'(A)') ' Correct use is:'
  write(stderr,*)
  write(stderr,'(A)') ' make_cart_blocks "valid_input_file_name" "valid_output_file_name"'
  stop
else
  call get_command_argument(1, ifile)
  call get_command_argument(2, ofile)
endif
inquire(file=trim(ifile),exist=is_file)
if (.not.is_file) then
  write(stderr,'(A)') ' Input file "'//trim(ifile)//'" not found'
  stop
endif
open(unit=Get_Unit(ui),file=trim(ifile),action="READ")
read(ui,*)domain(1)%x,domain(1)%y,domain(1)%z
read(ui,*)domain(2)%x,domain(2)%y,domain(2)%z
read(ui,*)Nbi,Nbj,Nbk
read(ui,*)(gc(c),c=1,6)
read(ui,*)Ni,Nj,Nk
read(ui,*)Ns
read(ui,*)(BC(c),c=1,6)
read(ui,*)
read(ui,*)Nzone
allocate(zone(1:Nzone))
do c=1,Nzone
  allocate(zone(c)%p(1:Ns+6))
  read(ui,*)
  read(ui,*)zone(c)%domain(1)%x,zone(c)%domain(1)%y,zone(c)%domain(1)%z
  read(ui,*)zone(c)%domain(2)%x,zone(c)%domain(2)%y,zone(c)%domain(2)%z
  do s=1,Ns+6
    read(ui,*)zone(c)%p(s)
  enddo
enddo
close(ui)
Nb = Nbi*Nbj*Nbk
Dx = (domain(2)%x-domain(1)%x)/real(Nbi,R8P)
Dy = (domain(2)%y-domain(1)%y)/real(Nbj,R8P)
Dz = (domain(2)%z-domain(1)%z)/real(Nbk,R8P)

open(unit=Get_Unit(uo),file=trim(ofile),action="WRITE")
write(uo,'(A)')repeat('-',132)
write(uo,'(A)')repeat(' ',34)//'BLOCKS BOUNDS'
write(uo,'(A)')trim(str(Ns))//repeat(' ',60)//'Ns = number of species'
write(uo,'(A)')trim(str(Nb))//repeat(' ',60)//'Nb = number of blocks'
do k=1,Nbk
  do j=1,Nbj
    do i=1,Nbi
      x = real((i-1),R8P)*Dx ; y = real((j-1),R8P)*Dy ; z = real((k-1),R8P)*Dz
      write(uo,*)
      write(uo,'(7A)')(trim(str(gc(c)))//' ',c=1,6),repeat(' ',49)//'gc(1:6)        = ghost cells of face boundary i,j,k'
      write(uo,'(A)') trim(str(Ni))//' '//trim(str(Nj))//' '//trim(str(Nk))//&
                      repeat(' ',54)//'Ni,Nj,Nk       = number of cells in i,j,k directions'
      write(uo,'(A)') trim(str(x))//' '//trim(str(y))//' '//trim(str(z))//' '//&
                      repeat(' ',45)//'xmin,ymin,zmin = minimum abscissa coordinates of block nodes'
      write(uo,'(A)') trim(str(x+Dx))//' '//trim(str(y+Dy))//' '//trim(str(z+Dz))//' '//&
                      repeat(' ',45)//'xmax,ymax,zmax = maximum abscissa coordinates of block nodes'
      if (i==1.AND.Nbi==1) then
        write(uo,'(A)')"'"//BC(1)//"' 0"//                 repeat(' ',54)//'face(1)        = left  i boundary condition'
        write(uo,'(A)')"'"//BC(2)//"' 0"//                 repeat(' ',54)//'face(2)        = right i boundary condition'
      elseif (i==1) then
        write(uo,'(A)')"'"//BC(1)//"' 0"//                 repeat(' ',54)//'face(1)        = left  i boundary condition'
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i+1,j,k)))//repeat(' ',54)//'face(2)        = right i boundary condition'
      elseif (i==Nbi) then
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i-1,j,k)))//repeat(' ',54)//'face(1)        = left  i boundary condition'
        write(uo,'(A)')"'"//BC(2)//"' 0"//                 repeat(' ',54)//'face(2)        = right i boundary condition'
      else
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i-1,j,k)))//repeat(' ',54)//'face(1)        = left  i boundary condition'
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i+1,j,k)))//repeat(' ',54)//'face(2)        = right i boundary condition'
      endif
      if (j==1.AND.Nbj==1) then
        write(uo,'(A)')"'"//BC(3)//"' 0"//                 repeat(' ',54)//'face(3)        = left  j boundary condition'
        write(uo,'(A)')"'"//BC(4)//"' 0"//                 repeat(' ',54)//'face(4)        = right j boundary condition'
      elseif (j==1) then
        write(uo,'(A)')"'"//BC(3)//"' 0"//                 repeat(' ',54)//'face(3)        = left  j boundary condition'
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i,j+1,k)))//repeat(' ',54)//'face(4)        = right j boundary condition'
      elseif (j==Nbj) then
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i,j-1,k)))//repeat(' ',54)//'face(3)        = left  j boundary condition'
        write(uo,'(A)')"'"//BC(4)//"' 0"//                 repeat(' ',54)//'face(4)        = right j boundary condition'
      else
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i,j-1,k)))//repeat(' ',54)//'face(3)        = left  j boundary condition'
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i,j+1,k)))//repeat(' ',54)//'face(4)        = right j boundary condition'
      endif
      if (k==1.AND.Nbk==1) then
        write(uo,'(A)')"'"//BC(5)//"' 0"//                 repeat(' ',54)//'face(5)        = left  k boundary condition'
        write(uo,'(A)')"'"//BC(6)//"' 0"//                 repeat(' ',54)//'face(6)        = right k boundary condition'
      elseif (k==1) then
        write(uo,'(A)')"'"//BC(5)//"' 0"//                 repeat(' ',54)//'face(5)        = left  k boundary condition'
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i,j,k+1)))//repeat(' ',54)//'face(6)        = right k boundary condition'
      elseif (k==Nbk) then
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i,j,k-1)))//repeat(' ',54)//'face(5)        = left  k boundary condition'
        write(uo,'(A)')"'"//BC(6)//"' 0"//                 repeat(' ',54)//'face(6)        = right k boundary condition'
      else
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i,j,k-1)))//repeat(' ',54)//'face(5)        = left  k boundary condition'
        write(uo,'(A)')"'ADJ' "//trim(str(bijk(i,j,k+1)))//repeat(' ',54)//'face(6)        = right k boundary condition'
      endif
      do c=1,Nzone
        if ((x+Dx*.5_R8P>zone(c)%domain(1)%x).AND.(x+Dx*.5_R8P<zone(c)%domain(2)%x).AND.&
            (y+Dy*.5_R8P>zone(c)%domain(1)%y).AND.(y+Dy*.5_R8P<zone(c)%domain(2)%y).AND.&
            (z+Dz*.5_R8P>zone(c)%domain(1)%z).AND.(z+Dz*.5_R8P<zone(c)%domain(2)%z)) then
          do s=1,Ns+6
            write(uo,'(A)')trim(str(zone(c)%p(s)))
          enddo
          exit
        endif
      enddo
    enddo
  enddo
enddo
write(uo,'(A)')repeat('-',132)
close(uo)
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  integer function Get_Unit(Free_Unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer, intent(OUT), optional:: Free_Unit !< Free logic unit.
  integer::                        n1        !< Counter.
  integer::                        ios       !< Inquiring flag.
  logical::                        lopen     !< Inquiring flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Get_Unit = -1
  n1=1
  do
    if ((n1/=stdout).AND.(n1/=stderr)) then
      inquire (unit=n1,opened=lopen,iostat=ios)
      if (ios==0) then
        if (.NOT.lopen) then
          Get_Unit = n1 ; if (present(Free_Unit)) Free_Unit = Get_Unit
          return
        endif
      endif
    endif
    n1=n1+1
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Get_Unit

  elemental function bijk(i,j,k) result(b)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: i,j,k !< Block coordinates.
  integer(I4P)::             b     !< Block index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  b = (k-1)*Nbj*Nbi + (j-1)*Nbi + i
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction bijk
endprogram make_cart_blocks
