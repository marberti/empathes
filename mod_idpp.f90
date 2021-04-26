! Copyright (C) 2020-2021  Marco Bertini
!
! This file is part of Empathes.
!
! Empathes is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Empathes is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Empathes.  If not, see <https://www.gnu.org/licenses/>.

module idpp

  !==================================================================
  !   IDPP Module
  !==================================================================
  ! Procedures contained in this module are meant to provide
  ! a better esteem of initial geometries than linear interpolation.
  !
  ! References:
  !    * Smidstrup 2014, doi: 10.1063/1.4878664
  !==================================================================

  use utility
  use geometry

  implicit none
  save
  private

  ! protected variables -----------------------------------
  public    :: flag_idpp,         &
               idpp_tol,          &
               idpp_energy,       &
               idpp_forces
  protected :: flag_idpp,         &
               idpp_tol,          &
               idpp_energy,       &
               idpp_forces
  ! public procedures -------------------------------------
  public    :: set_idpp,          &
               set_idpp_tol,      &
               init_idpp,         &
               compute_idpp_enfo, &
               driver_idpp

  !--------------------------------------------------------
  logical                                :: flag_idpp      = .false.
  logical                                :: flag_init_idpp = .false.
  integer                                :: d_idpp_len
  real(DBL)                              :: idpp_tol       = 1.0E-3_DBL
    ! default idpp convergence threshold ---^
  real(DBL), allocatable, dimension(:)   :: idpp_energy
  real(DBL), allocatable, dimension(:,:) :: idpp_forces
  real(DBL), allocatable, dimension(:,:) :: d_idpp ! contains
    ! the interpolated pair distances (0:images+1,d_idpp_len)

contains

!====================================================================
! Public
!====================================================================

subroutine set_idpp(flag)

  logical, intent(IN) :: flag

  logical, save       :: first_call = .true.

  if (first_call.eqv..false.) then
    call error("set_idpp: subroutine called more than once");
  end if

  flag_idpp = flag

  first_call = .false.

end subroutine set_idpp

!====================================================================

subroutine set_idpp_tol(str)

  character(*), intent(IN) :: str

  logical, save            :: first_call = .true.

  if (first_call.eqv..false.) then
    call error("set_idpp_tol: subroutine called more than once")
  end if

  if (isreal(trim(adjustl(str)))) then
    read(str,*) idpp_tol
  else
    call error("set_idpp_tol: argument must be a real")
  end if

  first_call = .false.

end subroutine set_idpp_tol

!====================================================================

subroutine init_idpp()

  !--------------------------------------------------------
  ! If flag_idpp is setted to .true., this subroutine
  ! initializes the distances array d_idpp used by
  ! other procedures of this module.
  !--------------------------------------------------------

  real(DBL), allocatable, dimension(:) :: delta
  integer                              :: i
  integer                              :: j
  integer                              :: k
  integer                              :: atoms
  integer                              :: indx
  integer                              :: err_n
  character(120)                       :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_idpp.eqv..false.) then
    return
  end if

  if (flag_init_images.eqv..false.) then
    call error("init_idpp: images are not initialized")
  end if

  if (flag_init_idpp.eqv..true.) then
    call error("init_idpp: idpp module already initialized")
  end if

  ! allocate the arrays -----------------------------------
  atoms      = geom_len/3
  d_idpp_len = triang_numb(atoms-1)

  allocate(d_idpp(0:image_n+1,d_idpp_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_idpp: "//trim(err_msg))
  end if

  allocate(delta(d_idpp_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_idpp: "//trim(err_msg))
  end if

  allocate(idpp_energy(0:image_n+1),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_idpp: "//trim(err_msg))
  end if

  allocate(idpp_forces(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_idpp: "//trim(err_msg))
  end if

  ! init d_idpp array -------------------------------------
  ! start and end geometries
  do k=0, image_n+1, +(image_n+1)
    do i=1, atoms
      do j=i+1, atoms
        indx           = get_d_idpp_index(i,j)
        d_idpp(k,indx) = atomic_distance(k,i,j)
      end do
    end do
  end do

  ! image geometries
  delta = d_idpp(image_n+1,:)-d_idpp(0,:)
  do k=1, image_n
    d_idpp(k,:) = d_idpp(0,:)+((k/real(image_n+1,DBL))*delta)
  end do

  ! init idpp_energy for start and end points -------------
  call get_idpp_energy(0)
  call get_idpp_energy(image_n+1)
  
  ! done --------------------------------------------------
  deallocate(delta,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_idpp: "//trim(err_msg))
  end if

  flag_init_idpp = .true.

end subroutine init_idpp

!====================================================================

subroutine compute_idpp_enfo()

  !--------------------------------------------------------
  ! Computes idpp energies and forces for all the images.
  ! Results are stored in idpp_energy and idpp_forces.
  !--------------------------------------------------------

  integer :: k

  ! preliminary checks ------------------------------------
  if (flag_init_idpp.eqv..false.) then
    call error("compute_idpp_enfo: idpp module not initialized")
  end if

  do k=1, image_n
    call get_idpp_energy(k)
    call get_idpp_forces(k)
  end do

end subroutine compute_idpp_enfo

!====================================================================
! Private
!====================================================================

subroutine get_idpp_energy(k)

  !--------------------------------------------------------
  ! Computes the idpp energy for the image k.
  ! Result is stored in idpp_energy array.
  !--------------------------------------------------------

  integer, intent(IN) :: k

  integer             :: i
  integer             :: j
  integer             :: atoms
  integer             :: indx
  character(8)        :: i_str
  real(DBL)           :: res
  real(DBL)           :: w
  real(DBL)           :: d

  ! preliminary checks ------------------------------------
  if ((k<0).or.(k>image_n+1)) then
    write(i_str,'(I8)') k
    i_str = adjustl(i_str)
    call error("get_idpp_energy: image "//trim(i_str)//" out of range")
  end if

  ! working section ---------------------------------------
  atoms = geom_len/3

  res = 0.0_DBL
  do i = 1, atoms
    do j = i+1, atoms
      indx = get_d_idpp_index(i,j)
      d    = atomic_distance(k,i,j)
      w    = 1.0_DBL/(d**4)
      res  = res + (w*((d_idpp(k,indx) - d)**2))
    end do
  end do

  idpp_energy(k) = res

end subroutine get_idpp_energy

!====================================================================

subroutine get_idpp_forces(k)

  !--------------------------------------------------------
  ! Computes the idpp forces for the image k.
  ! Results are stored in idpp_forces matrix.
  !--------------------------------------------------------

  integer, intent(IN) :: k     ! image number

  character(8)        :: i_str
  integer             :: atoms ! number of atoms
  integer             :: i     ! atom label
  integer             :: j     ! atom label
  integer             :: c     ! dummy coordinate, used for x, y and z
  integer             :: ci    ! dummy coordinate, used for xi, yi and zi
  integer             :: cj    ! dummy coordinate, used for xj, yj and zj
  real(DBL)           :: s     ! derivative of sij
  real(DBL)           :: d     ! atomic distance between i and j
  real(DBL)           :: dc    ! intermediate term
  real(DBL)           :: t     ! intermediate term

  ! preliminary checks ------------------------------------
  if ((k<1).or.(k>image_n)) then
    write(i_str,'(I8)') k
    i_str = adjustl(i_str)
    call error("get_idpp_forces: image "//trim(i_str)//" out of range")
  end if

  ! working section ---------------------------------------
  atoms            = geom_len/3
  idpp_forces(k,:) = 0.0_DBL

  do i=1, atoms
    do j=i+1, atoms
      d = atomic_distance(k,i,j)
      t = (d_idpp(k,get_d_idpp_index(i,j))-d)
      
      do c=1, 3
        ci = 3*(i-1)+c
        cj = 3*(j-1)+c
        dc = image_geom(k,ci)-image_geom(k,cj)
        s  = ((-4.0_DBL*dc*(t**2))/(d**6))+((-2.0_DBL*dc*t)/(d**5))

        idpp_forces(k,ci) = idpp_forces(k,ci)+s
        idpp_forces(k,cj) = idpp_forces(k,cj)-s
      end do
    end do
  end do

  ! convert from derivative to force ----------------------
  idpp_forces(k,:) = idpp_forces(k,:) * (-1)

end subroutine get_idpp_forces

!====================================================================

integer function get_d_idpp_index(r,c)

  !--------------------------------------------------------
  ! Given two atoms r and c,
  ! it returns the relative d_idpp array index.
  !--------------------------------------------------------

  integer, intent(IN) :: r
  integer, intent(IN) :: c

  integer             :: atoms
  integer             :: a

  atoms = geom_len/3

  ! preliminary checks ------------------------------------
  if (r>=c) then
    call error("get_d_idpp_index: row index >= column index")
  else if (r<=0) then
    call error("get_d_idpp_index: row index <= 0")
  else if (r>=atoms) then
    call error("get_d_idpp_index: row index >= atoms")
  else if (c<=0) then
    call error("get_d_idpp_index: column index <= 0")
  else if (c>atoms) then
    call error("get_d_idpp_index: column index > atoms")
  end if

  ! get the index -----------------------------------------
  if (r==1) then
    a = 0
  else
    a = (r-1)*atoms - triang_numb(r-1)
  end if

  get_d_idpp_index = c - r + a

end function get_d_idpp_index

!====================================================================

real(DBL) function atomic_distance(k,i,j)

  !--------------------------------------------------------
  ! Taken the image k and given two atoms i and j,
  ! this function computes the distance between them.
  ! If i and j are the same atom, distance = 0.0 is returned.
  !--------------------------------------------------------

  integer, intent(IN) :: k
  integer, intent(IN) :: i
  integer, intent(IN) :: j

  integer             :: atoms
  integer             :: coor
  real(DBL)           :: tot
  character(8)        :: str

  atoms = geom_len/3

  ! preliminary checks ------------------------------------
  if ((k<0).or.(k>image_n+1)) then
    write(str,'(I8)') k
    str = adjustl(str)
    call error("atomic_distance: image "//trim(str)//" out of range")
  else if ((i<1).or.(i>atoms)) then
    write(str,'(I8)') i
    str = adjustl(str)
    call error("atomic_distance: atom "//trim(str)//" out of range")
  else if ((j<1).or.(j>atoms)) then
    write(str,'(I8)') j
    str = adjustl(str)
    call error("atomic_distance: atom "//trim(str)//" out of range")
  end if

  ! compute the distance ----------------------------------
  if (i==j) then
    atomic_distance = 0.0_DBL
    return
  end if

  tot = 0.0_DBL
  do coor=1, 3
    tot = tot + (image_geom(k,3*(i-1)+coor)-image_geom(k,3*(j-1)+coor))**2
  end do

  atomic_distance = sqrt(tot)

end function atomic_distance

!====================================================================
! Driver
!====================================================================

subroutine driver_idpp()

  !--------------------------------------------------------
  ! For testing purposes only.
  !--------------------------------------------------------

  integer :: atoms
  integer :: indx
  integer :: i
  integer :: j
  integer :: k

  ! preliminary checks ------------------------------------
  if (flag_init_idpp.eqv..false.) then
    call error("driver_idpp: idpp module not initialized")
  end if

  write(FILEOUT,*) "DRV driver_idpp: start"
  atoms = geom_len/3

  ! print d_idpp matrix -----------------------------------
  write(FILEOUT,*) "DRV d_idpp matrix"
  do k=0, image_n+1
    write(FILEOUT,*) "DRV Image ", k
    do i=1, atoms
      do j=i+1, atoms
        indx = get_d_idpp_index(i,j)
        write(FILEOUT,'(5X,I4,2X,F15.6)') indx, d_idpp(k,indx)
      end do
    end do
  end do

  ! energy and forces -------------------------------------
  call compute_idpp_enfo()

  write(FILEOUT,*) "DRV energies and forces"
  do k=0, image_n+1
    write(FILEOUT,*) "DRV Image ", k
    write(FILEOUT,'(5X,"energy: ",5X,F15.6)') idpp_energy(k)
    
    if ((k==0).or.(k==image_n+1)) then
      cycle
    end if

    do i=1, geom_len
      j = ((i-1)/3)+1
      write(FILEOUT,'(5X,"force ",I4)',advance="no") j
      
      j = mod(i,3)
      select case (j)
      case (0)
        write(FILEOUT,'("z: ")',advance="no")
      case (1)
        write(FILEOUT,'("x: ")',advance="no")
      case (2)
        write(FILEOUT,'("y: ")',advance="no")
      end select

      write(FILEOUT,'(F18.9)') idpp_forces(k,i)
    end do
  end do

  ! ending execution --------------------------------------
  write(FILEOUT,*) "DRV driver_idpp: driver completed; ending execution"
  call end_main_exec() ! driver_idpp

end subroutine driver_idpp

!====================================================================

end module idpp

