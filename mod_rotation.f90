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

module rotation

  use utility

  implicit none
  save
  private

  ! protected variables -----------------------------------
  public    :: flag_rotation
  protected :: flag_rotation
  ! public procedures -------------------------------------
  public    :: set_rotation,   &
               rotate_geometry

  !--------------------------------------------------------
  logical :: flag_rotation = .false.

contains

!====================================================================
! Public
!====================================================================

subroutine set_rotation(flag)

  logical, intent(IN) :: flag

  logical, save       :: first_call = .true.

  if (first_call.eqv..false.) then
    call error("set_rotation: subroutine called more than once")
  end if

  flag_rotation = flag

  first_call = .false.

end subroutine set_rotation

!====================================================================

subroutine rotate_geometry(geom_in,geom_out)

  real(DBL), dimension(:),   intent(IN)  :: geom_in
  real(DBL), dimension(:),   intent(OUT) :: geom_out

  real(DBL), dimension(3,3), parameter   :: rx  =&
    &reshape([1.0_DBL,0.0_DBL,0.0_DBL,0.0_DBL,0.0_DBL,-1.0_DBL,0.0_DBL,1.0_DBL,0.0_DBL],[3,3])
  real(DBL), dimension(3,3), parameter   :: ry  =&
    &reshape([0.0_DBL,0.0_DBL,1.0_DBL,0.0_DBL,1.0_DBL,0.0_DBL,-1.0_DBL,0.0_DBL,0.0_DBL],[3,3])
  real(DBL), dimension(3,3), parameter   :: rzp =&
    &reshape([-1.0_DBL,0.0_DBL,0.0_DBL,0.0_DBL,-1.0_DBL,0.0_DBL,0.0_DBL,0.0_DBL,1.0_DBL],[3,3])
  real(DBL), dimension(3,3), parameter   :: rxp =&
    &reshape([1.0_DBL,0.0_DBL,0.0_DBL,0.0_DBL,-1.0_DBL,0.0_DBL,0.0_DBL,0.0_DBL,-1.0_DBL],[3,3])
  integer                                :: atoms
  integer                                :: geom_len
  integer                                :: i
  integer                                :: j
  integer                                :: indx
  real(DBL), allocatable, dimension(:,:) :: geom     ! contains input geom
  real(DBL), allocatable, dimension(:,:) :: geom_rot ! contains output geom
  real(DBL)                              :: aa
  real(DBL)                              :: bb
  real(DBL)                              :: cc
  real(DBL)                              :: ad
  real(DBL)                              :: bd
  real(DBL)                              :: cd
  real(DBL)                              :: d
  real(DBL), dimension(3,3)              :: rot
  real(DBL)                              :: xpr
  real(DBL)                              :: ypr
  real(DBL)                              :: theta
  integer                                :: err_n
  character(120)                         :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_rotation.eqv..false.) then
    call error("rotate_geometry: called with false flag_rotation")
  end if

  geom_len = size(geom_in,1)

  if (geom_len/=size(geom_out,1)) then
    call error("rotate_geometry: inconsistent argument lengths")
  end if

  if (mod(geom_len,3)/=0) then
    call error("rotate_geometry: argument lengths must be multiples of 3")
  end if

  atoms = geom_len/3

  !TODO implement rotation for atoms < 3
  if (atoms<3) then
    write(FILEOUT,*) "WAR rotate_geometry: not implemented yet for less than 3 atoms"
    return
  end if

  ! allocation section ------------------------------------
  allocate(geom(atoms,3),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("rotate_geometry: "//trim(err_msg))
  end if

  allocate(geom_rot(atoms,3),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("rotate_geometry: "//trim(err_msg))
  end if

  ! initialization section --------------------------------
  do i=1, atoms
    do j=1, 3
      indx      = 3*(i-1)+j
      geom(i,j) = geom_in(indx)
    end do
  end do

  ! working section ---------------------------------------
  ad = geom(1,1)-geom(2,1)
  bd = geom(1,2)-geom(2,2)
  cd = geom(1,3)-geom(2,3)
  d  = sqrt(ad*ad+bd*bd+cd*cd)

  if ((geom(1,3)==geom(2,3)).and.(geom(1,3)==geom(3,3))) then ! three points in XY plane
    rot = 1.0_DBL
  else if ((geom(1,2)==geom(2,2)).and.(geom(1,2)==geom(3,2))) then ! three points in XZ plane
    rot = rx
  else if ((geom(1,1)==geom(2,1)).and.(geom(1,1)==geom(3,1))) then ! three points in YZ plane
    rot = ry
  else
    call ppptp(geom(1:3,:),aa,bb,cc) ! find the plane passing for (at1,at2,at3)
    call rot_mat(rot,aa,bb,cc)
  end if

  ! translate at1 to (0,0,0)
  do i=1, 3
    geom(:,i) = geom(:,i)-geom(1,i) 
  end do

  ! put (at1,at2,at3) on XY plane
  geom_rot = matmul(geom,transpose(rot))

  theta = sign(1.0_DBL,geom_rot(2,1)) * acos(geom_rot(2,1)/d)

  ! put at2 on X axis
  do i=2, atoms
    xpr           = geom_rot(i,1)
    ypr           = geom_rot(i,2)
    geom_rot(i,1) = cos(theta)*xpr - sin(theta)*ypr
    geom_rot(i,2) = sin(theta)*xpr + cos(theta)*ypr
  end do

  ! put x(at2) on negative X axis 
  if (geom_rot(2,1)>0.0_DBL) then
    geom_rot = matmul(geom_rot,rzp)
  end if

  ! make y(at3) positive
  if (geom_rot(3,2)<0.0_DBL) then
    geom_rot = matmul(geom_rot,rxp)
  end if

  ! copy result -------------------------------------------
  do i=1, atoms
    do j=1, 3
      indx           = 3*(i-1)+j
      geom_out(indx) = geom_rot(i,j)
    end do
  end do

  ! deallocation section ----------------------------------
  deallocate(geom,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("rotate_geometry: "//trim(err_msg))
  end if

  deallocate(geom_rot,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("rotate_geometry: "//trim(err_msg))
  end if

end subroutine rotate_geometry

!====================================================================
! Private
!====================================================================

subroutine ppptp(p,aa,bb,cc)

  real(DBL), dimension(3,3), intent(IN)  :: p
  real(DBL),                 intent(OUT) :: aa
  real(DBL),                 intent(OUT) :: bb
  real(DBL),                 intent(OUT) :: cc

  real(DBL), dimension(5)                :: coef
  real(DBL)                              :: m
  real(DBL)                              :: mp
  real(DBL)                              :: n
  real(DBL)                              :: np
  real(DBL)                              :: o
  real(DBL)                              :: op
  real(DBL)                              :: A

  ! parametric equation of plane
  m  = p(2,1) - p(1,1)
  n  = p(2,2) - p(1,2)
  o  = p(2,3) - p(1,3)
  mp = p(3,1) - p(1,1)
  np = p(3,2) - p(1,2)
  op = p(3,3) - p(1,3)

  !TODO check denominators before divisions
  ! switch to cartesian equation
  a       = (n*op-np*o) / (m*np-mp*n)
  coef(1) = 1.0_DBL
  coef(2) = a*a
  coef(3) = (mp*mp*a*a) / (np*np)
  coef(4) = (op*op) / (np*np)
  coef(5) = (2.0_DBL*mp*op*a) / (np*np)

  cc = sqrt(1.0_DBL / sum(coef))
  aa = A*cc
  bb = -1.0_DBL*(mp*aa+op*cc) / np

end subroutine ppptp

!====================================================================

subroutine rot_mat(rot,aa,bb,cc)

  real(DBL), dimension(3,3), intent(OUT) :: rot
  real(DBL),                 intent(IN)  :: aa
  real(DBL),                 intent(IN)  :: bb
  real(DBL),                 intent(IN)  :: cc

  real(DBL), parameter                   :: PI  = 3.14159265358979323846264338327950288_DBL
  real(DBL)                              :: phi
  real(DBL)                              :: tet
  real(DBL)                              :: bbp

  phi = -PI/2.0_DBL
  tet = -PI/2.0_DBL
  
  if (bb/=0.0_DBL) then
    phi = atan(-1.0_DBL*aa/bb)
  end if
  
  bbp = -1.0_DBL*aa*sin(phi) + bb*cos(phi)
  
  if (cc/=0.0_DBL) then
    tet = atan(-1.0_DBL*bbp/cc)
  end if

  rot(1,1) = cos(phi)
  rot(1,2) = sin(phi)
  rot(1,3) = 0.0_DBL
  rot(2,1) = -1.0_DBL*sin(phi) * cos(tet)
  rot(2,2) = cos(phi) * cos(tet)
  rot(2,3) = sin(tet)
  rot(3,1) = sin(phi) * sin(tet)
  rot(3,2) = -1.0_DBL*cos(phi) * sin(tet)
  rot(3,3) = cos(tet)

end subroutine rot_mat

!====================================================================

end module rotation

